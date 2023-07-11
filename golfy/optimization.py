import random
import time

import numpy as np

from .solution import Solution
from .types import SwapCandidateList, ReplicateToNeighborDict
from .validity import count_violations, violations_per_replicate
from .util import pairs_to_dict


def find_violating_peptides(
    s: Solution,
) -> tuple[SwapCandidateList, ReplicateToNeighborDict]:
    replicate_to_pool_to_peptides = s.assignments
    # treat invalid pairs as if they've already been neighbors in a previous round
    peptide_to_neighbors = pairs_to_dict(s.invalid_neighbors)
    needs_swap = []
    replicate_to_neighbor_dict = {}
    for replicate_idx, pool_to_peptides in replicate_to_pool_to_peptides.items():
        for pool_idx, peptides in pool_to_peptides.items():
            for p1 in peptides:
                for p2 in peptides:
                    if p1 != p2:
                        if p2 in peptide_to_neighbors[p1]:
                            needs_swap.append((replicate_idx, pool_idx, p2))
                        else:
                            peptide_to_neighbors[p1].add(p2)
        # neighbor constaints at the end of this replicate
        replicate_to_neighbor_dict[replicate_idx] = {
            peptide: neighbors.copy()
            for (peptide, neighbors) in peptide_to_neighbors.items()
        }
    return needs_swap, replicate_to_neighbor_dict


def improve_solution(s: Solution, verbose: bool = False):
    replicate_to_pool_to_peptides = s.assignments

    needs_swap, replicate_to_neighbor_dict = find_violating_peptides(s)

    random.shuffle(needs_swap)

    swapped_pools = set()
    swapped_peptides = set()
    empty_pools_per_replicate = {
        replicate_idx: {
            pool_idx
            for (pool_idx, pool_peptides) in pool_to_peptides.items()
            if len(pool_peptides) == 0
        }
        for (replicate_idx, pool_to_peptides) in replicate_to_pool_to_peptides.items()
    }
    for replicate_idx, pool_idx_a, peptide_a in needs_swap:
        if pool_idx_a in swapped_pools or peptide_a in swapped_peptides:
            continue

        pool_to_peptides = replicate_to_pool_to_peptides[replicate_idx]
        pool_a = pool_to_peptides[pool_idx_a]
        empty_pools = empty_pools_per_replicate[replicate_idx]
        # if a pool is empty, just move the offending peptide there
        if empty_pools:
            pool_idx_b = random.choice(list(empty_pools))
            s.move_peptide(replicate_idx, pool_idx_a, peptide_a, pool_idx_b)
            swapped_pools.add(pool_idx_a)
            swapped_pools.add(pool_idx_b)
            swapped_peptides.add(peptide_a)
            empty_pools.remove(pool_idx_b)
        else:
            other_peptides = []
            peptide_to_pool_idx = {}
            for pool_idx_i, pool_peptides_i in pool_to_peptides.items():
                if pool_idx_i != pool_idx_a and pool_idx_i not in swapped_pools:
                    all_peptides_ok = True
                    for p in pool_peptides_i:
                        peptide_to_pool_idx[p] = pool_idx_i
                        all_peptides_ok = all_peptides_ok and (
                            peptide_a
                            not in replicate_to_neighbor_dict[replicate_idx][p]
                        )
                    if all_peptides_ok:
                        other_peptides.extend(
                            [p for p in pool_peptides_i if p not in swapped_peptides]
                        )

            if len(other_peptides) == 0:
                if verbose:
                    print(
                        "Not able to find a valid peptide to swap with for (%s, %s, %s)"
                        % (replicate_idx, pool_idx_a, peptide_a)
                    )
                continue
            peptide_b = random.choice(other_peptides)
            pool_idx_b = peptide_to_pool_idx[peptide_b]
            assert peptide_a != peptide_b
            assert pool_idx_a != pool_idx_b
            pool_b = pool_to_peptides[pool_idx_b]
            if (
                len(pool_a) > 1
                and len(pool_b) < s.max_peptides_per_pool
                and random.choice([False, True])
            ):
                # just move peptide a to the pool with fewer than max peptides
                if verbose:
                    print(
                        "Moving peptide %d from pool %d to pool %d"
                        % (peptide_a, pool_idx_a, pool_idx_b)
                    )
                s.move_peptide(replicate_idx, pool_idx_a, peptide_a, pool_idx_b)
                swapped_pools.add(pool_idx_a)
                swapped_pools.add(pool_idx_b)
                swapped_peptides.add(peptide_a)
            else:
                if verbose:
                    print("Before swap")
                    print(
                        "pool",
                        pool_idx_a,
                        "peptide",
                        peptide_a,
                        pool_to_peptides[pool_idx_a],
                    )
                    print(
                        "pool",
                        pool_idx_b,
                        "peptide",
                        peptide_b,
                        pool_to_peptides[pool_idx_b],
                    )

                # actually swap them
                s.swap_peptides(
                    replicate_idx, pool_idx_a, peptide_a, pool_idx_b, peptide_b
                )

                if verbose:
                    print("After")
                    print(
                        "pool",
                        pool_idx_a,
                        "peptide",
                        peptide_a,
                        pool_to_peptides[pool_idx_a],
                    )
                    print(
                        "pool",
                        pool_idx_b,
                        "peptide",
                        peptide_b,
                        pool_to_peptides[pool_idx_b],
                    )

                swapped_pools.add(pool_idx_a)
                swapped_pools.add(pool_idx_b)
                swapped_peptides.add(peptide_a)
                swapped_peptides.add(peptide_b)


def optimize(
    s: Solution,
    max_iters: int = 2000,
    verbose: bool = False,
    add_pool_if_stuck: bool = True,
    return_history: bool = False,
) -> bool:
    """
    Iteratively update solution by randomly swapping a violating peptide with a random other peptide

    Args
    ----
    solution
        Initial solution which will be modified in-place

    max_iters
        Maximum number of swaps to consider performing

    verbose
        print number of violations and pools for each iteration

    add_pool_if_stuck
        If no improvements have been made for 10 iters, add a pool
        to the last replicate

    return_history
        Return array of constraint validation counts per iteration


    Returns True if non-violating solution found, False if solution still has violations after
    max_iters
    """
    replicate_to_violation_count = violations_per_replicate(s)
    old_num_violations = sum(replicate_to_violation_count.values())

    history = [old_num_violations]
    if verbose:
        print("Initial solution has %s violations" % (old_num_violations,))

    if old_num_violations / s.num_replicates > 1000:
        # before first iteration, add a lot of empty pools to any replicates with
        # more than 1000 violations
        for replicate_idx, violation_count in replicate_to_violation_count.items():
            if violation_count > 1000:
                num_new_pools = int(np.ceil(violation_count / 1000))
                if verbose:
                    print(
                        "-- replicate %d has %d violations, adding %d pools"
                        % (replicate_idx + 1, violation_count, num_new_pools)
                    )
                for _ in range(num_new_pools):
                    s.add_empty_pool(replicate_idx)

    num_iters_without_improvement = 0
    for i in range(max_iters):
        t0 = time.time()
        improve_solution(s)

        replicate_to_violation_count = violations_per_replicate(s)

        new_num_violations = sum(replicate_to_violation_count.values())

        history.append(new_num_violations)

        if verbose:
            print(
                "%d) %d -> %d violations (%d pools, %0.2fs)"
                % (
                    i + 1,
                    old_num_violations,
                    new_num_violations,
                    s.num_pools(),
                    time.time() - t0,
                )
            )

        if old_num_violations <= new_num_violations:
            num_iters_without_improvement += 1
        else:
            num_iters_without_improvement = 0

        old_num_violations = new_num_violations

        if new_num_violations == 0:
            if verbose:
                print("Found valid solution after %d swaps" % (i + 1,))
            break

        if num_iters_without_improvement >= 3 and add_pool_if_stuck:
            replicates_in_need = [
                r for (r, v) in replicate_to_violation_count.items() if v > 0
            ]
            assert len(replicates_in_need) > 0
            for replicate_idx in replicates_in_need:
                if verbose:
                    print(
                        "-- adding pool %d to replicate %d"
                        % (len(s.assignments[replicate_idx]), replicate_idx + 1)
                    )
                s.add_empty_pool(replicate_idx)
            num_iters_without_improvement = 0

    result = old_num_violations == 0

    # just in case we ended up with any empty pools, remove them from the solution
    s.remove_empty_pools()

    if return_history:
        return result, np.array(history)
    else:
        return result
