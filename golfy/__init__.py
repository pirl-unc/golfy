from typing import Iterable
from collections import defaultdict

import numpy as np


# type aliases
Solution = dict[int, dict[int, Iterable[int]]]
ConstraintList = Iterable[tuple[int, int]]


def random_init(
    num_peptides: int = 100, peptides_per_pool: int = 5, num_replicates: int = 3
) -> Solution:
    num_pools = int(np.ceil(num_peptides / peptides_per_pool))
    replicate_to_pool_to_peptides = {}
    for i in range(num_replicates):
        peptide_array = np.arange(num_peptides)
        np.random.shuffle(peptide_array)
        pool_assignments = {}
        replicate_to_pool_to_peptides[i] = pool_assignments

        for j in range(num_pools):
            start_idx = peptides_per_pool * j
            end_idx = peptides_per_pool * (j + 1)
            pool_assignments[j] = peptide_array[start_idx:end_idx]
    return replicate_to_pool_to_peptides


def is_valid(
    replicate_to_pool_to_peptides: Solution, invalid_neighbors: ConstraintList = []
) -> bool:
    peptide_to_neighbors = defaultdict(set)

    if invalid_neighbors:
        # treat invalid pairs as if they've already been neighbors in a previous round
        for p1, p2 in invalid_neighbors:
            peptide_to_neighbors[p1].add(p2)
            peptide_to_neighbors[p2].add(p1)

    for replicate_idx, pool_to_peptides in replicate_to_pool_to_peptides.items():
        replicate_num = replicate_idx + 1

        # first check each peptide occurs once per replicate
        seen_peptides = set()
        for peptides in pool_to_peptides.values():
            for p in peptides:
                if p in seen_peptides:
                    print("Peptide %s twice in replicate %s" % (p, replicate_num))
                    return False
                seen_peptides.add(p)
        # next check to make sure that each peptides only paired with another at most once
        for peptides in pool_to_peptides.values():
            for p1 in peptides:
                for p2 in peptides:
                    if p1 != p2:
                        if p2 in peptide_to_neighbors[p1]:
                            print(
                                "Peptides %s and %s already together previous pool before replicate %s"
                                % (p1, p2, replicate_num)
                            )
                            return False
                        peptide_to_neighbors[p1].add(p2)
    return True


def count_violations(
    replicate_to_pool_to_peptides: Solution, invalid_neighbors: ConstraintList = []
) -> int:
    """
    Return the number of invalid peptide pairs in the solution
    """
    peptide_to_neighbors = defaultdict(set)

    if invalid_neighbors:
        # treat invalid pairs as if they've already been neighbors in a previous round
        for p1, p2 in invalid_neighbors:
            peptide_to_neighbors[p1].add(p2)
            peptide_to_neighbors[p2].add(p1)

    violations = 0
    for _, pool_to_peptides in replicate_to_pool_to_peptides.items():
        # first check each peptide occurs once per replicate
        seen_peptides = set()
        for peptides in pool_to_peptides.values():
            for p in peptides:
                if p in seen_peptides:
                    violations += 1
                else:
                    seen_peptides.add(p)
        # next check to make sure that each peptides only paired with another at most once
        for peptides in pool_to_peptides.values():
            for p1 in peptides:
                for p2 in peptides:
                    if p1 != p2:
                        if p2 in peptide_to_neighbors[p1]:
                            violations += 1
                        else:
                            peptide_to_neighbors[p1].add(p2)
    return violations


import random
from typing import Mapping

Replicate = int
Peptide = int
Pool = int
SwapCandidateList = Iterable[tuple[Replicate, Pool, Peptide]]
ReplicateToNeighborDict = Mapping[Replicate, Mapping[Peptide, set[Peptide]]]


def find_violating_peptides(
    replicate_to_pool_to_peptides: Solution, invalid_neighbors: ConstraintList = []
) -> tuple[SwapCandidateList, ReplicateToNeighborDict]:
    peptide_to_neighbors = defaultdict(set)

    if invalid_neighbors:
        # treat invalid pairs as if they've already been neighbors in a previous round
        for p1, p2 in invalid_neighbors:
            peptide_to_neighbors[p1].add(p2)
            peptide_to_neighbors[p2].add(p1)
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


def improve_solution(
    replicate_to_pool_to_peptides: Solution,
    invalid_neighbors: ConstraintList = [],
    verbose: bool = False,
):
    needs_swap, replicate_to_neighbor_dict = find_violating_peptides(
        replicate_to_pool_to_peptides, invalid_neighbors
    )

    random.shuffle(needs_swap)

    swapped_pools = set()
    swapped_peptides = set()
    for replicate_idx, pool_idx_a, peptide_a in needs_swap:
        if pool_idx_a in swapped_pools or peptide_a in swapped_peptides:
            continue

        pool_to_peptides = replicate_to_pool_to_peptides[replicate_idx]

        other_peptides = []
        peptide_to_pool_idx = {}
        for pool_idx_i, pool_peptides_i in pool_to_peptides.items():
            if pool_idx_i != pool_idx_a and pool_idx_i not in swapped_pools:
                all_peptides_ok = True
                for p in pool_peptides_i:
                    peptide_to_pool_idx[p] = pool_idx_i
                    all_peptides_ok = all_peptides_ok and (
                        peptide_a not in replicate_to_neighbor_dict[replicate_idx][p]
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
        if verbose:
            print("Before swap")
            print(
                "pool", pool_idx_a, "peptide", peptide_a, pool_to_peptides[pool_idx_a]
            )
            print(
                "pool", pool_idx_b, "peptide", peptide_b, pool_to_peptides[pool_idx_b]
            )

        # actually swap them
        pool_to_peptides[pool_idx_a] = np.array(
            [i for i in pool_to_peptides[pool_idx_a] if i != peptide_a] + [peptide_b]
        )
        pool_to_peptides[pool_idx_b] = np.array(
            [i for i in pool_to_peptides[pool_idx_b] if i != peptide_b] + [peptide_a]
        )

        if verbose:
            print("After")
            print(
                "pool", pool_idx_a, "peptide", peptide_a, pool_to_peptides[pool_idx_a]
            )
            print(
                "pool", pool_idx_b, "peptide", peptide_b, pool_to_peptides[pool_idx_b]
            )

        swapped_pools.add(pool_idx_a)
        swapped_pools.add(pool_idx_b)
        swapped_peptides.add(peptide_a)
        swapped_peptides.add(peptide_b)

    return replicate_to_pool_to_peptides


def optimize(
    replicate_to_pool_to_peptides: Solution,
    invalid_neighbors: ConstraintList = [],
    max_iters: int = 100,
    verbose: bool = False,
    return_history: bool = False,
) -> bool:
    """
    Iteratively update solution by randomly swapping a violating peptide with a random other peptide

    Parameters
    ----------
    replicate_to_pool_to_peptides:
        Initial solution which will be modified in-place

    invalid_neighbors:
        Constraint list of peptide pairs which cannot be placed together

    max_iters:
        Maximum number of swaps to consider performing

    verbose:
        Print number of violations for each iteration

    return_history:
        Return array of constraint validation counts per iteration


    Returns
    -------
    True if non-violating solution found, False if solution still has violations after
    max_iters
    """
    old_num_violations = count_violations(
        replicate_to_pool_to_peptides, invalid_neighbors
    )
    history = [old_num_violations]
    if verbose:
        print("Initial solution has %s violations" % (old_num_violations,))
    for i in range(max_iters):
        history.append(old_num_violations)
        improve_solution(replicate_to_pool_to_peptides, invalid_neighbors)
        new_num_violations = count_violations(
            replicate_to_pool_to_peptides, invalid_neighbors
        )
        history.append(new_num_violations)
        if verbose:
            print("%d) %d -> %d" % (i + 1, old_num_violations, new_num_violations))

        old_num_violations = new_num_violations

        if new_num_violations == 0:
            if verbose:
                print("Found valid solution after %d swaps" % (i + 1,))
            break

    result = old_num_violations == 0
    if return_history:
        return result, np.array(history)
    else:
        return result
