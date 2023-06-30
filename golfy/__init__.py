from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Optional, Mapping
import random

import numpy as np


Replicate = int
Peptide = int
Pool = int
SwapCandidateList = Iterable[tuple[Replicate, Pool, Peptide]]
ReplicateToNeighborDict = Mapping[Replicate, Mapping[Peptide, set[Peptide]]]
PeptidePairList = Iterable[tuple[Peptide, Peptide]]


@dataclass
class Spec:
    num_peptides: int
    max_peptides_per_pool: int
    num_replicates: int
    invalid_neighbors: PeptidePairList
    preferred_neighbors: PeptidePairList


@dataclass
class Solution(Spec):
    assignments: Mapping[Replicate, Mapping[Pool, Iterable[Peptide]]]

    def move_peptide(
        self,
        replicate_idx: Replicate,
        old_pool_idx: Pool,
        peptide: Peptide,
        new_pool_idx: Pool,
    ):
        """
        Move a peptide from its current pool to a new pool
        """
        pool_to_peptides = self.assignments[replicate_idx]
        old_pool = pool_to_peptides[old_pool_idx]
        new_pool = pool_to_peptides[new_pool_idx]
        assert len(new_pool) < self.max_peptides_per_pool
        pool_to_peptides[new_pool_idx] = np.array(list(new_pool) + [peptide])

        pool_to_peptides[old_pool_idx] = np.array([p for p in old_pool if p != peptide])

    def swap_peptides(
        self,
        replicate_idx: Replicate,
        pool_idx_a: Pool,
        peptide_a: Peptide,
        pool_idx_b: Pool,
        peptide_b: Peptide,
    ):
        """
        Move peptide A from its current pool to the pool of peptide B (and vice versa)
        """
        pool_to_peptides = self.assignments[replicate_idx]
        pool_a = pool_to_peptides[pool_idx_a]
        pool_b = pool_to_peptides[pool_idx_b]
        pool_to_peptides[pool_idx_a] = np.array(
            [i for i in pool_a if i != peptide_a] + [peptide_b]
        )
        pool_to_peptides[pool_idx_b] = np.array(
            [i for i in pool_b if i != peptide_b] + [peptide_a]
        )

    def remove_empty_pools(self):
        """
        Delete any empty pools and renumber the pools to be contiguous
        """
        for replicate_idx, pool_to_peptides in self.assignments.items():
            to_delete = [
                pool_idx
                for (pool_idx, pool) in pool_to_peptides.items()
                if len(pool) == 0
            ]
            if len(to_delete) > 0:
                for pool_idx in to_delete:
                    del pool_to_peptides[pool_idx]
                index_mapping = {
                    old_idx: new_idx
                    for (new_idx, old_idx) in enumerate(sorted(pool_to_peptides.keys()))
                }
                self.assignments[replicate_idx] = {
                    index_mapping[pool_idx]: pool
                    for (pool_idx, pool) in pool_to_peptides.items()
                }


def random_init(
    num_peptides: int = 100,
    peptides_per_pool: int = 5,
    num_replicates: int = 3,
    num_pools: Optional[int] = None,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
) -> Solution:
    if num_pools is None:
        num_pools = int(np.ceil(num_peptides / peptides_per_pool))

    preferred_neighbor_dict = defaultdict(set)
    for p1, p2 in preferred_neighbors:
        preferred_neighbor_dict[p1].add(p2)
        preferred_neighbor_dict[p2].add(p1)

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

    return Solution(
        num_peptides=num_peptides,
        max_peptides_per_pool=peptides_per_pool,
        num_replicates=num_replicates,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        assignments=replicate_to_pool_to_peptides,
    )


def is_valid(s: Solution) -> bool:
    peptide_to_neighbors = defaultdict(set)

    if s.invalid_neighbors:
        # treat invalid pairs as if they've already been neighbors in a previous round
        for p1, p2 in s.invalid_neighbors:
            peptide_to_neighbors[p1].add(p2)
            peptide_to_neighbors[p2].add(p1)

    for replicate_idx, pool_to_peptides in s.assignments.items():
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


def count_violations(s: Solution) -> int:
    """
    Return the number of invalid peptide pairs in the solution
    """
    invalid_neighbors = s.invalid_neighbors
    replicate_to_pool_to_peptides = s.assignments

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


def find_violating_peptides(
    s: Solution,
) -> tuple[SwapCandidateList, ReplicateToNeighborDict]:
    replicate_to_pool_to_peptides = s.assignments
    invalid_neighbors = s.invalid_neighbors

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


def improve_solution(s: Solution, verbose: bool = False):
    replicate_to_pool_to_peptides = s.assignments
    invalid_neighbors = s.invalid_neighbors

    needs_swap, replicate_to_neighbor_dict = find_violating_peptides(s)

    random.shuffle(needs_swap)

    swapped_pools = set()
    swapped_peptides = set()
    for replicate_idx, pool_idx_a, peptide_a in needs_swap:
        if pool_idx_a in swapped_pools or peptide_a in swapped_peptides:
            continue

        pool_to_peptides = replicate_to_pool_to_peptides[replicate_idx]
        pool_a = pool_to_peptides[pool_idx_a]

        # if a pool is empty, just move the offending peptide there
        empty_pools = {
            pool_idx
            for (pool_idx, pool_peptides) in pool_to_peptides.items()
            if len(pool_peptides) == 0
        }
        if empty_pools:
            pool_idx_b = random.choice(list(empty_pools))
            s.move_peptide(replicate_idx, pool_idx_a, peptide_a, pool_idx_b)
            swapped_pools.add(pool_idx_a)
            swapped_pools.add(pool_idx_b)
            swapped_peptides.add(peptide_a)
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
    max_iters: int = 100,
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
        Print number of violations for each iteration

    add_pool_if_stuck
        If no improvements have been made for 10 iters, add a pool
        to the last replicate

    return_history
        Return array of constraint validation counts per iteration


    Returns True if non-violating solution found, False if solution still has violations after
    max_iters
    """
    replicate_to_pool_to_peptides = s.assignments
    invalid_neighbors = s.invalid_neighbors

    old_num_violations = count_violations(s)
    history = [old_num_violations]
    if verbose:
        print("Initial solution has %s violations" % (old_num_violations,))
    num_iters_without_improvement = 0
    for i in range(max_iters):
        history.append(old_num_violations)
        improve_solution(s)
        new_num_violations = count_violations(s)

        history.append(new_num_violations)
        if verbose:
            print("%d) %d -> %d" % (i + 1, old_num_violations, new_num_violations))

        if old_num_violations <= new_num_violations:
            num_iters_without_improvement += 1
        else:
            num_iters_without_improvement = 0

        old_num_violations = new_num_violations

        if new_num_violations == 0:
            if verbose:
                print("Found valid solution after %d swaps" % (i + 1,))
            break

        if num_iters_without_improvement > 10 and add_pool_if_stuck:
            last_replicate_idx = s.num_replicates - 1
            last_replicate = s.assignments[last_replicate_idx]
            num_pools = len(last_replicate)
            last_replicate[num_pools] = np.array([])
            if verbose:
                print(
                    "Adding pool %d to replicate %d"
                    % (num_pools, last_replicate_idx + 1)
                )

    # just in case we ended up with any empty pools, remove them from the solution
    s.remove_empty_pools()

    result = old_num_violations == 0
    if return_history:
        return result, np.array(history)
    else:
        return result
