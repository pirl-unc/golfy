from collections import defaultdict
from typing import Optional, Literal

import numpy as np

from .solution import Solution
from .types import PeptidePairList
from .util import pairs_to_dict, transitive_closure
from .validity import count_violations


def _random_init(
    num_peptides: int = 100,
    peptides_per_pool: int = 5,
    num_replicates: int = 3,
    num_pools: Optional[int] = None,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    verbose: bool = False,
) -> Solution:
    if num_pools is None:
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

    return Solution(
        num_peptides=num_peptides,
        max_peptides_per_pool=peptides_per_pool,
        num_replicates=num_replicates,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        assignments=replicate_to_pool_to_peptides,
    )


def _greedy_init(
    num_peptides: int = 100,
    peptides_per_pool: int = 5,
    num_replicates: int = 3,
    num_pools: Optional[int] = None,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    verbose: bool = False,
) -> Solution:
    if num_pools is None:
        num_pools = int(np.ceil(num_peptides / peptides_per_pool))

    peptide_to_invalid = pairs_to_dict(invalid_neighbors)
    peptide_to_preferred = transitive_closure(pairs_to_dict(preferred_neighbors))
    if verbose:
        print("[_greedy_init] Invalid neighbors: %s" % (peptide_to_invalid,))
        print("[_greedy_init] Preferred neighbors: %s" % (peptide_to_preferred,))

    replicate_to_pool_to_peptides = {}

    random_peptide_order = np.arange(num_peptides)
    np.random.shuffle(random_peptide_order)
    assert len(set(random_peptide_order)) == num_peptides
    # assign all peptides with preferred neighbors first
    peptides_with_preferred_neighbors = [
        p for p in random_peptide_order if p in peptide_to_preferred
    ]
    peptides_without_preferred_neighbors = [
        p for p in random_peptide_order if p not in peptide_to_preferred
    ]
    peptide_list = (
        peptides_with_preferred_neighbors + peptides_without_preferred_neighbors
    )
    for i in range(num_replicates):
        peptide_to_pool = {}
        pool_to_peptides = defaultdict(set)

        def add_to_pool(peptide, pool_idx):
            pool = pool_to_peptides[pool_idx]
            pool.add(peptide)
            peptide_to_pool[peptide] = pool_idx
            for other_peptide in pool:
                peptide_to_invalid[peptide].add(other_peptide)
                peptide_to_invalid[other_peptide].add(peptide)

        def curr_num_pools():
            return len(pool_to_peptides)

        def make_new_pool(peptide):
            new_pool_idx = curr_num_pools()
            add_to_pool(peptide, new_pool_idx)
            return new_pool_idx

        for peptide in peptide_list:
            for preferred_neighbor in peptide_to_preferred.get(peptide, []):
                if preferred_neighbor in peptide_to_invalid[peptide]:
                    if verbose:
                        print(
                            "[_greedy_init] replicate %d, peptide %d already invalid with preferred neighbor %d"
                            % (i, peptide, preferred_neighbor)
                        )
                    continue
                preferred_pool_idx = peptide_to_pool.get(preferred_neighbor)
                if preferred_pool_idx is None:
                    if verbose:
                        print(
                            "[_greedy_init] replicate %d, peptide %d, preferred neighbor %d not in a pool yet"
                            % (i, peptide, preferred_neighbor)
                        )
                    continue
                preferred_pool = pool_to_peptides[preferred_pool_idx]
                if len(preferred_pool) >= peptides_per_pool:
                    if verbose:
                        print(
                            "[_greedy_init] replicate %d, preferred neighbor %d in pool %d which is already full"
                            % (i, preferred_neighbor, preferred_pool_idx)
                        )
                    continue

                add_to_pool(peptide, preferred_pool_idx)
                if verbose:
                    print(
                        "[_greedy_init] replicate %d, adding peptide %d to preferred pool %d to pair with peptide %d"
                        % (i, peptide, preferred_pool_idx, preferred_neighbor)
                    )
                break
            # if we didn't get a preferred peptide pool that's valid
            # and there's room for more pools, just make a singleton
            if peptide not in peptide_to_pool:
                if curr_num_pools() < num_pools:
                    assigned_pool_idx = make_new_pool(peptide)
                    if verbose:
                        print(
                            "[_greedy_init] replicate %d, making new pool %d for peptide %d"
                            % (i, assigned_pool_idx, peptide)
                        )

            # otherwise, try to find a valid pool
            if peptide not in peptide_to_pool:
                for pool_idx, pool in pool_to_peptides.items():
                    if len(pool) < peptides_per_pool:
                        disallowed_peptides = peptide_to_invalid[peptide]
                        valid = all(
                            [
                                other_peptide not in disallowed_peptides
                                for other_peptide in pool
                            ]
                        )
                        if valid:
                            add_to_pool(peptide, pool_idx)
                            break

            # otherwise, pick any pool less than the max size
            if peptide not in peptide_to_pool:
                for pool_idx, pool in pool_to_peptides.items():
                    if len(pool) < peptides_per_pool:
                        add_to_pool(peptide, pool_idx)
                        break

            # lastly, violate the num_pools constraint to make a
            # new singleton anyways
            if peptide not in peptide_to_pool:
                make_new_pool(peptide)

        replicate_to_pool_to_peptides[i] = {
            pool_idx: np.array(sorted(peptides))
            for (pool_idx, peptides) in pool_to_peptides.items()
        }

    return Solution(
        num_peptides=num_peptides,
        max_peptides_per_pool=peptides_per_pool,
        num_replicates=num_replicates,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        assignments=replicate_to_pool_to_peptides,
    )


def init(
    num_peptides: int = 100,
    peptides_per_pool: int = 5,
    num_replicates: int = 3,
    num_pools: Optional[int] = None,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    strategy: Literal["greedy", "random"] = "greedy",
    verbose=False,
) -> Solution:
    if strategy == "greedy":
        s = _greedy_init(
            num_peptides=num_peptides,
            peptides_per_pool=peptides_per_pool,
            num_replicates=num_replicates,
            num_pools=num_pools,
            invalid_neighbors=invalid_neighbors,
            preferred_neighbors=preferred_neighbors,
            verbose=verbose,
        )
    elif strategy == "random":
        s = _random_init(
            num_peptides=num_peptides,
            peptides_per_pool=peptides_per_pool,
            num_replicates=num_replicates,
            num_pools=num_pools,
            invalid_neighbors=invalid_neighbors,
            preferred_neighbors=preferred_neighbors,
            verbose=verbose,
        )
    else:
        raise ValueError("Unknown initialization strategy: '%s'" % strategy)

    violations = count_violations(s)
    if verbose:
        print("Generated solution with %d initial violations" % (violations))
    return s
