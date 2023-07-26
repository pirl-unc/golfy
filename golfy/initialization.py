from collections import defaultdict
import math
from typing import Optional, Literal

import numpy as np

from .design import Design
from .types import PeptidePairList, Replicate, Peptide
from .util import pairs_to_dict, transitive_closure
from .validity import count_violations


def _pools_per_replicate(
    num_peptides: int,
    max_peptides_per_pool: int,
    num_replicates: int,
    allow_extra_pools: bool,
    verbose: bool,
) -> dict[Replicate, int]:
    too_few_pools = num_peptides < max_peptides_per_pool**2
    if too_few_pools and not allow_extra_pools and verbose:
        print(
            "Warning: cannot satisfy constraints with %d peptides and %d peptides per pool (max allowed %d)"
            % (num_peptides, max_peptides_per_pool, int(np.sqrt(num_peptides)))
        )
    num_pools_first_replicate = int(np.ceil(num_peptides / max_peptides_per_pool))
    result = {0: num_pools_first_replicate}
    for replicate_idx in range(1, num_replicates):
        if too_few_pools and allow_extra_pools:
            max_peptides_per_pool = num_pools_first_replicate
            result[replicate_idx] = int(np.ceil(num_peptides / max_peptides_per_pool))
        else:
            result[replicate_idx] = num_pools_first_replicate
    return result


def _repeat_init(
    num_peptides: int,
    max_peptides_per_pool: int,
    num_replicates: int,
    num_pools_per_replicate: dict[Replicate, int],
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    allow_extra_pools: bool = False,
    verbose: bool = False,
):
    """
    Assign peptides to random blocks and then repeat those blocks r times
    """
    replicate_to_pool_to_peptides = {}
    for i in range(num_replicates):
        peptide_array = np.arange(num_peptides)
        np.random.shuffle(peptide_array)
        pool_assignments = {}
        num_pools = num_pools_per_replicate[i]
        peptides_per_pool = int(np.ceil(num_peptides / num_pools))
        for j in range(num_pools):
            start_idx = peptides_per_pool * j
            end_idx = peptides_per_pool * (j + 1)
            pool_assignments[j] = peptide_array[start_idx:end_idx]
        replicate_to_pool_to_peptides[i] = pool_assignments
    return Design(
        num_peptides=num_peptides,
        max_peptides_per_pool=max_peptides_per_pool,
        num_replicates=num_replicates,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        assignments=replicate_to_pool_to_peptides,
        allow_extra_pools=allow_extra_pools,
    )


def _random_peptide_order(
    num_peptides: int, peptide_to_preferred: dict[Peptide, set[Peptide]]
):
    """
    Return the peptide indices in a random order but make sure the ones with preferred neighbors are first
    """
    random_peptide_order = np.arange(num_peptides)
    np.random.shuffle(random_peptide_order)
    assert len(set(random_peptide_order)) == num_peptides
    # assign all peptides with preferred neighbors first
    peptides_with_preferred_neighbors = [
        p for p in random_peptide_order if peptide_to_preferred.get(p)
    ]
    peptides_without_preferred_neighbors = [
        p for p in random_peptide_order if not peptide_to_preferred.get(p)
    ]
    peptide_list = (
        peptides_with_preferred_neighbors + peptides_without_preferred_neighbors
    )
    return peptide_list


def _random_init(
    num_peptides: int,
    max_peptides_per_pool: int,
    num_replicates: int,
    num_pools_per_replicate: dict[Replicate, int],
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    allow_extra_pools: bool = False,
    verbose: bool = False,
) -> Design:
    replicate_to_pool_to_peptides = {}
    for i in range(num_replicates):
        peptide_array = np.arange(num_peptides)
        np.random.shuffle(peptide_array)
        pool_assignments = {}
        replicate_to_pool_to_peptides[i] = pool_assignments
        num_pools = num_pools_per_replicate[i]
        peptides_per_pool = int(np.ceil(num_peptides / num_pools))
        for j in range(num_pools):
            start_idx = peptides_per_pool * j
            end_idx = peptides_per_pool * (j + 1)
            pool_assignments[j] = peptide_array[start_idx:end_idx]

    return Design(
        num_peptides=num_peptides,
        max_peptides_per_pool=max_peptides_per_pool,
        num_replicates=num_replicates,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        assignments=replicate_to_pool_to_peptides,
        allow_extra_pools=allow_extra_pools,
    )


def _singleton_init(
    num_peptides: int,
    max_peptides_per_pool: int,
    num_replicates: int,
    num_pools_per_replicate: dict[Replicate, int],
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    allow_extra_pools: bool = False,
    verbose: bool = False,
) -> Design:
    """
    Initialize every peptide to be in its own pool
    """

    replicate_to_pool_to_peptides = {}
    for i in range(num_replicates):
        pool_to_peptides = {}
        for p in range(num_peptides):
            pool_to_peptides[p] = np.array([p])
        replicate_to_pool_to_peptides[i] = pool_to_peptides
    return Design(
        num_peptides=num_peptides,
        max_peptides_per_pool=max_peptides_per_pool,
        num_replicates=num_replicates,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        assignments=replicate_to_pool_to_peptides,
        allow_extra_pools=allow_extra_pools,
    )


def _valid_init(
    num_peptides: int,
    max_peptides_per_pool: int,
    num_replicates: int,
    num_pools_per_replicate: dict[Replicate, int],
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    allow_extra_pools: bool = False,
    verbose: bool = False,
) -> Design:
    peptide_to_invalid = pairs_to_dict(invalid_neighbors)
    peptide_to_preferred = transitive_closure(pairs_to_dict(preferred_neighbors))
    if verbose:
        print("[_valid_init] Invalid neighbors: %s" % (peptide_to_invalid,))
        print("[_valid_init] Preferred neighbors: %s" % (peptide_to_preferred,))

    replicate_to_pool_to_peptides = {}

    for i in range(num_replicates):
        num_pools = num_pools_per_replicate[i]
        peptides_per_pool = int(math.ceil(num_peptides / num_pools))

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

        for peptide in _random_peptide_order(
            num_peptides=num_peptides, peptide_to_preferred=peptide_to_preferred
        ):
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

            # otherwise, make a new pool
            if peptide not in peptide_to_pool:
                if not allow_extra_pools and verbose:
                    print(
                        "Warning: unable to create valid configuration without exceeding max pools for replicate %d (max %d, adding %d)"
                        % (i + 1, num_pools, curr_num_pools() + 1)
                    )
                make_new_pool(peptide)

        replicate_to_pool_to_peptides[i] = {
            pool_idx: np.array(sorted(peptides))
            for (pool_idx, peptides) in pool_to_peptides.items()
        }
    return Design(
        num_peptides=num_peptides,
        max_peptides_per_pool=max_peptides_per_pool,
        num_replicates=num_replicates,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        assignments=replicate_to_pool_to_peptides,
        allow_extra_pools=allow_extra_pools,
    )


def _greedy_init(
    num_peptides: int,
    max_peptides_per_pool: int,
    num_replicates: int,
    num_pools_per_replicate: dict[Replicate, int],
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    allow_extra_pools: bool = False,
    verbose: bool = False,
) -> Design:
    peptide_to_invalid = pairs_to_dict(invalid_neighbors)
    peptide_to_preferred = transitive_closure(pairs_to_dict(preferred_neighbors))
    if verbose:
        print("[_greedy_init] Invalid neighbors: %s" % (peptide_to_invalid,))
        print("[_greedy_init] Preferred neighbors: %s" % (peptide_to_preferred,))

    replicate_to_pool_to_peptides = {}

    for i in range(num_replicates):
        num_pools = num_pools_per_replicate[i]
        peptides_per_pool = int(math.ceil(num_peptides / num_pools))

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

        for peptide in _random_peptide_order(
            num_peptides=num_peptides, peptide_to_preferred=peptide_to_preferred
        ):
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
                if allow_extra_pools:
                    make_new_pool(peptide)
                else:
                    raise ValueError("Unable to create valid configuration")

        replicate_to_pool_to_peptides[i] = {
            pool_idx: np.array(sorted(peptides))
            for (pool_idx, peptides) in pool_to_peptides.items()
        }

    return Design(
        num_peptides=num_peptides,
        max_peptides_per_pool=max_peptides_per_pool,
        num_replicates=num_replicates,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        assignments=replicate_to_pool_to_peptides,
        allow_extra_pools=allow_extra_pools,
    )


def init(
    num_peptides: int = 100,
    max_peptides_per_pool: int = 5,
    num_replicates: int = 3,
    num_pools_per_replicate: Optional[int | dict[Replicate, int]] = None,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    strategy: Literal["greedy", "random", "valid", "singleton", "repeat"] = "greedy",
    allow_extra_pools: bool = False,
    verbose: bool = False,
) -> Design:
    """
    Initialize a Solution for a given configuration

    Args
    ----
    num_peptides
        number of peptides (default: 100)

    peptides_per_pool
        number of peptides per pool (default: 5)

    num_replicates
        number of replicates in the Solution (default: 3)

    num_pools_per_replicate
        number of pools in each replicate (default: ceil(num_peptides / peptides_per_pool))

    invalid_neighbors
        list of peptide pairs that cannot be in the same pool

    preferred_neighbors
        list of peptide pairs that should be in the same pool

    strategy
        initialization strategy, one of {"greedy", "random", "singleton", "valid", "repeat"} (default: "greedy")

    allow_extra_pools
        whether to allow extra pools to be created to satisfy constraints (default: False)

    verbose
        print diagnostic information during initialization (default: False)
    """
    fn = {
        "greedy": _greedy_init,
        "random": _random_init,
        "singleton": _singleton_init,
        "valid": _valid_init,
        "repeat": _repeat_init,
    }.get(strategy)
    if num_pools_per_replicate is None:
        num_pools_per_replicate = _pools_per_replicate(
            num_peptides=num_peptides,
            max_peptides_per_pool=max_peptides_per_pool,
            num_replicates=num_replicates,
            allow_extra_pools=allow_extra_pools,
            verbose=verbose,
        )
    elif isinstance(num_pools_per_replicate, int):
        num_pools_per_replicate = {
            i: num_pools_per_replicate for i in range(num_replicates)
        }
    if fn is None:
        raise ValueError("Unknown initialization strategy: '%s'" % (strategy,))

    if verbose:
        print(
            "Using '%s' initialization strategy, pools per replicate: %s"
            % (strategy, num_pools_per_replicate)
        )
    kwargs = dict(
        num_peptides=num_peptides,
        max_peptides_per_pool=max_peptides_per_pool,
        num_replicates=num_replicates,
        num_pools_per_replicate=num_pools_per_replicate,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        allow_extra_pools=allow_extra_pools,
        verbose=verbose,
    )
    s = fn(**kwargs)
    violations = count_violations(s)
    if verbose:
        print("Generated solution with %d initial violations" % (violations))
    return s
