import numpy as np

from .design import Design
from .types import Replicate, Pool, Peptide
from .util import pairs_to_dict, transitive_closure


def _make_peptide_to_invalid_excluding_replicate(
    s: Design,
) -> dict[Replicate, dict[Peptide, set[Peptide]]]:
    peptide_to_invalid = pairs_to_dict(s.invalid_neighbors)
    peptide_to_invalid_excluding_replicate = {
        replicate_idx: {
            peptide: peptide_to_invalid.get(peptide, set()).copy()
            for peptide in range(s.num_peptides)
        }
        for replicate_idx in range(s.num_replicates)
    }
    replicate_to_peptide_to_pool_idx = s.replicate_to_peptide_to_pool_dict()

    for replicate_idx, pool_to_peptides_1 in s.assignments.items():
        for other_replicate_idx, pool_to_peptides_2 in s.assignments.items():
            if replicate_idx == other_replicate_idx:
                continue
            for peptides in pool_to_peptides_1.values():
                for peptide in peptides:
                    other_pool_idx = replicate_to_peptide_to_pool_idx[
                        other_replicate_idx
                    ][peptide]

                    for other_peptide in pool_to_peptides_2[other_pool_idx]:
                        if other_peptide != peptide:
                            peptide_to_invalid_excluding_replicate[replicate_idx][
                                peptide
                            ].add(other_peptide)
                            peptide_to_invalid_excluding_replicate[replicate_idx][
                                other_peptide
                            ].add(peptide)
    return peptide_to_invalid_excluding_replicate


def merge_small_pools(s: Design, verbose: bool = False) -> int:
    num_merged = 0

    peptide_to_preferred = transitive_closure(pairs_to_dict(s.preferred_neighbors))
    peptide_to_invalid_excluding_replicate = (
        _make_peptide_to_invalid_excluding_replicate(s)
    )

    replicate_order = list(range(s.num_replicates))

    for i in replicate_order:
        pool_to_peptides = s.assignments[i]

        merged = set()
        # next, try to merge any pools that are smaller than the max size
        for pool_idx_1, peptides_1 in list(pool_to_peptides.items()):
            if pool_idx_1 in merged:
                continue

            if len(peptides_1) >= s.max_peptides_per_pool:
                continue

            n_peptides_1 = len(peptides_1)
            if n_peptides_1 == 0:
                # skip empty pools, they'll be removed at the end
                continue

            candidates = list(pool_to_peptides.items())
            # split the candidates into those with preferred peptide overlaps and those without
            # while also filtering out pools that are invalid or have already been merged
            candidates_with_preferred_neighbors = []
            candidates_without_preferred_neighbors = []
            for pool_idx_2, peptides_2 in candidates:
                if pool_idx_1 == pool_idx_2:
                    # can't merge a pool with itself
                    continue
                if len(peptides_2) == 0:
                    # skip empty pools, they'll be removed at the end
                    continue
                if n_peptides_1 + len(peptides_2) >= s.max_peptides_per_pool:
                    # can't exceed the max pool size
                    continue
                if pool_idx_2 in merged:
                    continue

                any_preferred = False
                peptides_2_set = set(peptides_2)

                for p1 in peptides_1:
                    preferred = peptide_to_preferred.get(p1)
                    if not preferred:
                        continue
                    intersection = peptides_2_set.intersection(preferred)
                    if len(intersection) > 0:
                        any_preferred = True
                        break
                if any_preferred:
                    candidates_with_preferred_neighbors.append((pool_idx_2, peptides_2))
                else:
                    candidates_without_preferred_neighbors.append(
                        (pool_idx_2, peptides_2)
                    )

            for pool_idx_2, peptides_2 in (
                candidates_with_preferred_neighbors
                + candidates_without_preferred_neighbors
            ):
                all_valid = True
                # check if any peptides in pool 1 are invalid with any peptides in pool 2
                combined_peptides = list(peptides_1) + list(peptides_2)
                for peptide in list(combined_peptides):
                    other_peptides = {p for p in combined_peptides if p != peptide}
                    all_valid = (
                        len(
                            peptide_to_invalid_excluding_replicate[i][
                                peptide
                            ].intersection(other_peptides)
                        )
                        == 0
                    )
                    if not all_valid:
                        break

                if all_valid:
                    num_merged += 1
                    if verbose:
                        print(
                            "-- merging pools %d and %d in replicate %d"
                            % (pool_idx_1, pool_idx_2, i + 1)
                        )
                    pool_to_peptides[pool_idx_1] = np.array(combined_peptides)
                    pool_to_peptides[pool_idx_2] = np.array([])
                    merged.update([pool_idx_1, pool_idx_2])
                    for p in combined_peptides:
                        other_peptides = {p for p in combined_peptides if p != peptide}
                        for other_replicate_idx in range(s.num_replicates):
                            if other_replicate_idx != i:
                                peptide_to_invalid_excluding_replicate[
                                    other_replicate_idx
                                ][p].update(other_peptides)
                    break

    # just in case we ended up with any empty pools, remove them from the solution
    s.remove_empty_pools()
    return num_merged


def cleanup(s: Design, verbose: bool = True, max_iters=1000) -> int:
    total_num_merged = 0
    prev_num_pools = s.num_pools()

    for iter_idx in range(max_iters):
        num_merged = merge_small_pools(s)
        total_num_merged += num_merged
        num_pools = s.num_pools()
        if verbose:
            print(
                "Merge iter %d: merged %d small pools, final number of pools: %d (prev %d) "
                % (iter_idx + 1, num_merged, s.num_pools(), prev_num_pools)
            )
        if num_pools >= prev_num_pools:
            break
        prev_num_pools = num_pools
    return total_num_merged
