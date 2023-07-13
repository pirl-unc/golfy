from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Mapping

import numpy as np

from .types import PeptidePairList, Replicate, Pool, Peptide


@dataclass
class Spec:
    num_peptides: int
    max_peptides_per_pool: int
    num_replicates: int
    allow_extra_pools: bool
    invalid_neighbors: PeptidePairList
    preferred_neighbors: PeptidePairList

    def to_tuple(self):
        return (
            self.num_peptides,
            self.max_peptides_per_pool,
            self.num_replicates,
            self.allow_extra_pools,
            tuple(sorted(self.invalid_neighbors)),
            tuple(sorted(self.preferred_neighbors)),
        )

    def __eq__(self, other):
        return isinstance(other, Spec) and self.to_tuple() == other.to_tuple()

    def __hash__(self):
        return hash(self.to_tuple())


@dataclass
class Design(Spec):
    assignments: Mapping[Replicate, Mapping[Pool, Iterable[Peptide]]]

    def to_spec(self) -> Spec:
        return Spec(
            num_peptides=self.num_peptides,
            max_peptides_per_pool=self.max_peptides_per_pool,
            num_replicates=self.num_replicates,
            allow_extra_pools=self.allow_extra_pools,
            invalid_neighbors=self.invalid_neighbors,
            preferred_neighbors=self.preferred_neighbors,
        )

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

    def add_empty_pool(self, replicate_idx: int):
        """
        Add an empty pool to the given replicate
        """
        replicate = self.assignments[replicate_idx]
        num_pools = len(replicate)
        replicate[num_pools] = np.array([])

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

    def num_pools(self):
        total = 0
        for pool_to_peptides in self.assignments.values():
            total += len(pool_to_peptides)
        return total

    def average_peptides_per_pool(self):
        return self.num_peptides * self.num_replicates / self.num_pools()

    def peptide_to_pool_dict_for_replicate(
        self, replicate_idx: Replicate
    ) -> Mapping[Peptide, Pool]:
        return {
            p: pool
            for (pool, peptides) in self.assignments[replicate_idx].items()
            for p in peptides
        }

    def replicate_to_peptide_to_pool_dict(
        self,
    ) -> Mapping[Replicate, Mapping[Peptide, Pool]]:
        return {
            replicate_idx: self.peptide_to_pool_dict_for_replicate(replicate_idx)
            for replicate_idx in range(self.num_replicates)
        }

    def peptide_to_replicate_pool_pairs_dict(
        self,
    ) -> Mapping[Peptide, list[tuple[Replicate, Pool]]]:
        result = defaultdict(list)
        for (
            replicate_idx,
            peptide_to_pool,
        ) in self.replicate_to_peptide_to_pool_dict().items():
            for peptide, pool in peptide_to_pool.items():
                result[peptide].append((replicate_idx, pool))
        return result
