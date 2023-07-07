from dataclasses import dataclass
import numpy as np
from typing import Iterable, Mapping

from .types import PeptidePairList, Replicate, Pool, Peptide


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

    def peptide_to_replicate_pool_pair_dict(
        self,
    ) -> Mapping[Peptide, tuple[Replicate, Pool]]:
        return {
            p: (replicate_idx, pool)
            for (
                replicate_idx,
                peptide_to_pool,
            ) in self.replicate_to_peptide_to_pool_dict().items()
            for (p, pool) in peptide_to_pool.items()
        }
