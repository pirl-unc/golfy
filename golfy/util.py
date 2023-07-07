from collections import defaultdict
from typing import Iterable, Mapping

from .types import Peptide


def pairs_to_dict(peptide_pairs: Iterable[Peptide]):
    peptide_to_set_dict = defaultdict(set)

    for p1, p2 in peptide_pairs:
        peptide_to_set_dict[p1].add(p2)
        peptide_to_set_dict[p2].add(p1)
    return peptide_to_set_dict


def peptide_to_transitive_neighbors(
    peptide: Peptide, peptide_to_neighbors: Mapping[Peptide, set[Peptide]]
) -> set[Peptide]:
    """
    Find all peptides that are neighbors of a given peptide, either directly or
    indirectly (via a chain of neighbors)
    """
    neighbors = set()
    to_visit = list(peptide_to_neighbors[peptide])
    while to_visit:
        p = to_visit.pop()
        if p not in neighbors and p != peptide:
            neighbors.add(p)
            to_visit.extend(peptide_to_neighbors[p])
    return neighbors


def transitive_closure(
    peptide_to_neighbors: Mapping[Peptide, set[Peptide]]
) -> Mapping[Peptide, set[Peptide]]:
    return {
        peptide: peptide_to_transitive_neighbors(peptide, peptide_to_neighbors)
        for peptide in peptide_to_neighbors
    }
