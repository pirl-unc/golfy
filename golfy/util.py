from collections import defaultdict
from typing import Iterable

from .types import Peptide


def pairs_to_dict(peptide_pairs: Iterable[Peptide]):
    peptide_to_set_dict = defaultdict(set)

    for p1, p2 in peptide_pairs:
        peptide_to_set_dict[p1].add(p2)
        peptide_to_set_dict[p2].add(p1)
    return peptide_to_set_dict
