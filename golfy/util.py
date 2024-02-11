# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations


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
