from typing import Iterable, Mapping

Replicate = int
Peptide = int
Pool = int
SwapCandidateList = Iterable[tuple[Replicate, Pool, Peptide]]
ReplicateToNeighborDict = Mapping[Replicate, Mapping[Peptide, set[Peptide]]]
PeptidePairList = Iterable[tuple[Peptide, Peptide]]
SpotCounts = Mapping[Replicate, Mapping[Pool, int]]
