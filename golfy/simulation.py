from dataclasses import dataclass
from typing import Optional

import numpy as np

from .design import Design
from .types import Replicate, Pool, Peptide, SpotCounts


def _split_peptides_by_hits(
    num_peptides: int,
    num_hits: int,
    hit_peptides: Optional[set[Peptide]] = None,
    verbose: bool = False,
) -> tuple[set[Peptide], set[Peptide]]:
    if hit_peptides is None:
        all_peptides = np.arange(num_peptides)
        np.random.shuffle(all_peptides)
        hit_peptides = set(all_peptides[:num_hits])
        not_hit_peptides = set(all_peptides[num_hits:])
        if verbose:
            print("Hits: %s" % (hit_peptides,))
    else:
        hit_peptides = set(hit_peptides)
        not_hit_peptides = set(range(num_peptides)).difference(hit_peptides)
    return (hit_peptides, not_hit_peptides)


def simulate_number_hits_per_pool(
    s: Design,
    num_hits: int = 5,
    activity_level: int = 1,
    hit_peptides: Optional[set[Peptide]] = None,
    verbose: bool = False,
) -> tuple[SpotCounts, set[Peptide]]:
    """
    Returns spot counts corresponding to the number of hit peptides per pool and set of hit peptides.
    """
    n = s.num_peptides
    hit_peptides, _ = _split_peptides_by_hits(
        n, num_hits, hit_peptides=hit_peptides, verbose=verbose
    )
    activity = np.zeros(n)
    activity[np.array(list(hit_peptides))] = activity_level
    spot_counts: SpotCounts = {
        r: {pool: sum(activity[peptides]) for (pool, peptides) in d.items()}
        for (r, d) in s.assignments.items()
    }
    if verbose:
        print("Spot counts: %s" % (spot_counts,))
    return spot_counts, hit_peptides


def simulate_any_hits_per_pool(
    s: Design,
    num_hits: int = 5,
    hit_peptides: Optional[set[Peptide]] = None,
    verbose: bool = False,
) -> tuple[SpotCounts, set[Peptide]]:
    """
    Returns a spot counts dictionary with a 1 if the pool contained a hit peptides and a 0 otherwise.
    Also returns the
    """
    (spot_counts, hit_peptides) = simulate_number_hits_per_pool(
        s, num_hits, 1, hit_peptides=hit_peptides, verbose=verbose
    )
    spot_counts_binary = {
        r: {pool: 1 if spots > 0 else 0 for (pool, spots) in d.items()}
        for (r, d) in spot_counts.items()
    }
    return (spot_counts_binary, hit_peptides)


def simulate_elispot_counts(
    s: Design,
    num_hits=5,
    max_activity_per_well=2000,
    min_background_peptide_activity=0,
    max_background_peptide_activity=20,
    min_hit_peptide_activity=50,
    max_hit_peptide_activity=500,
    hit_peptides: Optional[set[Peptide]] = None,
    verbose=False,
) -> tuple[SpotCounts, set[Peptide]]:
    num_peptides = s.num_peptides

    background_peptide_activity_range = (
        max_background_peptide_activity - min_background_peptide_activity
    )
    hit_peptide_activity_range = max_hit_peptide_activity - min_hit_peptide_activity

    hit_peptides, not_hit_peptides = _split_peptides_by_hits(
        s.num_peptides, num_hits, hit_peptides=hit_peptides, verbose=verbose
    )

    background = (
        np.random.rand(num_peptides) * background_peptide_activity_range
        + min_background_peptide_activity
    )
    if verbose:
        print("Background activity: %s" % (background,))

    hit_activity = (
        np.random.rand(num_peptides) * hit_peptide_activity_range
        + min_hit_peptide_activity
    )
    hit_activity[np.array(list(not_hit_peptides))] = 0
    if verbose:
        print("Hit activity: %s" % (hit_activity,))

    spot_counts: SpotCounts = {
        r: {
            pool: min(
                max_activity_per_well,
                int(sum([background[i] + hit_activity[i] for i in peptides])),
            )
            for (pool, peptides) in d.items()
        }
        for (r, d) in s.assignments.items()
    }
    if verbose:
        print("Spot counts: %s" % (spot_counts,))
    return (spot_counts, hit_peptides)
