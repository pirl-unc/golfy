from typing import Mapping
from dataclasses import dataclass

import numpy as np

from .types import Replicate, Pool, Peptide
from .solution import Solution

SpotCounts = Mapping[Replicate, Mapping[Pool, int]]


def simulate_elispot_counts(
    s: Solution,
    num_hits=5,
    max_activity_per_well=2000,
    min_background_peptide_activity=0,
    max_background_peptide_activity=20,
    min_hit_peptide_activity=50,
    max_hit_peptide_activity=500,
    verbose=False,
) -> tuple[SpotCounts, set[Peptide]]:
    num_peptides = s.num_peptides

    background_peptide_activity_range = (
        max_background_peptide_activity - min_background_peptide_activity
    )
    hit_peptide_activity_range = max_hit_peptide_activity - min_hit_peptide_activity

    all_peptides = np.arange(num_peptides)
    np.random.shuffle(all_peptides)
    hit_peptides = all_peptides[:num_hits]
    not_hit_peptides = all_peptides[num_hits:]
    if verbose:
        print("Hits: %s" % (hit_peptides,))

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
    hit_activity[not_hit_peptides] = 0
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
    return (spot_counts, set(hit_peptides))


@dataclass
class LinearSystem:
    A: np.ndarray
    b: np.ndarray
    pool_tuple_to_idx: Mapping[tuple[Replicate, Pool], int]
    idx_to_pool_tuple: Mapping[int, tuple[Replicate, Pool]]


def create_linear_system(
    s: Solution, spot_counts: SpotCounts, verbose=False
) -> LinearSystem:
    num_peptides = s.num_peptides
    num_pools = s.num_pools()

    A = np.zeros((num_pools, num_peptides + 1)).astype(float)
    b = np.zeros(num_pools).astype(float)

    pool_tuple_to_idx = {}
    idx_to_pool_tuple = {}
    i = 0
    for r, d in spot_counts.items():
        for pool, spots in d.items():
            b[i] = spots
            pool_tuple_to_idx[(r, pool)] = i
            idx_to_pool_tuple[i] = (r, pool)
            for p in s.assignments[r][pool]:
                A[i, p] = 1
            # add a ones column for a constant offset
            A[i, num_peptides] = 1
            i += 1
    if verbose:
        print("Ax = b")
        print("=======")
        print("A.shape: %s" % (A.shape,))
        print("b.shape: %s" % (b.shape,))
        print("A:\n%s" % (A,))
        print("A col sums: %s" % (A.sum(axis=0)))
        print("A row sums: %s" % (A.sum(axis=1)))
        print("b:\n%s" % (b,))
    return LinearSystem(
        A=A,
        b=b,
        pool_tuple_to_idx=pool_tuple_to_idx,
        idx_to_pool_tuple=idx_to_pool_tuple,
    )


@dataclass
class DeconvolutionResult:
    activity_per_peptide: np.ndarray
    prob_hit_per_peptide: np.ndarray
    high_confidence_hits: set[Peptide]


def solve_linear_system(
    linear_system: LinearSystem,
    min_peptide_activity: float = 1.0,
    leave_on_out=True,
    sparse_solution=True,
    verbose=False,
) -> DeconvolutionResult:
    A = linear_system.A
    b = linear_system.b
    num_pools, num_peptides_with_constant = A.shape
    num_peptides = num_peptides_with_constant - 1
    row_indices = list(range(num_pools))
    if leave_on_out:
        loo_indices = row_indices
    else:
        loo_indices = [None]

    avg_activity = np.zeros(num_peptides)
    frac_hit = np.zeros(num_peptides)

    for loo_idx in loo_indices:
        subset_indices = np.array([i for i in row_indices if i != loo_idx])
        A_subset = A[subset_indices, :]
        b_subset = b[subset_indices]
        if sparse_solution:
            # L1 minimizatin to get a small set of confident active peptides
            import sklearn.linear_model

            lasso = sklearn.linear_model.Lasso(fit_intercept=False)
            lasso.fit(A_subset, b_subset)
            x_with_offset = lasso.coef_
            x, c = x_with_offset[:-1], x_with_offset[-1]
            if verbose:
                print("x = %s" % (x,))
                print("c = %s" % (c,))
        else:
            # least squares non-sparse solution, works horribly, have fun
            x_with_offset, residuals, rank, singular_values = np.linalg.lstsq(
                A, b, rcond=None
            )
            x, c = x_with_offset[:-1], x_with_offset[-1]
            if verbose:
                print("x = %s" % (x[:-1],))
                print("c = %f" % (c,))
                print("residuals = %s" % (residuals,))
                print("rank = %s" % (rank,))
                print("singular_values = %s" % (singular_values,))

        avg_activity += x
        frac_hit += (x > min_peptide_activity).astype(float)

    avg_activity /= len(loo_indices)
    frac_hit /= len(loo_indices)
    high_confidence_hits = set(np.where(frac_hit > 0.5)[0])

    return DeconvolutionResult(
        activity_per_peptide=avg_activity,
        prob_hit_per_peptide=frac_hit,
        high_confidence_hits=high_confidence_hits,
    )
