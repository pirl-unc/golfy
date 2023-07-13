from typing import Mapping
from dataclasses import dataclass

import numpy as np

from .types import Peptide, SpotCounts, Pool, Replicate
from .design import Design


@dataclass
class LinearSystem:
    A: np.ndarray
    b: np.ndarray
    pool_tuple_to_idx: Mapping[tuple[Replicate, Pool], int]
    idx_to_pool_tuple: Mapping[int, tuple[Replicate, Pool]]


def create_linear_system(
    s: Design, spot_counts: SpotCounts, verbose=False
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
    alpha=None,
    verbose=False,
) -> DeconvolutionResult:
    from sklearn.linear_model import Lasso, Ridge, LassoCV

    A = linear_system.A
    b = linear_system.b

    num_pools, num_peptides_with_constant = A.shape
    if alpha is None:
        alpha = 1.0 / (2 * num_peptides_with_constant)
    if verbose:
        print(
            "[solve_linear_system] A.shape: %s, b.shape: %s, alpha = %0.4f"
            % (A.shape, b.shape, alpha)
        )

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
        b_min = np.percentile(b_subset, 1)
        b_max = np.percentile(b_subset, 99)
        scale = b_max - b_min

        assert scale > 0

        if sparse_solution:
            # L1 minimization to get a small set of confident active peptides
            model = Lasso(
                fit_intercept=False,
                positive=True,
                alpha=alpha,
                selection="random",
            )

        else:
            # this will work horribly, have fun
            model = Ridge(fit_intercept=False, positive=True, alpha=alpha)

        model.fit(A_subset, np.maximum(0, b_subset - b_min) / scale)
        x_with_offset = scale * model.coef_

        if verbose:
            print("x = %s" % (x,))
            print("c = %s" % (c,))
        x, c = x_with_offset[:-1], x_with_offset[-1]
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
