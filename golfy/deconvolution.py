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

from typing import Mapping, Literal
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
    background: float
    cutoff: float


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
    cutoff = 0
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

        x, c = x_with_offset[:-1], x_with_offset[-1]

        if verbose:
            print("x = %s" % (x,))
            print("c = %s" % (c,))

        avg_activity += x
        curr_cutoff = choose_cutoff(
            A=A_subset,
            b=b_subset,
            estimated_constant_offsset=c,
            min_peptide_activity=min_peptide_activity,
            sparse_estimate=sparse_solution,
        )

        frac_hit += (x > curr_cutoff).astype(float)
        cutoff += curr_cutoff

    avg_activity /= len(loo_indices)
    frac_hit /= len(loo_indices)
    high_confidence_hits = set(np.where(frac_hit > 0.5)[0])
    cutoff /= len(loo_indices)

    return DeconvolutionResult(
        activity_per_peptide=avg_activity,
        prob_hit_per_peptide=frac_hit,
        high_confidence_hits=high_confidence_hits,
        background=c,
        cutoff=cutoff,
    )


def em_step(
    linear_system: LinearSystem,
    beta: np.ndarray,
    # stability_factor=1e-6,
) -> np.ndarray:
    """
    Implements equation 5 from Strom et al (2016)
        (X.T @ (y / (X @ beta))) * (beta / X.sum(0))
    where X is the design matrix, y is the pool activities, and beta is the current estimate of the peptide activities
    """
    X = linear_system.A
    y = linear_system.b
    return (X.T @ (y / (X @ beta))) * (beta / X.sum(0))


def em_init(linear_system: LinearSystem) -> np.ndarray:
    """
    Initialize beta for E-M algorithm from Strom et al (2016)
    to have each peptide be the mean of its pool activities
    """
    X = linear_system.A
    y = linear_system.b
    return X.T @ y / X.sum(0)


def choose_cutoff(
    A: np.ndarray,
    b: np.ndarray,
    estimated_constant_offsset: float,
    min_peptide_activity: float,
    sparse_estimate: bool,
) -> float:
    empty_wells = np.where(A.sum(1) == 0)[0]
    any_empty_wells = len(empty_wells) > 0
    if any_empty_wells:
        empty_well_background = b[empty_wells].mean()
    else:
        empty_well_background = 0

    estimated_constant_offsset = max(0, estimated_constant_offsset)
    print(A, b, estimated_constant_offsset, empty_well_background, min_peptide_activity)
    if sparse_estimate:
        return (
            np.sqrt(estimated_constant_offsset)
            + empty_well_background
            + min_peptide_activity
        )
    else:
        return (
            3 * estimated_constant_offsset
            + empty_well_background
            + min_peptide_activity
        )


def em_deconvolve(
    linear_system: LinearSystem,
    min_peptide_activity=1.0,
    max_iters=500,
    verbose=False,
    tol=1e-3,
    sparsity_threshold=1e-4,
) -> DeconvolutionResult:
    """
    Implements the expectation-maximization algorithm from
    "A statistical approach to determining responses to individual peptides from pooled-peptide ELISpot data"
    """
    beta = em_init(linear_system=linear_system)
    last_beta = beta.copy()
    for i in range(max_iters):
        beta = em_step(linear_system=linear_system, beta=beta)
        mean_abs_diff = np.mean(np.abs(beta - last_beta))
        if verbose:
            print(
                "[em_deconvolve] Iter %d: beta = %s (mean abs change=%f)"
                % (i, beta.round(2), mean_abs_diff)
            )
        if mean_abs_diff < tol:
            if verbose:
                print("[em_deconvolve] Stopping at iter %d" % i)
            break
        beta[np.abs(beta) < sparsity_threshold] = 0
        last_beta = beta.copy()
    activity_per_peptide = beta[:-1]
    activity_per_peptide[activity_per_peptide < sparsity_threshold] = 0
    print(activity_per_peptide.round(1))
    background = beta[-1]
    cutoff = choose_cutoff(
        A=linear_system.A,
        b=linear_system.b,
        estimated_constant_offsset=background,
        min_peptide_activity=min_peptide_activity,
        sparse_estimate=False,
    )
    above_limit = activity_per_peptide > cutoff
    return DeconvolutionResult(
        activity_per_peptide=activity_per_peptide,
        prob_hit_per_peptide=above_limit.astype(float),
        high_confidence_hits=set(np.where(above_limit)[0]),
        background=background,
        cutoff=cutoff,
    )


def deconvolve(
    s: Design,
    spot_counts: SpotCounts,
    method: Literal["lasso", "ridge", "em"] = "lasso",
    min_peptide_activity=1.0,
    verbose=False,
):
    """
    Arguments
    ---------
    s
        Dictionary containing mapping of replicates -> pools -> peptides

    spot_counts
        Dictionary containing mapping of replicates -> pool -> spot counts

    min_peptide_activity
        Minimum estimated activity of a peptide to be considered for hit set

    method
        - "lasso": L1 regularized linear regression
        - "ridge": L2 regularized linear regression
        - "em": EM algorithm from Strom et al. 2016
                      "A statistical approach to determining responses to individual peptides from pooled-peptide ELISpot data"

    verbose
        Print diagnostic information

    Returns DeconvolutionResult
    """
    linear_system = create_linear_system(s, spot_counts, verbose=verbose)
    if method == "lasso":
        return solve_linear_system(
            linear_system,
            sparse_solution=True,
            leave_on_out=False,
            min_peptide_activity=min_peptide_activity,
        )
    elif method == "ridge":
        return solve_linear_system(
            linear_system,
            sparse_solution=False,
            leave_on_out=False,
            min_peptide_activity=min_peptide_activity,
        )
    elif method == "em":
        return em_deconvolve(
            linear_system, min_peptide_activity=min_peptide_activity, verbose=verbose
        )
    else:
        raise ValueError("Unknown method: %s" % (method,))
