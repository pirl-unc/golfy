import math
from typing import Optional, Literal

import numpy as np
import tqdm

from .design import Design
from .initialization import init
from .optimization import optimize
from .simulation import simulate_number_hits_per_pool
from .deconvolution import create_linear_system, solve_linear_system
from .types import Replicate, PeptidePairList
from .validity import count_violations


def find_best_design(
    num_peptides: int = 100,
    max_peptides_per_pool: int = 5,
    num_replicates: int = 3,
    num_pools_per_replicate: Optional[int | dict[Replicate, int]] = None,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    allow_extra_pools: bool = False,
    verbose: bool = False,
) -> Design:
    """
    Try several different initialization methods and return the solution which
    is optimized to have the fewest violations and the fewest pools.

    Args
    ----
    num_peptides
        The number total number of peptides in the experiment

    max_peptides_per_pool
        Maximum number of peptides which can be in any one pool. This is also the number of peptides
        which will be in most pools if allow_extra_pools is False

    num_replicates
        Number of replicates in the experiment (i.e. how many distinct pools each peptide occurs in)

    num_pools_per_replicate
        Number of pools in each replicate. If None, then this will be set to the minimum number of pools
        required to fit all peptides. If an int, then this will be the number of pools in each replicate.

    invalid_neighbors
        List of peptide pairs which cannot be in the same pool

    preferred_neighbors
        List of peptide pairs which should be in the same pool if possible

    allow_extra_pools
        If True, then the solution can have more than the minimum number of pools required to fit all peptides. If False,
        then the returned solution may not be valid (i.e. some peptide pairs will occur in more than one pool together)

    verbose
        If True, then print out information about the solution as it is being constructed
    """
    shared_kwargs = dict(
        num_peptides=num_peptides,
        max_peptides_per_pool=max_peptides_per_pool,
        num_replicates=num_replicates,
        num_pools_per_replicate=num_pools_per_replicate,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        allow_extra_pools=allow_extra_pools,
        verbose=verbose,
    )

    if len(preferred_neighbors) > 0:
        if allow_extra_pools:
            # these are the two methods which can group
            # preferred peptides in pools, either during
            # initialization or during merging
            init_strategies = ["greedy", "singleton"]
        else:
            # if we don't allow extra pools, then we can only
            # group preferred peptides during initialization
            init_strategies = ["greedy"]
    elif allow_extra_pools:
        # if there are no preferred peptides, then we can
        # only use all the initialization methods
        init_strategies = ["greedy", "random", "valid", "singleton"]
    else:
        # only greedy and random initialization methods
        # let us be strict about the number of pools
        init_strategies = ["greedy", "random"]

    designs = {
        strategy: init(strategy=strategy, **shared_kwargs)
        for strategy in init_strategies
    }
    best_design = None
    best_violations = None
    best_num_pools = None

    for strategy, s in designs.items():
        if verbose:
            print(
                "Initialized with strategy '%s': violations=%d, num_pools=%d"
                % (strategy, count_violations(s), s.num_pools())
            )
        optimize(s, allow_extra_pools=allow_extra_pools, verbose=verbose)
        violations = count_violations(s)
        num_pools = s.num_pools()
        if verbose:
            print(
                "-- after optimization of '%s' solution: violations=%d, num_pools=%d"
                % (strategy, violations, num_pools)
            )
        if (
            best_design is None
            or (violations < best_violations)
            or (violations == best_violations and num_pools < best_num_pools)
        ):
            best_design = s
            best_violations = violations
            best_num_pools = num_pools
            if verbose:
                print("^^ new best solution")
    return best_design


def evaluate_design(
    s: Design,
    num_simulation_iters: int = 10,
    min_hit_fraction: float = 0.01,
    max_hit_fraction: float = 0.05,
) -> tuple[float, float, float]:
    """
    Returns average precision, recall, F1 scores across multiple simulations of the given design for
    num_hits in range of 1% to 5% of the total number of peptides
    """
    ps = []
    rs = []
    f1s = []
    for _ in range(num_simulation_iters):
        for hit_fraction in np.arange(min_hit_fraction, max_hit_fraction + 0.001, 0.01):
            num_hits = int(np.ceil(hit_fraction * s.num_peptides))

            spot_counts, hit_peptides = simulate_number_hits_per_pool(
                s, num_hits=num_hits
            )
            linear_system = create_linear_system(s, spot_counts)
            result = solve_linear_system(
                linear_system, leave_on_out=False, min_peptide_activity=0.2
            )

            predicted_hits = result.high_confidence_hits

            tp = len(predicted_hits.intersection(hit_peptides))
            fp = len({p for p in predicted_hits if p not in hit_peptides})
            fn = len({p for p in hit_peptides if p not in predicted_hits})

            precision = tp / (tp + fp) if tp + fp > 0 else 0
            recall = tp / (tp + fn) if tp + fn > 0 else 0
            f1 = (
                2 * precision * recall / (precision + recall)
                if precision + recall > 0
                else 0
            )

            ps.append(precision)
            rs.append(recall)
            f1s.append(f1)
    return (np.mean(ps), np.mean(rs), np.mean(f1s))


def best_design_for_pool_budget(
    num_peptides: int = 100,
    max_pools: int = 96,
    num_simulation_iters: int = 2,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    verbose: bool = False,
):
    assert num_peptides > 1, "No need to pool if there's only one peptide"
    assert max_pools > 1, "Must have more than one pool"
    assert max_pools <= num_peptides, "Can't have more pools than peptides"

    shared_kwargs = dict(
        num_peptides=num_peptides,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        verbose=verbose,
    )

    designs = {}
    for max_peptides_per_pool in range(2, max(3, num_peptides // 5)):
        if max_peptides_per_pool * max_pools < num_peptides:
            # not enough pools to fit all peptides without replicates
            continue

        for num_replicates in range(2, 6):
            min_pools_for_spec = math.ceil(
                num_peptides * num_replicates / max_peptides_per_pool
            )
            if min_pools_for_spec > max_pools:
                # not enough pools to fit all peptides with given number of replicates
                continue

            for allow_extra_pools in [True, False]:
                if min_pools_for_spec == max_pools and allow_extra_pools:
                    # no need to allow extra pools if we're already at the maximum
                    continue

                s = find_best_design(
                    max_peptides_per_pool=max_peptides_per_pool,
                    num_replicates=num_replicates,
                    allow_extra_pools=allow_extra_pools,
                    **shared_kwargs,
                )
                if s is None:
                    continue
                num_pools = s.num_pools()
                num_violations = count_violations(s)
                key = s.to_spec()
                if num_pools <= max_pools:
                    p, r, f1 = evaluate_design(s, num_simulation_iters)

                    print(
                        "%s: %d pools, %d violations, precision=%0.2f, recall=%0.2f, f1=%0.2f"
                        % (key, num_pools, num_violations, p, r, f1)
                    )
                    # maximize f1, minimize violations,  maximize (precision, recall), minimize pools, maximize replicates
                    sort_key = (
                        -round(f1, 2),
                        num_violations,
                        -round(p, 2),
                        -round(r, 2),
                        num_pools,
                        -num_replicates,
                    )
                    designs[sort_key] = s
    best_key = sorted(designs.keys())[0]
    return designs[best_key]
