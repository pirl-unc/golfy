import math
from typing import Optional, Literal

import numpy as np

from .design import Design, Spec
from .evaluation import evaluate_design, EvaluationResult
from .initialization import init
from .optimization import optimize
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


def generate_candidate_specs_for_pool_budget(
    num_peptides: int = 100,
    max_pools: int = 96,
    min_num_replicates: int = 2,
    max_num_replicates: int = 5,
    min_max_peptides_per_pool: int = 2,
    max_max_peptides_per_pool: Optional[int] = None,
    allow_extra_pools: Optional[bool] = None,
) -> list[Spec]:
    """
    Generator of candidate configurations from which we might want to try
    to make Designs that fit into the given pool budget.
    """
    if allow_extra_pools is None:
        allow_extra_pools = [True, False]
    elif type(allow_extra_pools) is bool:
        allow_extra_pools = [allow_extra_pools]

    if max_max_peptides_per_pool is None:
        max_max_peptides_per_pool = max(3, num_peptides // 5)

    if min_num_replicates > max_num_replicates:
        max_num_replicates = min_num_replicates

    if min_max_peptides_per_pool > max_max_peptides_per_pool:
        max_max_peptides_per_pool = min_max_peptides_per_pool

    for max_peptides_per_pool in range(
        min_max_peptides_per_pool, max_max_peptides_per_pool
    ):
        if max_peptides_per_pool * max_pools < num_peptides:
            # not enough pools to fit all peptides without replicates
            continue

        for num_replicates in range(min_num_replicates, max_num_replicates + 1):
            min_pools_for_spec = math.ceil(
                num_peptides * num_replicates / max_peptides_per_pool
            )
            if min_pools_for_spec > max_pools:
                # not enough pools to fit all peptides with given number of replicates
                continue

            for curr_allow_extra_pools in allow_extra_pools:
                yield Spec(
                    num_peptides=num_peptides,
                    max_peptides_per_pool=max_peptides_per_pool,
                    num_replicates=num_replicates,
                    allow_extra_pools=curr_allow_extra_pools,
                    invalid_neighbors=[],
                    preferred_neighbors=[],
                )


def score_designs_for_pool_budget(
    num_peptides: int = 100,
    max_pools: int = 96,
    num_simulation_iters: int = 3,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    min_pools: Optional[int] = None,
    min_num_replicates: int = 2,
    max_num_replicates: int = 5,
    min_max_peptides_per_pool: int = 2,
    max_max_peptides_per_pool: Optional[int] = None,
    allow_extra_pools: Optional[bool] = None,
    verbose: bool = False,
) -> list[tuple[Design, EvaluationResult]]:
    assert num_peptides > 1, "No need to pool if there's only one peptide"
    assert max_pools > 1, "Must have more than one pool"
    assert max_pools <= num_peptides, "Can't have more pools than peptides"

    designs = []
    specs = list(
        generate_candidate_specs_for_pool_budget(
            num_peptides=num_peptides,
            max_pools=max_pools,
            min_num_replicates=min_num_replicates,
            max_num_replicates=max_num_replicates,
            min_max_peptides_per_pool=min_max_peptides_per_pool,
            max_max_peptides_per_pool=max_max_peptides_per_pool,
            allow_extra_pools=allow_extra_pools,
        )
    )
    print("Generated %d candidate specs" % (len(specs),))

    shared_kwargs = dict(
        num_peptides=num_peptides,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        verbose=verbose,
    )

    for spec in specs:
        s = find_best_design(
            max_peptides_per_pool=spec.max_peptides_per_pool,
            num_replicates=spec.num_replicates,
            allow_extra_pools=spec.allow_extra_pools,
            **shared_kwargs,
        )
        num_pools = s.num_pools()
        if min_pools is not None and num_pools < min_pools:
            continue

        if num_pools > max_pools:
            print("%s: %d pools > %d max pools" % (spec, num_pools, max_pools))
        else:
            scores = evaluate_design(s, num_simulation_iters)
            print("%s: %s" % (spec, scores))
            designs.append((s, scores))
    return designs


def best_design_for_pool_budget(
    num_peptides: int = 100,
    max_pools: int = 96,
    num_simulation_iters: int = 3,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    min_pools: Optional[int] = None,
    min_num_replicates: int = 3,
    max_num_replicates: int = 8,
    min_max_peptides_per_pool: int = 2,
    max_max_peptides_per_pool: Optional[int] = None,
    verbose: bool = False,
):
    assert num_peptides > 1, "No need to pool if there's only one peptide"
    assert max_pools > 1, "Must have more than one pool"
    assert max_pools <= num_peptides, "Can't have more pools than peptides"
    designs_with_scores = score_designs_for_pool_budget(
        num_peptides=num_peptides,
        max_pools=max_pools,
        num_simulation_iters=num_simulation_iters,
        invalid_neighbors=invalid_neighbors,
        preferred_neighbors=preferred_neighbors,
        min_pools=min_pools,
        min_num_replicates=min_num_replicates,
        max_num_replicates=max_num_replicates,
        min_max_peptides_per_pool=min_max_peptides_per_pool,
        max_max_peptides_per_pool=max_max_peptides_per_pool,
        verbose=verbose,
    )
    best_pair = sorted(designs_with_scores, key=lambda x: x[1].sort_key())[0]
    return best_pair[0]
