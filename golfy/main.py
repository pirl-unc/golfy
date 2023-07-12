from typing import Optional, Literal

from .initialization import init
from .optimization import optimize
from .solution import Solution
from .types import Replicate, Pool, Peptide, PeptidePairList
from .validity import count_violations


def find_best_solution(
    num_peptides: int = 100,
    max_peptides_per_pool: int = 5,
    num_replicates: int = 3,
    num_pools_per_replicate: Optional[int | dict[Replicate, int]] = None,
    invalid_neighbors: PeptidePairList = [],
    preferred_neighbors: PeptidePairList = [],
    allow_extra_pools: bool = False,
    verbose: bool = False,
) -> Solution:
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

    solutions = {
        strategy: init(strategy=strategy, **shared_kwargs)
        for strategy in init_strategies
    }
    best_solution = None
    best_violations = None
    best_num_pools = None

    for strategy, s in solutions.items():
        print(
            "Initial solution for init strategy '%s', violations=%d, num_pools=%d"
            % (strategy, count_violations(s), s.num_pools())
        )
        optimize(s, allow_extra_pools=allow_extra_pools, verbose=verbose)
        violations = count_violations(s)
        num_pools = s.num_pools()
        print(
            "After optimization strategy '%s', violations=%d, num_pools=%d"
            % (strategy, violations, num_pools)
        )
        if (
            best_solution is None
            or (violations < best_violations)
            or (violations == best_violations and num_pools < best_num_pools)
        ):
            best_solution = s
            best_violations = violations
            best_num_pools = num_pools
            print("^ new best solution")
    return best_solution
