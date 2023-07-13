from .design import Design
from .types import Replicate, Pool, Peptide, Mapping
from .util import pairs_to_dict


def violations_per_replicate(
    s: Design,
    error_on_duplicate=True,
    error_on_extra=True,
    error_on_missing=True,
    verbose=False,
) -> Mapping[Replicate, int]:
    replicate_to_violations = {}

    # treat invalid pairs as if they've already been neighbors in a previous round
    peptide_to_neighbors = pairs_to_dict(s.invalid_neighbors)

    for replicate_idx, pool_to_peptides in s.assignments.items():
        violations = 0
        replicate_num = replicate_idx + 1

        # first check each peptide occurs once per replicate
        seen_peptides = set()
        for peptides in pool_to_peptides.values():
            for p in peptides:
                if p in seen_peptides:
                    if error_on_duplicate:
                        raise ValueError(
                            "Peptide %s twice in replicate %s" % (p, replicate_num)
                        )
                    violations += 1
                seen_peptides.add(p)
        expected_peptides = set(range(s.num_peptides))
        if seen_peptides != expected_peptides:
            extra_peptides = seen_peptides.difference(expected_peptides)
            if extra_peptides:
                if error_on_extra:
                    raise ValueError(
                        "Unexpected extra peptides in replicate %d: %s"
                        % (
                            replicate_num,
                            extra_peptides,
                        )
                    )
                violations += len(extra_peptides)
            missing_peptides = expected_peptides.difference(seen_peptides)
            if missing_peptides:
                if error_on_missing:
                    raise ValueError(
                        "Unexpected missing peptides in replicate %d: %s"
                        % (replicate_num, missing_peptides)
                    )
                violations += len(missing_peptides)

        # next check to make sure that each peptides only paired with another at most once
        for peptides in pool_to_peptides.values():
            for p1 in peptides:
                for p2 in peptides:
                    if p1 != p2:
                        if p2 in peptide_to_neighbors[p1]:
                            if verbose:
                                print(
                                    "Peptides %s and %s already together previous pool before replicate %s"
                                    % (p1, p2, replicate_num)
                                )
                            violations += 1
                        peptide_to_neighbors[p1].add(p2)
        replicate_to_violations[replicate_idx] = violations
    return replicate_to_violations


def count_violations(
    s: Design,
    error_on_duplicate=True,
    error_on_extra=True,
    error_on_missing=True,
    verbose=False,
) -> int:
    return sum(
        violations_per_replicate(
            s,
            error_on_duplicate=error_on_duplicate,
            error_on_extra=error_on_extra,
            error_on_missing=error_on_missing,
            verbose=verbose,
        ).values()
    )


def is_valid(
    s: Design,
    error_on_duplicate=True,
    error_on_extra=True,
    error_on_missing=True,
    verbose=False,
) -> bool:
    return (
        count_violations(
            s,
            error_on_duplicate=error_on_duplicate,
            error_on_extra=error_on_extra,
            error_on_missing=error_on_missing,
            verbose=verbose,
        )
        == 0
    )
