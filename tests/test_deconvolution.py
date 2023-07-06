from golfy import init, optimize, is_valid
from golfy.deconvolution import (
    simulate_elispot_counts,
    create_linear_system,
    solve_linear_system,
)


def test_deconvolution():
    s = init(100, 3, 5)
    optimize(s)
    assert is_valid(s)
    # make sure the margin between background and hit activities
    # is so large that we never falsely identify a background peptide
    counts, hit_peptides = simulate_elispot_counts(
        s,
        num_hits=3,
        max_background_peptide_activity=1,
        min_hit_peptide_activity=100,
    )

    linear_system = create_linear_system(s, counts)
    for loo in [False, True]:
        # solver should work regardless of whether we use LOO to estimate
        # probabilities of hits or not
        solution = solve_linear_system(linear_system, leave_on_out=loo)
        assert solution.high_confidence_hits == hit_peptides
