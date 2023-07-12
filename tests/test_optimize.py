from golfy import init, optimize, is_valid


def test_optimize():
    s = init(num_peptides=200, max_peptides_per_pool=10, num_replicates=3)
    assert not is_valid(s)
    optimize(s)
    assert is_valid(s)
