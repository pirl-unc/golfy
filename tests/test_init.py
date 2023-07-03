from golfy import init


def test_random_init():
    s = init(num_peptides=100, peptides_per_pool=5, num_replicates=3)
    assert len(s.assignments) == 3
    assert len(s.assignments[0]) == 20
    assert len(s.assignments[0][0]) == 5
