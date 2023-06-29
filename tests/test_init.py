from golfy import random_init


def test_random_init():
    s = random_init(num_peptides=100, peptides_per_pool=5, num_replicates=3)
    assert len(s) == 3
    assert len(s[0]) == 20
    assert len(s[0][0]) == 5
