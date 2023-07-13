from golfy import find_best_design


def test_find_best_design_allow_extra_pools():
    s = find_best_design(100, 20, 3, allow_extra_pools=True)
    assert s.num_peptides == 100
    assert s.num_replicates == 3
    assert s.max_peptides_per_pool == 20
    assert s.num_pools() > 15


def test_find_best_design_no_extra_pools():
    s = find_best_design(100, 20, 3, allow_extra_pools=False)
    assert s.num_peptides == 100
    assert s.num_replicates == 3
    assert s.max_peptides_per_pool == 20
    assert s.num_pools() == 15
