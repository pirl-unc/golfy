from golfy import best_design_for_pool_budget


def test_best_design_for_pool_budget():
    s = best_design_for_pool_budget(num_peptides=50, max_pools=20)
    assert s.num_peptides == 50
    assert s.num_replicates > 1
    assert s.num_pools() <= 20
