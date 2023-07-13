from golfy.simulation import (
    simulate_elispot_counts,
    simulate_any_hits_per_pool,
    simulate_number_hits_per_pool,
)
from golfy import init


def test_simulate_elispot_counts_known_hits():
    # create 100 blocks of one peptide each
    s = init(20, 1, 5, strategy="singleton")
    # make the first two peptides hits
    hits = {0, 1}
    spot_counts, hits2 = simulate_elispot_counts(s, hit_peptides=hits)
    assert hits == hits2


def test_simulate_any_hits_per_pool_known_hits():
    # create 100 blocks of one peptide each
    s = init(20, 1, 5, strategy="singleton")
    hits = {0, 1}
    spot_counts, hits2 = simulate_any_hits_per_pool(s, hit_peptides=hits)
    assert hits == hits2


def test_simulate_number_hits_per_pool_known_hits():
    # create 100 blocks of one peptide each
    s = init(20, 1, 5, strategy="singleton")
    hits = {0, 1}
    spot_counts, hits2 = simulate_number_hits_per_pool(s, hit_peptides=hits)
    assert hits == hits2
