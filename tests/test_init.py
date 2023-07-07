from golfy import init


def test_random_init():
    s = init(num_peptides=100, peptides_per_pool=5, num_replicates=3, strategy="random")
    assert len(s.assignments) == 3
    assert len(s.assignments[0]) == 20
    assert len(s.assignments[0][0]) == 5


def test_greedy_init():
    s = init(num_peptides=100, peptides_per_pool=5, num_replicates=3, strategy="greedy")
    assert len(s.assignments) == 3
    assert len(s.assignments[0]) == 20
    assert len(s.assignments[0][0]) == 5


def test_greedy_init_with_preferred_neighbors():
    """
    Test case from Andy:
    25 peptides / 5 peptides per pool / 3x
        2023-07-06 17:57:30 INFO     peptide_7 and peptide_10
        2023-07-06 17:57:30 INFO     peptide_12 and peptide_17
        2023-07-06 17:57:30 INFO     peptide_14 and peptide_15
        2023-07-06 17:57:30 INFO     peptide_15 and peptide_17
    """
    s = init(
        num_peptides=25,
        peptides_per_pool=5,
        num_replicates=3,
        strategy="greedy",
        preferred_neighbors=[(7, 10), (12, 17), (14, 15), (15, 17)],
        verbose=True,
    )
    assert len(s.assignments) == 3
    assert len(s.assignments[0]) == 5
    assert len(s.assignments[0][0]) == 5
    # expect that greedy fill of first replicate groups {7, 10} and {12, 17, 14, 15}
    pool_to_peptides = s.assignments[0]
    print(pool_to_peptides)

    peptide_to_pool = s.peptide_to_pool_dict_for_replicate(0)
    print(peptide_to_pool)

    assert peptide_to_pool[7] == peptide_to_pool[10]
    assert peptide_to_pool[12] == peptide_to_pool[17]
    assert peptide_to_pool[14] == peptide_to_pool[15]
    assert peptide_to_pool[15] == peptide_to_pool[17]
