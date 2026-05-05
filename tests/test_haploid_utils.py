from genomatch.haploid_utils import expected_ploidy_pair


def test_expected_ploidy_pair_uses_t2t_boundaries():
    build = "T2T-CHM13v2.0"

    assert expected_ploidy_pair("X", "1", genome_build=build) == (2, 2)
    assert expected_ploidy_pair("X", "2394410", genome_build=build) == (2, 2)
    assert expected_ploidy_pair("X", "2394411", genome_build=build) == (1, 2)
    assert expected_ploidy_pair("X", "153925834", genome_build=build) == (1, 2)
    assert expected_ploidy_pair("X", "153925835", genome_build=build) == (2, 2)
    assert expected_ploidy_pair("X", "154259566", genome_build=build) == (2, 2)
    assert expected_ploidy_pair("X", "154259567", genome_build=build) == (1, 2)

    assert expected_ploidy_pair("Y", "2458320", genome_build=build) == (2, 2)
    assert expected_ploidy_pair("Y", "2458321", genome_build=build) == (1, 0)
    assert expected_ploidy_pair("Y", "62122809", genome_build=build) == (1, 0)
    assert expected_ploidy_pair("Y", "62122810", genome_build=build) == (2, 2)
    assert expected_ploidy_pair("Y", "62460029", genome_build=build) == (2, 2)
    assert expected_ploidy_pair("Y", "62460030", genome_build=build) == (1, 0)

    assert expected_ploidy_pair("MT", "1", genome_build=build) == (1, 1)
