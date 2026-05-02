import numpy as np

from genomatch.bfile_utils import (
    BimRow,
    build_packed_ploidy_validation_plan,
    count_target_ploidy_genotype_issues,
    count_target_ploidy_genotype_issues_packed,
    decode_bed_chunk,
    encode_bed_row,
)


def test_packed_validation_parity_randomized():
    rng = np.random.default_rng(12345)
    row = BimRow("X", "rs1", "0", "5000000", "A", "G")
    ploidy_pairs = [(2, 2), (1, 2), (0, 2), (1, 1), (0, 0)]
    sample_sizes = [1, 2, 3, 4, 5, 7, 31, 64]

    for n_samples in sample_sizes:
        for _ in range(8):
            sexes = rng.integers(0, 3, size=n_samples, endpoint=False).tolist()
            genos = rng.integers(0, 4, size=n_samples, endpoint=False).tolist()
            packed = encode_bed_row(genos, n_samples)
            decoded = decode_bed_chunk(packed, n_samples)
            for ploidy_male, ploidy_female in ploidy_pairs:
                expected = count_target_ploidy_genotype_issues(
                    decoded,
                    sexes,
                    ploidy_male,
                    ploidy_female,
                    row,
                )
                plan = build_packed_ploidy_validation_plan(
                    sexes,
                    ploidy_male,
                    ploidy_female,
                )
                observed = count_target_ploidy_genotype_issues_packed(packed, plan)
                assert observed == expected


def test_packed_validation_ignores_padding_bits_in_partial_last_byte():
    # n_samples=5 uses two BED bytes; only the first slot in byte 2 is real data.
    n_samples = 5
    sexes = [1, 1, 1, 1, 1]
    ploidy_male, ploidy_female = (1, 1)
    genos = [1, 1, 1, 1, 1]
    packed = bytearray(encode_bed_row(genos, n_samples))
    # Keep slot0 of last byte missing (=01), poison padding slots 1..3 with het (=10).
    packed[-1] = 0b10101001
    packed_bytes = bytes(packed)
    decoded = decode_bed_chunk(packed_bytes, n_samples)
    expected = count_target_ploidy_genotype_issues(
        decoded,
        sexes,
        ploidy_male,
        ploidy_female,
        BimRow("MT", "rs1", "0", "10", "A", "G"),
    )
    plan = build_packed_ploidy_validation_plan(sexes, ploidy_male, ploidy_female)
    observed = count_target_ploidy_genotype_issues_packed(packed_bytes, plan)
    assert observed == expected
    assert observed == (0, 0, 0)


def test_packed_validation_unknown_sex_accounting_matches_reference():
    n_samples = 6
    sexes = [0, 0, 1, 2, 0, 1]
    ploidy_male, ploidy_female = (1, 2)
    genos = [2, 0, 2, 3, 1, 0]
    packed = encode_bed_row(genos, n_samples)
    decoded = decode_bed_chunk(packed, n_samples)
    expected = count_target_ploidy_genotype_issues(
        decoded,
        sexes,
        ploidy_male,
        ploidy_female,
        BimRow("X", "rs1", "0", "5000000", "A", "G"),
    )
    plan = build_packed_ploidy_validation_plan(sexes, ploidy_male, ploidy_female)
    observed = count_target_ploidy_genotype_issues_packed(packed, plan)
    assert observed == expected
    assert observed[2] == 3


def test_packed_validation_haploid_and_absent_issue_counts_match_reference():
    n_samples = 4
    genos = [2, 0, 1, 3]
    packed = encode_bed_row(genos, n_samples)
    decoded = decode_bed_chunk(packed, n_samples)
    sexes = [1, 1, 1, 1]
    row = BimRow("Y", "rs1", "0", "100", "A", "G")

    hap_expected = count_target_ploidy_genotype_issues(decoded, sexes, 1, 1, row)
    hap_plan = build_packed_ploidy_validation_plan(sexes, 1, 1)
    hap_observed = count_target_ploidy_genotype_issues_packed(packed, hap_plan)
    assert hap_observed == hap_expected
    assert hap_observed[0] == 1

    absent_expected = count_target_ploidy_genotype_issues(decoded, sexes, 0, 0, row)
    absent_plan = build_packed_ploidy_validation_plan(sexes, 0, 0)
    absent_observed = count_target_ploidy_genotype_issues_packed(packed, absent_plan)
    assert absent_observed == absent_expected
    assert absent_observed[1] == 3


def test_packed_validation_y_msy_haploid_and_absent_simultaneously():
    # Y_MSY: males haploid, females absent.
    sexes = [1, 2, 1, 2]
    genos = [2, 3, 0, 1]
    packed = encode_bed_row(genos, 4)
    decoded = decode_bed_chunk(packed, 4)
    row = BimRow("Y", "rs1", "0", "5000000", "A", "G")
    expected = count_target_ploidy_genotype_issues(decoded, sexes, 1, 0, row)
    plan = build_packed_ploidy_validation_plan(sexes, 1, 0)
    observed = count_target_ploidy_genotype_issues_packed(packed, plan)
    assert observed == expected
    assert observed == (1, 1, 0)
