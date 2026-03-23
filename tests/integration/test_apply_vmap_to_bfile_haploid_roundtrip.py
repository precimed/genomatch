from __future__ import annotations

import random
from pathlib import Path

import pytest

from plink_vcf_helpers import (
    HaploidVariant,
    apply_vmap_to_bfile,
    assert_plink_ok,
    build_source_vmap,
    build_vmap,
    gt_to_dosage,
    load_haploid_variants,
    parse_vcf_genotypes,
    ploidy_for_region,
    read_variant_object_target_metadata,
    read_ploidy,
    recode_bfile_to_vcf,
    write_fam_sexes,
    write_target_vtable,
    write_vcf,
)
from utils import REPO_ROOT, read_bed_matrix, read_bim, resolve_plink, run_cmd


COUNTS = {
    "22": 2,
    "X_par1": 2,
    "X_par2": 2,
    "X_nonpar": 2,
    "Y_msy": 2,
    "MT": 2,
}


@pytest.fixture(params=["GRCh37", "GRCh38"])
def haploid_variants(request) -> tuple[str, list[HaploidVariant]]:
    bim_path = REPO_ROOT / "match-test" / "data" / "generated_reference" / "haploid_regions_target_500.bim"
    build_name = request.param
    return build_name, load_haploid_variants(bim_path, build_name=build_name, counts=COUNTS)


def make_source_variants(variants: list[HaploidVariant], build_name: str) -> list[HaploidVariant]:
    source_rows = []
    for variant in variants:
        source_rows.append(
            HaploidVariant(
                chrom="23" if variant.chrom == "X" else "24" if variant.chrom == "Y" else "26" if variant.chrom == "MT" else "22",
                snp=variant.snp,
                pos=variant.pos,
                a1=variant.a1,
                a2=variant.a2,
                region=variant.region,
            )
        )
    return source_rows


def make_target_rows(variants: list[HaploidVariant], source_chroms: list[str]):
    rows = []
    swap_mask = []
    for idx, (variant, source_chrom) in enumerate(zip(variants, source_chroms)):
        pos = str(variant.pos)
        swap = variant.region in ("X_nonpar", "MT") and idx % 2 == 0
        swap_mask.append(swap)
        if swap:
            rows.append([source_chrom, pos, variant.snp, variant.a2, variant.a1])
        else:
            rows.append([source_chrom, pos, variant.snp, variant.a1, variant.a2])
    return rows, swap_mask


def make_genotypes(variants: list[HaploidVariant], sexes: list[int], seed: int) -> list[list[str]]:
    rng = random.Random(seed)
    genotypes = []
    for variant in variants:
        row = []
        for sex in sexes:
            ploidy = ploidy_for_region(variant.region, sex)
            if ploidy == 0:
                row.append("./.")
                continue
            if rng.random() < 0.10:
                row.append("./.")
                continue
            if ploidy == 2:
                dose = rng.choice([0, 1, 2])
                row.append("0/0" if dose == 0 else "0/1" if dose == 1 else "1/1")
            else:
                allele = rng.choice([0, 1])
                row.append("0/0" if allele == 0 else "1/1")
        genotypes.append(row)
    return genotypes


def index_for_region(variants: list[HaploidVariant], region: str) -> int:
    for idx, variant in enumerate(variants):
        if variant.region == region:
            return idx
    raise AssertionError(f"region {region} not found")


def build_pipeline(
    tmp_path: Path,
    *,
    build_name: str,
    variants: list[HaploidVariant],
    sexes: list[int],
    genotypes: list[list[str]],
):
    plink_cmd = resolve_plink()
    samples = [f"S{idx:02d}" for idx in range(1, len(sexes) + 1)]
    source_variants = make_source_variants(variants, build_name)

    input_vcf = tmp_path / "input.vcf"
    source_prefix = tmp_path / "source"
    source_vmap = tmp_path / "source.vmap"
    target_vtable = tmp_path / "target.vtable"
    vmap_path = tmp_path / "aligned.vmap"
    out_prefix = tmp_path / "out"

    write_vcf(input_vcf, samples, source_variants, genotypes)
    result = run_cmd(
        plink_cmd
        + ["--vcf", str(input_vcf), "--double-id", "--keep-allele-order", "--make-bed", "--out", str(source_prefix)]
    )
    assert_plink_ok(
        result,
        allow_warnings=True,
        allowed_warning_substrings=["nonmale y chromosome genotype"],
    )

    write_fam_sexes(source_prefix.with_suffix(".fam"), sexes)
    build_source_vmap(source_prefix, genome_build=build_name, output_path=source_vmap)
    source_bim_rows = read_bim(source_prefix.with_suffix(".bim"))
    target_rows, swap_mask = make_target_rows(variants, [row.chrom for row in source_bim_rows])
    source_meta = read_variant_object_target_metadata(source_vmap)
    write_target_vtable(
        target_vtable,
        target_rows,
        genome_build=build_name,
        contig_naming=str(source_meta["contig_naming"]),
    )
    build_vmap(source_vmap, target_vtable, vmap_path)
    return plink_cmd, out_prefix, vmap_path, swap_mask, source_prefix


def test_apply_vmap_to_bfile_haploid_roundtrip_ploidy(tmp_path: Path, haploid_variants) -> None:
    build_name, variants = haploid_variants
    sexes = [1, 1, 1, 2, 2, 2]
    genotypes = make_genotypes(variants, sexes, seed=20240219)

    plink_cmd, out_prefix, _vmap_path, swap_mask, _source_prefix = build_pipeline(
        tmp_path,
        build_name=build_name,
        variants=variants,
        sexes=sexes,
        genotypes=genotypes,
    )

    result = apply_vmap_to_bfile(tmp_path / "source", tmp_path / "aligned.vmap", out_prefix)
    assert result.returncode == 0, result.stderr

    ploidy_rows = read_ploidy(out_prefix.with_suffix(".ploidy"))
    expected_ploidy = [(ploidy_for_region(variant.region, 1), ploidy_for_region(variant.region, 2)) for variant in variants]
    assert ploidy_rows == expected_ploidy

    recode_bfile_to_vcf(plink_cmd, out_prefix)
    _, out_alleles, out_genotypes = parse_vcf_genotypes(out_prefix.with_suffix(".vcf"))
    assert len(out_genotypes) == len(genotypes)
    for idx, (row_gts, swap) in enumerate(zip(out_genotypes, swap_mask)):
        expected_row = genotypes[idx]
        ref, alt = out_alleles[idx]
        target = variants[idx]
        target_a1 = target.a2 if swap else target.a1
        target_a2 = target.a1 if swap else target.a2
        assert {ref, alt} == {target_a1, target_a2}
        alt_is_a1 = alt == target_a1
        alt_is_a2 = alt == target_a2
        assert alt_is_a1 or alt_is_a2
        for obs_gt, exp_gt in zip(row_gts, expected_row):
            obs_dosage = gt_to_dosage(obs_gt)
            exp_dosage = gt_to_dosage(exp_gt)
            if exp_dosage is None:
                assert obs_dosage is None
                continue
            if alt_is_a1:
                expected_dosage = exp_dosage if not swap else 2 - exp_dosage
            else:
                expected_dosage = 2 - exp_dosage if not swap else exp_dosage
            assert obs_dosage == expected_dosage


@pytest.mark.parametrize(
    "case_name,expected_warning",
    [
        ("male_x_nonpar_het", "heterozygous haploid-target"),
        ("male_y_msy_het", "heterozygous haploid-target"),
        ("mt_male_het", "heterozygous haploid-target"),
        ("mt_female_het", "heterozygous haploid-target"),
        ("female_y_call", "absent-region nonmissing"),
    ],
)
def test_apply_vmap_to_bfile_haploid_incompatible_target_calls_are_reported_without_rewriting(
    tmp_path: Path,
    haploid_variants,
    case_name: str,
    expected_warning: str,
) -> None:
    build_name, variants = haploid_variants
    sexes = [1, 1, 2, 2]
    genotypes = make_genotypes(variants, sexes, seed=4242)

    if case_name == "male_x_nonpar_het":
        genotypes[index_for_region(variants, "X_nonpar")][0] = "0/1"
    elif case_name == "male_y_msy_het":
        genotypes[index_for_region(variants, "Y_msy")][0] = "0/1"
    elif case_name == "mt_male_het":
        genotypes[index_for_region(variants, "MT")][0] = "0/1"
    elif case_name == "mt_female_het":
        genotypes[index_for_region(variants, "MT")][2] = "0/1"
    elif case_name == "female_y_call":
        genotypes[index_for_region(variants, "Y_msy")][2] = "0/0"
    else:
        raise AssertionError(f"unknown case {case_name}")

    _plink_cmd, out_prefix, vmap_path, _swap_mask, source_prefix = build_pipeline(
        tmp_path,
        build_name=build_name,
        variants=variants,
        sexes=sexes,
        genotypes=genotypes,
    )

    result = apply_vmap_to_bfile(source_prefix, vmap_path, out_prefix)
    assert result.returncode == 0, result.stderr
    assert expected_warning in result.stderr
    matrix = read_bed_matrix(out_prefix.with_suffix(".bed"), len(sexes), len(variants))
    if case_name == "male_x_nonpar_het":
        assert matrix[index_for_region(variants, "X_nonpar")][0] == 2
    elif case_name == "male_y_msy_het":
        assert matrix[index_for_region(variants, "Y_msy")][0] == 2
    elif case_name == "mt_male_het":
        assert matrix[index_for_region(variants, "MT")][0] == 2
    elif case_name == "mt_female_het":
        assert matrix[index_for_region(variants, "MT")][2] == 2
    elif case_name == "female_y_call":
        assert matrix[index_for_region(variants, "Y_msy")][2] != 1
    else:
        raise AssertionError(f"unknown case {case_name}")


def test_apply_vmap_to_bfile_haploid_missing_sex_warns_and_skips_ploidy(tmp_path: Path, haploid_variants) -> None:
    build_name, variants = haploid_variants
    sexes = [0, 1, 2, 2]
    genotypes = make_genotypes(variants, sexes, seed=4242)

    _plink_cmd, out_prefix, vmap_path, _swap_mask, source_prefix = build_pipeline(
        tmp_path,
        build_name=build_name,
        variants=variants,
        sexes=sexes,
        genotypes=genotypes,
    )

    result = apply_vmap_to_bfile(source_prefix, vmap_path, out_prefix)
    assert result.returncode == 0, result.stderr
    assert "skipped ploidy validation" in result.stderr
    assert out_prefix.with_suffix(".bed").exists()
    assert out_prefix.with_suffix(".ploidy").exists()
