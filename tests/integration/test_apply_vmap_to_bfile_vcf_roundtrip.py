from __future__ import annotations

from pathlib import Path

import pytest

from plink_vcf_helpers import (
    apply_vmap_to_bfile,
    assert_plink_ok,
    build_source_vmap,
    build_vmap,
    gt_to_dosage,
    parse_vcf,
    read_variant_object_target_metadata,
    recode_bfile_to_vcf,
    write_random_vcf,
    write_target_vtable,
)
from utils import read_bim, read_tsv, resolve_plink, run_cmd


@pytest.mark.parametrize("n_samples", [20, 21, 22, 23])
def test_apply_vmap_to_bfile_vcf_roundtrip_order_and_swaps(tmp_path: Path, n_samples: int) -> None:
    plink_cmd = resolve_plink()
    input_vcf = tmp_path / "input.vcf"
    source_prefix = tmp_path / "source"
    source_vmap = tmp_path / "source.vmap"
    target_vtable = tmp_path / "target.vtable"
    vmap_path = tmp_path / "aligned.vmap"
    out_prefix = tmp_path / "out"

    write_random_vcf(input_vcf, n_samples=n_samples, n_snps=12, seed=12345 + n_samples, missing_rate=0.10)
    result = run_cmd(
        plink_cmd
        + ["--vcf", str(input_vcf), "--double-id", "--keep-allele-order", "--make-bed", "--out", str(source_prefix)]
    )
    assert_plink_ok(result)

    build_source_vmap(source_prefix, genome_build="GRCh37", output_path=source_vmap)
    source_meta = read_variant_object_target_metadata(source_vmap)
    source_rows = read_bim(source_prefix.with_suffix(".bim"))
    target_rows = []
    for idx, row in enumerate(reversed(source_rows)):
        if idx % 2 == 1:
            target_rows.append([row.chrom, row.bp, row.snp, row.a2, row.a1])
        else:
            target_rows.append([row.chrom, row.bp, row.snp, row.a1, row.a2])
    write_target_vtable(
        target_vtable,
        target_rows,
        genome_build="GRCh37",
        contig_naming=str(source_meta["contig_naming"]),
    )

    build_vmap(source_vmap, target_vtable, vmap_path)
    assert read_tsv(vmap_path) == [
        [
            row[0],
            row[1],
            row[2],
            row[3],
            row[4],
            ".",
            str(len(target_rows) - 1 - idx),
            "swap" if idx % 2 == 1 else "identity",
        ]
        for idx, row in enumerate(target_rows)
    ]

    result = apply_vmap_to_bfile(source_prefix, vmap_path, out_prefix)
    assert result.returncode == 0, result.stderr

    recode_bfile_to_vcf(plink_cmd, out_prefix)
    in_samples, in_variants = parse_vcf(input_vcf)
    out_samples, out_variants = parse_vcf(out_prefix.with_suffix(".vcf"))
    assert out_samples == [f"{sample}_{sample}" for sample in in_samples]

    reversed_in = list(reversed(in_variants))
    assert len(out_variants) == len(target_rows)
    for idx, observed in enumerate(out_variants):
        source_variant = reversed_in[idx]
        target = target_rows[idx]
        source_row = source_rows[len(source_rows) - 1 - idx]
        swap = idx % 2 == 1

        assert (observed.chrom, observed.pos, observed.vid) == (target[0], int(target[1]), target[2])
        assert {observed.ref, observed.alt} == {target[3], target[4]}
        alt_is_a1 = observed.alt == target[3]
        alt_is_a2 = observed.alt == target[4]
        assert alt_is_a1 or alt_is_a2

        for obs_gt, exp_gt in zip(observed.genotypes, source_variant.genotypes):
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
