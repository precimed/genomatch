from __future__ import annotations

import json
from pathlib import Path

from real_liftover_helpers import resolve_bcftools_with_liftover, select_mapping_entries, write_real_match_config
from utils import read_tsv, run_py, run_py_with_env, write_json, write_lines


def write_sumstats_metadata(path: Path) -> None:
    write_lines(
        path,
        [
            "cleansumstats_metafile_kind: minimal",
            "path_sumStats: input.tsv",
            "stats_Model: linear",
            "col_CHR: chromosome",
            "col_POS: position",
            "col_SNP: snp_id",
            "col_EffectAllele: effect_allele",
            "col_OtherAllele: other_allele",
            "col_BETA: beta",
            "col_OR: odds_ratio",
            "col_EAF: effect_allele_frequency",
            "col_P: pvalue",
            "col_N: sample_size",
        ],
    )


def build_expected_output(entries, source_rows, target_rows):
    index_by_rsid = {entry.rsid: idx for idx, entry in enumerate(entries)}
    expected = [["chromosome", "position", "snp_id", "effect_allele", "other_allele", "beta", "odds_ratio", "effect_allele_frequency", "pvalue", "sample_size"]]
    for target in target_rows:
        rsid = target[2]
        if rsid == "missing_rs":
            expected.append(target + ["n/a", "n/a", "n/a", "n/a", "n/a"])
            continue
        source = list(source_rows[index_by_rsid[rsid]])
        source[0] = target[0]
        source[1] = target[1]
        source[2] = target[2]
        swap = (source[3], source[4]) != (target[3], target[4])
        source[3] = target[3]
        source[4] = target[4]
        if swap:
            source[5] = str(-float(source[5]))
            source[6] = str(1.0 / float(source[6]))
            source[7] = str(1.0 - float(source[7]))
        expected.append(source)
    return expected


def test_end_to_end_sumstats_liftover_match_apply_real(tmp_path: Path) -> None:
    entries = select_mapping_entries(max_rows=5)
    bcftools = resolve_bcftools_with_liftover()
    config = write_real_match_config(tmp_path, chroms={entry.chrom for entry in entries})
    env = {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": bcftools}

    raw = tmp_path / "input.tsv"
    metadata = tmp_path / "metadata.yaml"
    source_vmap = tmp_path / "source.vmap"
    normalized_vmap = tmp_path / "source.normalized.vmap"
    lifted = tmp_path / "lifted.vmap"
    target_vtable = tmp_path / "target.vtable"
    final_vmap = tmp_path / "final.vmap"
    output = tmp_path / "output.tsv"

    source_rows = []
    for idx, entry in enumerate(entries):
        beta = ["0.10", "-0.20", "0.30", "-0.40", "0.50"][idx]
        odds_ratio = ["1.10", "0.90", "1.20", "0.80", "1.30"][idx]
        eaf = ["0.20", "0.30", "0.40", "0.45", "0.15"][idx]
        source_rows.append(
            [
                entry.chrom,
                entry.bp37,
                entry.rsid,
                entry.a1,
                entry.a2,
                beta,
                odds_ratio,
                eaf,
                f"0.00{idx + 1}",
                str(10000 + idx),
            ]
        )
    write_lines(raw, ["\t".join(["chromosome", "position", "snp_id", "effect_allele", "other_allele", "beta", "odds_ratio", "effect_allele_frequency", "pvalue", "sample_size"])] + ["\t".join(row) for row in source_rows])
    write_sumstats_metadata(metadata)

    result = run_py("import_sumstats.py", "--input", raw, "--output", source_vmap, "--sumstats-metadata", metadata, "--genome-build", "unknown")
    assert result.returncode == 0, result.stderr
    assert all(row[5] == "." and row[6] == str(idx) for idx, row in enumerate(read_tsv(source_vmap)))

    result = run_py("normalize_contigs.py", "--input", source_vmap, "--output", normalized_vmap, "--to", "ncbi")
    assert result.returncode == 0, result.stderr

    result = run_py_with_env("guess_build.py", env, "--input", normalized_vmap, "--write", "--force")
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "GRCh37"

    result = run_py_with_env("liftover_build.py", env, "--input", normalized_vmap, "--output", lifted, "--target-build", "GRCh38")
    assert result.returncode == 0, result.stderr

    target_rows = []
    for idx, entry in enumerate(reversed(entries)):
        if idx % 2 == 1:
            target_rows.append([entry.chrom, entry.bp38, entry.rsid, entry.a2, entry.a1])
        else:
            target_rows.append([entry.chrom, entry.bp38, entry.rsid, entry.a1, entry.a2])
    target_rows.append([entries[0].chrom, str(int(entries[-1].bp38) + 1000), "missing_rs", "A", "C"])
    write_lines(target_vtable, ["\t".join(row) for row in target_rows])
    write_json(
        target_vtable.with_name(target_vtable.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"},
    )

    result = run_py("match_vmap_to_target.py", "--source", lifted, "--target", target_vtable, "--output", final_vmap)
    assert result.returncode == 0, result.stderr
    expected_vmap = []
    source_index_by_rsid = {entry.rsid: idx for idx, entry in enumerate(entries)}
    for idx, target in enumerate(target_rows[:-1]):
        source_index = source_index_by_rsid[target[2]]
        op = "swap" if idx % 2 == 1 else "identity"
        expected_vmap.append(target + [".", str(source_index), op])
    expected_vmap.append(target_rows[-1] + [".", "-1", "missing"])
    assert read_tsv(final_vmap) == expected_vmap

    result = run_py("apply_vmap_to_sumstats.py", "--input", raw, "--sumstats-metadata", metadata, "--vmap", final_vmap, "--output", output)
    assert result.returncode == 0, result.stderr
    assert read_tsv(output) == build_expected_output(entries, source_rows, target_rows)
