from __future__ import annotations

import json
from pathlib import Path

from real_liftover_helpers import resolve_bcftools_with_liftover, select_mapping_entries, write_real_match_config
from utils import read_tsv, run_py, run_py_with_env, write_json, write_lines


def format_output_float(value: float) -> str:
    return format(value, ".15g")


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
    expected = [["CHR", "POS", "SNP", "EffectAllele", "OtherAllele", "beta", "odds_ratio", "effect_allele_frequency", "pvalue", "sample_size"]]
    for target in target_rows:
        rsid = target[2]
        source = list(source_rows[index_by_rsid[rsid]])
        swap = (source[3], source[4]) != (target[3], target[4])
        output_row = target[:5] + source[5:]
        if swap:
            output_row[5] = format_output_float(-float(output_row[5]))
            output_row[6] = format_output_float(1.0 / float(output_row[6]))
            output_row[7] = format_output_float(1.0 - float(output_row[7]))
        expected.append(output_row)
    return expected


def test_end_to_end_sumstats_liftover_restrict_apply_real(tmp_path: Path) -> None:
    entries = select_mapping_entries(max_rows=5)
    bcftools = resolve_bcftools_with_liftover()
    config = write_real_match_config(tmp_path, chroms={entry.chrom for entry in entries})
    env = {
        "MATCH_CONFIG": str(config),
        "MATCH_BCFTOOLS": bcftools,
        "MATCH_REFERENCE_ACCESS_MODE": "LEGACY",
    }

    raw = tmp_path / "input.tsv"
    metadata = tmp_path / "metadata.yaml"
    source_vmap = tmp_path / "source.vmap"
    normalized_vmap = tmp_path / "source.normalized.vmap"
    lifted = tmp_path / "lifted.vmap"
    final_vmap = tmp_path / "final.vmap"
    output = tmp_path / "output.tsv"
    membership_vtable = tmp_path / "membership.vtable"

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

    lifted_rows = read_tsv(lifted)
    expected_vmap = [row for idx, row in enumerate(lifted_rows) if idx % 2 == 0]
    write_lines(membership_vtable, ["\t".join(row[:5]) for row in reversed(expected_vmap)])
    write_json(
        membership_vtable.with_name(membership_vtable.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"},
    )
    result = run_py("restrict_vmap.py", lifted, membership_vtable, "--output", final_vmap)
    assert result.returncode == 0, result.stderr
    assert read_tsv(final_vmap) == expected_vmap

    result = run_py(
        "project_payload.py",
        "--input-format",
        "sumstats",
        "--input",
        raw,
        "--sumstats-metadata",
        metadata,
        "--vmap",
        final_vmap,
        "--output",
        output,
        "--retain-snp-id",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(output) == build_expected_output(entries, source_rows, [row[:5] for row in expected_vmap])
