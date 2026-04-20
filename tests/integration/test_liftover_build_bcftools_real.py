from __future__ import annotations

import json
from pathlib import Path

from real_liftover_helpers import (
    canonical_alleles_for_build,
    resolve_bcftools_with_liftover,
    select_mapping_entries,
    write_real_match_config,
)
from utils import read_tsv, run_py_with_env, write_json, write_lines


def expected_lifted_table_rows(entries):
    return [
        [entry.chrom, entry.bp38, entry.rsid, *canonical_alleles_for_build(entry, build_name="GRCh38")]
        for entry in sorted(entries, key=lambda item: (int(item.chrom), int(item.bp38)))
    ]


def expected_lifted_vmap_rows(entries):
    return [
        [
            entry.chrom,
            entry.bp38,
            entry.rsid,
            *canonical_alleles_for_build(entry, build_name="GRCh38"),
            f"orig{entry.chrom}",
            str(100 + idx),
            "identity",
        ]
        for idx, entry in sorted(enumerate(entries), key=lambda item: (int(item[1].chrom), int(item[1].bp38)))
    ]


def test_real_guess_liftover_and_validate_pipeline(tmp_path: Path) -> None:
    entries = select_mapping_entries(max_rows=4)
    bcftools = resolve_bcftools_with_liftover()
    config = write_real_match_config(tmp_path, chroms={entry.chrom for entry in entries})
    env = {
        "MATCH_CONFIG": str(config),
        "MATCH_BCFTOOLS": bcftools,
        "MATCH_REFERENCE_ACCESS_MODE": "BULK",
    }

    source = tmp_path / "source.vtable"
    lifted = tmp_path / "lifted.vtable"
    validated = tmp_path / "validated.vtable"
    write_lines(source, [f"{entry.chrom}\t{entry.bp37}\t{entry.rsid}\t{entry.a1}\t{entry.a2}" for entry in entries])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env("guess_build.py", env, "--input", source, "--write")
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "GRCh37"
    meta = json.loads(source.with_name(source.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["genome_build"] == "GRCh37"
    assert meta["contig_naming"] == "ncbi"

    result = run_py_with_env("liftover_build.py", env, "--input", source, "--output", lifted, "--target-build", "GRCh38")
    assert result.returncode == 0, result.stderr
    assert read_tsv(lifted) == expected_lifted_table_rows(entries)
    lifted_meta = json.loads(lifted.with_name(lifted.name + ".meta.json").read_text(encoding="utf-8"))
    assert lifted_meta == {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"}
    expected_qc = [["source_shard", "source_index", "source_id", "status"]]
    for idx, entry in enumerate(entries):
        expected_qc.append([".", str(idx), entry.rsid, "lifted"])
    assert read_tsv(lifted.with_name(lifted.name + ".qc.tsv")) == expected_qc

    result = run_py_with_env("restrict_build_compatible.py", env, "--source", lifted, "--output", validated)
    assert result.returncode == 0, result.stderr
    assert read_tsv(validated) == expected_lifted_table_rows(entries)


def test_real_liftover_composes_existing_vmap_provenance(tmp_path: Path) -> None:
    entries = select_mapping_entries(max_rows=4)
    bcftools = resolve_bcftools_with_liftover()
    config = write_real_match_config(tmp_path, chroms={entry.chrom for entry in entries})
    env = {
        "MATCH_CONFIG": str(config),
        "MATCH_BCFTOOLS": bcftools,
        "MATCH_REFERENCE_ACCESS_MODE": "LEGACY",
    }

    source = tmp_path / "base.vmap"
    lifted = tmp_path / "lifted.vmap"
    rows = []
    chrom_counts = {}
    for idx, entry in enumerate(entries):
        local_index = chrom_counts.get(entry.chrom, 0)
        chrom_counts[entry.chrom] = local_index + 1
        canonical_a1, canonical_a2 = canonical_alleles_for_build(entry, build_name="GRCh37")
        if idx % 2 == 1:
            rows.append(
                f"{entry.chrom}\t{entry.bp37}\t{entry.rsid}\t{canonical_a2}\t{canonical_a1}\torig{entry.chrom}\t{100 + idx}\tswap"
            )
        else:
            rows.append(
                f"{entry.chrom}\t{entry.bp37}\t{entry.rsid}\t{canonical_a1}\t{canonical_a2}\torig{entry.chrom}\t{100 + idx}\tidentity"
            )
    write_lines(source, rows)
    write_json(
        source.with_name(source.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"},
        },
    )

    result = run_py_with_env("liftover_build.py", env, "--input", source, "--output", lifted, "--target-build", "GRCh38")
    assert result.returncode == 0, result.stderr

    assert read_tsv(lifted) == expected_lifted_vmap_rows(entries)

    meta = json.loads(lifted.with_name(lifted.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["target"] == {"genome_build": "GRCh38", "contig_naming": "ncbi"}
