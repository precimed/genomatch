from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List

import pytest

from genomatch.vtable_utils import classify_allele_operation
from utils import REPO_ROOT, read_tsv, run_py, write_json, write_lines


REGRESSION_DIR = REPO_ROOT / "tests" / "data" / "regressions"

CASE_CONFIGS = {
    "test_regression_280": {
        "source_build": "GRCh37",
        "target_build": "GRCh38",
        "allow_strand_flips": True,
        "column_map": {
            "SNP": "RSID",
            "A1": "EffectAllele",
            "A2": "OtherAllele",
            "EFFECT_A1": "B",
            "P": "P",
            "CaEAF": "CaseEAF",
        },
    },
    "test_regression_393": {
        "source_build": "GRCh37",
        "target_build": "GRCh38",
        "allow_strand_flips": True,
        "column_map": {
            "SNP": "RSID",
            "A1": "EffectAllele",
            "A2": "OtherAllele",
            "EFFECT_A1": "B",
            "P": "P",
            "Ncase": "CaseN",
            "Ncont": "ControlN",
            "CaseAF": "CaseEAF",
            "ControlAF": "ControlEAF",
        },
    },
    "test_regression_438": {
        "source_build": "GRCh38",
        "target_build": "GRCh38",
        "allow_strand_flips": False,
        "column_map": {
            "ID": "RSID",
            "ALLELE1": "EffectAllele",
            "ALLELE0": "OtherAllele",
            "BETA": "B",
            "SE": "SE",
            "N": "N",
            "A1FREQ": "EAF",
        },
    },
    "test_regression_missing_variants": {
        "source_build": "GRCh38",
        "target_build": "GRCh38",
        "allow_strand_flips": True,
        "column_map": {
            "A1": "EffectAllele",
            "A2": "OtherAllele",
            "EFFECT_A1": "B",
            "P": "P",
            "Ncase": "CaseN",
            "Ncont": "ControlN",
            "CaseAF": "CaseEAF",
        },
    },
}


def read_table(path: Path) -> List[List[str]]:
    return [line.rstrip("\n").split("\t") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def column_index(header: List[str], name: str) -> int:
    lowered = {col.lower(): idx for idx, col in enumerate(header)}
    idx = lowered.get(name.lower())
    if idx is None:
        raise AssertionError(f"missing column {name!r} in header {header}")
    return idx


def normalize_chrom_token(value: str) -> str:
    token = value.strip()
    return token[3:] if token.lower().startswith("chr") else token


def case_chroms(raw_rows: List[List[str]], raw_header: List[str], expected_rows: List[List[str]], expected_header: List[str]) -> set[str]:
    raw_chrom_idx = column_index(raw_header, "CHR") if "CHR" in raw_header else column_index(raw_header, "CHROM")
    expected_chrom_idx = column_index(expected_header, "CHR")
    chroms = {normalize_chrom_token(row[raw_chrom_idx]) for row in raw_rows}
    chroms.update({normalize_chrom_token(row[expected_chrom_idx]) for row in expected_rows})
    return chroms


def write_target_vtable(expected_header: List[str], expected_rows: List[List[str]], path: Path, *, genome_build: str) -> None:
    idx_chr = column_index(expected_header, "CHR")
    idx_pos = column_index(expected_header, "POS")
    idx_id = column_index(expected_header, "RSID")
    idx_a1 = column_index(expected_header, "EffectAllele")
    idx_a2 = column_index(expected_header, "OtherAllele")
    write_lines(
        path,
        [
            "\t".join([row[idx_chr], row[idx_pos], row[idx_id], row[idx_a1], row[idx_a2]])
            for row in expected_rows
        ],
    )
    write_json(
        path.with_name(path.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": genome_build, "contig_naming": "ncbi"},
    )


def read_metadata(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_source_provenance(
    case_name: str,
    raw_header: List[str],
    raw_rows: List[List[str]],
    expected_header: List[str],
    expected_row: List[str],
) -> tuple[str, int, str]:
    expected_a1 = expected_row[column_index(expected_header, "EffectAllele")]
    expected_a2 = expected_row[column_index(expected_header, "OtherAllele")]
    expected_chrom = normalize_chrom_token(expected_row[column_index(expected_header, "CHR")])
    expected_pos = expected_row[column_index(expected_header, "POS")]

    idx_a1 = column_index(raw_header, "A1") if "A1" in raw_header else column_index(raw_header, "ALLELE1")
    idx_a2 = column_index(raw_header, "A2") if "A2" in raw_header else column_index(raw_header, "ALLELE0")
    idx_chr = column_index(raw_header, "CHR") if "CHR" in raw_header else column_index(raw_header, "CHROM")
    idx_pos = column_index(raw_header, "BP") if "BP" in raw_header else column_index(raw_header, "GENPOS")

    idx_snp = None
    if case_name != "test_regression_missing_variants":
        idx_snp = column_index(raw_header, "SNP") if "SNP" in raw_header else column_index(raw_header, "ID")
        expected_rsid = expected_row[column_index(expected_header, "RSID")]
        for idx, raw_row in enumerate(raw_rows):
            if raw_row[idx_snp] != expected_rsid:
                continue
            op = classify_allele_operation(raw_row[idx_a1], raw_row[idx_a2], expected_a1, expected_a2, allow_strand_flips=True)
            if op != "missing":
                return ".", idx, op
        raise AssertionError(f"no source row found for expected rsid {expected_rsid!r}")

    for idx, raw_row in enumerate(raw_rows):
        raw_chrom = normalize_chrom_token(raw_row[idx_chr])
        if raw_chrom != expected_chrom:
            continue
        if raw_row[idx_pos] != expected_pos:
            continue
        op = classify_allele_operation(raw_row[idx_a1], raw_row[idx_a2], expected_a1, expected_a2, allow_strand_flips=True)
        if op != "missing":
            return ".", idx, op
    raise AssertionError(f"no source row found for expected row {expected_row}")


def write_bridge_vmap(
    case_name: str,
    raw_header: List[str],
    raw_rows: List[List[str]],
    expected_header: List[str],
    expected_rows: List[List[str]],
    path: Path,
    *,
    source_meta: Dict[str, object],
    target_build: str,
) -> None:
    idx_chr = column_index(expected_header, "CHR")
    idx_pos = column_index(expected_header, "POS")
    idx_rsid = column_index(expected_header, "RSID")
    idx_a1 = column_index(expected_header, "EffectAllele")
    idx_a2 = column_index(expected_header, "OtherAllele")

    rows = []
    for expected_row in expected_rows:
        source_shard, source_index, allele_op = resolve_source_provenance(
            case_name, raw_header, raw_rows, expected_header, expected_row
        )
        rows.append(
            "\t".join(
                [
                    expected_row[idx_chr],
                    expected_row[idx_pos],
                    expected_row[idx_rsid],
                    expected_row[idx_a1],
                    expected_row[idx_a2],
                    source_shard,
                    str(source_index),
                    allele_op,
                ]
            )
        )
    write_lines(path, rows)
    write_json(
        path.with_name(path.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": target_build, "contig_naming": "ncbi"},
        },
    )


def assert_equivalent(observed: str, expected: str) -> None:
    try:
        assert abs(float(observed) - float(expected)) < 1e-6
    except ValueError:
        assert observed == expected


def assert_selected_columns_match(
    output_rows: List[List[str]],
    output_header: List[str],
    expected_rows: List[List[str]],
    expected_header: List[str],
    column_map: Dict[str, str],
) -> None:
    assert len(output_rows) == len(expected_rows)
    for output_row, expected_row in zip(output_rows, expected_rows):
        for output_name, expected_name in column_map.items():
            assert_equivalent(
                output_row[column_index(output_header, output_name)],
                expected_row[column_index(expected_header, expected_name)],
            )


@pytest.mark.parametrize("case_name", list(CASE_CONFIGS))
def test_apply_vmap_to_sumstats_selected_regressions(tmp_path: Path, case_name: str) -> None:
    case = CASE_CONFIGS[case_name]
    raw = REGRESSION_DIR / case_name / "input.txt"
    metadata = REGRESSION_DIR / case_name / "metadata.yaml"
    expected = REGRESSION_DIR / case_name / "expected-result1.tsv"

    raw_table = read_table(raw)
    expected_table = read_table(expected)
    raw_header, raw_rows = raw_table[0], raw_table[1:]
    expected_header, expected_rows = expected_table[0], expected_table[1:]

    source_vmap = tmp_path / f"{case_name}.source.vmap"
    bridge_vmap = tmp_path / f"{case_name}.bridge.vmap"
    target_vtable = tmp_path / f"{case_name}.target.vtable"
    final_vmap = tmp_path / f"{case_name}.final.vmap"
    output = tmp_path / f"{case_name}.output.tsv"

    result = run_py(
        "import_sumstats.py",
        "--input",
        raw,
        "--output",
        source_vmap,
        "--sumstats-metadata",
        metadata,
        "--genome-build",
        case["source_build"],
    )
    assert result.returncode == 0, result.stderr
    assert all(row[5] == "." and row[6] == str(idx) for idx, row in enumerate(read_tsv(source_vmap)))
    source_meta = read_metadata(source_vmap.with_name(source_vmap.name + ".meta.json"))
    assert source_meta["object_type"] == "variant_map"
    assert source_meta["target"]["genome_build"] == case["source_build"]
    write_bridge_vmap(
        case_name,
        raw_header,
        raw_rows,
        expected_header,
        expected_rows,
        bridge_vmap,
        source_meta=source_meta,
        target_build=case["target_build"],
    )
    write_target_vtable(expected_header, expected_rows, target_vtable, genome_build=case["target_build"])
    result = run_py("match_vmap_to_target.py", "--source", bridge_vmap, "--target", target_vtable, "--output", final_vmap)
    assert result.returncode == 0, result.stderr
    assert sum(1 for row in read_tsv(final_vmap) if row[6] != "-1") == len(expected_rows)

    result = run_py("apply_vmap_to_sumstats.py", "--input", raw, "--sumstats-metadata", metadata, "--vmap", final_vmap, "--output", output)
    assert result.returncode == 0, result.stderr

    output_table = read_tsv(output)
    output_header, output_rows = output_table[0], output_table[1:]
    assert_selected_columns_match(output_rows, output_header, expected_rows, expected_header, case["column_map"])
