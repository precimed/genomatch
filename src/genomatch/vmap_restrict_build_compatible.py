#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

import pandas as pd

from .reference_utils import fetch_reference_bases, resolve_bcftools_binary, resolve_internal_reference_fasta
from .haploid_utils import expected_ploidy_pair
from .tabular_rows import VMapRowsTable, VariantRowsTable
from .vtable_utils import (
    VariantRow,
    complement_allele,
    compose_allele_ops_series,
    convert_contig_label,
    duplicate_target_rows_mask_table,
    load_variant_object_tables,
    normalize_contig_for_reference,
    require_contig_naming,
    require_table_matches_contig_naming,
    sort_target_table_by_declared_coordinate,
    validate_allele_values,
    write_metadata,
    write_vmap_status_qc,
    write_vmap_table,
    write_vtable_table,
)


NORM_STATUS_PRIORITY = (
    "norm_multiple_output_records",
    "norm_multiallelic",
    "norm_not_atcg_alleles",
    "norm_identical_ref_alt_alleles",
    "norm_invalid_position",
    "unsupported_target_contig",
    "ploidy_class_changed",
    "norm_unsupported_complex_indel",
    "norm_ref_mismatch",
)

ALLELE_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class NormalizationCandidate:
    candidate_id: str
    input_index: int
    row: VariantRow
    local_op: str


@dataclass(frozen=True)
class NormalizationOutcome:
    row: VariantRow | None
    local_op: str
    status: str | None = None


def validate_reference_aware_table(table: VariantRowsTable, *, label: str) -> None:
    """Validate input alleles once at boundary for restrict-build flows."""
    frame = table.to_frame(copy=False)
    validate_allele_values(frame["a1"], label=label)
    validate_allele_values(frame["a2"], label=label)


def restrict_rows_table(
    table: VariantRowsTable,
    fasta_path: Path,
    contig_naming: str,
    allow_strand_flips: bool,
) -> pd.DataFrame:
    # Assumes: input table is validated for allele token domain and required schema columns.
    # Performs: SV(reference-anchored restriction + optional strand-flip matching) in vectorized form.
    # Guarantees: one outcome row per input row with keep/local_op and canonical output allele columns.
    frame = table.to_frame(copy=False).loc[:, ["chrom", "pos", "id", "a1", "a2"]].copy()
    frame["pos_int"] = pd.to_numeric(frame["pos"], errors="raise").astype("int64")

    chrom_tokens = frame["chrom"].astype(str)
    contig_map: Dict[str, str | None] = {}
    # PERF: loop retained over unique contigs only (small cardinality), not per-row.
    for chrom in chrom_tokens.drop_duplicates().tolist():
        try:
            contig_map[chrom] = normalize_contig_for_reference(str(chrom), contig_naming, "ucsc")
        except ValueError:
            contig_map[chrom] = None
    frame["ucsc_contig"] = chrom_tokens.map(contig_map)

    valid_contig_mask = frame["ucsc_contig"].notna()
    query_pairs = list(
        dict.fromkeys(
            zip(
                frame.loc[valid_contig_mask, "ucsc_contig"].astype(str).tolist(),
                frame.loc[valid_contig_mask, "pos_int"].astype(int).tolist(),
            )
        )
    )
    reference_bases = fetch_reference_bases(fasta_path, query_pairs)
    ref_frame = pd.DataFrame(
        [
            {"ucsc_contig": contig, "pos_int": pos_int, "ref_base": base}
            for (contig, pos_int), base in reference_bases.items()
        ]
    )
    if ref_frame.empty:
        frame["ref_base"] = ""
    else:
        frame = frame.merge(ref_frame, on=["ucsc_contig", "pos_int"], how="left")
        frame["ref_base"] = frame["ref_base"].fillna("")

    frame["keep"] = False
    frame["local_op"] = "identity"
    frame["out_chrom"] = frame["chrom"]
    frame["out_pos"] = frame["pos"]
    frame["out_id"] = frame["id"]
    frame["out_a1"] = frame["a1"]
    frame["out_a2"] = frame["a2"]

    ref_ok = frame["ref_base"].isin(["A", "C", "G", "T"])
    a1_len1 = frame["a1"].str.len().eq(1)
    a2_len1 = frame["a2"].str.len().eq(1)

    identity_mask = ref_ok & a2_len1 & frame["a2"].eq(frame["ref_base"])
    swap_mask = ref_ok & ~identity_mask & a1_len1 & frame["a1"].eq(frame["ref_base"])

    frame.loc[identity_mask, "keep"] = True
    frame.loc[identity_mask, "local_op"] = "identity"

    frame.loc[swap_mask, "keep"] = True
    frame.loc[swap_mask, "local_op"] = "swap"
    frame.loc[swap_mask, ["out_a1", "out_a2"]] = frame.loc[swap_mask, ["a2", "a1"]].to_numpy()

    if allow_strand_flips:
        unresolved_mask = ref_ok & ~(identity_mask | swap_mask)
        full_acgt_mask = frame["a1"].str.fullmatch(r"[ACGT]+") & frame["a2"].str.fullmatch(r"[ACGT]+")
        complement_mask = unresolved_mask & full_acgt_mask
        comp_a1 = frame.loc[complement_mask, "a1"].str.translate(ALLELE_COMPLEMENT_TABLE).str[::-1]
        comp_a2 = frame.loc[complement_mask, "a2"].str.translate(ALLELE_COMPLEMENT_TABLE).str[::-1]
        comp_ref = frame.loc[complement_mask, "ref_base"]
        comp_identity_mask = comp_a2.str.len().eq(1) & comp_a2.eq(comp_ref)
        comp_swap_mask = (~comp_identity_mask) & comp_a1.str.len().eq(1) & comp_a1.eq(comp_ref)

        comp_identity_idx = comp_identity_mask[comp_identity_mask].index
        if not comp_identity_idx.empty:
            frame.loc[comp_identity_idx, "keep"] = True
            frame.loc[comp_identity_idx, "local_op"] = "flip"
            frame.loc[comp_identity_idx, "out_a1"] = comp_a1.loc[comp_identity_idx].to_numpy()
            frame.loc[comp_identity_idx, "out_a2"] = comp_a2.loc[comp_identity_idx].to_numpy()

        comp_swap_idx = comp_swap_mask[comp_swap_mask].index
        if not comp_swap_idx.empty:
            frame.loc[comp_swap_idx, "keep"] = True
            frame.loc[comp_swap_idx, "local_op"] = "flip_swap"
            frame.loc[comp_swap_idx, "out_a1"] = comp_a2.loc[comp_swap_idx].to_numpy()
            frame.loc[comp_swap_idx, "out_a2"] = comp_a1.loc[comp_swap_idx].to_numpy()

    return frame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Restrict a .vtable/.vmap to same-build variants that are reference-compatible "
            "with the declared build, using UCSC-style internal reference FASTA with contig "
            "normalization as needed."
        )
    )
    parser.add_argument("--source", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Filtered output .vtable or .vmap")
    parser.add_argument(
        "--allow-strand-flips",
        action="store_true",
        help="Allow strand-complemented allele restriction against the declared build",
    )
    parser.add_argument(
        "--norm-indels",
        action="store_true",
        help="Normalize retained indels with bcftools norm while keeping canonical target ordering",
    )
    parser.add_argument(
        "--sort",
        action="store_true",
        help="Sort retained output rows into declared coordinate order",
    )
    return parser.parse_args()



def canonicalize_reference_anchored_row(row: VariantRow, ref_base: str) -> Tuple[VariantRow, str] | None:
    if len(row.a2) == 1 and row.a2 == ref_base:
        return row, "identity"
    if len(row.a1) == 1 and row.a1 == ref_base:
        return VariantRow(row.chrom, row.pos, row.id, row.a2, row.a1), "swap"
    return None


def restrict_row_against_reference(
    row: VariantRow,
    ref_base: str,
    *,
    allow_strand_flips: bool,
) -> Tuple[VariantRow, str] | None:
    canonical = canonicalize_reference_anchored_row(row, ref_base)
    if canonical is not None:
        return canonical
    if not allow_strand_flips:
        return None
    try:
        complemented = VariantRow(
            row.chrom,
            row.pos,
            row.id,
            complement_allele(row.a1),
            complement_allele(row.a2),
        )
    except ValueError:
        return None
    canonical = canonicalize_reference_anchored_row(complemented, ref_base)
    if canonical is None:
        return None
    restricted_row, local_op = canonical
    if local_op == "identity":
        return restricted_row, "flip"
    return restricted_row, "flip_swap"


def restrict_rows(
    rows: Sequence[VariantRow],
    fasta_path: Path,
    contig_naming: str,
    allow_strand_flips: bool,
) -> List[Tuple[VariantRow | None, str]]:
    # Assumes: input rows are already validated for allele token domain.
    # Performs: SV(reference-anchored restriction + optional strand-flip matching).
    # Guarantees: row-aligned restriction outcomes with local allele operation tags.
    row_list = list(rows)
    frame = pd.DataFrame(
        {
            "chrom": [row.chrom for row in row_list],
            "pos": [row.pos for row in row_list],
        }
    )

    positions = [int(row.pos) for row in row_list]
    frame["pos_int"] = positions

    ucsc_contigs: List[str | None] = []
    for chrom in frame["chrom"].tolist():
        try:
            ucsc_contigs.append(normalize_contig_for_reference(str(chrom), contig_naming, "ucsc"))
        except ValueError:
            ucsc_contigs.append(None)
    frame["ucsc_contig"] = ucsc_contigs

    query_pairs = [
        (ucsc_contig, pos)
        for ucsc_contig, pos in zip(frame["ucsc_contig"].tolist(), frame["pos_int"].tolist())
        if ucsc_contig is not None
    ]
    reference_bases = fetch_reference_bases(fasta_path, query_pairs)
    frame["ref_base"] = [
        reference_bases.get((ucsc_contig, pos), "") if ucsc_contig is not None else ""
        for ucsc_contig, pos in zip(frame["ucsc_contig"].tolist(), frame["pos_int"].tolist())
    ]

    out_rows: List[Tuple[VariantRow | None, str]] = []
    # PERF: loop retained; allele reconciliation logic is row-local with no cross-row state, so
    # vectorization would require encoding the decision tree in masks without clarity benefit;
    # expected cardinality matches caller's row count (SNV or indel subset, not full table).
    for row, ref_base in zip(row_list, frame["ref_base"].tolist()):
        if ref_base not in {"A", "C", "G", "T"}:
            out_rows.append((None, "identity"))
            continue
        restricted = restrict_row_against_reference(row, ref_base, allow_strand_flips=allow_strand_flips)
        out_rows.append(restricted if restricted is not None else (None, "identity"))
    return out_rows


def normalize_candidate_to_ucsc_vcf(candidate: NormalizationCandidate, contig_naming: str) -> Tuple[str, str, str, str, str]:
    chrom = convert_contig_label(candidate.row.chrom, contig_naming, "ucsc")
    return chrom, candidate.row.pos, candidate.candidate_id, candidate.row.a2, candidate.row.a1


def write_normalization_vcf(
    path: Path,
    candidates: Sequence[NormalizationCandidate],
    contig_naming: str,
) -> None:
    # Assumes: candidates are prevalidated normalization candidates.
    # Performs: CV(temporary VCF row-shape assembly for bcftools norm).
    # Guarantees: debug/input VCF consumable by bcftools norm.
    ucsc_contigs = sorted(
        {normalize_candidate_to_ucsc_vcf(candidate, contig_naming)[0] for candidate in candidates}
    )
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom in ucsc_contigs:
            handle.write(f"##contig=<ID={chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for candidate in candidates:
            chrom, pos, candidate_id, ref, alt = normalize_candidate_to_ucsc_vcf(candidate, contig_naming)
            handle.write("\t".join([chrom, pos, candidate_id, ref, alt, ".", "PASS", "."]) + "\n")


def write_bcftools_norm_log(log_path: Path, cmd: Sequence[str], result: subprocess.CompletedProcess[str]) -> None:
    log_path.write_text(
        "\n".join(
            [
                f"command: {shlex.join(cmd)}",
                f"returncode: {result.returncode}",
                "stdout:",
                result.stdout.rstrip("\n"),
                "stderr:",
                result.stderr.rstrip("\n"),
                "",
            ]
        ),
        encoding="utf-8",
    )


def append_bcftools_norm_log(
    log_path: Path,
    cmd: Sequence[str],
    result: subprocess.CompletedProcess[str],
    *,
    label: str | None = None,
    mode: str,
) -> None:
    lines: List[str] = []
    if label:
        lines.append(f"[{label}]")
    lines.extend(
        [
            f"command: {shlex.join(cmd)}",
            f"returncode: {result.returncode}",
            "stdout:",
            result.stdout.rstrip("\n"),
            "stderr:",
            result.stderr.rstrip("\n"),
            "",
        ]
    )
    with open(log_path, mode, encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))


def run_bcftools_norm_once(
    input_vcf: Path,
    output_vcf: Path,
    fasta_path: Path,
    *,
    check_mode: str,
) -> Tuple[List[str], subprocess.CompletedProcess[str]]:
    bcftools = resolve_bcftools_binary()
    cmd = [
        bcftools,
        "norm",
        "-f",
        str(fasta_path),
        "-c",
        check_mode,
        "-o",
        str(output_vcf),
        "-Ov",
        str(input_vcf),
    ]
    return cmd, subprocess.run(cmd, capture_output=True, text=True, check=False)


def parse_bcftools_check_ref_warnings(stderr_text: str) -> Set[Tuple[str, str, str, str]]:
    warned_records: Set[Tuple[str, str, str, str]] = set()
    for line in stderr_text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 5 and parts[0] == "REF_MISMATCH":
            warned_records.add((parts[1], parts[2], parts[3], parts[4]))
    return warned_records


def write_vcf_excluding_ref_mismatch_records(input_vcf: Path, output_vcf: Path, warned_records: Set[Tuple[str, str, str, str]]) -> None:
    with open(input_vcf, "r", encoding="utf-8", newline="") as src, open(
        output_vcf, "w", encoding="utf-8", newline="\n"
    ) as dst:
        for line in src:
            if not line.strip() or line.startswith("#"):
                dst.write(line)
                continue
            chrom, pos, _row_id, ref, alt, *_rest = line.rstrip("\n").split("\t")
            if (chrom, pos, ref, alt) in warned_records:
                continue
            dst.write(line)


def run_bcftools_norm_check_ref_x_workaround(input_vcf: Path, output_vcf: Path, fasta_path: Path, *, log_path: Path) -> None:
    # Work around upstream bcftools issue #2427: `bcftools norm -c x` may segfault when the
    # first non-header record has a REF mismatch. Emulate strict `-c x` filtering via two
    # `-c w` passes by warning on REF mismatches first, then excluding exactly those records.
    first_cmd, first_result = run_bcftools_norm_once(input_vcf, output_vcf, fasta_path, check_mode="w")
    append_bcftools_norm_log(log_path, first_cmd, first_result, label="check_ref_w_pass1", mode="w")
    if first_result.returncode != 0:
        text = (first_result.stderr or first_result.stdout).strip()
        if text:
            raise ValueError(f"bcftools norm failed: {text} (see {log_path})")
        raise ValueError(f"bcftools norm failed with exit code {first_result.returncode}; see {log_path}")
    warned_records = parse_bcftools_check_ref_warnings(first_result.stderr)
    with tempfile.TemporaryDirectory(prefix="restrict_build_compatible.issue2427.", dir=output_vcf.parent) as temp_dir_name:
        filtered_input = Path(temp_dir_name) / "bcftools_norm_input_x_filtered.vcf"
        write_vcf_excluding_ref_mismatch_records(input_vcf, filtered_input, warned_records)
        second_cmd, second_result = run_bcftools_norm_once(filtered_input, output_vcf, fasta_path, check_mode="w")
        append_bcftools_norm_log(log_path, second_cmd, second_result, label="check_ref_w_pass2", mode="a")
    if second_result.returncode != 0:
        text = (second_result.stderr or second_result.stdout).strip()
        if text:
            raise ValueError(f"bcftools norm failed: {text} (see {log_path})")
        raise ValueError(f"bcftools norm failed with exit code {second_result.returncode}; see {log_path}")


def run_bcftools_norm(input_vcf: Path, output_vcf: Path, fasta_path: Path, *, check_mode: str, log_path: Path) -> None:
    if check_mode == "x":
        run_bcftools_norm_check_ref_x_workaround(input_vcf, output_vcf, fasta_path, log_path=log_path)
        return
    cmd, result = run_bcftools_norm_once(input_vcf, output_vcf, fasta_path, check_mode=check_mode)
    write_bcftools_norm_log(log_path, cmd, result)
    if result.returncode != 0:
        text = (result.stderr or result.stdout).strip()
        if text:
            raise ValueError(f"bcftools norm failed: {text} (see {log_path})")
        raise ValueError(f"bcftools norm failed with exit code {result.returncode}; see {log_path}")


def normalization_debug_vcf_path(output_path: Path, *, check_mode: str, kind: str) -> Path:
    return output_path.with_name(output_path.name + f".bcftools_norm_{kind}_{check_mode}.vcf")


def normalization_debug_log_path(output_path: Path, *, check_mode: str) -> Path:
    return output_path.with_name(output_path.name + f".bcftools_norm_output_{check_mode}.log")


def is_normalized_acgt_allele(value: str) -> bool:
    token = value.strip()
    return bool(token) and set(token) <= {"A", "C", "G", "T"}


def parse_normalized_candidates(
    output_vcf: Path,
    candidate_lookup: Dict[str, NormalizationCandidate],
    final_contig_naming: str,
) -> Dict[str, NormalizationOutcome]:
    # Assumes: output_vcf is produced by bcftools norm for candidate_lookup inputs.
    # Performs: PV(VCF row-shape checks), SV(normalization outcome classification and filtering).
    # Guarantees: one outcome per candidate_id with row or explicit failure status.
    records_by_id: Dict[str, List[Tuple[str, str, str, str]]] = {candidate_id: [] for candidate_id in candidate_lookup}
    # PERF: loop retained to preserve boundary-local parsing semantics from bcftools VCF output.
    with open(output_vcf, "r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 5:
                raise ValueError(f"invalid normalized VCF row: {raw_line.strip()}")
            chrom, pos, candidate_id, ref, alt = parts[:5]
            if candidate_id not in candidate_lookup:
                raise ValueError(f"unexpected normalized row ID: {candidate_id}")
            records_by_id[candidate_id].append((chrom, pos, ref, alt))

    outcomes: Dict[str, NormalizationOutcome] = {}
    for candidate_id, candidate in candidate_lookup.items():
        records = records_by_id.get(candidate_id, [])
        if not records:
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "norm_ref_mismatch")
            continue
        if len(records) != 1:
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "norm_multiple_output_records")
            continue
        chrom, pos, ref, alt = records[0]
        if "," in alt:
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "norm_multiallelic")
            continue
        try:
            final_pos = int(pos)
        except ValueError:
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "norm_invalid_position")
            continue
        if final_pos <= 0:
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "norm_invalid_position")
            continue
        if ref == alt:
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "norm_identical_ref_alt_alleles")
            continue
        if not is_normalized_acgt_allele(ref) or not is_normalized_acgt_allele(alt):
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "norm_not_atcg_alleles")
            continue
        try:
            final_chrom = convert_contig_label(chrom, "ucsc", final_contig_naming)
        except ValueError:
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "unsupported_target_contig")
            continue
        canonical = VariantRow(final_chrom, str(final_pos), candidate.row.id, alt, ref)
        if len(canonical.a1) > 1 and len(canonical.a2) > 1:
            outcomes[candidate_id] = NormalizationOutcome(None, candidate.local_op, "norm_unsupported_complex_indel")
            continue
        outcomes[candidate_id] = NormalizationOutcome(canonical, candidate.local_op)
    return outcomes


def normalize_candidates(
    candidates: Sequence[NormalizationCandidate],
    *,
    fasta_path: Path,
    contig_naming: str,
    check_mode: str,
    output_path: Path,
) -> Dict[str, NormalizationOutcome]:
    if not candidates:
        return {}
    input_vcf = normalization_debug_vcf_path(output_path, check_mode=check_mode, kind="input")
    output_vcf = normalization_debug_vcf_path(output_path, check_mode=check_mode, kind="output")
    log_path = normalization_debug_log_path(output_path, check_mode=check_mode)
    write_normalization_vcf(input_vcf, candidates, contig_naming)
    run_bcftools_norm(input_vcf, output_vcf, fasta_path, check_mode=check_mode, log_path=log_path)
    return parse_normalized_candidates(
        output_vcf,
        {candidate.candidate_id: candidate for candidate in candidates},
        contig_naming,
    )


def pick_branch3_failure_status(statuses: Sequence[str]) -> str:
    if statuses and all(status == "norm_ref_mismatch" for status in statuses):
        return "norm_ref_mismatch"
    for status in NORM_STATUS_PRIORITY:
        if status in statuses:
            return status
    return "norm_ref_mismatch"


def apply_normalization_ploidy_filter(
    source_row: VariantRow,
    outcome: NormalizationOutcome,
    *,
    genome_build: str,
) -> NormalizationOutcome:
    if outcome.row is None:
        return outcome
    if source_row.chrom == outcome.row.chrom and source_row.pos == outcome.row.pos:
        return outcome
    source_ploidy_pair = expected_ploidy_pair(
        source_row.chrom,
        source_row.pos,
        genome_build=genome_build,
    )
    target_ploidy_pair = expected_ploidy_pair(
        outcome.row.chrom,
        outcome.row.pos,
        genome_build=genome_build,
    )
    if source_ploidy_pair != target_ploidy_pair:
        return NormalizationOutcome(None, outcome.local_op, "ploidy_class_changed")
    return outcome


def restrict_rows_with_indel_normalization(
    rows: Sequence[VariantRow],
    fasta_path: Path,
    contig_naming: str,
    *,
    allow_strand_flips: bool,
    genome_build: str,
    output_path: Path,
) -> List[NormalizationOutcome]:
    # Assumes: input rows are already validated for allele token domain.
    # Performs: SV(branching by allele lengths + normalization/restriction + ploidy filtering).
    # Guarantees: one normalization outcome per input row with canonical local-op/status encoding.
    outcomes: List[NormalizationOutcome] = [NormalizationOutcome(None, "identity") for _ in rows]
    snv_indices: List[int] = []
    one_multibase_indices: List[int] = []
    both_multibase_indices: List[int] = []
    for idx, row in enumerate(rows):
        a1_len = len(row.a1.strip())
        a2_len = len(row.a2.strip())
        if a1_len == 1 and a2_len == 1:
            snv_indices.append(idx)
        elif (a1_len == 1) != (a2_len == 1):
            one_multibase_indices.append(idx)
        else:
            both_multibase_indices.append(idx)

    if snv_indices:
        snv_results = restrict_rows(
            [rows[idx] for idx in snv_indices],
            fasta_path,
            contig_naming,
            allow_strand_flips,
        )
        for input_index, (restricted_row, local_op) in zip(snv_indices, snv_results):
            if restricted_row is None:
                outcomes[input_index] = NormalizationOutcome(None, local_op)
            else:
                outcomes[input_index] = NormalizationOutcome(restricted_row, local_op)

    branch2_candidates: List[NormalizationCandidate] = []
    if one_multibase_indices:
        branch2_results = restrict_rows(
            [rows[idx] for idx in one_multibase_indices],
            fasta_path,
            contig_naming,
            allow_strand_flips,
        )
        for input_index, (restricted_row, local_op) in zip(one_multibase_indices, branch2_results):
            if restricted_row is None:
                outcomes[input_index] = NormalizationOutcome(None, local_op)
                continue
            branch2_candidates.append(
                NormalizationCandidate(
                    candidate_id=f"branch2_row{input_index}",
                    input_index=input_index,
                    row=restricted_row,
                    local_op=local_op,
                )
            )
        branch2_outcomes = normalize_candidates(
            branch2_candidates,
            fasta_path=fasta_path,
            contig_naming=contig_naming,
            check_mode="e",
            output_path=output_path,
        )
        for candidate in branch2_candidates:
            parsed = branch2_outcomes[candidate.candidate_id]
            outcomes[candidate.input_index] = apply_normalization_ploidy_filter(
                rows[candidate.input_index],
                parsed,
                genome_build=genome_build,
            )

    branch3_candidates: List[NormalizationCandidate] = []
    if both_multibase_indices:
        for input_index in both_multibase_indices:
            row = rows[input_index]
            branch3_candidates.extend(
                [
                    NormalizationCandidate(
                        candidate_id=f"branch3_row{input_index}_identity",
                        input_index=input_index,
                        row=row,
                        local_op="identity",
                    ),
                    NormalizationCandidate(
                        candidate_id=f"branch3_row{input_index}_swap",
                        input_index=input_index,
                        row=VariantRow(row.chrom, row.pos, row.id, row.a2, row.a1),
                        local_op="swap",
                    ),
                ]
            )
        branch3_outcomes = normalize_candidates(
            branch3_candidates,
            fasta_path=fasta_path,
            contig_naming=contig_naming,
            check_mode="x",
            output_path=output_path,
        )
        grouped: Dict[int, List[NormalizationOutcome]] = {}
        for candidate in branch3_candidates:
            grouped.setdefault(candidate.input_index, []).append(
                apply_normalization_ploidy_filter(
                    rows[candidate.input_index],
                    branch3_outcomes[candidate.candidate_id],
                    genome_build=genome_build,
                )
            )
        for input_index, candidate_outcomes in grouped.items():
            survivors = [candidate for candidate in candidate_outcomes if candidate.row is not None]
            if len(survivors) == 1:
                outcomes[input_index] = survivors[0]
                continue
            if len(survivors) == 2:
                outcomes[input_index] = NormalizationOutcome(None, "identity", "norm_ambiguous_orientation")
                continue
            outcomes[input_index] = NormalizationOutcome(
                None,
                "identity",
                pick_branch3_failure_status(
                    [candidate.status for candidate in candidate_outcomes if candidate.status is not None]
                ),
            )
    return outcomes


def build_restriction_outcomes_table(
    source_rows_table: VariantRowsTable,
    *,
    fasta_path: Path,
    contig_naming: str,
    allow_strand_flips: bool,
    norm_indels: bool,
    genome_build: str,
    output_path: Path,
) -> pd.DataFrame:
    # Assumes: source rows already pass boundary contig and allele validation checks.
    # Performs: SV(split by SNV vs indel; vectorized SNV restriction; row-based indel restriction/normalization).
    # Guarantees: one aligned outcome record per input row with keep/local_op/status and canonical out_* fields.
    source_frame = (
        source_rows_table.to_frame(copy=False)
        .loc[:, ["chrom", "pos", "id", "a1", "a2"]]
        .reset_index(drop=True)
    )
    outcomes = pd.DataFrame(
        {
            "keep": pd.Series(False, index=source_frame.index, dtype=bool),
            "local_op": pd.Series("identity", index=source_frame.index, dtype="object"),
            "status": pd.Series([None] * len(source_frame), dtype="object"),
            "out_chrom": source_frame["chrom"],
            "out_pos": source_frame["pos"],
            "out_id": source_frame["id"],
            "out_a1": source_frame["a1"],
            "out_a2": source_frame["a2"],
        }
    )

    snv_mask = source_frame["a1"].str.len().eq(1) & source_frame["a2"].str.len().eq(1)
    snv_idx = source_frame.index[snv_mask]
    indel_idx = source_frame.index[~snv_mask]

    if not snv_idx.empty:
        snv_input = source_frame.loc[snv_idx, ["chrom", "pos", "id", "a1", "a2"]].reset_index(drop=True)
        snv_restricted = restrict_rows_table(
            VariantRowsTable.from_frame(snv_input, copy=False),
            fasta_path,
            contig_naming,
            allow_strand_flips,
        ).reset_index(drop=True)
        outcomes.loc[snv_idx, "keep"] = snv_restricted["keep"].to_numpy(dtype=bool)
        outcomes.loc[snv_idx, "local_op"] = snv_restricted["local_op"].to_numpy(dtype=object)
        outcomes.loc[snv_idx, "out_chrom"] = snv_restricted["out_chrom"].to_numpy(dtype=object)
        outcomes.loc[snv_idx, "out_pos"] = snv_restricted["out_pos"].to_numpy(dtype=object)
        outcomes.loc[snv_idx, "out_id"] = snv_restricted["out_id"].to_numpy(dtype=object)
        outcomes.loc[snv_idx, "out_a1"] = snv_restricted["out_a1"].to_numpy(dtype=object)
        outcomes.loc[snv_idx, "out_a2"] = snv_restricted["out_a2"].to_numpy(dtype=object)

    if not indel_idx.empty:
        indel_frame = source_frame.loc[indel_idx, ["chrom", "pos", "id", "a1", "a2"]].reset_index(drop=True)
        indel_rows: List[VariantRow] = []
        for chrom, pos, row_id, a1, a2 in indel_frame.itertuples(index=False, name=None):
            indel_rows.append(VariantRow(str(chrom), str(pos), str(row_id), str(a1), str(a2)))
        if norm_indels:
            indel_outcomes = restrict_rows_with_indel_normalization(
                indel_rows,
                fasta_path,
                contig_naming,
                allow_strand_flips=allow_strand_flips,
                genome_build=genome_build,
                output_path=output_path,
            )
        else:
            restricted = restrict_rows(indel_rows, fasta_path, contig_naming, allow_strand_flips)
            indel_outcomes = [NormalizationOutcome(row, local_op) for row, local_op in restricted]

        keep_values: List[bool] = []
        local_ops: List[str] = []
        statuses: List[str | None] = []
        out_chrom: List[str] = []
        out_pos: List[str] = []
        out_id: List[str] = []
        out_a1: List[str] = []
        out_a2: List[str] = []
        # PERF: loop retained for indel outcome projection; indel subset is expected smaller than SNV subset.
        for fallback_row, outcome in zip(indel_rows, indel_outcomes):
            local_ops.append(outcome.local_op)
            statuses.append(outcome.status)
            if outcome.row is None:
                keep_values.append(False)
                out_chrom.append(fallback_row.chrom)
                out_pos.append(fallback_row.pos)
                out_id.append(fallback_row.id)
                out_a1.append(fallback_row.a1)
                out_a2.append(fallback_row.a2)
                continue
            keep_values.append(True)
            out_chrom.append(outcome.row.chrom)
            out_pos.append(outcome.row.pos)
            out_id.append(outcome.row.id)
            out_a1.append(outcome.row.a1)
            out_a2.append(outcome.row.a2)

        outcomes.loc[indel_idx, "keep"] = keep_values
        outcomes.loc[indel_idx, "local_op"] = local_ops
        outcomes.loc[indel_idx, "status"] = statuses
        outcomes.loc[indel_idx, "out_chrom"] = out_chrom
        outcomes.loc[indel_idx, "out_pos"] = out_pos
        outcomes.loc[indel_idx, "out_id"] = out_id
        outcomes.loc[indel_idx, "out_a1"] = out_a1
        outcomes.loc[indel_idx, "out_a2"] = out_a2

    return outcomes


def main() -> int:
    # Assumes: source variant object path and metadata sidecar are valid and loadable.
    # Performs: orchestration of restriction/normalization SV kernels and output CV checks.
    # Guarantees: output object type matches input type with preserved metadata and QC sidecar behavior.
    args = parse_args()
    source_path = Path(args.source)
    output_path = Path(args.output)
    if not source_path.exists():
        raise ValueError(f"source not found: {source_path}")

    loaded = load_variant_object_tables(source_path)
    source_rows_table = loaded.target_rows_table
    input_row_count = len(source_rows_table)
    source_meta = dict(loaded.target_metadata)
    genome_build = source_meta.get("genome_build")
    contig_naming = require_contig_naming(source_meta, label="variant object")
    if genome_build == "unknown":
        raise ValueError(
            "source metadata genome_build must be known for build-compatible restriction "
            f"(declared genome_build={genome_build!r}, declared contig_naming={contig_naming!r}, "
            "internal_reference_naming='ucsc')"
        )
    if args.norm_indels and contig_naming == "plink_splitx":
        raise ValueError(
            "restrict_build_compatible.py --norm-indels does not support contig_naming=plink_splitx; "
            "normalize to build-independent ncbi/ucsc/plink first, then rerun --norm-indels"
        )
    require_table_matches_contig_naming(source_rows_table, contig_naming, label="variant object")
    validate_reference_aware_table(source_rows_table, label="restrict_build_compatible.py input")
    fasta_path = resolve_internal_reference_fasta(str(genome_build))
    duplicate_target_count = 0
    qc_output_path = output_path.with_name(output_path.name + ".qc.tsv")
    if args.norm_indels:
        for check_mode in ("e", "x"):
            for kind in ("input", "output"):
                debug_vcf = normalization_debug_vcf_path(output_path, check_mode=check_mode, kind=kind)
                if debug_vcf.exists():
                    debug_vcf.unlink()
            debug_log = normalization_debug_log_path(output_path, check_mode=check_mode)
            if debug_log.exists():
                debug_log.unlink()

    restricted_outcomes = build_restriction_outcomes_table(
        source_rows_table,
        fasta_path=fasta_path,
        contig_naming=contig_naming,
        allow_strand_flips=args.allow_strand_flips,
        norm_indels=bool(args.norm_indels),
        genome_build=str(genome_build),
        output_path=output_path,
    )
    keep_mask = restricted_outcomes["keep"].astype(bool)
    flip_count = int((keep_mask & restricted_outcomes["local_op"].isin(["flip", "flip_swap"])).sum())

    if loaded.base_vmap_table is not None:
        base_frame = loaded.base_vmap_table.to_frame(copy=False).reset_index(drop=True)
        if len(base_frame) != len(restricted_outcomes):
            raise ValueError("input variant_map source and target row counts are inconsistent")
        kept_idx = restricted_outcomes.index[keep_mask]
        out_vmap_frame = pd.DataFrame(
            {
                "chrom": restricted_outcomes.loc[kept_idx, "out_chrom"].reset_index(drop=True),
                "pos": restricted_outcomes.loc[kept_idx, "out_pos"].reset_index(drop=True),
                "id": restricted_outcomes.loc[kept_idx, "out_id"].reset_index(drop=True),
                "a1": restricted_outcomes.loc[kept_idx, "out_a1"].reset_index(drop=True),
                "a2": restricted_outcomes.loc[kept_idx, "out_a2"].reset_index(drop=True),
                "source_shard": base_frame.loc[kept_idx, "source_shard"].reset_index(drop=True),
                "source_index": pd.to_numeric(base_frame.loc[kept_idx, "source_index"], errors="raise")
                .astype("int64")
                .reset_index(drop=True),
                "allele_op": compose_allele_ops_series(
                    base_frame.loc[kept_idx, "allele_op"].reset_index(drop=True),
                    restricted_outcomes.loc[kept_idx, "local_op"].reset_index(drop=True),
                ),
            }
        )

        qc_frames: List[pd.DataFrame] = []
        status_mask = (~keep_mask) & restricted_outcomes["status"].notna()
        if bool(status_mask.any()):
            dropped_idx = restricted_outcomes.index[status_mask]
            qc_status = pd.DataFrame(
                {
                    "source_shard": base_frame.loc[dropped_idx, "source_shard"].reset_index(drop=True),
                    "source_index": pd.to_numeric(base_frame.loc[dropped_idx, "source_index"], errors="raise")
                    .astype("int64")
                    .reset_index(drop=True),
                    "id": base_frame.loc[dropped_idx, "id"].reset_index(drop=True),
                    "status": restricted_outcomes.loc[dropped_idx, "status"].reset_index(drop=True),
                }
            )
            qc_frames.append(qc_status)

        duplicate_mask = duplicate_target_rows_mask_table(out_vmap_frame)
        duplicate_target_count = int(duplicate_mask.sum())
        if bool(duplicate_mask.any()):
            qc_duplicate = out_vmap_frame.loc[duplicate_mask, ["source_shard", "source_index", "id"]].copy()
            qc_duplicate["status"] = "duplicate_target"
            qc_frames.append(qc_duplicate)
            out_vmap_frame = out_vmap_frame.loc[~duplicate_mask].reset_index(drop=True)

        if args.sort and not out_vmap_frame.empty:
            out_vmap_frame = sort_target_table_by_declared_coordinate(
                out_vmap_frame,
                contig_naming,
                label="restrict_build_compatible.py output",
            )

        if qc_frames:
            qc_all = pd.concat(qc_frames, axis=0, ignore_index=True)
            write_vmap_status_qc(
                qc_output_path,
                qc_all.loc[:, ["source_shard", "source_index", "id", "status"]].itertuples(index=False, name=None),
            )
        elif qc_output_path.exists():
            qc_output_path.unlink()

        write_vmap_table(output_path, VMapRowsTable.from_frame(out_vmap_frame, copy=False), assume_validated=True)
        output_row_count = len(out_vmap_frame)
    else:
        out_vtable_frame = (
            restricted_outcomes.loc[keep_mask, ["out_chrom", "out_pos", "out_id", "out_a1", "out_a2"]]
            .rename(
                columns={
                    "out_chrom": "chrom",
                    "out_pos": "pos",
                    "out_id": "id",
                    "out_a1": "a1",
                    "out_a2": "a2",
                }
            )
            .reset_index(drop=True)
        )
        if args.sort and not out_vtable_frame.empty:
            out_vtable_frame = sort_target_table_by_declared_coordinate(
                out_vtable_frame,
                contig_naming,
                label="restrict_build_compatible.py output",
            )
        if qc_output_path.exists():
            qc_output_path.unlink()
        write_vtable_table(output_path, VariantRowsTable.from_frame(out_vtable_frame, copy=False), assume_validated=True)
        output_row_count = len(out_vtable_frame)
    write_metadata(output_path, dict(loaded.raw_metadata))
    summary = {
        "input": str(source_path),
        "output": str(output_path),
        "object_type": loaded.object_type,
        "genome_build": str(genome_build),
        "contig_naming": contig_naming,
        "internal_reference_naming": "ucsc",
        "normalization": "none" if contig_naming == "ucsc" else f"{contig_naming}->ucsc",
        "allow_strand_flips": bool(args.allow_strand_flips),
        "norm_indels": bool(args.norm_indels),
        "input_rows": input_row_count,
        "output_rows": output_row_count,
        "dropped_rows": input_row_count - output_row_count,
        "duplicate_target_rows": duplicate_target_count,
        "flip_rows": flip_count,
        "reference_fasta": str(fasta_path),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
