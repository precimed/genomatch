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

from .reference_utils import fetch_reference_base, resolve_bcftools_binary, resolve_internal_reference_fasta
from .haploid_utils import expected_ploidy_pair
from .vtable_utils import (
    VMapRow,
    VariantRow,
    complement_allele,
    compose_allele_ops,
    convert_contig_label,
    duplicate_target_row_keys,
    load_variant_object,
    normalize_contig_for_reference,
    require_contig_naming,
    require_rows_match_contig_naming,
    sort_target_rows_by_declared_coordinate,
    target_row_key,
    validate_allele_value,
    write_metadata,
    write_vmap_status_qc,
    write_vmap,
    write_vtable,
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


def candidate_ref_base(row: VariantRow, fasta_path: Path, contig_naming: str) -> str:
    pos = int(row.pos)
    try:
        ucsc_contig = normalize_contig_for_reference(row.chrom, contig_naming, "ucsc")
    except ValueError:
        return ""
    return fetch_reference_base(fasta_path, ucsc_contig, pos)


def validate_reference_aware_alleles(row: VariantRow, *, label: str) -> None:
    validate_allele_value(row.a1, label=label)
    validate_allele_value(row.a2, label=label)


def is_reference_anchored(a1: str, a2: str, ref_base: str) -> bool:
    alleles = [a1.upper(), a2.upper()]
    return any(len(allele) == 1 and allele == ref_base for allele in alleles)


def canonicalize_reference_anchored_row(row: VariantRow, ref_base: str) -> Tuple[VariantRow, str] | None:
    if len(row.a2) == 1 and row.a2.upper() == ref_base:
        return row, "identity"
    if len(row.a1) == 1 and row.a1.upper() == ref_base:
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
    out_rows: List[Tuple[VariantRow | None, str]] = []
    for row in rows:
        validate_reference_aware_alleles(row, label="restrict_build_compatible.py input")
        ref_base = candidate_ref_base(row, fasta_path, contig_naming)
        if ref_base not in {"A", "C", "G", "T"}:
            out_rows.append((None, "identity"))
            continue
        restricted = restrict_row_against_reference(row, ref_base, allow_strand_flips=allow_strand_flips)
        if restricted is None:
            out_rows.append((None, "identity"))
            continue
        out_rows.append(restricted)
    return out_rows


def normalize_candidate_to_ucsc_vcf(candidate: NormalizationCandidate, contig_naming: str) -> Tuple[str, str, str, str, str]:
    chrom = convert_contig_label(candidate.row.chrom, contig_naming, "ucsc")
    return chrom, candidate.row.pos, candidate.candidate_id, candidate.row.a2.upper(), candidate.row.a1.upper()


def write_normalization_vcf(
    path: Path,
    candidates: Sequence[NormalizationCandidate],
    contig_naming: str,
) -> None:
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
    token = value.strip().upper()
    return bool(token) and set(token) <= {"A", "C", "G", "T"}


def parse_normalized_candidates(
    output_vcf: Path,
    candidate_lookup: Dict[str, NormalizationCandidate],
    final_contig_naming: str,
) -> Dict[str, NormalizationOutcome]:
    records_by_id: Dict[str, List[Tuple[str, str, str, str]]] = {candidate_id: [] for candidate_id in candidate_lookup}
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
        if ref.upper() == alt.upper():
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
        canonical = VariantRow(final_chrom, str(final_pos), candidate.row.id, alt.upper(), ref.upper())
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
    outcomes: List[NormalizationOutcome] = [NormalizationOutcome(None, "identity") for _ in rows]
    snv_indices: List[int] = []
    one_multibase_indices: List[int] = []
    both_multibase_indices: List[int] = []
    for idx, row in enumerate(rows):
        validate_reference_aware_alleles(row, label="restrict_build_compatible.py input")
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


def main() -> int:
    args = parse_args()
    source_path = Path(args.source)
    output_path = Path(args.output)
    if not source_path.exists():
        raise ValueError(f"source not found: {source_path}")

    loaded = load_variant_object(source_path)
    source_rows = loaded.target_rows
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
    require_rows_match_contig_naming(source_rows, contig_naming, label="variant object")
    fasta_path = resolve_internal_reference_fasta(str(genome_build))
    if args.norm_indels:
        for check_mode in ("e", "x"):
            for kind in ("input", "output"):
                debug_vcf = normalization_debug_vcf_path(output_path, check_mode=check_mode, kind=kind)
                if debug_vcf.exists():
                    debug_vcf.unlink()
            debug_log = normalization_debug_log_path(output_path, check_mode=check_mode)
            if debug_log.exists():
                debug_log.unlink()
        restricted_rows = restrict_rows_with_indel_normalization(
            source_rows,
            fasta_path,
            contig_naming,
            allow_strand_flips=args.allow_strand_flips,
            genome_build=str(genome_build),
            output_path=output_path,
        )
    else:
        restricted_rows = [
            NormalizationOutcome(restricted_row, local_op)
            for restricted_row, local_op in restrict_rows(source_rows, fasta_path, contig_naming, args.allow_strand_flips)
        ]
    flip_count = sum(
        1
        for outcome in restricted_rows
        if outcome.row is not None and outcome.local_op in {"flip", "flip_swap"}
    )
    duplicate_target_count = 0
    qc_output_path = output_path.with_name(output_path.name + ".qc.tsv")
    if loaded.base_vmap_rows is not None:
        base_rows = loaded.base_vmap_rows
        out_vmap_rows: List[VMapRow] = []
        qc_rows: List[tuple[str, int, str, str]] = []
        for base_row, outcome in zip(base_rows, restricted_rows):
            if outcome.row is None:
                if outcome.status is not None:
                    qc_rows.append((base_row.source_shard, base_row.source_index, base_row.id, outcome.status))
                continue
            out_vmap_rows.append(
                VMapRow(
                    outcome.row.chrom,
                    outcome.row.pos,
                    outcome.row.id,
                    outcome.row.a1,
                    outcome.row.a2,
                    base_row.source_shard,
                    base_row.source_index,
                    compose_allele_ops(base_row.allele_op, outcome.local_op),
                )
            )
        duplicate_keys = duplicate_target_row_keys(out_vmap_rows)
        duplicate_target_count = sum(
            1
            for row in out_vmap_rows
            if target_row_key(row) in duplicate_keys
        )
        if duplicate_keys:
            qc_rows.extend(
                [
                (row.source_shard, row.source_index, row.id, "duplicate_target")
                for row in out_vmap_rows
                if target_row_key(row) in duplicate_keys
                ]
            )
            out_vmap_rows = [
                row
                for row in out_vmap_rows
                if target_row_key(row) not in duplicate_keys
            ]
        if args.sort:
            out_vmap_rows = sort_target_rows_by_declared_coordinate(
                out_vmap_rows,
                contig_naming,
                label="restrict_build_compatible.py output",
            )
        if qc_rows:
            write_vmap_status_qc(qc_output_path, qc_rows)
        elif qc_output_path.exists():
            qc_output_path.unlink()
        write_vmap(output_path, out_vmap_rows)
        output_row_count = len(out_vmap_rows)
    else:
        out_rows = [outcome.row for outcome in restricted_rows if outcome.row is not None]
        if args.sort:
            out_rows = sort_target_rows_by_declared_coordinate(
                out_rows,
                contig_naming,
                label="restrict_build_compatible.py output",
            )
        write_vtable(output_path, out_rows)
        output_row_count = len(out_rows)
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
        "input_rows": len(source_rows),
        "output_rows": output_row_count,
        "dropped_rows": len(source_rows) - output_row_count,
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
