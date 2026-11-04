#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import pandas as pd

from ._cli_utils import run_cli
from .reference_utils import fetch_reference_bases, reference_alleles_match, resolve_bcftools_binary, resolve_internal_reference_fasta
from .haploid_utils import expected_ploidy_pair
from .tabular_rows import VMapRowsTable, VariantRowsTable
from .vtable_utils import (
    VariantRow,
    compose_allele_ops_series,
    convert_contig_label,
    duplicate_target_rows_mask_table,
    load_variant_object_tables,
    metadata_with_variants_count,
    normalize_contig_for_reference,
    require_contig_naming,
    require_table_matches_contig_naming,
    validate_allele_values,
    variant_row_qc_rows,
    write_metadata,
    write_variant_status_qc,
    write_vmap_table,
    write_vtable_table,
)


NORMALIZATION_BATCH_ROWS = 100_000

ALLELE_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")
logger = logging.getLogger(__name__)


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
    frame.loc[swap_mask, ["out_a1", "out_a2"]] = frame.loc[swap_mask, ["a2", "a1"]].to_numpy(dtype=object)

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
            frame.loc[comp_identity_idx, "out_a1"] = comp_a1.loc[comp_identity_idx].to_numpy(dtype=object)
            frame.loc[comp_identity_idx, "out_a2"] = comp_a2.loc[comp_identity_idx].to_numpy(dtype=object)

        comp_swap_idx = comp_swap_mask[comp_swap_mask].index
        if not comp_swap_idx.empty:
            frame.loc[comp_swap_idx, "keep"] = True
            frame.loc[comp_swap_idx, "local_op"] = "flip_swap"
            frame.loc[comp_swap_idx, "out_a1"] = comp_a2.loc[comp_swap_idx].to_numpy(dtype=object)
            frame.loc[comp_swap_idx, "out_a2"] = comp_a1.loc[comp_swap_idx].to_numpy(dtype=object)

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
    parser.add_argument("--include-indels", action="store_true", help="Admit multibase input (default: SNV-only)")
    parser.add_argument("--allow-indel-swaps", action="store_true", help="Allow REF/ALT exchange for multibase input")
    parser.add_argument(
        "--norm-indels",
        action="store_true",
        help="Normalize retained indels with bcftools norm while keeping canonical target ordering",
    )
    return parser.parse_args()



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


def run_bcftools_norm(input_vcf: Path, output_vcf: Path, fasta_path: Path, *, log_path: Path) -> None:
    cmd = [resolve_bcftools_binary(), "norm", "-f", str(fasta_path), "-c", "e",
           "-o", str(output_vcf), "-Ov", str(input_vcf)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    write_bcftools_norm_log(log_path, cmd, result)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip()
        raise ValueError(f"bcftools norm failed: {detail or result.returncode} (see {log_path})")


def is_normalized_acgt_allele(value: str) -> bool:
    return bool(value) and set(value) <= {"A", "C", "G", "T"}


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
            if len(parts) < 8:
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
        ref, alt = ref.upper(), alt.upper()
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
        outcomes[candidate_id] = NormalizationOutcome(canonical, candidate.local_op)
    return outcomes


def require_reference_matches(rows: Sequence[VariantRow], fasta_path: Path, contig_naming: str) -> None:
    # Assumes: canonical sequence rows from assignment or normalized-output parsing.
    # Performs: SV(full REF invariant at the external-tool boundary).
    # Guarantees: every supplied a2 matches the configured FASTA or the run fails.
    matches = reference_alleles_match(
        fasta_path,
        [normalize_contig_for_reference(row.chrom, contig_naming, "ucsc") for row in rows],
        [int(row.pos) for row in rows], [row.a2 for row in rows],
    )
    if not matches.all():
        raise ValueError("unexpected REF mismatch after reference-allele assignment")


def normalize_candidates(
    candidates: Sequence[NormalizationCandidate], *, fasta_path: Path,
    contig_naming: str, output_path: Path,
) -> Dict[str, NormalizationOutcome]:
    # Assumes: each record carries one selected REF/ALT assignment and a stable ID.
    # Performs: SV(reference invariants), batched external normalization and boundary PV/SV.
    # Guarantees: one outcome per selected ID, retaining each batch's diagnostics.
    outcomes: Dict[str, NormalizationOutcome] = {}
    # PERF: bounded batches limit VCF/record bookkeeping; never spawn per variant.
    for start in range(0, len(candidates), NORMALIZATION_BATCH_ROWS):
        batch = candidates[start:start + NORMALIZATION_BATCH_ROWS]
        require_reference_matches([candidate.row for candidate in batch], fasta_path, contig_naming)
        batch_path = output_path.with_name(output_path.name + f".bcftools_norm_batch{start // NORMALIZATION_BATCH_ROWS}")
        input_vcf = batch_path.with_name(batch_path.name + ".input.vcf")
        output_vcf = batch_path.with_name(batch_path.name + ".output.vcf")
        log_path = batch_path.with_name(batch_path.name + ".log")
        write_normalization_vcf(input_vcf, batch, contig_naming)
        run_bcftools_norm(input_vcf, output_vcf, fasta_path, log_path=log_path)
        parsed = parse_normalized_candidates(output_vcf, {c.candidate_id: c for c in batch}, contig_naming)
        require_reference_matches([o.row for o in parsed.values() if o.row is not None], fasta_path, contig_naming)
        outcomes.update(parsed)
    return outcomes


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


def build_restriction_outcomes_table(
    source_rows_table: VariantRowsTable, *, fasta_path: Path, contig_naming: str,
    allow_strand_flips: bool, norm_indels: bool, genome_build: str,
    output_path: Path, include_indels: bool = False, allow_indel_swaps: bool = False,
) -> pd.DataFrame:
    # Assumes: loader-owned PN/PV and boundary allele/contig validation have completed.
    # Performs: SV(input filtering, ordered full-REF assignment, normalization and ploidy).
    # Guarantees: aligned outcomes; structural indexes never replace stored QC locators.
    source = source_rows_table.frame.loc[:, ["chrom", "pos", "id", "a1", "a2"]].reset_index(drop=True)
    outcomes = source.add_prefix("out_").copy()
    outcomes["keep"] = False
    outcomes["local_op"] = "identity"
    outcomes["status"] = None
    identical = source["a1"].eq(source["a2"])
    snv = source["a1"].str.len().eq(1) & source["a2"].str.len().eq(1)
    outcomes.loc[identical, "status"] = "identical_alleles"
    outcomes.loc[~identical & ~snv & (not include_indels), "status"] = "non_snv_input"
    snv_idx = source.index[snv & ~identical]
    if len(snv_idx):
        selected = restrict_rows_table(
            VariantRowsTable.from_frame(source.loc[snv_idx].reset_index(drop=True), copy=False),
            fasta_path, contig_naming, allow_strand_flips,
        )
        for column in ("keep", "local_op", "out_a1", "out_a2"):
            outcomes.loc[snv_idx, column] = selected[column].to_numpy(dtype=bool if column == "keep" else object)

    multibase = source.loc[~snv & ~identical] if include_indels else source.iloc[:0]
    contigs = [normalize_contig_for_reference(chrom, contig_naming, "ucsc") for chrom in multibase["chrom"]]
    positions = multibase["pos"].astype(int).tolist()
    identity = reference_alleles_match(fasta_path, contigs, positions, multibase["a2"].tolist())
    swap = reference_alleles_match(fasta_path, contigs, positions, multibase["a1"].tolist()) if allow_indel_swaps else identity & False
    swapped_idx = multibase.index[~identity & swap]
    selected_idx = multibase.index[identity | swap]
    outcomes.loc[selected_idx, "keep"] = True
    outcomes.loc[swapped_idx, "local_op"] = "swap"
    outcomes.loc[swapped_idx, ["out_a1", "out_a2"]] = source.loc[
        swapped_idx, ["a2", "a1"]
    ].to_numpy(dtype=object)
    selected_candidates: List[NormalizationCandidate] = []
    if norm_indels:
        # PERF: row objects are needed only for multibase VCF emission and round-trip bookkeeping.
        for idx, values in zip(selected_idx, outcomes.loc[selected_idx].itertuples(index=False, name=None)):
            chrom, pos, row_id, a1, a2, _keep, local_op, _status = values
            selected_candidates.append(NormalizationCandidate(
                f"row{idx}", idx, VariantRow(chrom, pos, row_id, a1, a2), local_op,
            ))
    if selected_candidates:
        normalized = normalize_candidates(selected_candidates, fasta_path=fasta_path, contig_naming=contig_naming, output_path=output_path)
        # PERF: normalization bookkeeping is multibase-only and bounded by candidate batches.
        for candidate in selected_candidates:
            idx = candidate.input_index
            result = apply_normalization_ploidy_filter(candidate.row, normalized[candidate.candidate_id], genome_build=genome_build)
            outcomes.at[idx, "keep"] = result.row is not None
            outcomes.at[idx, "status"] = result.status
            if result.row is not None:
                outcomes.loc[idx, ["out_chrom", "out_pos", "out_a1", "out_a2"]] = [result.row.chrom, result.row.pos, result.row.a1, result.row.a2]
    unmatched = ~outcomes["keep"] & outcomes["status"].isna()
    outcomes.loc[unmatched, "status"] = "no_reference_match"
    return outcomes


def main() -> int:
    # Assumes: source variant object path and metadata sidecar are valid and loadable.
    # Performs: orchestration of restriction/normalization SV kernels and output CV checks.
    # Guarantees: output object type matches input type with preserved metadata and QC sidecar behavior.
    args = parse_args()
    source_path = Path(args.source)
    output_path = Path(args.output)
    logger.info("restrict_build_compatible.py: restricting %s -> %s", source_path, output_path)
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
    output_path.parent.mkdir(parents=True, exist_ok=True)
    for diagnostic in output_path.parent.iterdir():
        if diagnostic.name.startswith(output_path.name + ".bcftools_norm") and diagnostic.is_file():
            diagnostic.unlink()
    if not args.include_indels and (args.allow_indel_swaps or args.norm_indels):
        logger.warning("multibase operations are no-ops without --include-indels because multibase input is excluded")

    restricted_outcomes = build_restriction_outcomes_table(
        source_rows_table,
        fasta_path=fasta_path,
        contig_naming=contig_naming,
        allow_strand_flips=args.allow_strand_flips,
        norm_indels=bool(args.norm_indels),
        include_indels=args.include_indels,
        allow_indel_swaps=args.allow_indel_swaps,
        genome_build=str(genome_build),
        output_path=output_path,
    )
    keep_mask = restricted_outcomes["keep"].astype(bool)
    flip_count = int((keep_mask & restricted_outcomes["local_op"].isin(["flip", "flip_swap"])).sum())

    # Both object types share target selection, duplicate handling, and QC.
    # Only .vmap carries allele-operation history into the output object.
    base_frame = (
        loaded.base_vmap_table.frame if loaded.base_vmap_table is not None else source_rows_table.frame
    )
    base_qc = variant_row_qc_rows(base_frame, is_vmap=loaded.base_vmap_table is not None)
    kept_idx = restricted_outcomes.index[keep_mask]
    out_frame = restricted_outcomes.loc[
        kept_idx, ["out_chrom", "out_pos", "out_id", "out_a1", "out_a2"]
    ].rename(columns={
        "out_chrom": "chrom", "out_pos": "pos", "out_id": "id",
        "out_a1": "a1", "out_a2": "a2",
    }).copy()
    if loaded.base_vmap_table is not None:
        if len(base_frame) != len(restricted_outcomes):
            raise ValueError("input variant_map source and target row counts are inconsistent")
        out_frame["source_shard"] = base_frame.loc[kept_idx, "source_shard"]
        out_frame["source_index"] = base_frame.loc[kept_idx, "source_index"]
        out_frame["allele_op"] = compose_allele_ops_series(
            base_frame.loc[kept_idx, "allele_op"],
            restricted_outcomes.loc[kept_idx, "local_op"],
        )
    else:
        out_frame["row_idx"] = source_rows_table.frame.loc[kept_idx, "row_idx"]

    qc_frames: List[pd.DataFrame] = []
    status_mask = (~keep_mask) & restricted_outcomes["status"].notna()
    if bool(status_mask.any()):
        dropped_idx = restricted_outcomes.index[status_mask]
        qc_status = base_qc.loc[dropped_idx].copy()
        qc_status["status"] = restricted_outcomes.loc[dropped_idx, "status"]
        qc_frames.append(qc_status)

    duplicate_mask = duplicate_target_rows_mask_table(out_frame)
    duplicate_target_count = int(duplicate_mask.sum())
    if bool(duplicate_mask.any()):
        qc_duplicate = base_qc.loc[out_frame.index[duplicate_mask]].copy()
        qc_duplicate["status"] = "duplicate_target"
        qc_frames.append(qc_duplicate)
        out_frame = out_frame.loc[~duplicate_mask]

    if qc_frames:
        qc_all = pd.concat(qc_frames, axis=0, ignore_index=True)
        write_variant_status_qc(qc_output_path, qc_all, is_vmap=loaded.base_vmap_table is not None)
    elif qc_output_path.exists():
        qc_output_path.unlink()

    out_frame = out_frame.reset_index(drop=True)
    if loaded.base_vmap_table is not None:
        write_vmap_table(output_path, VMapRowsTable.from_frame(out_frame, copy=False), assume_validated=True)
    else:
        write_vtable_table(output_path, VariantRowsTable.from_frame(out_frame, copy=False), assume_validated=True)
    output_row_count = len(out_frame)
    write_metadata(output_path, metadata_with_variants_count(loaded.raw_metadata, output_row_count))
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
        "include_indels": args.include_indels,
        "allow_indel_swaps": args.allow_indel_swaps,
        "assignment_counts": restricted_outcomes.loc[
            ~restricted_outcomes["status"].isin(["identical_alleles", "non_snv_input", "no_reference_match"]),
            "local_op",
        ].value_counts().to_dict(),
        "exclusion_counts": {**restricted_outcomes["status"].dropna().value_counts().to_dict(), "duplicate_target": duplicate_target_count},
        "input_rows": input_row_count,
        "output_rows": output_row_count,
        "dropped_rows": input_row_count - output_row_count,
        "duplicate_target_rows": duplicate_target_count,
        "flip_rows": flip_count,
        "reference_fasta": str(fasta_path),
    }
    logger.info(
        "restrict_build_compatible.py: completed with %s input rows, %s output rows, %s dropped rows",
        input_row_count,
        output_row_count,
        input_row_count - output_row_count,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli(main))
