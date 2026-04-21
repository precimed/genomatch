#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from ._cli_utils import run_cli
from .contig_utils import UNKNOWN_CONTIG, canonical_contig_from_label
from .reference_utils import fetch_reference_bases, resolve_bcftools_binary, resolve_liftover_assets
from .haploid_utils import expected_ploidy_pair
from .tabular_rows import VMapRowsTable, VariantRowsTable
from .vectorization_utils import (
    first_true_index,
    map_unique_values,
    require_columns,
    strict_int_series,
)
from .vtable_utils import (
    CANONICAL_CONTIG_RANK,
    compose_allele_ops_series,
    convert_contig_label,
    load_variant_object_tables,
    normalize_chrom_label,
    require_contig_naming,
    require_table_matches_contig_naming,
    write_metadata,
    write_vmap_status_qc,
    write_vmap_table,
    write_vtable_table,
)

PREPARE_INPUT_COLUMNS = (
    "row_id",
    "chrom",
    "pos",
    "id",
    "a1",
    "a2",
    "source_shard",
    "source_index",
    "upstream_allele_op",
)
PREPARED_COLUMNS = (
    "row_id",
    "chrom",
    "pos",
    "id",
    "a1",
    "a2",
    "source_shard",
    "source_index",
    "upstream_allele_op",
    "status",
    "input_allele_op",
    "ref",
    "alt",
)
ROW_LOOKUP_COLUMNS = (
    "row_id",
    "chrom",
    "pos",
    "id",
    "source_shard",
    "source_index",
    "upstream_allele_op",
    "input_allele_op",
)

WRITE_CHUNK_ROWS = 100_000
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Liftover a .vtable/.vmap across genome builds, preserving the input object type. "
            "Input contigs may be UCSC or numeric/NCBI-style, but the internal reference path is UCSC-only."
        )
    )
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable or .vmap")
    parser.add_argument("--target-build", required=True, help="Target genome build")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse retained bcftools +liftover output if it already exists, then rerun parsing/QC/output steps",
    )
    return parser.parse_args()


def convert_rows_to_ucsc(rows: VariantRowsTable, contig_naming: str) -> VariantRowsTable:
    frame = rows.to_frame(copy=True)
    def to_ucsc(chrom_str: str) -> str:
        try:
            return convert_contig_label(chrom_str, contig_naming, "ucsc")
        except ValueError as exc:
            raise ValueError(
                f"unable to normalize input contig {chrom_str!r} from declared contig_naming={contig_naming!r} "
                "to internal UCSC naming for liftover"
            ) from exc
    frame["chrom"] = map_unique_values(frame["chrom"], to_ucsc)
    return VariantRowsTable.from_frame(frame, copy=False)


def prepare_liftover_table(frame: pd.DataFrame, source_fasta: Path) -> pd.DataFrame:
    # Assumes: liftover_input columns exist and upstream loaders already canonicalized alleles.
    # Performs: SV(position/allele/reference-anchoring checks), query construction for source reference bases.
    # Guarantees: PREPARED_COLUMNS with status/ref/alt fields ready for VCF emission.
    prepared = frame.copy(deep=True)
    require_columns(prepared, PREPARE_INPUT_COLUMNS, label="liftover_input")

    pos_int, valid_pos = strict_int_series(prepared["pos"])
    if not bool(valid_pos.all()):
        first_bad_idx = first_true_index(~valid_pos)
        raise ValueError(
            f"invalid source position for liftover input: "
            f"{prepared.at[first_bad_idx, 'chrom']}:{prepared.at[first_bad_idx, 'pos']}"
        )
    prepared["pos_int"] = pos_int.astype("int64")
    a1_tokens = prepared["a1"]
    a2_tokens = prepared["a2"]
    invalid_a1 = ~a1_tokens.str.fullmatch(r"[ACGT]+")
    invalid_a2 = ~a2_tokens.str.fullmatch(r"[ACGT]+")
    any_invalid = invalid_a1 | invalid_a2
    if any_invalid.any():
        first_idx = int(any_invalid.idxmax())
        bad_value = prepared.at[first_idx, "a1"] if bool(invalid_a1.iloc[first_idx]) else prepared.at[first_idx, "a2"]
        raise ValueError(f"invalid allele code in liftover_build.py input: {bad_value!r}")

    query_frame = prepared.loc[:, ["chrom", "pos_int"]].drop_duplicates(ignore_index=True)
    query_pairs = list(query_frame.itertuples(index=False, name=None))
    reference_bases = fetch_reference_bases(source_fasta, query_pairs)
    ref_frame = pd.DataFrame(
        [{"chrom": chrom, "pos_int": pos, "ref_base": ref_base} for (chrom, pos), ref_base in reference_bases.items()]
    )
    if ref_frame.empty:
        prepared["ref_base"] = ""
    else:
        prepared = prepared.merge(ref_frame, on=["chrom", "pos_int"], how="left")
        prepared["ref_base"] = prepared["ref_base"].fillna("")
    valid_ref_mask = prepared["ref_base"].isin({"A", "C", "G", "T"})
    if not bool(valid_ref_mask.all()):
        first_bad_idx = first_true_index(~valid_ref_mask)
        raise ValueError(
            f"source reference base not found for liftover input: {prepared.at[first_bad_idx, 'chrom']}:{prepared.at[first_bad_idx, 'pos']}"
        )

    is_snv = (a1_tokens.str.len() == 1) & (a2_tokens.str.len() == 1)
    ready_swap = (a2_tokens == prepared["ref_base"]) & (a1_tokens != prepared["ref_base"])
    ready_identity = (a1_tokens == prepared["ref_base"]) & (a2_tokens != prepared["ref_base"])
    ready_mask = ready_swap | ready_identity
    bad_snv_mask = (~ready_mask) & is_snv
    if bool(bad_snv_mask.any()):
        first_bad_idx = first_true_index(bad_snv_mask)
        raise ValueError(
            f"source alleles do not match reference for liftover input: {prepared.at[first_bad_idx, 'chrom']}:{prepared.at[first_bad_idx, 'pos']} "
            f"({prepared.at[first_bad_idx, 'a1']}/{prepared.at[first_bad_idx, 'a2']}, ref={prepared.at[first_bad_idx, 'ref_base']})"
        )

    prepared["status"] = ready_mask.map({True: "ready", False: "unsupported_non_snv"})
    prepared["input_allele_op"] = ready_swap.map({True: "swap", False: "identity"})
    prepared["ref"] = prepared["ref_base"].where(ready_mask, None)
    prepared["alt"] = a2_tokens.where(~ready_swap, a1_tokens).where(ready_mask, None)
    return prepared.loc[:, list(PREPARED_COLUMNS)]


def write_temp_vcf(
    path: Path,
    prepared: pd.DataFrame,
) -> pd.DataFrame:
    # Assumes: prepared follows PREPARED_COLUMNS contract from prepare_liftover_table.
    # Performs: CV(VCF row-shape construction for ready rows), VCF text emission.
    # Guarantees: VCF file for bcftools liftover and a row_id-indexed lookup table.
    require_columns(prepared, PREPARED_COLUMNS, label="prepared_liftover")
    ready_frame = prepared.loc[prepared["status"] == "ready"].copy()
    row_lookup = ready_frame.loc[:, list(ROW_LOOKUP_COLUMNS)].set_index("row_id", drop=False)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom in sorted(set(ready_frame["chrom"].tolist()), key=chrom_sort_key):
            handle.write(f"##contig=<ID={chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        lines = (
            ready_frame["chrom"].astype(str) + "\t"
            + ready_frame["pos"].astype(str) + "\t"
            + ready_frame["row_id"].astype(str) + "\t"
            + ready_frame["ref"].fillna("").astype(str) + "\t"
            + ready_frame["alt"].fillna("").astype(str)
            + "\t.\tPASS\t.\n"
        )
        # PERF: loop retained for bounded-memory emission; full-string materialization scales poorly.
        for start in range(0, len(lines), WRITE_CHUNK_ROWS):
            stop = start + WRITE_CHUNK_ROWS
            handle.write("".join(lines.iloc[start:stop].tolist()))
    return row_lookup


def run_bcftools_liftover(
    input_vcf: Path,
    output_vcf: Path,
    source_fasta: Path,
    target_fasta: Path,
    chain_path: Path,
) -> None:
    bcftools = resolve_bcftools_binary()
    cmd = [
        bcftools,
        "+liftover",
        str(input_vcf),
        "-o",
        str(output_vcf),
        "--",
        "-s",
        str(source_fasta),
        "-f",
        str(target_fasta),
        "-c",
        str(chain_path),
        "--flip-tag",
        "FLIP",
        "--swap-tag",
        "SWAP",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        text = (result.stderr or result.stdout).strip()
        raise ValueError(f"bcftools +liftover failed: {text}")
def chrom_sort_key(label: str) -> Tuple[int, str]:
    canonical = normalize_chrom_label(label)
    if canonical.isdigit():
        return (int(canonical), "")
    order = {"X": 23, "Y": 24, "MT": 25}
    if canonical in order:
        return (order[canonical], "")
    return (999, label)


def parse_lifted_vcf(
    path: Path,
    row_lookup: pd.DataFrame,
    final_contig_naming: str,
    source_build: str,
    target_build: str,
) -> Tuple[VMapRowsTable, Dict[str, str]]:
    # Assumes: lifted VCF is produced by bcftools +liftover and row_lookup is keyed by row_id.
    # Performs: PV(lifted VCF row-shape checks), SV(liftover-specific filtering, contig/ploidy checks, allele-op composition).
    # Guarantees: deduplicated VMapRowsTable of lifted rows plus per-row QC status map.
    status_by_row_id: Dict[str, str] = {}

    try:
        frame = pd.read_table(
            path,
            sep="\t",
            header=None,
            names=["chrom", "pos", "row_id", "ref", "alt", "qual", "filter", "info"],
            usecols=range(8),
            comment="#",
            dtype="object",
            keep_default_na=False,
            na_values=[],
            skip_blank_lines=True,
        )
    except Exception as exc:
        raise ValueError(f"invalid lifted VCF content: {path}") from exc

    if frame.empty:
        return VMapRowsTable.from_rows([]), status_by_row_id

    if frame.isna().any(axis=None):
        raise ValueError(f"invalid lifted VCF content: {path}")

    unknown_row_mask = ~frame["row_id"].isin(row_lookup.index)
    if bool(unknown_row_mask.any()):
        first_bad_idx = first_true_index(unknown_row_mask)
        raise ValueError(f"unexpected liftover row ID: {frame.at[first_bad_idx, 'row_id']}")

    frame["pass_filter"] = frame["filter"].isin({".", "PASS"})
    frame["ref_token"] = frame["ref"].str.strip().str.upper()
    frame["alt_token"] = frame["alt"].str.strip().str.upper()
    frame["status"] = None

    active_mask = frame["pass_filter"]
    multiallelic_mask = active_mask & frame["alt"].str.contains(",", regex=False)
    frame.loc[multiallelic_mask, "status"] = "unsupported_non_snv"

    token_valid_mask = (
        frame["ref_token"].ne("")
        & frame["alt_token"].ne("")
        & frame["ref_token"].str.fullmatch(r"[ACGT]+")
        & frame["alt_token"].str.fullmatch(r"[ACGT]+")
    )
    invalid_token_mask = active_mask & frame["status"].isna() & ~token_valid_mask
    frame.loc[invalid_token_mask, "status"] = "unsupported_non_snv"

    both_multibase_mask = (
        active_mask
        & frame["status"].isna()
        & (frame["ref_token"].str.len() > 1)
        & (frame["alt_token"].str.len() > 1)
    )
    frame.loc[both_multibase_mask, "status"] = "unsupported_non_snv"

    frame["final_chrom"] = None
    contig_candidates = frame.loc[active_mask & frame["status"].isna(), "chrom"]
    contig_map: Dict[str, str] = {}
    invalid_contigs: set[str] = set()
    # PERF: retained loop is over unique contigs (small cardinality), not per-row.
    for chrom in pd.unique(contig_candidates):
        chrom_str = str(chrom)
        try:
            contig_map[chrom_str] = convert_contig_label(chrom_str, "ucsc", final_contig_naming)
        except ValueError:
            invalid_contigs.add(chrom_str)
    contig_candidate_mask = active_mask & frame["status"].isna()
    invalid_contig_mask = contig_candidate_mask & frame["chrom"].isin(invalid_contigs)
    frame.loc[invalid_contig_mask, "status"] = "unsupported_target_contig"
    valid_contig_mask = contig_candidate_mask & ~invalid_contig_mask
    frame.loc[valid_contig_mask, "final_chrom"] = frame.loc[valid_contig_mask, "chrom"].map(contig_map)

    row_lookup_merge = row_lookup.reset_index(drop=True)
    frame = frame.merge(
        row_lookup_merge.loc[
            :, ["row_id", "id", "source_shard", "source_index", "upstream_allele_op", "input_allele_op"]
        ].rename(
            columns={
                "id": "lookup_id",
                "source_shard": "lookup_source_shard",
                "source_index": "lookup_source_index",
                "upstream_allele_op": "lookup_upstream_allele_op",
                "input_allele_op": "lookup_input_allele_op",
            }
        ),
        on="row_id",
        how="left",
    )

    ploidy_candidate_mask = active_mask & frame["status"].isna()
    canonical_target = frame["final_chrom"].map(lambda chrom: normalize_chrom_label(str(chrom)))
    non_autosome_mask = ~canonical_target.isin({str(i) for i in range(1, 23)})
    ploidy_check_mask = ploidy_candidate_mask & non_autosome_mask
    if bool(ploidy_check_mask.any()):
        ploidy_frame = frame.loc[ploidy_check_mask, ["row_id", "final_chrom", "pos"]].copy()
        source_lookup = row_lookup.loc[ploidy_frame["row_id"], ["chrom", "pos"]].reset_index(drop=True)
        ploidy_frame["source_key"] = list(zip(source_lookup["chrom"], source_lookup["pos"]))
        ploidy_frame["target_key"] = list(
            zip(ploidy_frame["final_chrom"], ploidy_frame["pos"])
        )

        source_ploidy_by_key = {
            key: expected_ploidy_pair(key[0], key[1], genome_build=source_build)
            for key in pd.unique(ploidy_frame["source_key"])
        }
        target_ploidy_by_key = {
            key: expected_ploidy_pair(key[0], key[1], genome_build=target_build)
            for key in pd.unique(ploidy_frame["target_key"])
        }

        source_pairs = ploidy_frame["source_key"].map(source_ploidy_by_key)
        target_pairs = ploidy_frame["target_key"].map(target_ploidy_by_key)
        mismatch_mask = target_pairs != source_pairs
        if bool(mismatch_mask.any()):
            mismatch_indices = ploidy_frame.index[mismatch_mask]
            frame.loc[mismatch_indices, "status"] = "ploidy_class_changed"

    success_mask = active_mask & frame["status"].isna()
    success_frame = frame.loc[success_mask].copy()
    success_frame["allele_op"] = pd.Series(index=success_frame.index, dtype="object")
    if not success_frame.empty:
        success_frame["liftover_op"] = parse_liftover_info_ops(success_frame["info"])
        success_frame["op_after_output"] = compose_allele_ops_series(
            success_frame["liftover_op"],
            pd.Series(["swap"] * len(success_frame), index=success_frame.index, dtype="object"),
        )
        success_frame["op_after_input"] = compose_allele_ops_series(
            success_frame["lookup_input_allele_op"],
            success_frame["op_after_output"],
        )
        success_frame["allele_op"] = compose_allele_ops_series(
            success_frame["lookup_upstream_allele_op"],
            success_frame["op_after_input"],
        )

    rejected_rows = frame.loc[frame["status"].notna(), ["row_id", "status"]]
    status_by_row_id.update(dict(zip(rejected_rows["row_id"], rejected_rows["status"])))
    status_by_row_id.update(dict(zip(success_frame["row_id"], ["lifted"] * len(success_frame))))

    parsed_frame = (
        success_frame.loc[
            :,
            [
                "row_id",
                "final_chrom",
                "pos",
                "lookup_id",
                "alt_token",
                "ref_token",
                "lookup_source_shard",
                "lookup_source_index",
                "allele_op",
            ],
        ]
        .rename(
            columns={
                "final_chrom": "chrom",
                "lookup_id": "id",
                "alt_token": "a1",
                "ref_token": "a2",
                "lookup_source_shard": "source_shard",
                "lookup_source_index": "source_index",
            }
        )
        .copy()
    )
    if parsed_frame.empty:
        return VMapRowsTable.from_rows([]), status_by_row_id

    canonical = parsed_frame["chrom"].map(lambda chrom: canonical_contig_from_label(str(chrom), final_contig_naming))
    invalid_contig = canonical.isna()
    if bool(invalid_contig.any()):
        first_bad_idx = first_true_index(invalid_contig)
        bad_chrom = parsed_frame.at[first_bad_idx, "chrom"]
        raise ValueError(
            f"liftover_build.py output has contigs inconsistent with declared contig_naming={final_contig_naming!r}; "
            f"first invalid label: {bad_chrom!r}. Run normalize_contigs.py first"
        )
    unknown_contig = canonical.eq(UNKNOWN_CONTIG)
    if bool(unknown_contig.any()):
        raise ValueError("liftover_build.py output contains chrom=unknown; run normalize_contigs.py first")
    parsed_frame["contig_rank"] = canonical.map(CANONICAL_CONTIG_RANK)
    pos_int, valid_pos = strict_int_series(parsed_frame["pos"])
    invalid_pos = (~valid_pos) | (pos_int <= 0)
    if bool(invalid_pos.any()):
        first_bad_idx = first_true_index(invalid_pos)
        raise ValueError(f"liftover_build.py output row has invalid pos: {parsed_frame.at[first_bad_idx, 'pos']!r}")
    parsed_frame["pos_int"] = pos_int.astype("int64")
    parsed_frame = parsed_frame.sort_values(["contig_rank", "pos_int"], kind="stable")

    duplicate_mask = parsed_frame.duplicated(subset=["chrom", "pos", "a1", "a2"], keep=False)
    status_by_row_id.update(dict(zip(parsed_frame.loc[duplicate_mask, "row_id"], ["duplicate_target"] * int(duplicate_mask.sum()))))

    out_frame = parsed_frame.loc[
        ~duplicate_mask, ["chrom", "pos", "id", "a1", "a2", "source_shard", "source_index", "allele_op"]
    ]
    return VMapRowsTable.from_frame(out_frame.reset_index(drop=True), copy=False), status_by_row_id


def parse_liftover_info_ops(info: pd.Series) -> pd.Series:
    # Assumes: `info` holds VCF INFO tokens from bcftools +liftover.
    # Performs: vectorized FLIP/SWAP token parsing with strict SWAP-value validation.
    # Guarantees: allele-op series aligned to input index with canonical v1 operation tokens.
    info_tokens = info.fillna("").astype(str)
    has_flip = info_tokens.str.contains(r"(?:^|;)FLIP(?:;|$)", regex=True)
    swap_value = info_tokens.str.extract(r"(?:^|;)SWAP=([^;]*)(?:;|$)", expand=False)
    invalid_swap = swap_value.notna() & swap_value.ne("1")
    if bool(invalid_swap.any()):
        first_bad_idx = first_true_index(invalid_swap)
        bad_value = swap_value.at[first_bad_idx]
        raise ValueError(
            f"liftover output contains unsupported SWAP annotation {bad_value!r}; "
            "cannot represent this variant in canonical v1 allele_op"
        )
    has_swap = swap_value.eq("1").fillna(False)

    out = pd.Series("identity", index=info_tokens.index, dtype="object")
    out.loc[has_swap] = "swap"
    out.loc[has_flip] = "flip"
    out.loc[has_flip & has_swap] = "flip_swap"
    return out


def liftover_debug_vcf_path(output_path: Path, *, kind: str) -> Path:
    return output_path.with_name(output_path.name + f".bcftools_{kind}.vcf")


def print_resume_skip(output_vcf: Path) -> None:
    logger.info("liftover_build.py: skipping bcftools +liftover; retained output exists at %s", output_vcf)


def main() -> int:
    # Assumes: input object and metadata sidecar are present and valid for loading.
    # Performs: orchestration of PN/PV/SV/CV stages across load, liftover prep, parse, and output writing.
    # Guarantees: output object type matches input type with updated build metadata and QC sidecar.
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    logger.info("liftover_build.py: lifting %s -> %s (target_build=%s)", input_path, output_path, args.target_build)
    if not input_path.exists():
        raise ValueError(f"input not found: {input_path}")

    loaded = load_variant_object_tables(input_path)
    source_rows_table = loaded.target_rows_table
    source_meta = dict(loaded.target_metadata)
    source_build = str(source_meta.get("genome_build"))
    source_contig_naming = require_contig_naming(source_meta, label="variant object")
    if source_contig_naming == "plink_splitx":
        raise ValueError(
            "liftover_build.py does not support contig_naming=plink_splitx; "
            "normalize to build-independent ncbi/ucsc/plink first, liftover, then normalize back if needed"
        )
    target_build = str(args.target_build)
    if source_build == "unknown":
        raise ValueError(
            "source metadata genome_build must be known for liftover "
            f"(declared genome_build={source_build!r}, declared contig_naming={source_contig_naming!r}, "
            "internal_reference_naming='ucsc')"
        )
    if source_build == target_build:
        raise ValueError("source and target builds are identical; liftover is not required")
    require_table_matches_contig_naming(source_rows_table, source_contig_naming, label="variant object")

    source_fasta, target_fasta, chain_path = resolve_liftover_assets(source_build, target_build)
    ucsc_rows = convert_rows_to_ucsc(source_rows_table, source_contig_naming)
    input_frame = ucsc_rows.to_frame(copy=True)
    input_frame["row_id"] = [f"row{idx}" for idx in range(len(input_frame))]
    if loaded.base_vmap_table is not None:
        base_frame = loaded.base_vmap_table.to_frame(copy=False)
        input_frame["source_shard"] = base_frame["source_shard"].tolist()
        input_frame["source_index"] = base_frame["source_index"].tolist()
        input_frame["upstream_allele_op"] = base_frame["allele_op"].tolist()
    else:
        input_frame["source_shard"] = "."
        input_frame["source_index"] = list(range(len(input_frame)))
        input_frame["upstream_allele_op"] = "identity"

    input_vcf = liftover_debug_vcf_path(output_path, kind="input")
    output_vcf = liftover_debug_vcf_path(output_path, kind="output")
    input_vcf.parent.mkdir(parents=True, exist_ok=True)
    prepared = prepare_liftover_table(input_frame, source_fasta)
    row_lookup = write_temp_vcf(input_vcf, prepared)
    if args.resume and output_vcf.exists():
        print_resume_skip(output_vcf)
    else:
        run_bcftools_liftover(input_vcf, output_vcf, source_fasta, target_fasta, chain_path)
    raw_rows_table, parse_status_by_row_id = parse_lifted_vcf(
        output_vcf,
        row_lookup,
        source_contig_naming,
        source_build,
        target_build,
    )

    prepared["qc_status"] = prepared["status"].astype(str)
    ready_mask = prepared["qc_status"] == "ready"
    prepared.loc[ready_mask, "qc_status"] = prepared.loc[ready_mask, "row_id"].map(parse_status_by_row_id).fillna("unmapped")
    qc_rows = list(
        prepared.loc[:, ["source_shard", "source_index", "id", "qc_status"]]
        .assign(source_index=lambda df: df["source_index"].astype(int))
        .itertuples(index=False, name=None)
    )

    if loaded.base_vmap_table is not None:
        out_meta = dict(loaded.raw_metadata)
        out_meta["target"] = dict(out_meta["target"])
        out_meta["target"]["genome_build"] = target_build
        write_vmap_table(output_path, raw_rows_table, assume_validated=True)
    else:
        out_meta = dict(loaded.raw_metadata)
        out_meta["genome_build"] = target_build
        variant_table = VariantRowsTable.from_frame(
            raw_rows_table.to_frame(copy=False).loc[:, ["chrom", "pos", "id", "a1", "a2"]],
            copy=False,
        )
        write_vtable_table(output_path, variant_table, assume_validated=True)
    write_metadata(output_path, out_meta)
    write_vmap_status_qc(output_path.with_name(output_path.name + ".qc.tsv"), qc_rows)
    logger.info(
        "liftover_build.py: wrote %s with %s output rows and %s QC rows",
        output_path,
        len(raw_rows_table),
        len(qc_rows),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli(main))
