#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd

from ._cli_utils import run_cli
from .exact_set_utils import load_target_object_info
from .importer_utils import (
    finalize_imported_vmap_vectorized,
    is_canonical_allele_token,
    is_valid_import_position,
    reject_template_argument,
)
from .contig_utils import canonical_contig_from_any_supported_label
from .sumstats_utils import (
    extract_variant_field,
    find_metadata_value,
    iter_sumstats_table_chunks,
    is_missing_token_series,
    load_metadata,
    read_sumstats_header,
    resolve_sumstats_input_path,
    resolve_column,
    resolve_variant_columns,
    VariantColumnMapping,
)
from .vtable_utils import (
    normalize_allele_token,
    open_text,
    parse_chr2use,
    read_vmap_table,
    VariantRow,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract a .vmap from summary statistics.")
    parser.add_argument("--input", help="Input summary statistics file (optional when metadata defines path_sumStats)")
    parser.add_argument("--output", required=True, help="Output .vmap file")
    parser.add_argument("--sumstats-metadata", required=True, help="Cleansumstats-style metadata YAML")
    parser.add_argument("--genome-build", default="unknown", help="Genome build for metadata")
    parser.add_argument("--id-lookup", help="Optional .vtable or .vmap for ID-based coordinate enrichment")
    parser.add_argument("--chr2use", "--contigs", dest="chr2use", help="Comma-separated chromosome list or ranges")
    parser.add_argument("--max-allele-length", type=int, default=150, help="Maximum allele length; rows exceeding this are dropped (default: 150)")
    return parser.parse_args()


def load_id_lookup_object(path: Path) -> tuple[Dict[str, VariantRow], Set[str], Dict[str, object]]:
    if not (path.name.endswith(".vtable") or path.name.endswith(".vmap")):
        raise ValueError("--id-lookup must point to a .vtable or .vmap")
    info = load_target_object_info(path)
    if info.object_type == "variant_map":
        vmap_frame = read_vmap_table(path, check_duplicates=False).to_frame(copy=False)
        rows = [
            VariantRow(chrom, pos, row_id, a1, a2)
            for chrom, pos, row_id, a1, a2 in vmap_frame.loc[:, ["chrom", "pos", "id", "a1", "a2"]].itertuples(
                index=False,
                name=None,
            )
        ]
    else:
        rows = []
        with open_text(path, "rt") as handle:
            for line in handle:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 5:
                    raise ValueError(f"invalid vtable row in {path}: {line.strip()}")
                chrom, pos, row_id, a1, a2 = parts
                a1 = normalize_allele_token(a1)
                a2 = normalize_allele_token(a2)
                if not chrom or not is_valid_import_position(pos):
                    raise ValueError(f"invalid vtable row in {path}: {line.strip()}")
                if not is_canonical_allele_token(a1) or not is_canonical_allele_token(a2):
                    raise ValueError(f"invalid vtable row in {path}: {line.strip()}")
                rows.append(VariantRow(chrom, pos, row_id, a1, a2))
    unique_matches: Dict[str, VariantRow] = {}
    ambiguous_ids: Set[str] = set()
    ignored_ids = 0
    for row in rows:
        lookup_id = row.id.strip()
        if not lookup_id or lookup_id == ".":
            ignored_ids += 1
            continue
        if lookup_id in ambiguous_ids:
            continue
        if lookup_id in unique_matches:
            unique_matches.pop(lookup_id, None)
            ambiguous_ids.add(lookup_id)
            continue
        unique_matches[lookup_id] = row
    if ignored_ids:
        logger.warning("ignored %s --id-lookup rows whose id is missing, empty, or '.'", ignored_ids)
    return unique_matches, ambiguous_ids, {
        "genome_build": info.target_metadata["genome_build"],
        "contig_naming": info.target_metadata.get("contig_naming"),
    }


def resolve_id_enrichment_columns(header: List[str], metadata: Dict[str, object]) -> tuple[int, int, int]:
    if find_metadata_value(metadata, "col_CHR") is not None or find_metadata_value(metadata, "col_POS") is not None:
        raise ValueError("--id-lookup requires metadata to omit col_CHR and col_POS")
    snp_idx = resolve_column(header, find_metadata_value(metadata, "col_SNP"), "col_SNP", required=True)
    effect_idx = resolve_column(
        header,
        find_metadata_value(metadata, "col_EffectAllele"),
        "col_EffectAllele",
        required=True,
    )
    other_idx = resolve_column(
        header,
        find_metadata_value(metadata, "col_OtherAllele"),
        "col_OtherAllele",
        required=True,
    )
    return snp_idx, effect_idx, other_idx


def apply_common_allele_validation(
    reason: pd.Series,
    a1_series: pd.Series,
    a2_series: pd.Series,
    *,
    max_allele_length: int,
) -> None:
    missing_required_mask = reason.isna() & (
        is_missing_token_series(a1_series)
        | is_missing_token_series(a2_series)
    )
    reason.loc[missing_required_mask] = "malformed_row"

    non_actg_mask = reason.isna() & (
        ~a1_series.map(is_canonical_allele_token).fillna(False)
        | ~a2_series.map(is_canonical_allele_token).fillna(False)
    )
    reason.loc[non_actg_mask] = "non_actg_allele"

    allele_too_long_mask = reason.isna() & (
        (a1_series.map(len).fillna(0) > max_allele_length)
        | (a2_series.map(len).fillna(0) > max_allele_length)
    )
    reason.loc[allele_too_long_mask] = "allele_too_long"


def qc_frame_from_reason(reason: pd.Series, source_index_series: pd.Series) -> pd.DataFrame:
    qc_positions = reason[reason.notna()].index
    return pd.DataFrame(
        {
            "source_shard": ".",
            "source_index": source_index_series.iloc[qc_positions].values,
            "reason": reason.iloc[qc_positions].astype(str).values,
        },
        columns=["source_shard", "source_index", "reason"],
    )


def import_sumstats_variant_chunk(
    frame: pd.DataFrame,
    source_index_series: pd.Series,
    variant_columns: VariantColumnMapping,
    *,
    max_allele_length: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    n_rows = len(frame)
    reason = pd.Series([None] * n_rows, dtype="object")

    valid_mask = reason.isna()
    chrom_series = pd.Series([""] * n_rows, dtype="object")
    pos_series = pd.Series([""] * n_rows, dtype="object")
    row_id_series = pd.Series(["."] * n_rows, dtype="object")
    a1_series = pd.Series([""] * n_rows, dtype="object")
    a2_series = pd.Series([""] * n_rows, dtype="object")

    if bool(valid_mask.any()):
        active_idx = valid_mask[valid_mask].index
        chrom_raw = frame.iloc[active_idx, variant_columns.chr].astype(str)
        pos_raw = frame.iloc[active_idx, variant_columns.pos].astype(str)
        a1_raw = frame.iloc[active_idx, variant_columns.effect_allele].astype(str)
        a2_raw = frame.iloc[active_idx, variant_columns.other_allele].astype(str)

        chrom_series.iloc[active_idx] = chrom_raw.map(lambda value: extract_variant_field(value, "CHR"))
        pos_series.iloc[active_idx] = pos_raw.map(lambda value: extract_variant_field(value, "POS"))
        a1_series.iloc[active_idx] = a1_raw.map(lambda value: normalize_allele_token(extract_variant_field(value, "EffectAllele")))
        a2_series.iloc[active_idx] = a2_raw.map(lambda value: normalize_allele_token(extract_variant_field(value, "OtherAllele")))

        missing_position_mask = valid_mask & (
            is_missing_token_series(chrom_series)
            | is_missing_token_series(pos_series)
        )
        reason.loc[missing_position_mask] = "malformed_row"

        if variant_columns.snp is not None:
            snp_raw = frame.iloc[active_idx, variant_columns.snp].astype(str)
            row_id_series.iloc[active_idx] = snp_raw.where(snp_raw != "", ".")

        apply_common_allele_validation(
            reason,
            a1_series,
            a2_series,
            max_allele_length=max_allele_length,
        )

        malformed_mask = reason.isna() & (
            chrom_series.astype(str).str.strip().eq("")
            | ~pos_series.map(is_valid_import_position).fillna(False)
        )
        reason.loc[malformed_mask] = "malformed_row"

    retained_mask = reason.isna()
    retained_positions = retained_mask[retained_mask].index
    rows_frame = pd.DataFrame(
        {
            "chrom": chrom_series.iloc[retained_positions].astype(str).values,
            "pos": pos_series.iloc[retained_positions].astype(str).values,
            "id": row_id_series.iloc[retained_positions].astype(str).values,
            "a1": a1_series.iloc[retained_positions].astype(str).values,
            "a2": a2_series.iloc[retained_positions].astype(str).values,
            "source_shard": ".",
            "source_index": source_index_series.iloc[retained_positions].values,
        },
        columns=["chrom", "pos", "id", "a1", "a2", "source_shard", "source_index"],
    )
    qc_rows_frame = qc_frame_from_reason(reason, source_index_series)
    return rows_frame, qc_rows_frame


def import_sumstats_id_enrichment_chunk(
    frame: pd.DataFrame,
    source_index_series: pd.Series,
    *,
    snp_idx: int,
    effect_idx: int,
    other_idx: int,
    id_lookup_rows: Dict[str, VariantRow],
    ambiguous_lookup_ids: Set[str],
    max_allele_length: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    n_rows = len(frame)
    reason = pd.Series([None] * n_rows, dtype="object")

    valid_mask = reason.isna()
    raw_id_series = pd.Series([""] * n_rows, dtype="object")
    a1_series = pd.Series([""] * n_rows, dtype="object")
    a2_series = pd.Series([""] * n_rows, dtype="object")

    if bool(valid_mask.any()):
        active_idx = valid_mask[valid_mask].index
        raw_id_series.iloc[active_idx] = frame.iloc[active_idx, snp_idx].astype(str).str.strip()
        a1_raw = frame.iloc[active_idx, effect_idx].astype(str)
        a2_raw = frame.iloc[active_idx, other_idx].astype(str)
        a1_series.iloc[active_idx] = a1_raw.map(lambda value: normalize_allele_token(extract_variant_field(value, "EffectAllele")))
        a2_series.iloc[active_idx] = a2_raw.map(lambda value: normalize_allele_token(extract_variant_field(value, "OtherAllele")))

        apply_common_allele_validation(
            reason,
            a1_series,
            a2_series,
            max_allele_length=max_allele_length,
        )

        invalid_id_mask = reason.isna() & (
            raw_id_series.astype(str).str.strip().eq("")
            | raw_id_series.astype(str).str.strip().eq(".")
        )
        reason.loc[invalid_id_mask] = "invalid_id"

        ambiguous_id_mask = reason.isna() & raw_id_series.isin(ambiguous_lookup_ids)
        reason.loc[ambiguous_id_mask] = "ambiguous_id_match"

        missing_id_mask = reason.isna() & ~raw_id_series.isin(set(id_lookup_rows))
        reason.loc[missing_id_mask] = "id_not_found"

    retained_mask = reason.isna()
    retained_positions = retained_mask[retained_mask].index
    raw_id_retained = raw_id_series.iloc[retained_positions].astype(str)
    rows_frame = pd.DataFrame(
        {
            "chrom": raw_id_retained.map(lambda value: id_lookup_rows[value].chrom).values,
            "pos": raw_id_retained.map(lambda value: id_lookup_rows[value].pos).values,
            "id": raw_id_retained.values,
            "a1": a1_series.iloc[retained_positions].astype(str).values,
            "a2": a2_series.iloc[retained_positions].astype(str).values,
            "source_shard": ".",
            "source_index": source_index_series.iloc[retained_positions].values,
        },
        columns=["chrom", "pos", "id", "a1", "a2", "source_shard", "source_index"],
    )
    qc_rows_frame = qc_frame_from_reason(reason, source_index_series)
    return rows_frame, qc_rows_frame


def main() -> int:
    args = parse_args()
    meta_path = Path(args.sumstats_metadata)
    output_path = Path(args.output)
    logger.info("import_sumstats.py: importing sumstats -> %s", output_path)
    if args.input:
        reject_template_argument(args.input, label="import_sumstats.py --input")
    if args.id_lookup:
        reject_template_argument(args.id_lookup, label="import_sumstats.py --id-lookup")
    if not meta_path.exists():
        raise ValueError(f"metadata file not found: {meta_path}")
    metadata: Dict[str, object] = load_metadata(meta_path)
    input_path = resolve_sumstats_input_path(
        args.input,
        metadata_path=meta_path,
        metadata=metadata,
        consumer_label="import_sumstats.py",
    )
    reject_template_argument(str(input_path), label="import_sumstats.py --input")
    if not input_path.exists():
        raise ValueError(f"sumstats file not found: {input_path}")
    id_lookup_path: Optional[Path] = Path(args.id_lookup) if args.id_lookup else None
    if id_lookup_path is not None and not id_lookup_path.exists():
        raise ValueError(f"id lookup not found: {id_lookup_path}")
    id_lookup_rows: Dict[str, VariantRow] = {}
    ambiguous_lookup_ids: Set[str] = set()
    inherited_target_meta: Optional[Dict[str, object]] = None
    if id_lookup_path is not None:
        id_lookup_rows, ambiguous_lookup_ids, inherited_target_meta = load_id_lookup_object(id_lookup_path)
    _header_line, header, _delimiter, _header_line_number = read_sumstats_header(input_path)
    if id_lookup_path is None:
        variant_columns = resolve_variant_columns(header, metadata, require_pos=True)
    else:
        snp_idx, effect_idx, other_idx = resolve_id_enrichment_columns(header, metadata)

    rows_chunks: List[pd.DataFrame] = []
    qc_chunks: List[pd.DataFrame] = []
    for sumstats_chunk in iter_sumstats_table_chunks(input_path):
        frame = sumstats_chunk.frame
        source_index_series = sumstats_chunk.source_index
        if id_lookup_path is None:
            rows_chunk, qc_chunk = import_sumstats_variant_chunk(
                frame,
                source_index_series,
                variant_columns,
                max_allele_length=args.max_allele_length,
            )
        else:
            rows_chunk, qc_chunk = import_sumstats_id_enrichment_chunk(
                frame,
                source_index_series,
                snp_idx=snp_idx,
                effect_idx=effect_idx,
                other_idx=other_idx,
                id_lookup_rows=id_lookup_rows,
                ambiguous_lookup_ids=ambiguous_lookup_ids,
                max_allele_length=args.max_allele_length,
            )
        if not rows_chunk.empty:
            rows_chunks.append(rows_chunk)
        if not qc_chunk.empty:
            qc_chunks.append(qc_chunk)

    rows_frame = (
        pd.concat(rows_chunks, ignore_index=True)
        if rows_chunks
        else pd.DataFrame(columns=["chrom", "pos", "id", "a1", "a2", "source_shard", "source_index"])
    )
    qc_rows_frame = (
        pd.concat(qc_chunks, ignore_index=True)
        if qc_chunks
        else pd.DataFrame(columns=["source_shard", "source_index", "reason"])
    )

    allowed_chr, chr_filter_enabled = parse_chr2use(args.chr2use)
    if chr_filter_enabled and not rows_frame.empty:
        allowed_set = set(allowed_chr)
        canonical = rows_frame["chrom"].astype(str).map(canonical_contig_from_any_supported_label)
        keep_mask = canonical.isin(allowed_set)
        dropped = rows_frame.loc[~keep_mask, ["source_shard", "source_index"]].copy()
        if not dropped.empty:
            dropped["reason"] = "filtered_by_chr2use"
            qc_rows_frame = pd.concat([qc_rows_frame, dropped], ignore_index=True)
        rows_frame = rows_frame.loc[keep_mask].reset_index(drop=True)
    finalize_imported_vmap_vectorized(
        output_path=output_path,
        rows_frame=rows_frame,
        genome_build=args.genome_build if inherited_target_meta is None else str(inherited_target_meta["genome_build"]),
        target_contig_naming=None if inherited_target_meta is None else inherited_target_meta.get("contig_naming"),
        infer_target_contig_naming=inherited_target_meta is None,
        created_by="import_sumstats.py",
        derived_from=input_path,
        qc_rows_frame=qc_rows_frame,
    )
    logger.info(
        "import_sumstats.py: wrote %s with %s retained rows and %s QC rows",
        output_path,
        len(rows_frame),
        len(qc_rows_frame),
    )
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
