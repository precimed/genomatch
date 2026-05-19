#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List, Optional

import pandas as pd

from ._cli_utils import run_cli
from .importer_utils import (
    finalize_imported_vmap_vectorized,
    is_canonical_allele_token,
    is_valid_import_position,
    reject_template_argument,
)
from .contig_utils import canonical_contig_from_any_supported_label
from .sumstats_utils import (
    extract_variant_field,
    iter_sumstats_table_chunks,
    is_missing_token_series,
    load_metadata,
    read_sumstats_header,
    resolve_sumstats_input_path,
    resolve_variant_columns,
    VariantColumnMapping,
)
from .sumstats_id_lookup import (
    augment_rows_frame_with_id_lookup,
    load_id_lookup_info,
    resolve_id_enrichment_columns,
)
from .vtable_utils import (
    normalize_allele_token,
    parse_chr2use,
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
    allow_missing_coordinates: bool = False,
    require_valid_id: bool = False,
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
        a1_raw = frame.iloc[active_idx, variant_columns.effect_allele].astype(str)
        a2_raw = frame.iloc[active_idx, variant_columns.other_allele].astype(str)

        if variant_columns.chr is not None:
            chrom_raw = frame.iloc[active_idx, variant_columns.chr].astype(str)
            chrom_series.iloc[active_idx] = chrom_raw.map(lambda value: extract_variant_field(value, "CHR"))
        if variant_columns.pos is not None:
            pos_raw = frame.iloc[active_idx, variant_columns.pos].astype(str)
            pos_series.iloc[active_idx] = pos_raw.map(lambda value: extract_variant_field(value, "POS"))
        a1_series.iloc[active_idx] = a1_raw.map(lambda value: normalize_allele_token(extract_variant_field(value, "EffectAllele")))
        a2_series.iloc[active_idx] = a2_raw.map(lambda value: normalize_allele_token(extract_variant_field(value, "OtherAllele")))

        if variant_columns.snp is not None:
            snp_raw = frame.iloc[active_idx, variant_columns.snp].astype(str)
            row_id_series.iloc[active_idx] = snp_raw.astype(str).str.strip().where(snp_raw.astype(str).str.strip() != "", ".")
        if require_valid_id:
            invalid_id_mask = reason.isna() & (
                row_id_series.astype(str).str.strip().eq("")
                | row_id_series.astype(str).str.strip().eq(".")
            )
            reason.loc[invalid_id_mask] = "invalid_id"

        apply_common_allele_validation(
            reason,
            a1_series,
            a2_series,
            max_allele_length=max_allele_length,
        )

        missing_chrom_mask = chrom_series.astype(str).str.strip().eq("")
        invalid_pos_mask = ~pos_series.map(is_valid_import_position).fillna(False)
        malformed_position_mask = missing_chrom_mask | invalid_pos_mask
        if allow_missing_coordinates:
            missing_pos_mask = pos_series.astype(str).str.strip().eq("")
            any_coordinate_present_mask = ~missing_chrom_mask | ~missing_pos_mask
            malformed_position_mask = any_coordinate_present_mask & malformed_position_mask
        reason.loc[reason.isna() & malformed_position_mask] = "malformed_row"

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
    metadata: dict[str, object] = load_metadata(meta_path)
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
    id_lookup_info = None
    inherited_target_meta: Optional[dict[str, object]] = None
    if id_lookup_path is not None:
        id_lookup_info, inherited_target_meta = load_id_lookup_info(id_lookup_path)
    _header_line, header, _delimiter, _header_line_number = read_sumstats_header(input_path)
    if id_lookup_path is None:
        variant_columns = resolve_variant_columns(header, metadata, require_pos=True)
    else:
        variant_columns = resolve_id_enrichment_columns(header, metadata)

    rows_chunks: List[pd.DataFrame] = []
    qc_chunks: List[pd.DataFrame] = []
    for sumstats_chunk in iter_sumstats_table_chunks(input_path):
        frame = sumstats_chunk.frame
        source_index_series = sumstats_chunk.source_index
        rows_chunk, qc_chunk = import_sumstats_variant_chunk(
            frame,
            source_index_series,
            variant_columns,
            allow_missing_coordinates=id_lookup_path is not None,
            require_valid_id=id_lookup_path is not None,
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

    if id_lookup_path is not None:
        if id_lookup_info is None:
            raise AssertionError("id_lookup_info missing")
        rows_frame, id_lookup_qc_frame = augment_rows_frame_with_id_lookup(rows_frame, id_lookup_path, id_lookup_info)
        if not id_lookup_qc_frame.empty:
            qc_rows_frame = pd.concat([qc_rows_frame, id_lookup_qc_frame], ignore_index=True)

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
