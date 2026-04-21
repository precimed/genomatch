#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

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
    find_metadata_value,
    is_missing_token_series,
    load_metadata,
    resolve_sumstats_input_path,
    read_sumstats_table,
    resolve_column,
    resolve_variant_columns,
)
from .vtable_utils import (
    load_metadata as load_variant_metadata,
    normalize_allele_token,
    open_text,
    parse_chr2use,
    validate_vtable_metadata,
    VariantRow,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract a .vmap from summary statistics.")
    parser.add_argument("--input", help="Input summary statistics file (optional when metadata defines path_sumStats)")
    parser.add_argument("--output", required=True, help="Output .vmap file")
    parser.add_argument("--sumstats-metadata", required=True, help="Cleansumstats-style metadata YAML")
    parser.add_argument("--genome-build", default="unknown", help="Genome build for metadata")
    parser.add_argument("--id-vtable", help="Optional .vtable for ID-based coordinate enrichment")
    parser.add_argument("--chr2use", "--contigs", dest="chr2use", help="Comma-separated chromosome list or ranges")
    return parser.parse_args()


def load_id_lookup_vtable(path: Path) -> tuple[Dict[str, VariantRow], Set[str], Dict[str, object]]:
    if path.suffix != ".vtable":
        raise ValueError("--id-vtable must point to a .vtable")
    metadata = load_variant_metadata(path)
    validate_vtable_metadata(metadata)
    unique_matches: Dict[str, VariantRow] = {}
    ambiguous_ids: Set[str] = set()
    ignored_ids = 0
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
            lookup_id = row_id.strip()
            if not lookup_id or lookup_id == ".":
                ignored_ids += 1
                continue
            row = VariantRow(chrom, pos, row_id, a1, a2)
            if lookup_id in ambiguous_ids:
                continue
            if lookup_id in unique_matches:
                unique_matches.pop(lookup_id, None)
                ambiguous_ids.add(lookup_id)
                continue
            unique_matches[lookup_id] = row
    if ignored_ids:
        logger.warning("ignored %s --id-vtable rows whose id is missing, empty, or '.'", ignored_ids)
    return unique_matches, ambiguous_ids, {
        "genome_build": metadata["genome_build"],
        "contig_naming": metadata.get("contig_naming"),
    }


def resolve_id_enrichment_columns(header: List[str], metadata: Dict[str, object]) -> tuple[int, int, int]:
    if find_metadata_value(metadata, "col_CHR") is not None or find_metadata_value(metadata, "col_POS") is not None:
        raise ValueError("--id-vtable requires metadata to omit col_CHR and col_POS")
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


def main() -> int:
    args = parse_args()
    meta_path = Path(args.sumstats_metadata)
    output_path = Path(args.output)
    logger.info("import_sumstats.py: importing sumstats -> %s", output_path)
    if args.input:
        reject_template_argument(args.input, label="import_sumstats.py --input")
    if args.id_vtable:
        reject_template_argument(args.id_vtable, label="import_sumstats.py --id-vtable")
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
    id_vtable_path: Optional[Path] = Path(args.id_vtable) if args.id_vtable else None
    if id_vtable_path is not None and not id_vtable_path.exists():
        raise ValueError(f"id-vtable not found: {id_vtable_path}")
    id_lookup_rows: Dict[str, VariantRow] = {}
    ambiguous_lookup_ids: Set[str] = set()
    inherited_target_meta: Optional[Dict[str, object]] = None
    if id_vtable_path is not None:
        id_lookup_rows, ambiguous_lookup_ids, inherited_target_meta = load_id_lookup_vtable(id_vtable_path)
    sumstats_table = read_sumstats_table(input_path)
    header = list(sumstats_table.header)
    frame = sumstats_table.frame
    source_index_series = sumstats_table.source_index
    if id_vtable_path is None:
        variant_columns = resolve_variant_columns(header, metadata, require_pos=True)
    else:
        snp_idx, effect_idx, other_idx = resolve_id_enrichment_columns(header, metadata)

    rows_frame = pd.DataFrame(columns=["chrom", "pos", "id", "a1", "a2", "source_shard", "source_index"])
    qc_rows_frame = pd.DataFrame(columns=["source_shard", "source_index", "reason"])
    if id_vtable_path is None:
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

            missing_required_mask = valid_mask & (
                is_missing_token_series(chrom_series)
                | is_missing_token_series(pos_series)
                | is_missing_token_series(a1_series)
                | is_missing_token_series(a2_series)
            )
            reason.loc[missing_required_mask] = "malformed_row"

            if variant_columns.snp is not None:
                snp_raw = frame.iloc[active_idx, variant_columns.snp].astype(str)
                row_id_series.iloc[active_idx] = snp_raw.where(snp_raw != "", ".")

            non_actg_mask = reason.isna() & (
                ~a1_series.map(is_canonical_allele_token).fillna(False)
                | ~a2_series.map(is_canonical_allele_token).fillna(False)
            )
            reason.loc[non_actg_mask] = "non_actg_allele"

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
        qc_positions = reason[reason.notna()].index
        qc_rows_frame = pd.DataFrame(
            {
                "source_shard": ".",
                "source_index": source_index_series.iloc[qc_positions].values,
                "reason": reason.iloc[qc_positions].astype(str).values,
            },
            columns=["source_shard", "source_index", "reason"],
        )
    else:
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

            missing_required_mask = valid_mask & (
                is_missing_token_series(a1_series)
                | is_missing_token_series(a2_series)
            )
            reason.loc[missing_required_mask] = "malformed_row"

            non_actg_mask = reason.isna() & (
                ~a1_series.map(is_canonical_allele_token).fillna(False)
                | ~a2_series.map(is_canonical_allele_token).fillna(False)
            )
            reason.loc[non_actg_mask] = "non_actg_allele"

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
        qc_positions = reason[reason.notna()].index
        qc_rows_frame = pd.DataFrame(
            {
                "source_shard": ".",
                "source_index": source_index_series.iloc[qc_positions].values,
                "reason": reason.iloc[qc_positions].astype(str).values,
            },
            columns=["source_shard", "source_index", "reason"],
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
