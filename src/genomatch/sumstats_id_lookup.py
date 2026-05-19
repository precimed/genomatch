#!/usr/bin/env python3
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterator, List

import pandas as pd

from .exact_set_utils import TargetObjectInfo, load_target_object_info, validate_loaded_row_count
from .sumstats_utils import find_metadata_value, resolve_column, VariantColumnMapping
from .vtable_utils import iter_vmap_table_chunks, iter_vtable_table_chunks

logger = logging.getLogger(__name__)


def load_id_lookup_info(path: Path) -> tuple[TargetObjectInfo, dict[str, object]]:
    if not (path.name.endswith(".vtable") or path.name.endswith(".vmap")):
        raise ValueError("--id-lookup must point to a .vtable or .vmap")
    info = load_target_object_info(path)
    return info, {
        "genome_build": info.target_metadata["genome_build"],
        "contig_naming": info.target_metadata.get("contig_naming"),
    }


def iter_id_lookup_target_chunks(path: Path) -> Iterator[pd.DataFrame]:
    if path.name.endswith(".vmap"):
        for table in iter_vmap_table_chunks(path, check_duplicates=False):
            yield table.to_frame(copy=False).loc[:, ["chrom", "pos", "id", "a1", "a2"]]
    else:
        for table in iter_vtable_table_chunks(path):
            yield table.to_frame(copy=False)


def raise_ambiguous_id_lookup_error(lookup_id: object, first: pd.Series, second: pd.Series) -> None:
    raise ValueError(
        f"--id-lookup is ambiguous for source ID {str(lookup_id)!r}: "
        f"{first['lookup_chrom']}:{first['lookup_pos']} and {second['lookup_chrom']}:{second['lookup_pos']}"
    )


def require_unambiguous_lookup_coordinates(coordinates: pd.DataFrame) -> None:
    duplicate_id_mask = coordinates["id"].duplicated(keep=False)
    if not bool(duplicate_id_mask.any()):
        return
    first_id = coordinates.loc[duplicate_id_mask, "id"].iloc[0]
    matches = coordinates.loc[coordinates["id"].eq(first_id)]
    raise_ambiguous_id_lookup_error(first_id, matches.iloc[0], matches.iloc[1])


def require_compatible_filled_lookup_coordinates(coordinates: pd.DataFrame, rows_frame: pd.DataFrame) -> None:
    if coordinates.empty:
        return
    filled_mask = rows_frame["chrom"].astype(str).str.strip().ne("") & rows_frame["pos"].astype(str).str.strip().ne("")
    if not bool(filled_mask.any()):
        return
    existing_coordinates = rows_frame.loc[
        filled_mask & rows_frame["id"].isin(pd.Index(coordinates["id"])),
        ["id", "chrom", "pos"],
    ].rename(columns={"chrom": "lookup_chrom", "pos": "lookup_pos"})
    if existing_coordinates.empty:
        return
    existing_coordinates = existing_coordinates.drop_duplicates(
        subset=["id", "lookup_chrom", "lookup_pos"],
        keep="first",
    )
    overlap = coordinates.merge(
        existing_coordinates,
        on="id",
        how="inner",
        suffixes=("_new", "_seen"),
        sort=False,
    )
    if overlap.empty:
        return
    conflict_mask = (
        overlap["lookup_chrom_new"].ne(overlap["lookup_chrom_seen"])
        | overlap["lookup_pos_new"].ne(overlap["lookup_pos_seen"])
    )
    if not bool(conflict_mask.any()):
        return
    conflict = overlap.loc[conflict_mask].iloc[0]
    raise ValueError(
        f"--id-lookup is ambiguous for source ID {str(conflict['id'])!r}: "
        f"{conflict['lookup_chrom_seen']}:{conflict['lookup_pos_seen']} and "
        f"{conflict['lookup_chrom_new']}:{conflict['lookup_pos_new']}"
    )


def augment_rows_frame_with_id_lookup(
    rows_frame: pd.DataFrame,
    path: Path,
    info: TargetObjectInfo,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Assumes: rows_frame is the retained imported sumstats frame with valid non-missing IDs
    # and placeholder chrom/pos columns in --id-lookup mode.
    # Performs: SV(streamed lookup ID filtering, duplicate-coordinate ambiguity checks, vectorized merge assignment).
    # Guarantees: chrom/pos filled for matched IDs; unmatched rows returned as id_not_found QC rows.
    if rows_frame.empty:
        return rows_frame, pd.DataFrame(columns=["source_shard", "source_index", "reason"])

    needed_ids = pd.Index(rows_frame["id"].astype(str).drop_duplicates())
    row_idx_column = "_row_idx"
    rows_frame[row_idx_column] = pd.RangeIndex(len(rows_frame))
    ignored_ids = 0
    observed_rows = 0

    for chunk in iter_id_lookup_target_chunks(path):
        observed_rows += len(chunk)
        lookup_id_series = chunk["id"].astype(str).str.strip()
        invalid_id_mask = lookup_id_series.eq("") | lookup_id_series.eq(".")
        ignored_ids += int(invalid_id_mask.sum())
        relevant_mask = ~invalid_id_mask & lookup_id_series.isin(needed_ids)
        if not bool(relevant_mask.any()):
            continue

        coordinates = pd.DataFrame(
            {
                "id": lookup_id_series.loc[relevant_mask].to_numpy(dtype=object),
                "lookup_chrom": chunk.loc[relevant_mask, "chrom"].astype(str).to_numpy(dtype=object),
                "lookup_pos": chunk.loc[relevant_mask, "pos"].astype(str).to_numpy(dtype=object),
            }
        ).drop_duplicates(subset=["id", "lookup_chrom", "lookup_pos"], keep="first")
        require_unambiguous_lookup_coordinates(coordinates)
        require_compatible_filled_lookup_coordinates(coordinates, rows_frame)

        assignments = rows_frame.loc[:, [row_idx_column, "id"]].merge(coordinates, on="id", how="inner", sort=False)
        if assignments.empty:
            continue
        row_positions = assignments[row_idx_column].to_numpy()
        rows_frame.loc[row_positions, "chrom"] = assignments["lookup_chrom"].to_numpy(dtype=object)
        rows_frame.loc[row_positions, "pos"] = assignments["lookup_pos"].to_numpy(dtype=object)

    validate_loaded_row_count(info, observed_rows)
    if ignored_ids:
        logger.warning("ignored %s --id-lookup rows whose id is missing, empty, or '.'", ignored_ids)

    matched_mask = rows_frame["chrom"].astype(str).str.strip().ne("") & rows_frame["pos"].astype(str).str.strip().ne("")
    dropped = rows_frame.loc[~matched_mask, ["source_shard", "source_index"]].copy()
    if not dropped.empty:
        dropped["reason"] = "id_not_found"
    return rows_frame.loc[matched_mask].drop(columns=[row_idx_column]).reset_index(drop=True), dropped


def resolve_id_enrichment_columns(header: List[str], metadata: dict[str, object]) -> VariantColumnMapping:
    if find_metadata_value(metadata, "col_CHR") is not None or find_metadata_value(metadata, "col_POS") is not None:
        raise ValueError("--id-lookup requires metadata to omit col_CHR and col_POS")
    return VariantColumnMapping(
        chr=None,
        pos=None,
        snp=resolve_column(header, find_metadata_value(metadata, "col_SNP"), "col_SNP", required=True),
        effect_allele=resolve_column(
            header,
            find_metadata_value(metadata, "col_EffectAllele"),
            "col_EffectAllele",
            required=True,
        ),
        other_allele=resolve_column(
            header,
            find_metadata_value(metadata, "col_OtherAllele"),
            "col_OtherAllele",
            required=True,
        ),
    )
