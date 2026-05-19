#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterator

import pandas as pd

from ._cli_utils import run_cli
from .exact_set_utils import (
    exact_key_index,
    format_target_key_ids,
    load_target_object_info,
    require_shared_target_metadata,
    TARGET_KEY_COLUMNS,
    validate_loaded_row_count,
)
from .tabular_rows import VMapRowsTable
from .vtable_utils import (
    iter_vmap_table_chunks,
    iter_vtable_table_chunks,
    metadata_with_variants_count,
    read_vmap_table,
    write_metadata,
    write_vmap_status_qc,
    write_vmap_table,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Assign IDs in a .vmap from a .vmap/.vtable ID source by exact chrom:pos:a1:a2 match."
        )
    )
    parser.add_argument("--vmap", required=True, help="Source .vmap to update")
    parser.add_argument("--id-source", required=True, help="ID source .vmap or .vtable")
    parser.add_argument("--output", required=True, help="Output .vmap")
    parser.add_argument(
        "--unmatched-id-policy",
        choices=["drop", "variant-key", "missing"],
        default="drop",
        help=(
            "How to handle source .vmap rows absent from --id-source: drop and audit them (default), "
            "keep them with id=chrom:pos:a1:a2, or keep them with id='.'. "
            "Rows matched to missing ID-source IDs are always dropped and audited."
        ),
    )
    parser.add_argument(
        "--duplicate-id-policy",
        choices=["allow", "fail", "drop-all"],
        default="fail",
        help=(
            "How to handle duplicate non-missing IDs in the retained output .vmap: "
            "fail clearly (default), allow them, or drop all rows carrying duplicated IDs to QC."
        ),
    )
    return parser.parse_args()


def iter_id_source_target_chunks(path: Path) -> Iterator[pd.DataFrame]:
    if path.name.endswith(".vmap"):
        for table in iter_vmap_table_chunks(path, check_duplicates=False):
            yield table.to_frame(copy=False).loc[:, ["chrom", "pos", "id", "a1", "a2"]]
    else:
        for table in iter_vtable_table_chunks(path):
            yield table.to_frame(copy=False)


def raise_duplicate_used_key(path: Path, row: pd.Series) -> None:
    key = ":".join(str(row[column]) for column in TARGET_KEY_COLUMNS)
    raise ValueError(f"{path} contains duplicate used chrom:pos:a1:a2 key: {key}")


def require_no_duplicate_candidate_keys(coordinates: pd.DataFrame, *, path: Path) -> None:
    duplicate_mask = coordinates.duplicated(subset=TARGET_KEY_COLUMNS, keep=False)
    if bool(duplicate_mask.any()):
        raise_duplicate_used_key(path, coordinates.loc[duplicate_mask].iloc[0])


def duplicated_output_id_mask(frame: pd.DataFrame) -> pd.Series:
    ids = frame["id"].astype(str).str.strip()
    non_missing_id_mask = ~ids.isin({"", "."})
    return non_missing_id_mask & ids.duplicated(keep=False)


def raise_duplicate_output_id(frame: pd.DataFrame) -> None:
    duplicate_id = str(frame.loc[duplicated_output_id_mask(frame), "id"].iloc[0])
    raise ValueError(
        f"output .vmap would contain duplicate non-missing ID {duplicate_id!r}; "
        "use --duplicate-id-policy allow or --duplicate-id-policy drop-all to override"
    )


def assign_ids_from_id_source(
    source_frame: pd.DataFrame,
    id_source_path: Path,
    id_source_info,
    *,
    unmatched_id_policy: str,
    duplicate_id_policy: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Assumes: source_frame is a loaded, unique-target source .vmap frame.
    # Performs: SV(streamed ID-source filtering, duplicate-used-key checks, vectorized ID assignment).
    # Guarantees: retained frame has assigned IDs; dropped frame has id_not_found or missing_id status.
    source_keys = exact_key_index(source_frame)
    row_idx_column = "_row_idx"
    assigned_id_column = "_assigned_id"
    source_id_column = "_source_id"
    source_frame[row_idx_column] = pd.RangeIndex(len(source_frame))
    source_frame[assigned_id_column] = pd.Series([None] * len(source_frame), dtype="object")
    observed_rows = 0

    for chunk in iter_id_source_target_chunks(id_source_path):
        observed_rows += len(chunk)
        chunk_keys = exact_key_index(chunk)
        used_mask = chunk_keys.isin(source_keys)
        if not bool(used_mask.any()):
            continue

        coordinates = chunk.loc[used_mask, [*TARGET_KEY_COLUMNS, "id"]].rename(columns={"id": "_candidate_id"})
        require_no_duplicate_candidate_keys(coordinates, path=id_source_path)

        assignments = source_frame.loc[:, [row_idx_column, assigned_id_column, *TARGET_KEY_COLUMNS]].merge(
            coordinates,
            on=TARGET_KEY_COLUMNS,
            how="inner",
            sort=False,
        )
        if assignments.empty:
            continue
        duplicate_assigned_mask = assignments[assigned_id_column].notna()
        if bool(duplicate_assigned_mask.any()):
            raise_duplicate_used_key(id_source_path, assignments.loc[duplicate_assigned_mask].iloc[0])

        row_positions = assignments[row_idx_column].to_numpy()
        source_frame.loc[row_positions, assigned_id_column] = assignments["_candidate_id"].astype(str).to_numpy(
            dtype=object
        )

    validate_loaded_row_count(id_source_info, observed_rows)

    assigned_ids = source_frame[assigned_id_column]
    matched_mask = assigned_ids.notna()
    missing_id_mask = matched_mask & assigned_ids.astype(str).str.strip().isin({"", "."})
    unmatched_mask = ~matched_mask
    retain_mask = (matched_mask & ~missing_id_mask) | (unmatched_mask & (unmatched_id_policy != "drop"))

    out_frame = source_frame.loc[retain_mask].copy()
    out_frame[source_id_column] = out_frame["id"].astype(str)
    out_assigned_ids = assigned_ids.loc[retain_mask].copy()
    out_unmatched_mask = out_assigned_ids.isna()
    if bool(out_unmatched_mask.any()):
        if unmatched_id_policy == "variant-key":
            unmatched_frame = out_frame.loc[out_unmatched_mask, TARGET_KEY_COLUMNS]
            out_assigned_ids.loc[out_unmatched_mask] = format_target_key_ids(unmatched_frame)
        elif unmatched_id_policy == "missing":
            out_assigned_ids.loc[out_unmatched_mask] = "."
    out_frame["id"] = out_assigned_ids.astype(str).to_numpy()

    dropped_frame = source_frame.loc[~retain_mask].copy()
    if not dropped_frame.empty:
        dropped_frame["status"] = "missing_id"
        dropped_frame.loc[unmatched_mask.loc[dropped_frame.index], "status"] = "id_not_found"
        dropped_frame.loc[missing_id_mask.loc[dropped_frame.index], "status"] = "missing_id"

    if duplicate_id_policy != "allow":
        duplicate_id_mask = duplicated_output_id_mask(out_frame)
        if bool(duplicate_id_mask.any()):
            if duplicate_id_policy == "fail":
                raise_duplicate_output_id(out_frame)
            duplicate_dropped_frame = out_frame.loc[duplicate_id_mask].copy()
            duplicate_dropped_frame["status"] = "duplicate_id"
            duplicate_dropped_frame["id"] = duplicate_dropped_frame[source_id_column].astype(str)
            dropped_frame = pd.concat([dropped_frame, duplicate_dropped_frame], axis=0)
            out_frame = out_frame.loc[~duplicate_id_mask].copy()

    out_frame = out_frame.drop(columns=[row_idx_column, assigned_id_column, source_id_column])
    if not dropped_frame.empty:
        dropped_frame = dropped_frame.sort_values(row_idx_column, kind="stable")
        dropped_frame = dropped_frame.drop(
            columns=[row_idx_column, assigned_id_column, source_id_column],
            errors="ignore",
        )
    return out_frame.reset_index(drop=True), dropped_frame.reset_index(drop=True)


def main() -> int:
    args = parse_args()
    vmap_path = Path(args.vmap)
    id_source_path = Path(args.id_source)
    output_path = Path(args.output)
    logger.info("assign_vmap_ids.py: assigning IDs for %s from %s -> %s", vmap_path, id_source_path, output_path)

    if "@" in args.vmap or "@" in args.id_source or "@" in args.output:
        raise ValueError("assign_vmap_ids.py does not accept '@' paths")
    if not vmap_path.name.endswith(".vmap"):
        raise ValueError("assign_vmap_ids.py requires --vmap to be a .vmap")
    if not id_source_path.name.endswith((".vmap", ".vtable")):
        raise ValueError("assign_vmap_ids.py requires --id-source to be a .vmap or .vtable")
    if not output_path.name.endswith(".vmap"):
        raise ValueError("assign_vmap_ids.py requires --output to end with .vmap")

    vmap_info = load_target_object_info(vmap_path)
    id_source_info = load_target_object_info(id_source_path)
    require_shared_target_metadata([vmap_info, id_source_info])

    source_table = read_vmap_table(vmap_path)
    source_frame = source_table.to_frame(copy=False)
    validate_loaded_row_count(vmap_info, len(source_frame))

    out_frame, dropped_frame = assign_ids_from_id_source(
        source_frame,
        id_source_path,
        id_source_info,
        unmatched_id_policy=args.unmatched_id_policy,
        duplicate_id_policy=args.duplicate_id_policy,
    )
    qc_path = output_path.with_name(output_path.name + ".qc.tsv")
    if dropped_frame.empty:
        if qc_path.exists():
            qc_path.unlink()
    else:
        write_vmap_status_qc(
            qc_path,
            zip(
                dropped_frame["source_shard"].astype(str),
                dropped_frame["source_index"].astype(int),
                dropped_frame["id"].astype(str),
                dropped_frame["status"].astype(str),
            ),
        )

    write_vmap_table(output_path, VMapRowsTable.from_frame(out_frame, copy=False))
    write_metadata(output_path, metadata_with_variants_count(vmap_info.raw_metadata, len(out_frame)))
    logger.info(
        "assign_vmap_ids.py: wrote %s with %s assigned rows; dropped %s rows",
        output_path,
        len(out_frame),
        len(dropped_frame),
    )
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
