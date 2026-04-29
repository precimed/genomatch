#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from ._cli_utils import run_cli
from .tabular_rows import VariantRowsTable, VMapRowsTable
from .vtable_utils import (
    load_metadata,
    read_vmap_table,
    read_vtable_table,
    require_contig_naming,
    sort_target_table_by_declared_coordinate,
    validate_vmap_metadata,
    validate_vtable_metadata,
    write_metadata,
    write_vmap_table,
    write_vtable_table,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sort target rows in a .vtable/.vmap by declared contig order, then numeric position. "
            "The input object type and any .vmap provenance are preserved."
        )
    )
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable or .vmap")
    parser.add_argument(
        "--drop-duplicates",
        action="store_true",
        help="Drop duplicate target identities (chrom,pos,a1,a2), retaining first in stable sorted order",
    )
    return parser.parse_args()


def object_type_for_path(path: Path) -> str:
    if path.name.endswith(".vtable"):
        return "variant_table"
    if path.name.endswith(".vmap"):
        return "variant_map"
    raise ValueError(f"unsupported input: {path}")


def load_sort_input(path: Path, *, allow_vmap_duplicates: bool) -> tuple[str, dict, str, pd.DataFrame]:
    object_type = object_type_for_path(path)
    metadata = load_metadata(path)
    if object_type == "variant_table":
        validate_vtable_metadata(metadata)
        contig_naming = require_contig_naming(metadata, label="variant object")
        frame = read_vtable_table(path).to_frame(copy=False)
    else:
        validate_vmap_metadata(metadata)
        target_metadata = dict(metadata["target"])
        contig_naming = require_contig_naming(target_metadata, label="variant object")
        frame = read_vmap_table(path, check_duplicates=not allow_vmap_duplicates).to_frame(copy=False)
    return object_type, metadata, contig_naming, frame


def drop_duplicate_target_identities(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    duplicate_mask = frame.duplicated(subset=["chrom", "pos", "a1", "a2"], keep="first")
    dropped = int(duplicate_mask.sum())
    if not dropped:
        return frame, 0
    return frame.loc[~duplicate_mask].reset_index(drop=True), dropped


def write_sorted_object(path: Path, *, object_type: str, frame: pd.DataFrame) -> None:
    if object_type == "variant_table":
        write_vtable_table(path, VariantRowsTable.from_frame(frame, copy=False), assume_validated=True)
        return
    write_vmap_table(path, VMapRowsTable.from_frame(frame, copy=False), assume_validated=True)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    logger.info("sort_variants.py: sorting %s -> %s", input_path, output_path)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")
    if object_type_for_path(input_path) != object_type_for_path(output_path):
        raise ValueError("sort_variants.py output suffix must match input suffix")

    object_type, metadata, contig_naming, frame = load_sort_input(
        input_path,
        allow_vmap_duplicates=args.drop_duplicates,
    )
    sorted_frame = sort_target_table_by_declared_coordinate(frame, contig_naming, label="variant object")
    if args.drop_duplicates:
        sorted_frame, dropped = drop_duplicate_target_identities(sorted_frame)
    else:
        dropped = 0

    write_sorted_object(output_path, object_type=object_type, frame=sorted_frame)
    write_metadata(output_path, dict(metadata))
    qc_path = output_path.with_name(output_path.name + ".qc.tsv")
    if qc_path.exists():
        qc_path.unlink()

    if args.drop_duplicates:
        logger.info("sort_variants.py: dropped %s duplicate target rows.", dropped)
    logger.info("sort_variants.py: retained %s rows in declared coordinate order.", len(sorted_frame))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
