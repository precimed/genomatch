#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from ._cli_utils import run_cli
from .exact_set_utils import (
    exact_key_index,
    load_target_object_info,
    read_target_table,
    require_shared_target_metadata,
    vtable_metadata_from_target_info,
)
from .tabular_rows import VariantRowsTable
from .vtable_utils import (
    require_contig_naming,
    sort_target_table_by_declared_coordinate,
    write_metadata,
    write_vtable_table,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Union 2+ .vtable/.vmap inputs on exact variant rows, keep the first occurrence of each row, "
            "and sort the result into declared coordinate order."
        )
    )
    parser.add_argument("inputs", nargs="+", help="Input .vtable/.vmap files")
    parser.add_argument("--output", required=True, help="Output .vtable")
    return parser.parse_args()


def print_loaded(path: Path, row_count: int) -> None:
    logger.info("union_variants.py: loaded %s variants from %s", row_count, path)


def print_accumulated(path: Path, row_count: int) -> None:
    logger.info("union_variants.py: after unioning %s, %s variants accumulated", path, row_count)


def drop_duplicate_target_identities(frame: pd.DataFrame) -> pd.DataFrame:
    duplicate_mask = exact_key_index(frame).duplicated(keep="first")
    if not duplicate_mask.any():
        return frame
    return frame.loc[~duplicate_mask].reset_index(drop=True)


def main() -> int:
    args = parse_args()
    if len(args.inputs) < 2:
        raise ValueError("union_variants.py requires at least two inputs")
    if "@" in args.output or any("@" in item for item in args.inputs):
        raise ValueError("union_variants.py does not accept '@' paths")
    output_path = Path(args.output)
    if not output_path.name.endswith(".vtable"):
        raise ValueError("union_variants.py requires --output to end with .vtable")

    input_paths = [Path(item) for item in args.inputs]
    logger.info("union_variants.py: unioning %s inputs -> %s", len(input_paths), output_path)
    infos = [load_target_object_info(path) for path in input_paths]
    require_shared_target_metadata(infos)
    contig_naming = require_contig_naming(infos[0].target_metadata, label="variant table")

    accumulated = pd.DataFrame(columns=["chrom", "pos", "id", "a1", "a2"], dtype="object")
    for path in input_paths:
        table = read_target_table(path)
        frame = table.to_frame(copy=False)
        print_loaded(path, len(frame))
        accumulated = pd.concat([accumulated, frame], ignore_index=True, copy=False)
        print_accumulated(path, len(accumulated))

    sorted_frame = sort_target_table_by_declared_coordinate(accumulated, contig_naming, label="variant table")
    sorted_frame = drop_duplicate_target_identities(sorted_frame)
    write_vtable_table(output_path, VariantRowsTable.from_frame(sorted_frame, copy=False), assume_validated=True)
    write_metadata(output_path, vtable_metadata_from_target_info(infos[0]))
    logger.info("union_variants.py: wrote %s with %s union rows", output_path, len(sorted_frame))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
