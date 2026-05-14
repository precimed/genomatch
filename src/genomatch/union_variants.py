#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from ._cli_utils import run_cli
from .exact_set_utils import (
    TargetObjectInfo,
    assign_target_derived_ids,
    drop_duplicate_target_identities,
    load_target_object_info,
    read_target_table,
    require_shared_target_metadata,
    validate_loaded_row_count,
    vtable_metadata_from_target_info,
)
from .tabular_rows import VariantRowsTable
from .vtable_utils import (
    metadata_with_variants_count,
    require_contig_naming,
    sort_target_table_by_declared_coordinate,
    write_metadata,
    write_vtable_table,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Union 2+ .vtable/.vmap inputs on exact variant rows, "
            "write target-derived IDs, and sort the result into declared coordinate order."
        )
    )
    parser.add_argument("inputs", nargs="+", help="Input .vtable/.vmap files")
    parser.add_argument("--output", required=True, help="Output .vtable")
    return parser.parse_args()


def print_loaded(path: Path, row_count: int) -> None:
    logger.info("union_variants.py: loaded %s variants from %s", row_count, path)


def print_accumulated(path: Path, row_count: int) -> None:
    logger.info("union_variants.py: after unioning %s, %s variants accumulated", path, row_count)


def missing_count_paths(infos: list[TargetObjectInfo]) -> list[Path]:
    return [info.path for info in infos if info.variants_count is None]


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
    missing_counts = missing_count_paths(infos)
    if missing_counts:
        logger.warning(
            "union_variants.py: variants_count missing from metadata for %s; using CLI input order",
            ", ".join(str(path) for path in missing_counts),
        )
        ordered_infos = infos
    else:
        declared_rank = {path: index for index, path in enumerate(input_paths)}
        ordered_infos = sorted(infos, key=lambda info: (info.variants_count, declared_rank[info.path]))

    accumulated = read_target_table(ordered_infos[0].path).to_frame(copy=False)
    validate_loaded_row_count(ordered_infos[0], len(accumulated))
    print_loaded(ordered_infos[0].path, len(accumulated))
    print_accumulated(ordered_infos[0].path, len(accumulated))
    for info in ordered_infos[1:]:
        table = read_target_table(info.path)
        frame = table.to_frame(copy=False)
        validate_loaded_row_count(info, len(frame))
        print_loaded(info.path, len(frame))
        accumulated = pd.concat([accumulated, frame], ignore_index=True, copy=False)
        print_accumulated(info.path, len(accumulated))

    sorted_frame = sort_target_table_by_declared_coordinate(accumulated, contig_naming, label="variant table")
    sorted_frame = assign_target_derived_ids(sorted_frame)
    sorted_frame = drop_duplicate_target_identities(sorted_frame)
    write_vtable_table(output_path, VariantRowsTable.from_frame(sorted_frame, copy=False), assume_validated=True)
    write_metadata(output_path, metadata_with_variants_count(vtable_metadata_from_target_info(infos[0]), len(sorted_frame)))
    logger.info("union_variants.py: wrote %s with %s union rows", output_path, len(sorted_frame))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
