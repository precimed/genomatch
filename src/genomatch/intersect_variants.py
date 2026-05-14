#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ._cli_utils import run_cli
from .exact_set_utils import (
    load_target_object_info,
    read_target_table,
    require_shared_target_metadata,
    restrict_frame_to_membership_chunks,
    vtable_metadata_from_target_info,
)
from .tabular_rows import VariantRowsTable
from .vtable_utils import write_metadata, write_vtable_table

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Intersect 2+ .vtable/.vmap inputs on exact variant rows.")
    parser.add_argument("inputs", nargs="+", help="Input .vtable/.vmap files")
    parser.add_argument("--output", required=True, help="Output .vtable")
    return parser.parse_args()


def print_loaded(path: Path, row_count: int) -> None:
    logger.info("intersect_variants.py: loaded %s variants from %s", row_count, path)


def print_remaining(path: Path, row_count: int) -> None:
    logger.info("intersect_variants.py: after intersecting %s, %s variants remain", path, row_count)


def main() -> int:
    args = parse_args()
    if len(args.inputs) < 2:
        raise ValueError("intersect_variants.py requires at least two inputs")
    if "@" in args.output or any("@" in item for item in args.inputs):
        raise ValueError("intersect_variants.py does not accept '@' paths")
    output_path = Path(args.output)
    if not output_path.name.endswith(".vtable"):
        raise ValueError("intersect_variants.py requires --output to end with .vtable")

    input_paths = [Path(item) for item in args.inputs]
    logger.info("intersect_variants.py: intersecting %s inputs -> %s", len(input_paths), output_path)
    infos = [load_target_object_info(path) for path in input_paths]
    require_shared_target_metadata(infos)

    first_table = read_target_table(input_paths[0])
    retained_frame = first_table.to_frame(copy=False)
    print_loaded(input_paths[0], len(retained_frame))
    for path in input_paths[1:]:
        retained_frame, row_count = restrict_frame_to_membership_chunks(retained_frame, path)
        print_loaded(path, row_count)
        print_remaining(path, len(retained_frame))

    write_vtable_table(output_path, VariantRowsTable.from_frame(retained_frame, copy=False))
    write_metadata(output_path, vtable_metadata_from_target_info(infos[0]))
    logger.info("intersect_variants.py: wrote %s with %s intersected rows", output_path, len(retained_frame))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
