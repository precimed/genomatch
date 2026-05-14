#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ._cli_utils import run_cli
from .exact_set_utils import (
    load_target_object_info,
    read_mapped_vmap_table,
    require_shared_target_metadata,
    restrict_frame_to_membership_chunks,
)
from .tabular_rows import VMapRowsTable
from .vtable_utils import write_metadata, write_vmap_table

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Restrict a mapped-only .vmap by exact membership in one or more .vtable/.vmap inputs."
    )
    parser.add_argument("source", help="Source mapped-only .vmap")
    parser.add_argument("restrictions", nargs="+", help="Restriction .vtable/.vmap files")
    parser.add_argument("--output", required=True, help="Output .vmap")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_path = Path(args.source)
    restriction_paths = [Path(item) for item in args.restrictions]
    output_path = Path(args.output)
    logger.info(
        "restrict_vmap.py: restricting %s by %s membership inputs -> %s",
        source_path,
        len(restriction_paths),
        output_path,
    )
    if "@" in args.source or "@" in args.output or any("@" in item for item in args.restrictions):
        raise ValueError("restrict_vmap.py does not accept '@' paths")
    if not source_path.name.endswith(".vmap"):
        raise ValueError("restrict_vmap.py requires source to be a .vmap")
    if not output_path.name.endswith(".vmap"):
        raise ValueError("restrict_vmap.py requires --output to end with .vmap")
    for path in restriction_paths:
        if not path.name.endswith((".vtable", ".vmap")):
            raise ValueError("restrict_vmap.py restriction inputs must be .vtable or .vmap")

    infos = [load_target_object_info(source_path), *(load_target_object_info(path) for path in restriction_paths)]
    require_shared_target_metadata(infos)
    source_table = read_mapped_vmap_table(source_path)
    out_frame = source_table.to_frame(copy=False)
    logger.info("restrict_vmap.py: loaded %s source variants from %s", len(out_frame), source_path)
    for path in restriction_paths:
        out_frame, row_count = restrict_frame_to_membership_chunks(out_frame, path)
        logger.info("restrict_vmap.py: loaded %s restriction variants from %s", row_count, path)
    write_vmap_table(output_path, VMapRowsTable.from_frame(out_frame, copy=False))
    write_metadata(output_path, dict(infos[0].raw_metadata))
    logger.info("restrict_vmap.py: wrote %s with %s restricted rows", output_path, len(out_frame))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
