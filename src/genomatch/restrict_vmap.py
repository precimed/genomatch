#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ._cli_utils import run_cli
from .exact_set_utils import (
    TargetObjectInfo,
    load_target_object_info,
    read_target_table,
    read_mapped_vmap_table,
    require_shared_target_metadata,
    restrict_frame_to_membership_chunks,
    restrict_frame_to_membership_frame,
    validate_loaded_row_count,
)
from .tabular_rows import VMapRowsTable
from .vtable_utils import metadata_with_variants_count, write_metadata, write_vmap_table

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Restrict a mapped-only .vmap by exact membership in one or more .vtable/.vmap inputs."
    )
    parser.add_argument("source", help="Source mapped-only .vmap")
    parser.add_argument("restrictions", nargs="+", help="Restriction .vtable/.vmap files")
    parser.add_argument("--output", required=True, help="Output .vmap")
    return parser.parse_args()


def missing_count_paths(infos: list[TargetObjectInfo]) -> list[Path]:
    return [info.path for info in infos if info.variants_count is None]


def restrict_declared_order(source_path: Path, restriction_paths: list[Path], infos: list[TargetObjectInfo]):
    source_table = read_mapped_vmap_table(source_path)
    out_frame = source_table.to_frame(copy=False)
    validate_loaded_row_count(infos[0], len(out_frame))
    logger.info("restrict_vmap.py: loaded %s source variants from %s", len(out_frame), source_path)
    info_by_path = {info.path: info for info in infos}
    for path in restriction_paths:
        out_frame, row_count = restrict_frame_to_membership_chunks(out_frame, path)
        validate_loaded_row_count(info_by_path[path], row_count)
        logger.info("restrict_vmap.py: loaded %s restriction variants from %s", row_count, path)
    return out_frame


def read_target_frame_for_info(info: TargetObjectInfo):
    frame = read_target_table(info.path).to_frame(copy=False)
    validate_loaded_row_count(info, len(frame))
    return frame


def restrict_count_order(source_path: Path, restriction_paths: list[Path], infos: list[TargetObjectInfo]):
    source_info = infos[0]
    ordered_items = sorted(
        enumerate(infos),
        key=lambda item: (item[1].variants_count, item[0]),
    )
    candidate_frame = None
    out_frame = None
    for index, info in ordered_items:
        if index == 0:
            source_table = read_mapped_vmap_table(source_path)
            source_frame = source_table.to_frame(copy=False)
            validate_loaded_row_count(source_info, len(source_frame))
            logger.info("restrict_vmap.py: loaded %s source variants from %s", len(source_frame), source_path)
            if candidate_frame is None:
                out_frame = source_frame
            else:
                out_frame = restrict_frame_to_membership_frame(source_frame, candidate_frame)
            continue

        if out_frame is None:
            if candidate_frame is None:
                candidate_frame = read_target_frame_for_info(info)
                row_count = len(candidate_frame)
            else:
                candidate_frame, row_count = restrict_frame_to_membership_chunks(candidate_frame, info.path)
                validate_loaded_row_count(info, row_count)
            logger.info("restrict_vmap.py: loaded %s restriction variants from %s", row_count, info.path)
        else:
            out_frame, row_count = restrict_frame_to_membership_chunks(out_frame, info.path)
            validate_loaded_row_count(info, row_count)
            logger.info("restrict_vmap.py: loaded %s restriction variants from %s", row_count, info.path)

    assert out_frame is not None
    return out_frame


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
    missing_counts = missing_count_paths(infos)
    if missing_counts:
        logger.warning(
            "restrict_vmap.py: variants_count missing from metadata for %s; using CLI input order",
            ", ".join(str(path) for path in missing_counts),
        )
        out_frame = restrict_declared_order(source_path, restriction_paths, infos)
    else:
        out_frame = restrict_count_order(source_path, restriction_paths, infos)
    write_vmap_table(output_path, VMapRowsTable.from_frame(out_frame, copy=False))
    write_metadata(output_path, metadata_with_variants_count(infos[0].raw_metadata, len(out_frame)))
    logger.info("restrict_vmap.py: wrote %s with %s restricted rows", output_path, len(out_frame))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
