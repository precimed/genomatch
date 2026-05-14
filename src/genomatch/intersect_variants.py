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
    require_shared_target_metadata,
    restrict_frame_to_membership_chunks,
    restrict_frame_to_membership_frame,
    validate_loaded_row_count,
    vtable_metadata_from_target_info,
)
from .tabular_rows import VariantRowsTable
from .vtable_utils import metadata_with_variants_count, write_metadata, write_vtable_table

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


def missing_count_paths(infos: list[TargetObjectInfo]) -> list[Path]:
    return [info.path for info in infos if info.variants_count is None]


def intersect_declared_order(input_paths: list[Path], infos: list[TargetObjectInfo]):
    first_table = read_target_table(input_paths[0])
    retained_frame = first_table.to_frame(copy=False)
    validate_loaded_row_count(infos[0], len(retained_frame))
    print_loaded(input_paths[0], len(retained_frame))
    info_by_path = {info.path: info for info in infos}
    for path in input_paths[1:]:
        retained_frame, row_count = restrict_frame_to_membership_chunks(retained_frame, path)
        validate_loaded_row_count(info_by_path[path], row_count)
        print_loaded(path, row_count)
        print_remaining(path, len(retained_frame))
    return retained_frame


def intersect_count_order(input_paths: list[Path], infos: list[TargetObjectInfo]):
    declared_rank = {path: index for index, path in enumerate(input_paths)}
    ordered_infos = sorted(infos, key=lambda info: (info.variants_count, declared_rank[info.path]))
    driver_info = ordered_infos[0]
    driver_table = read_target_table(driver_info.path)
    candidate_frame = driver_table.to_frame(copy=False)
    validate_loaded_row_count(driver_info, len(candidate_frame))
    print_loaded(driver_info.path, len(candidate_frame))
    for info in ordered_infos[1:]:
        candidate_frame, row_count = restrict_frame_to_membership_chunks(candidate_frame, info.path)
        validate_loaded_row_count(info, row_count)
        print_loaded(info.path, row_count)
        print_remaining(info.path, len(candidate_frame))

    if driver_info is infos[0]:
        return candidate_frame

    first_info = infos[0]
    first_frame = read_target_table(input_paths[0]).to_frame(copy=False)
    validate_loaded_row_count(first_info, len(first_frame))
    return restrict_frame_to_membership_frame(first_frame, candidate_frame)


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
    missing_counts = missing_count_paths(infos)
    if missing_counts:
        logger.warning(
            "intersect_variants.py: variants_count missing from metadata for %s; using CLI input order",
            ", ".join(str(path) for path in missing_counts),
        )
        retained_frame = intersect_declared_order(input_paths, infos)
    else:
        retained_frame = intersect_count_order(input_paths, infos)

    write_vtable_table(output_path, VariantRowsTable.from_frame(retained_frame, copy=False))
    write_metadata(output_path, metadata_with_variants_count(vtable_metadata_from_target_info(infos[0]), len(retained_frame)))
    logger.info("intersect_variants.py: wrote %s with %s intersected rows", output_path, len(retained_frame))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
