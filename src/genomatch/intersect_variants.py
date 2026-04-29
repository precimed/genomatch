#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterator, List, Tuple

from ._cli_utils import run_cli
from .vtable_utils import (
    VariantRow,
    iter_vmap_rows,
    iter_vtable_rows,
    load_metadata,
    read_vmap,
    read_vtable,
    require_contig_naming,
    validate_vmap_metadata,
    validate_vtable_metadata,
    variant_row_identity,
    variant_row_from_vmap_row,
    variant_rows_from_vmap_rows,
    write_metadata,
    write_vtable,
)

logger = logging.getLogger(__name__)


def variant_key(row: VariantRow) -> Tuple[str, str, str, str]:
    return variant_row_identity(row)


def load_intersect_metadata(path: Path) -> dict:
    metadata = load_metadata(path)
    if path.name.endswith(".vtable"):
        validate_vtable_metadata(metadata)
        return metadata
    if path.name.endswith(".vmap"):
        validate_vmap_metadata(metadata)
        target_meta = dict(metadata["target"])
        target_meta["object_type"] = "variant_table"
        return target_meta
    raise ValueError(f"unsupported input: {path}")


def load_intersect_rows(path: Path) -> list[VariantRow]:
    if path.name.endswith(".vtable"):
        return read_vtable(path)
    if path.name.endswith(".vmap"):
        return variant_rows_from_vmap_rows(read_vmap(path))
    raise ValueError(f"unsupported input: {path}")


def iter_intersect_rows(path: Path) -> Iterator[VariantRow]:
    if path.name.endswith(".vtable"):
        yield from iter_vtable_rows(path, label="variant table")
        return
    elif path.name.endswith(".vmap"):
        for vmap_row in iter_vmap_rows(path):
            yield variant_row_from_vmap_row(vmap_row)
        return
    else:
        raise ValueError(f"unsupported input: {path}")


def stream_intersect_keys(path: Path, *, common: set[Tuple[str, str, str, str]]) -> tuple[set[Tuple[str, str, str, str]], int]:
    # Assumes: input metadata was prevalidated against all inputs.
    # Performs: streamed row PN/PV/SV sufficient for exact intersection membership.
    # Guarantees: returns keys from this input that are still candidates in `common`, plus row count.
    seen_in_this_input: set[Tuple[str, str, str, str]] = set()
    row_count = 0
    # PERF: line loop retained for streaming exact-key extraction; keeping row objects for later
    # inputs would scale memory with those inputs, which this function intentionally avoids.
    for row in iter_intersect_rows(path):
        row_count += 1
        key = variant_key(row)
        if key in common:
            seen_in_this_input.add(key)
    return seen_in_this_input, row_count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Intersect 2+ .vtable/.vmap inputs on exact variant rows.")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input .vtable/.vmap files")
    parser.add_argument("--output", required=True, help="Output .vtable")
    return parser.parse_args()


def print_loaded(path: Path, row_count: int) -> None:
    logger.info("intersect_variants.py: loaded %s variants from %s", row_count, path)


def print_remaining(path: Path, row_count: int) -> None:
    logger.info("intersect_variants.py: after intersecting %s, %s variants remain", path, row_count)


def format_metadata_summary(input_paths: list[Path], loaded_meta: list[dict]) -> str:
    lines = ["input metadata:"]
    for path, meta in zip(input_paths, loaded_meta):
        contig_naming = meta.get("contig_naming", "<missing>")
        lines.append(
            f"- {path}: genome_build={meta['genome_build']}, contig_naming={contig_naming}"
        )
    return "\n".join(lines)


def require_shared_metadata(input_paths: list[Path], loaded_meta: list[dict]) -> None:
    builds = {str(meta["genome_build"]) for meta in loaded_meta}
    metadata_summary = format_metadata_summary(input_paths, loaded_meta)
    contig_namings = set()
    for meta in loaded_meta:
        try:
            contig_namings.add(require_contig_naming(meta, label="variant object"))
        except ValueError as exc:
            raise ValueError(f"{exc}\n{metadata_summary}") from exc
    if len(builds) != 1:
        raise ValueError(f"all inputs must have the same genome_build\n{metadata_summary}")
    if len(contig_namings) != 1:
        raise ValueError(f"all inputs must have the same contig_naming\n{metadata_summary}")


def main() -> int:
    args = parse_args()
    if len(args.inputs) < 2:
        raise ValueError("intersect_variants.py requires at least two inputs")
    input_paths = [Path(item) for item in args.inputs]
    logger.info("intersect_variants.py: intersecting %s inputs -> %s", len(input_paths), args.output)
    loaded_meta = [load_intersect_metadata(path) for path in input_paths]
    require_shared_metadata(input_paths, loaded_meta)
    first_meta = loaded_meta[0]
    first_rows = load_intersect_rows(input_paths[0])
    print_loaded(input_paths[0], len(first_rows))
    common = set(variant_key(row) for row in first_rows)
    for path in input_paths[1:]:
        seen_in_this_input, row_count = stream_intersect_keys(path, common=common)
        print_loaded(path, row_count)
        common &= seen_in_this_input
        print_remaining(path, len(common))
    out_rows = [row for row in first_rows if variant_key(row) in common]
    write_vtable(Path(args.output), out_rows)
    write_metadata(Path(args.output), dict(first_meta))
    logger.info("intersect_variants.py: wrote %s with %s intersected rows", args.output, len(out_rows))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
