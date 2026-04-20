#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

from ._cli_utils import run_cli
from .vtable_utils import (
    VariantRow,
    load_metadata,
    read_vmap,
    read_vtable,
    require_contig_naming,
    require_rows_match_contig_naming,
    validate_vmap_metadata,
    validate_vtable_metadata,
    variant_rows_from_vmap_rows,
    write_metadata,
    write_vtable,
)


def variant_key(row: VariantRow) -> Tuple[str, str, str, str]:
    return (row.chrom, row.pos, row.a1, row.a2)


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


def load_intersect_rows(path: Path, metadata: dict) -> list[VariantRow]:
    if path.name.endswith(".vtable"):
        rows = read_vtable(path)
        require_rows_match_contig_naming(rows, require_contig_naming(metadata, label="variant table"), label="variant table")
        return rows
    if path.name.endswith(".vmap"):
        rows = variant_rows_from_vmap_rows(read_vmap(path))
        require_rows_match_contig_naming(rows, require_contig_naming(metadata, label="variant map target"), label="variant map target")
        return rows
    raise ValueError(f"unsupported input: {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Intersect 2+ .vtable/.vmap inputs on exact variant rows.")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input .vtable/.vmap files")
    parser.add_argument("--output", required=True, help="Output .vtable")
    return parser.parse_args()


def print_loaded(path: Path, row_count: int) -> None:
    print(f"intersect_variants.py: loaded {row_count} variants from {path}", file=sys.stderr)


def print_remaining(path: Path, row_count: int) -> None:
    print(f"intersect_variants.py: after intersecting {path}, {row_count} variants remain", file=sys.stderr)


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
    loaded_meta = [load_intersect_metadata(path) for path in input_paths]
    require_shared_metadata(input_paths, loaded_meta)
    first_meta = loaded_meta[0]
    loaded_rows = []
    for path, metadata in zip(input_paths, loaded_meta):
        rows = load_intersect_rows(path, metadata)
        loaded_rows.append(rows)
        print_loaded(path, len(rows))
    first_rows = loaded_rows[0]
    common = set(variant_key(row) for row in first_rows)
    for path, rows in zip(input_paths[1:], loaded_rows[1:]):
        common &= set(variant_key(row) for row in rows)
        print_remaining(path, len(common))
    out_rows = [row for row in first_rows if variant_key(row) in common]
    write_vtable(Path(args.output), out_rows)
    write_metadata(Path(args.output), dict(first_meta))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
