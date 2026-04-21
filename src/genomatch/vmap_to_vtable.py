#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from ._cli_utils import run_cli
from .vtable_utils import (
    load_metadata,
    make_vtable_metadata,
    read_vmap,
    require_contig_naming,
    require_rows_match_contig_naming,
    validate_vmap_metadata,
    variant_rows_from_vmap_rows,
    write_metadata,
    write_vtable,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Materialize the target side of a .vmap as a standalone .vtable.")
    parser.add_argument("--source", required=True, help="Input .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_path = Path(args.source)
    output_path = Path(args.output)
    if not source_path.exists():
        raise ValueError(f"source .vmap not found: {source_path}")
    metadata = load_metadata(source_path)
    validate_vmap_metadata(metadata)
    vmap_rows = read_vmap(source_path)
    target_meta = dict(metadata["target"])
    require_rows_match_contig_naming(vmap_rows, require_contig_naming(target_meta, label="variant map target"), label="variant map target")
    rows = variant_rows_from_vmap_rows(vmap_rows)
    write_vtable(output_path, rows)
    write_metadata(
        output_path,
        make_vtable_metadata(
            genome_build=target_meta["genome_build"],
            contig_naming=target_meta["contig_naming"],
            provenance={"created_by": "convert_vmap_to_target.py", "derived_from": str(source_path)},
        ),
    )
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
