#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from contig_cleanup_utils import load_target_variant_object, write_variant_object_like_input
from vtable_utils import require_contig_naming, sort_target_rows_by_declared_coordinate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sort target rows in a .vtable/.vmap by declared contig order, then numeric position. "
            "The input object type and any .vmap provenance are preserved."
        )
    )
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable or .vmap")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")

    loaded = load_target_variant_object(input_path)
    contig_naming = require_contig_naming(loaded.target_metadata, label="variant object")
    input_rows = loaded.base_vmap_rows if loaded.base_vmap_rows is not None else loaded.target_rows
    out_rows = sort_target_rows_by_declared_coordinate(input_rows, contig_naming, label="variant object")
    write_variant_object_like_input(loaded, output_path, out_rows, dict(loaded.raw_metadata))
    print(f"sort_variants.py: retained {len(out_rows)} rows in declared coordinate order.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
