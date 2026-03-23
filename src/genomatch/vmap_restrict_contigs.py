#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .contig_cleanup_utils import (
    build_target_restriction_selection_for_contigs,
    load_target_variant_object,
    restriction_output_metadata,
    write_variant_object_like_input,
)
from .vtable_utils import parse_chr2use


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Filter retained target rows by --chr2use / --contigs.")
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable or .vmap")
    parser.add_argument("--chr2use", "--contigs", dest="chr2use", help="Comma-separated chromosome list or ranges")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")
    loaded = load_target_variant_object(input_path)
    input_rows = loaded.base_vmap_rows if loaded.base_vmap_rows is not None else loaded.target_rows
    target_meta = dict(loaded.target_metadata)
    restriction_output_metadata(loaded)
    allowed_contigs, explicit = parse_chr2use(args.chr2use)
    selection_rows = build_target_restriction_selection_for_contigs(
        input_rows,
        allowed_canonical=set(allowed_contigs) if explicit else None,
        contig_naming=str(target_meta["contig_naming"]),
    )
    out_rows = list(selection_rows)
    out_meta = restriction_output_metadata(loaded)
    write_variant_object_like_input(loaded, output_path, out_rows, out_meta)
    dropped_chr2use = len(input_rows) - len(selection_rows)
    message = f"restrict_contigs.py: retained {len(selection_rows)} rows"
    if explicit:
        message += f" and dropped {dropped_chr2use} rows outside --chr2use"
    else:
        message += f" and dropped {dropped_chr2use} rows"
    print(message + ".", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
