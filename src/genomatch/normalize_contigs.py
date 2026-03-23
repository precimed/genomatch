#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ._cli_utils import run_cli
from .contig_cleanup_utils import (
    load_target_variant_object,
    normalize_target_rows,
    normalized_output_metadata,
    write_variant_object_like_input,
)
from .vtable_utils import SUPPORTED_CONTIG_NAMINGS, duplicate_target_row_keys, target_row_key, write_vmap_status_qc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Normalize contig labels in a .vtable or .vmap.")
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable or .vmap")
    parser.add_argument("--to", required=True, choices=sorted(SUPPORTED_CONTIG_NAMINGS))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")
    loaded = load_target_variant_object(input_path)
    input_rows = loaded.base_vmap_rows if loaded.base_vmap_rows is not None else loaded.target_rows
    result = normalize_target_rows(input_rows, args.to, genome_build=str(loaded.target_metadata.get("genome_build")))
    output_rows = list(result.rows)
    duplicate_target_count = 0
    qc_path = output_path.with_name(output_path.name + ".qc.tsv")
    if loaded.base_vmap_rows is not None:
        duplicate_keys = duplicate_target_row_keys(output_rows)
        duplicate_target_count = sum(1 for row in output_rows if target_row_key(row) in duplicate_keys)
        if duplicate_keys:
            qc_rows = [
                (row.source_shard, row.source_index, row.id, "duplicate_target")
                for row in output_rows
                if target_row_key(row) in duplicate_keys
            ]
            output_rows = [row for row in output_rows if target_row_key(row) not in duplicate_keys]
            write_vmap_status_qc(qc_path, qc_rows)
        elif qc_path.exists():
            qc_path.unlink()
    write_variant_object_like_input(
        loaded,
        output_path,
        output_rows,
        normalized_output_metadata(loaded, args.to),
        preserve_qc=True,
    )
    print(
        f"normalize_contigs.py: normalized {result.normalized_count} rows to {args.to}; "
        f"dropped {result.unknown_count} unresolved rows and {duplicate_target_count} duplicate target rows.",
        file=sys.stderr,
    )
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
