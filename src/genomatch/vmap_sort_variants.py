#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ._cli_utils import run_cli
from .contig_cleanup_utils import load_target_variant_object, write_variant_object_like_input
from .vtable_utils import require_contig_naming, sort_target_rows_by_declared_coordinate

logger = logging.getLogger(__name__)


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
    logger.info("sort_variants.py: sorting %s -> %s", input_path, output_path)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")

    loaded = load_target_variant_object(input_path)
    contig_naming = require_contig_naming(loaded.target_metadata, label="variant object")
    input_rows = loaded.base_vmap_rows if loaded.base_vmap_rows is not None else loaded.target_rows
    out_rows = sort_target_rows_by_declared_coordinate(input_rows, contig_naming, label="variant object")
    write_variant_object_like_input(loaded, output_path, out_rows, dict(loaded.raw_metadata))
    logger.info("sort_variants.py: retained %s rows in declared coordinate order.", len(out_rows))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
