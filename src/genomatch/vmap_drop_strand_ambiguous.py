#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ._cli_utils import run_cli
from .contig_cleanup_utils import load_target_variant_object, write_variant_object_like_input
from .vtable_utils import (
    complement_allele,
    require_contig_naming,
    require_rows_match_contig_naming,
    validate_allele_value,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Drop target-side strand-ambiguous biallelic rows.")
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable or .vmap")
    return parser.parse_args()


def is_strand_ambiguous(a1: str, a2: str) -> bool:
    validate_allele_value(a1, label="drop_strand_ambiguous.py input")
    validate_allele_value(a2, label="drop_strand_ambiguous.py input")
    return complement_allele(a1) == a2


def build_selection_rows(rows):
    return [row for row in rows if not is_strand_ambiguous(row.a1, row.a2)]


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    logger.info("drop_strand_ambiguous.py: filtering %s -> %s", input_path, output_path)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")

    loaded = load_target_variant_object(input_path)
    target_meta = dict(loaded.target_metadata)
    contig_naming = require_contig_naming(target_meta, label="variant object")
    input_rows = loaded.base_vmap_rows if loaded.base_vmap_rows is not None else loaded.target_rows
    require_rows_match_contig_naming(input_rows, contig_naming, label="variant object")

    selection_rows = build_selection_rows(input_rows)
    out_rows = list(selection_rows)
    out_meta = dict(loaded.raw_metadata)
    write_variant_object_like_input(loaded, output_path, out_rows, out_meta)
    dropped = len(input_rows) - len(selection_rows)
    logger.info(
        "drop_strand_ambiguous.py: retained %s rows and dropped %s strand-ambiguous rows.",
        len(selection_rows),
        dropped,
    )
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
