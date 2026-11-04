#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ._cli_utils import run_cli
from .contig_cleanup_utils import (
    normalize_target_table,
    normalized_output_metadata,
)
from .vtable_utils import (
    SUPPORTED_CONTIG_NAMINGS,
    VMapRowsTable,
    VariantRowsTable,
    duplicate_target_rows_mask_table,
    load_variant_object_tables,
    metadata_with_variants_count,
    variant_row_qc_rows,
    write_metadata,
    write_variant_status_qc,
    write_vmap_table,
    write_vtable_table,
)

logger = logging.getLogger(__name__)

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
    logger.info("normalize_contigs.py: normalizing %s -> %s (to=%s)", input_path, output_path, args.to)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")
    loaded = load_variant_object_tables(input_path)
    genome_build = str(loaded.target_metadata.get("genome_build"))
    is_vmap = loaded.base_vmap_table is not None
    source_frame = (
        loaded.base_vmap_table.to_frame(copy=False) if is_vmap else loaded.target_rows_table.to_frame(copy=False)
    )
    # Assumes: source_frame has canonical vmap/vtable columns; args.to is a valid contig naming token.
    # Guarantees: result.frame has contig-normalized target rows; result.duplicate_mask flags chrom:pos:a1:a2 duplicates.
    result = normalize_target_table(source_frame, args.to, genome_build=genome_build)
    out_frame = result.frame
    qc_path = output_path.with_name(output_path.name + ".qc.tsv")
    duplicate_mask = duplicate_target_rows_mask_table(out_frame)
    duplicate_target_count = int(duplicate_mask.sum())
    if duplicate_target_count > 0:
        # Normalization may discard unresolved rows and reset the structural index.
        # Match QC by stored identifiers, never by positions in the filtered frame.
        dup_frame = out_frame.loc[duplicate_mask]
        qc_rows = variant_row_qc_rows(dup_frame, is_vmap=is_vmap)
        write_variant_status_qc(qc_path, qc_rows.assign(status="duplicate_target"), is_vmap=is_vmap)
        out_frame = out_frame.loc[~duplicate_mask].reset_index(drop=True)
    elif qc_path.exists():
        qc_path.unlink()
    if is_vmap:
        write_vmap_table(output_path, VMapRowsTable.from_frame(out_frame, copy=False), assume_validated=True)
    else:
        write_vtable_table(output_path, VariantRowsTable.from_frame(out_frame, copy=False), assume_validated=True)
    write_metadata(output_path, metadata_with_variants_count(normalized_output_metadata(loaded, args.to), len(out_frame)))
    logger.info(
        "normalize_contigs.py: normalized %s rows to %s; dropped %s unresolved rows and %s duplicate target rows.",
        result.normalized_count,
        args.to,
        result.unknown_count,
        duplicate_target_count,
    )
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
