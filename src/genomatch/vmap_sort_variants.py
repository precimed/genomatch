#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path

from ._cli_utils import run_cli
from .contig_utils import UNKNOWN_CONTIG, canonical_contig_from_label
from .tabular_rows import VMapRow, VariantRow
from .vtable_utils import (
    CANONICAL_CONTIG_RANK,
    ensure_parent_dir,
    iter_vmap_rows,
    iter_vtable_rows,
    load_metadata,
    require_contig_naming,
    validate_vmap_metadata,
    validate_vtable_metadata,
    variant_row_identity,
    write_metadata,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SortableRow:
    row_index: int
    rank: int
    pos_int: int
    row: VariantRow | VMapRow


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sort target rows in a .vtable/.vmap by declared contig order, then numeric position. "
            "The input object type and any .vmap provenance are preserved."
        )
    )
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable or .vmap")
    parser.add_argument(
        "--drop-duplicates",
        action="store_true",
        help="Drop duplicate target identities (chrom,pos,a1,a2), retaining first in stable sorted order",
    )
    return parser.parse_args()


def object_type_for_path(path: Path) -> str:
    if path.name.endswith(".vtable"):
        return "variant_table"
    if path.name.endswith(".vmap"):
        return "variant_map"
    raise ValueError(f"unsupported input: {path}")


def load_sort_metadata(path: Path) -> tuple[dict, str]:
    metadata = load_metadata(path)
    object_type = object_type_for_path(path)
    if object_type == "variant_table":
        validate_vtable_metadata(metadata)
        target_metadata = metadata
    else:
        validate_vmap_metadata(metadata)
        target_metadata = dict(metadata["target"])
    contig_naming = require_contig_naming(target_metadata, label="variant object")
    return metadata, contig_naming


def build_sort_rows(path: Path, *, object_type: str, contig_naming: str) -> list[SortableRow]:
    row_iter = iter_vtable_rows(path) if object_type == "variant_table" else iter_vmap_rows(path)
    out: list[SortableRow] = []
    for row_index, row in enumerate(row_iter):
        canonical = canonical_contig_from_label(row.chrom, contig_naming)
        if canonical is None or canonical == UNKNOWN_CONTIG:
            raise ValueError(
                f"variant object has contigs inconsistent with declared contig_naming={contig_naming!r}: {row.chrom!r}. "
                "Run normalize_contigs.py first"
            )
        out.append(SortableRow(row_index=row_index, rank=CANONICAL_CONTIG_RANK[canonical], pos_int=int(row.pos), row=row))
    return out


def write_sorted_rows(path: Path, *, rows: list[VariantRow | VMapRow], object_type: str) -> None:
    ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        if object_type == "variant_table":
            for row in rows:
                handle.write(f"{row.chrom}\t{row.pos}\t{row.id}\t{row.a1}\t{row.a2}\n")
            return
        for row in rows:
            if not isinstance(row, VMapRow):
                raise ValueError("internal error: expected VMapRow when writing variant_map output")
            handle.write(
                f"{row.chrom}\t{row.pos}\t{row.id}\t{row.a1}\t{row.a2}\t{row.source_shard}\t{row.source_index}\t{row.allele_op}\n"
            )


def sort_rows(rows: list[SortableRow], *, drop_duplicates: bool) -> tuple[list[VariantRow | VMapRow], int]:
    sorted_rows = sorted(rows, key=lambda item: (item.rank, item.pos_int, item.row_index))
    if not drop_duplicates:
        return [item.row for item in sorted_rows], 0
    dropped = 0
    out: list[VariantRow | VMapRow] = []
    seen: set[tuple[str, str, str, str]] = set()
    for item in sorted_rows:
        identity = variant_row_identity(item.row)
        if identity in seen:
            dropped += 1
            continue
        seen.add(identity)
        out.append(item.row)
    return out, dropped


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    logger.info("sort_variants.py: sorting %s -> %s", input_path, output_path)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")
    if object_type_for_path(input_path) != object_type_for_path(output_path):
        raise ValueError("sort_variants.py output suffix must match input suffix")

    object_type = object_type_for_path(input_path)
    metadata, contig_naming = load_sort_metadata(input_path)
    sortable_rows = build_sort_rows(input_path, object_type=object_type, contig_naming=contig_naming)
    rows, dropped = sort_rows(sortable_rows, drop_duplicates=args.drop_duplicates)

    temp_output = output_path.with_name(output_path.name + ".tmp")
    if temp_output.exists():
        temp_output.unlink()
    try:
        write_sorted_rows(temp_output, rows=rows, object_type=object_type)
        temp_output.replace(output_path)
    except Exception:
        if temp_output.exists():
            temp_output.unlink()
        raise

    write_metadata(output_path, dict(metadata))
    qc_path = output_path.with_name(output_path.name + ".qc.tsv")
    if qc_path.exists():
        qc_path.unlink()

    if args.drop_duplicates:
        logger.info("sort_variants.py: dropped %s duplicate target rows.", dropped)
    logger.info("sort_variants.py: retained %s rows in declared coordinate order.", len(rows))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
