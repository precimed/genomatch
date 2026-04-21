#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from ._cli_utils import run_cli
from .vtable_utils import (
    MISSING_SOURCE_SHARD,
    VMapRow,
    VariantRow,
    canonical_contig_from_row,
    classify_allele_operation,
    compose_allele_ops,
    load_variant_object,
    load_metadata,
    read_vmap,
    read_vtable,
    require_contig_naming,
    require_rows_match_contig_naming,
    validate_vmap_metadata,
    validate_vtable_metadata,
    write_metadata,
    write_vmap,
)

logger = logging.getLogger(__name__)


def target_key(row: VariantRow, contig_naming: str, *, label: str) -> Tuple[str, str, str, str]:
    return (canonical_contig_from_row(row, contig_naming, label=label), row.pos, row.a1, row.a2)


def detect_target_duplicates(rows: Sequence[VariantRow], contig_naming: str, *, label: str) -> None:
    seen = set()
    for row in rows:
        key = target_key(row, contig_naming, label=label)
        if key in seen:
            raise ValueError("duplicate chrom:pos:a1:a2 in target rows")
        seen.add(key)


def build_lookup(
    source_rows: Sequence[VariantRow],
    source_vmap_rows: Sequence[VMapRow],
    contig_naming: str,
) -> Dict[Tuple[str, str], List[Tuple[VMapRow, VariantRow]]]:
    lookup: Dict[Tuple[str, str], List[Tuple[VMapRow, VariantRow]]] = {}
    for base_row, source_row in zip(source_vmap_rows, source_rows):
        chrom = canonical_contig_from_row(source_row, contig_naming, label="source variant map target")
        lookup.setdefault((chrom, source_row.pos), []).append((base_row, source_row))
    return lookup


def match_rows(
    source_rows: Sequence[VariantRow],
    source_vmap_rows: Sequence[VMapRow],
    source_contig_naming: str,
    target_rows: Sequence[VariantRow],
    target_contig_naming: str,
) -> List[VMapRow]:
    lookup = build_lookup(source_rows, source_vmap_rows, source_contig_naming)
    out: List[VMapRow] = []
    for row in target_rows:
        chrom = canonical_contig_from_row(row, target_contig_naming, label="target variant table")
        matched_row = VMapRow(row.chrom, row.pos, row.id, row.a1, row.a2, MISSING_SOURCE_SHARD, -1, "missing")
        for base_row, source_row in lookup.get((chrom, row.pos), []):
            local_op = classify_allele_operation(
                source_row.a1,
                source_row.a2,
                row.a1,
                row.a2,
                allow_strand_flips=False,
            )
            if local_op == "missing":
                continue
            if base_row.source_index == -1:
                continue
            else:
                matched_row = VMapRow(
                    row.chrom,
                    row.pos,
                    row.id,
                    row.a1,
                    row.a2,
                    base_row.source_shard,
                    base_row.source_index,
                    compose_allele_ops(base_row.allele_op, local_op),
                )
            break
        out.append(matched_row)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a .vmap from a source .vmap to target rows from a .vtable or .vmap.")
    parser.add_argument("--source", required=True, help="Source .vmap")
    parser.add_argument("--target", required=True, help="Target .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vmap")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_path = Path(args.source)
    target_path = Path(args.target)
    output_path = Path(args.output)
    logger.info("match_vmap_to_target.py: matching %s to %s -> %s", source_path, target_path, output_path)
    if "@" in args.source or "@" in args.target or "@" in args.output:
        raise ValueError("match_vmap_to_target.py does not accept '@' paths")
    if not source_path.exists():
        raise ValueError(f"source file not found: {source_path}")
    if not target_path.exists():
        raise ValueError(f"target file not found: {target_path}")
    if not source_path.name.endswith(".vmap"):
        raise ValueError("match_vmap_to_target.py requires --source to be a .vmap")
    if not target_path.name.endswith((".vtable", ".vmap")):
        raise ValueError("match_vmap_to_target.py requires --target to be a .vtable or .vmap")

    source_meta = load_metadata(source_path)
    validate_vmap_metadata(source_meta)
    source_contig_naming = require_contig_naming(dict(source_meta["target"]), label="source variant map target")
    source_vmap_rows = read_vmap(source_path)
    source_rows = [VariantRow(row.chrom, row.pos, row.id, row.a1, row.a2) for row in source_vmap_rows]
    require_rows_match_contig_naming(source_rows, source_contig_naming, label="source variant map target")

    target_loaded = load_variant_object(target_path)
    if target_loaded.object_type == "variant_map":
        logger.warning(
            "target .vmap provenance is ignored by match_vmap_to_target.py; "
            "use convert_vmap_to_target.py first to materialize a provenance-free .vtable and silence this warning."
        )
    target_meta = dict(target_loaded.target_metadata)
    target_contig_naming = require_contig_naming(target_meta, label="target variant object")
    target_rows = list(target_loaded.target_rows)
    require_rows_match_contig_naming(target_rows, target_contig_naming, label="target variant object")
    if source_meta["target"]["genome_build"] != target_meta["genome_build"]:
        raise ValueError("source and target genome_build do not match")
    if source_contig_naming != target_contig_naming:
        raise ValueError("source and target contig_naming do not match")
    detect_target_duplicates(target_rows, target_contig_naming, label="target variant object")

    out_rows = match_rows(source_rows, source_vmap_rows, source_contig_naming, target_rows, target_contig_naming)
    write_vmap(output_path, out_rows)
    write_metadata(
        output_path,
        {
            "object_type": "variant_map",
            "target": {
                "genome_build": target_meta["genome_build"],
                "contig_naming": target_meta["contig_naming"],
            },
        },
    )
    missing = sum(1 for row in out_rows if row.source_index == -1)
    logger.info(
        "match_vmap_to_target.py: matched %s/%s target rows; missing %s.",
        len(out_rows) - missing,
        len(out_rows),
        missing,
    )
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
