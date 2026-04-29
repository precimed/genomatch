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
    classify_allele_operation,
    compose_allele_ops,
    iter_vmap_rows,
    load_variant_object,
    load_metadata,
    require_contig_naming,
    require_unique_target_rows,
    validate_vmap_metadata,
    variant_row_from_vmap_row,
    write_metadata,
    write_vmap,
)

logger = logging.getLogger(__name__)


def build_target_lookup(target_rows: Sequence[VariantRow]) -> Dict[Tuple[str, str], List[int]]:
    lookup: Dict[Tuple[str, str], List[int]] = {}
    for idx, row in enumerate(target_rows):
        lookup.setdefault((row.chrom, row.pos), []).append(idx)
    return lookup


def stream_source_matches(
    *,
    source_path: Path,
    target_rows: Sequence[VariantRow],
    target_lookup: Dict[Tuple[str, str], List[int]],
) -> dict[int, VMapRow]:
    # Assumes: source metadata and target metadata have already been checked for same build/naming.
    # Performs: streamed source row parsing plus coordinate lookup and swap-aware allele classification.
    # Guarantees: sparse target-index map preserving target order, source provenance, and first-match wins.
    output_map: dict[int, VMapRow] = {}
    # PERF: line loop retained for streaming source matching; target candidates are bounded by
    # target coordinate multiplicity, and retaining all source rows would dominate memory.
    for base_row in iter_vmap_rows(source_path):
        source_row = variant_row_from_vmap_row(base_row)
        for target_idx in target_lookup.get((source_row.chrom, source_row.pos), []):
            if target_idx in output_map:
                continue
            target_row = target_rows[target_idx]
            local_op = classify_allele_operation(
                source_row.a1,
                source_row.a2,
                target_row.a1,
                target_row.a2,
                allow_strand_flips=False,
            )
            if local_op == "missing" or base_row.source_index == -1:
                continue
            output_map[target_idx] = VMapRow(
                target_row.chrom,
                target_row.pos,
                target_row.id,
                target_row.a1,
                target_row.a2,
                base_row.source_shard,
                base_row.source_index,
                compose_allele_ops(base_row.allele_op, local_op),
            )
    return output_map


def materialize_matched_rows(target_rows: Sequence[VariantRow], output_map: dict[int, VMapRow]) -> list[VMapRow]:
    out: list[VMapRow] = []
    for idx, row in enumerate(target_rows):
        matched = output_map.get(idx)
        if matched is not None:
            out.append(matched)
        else:
            out.append(VMapRow(row.chrom, row.pos, row.id, row.a1, row.a2, MISSING_SOURCE_SHARD, -1, "missing"))
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
    target_loaded = load_variant_object(target_path)
    if target_loaded.object_type == "variant_map":
        logger.warning(
            "target .vmap provenance is ignored by match_vmap_to_target.py; "
            "use convert_vmap_to_target.py first to materialize a provenance-free .vtable and silence this warning."
        )
    target_meta = dict(target_loaded.target_metadata)
    target_contig_naming = require_contig_naming(target_meta, label="target variant object")
    target_rows = list(target_loaded.target_rows)
    if source_meta["target"]["genome_build"] != target_meta["genome_build"]:
        raise ValueError("source and target genome_build do not match")
    if source_contig_naming != target_contig_naming:
        raise ValueError("source and target contig_naming do not match")
    require_unique_target_rows(target_rows)

    target_lookup = build_target_lookup(target_rows)
    output_map = stream_source_matches(
        source_path=source_path,
        target_rows=target_rows,
        target_lookup=target_lookup,
    )
    out_rows = materialize_matched_rows(target_rows, output_map)
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
