#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from ._cli_utils import run_cli
from .vtable_utils import (
    MISSING_SOURCE_SHARD,
    VALID_ALLELE_OPS,
    VMapRow,
    VariantRow,
    canonical_contig_from_row,
    classify_allele_operation,
    compose_allele_ops,
    load_variant_object,
    load_metadata,
    require_contig_naming,
    require_rows_match_contig_naming,
    validate_vmap_metadata,
    write_metadata,
    write_vmap,
)

logger = logging.getLogger(__name__)
MATCH_CHUNK_LINES = 100_000
POSITIVE_INT_TOKEN_RE = re.compile(r"\+?0*[1-9]\d*\Z")
VALID_NONMISSING_ALLELE_OPS = VALID_ALLELE_OPS - {"missing"}


def target_key(row: VariantRow, contig_naming: str, *, label: str) -> Tuple[str, str, str, str]:
    return (canonical_contig_from_row(row, contig_naming, label=label), row.pos, row.a1, row.a2)


def detect_target_duplicates(rows: Sequence[VariantRow], contig_naming: str, *, label: str) -> None:
    seen = set()
    for row in rows:
        key = target_key(row, contig_naming, label=label)
        if key in seen:
            raise ValueError("duplicate chrom:pos:a1:a2 in target rows")
        seen.add(key)


def build_target_lookup(
    target_rows: Sequence[VariantRow],
    target_contig_naming: str,
) -> Dict[Tuple[str, str], List[int]]:
    lookup: Dict[Tuple[str, str], List[int]] = {}
    for idx, row in enumerate(target_rows):
        chrom = canonical_contig_from_row(row, target_contig_naming, label="target variant object")
        lookup.setdefault((chrom, row.pos), []).append(idx)
    return lookup


def parse_source_vmap_line(path: Path, line: str) -> VMapRow:
    parts = line.rstrip("\n").split("\t")
    if len(parts) != 8:
        raise ValueError(f"invalid vmap row in {path}")
    chrom, pos, row_id, a1, a2, source_shard, source_index_raw, allele_op = parts
    if not POSITIVE_INT_TOKEN_RE.fullmatch(pos):
        raise ValueError(f"vmap row has invalid pos: {pos!r}")
    source_index_token = source_index_raw.strip()
    if not re.fullmatch(r"[+-]?\d+", source_index_token):
        raise ValueError(f"invalid vmap row in {path}")
    source_index = int(source_index_token)
    a1 = a1.strip().upper()
    a2 = a2.strip().upper()
    if re.fullmatch(r"[ACGT]+", a1) is None:
        raise ValueError(f"invalid allele code in vmap: {parts[3]!r}")
    if re.fullmatch(r"[ACGT]+", a2) is None:
        raise ValueError(f"invalid allele code in vmap: {parts[4]!r}")
    if source_index == -1:
        if source_shard != MISSING_SOURCE_SHARD:
            raise ValueError("vmap row with source_index=-1 must have source_shard='.'")
        if allele_op != "missing":
            raise ValueError("vmap row with source_index=-1 must have allele_op=missing")
    elif source_index < 0:
        raise ValueError("vmap row source_index out of range")
    elif not source_shard:
        raise ValueError("vmap row with source_index>=0 must define source_shard")
    elif allele_op not in VALID_NONMISSING_ALLELE_OPS:
        raise ValueError(f"invalid allele_op in vmap row: {allele_op!r}")
    return VMapRow(chrom, pos, row_id, a1, a2, source_shard, source_index, allele_op)


def stream_source_matches(
    *,
    source_path: Path,
    source_contig_naming: str,
    target_rows: Sequence[VariantRow],
    target_lookup: Dict[Tuple[str, str], List[int]],
) -> dict[int, VMapRow]:
    # Assumes: source metadata and target metadata have already been checked for same build/naming.
    # Performs: streamed source row parsing plus coordinate lookup and swap-aware allele classification.
    # Guarantees: sparse target-index map preserving target order, source provenance, and first-match wins.
    output_map: dict[int, VMapRow] = {}
    with open(source_path, "r", encoding="utf-8", newline="") as handle:
        chunk: list[str] = []
        for line in handle:
            if not line.strip():
                continue
            chunk.append(line)
            if len(chunk) >= MATCH_CHUNK_LINES:
                process_source_match_chunk(
                    source_path=source_path,
                    source_contig_naming=source_contig_naming,
                    target_rows=target_rows,
                    target_lookup=target_lookup,
                    output_map=output_map,
                    chunk=chunk,
                )
                chunk = []
        if chunk:
            process_source_match_chunk(
                source_path=source_path,
                source_contig_naming=source_contig_naming,
                target_rows=target_rows,
                target_lookup=target_lookup,
                output_map=output_map,
                chunk=chunk,
            )
    return output_map


def process_source_match_chunk(
    *,
    source_path: Path,
    source_contig_naming: str,
    target_rows: Sequence[VariantRow],
    target_lookup: Dict[Tuple[str, str], List[int]],
    output_map: dict[int, VMapRow],
    chunk: list[str],
) -> None:
    # PERF: line loop retained for streaming source matching; target candidates are bounded by
    # target coordinate multiplicity, and retaining all source rows would dominate memory.
    for line in chunk:
        base_row = parse_source_vmap_line(source_path, line)
        source_row = VariantRow(base_row.chrom, base_row.pos, base_row.id, base_row.a1, base_row.a2)
        chrom = canonical_contig_from_row(source_row, source_contig_naming, label="source variant map target")
        for target_idx in target_lookup.get((chrom, source_row.pos), []):
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
    require_rows_match_contig_naming(target_rows, target_contig_naming, label="target variant object")
    if source_meta["target"]["genome_build"] != target_meta["genome_build"]:
        raise ValueError("source and target genome_build do not match")
    if source_contig_naming != target_contig_naming:
        raise ValueError("source and target contig_naming do not match")
    detect_target_duplicates(target_rows, target_contig_naming, label="target variant object")

    target_lookup = build_target_lookup(target_rows, target_contig_naming)
    output_map = stream_source_matches(
        source_path=source_path,
        source_contig_naming=source_contig_naming,
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
