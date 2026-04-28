#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import List, Tuple

from ._cli_utils import run_cli
from .vtable_utils import (
    MISSING_SOURCE_SHARD,
    SIGNED_INT_TOKEN_PATTERN,
    VALID_ALLELE_OPS,
    VariantRow,
    load_metadata,
    read_vmap,
    read_vtable,
    require_contig_naming,
    require_rows_match_contig_naming,
    validate_vmap_metadata,
    validate_vtable_metadata,
    variant_rows_from_vmap_rows,
    write_metadata,
    write_vtable,
)

logger = logging.getLogger(__name__)
INTERSECT_CHUNK_LINES = 100_000
POSITIVE_INT_TOKEN_RE = re.compile(r"\+?0*[1-9]\d*\Z")
SIGNED_INT_TOKEN_RE = re.compile(SIGNED_INT_TOKEN_PATTERN + r"\Z")
VALID_NONMISSING_ALLELE_OPS = VALID_ALLELE_OPS - {"missing"}


def variant_key(row: VariantRow) -> Tuple[str, str, str, str]:
    return (row.chrom, row.pos, row.a1, row.a2)


def load_intersect_metadata(path: Path) -> dict:
    metadata = load_metadata(path)
    if path.name.endswith(".vtable"):
        validate_vtable_metadata(metadata)
        return metadata
    if path.name.endswith(".vmap"):
        validate_vmap_metadata(metadata)
        target_meta = dict(metadata["target"])
        target_meta["object_type"] = "variant_table"
        return target_meta
    raise ValueError(f"unsupported input: {path}")


def load_intersect_rows(path: Path, metadata: dict) -> list[VariantRow]:
    if path.name.endswith(".vtable"):
        rows = read_vtable(path)
        require_rows_match_contig_naming(rows, require_contig_naming(metadata, label="variant table"), label="variant table")
        return rows
    if path.name.endswith(".vmap"):
        rows = variant_rows_from_vmap_rows(read_vmap(path))
        require_rows_match_contig_naming(rows, require_contig_naming(metadata, label="variant map target"), label="variant map target")
        return rows
    raise ValueError(f"unsupported input: {path}")


def parse_intersect_line(path: Path, line: str, metadata: dict) -> VariantRow:
    if path.name.endswith(".vtable"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 5:
            raise ValueError(f"invalid vtable row in {path}")
        chrom, pos, row_id, a1, a2 = parts
        label = "variant table"
    elif path.name.endswith(".vmap"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 8:
            raise ValueError(f"invalid vmap row in {path}")
        chrom, pos, row_id, a1, a2 = parts[:5]
        source_shard, source_index_raw, allele_op = parts[5:]
        source_index_token = source_index_raw.strip()
        if not SIGNED_INT_TOKEN_RE.fullmatch(source_index_token):
            raise ValueError(f"invalid vmap row in {path}")
        source_index = int(source_index_token)
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
        label = "variant map target"
    else:
        raise ValueError(f"unsupported input: {path}")
    if not POSITIVE_INT_TOKEN_RE.fullmatch(pos):
        raise ValueError(f"{label} row has invalid pos: {pos!r}")
    a1 = a1.strip().upper()
    a2 = a2.strip().upper()
    if re.fullmatch(r"[ACGT]+", a1) is None:
        raise ValueError(f"invalid allele code in {label}: {parts[3]!r}")
    if re.fullmatch(r"[ACGT]+", a2) is None:
        raise ValueError(f"invalid allele code in {label}: {parts[4]!r}")
    row = VariantRow(chrom, pos, row_id, a1, a2)
    require_rows_match_contig_naming([row], require_contig_naming(metadata, label=label), label=label)
    return row


def stream_intersect_keys(path: Path, metadata: dict, *, common: set[Tuple[str, str, str, str]]) -> tuple[set[Tuple[str, str, str, str]], int]:
    # Assumes: metadata was prevalidated against all inputs.
    # Performs: streamed row PN/PV/SV sufficient for exact intersection membership.
    # Guarantees: returns keys from this input that are still candidates in `common`, plus row count.
    seen_in_this_input: set[Tuple[str, str, str, str]] = set()
    row_count = 0
    with open(path, "r", encoding="utf-8", newline="") as handle:
        chunk: list[str] = []
        for line in handle:
            if not line.strip():
                continue
            chunk.append(line)
            if len(chunk) >= INTERSECT_CHUNK_LINES:
                row_count += process_intersect_chunk(path, metadata, chunk, common, seen_in_this_input)
                chunk = []
        if chunk:
            row_count += process_intersect_chunk(path, metadata, chunk, common, seen_in_this_input)
    return seen_in_this_input, row_count


def process_intersect_chunk(
    path: Path,
    metadata: dict,
    chunk: list[str],
    common: set[Tuple[str, str, str, str]],
    seen_in_this_input: set[Tuple[str, str, str, str]],
) -> int:
    # PERF: line loop retained for streaming exact-key extraction; keeping row objects for later
    # inputs would scale memory with those inputs, which this function intentionally avoids.
    for line in chunk:
        key = variant_key(parse_intersect_line(path, line, metadata))
        if key in common:
            seen_in_this_input.add(key)
    return len(chunk)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Intersect 2+ .vtable/.vmap inputs on exact variant rows.")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input .vtable/.vmap files")
    parser.add_argument("--output", required=True, help="Output .vtable")
    return parser.parse_args()


def print_loaded(path: Path, row_count: int) -> None:
    logger.info("intersect_variants.py: loaded %s variants from %s", row_count, path)


def print_remaining(path: Path, row_count: int) -> None:
    logger.info("intersect_variants.py: after intersecting %s, %s variants remain", path, row_count)


def format_metadata_summary(input_paths: list[Path], loaded_meta: list[dict]) -> str:
    lines = ["input metadata:"]
    for path, meta in zip(input_paths, loaded_meta):
        contig_naming = meta.get("contig_naming", "<missing>")
        lines.append(
            f"- {path}: genome_build={meta['genome_build']}, contig_naming={contig_naming}"
        )
    return "\n".join(lines)


def require_shared_metadata(input_paths: list[Path], loaded_meta: list[dict]) -> None:
    builds = {str(meta["genome_build"]) for meta in loaded_meta}
    metadata_summary = format_metadata_summary(input_paths, loaded_meta)
    contig_namings = set()
    for meta in loaded_meta:
        try:
            contig_namings.add(require_contig_naming(meta, label="variant object"))
        except ValueError as exc:
            raise ValueError(f"{exc}\n{metadata_summary}") from exc
    if len(builds) != 1:
        raise ValueError(f"all inputs must have the same genome_build\n{metadata_summary}")
    if len(contig_namings) != 1:
        raise ValueError(f"all inputs must have the same contig_naming\n{metadata_summary}")


def main() -> int:
    args = parse_args()
    if len(args.inputs) < 2:
        raise ValueError("intersect_variants.py requires at least two inputs")
    input_paths = [Path(item) for item in args.inputs]
    logger.info("intersect_variants.py: intersecting %s inputs -> %s", len(input_paths), args.output)
    loaded_meta = [load_intersect_metadata(path) for path in input_paths]
    require_shared_metadata(input_paths, loaded_meta)
    first_meta = loaded_meta[0]
    first_rows = load_intersect_rows(input_paths[0], loaded_meta[0])
    print_loaded(input_paths[0], len(first_rows))
    common = set(variant_key(row) for row in first_rows)
    for path, metadata in zip(input_paths[1:], loaded_meta[1:]):
        seen_in_this_input, row_count = stream_intersect_keys(path, metadata, common=common)
        print_loaded(path, row_count)
        common &= seen_in_this_input
        print_remaining(path, len(common))
    out_rows = [row for row in first_rows if variant_key(row) in common]
    write_vtable(Path(args.output), out_rows)
    write_metadata(Path(args.output), dict(first_meta))
    logger.info("intersect_variants.py: wrote %s with %s intersected rows", args.output, len(out_rows))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
