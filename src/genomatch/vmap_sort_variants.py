#!/usr/bin/env python3
from __future__ import annotations

import argparse
import heapq
import logging
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO

from ._cli_utils import run_cli
from .contig_utils import UNKNOWN_CONTIG, canonical_contig_from_label
from .vtable_utils import (
    CANONICAL_CONTIG_RANK,
    MISSING_SOURCE_SHARD,
    SIGNED_INT_TOKEN_PATTERN,
    VALID_ALLELE_OPS,
    ensure_parent_dir,
    load_metadata,
    require_contig_naming,
    validate_vmap_metadata,
    validate_vtable_metadata,
    write_metadata,
)

logger = logging.getLogger(__name__)

SORT_CHUNK_LINES = 1_000_000
SORT_MERGE_FAN_IN = 64
SORT_CHUNK_LINES_ENV = "GENOMATCH_SORT_CHUNK_LINES"
POSITIVE_INT_TOKEN_RE = re.compile(r"\+?0*[1-9]\d*\Z")
SIGNED_INT_TOKEN_RE = re.compile(SIGNED_INT_TOKEN_PATTERN + r"\Z")
VALID_NONMISSING_ALLELE_OPS = VALID_ALLELE_OPS - {"missing"}


@dataclass(frozen=True)
class ParsedSortRow:
    line: str
    rank: int
    pos_int: int
    identity: tuple[str, str, str, str]


@dataclass
class RunCursor:
    run_index: int
    handle: TextIO
    row_index: int = 0


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
        "--prefix",
        help="Optional scratch prefix for external-sort temporary files; not a variant-object output stem",
    )
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


def chunk_size_from_env() -> int:
    raw = os.environ.get(SORT_CHUNK_LINES_ENV)
    if raw is None:
        return SORT_CHUNK_LINES
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"{SORT_CHUNK_LINES_ENV} must be a positive integer") from exc
    if value <= 0:
        raise ValueError(f"{SORT_CHUNK_LINES_ENV} must be a positive integer")
    return value


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


def require_positive_pos(token: str, *, label: str) -> int:
    if not POSITIVE_INT_TOKEN_RE.fullmatch(token):
        raise ValueError(f"{label} row has invalid pos: {token!r}")
    value = int(token)
    if value <= 0:
        raise ValueError(f"{label} row has invalid pos: {token!r}")
    return value


def normalize_allele(value: str, *, label: str) -> str:
    allele = value.strip().upper()
    if not allele or re.fullmatch(r"[ACGT]+", allele) is None:
        raise ValueError(f"invalid allele code in {label}: {value!r}")
    return allele


def parse_source_index(token: str, *, label: str) -> int:
    stripped = token.strip()
    if not SIGNED_INT_TOKEN_RE.fullmatch(stripped):
        raise ValueError(f"invalid {label} row")
    return int(stripped)


def parse_variant_line(
    line: str,
    *,
    object_type: str,
    contig_naming: str,
    path: Path,
) -> ParsedSortRow:
    label = "vtable" if object_type == "variant_table" else "vmap"
    parts = line.rstrip("\n").split("\t")
    expected_columns = 5 if object_type == "variant_table" else 8
    if len(parts) != expected_columns:
        raise ValueError(f"invalid {label} row in {path}")

    chrom, pos, row_id = parts[0], parts[1], parts[2]
    if not chrom:
        raise ValueError(f"{label} row is missing chrom")
    if object_type == "variant_table" and not row_id:
        raise ValueError(f"{label} row is missing id")
    canonical = canonical_contig_from_label(chrom, contig_naming)
    if canonical is None or canonical == UNKNOWN_CONTIG:
        raise ValueError(
            f"variant object has contigs inconsistent with declared contig_naming={contig_naming!r}: {chrom!r}. "
            "Run normalize_contigs.py first"
        )
    pos_int = require_positive_pos(pos, label=label)
    a1 = normalize_allele(parts[3], label=label)
    a2 = normalize_allele(parts[4], label=label)
    normalized_parts = [chrom, pos, row_id, a1, a2]

    if object_type == "variant_map":
        source_shard = parts[5]
        source_index = parse_source_index(parts[6], label=label)
        allele_op = parts[7]
        if source_index == -1:
            if source_shard != MISSING_SOURCE_SHARD:
                raise ValueError("vmap row with source_index=-1 must have source_shard='.'")
            if allele_op != "missing":
                raise ValueError("vmap row with source_index=-1 must have allele_op=missing")
        elif source_index < 0:
            raise ValueError("vmap row source_index out of range")
        else:
            if not source_shard:
                raise ValueError("vmap row with source_index>=0 must define source_shard")
            if allele_op not in VALID_NONMISSING_ALLELE_OPS:
                raise ValueError(f"invalid allele_op in vmap row: {allele_op!r}")
        normalized_parts.extend([source_shard, str(source_index), allele_op])

    return ParsedSortRow(
        line="\t".join(normalized_parts) + "\n",
        rank=CANONICAL_CONTIG_RANK[canonical],
        pos_int=pos_int,
        identity=(chrom, pos, a1, a2),
    )


def iter_parsed_rows(path: Path, *, object_type: str, contig_naming: str) -> Iterator[ParsedSortRow]:
    # Assumes: metadata has already been validated and contig_naming is declared.
    # Performs: PN(allele uppercasing, source_index integer rendering) and PV/SV(row shape,
    # coordinate sortability, declared contig compatibility, vmap provenance fields).
    # Guarantees: parsed rows have stable line payloads and declared-coordinate sort keys.
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for line in handle:
            if not line.strip():
                continue
            yield parse_variant_line(line, object_type=object_type, contig_naming=contig_naming, path=path)


def scratch_directory(prefix: str | None) -> tempfile.TemporaryDirectory[str]:
    if prefix is None:
        return tempfile.TemporaryDirectory(prefix="genomatch-sort.")
    scratch_prefix = Path(prefix)
    scratch_dir = scratch_prefix.parent if str(scratch_prefix.parent) else Path(".")
    scratch_dir.mkdir(parents=True, exist_ok=True)
    return tempfile.TemporaryDirectory(prefix=scratch_prefix.name + ".", dir=str(scratch_dir))


def write_run(path: Path, rows: list[ParsedSortRow]) -> None:
    rows.sort(key=lambda row: (row.rank, row.pos_int))
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        # PERF: row loop retained for bounded-memory run emission; sorting requires per-row keys,
        # and output cost is I/O-bound over the configured chunk size.
        for row in rows:
            handle.write(row.line)


def create_sorted_runs(
    *,
    input_path: Path,
    object_type: str,
    contig_naming: str,
    chunk_lines: int,
    scratch_path: Path,
) -> tuple[list[Path], int]:
    runs: list[Path] = []
    chunk: list[ParsedSortRow] = []
    row_count = 0
    for row in iter_parsed_rows(input_path, object_type=object_type, contig_naming=contig_naming):
        chunk.append(row)
        row_count += 1
        if len(chunk) >= chunk_lines:
            run_path = scratch_path / f"run.{len(runs):06d}.tsv"
            write_run(run_path, chunk)
            runs.append(run_path)
            chunk = []
    if chunk:
        run_path = scratch_path / f"run.{len(runs):06d}.tsv"
        write_run(run_path, chunk)
        runs.append(run_path)
    return runs, row_count


def read_next_from_run(cursor: RunCursor, *, object_type: str, contig_naming: str, path: Path) -> tuple | None:
    while True:
        line = cursor.handle.readline()
        if not line:
            return None
        if not line.strip():
            continue
        row = parse_variant_line(line, object_type=object_type, contig_naming=contig_naming, path=path)
        heap_item = (row.rank, row.pos_int, cursor.run_index, cursor.row_index, row.line, row.identity, cursor)
        cursor.row_index += 1
        return heap_item


def merge_runs_once(
    *,
    run_paths: list[Path],
    output_path: Path,
    object_type: str,
    contig_naming: str,
    drop_duplicates: bool,
    enforce_vmap_uniqueness: bool,
) -> tuple[int, int]:
    cursors: list[RunCursor] = []
    heap: list[tuple] = []
    for run_index, run_path in enumerate(run_paths):
        handle = open(run_path, "r", encoding="utf-8", newline="")
        cursor = RunCursor(run_index=run_index, handle=handle)
        cursors.append(cursor)
        first = read_next_from_run(cursor, object_type=object_type, contig_naming=contig_naming, path=run_path)
        if first is not None:
            heapq.heappush(heap, first)

    emitted = 0
    dropped = 0
    current_coordinate: tuple[str, str] | None = None
    seen_identities_for_coordinate: set[tuple[str, str, str, str]] = set()
    try:
        with open(output_path, "w", encoding="utf-8", newline="\n") as output:
            while heap:
                _rank, _pos_int, _run_index, _row_index, line, identity, cursor = heapq.heappop(heap)
                coordinate = (identity[0], identity[1])
                if coordinate != current_coordinate:
                    current_coordinate = coordinate
                    seen_identities_for_coordinate = set()

                duplicate = identity in seen_identities_for_coordinate
                if duplicate and drop_duplicates:
                    dropped += 1
                elif duplicate and enforce_vmap_uniqueness:
                    raise ValueError("duplicate chrom:pos:a1:a2 in vmap target rows")
                else:
                    output.write(line)
                    emitted += 1
                    seen_identities_for_coordinate.add(identity)

                next_item = read_next_from_run(cursor, object_type=object_type, contig_naming=contig_naming, path=run_paths[cursor.run_index])
                if next_item is not None:
                    heapq.heappush(heap, next_item)
    finally:
        for cursor in cursors:
            cursor.handle.close()
    return emitted, dropped


def merge_runs_bounded(
    *,
    run_paths: list[Path],
    final_output_path: Path,
    object_type: str,
    contig_naming: str,
    scratch_path: Path,
    drop_duplicates: bool,
) -> tuple[int, int]:
    if not run_paths:
        with open(final_output_path, "w", encoding="utf-8", newline="\n"):
            pass
        return 0, 0

    current_runs = run_paths
    round_index = 0
    while len(current_runs) > SORT_MERGE_FAN_IN:
        next_runs: list[Path] = []
        for start in range(0, len(current_runs), SORT_MERGE_FAN_IN):
            group = current_runs[start : start + SORT_MERGE_FAN_IN]
            merged_path = scratch_path / f"merge.{round_index:03d}.{len(next_runs):06d}.tsv"
            merge_runs_once(
                run_paths=group,
                output_path=merged_path,
                object_type=object_type,
                contig_naming=contig_naming,
                drop_duplicates=False,
                enforce_vmap_uniqueness=False,
            )
            next_runs.append(merged_path)
        current_runs = next_runs
        round_index += 1

    return merge_runs_once(
        run_paths=current_runs,
        output_path=final_output_path,
        object_type=object_type,
        contig_naming=contig_naming,
        drop_duplicates=drop_duplicates,
        enforce_vmap_uniqueness=object_type == "variant_map" and not drop_duplicates,
    )


def write_sorted_output(
    *,
    input_path: Path,
    output_path: Path,
    metadata: dict,
    object_type: str,
    contig_naming: str,
    prefix: str | None,
    drop_duplicates: bool,
    chunk_lines: int,
) -> tuple[int, int]:
    ensure_parent_dir(output_path)
    temp_output = output_path.with_name(output_path.name + ".tmp")
    if temp_output.exists():
        temp_output.unlink()
    try:
        with scratch_directory(prefix) as scratch_dir_raw:
            scratch_path = Path(scratch_dir_raw)
            runs, _row_count = create_sorted_runs(
                input_path=input_path,
                object_type=object_type,
                contig_naming=contig_naming,
                chunk_lines=chunk_lines,
                scratch_path=scratch_path,
            )
            emitted, dropped = merge_runs_bounded(
                run_paths=runs,
                final_output_path=temp_output,
                object_type=object_type,
                contig_naming=contig_naming,
                scratch_path=scratch_path,
                drop_duplicates=drop_duplicates,
            )
        temp_output.replace(output_path)
    except Exception:
        if temp_output.exists():
            temp_output.unlink()
        raise
    write_metadata(output_path, metadata)
    qc_path = output_path.with_name(output_path.name + ".qc.tsv")
    if qc_path.exists():
        qc_path.unlink()
    return emitted, dropped


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    logger.info("sort_variants.py: sorting %s -> %s", input_path, output_path)
    if not input_path.exists():
        raise ValueError(f"input file not found: {input_path}")
    if object_type_for_path(input_path) != object_type_for_path(output_path):
        raise ValueError("sort_variants.py output suffix must match input suffix")

    metadata, contig_naming = load_sort_metadata(input_path)
    chunk_lines = chunk_size_from_env()
    emitted, dropped = write_sorted_output(
        input_path=input_path,
        output_path=output_path,
        metadata=dict(metadata),
        object_type=object_type_for_path(input_path),
        contig_naming=contig_naming,
        prefix=args.prefix,
        drop_duplicates=args.drop_duplicates,
        chunk_lines=chunk_lines,
    )
    if args.drop_duplicates:
        logger.info("sort_variants.py: dropped %s duplicate target rows.", dropped)
    logger.info("sort_variants.py: retained %s rows in declared coordinate order.", emitted)
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
