from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

from .contig_utils import (
    canonical_contig_from_any_supported_label,
    canonical_contig_from_label,
    SUPPORTED_CONTIG_NAMINGS,
    supported_exact_contig_tokens,
)
from .vtable_utils import (
    COMPLEMENT,
    ensure_parent_dir,
    MISSING_SOURCE_SHARD,
    VMapRow,
    VariantRow,
    infer_contig_naming,
    make_vmap_metadata,
    parse_chr2use,
    write_metadata,
    write_vmap,
)


@dataclass(frozen=True)
class DiscoveredInputShard:
    path: Path
    source_shard: str


@dataclass(frozen=True)
class ImportedVariantRow:
    row: VariantRow
    source_shard: str
    source_index: int


@dataclass(frozen=True)
class ImportQcRow:
    source_shard: str
    source_index: int
    reason: str


def finalize_imported_vmap(
    *,
    output_path: Path,
    rows: Iterable[ImportedVariantRow],
    genome_build: str,
    target_contig_naming: str | None = None,
    infer_target_contig_naming: bool = True,
    created_by: str,
    derived_from: Path,
    qc_rows: Sequence[ImportQcRow],
) -> None:
    rows_list = list(rows)
    chroms = [item.row.chrom for item in rows_list]
    if infer_target_contig_naming:
        contig_naming = infer_contig_naming(chroms)
        if contig_naming is None and importer_should_warn_for_contigs(chroms):
            print(
                "Warning: retained rows include mixed or invalid contig labels; preserving contigs as-is, "
                "omitting contig_naming from metadata, and requiring normalize_contigs.py before downstream use.",
                file=sys.stderr,
            )
    else:
        contig_naming = target_contig_naming
    write_vmap(
        output_path,
        [
            VMapRow(
                item.row.chrom,
                item.row.pos,
                item.row.id,
                item.row.a1,
                item.row.a2,
                item.source_shard,
                item.source_index,
                "identity",
            )
            for item in rows_list
        ],
    )
    metadata = make_vmap_metadata(
        {"genome_build": genome_build, "contig_naming": contig_naming},
        provenance={"created_by": created_by, "derived_from": str(derived_from)},
    )
    write_metadata(output_path, metadata)
    if qc_rows:
        write_import_qc(output_path.with_name(output_path.name + ".qc.tsv"), qc_rows)


def write_import_qc(path: Path, rows: Sequence[ImportQcRow]) -> None:
    ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("source_shard\tsource_index\treason\n")
        for row in rows:
            handle.write(f"{row.source_shard}\t{row.source_index}\t{row.reason}\n")


def resolve_import_input_paths(input_arg: str, *, kind_label: str) -> List[DiscoveredInputShard]:
    input_path = Path(input_arg)
    if "@" not in input_arg:
        if not input_path.exists():
            raise ValueError(f"input {kind_label} not found: {input_path}")
        return [DiscoveredInputShard(path=input_path, source_shard=MISSING_SOURCE_SHARD)]
    name_prefix, name_suffix = input_path.name.split("@", 1)
    discovered: List[DiscoveredInputShard] = []
    for token in supported_exact_contig_tokens():
        candidate = input_path.parent / f"{name_prefix}{token}{name_suffix}"
        if candidate.is_file():
            discovered.append(DiscoveredInputShard(path=candidate, source_shard=token))
    if not discovered:
        raise ValueError(f"no input shards found for template: {input_path}")
    return sorted(discovered, key=lambda shard: str(shard.path))


def importer_should_warn_for_contigs(labels: Sequence[str]) -> bool:
    seen = [label.strip() for label in labels if label and label.strip()]
    if not seen:
        return False
    if infer_contig_naming(seen) is not None:
        return False
    if any(
        not any(canonical_contig_from_label(label, naming) is not None for naming in SUPPORTED_CONTIG_NAMINGS)
        for label in seen
    ):
        return True
    return not any(
        all(canonical_contig_from_label(label, naming) is not None for label in seen)
        for naming in SUPPORTED_CONTIG_NAMINGS
    )


def is_canonical_import_allele(value: str) -> bool:
    return bool(value) and all(base in COMPLEMENT for base in value)


def is_valid_import_position(value: str) -> bool:
    try:
        return int(value) > 0
    except Exception:
        return False


def filter_import_rows_by_chr2use(
    rows: Sequence[ImportedVariantRow],
    chr2use_raw: str | None,
) -> tuple[List[ImportedVariantRow], List[ImportQcRow]]:
    allowed, explicit = parse_chr2use(chr2use_raw)
    if not explicit:
        return list(rows), []
    allowed_set = set(allowed)
    filtered_rows: List[ImportedVariantRow] = []
    qc_rows: List[ImportQcRow] = []
    for item in rows:
        canonical = canonical_contig_from_any_supported_label(item.row.chrom)
        if canonical in allowed_set:
            filtered_rows.append(item)
        else:
            qc_rows.append(ImportQcRow(item.source_shard, item.source_index, "filtered_by_chr2use"))
    return filtered_rows, qc_rows


def reject_template_argument(raw: str, *, label: str) -> None:
    if "@" in raw:
        raise ValueError(f"{label} does not accept '@' paths")
