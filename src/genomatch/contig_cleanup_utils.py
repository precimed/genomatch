from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Union

from .vtable_utils import (
    LoadedVariantObject,
    UNKNOWN_CONTIG,
    VMapRow,
    VariantRow,
    canonical_contig_from_label,
    load_variant_object,
    repair_contig_label,
    require_contig_naming,
    require_rows_match_contig_naming,
    write_metadata,
    write_vmap,
    write_vtable,
)

TargetRow = Union[VariantRow, VMapRow]

@dataclass(frozen=True)
class NormalizedTargetRows:
    rows: List[TargetRow]
    normalized_count: int
    unknown_count: int


def load_target_variant_object(path: Path) -> LoadedVariantObject:
    return load_variant_object(path)


def _replace_row_chrom(row: TargetRow, chrom: str) -> TargetRow:
    if isinstance(row, VariantRow):
        return VariantRow(chrom, row.pos, row.id, row.a1, row.a2)
    return type(row)(chrom, row.pos, row.id, row.a1, row.a2, row.source_shard, row.source_index, row.allele_op)


def normalize_target_rows(rows: Sequence[TargetRow], to_naming: str, *, genome_build: str | None = None) -> NormalizedTargetRows:
    normalized_rows: List[TargetRow] = []
    normalized_count = 0
    unknown_count = 0
    for row in rows:
        chrom = repair_contig_label(row.chrom, to_naming, pos=row.pos, genome_build=genome_build)
        if chrom == UNKNOWN_CONTIG:
            unknown_count += 1
            continue
        if chrom != row.chrom:
            normalized_count += 1
        normalized_rows.append(_replace_row_chrom(row, chrom))
    return NormalizedTargetRows(normalized_rows, normalized_count, unknown_count)


def require_target_contig_cleanup_contract(loaded: LoadedVariantObject) -> Dict[str, object]:
    target_meta = dict(loaded.target_metadata)
    contig_naming = require_contig_naming(target_meta, label="variant object")
    rows = loaded.base_vmap_rows if loaded.base_vmap_rows is not None else loaded.target_rows
    require_rows_match_contig_naming(rows, contig_naming, label="variant object")
    return {
        "genome_build": target_meta["genome_build"],
        "contig_naming": contig_naming,
    }


def build_target_restriction_selection_for_contigs(
    rows: Sequence[TargetRow],
    *,
    allowed_canonical: set[str] | None,
    contig_naming: str | None,
) -> List[TargetRow]:
    selected: List[TargetRow] = []
    for row in rows:
        if allowed_canonical is not None:
            if contig_naming is None:
                raise ValueError("contig_naming is required when restricting target rows by contig")
            canonical = canonical_contig_from_label(row.chrom, contig_naming)
            if canonical not in allowed_canonical:
                continue
        selected.append(row)
    return selected


def normalized_output_metadata(loaded: LoadedVariantObject, to_naming: str) -> Dict[str, object]:
    if loaded.object_type == "variant_table":
        out_meta = dict(loaded.raw_metadata)
        out_meta["contig_naming"] = to_naming
        return out_meta
    out_meta = dict(loaded.raw_metadata)
    out_meta["target"] = dict(loaded.target_metadata)
    out_meta["target"]["contig_naming"] = to_naming
    return out_meta


def write_variant_object_like_input(
    loaded: LoadedVariantObject,
    output_path: Path,
    rows: Sequence[TargetRow],
    metadata: Dict[str, object],
    *,
    preserve_qc: bool = False,
) -> None:
    qc_path = output_path.with_name(output_path.name + ".qc.tsv")
    if loaded.object_type == "variant_table":
        write_vtable(output_path, rows)  # type: ignore[arg-type]
    else:
        write_vmap(output_path, rows)  # type: ignore[arg-type]
    write_metadata(output_path, metadata)
    if not preserve_qc and qc_path.exists():
        qc_path.unlink()


def restriction_output_metadata(loaded: LoadedVariantObject) -> Dict[str, object]:
    require_target_contig_cleanup_contract(loaded)
    return dict(loaded.raw_metadata)
