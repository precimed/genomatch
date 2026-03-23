#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set

from importer_utils import (
    filter_import_rows_by_chr2use,
    finalize_imported_vmap,
    ImportQcRow,
    ImportedVariantRow,
    is_canonical_import_allele,
    is_valid_import_position,
    reject_template_argument,
)
from sumstats_utils import (
    extract_variant_field,
    find_metadata_value,
    load_metadata,
    open_sumstats_data,
    resolve_column,
    resolve_variant_columns,
    split_line,
)
from vtable_utils import (
    load_metadata as load_variant_metadata,
    open_text,
    validate_vtable_metadata,
    VariantRow,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract a .vmap from summary statistics.")
    parser.add_argument("--input", required=True, help="Input summary statistics file")
    parser.add_argument("--output", required=True, help="Output .vmap file")
    parser.add_argument("--sumstats-metadata", required=True, help="Cleansumstats-style metadata YAML")
    parser.add_argument("--genome-build", default="unknown", help="Genome build for metadata")
    parser.add_argument("--id-vtable", help="Optional .vtable for ID-based coordinate enrichment")
    parser.add_argument("--chr2use", "--contigs", dest="chr2use", help="Comma-separated chromosome list or ranges")
    return parser.parse_args()


def load_id_lookup_vtable(path: Path) -> tuple[Dict[str, VariantRow], Set[str], Dict[str, object]]:
    if path.suffix != ".vtable":
        raise ValueError("--id-vtable must point to a .vtable")
    metadata = load_variant_metadata(path)
    validate_vtable_metadata(metadata)
    unique_matches: Dict[str, VariantRow] = {}
    ambiguous_ids: Set[str] = set()
    ignored_ids = 0
    with open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 5:
                raise ValueError(f"invalid vtable row in {path}: {line.strip()}")
            chrom, pos, row_id, a1, a2 = parts
            if not chrom or not is_valid_import_position(pos):
                raise ValueError(f"invalid vtable row in {path}: {line.strip()}")
            if not is_canonical_import_allele(a1) or not is_canonical_import_allele(a2):
                raise ValueError(f"invalid vtable row in {path}: {line.strip()}")
            lookup_id = row_id.strip()
            if not lookup_id or lookup_id == ".":
                ignored_ids += 1
                continue
            row = VariantRow(chrom, pos, row_id, a1, a2)
            if lookup_id in ambiguous_ids:
                continue
            if lookup_id in unique_matches:
                unique_matches.pop(lookup_id, None)
                ambiguous_ids.add(lookup_id)
                continue
            unique_matches[lookup_id] = row
    if ignored_ids:
        print(
            f"Warning: ignored {ignored_ids} --id-vtable rows whose id is missing, empty, or '.'",
            file=sys.stderr,
        )
    return unique_matches, ambiguous_ids, {
        "genome_build": metadata["genome_build"],
        "contig_naming": metadata.get("contig_naming"),
    }


def resolve_id_enrichment_columns(header: List[str], metadata: Dict[str, object]) -> tuple[int, int, int]:
    if find_metadata_value(metadata, "col_CHR") is not None or find_metadata_value(metadata, "col_POS") is not None:
        raise ValueError("--id-vtable requires metadata to omit col_CHR and col_POS")
    snp_idx = resolve_column(header, find_metadata_value(metadata, "col_SNP"), "col_SNP", required=True)
    effect_idx = resolve_column(
        header,
        find_metadata_value(metadata, "col_EffectAllele"),
        "col_EffectAllele",
        required=True,
    )
    other_idx = resolve_column(
        header,
        find_metadata_value(metadata, "col_OtherAllele"),
        "col_OtherAllele",
        required=True,
    )
    return snp_idx, effect_idx, other_idx


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    meta_path = Path(args.sumstats_metadata)
    output_path = Path(args.output)
    reject_template_argument(args.input, label="import_sumstats.py --input")
    if args.id_vtable:
        reject_template_argument(args.id_vtable, label="import_sumstats.py --id-vtable")
    if not input_path.exists():
        raise ValueError(f"sumstats file not found: {input_path}")
    if not meta_path.exists():
        raise ValueError(f"metadata file not found: {meta_path}")
    id_vtable_path: Optional[Path] = Path(args.id_vtable) if args.id_vtable else None
    if id_vtable_path is not None and not id_vtable_path.exists():
        raise ValueError(f"id-vtable not found: {id_vtable_path}")
    metadata: Dict[str, object] = load_metadata(meta_path)
    id_lookup_rows: Dict[str, VariantRow] = {}
    ambiguous_lookup_ids: Set[str] = set()
    inherited_target_meta: Optional[Dict[str, object]] = None
    if id_vtable_path is not None:
        id_lookup_rows, ambiguous_lookup_ids, inherited_target_meta = load_id_lookup_vtable(id_vtable_path)
    with open_sumstats_data(input_path) as (handle, _header_line, header, delimiter):
        if id_vtable_path is None:
            variant_columns = resolve_variant_columns(header, metadata, require_pos=True)
            max_idx = max(
                variant_columns.chr,
                variant_columns.pos or 0,
                variant_columns.effect_allele,
                variant_columns.other_allele,
                variant_columns.snp or 0,
            )
        else:
            snp_idx, effect_idx, other_idx = resolve_id_enrichment_columns(header, metadata)
            max_idx = max(snp_idx, effect_idx, other_idx)
        rows: List[ImportedVariantRow] = []
        qc_rows: List[ImportQcRow] = []
        source_index = 0
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = split_line(line, delimiter)
            if len(cols) <= max_idx:
                qc_rows.append(ImportQcRow(".", source_index, "malformed_row"))
                source_index += 1
                continue
            try:
                if id_vtable_path is None:
                    chrom = extract_variant_field(cols[variant_columns.chr], "CHR")
                    pos = extract_variant_field(cols[variant_columns.pos], "POS")
                    row_id = cols[variant_columns.snp] if variant_columns.snp is not None and cols[variant_columns.snp] else "."
                    a1 = extract_variant_field(cols[variant_columns.effect_allele], "EffectAllele")
                    a2 = extract_variant_field(cols[variant_columns.other_allele], "OtherAllele")
                else:
                    raw_id = cols[snp_idx].strip()
                    a1 = extract_variant_field(cols[effect_idx], "EffectAllele")
                    a2 = extract_variant_field(cols[other_idx], "OtherAllele")
            except Exception:
                qc_rows.append(ImportQcRow(".", source_index, "malformed_row"))
                source_index += 1
                continue
            if not is_canonical_import_allele(a1) or not is_canonical_import_allele(a2):
                qc_rows.append(ImportQcRow(".", source_index, "non_actg_allele"))
                source_index += 1
                continue
            if id_vtable_path is None:
                if not chrom or not is_valid_import_position(pos):
                    qc_rows.append(ImportQcRow(".", source_index, "malformed_row"))
                    source_index += 1
                    continue
            else:
                if not raw_id or raw_id == ".":
                    qc_rows.append(ImportQcRow(".", source_index, "invalid_id"))
                    source_index += 1
                    continue
                if raw_id in ambiguous_lookup_ids:
                    qc_rows.append(ImportQcRow(".", source_index, "ambiguous_id_match"))
                    source_index += 1
                    continue
                matched_row = id_lookup_rows.get(raw_id)
                if matched_row is None:
                    qc_rows.append(ImportQcRow(".", source_index, "id_not_found"))
                    source_index += 1
                    continue
                chrom = matched_row.chrom
                pos = matched_row.pos
                row_id = raw_id
            rows.append(
                ImportedVariantRow(
                    VariantRow(
                        chrom,
                        pos,
                        row_id,
                        a1,
                        a2,
                    ),
                    ".",
                    source_index,
                )
            )
            source_index += 1
    rows, chr_qc_rows = filter_import_rows_by_chr2use(rows, args.chr2use)
    qc_rows.extend(chr_qc_rows)
    finalize_imported_vmap(
        output_path=output_path,
        rows=rows,
        genome_build=args.genome_build if inherited_target_meta is None else str(inherited_target_meta["genome_build"]),
        target_contig_naming=None if inherited_target_meta is None else inherited_target_meta.get("contig_naming"),
        infer_target_contig_naming=inherited_target_meta is None,
        created_by="import_sumstats.py",
        derived_from=input_path,
        qc_rows=qc_rows,
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
