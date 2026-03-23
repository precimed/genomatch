#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from contig_utils import (
    CANONICAL_CONTIG_ORDER,
    NUMERIC_CONTIG_NAMINGS,
    SUPPORTED_CONTIG_NAMINGS,
    UNKNOWN_CONTIG,
    canonical_contig_from_label,
    contig_label_for_naming,
    convert_contig_label,
    normalize_chrom_label,
    normalize_contig_for_reference,
    repair_contig_label,
)

SUPPORTED_GENOME_BUILDS = {"GRCh37", "GRCh38", "unknown"}
VALID_ALLELE_OPS = {"identity", "swap", "flip", "flip_swap", "missing"}
ALLELE_OP_COMPOSITION = {
    ("identity", "identity"): "identity",
    ("identity", "swap"): "swap",
    ("identity", "flip"): "flip",
    ("identity", "flip_swap"): "flip_swap",
    ("swap", "identity"): "swap",
    ("swap", "swap"): "identity",
    ("swap", "flip"): "flip_swap",
    ("swap", "flip_swap"): "flip",
    ("flip", "identity"): "flip",
    ("flip", "swap"): "flip_swap",
    ("flip", "flip"): "identity",
    ("flip", "flip_swap"): "swap",
    ("flip_swap", "identity"): "flip_swap",
    ("flip_swap", "swap"): "flip",
    ("flip_swap", "flip"): "swap",
    ("flip_swap", "flip_swap"): "identity",
}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
MISSING_SOURCE_SHARD = "."
CANONICAL_CONTIG_RANK = {label: idx for idx, label in enumerate(CANONICAL_CONTIG_ORDER, start=1)}


@dataclass(frozen=True)
class VariantRow:
    chrom: str
    pos: str
    id: str
    a1: str
    a2: str


@dataclass(frozen=True)
class VMapRow:
    chrom: str
    pos: str
    id: str
    a1: str
    a2: str
    source_shard: str
    source_index: int
    allele_op: str


@dataclass(frozen=True)
class LoadedVariantObject:
    path: Path
    object_type: str
    target_rows: List[VariantRow]
    target_metadata: Dict[str, object]
    raw_metadata: Dict[str, object]
    base_vmap_rows: Optional[List[VMapRow]]


def open_text(path: Path, mode: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def metadata_path_for(path: Path) -> Path:
    suffix = "".join(path.suffixes)
    if suffix.endswith(".vtable"):
        return path.with_name(path.name + ".meta.json")
    if suffix.endswith(".vmap"):
        return path.with_name(path.name + ".meta.json")
    raise ValueError(f"unsupported variant object path: {path}")


def ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _load_json(path: Path) -> Dict[str, object]:
    with open(path, "r", encoding="utf-8", newline="\n") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"metadata must be a JSON object: {path}")
    return data


def load_metadata(path: Path) -> Dict[str, object]:
    meta_path = metadata_path_for(path)
    if not meta_path.exists():
        raise ValueError(f"metadata sidecar not found: {meta_path}")
    return _load_json(meta_path)


def write_metadata(path: Path, metadata: Dict[str, object]) -> None:
    validate_metadata_for_path(path, metadata)
    meta_path = metadata_path_for(path)
    ensure_parent_dir(meta_path)
    with open(meta_path, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def validate_metadata_for_path(path: Path, metadata: Dict[str, object]) -> None:
    if path.name.endswith(".vtable"):
        validate_vtable_metadata(metadata)
        return
    if path.name.endswith(".vmap"):
        validate_vmap_metadata(metadata)
        return
    raise ValueError(f"unsupported variant object path: {path}")


def validate_vtable_metadata(metadata: Dict[str, object]) -> None:
    if metadata.get("object_type") != "variant_table":
        raise ValueError("vtable metadata must have object_type=variant_table")
    genome_build = metadata.get("genome_build")
    contig_naming = metadata.get("contig_naming")
    if genome_build not in SUPPORTED_GENOME_BUILDS:
        raise ValueError(f"unsupported genome_build in metadata: {genome_build!r}")
    if contig_naming is not None and contig_naming not in SUPPORTED_CONTIG_NAMINGS:
        raise ValueError(f"unsupported contig_naming in metadata: {contig_naming!r}")


def validate_vmap_metadata(metadata: Dict[str, object]) -> None:
    if metadata.get("object_type") != "variant_map":
        raise ValueError("vmap metadata must have object_type=variant_map")
    target = metadata.get("target")
    if not isinstance(target, dict):
        raise ValueError("vmap metadata must define a target object")
    genome_build = target.get("genome_build")
    contig_naming = target.get("contig_naming")
    if genome_build not in SUPPORTED_GENOME_BUILDS:
        raise ValueError(f"unsupported target.genome_build in metadata: {genome_build!r}")
    if contig_naming is not None and contig_naming not in SUPPORTED_CONTIG_NAMINGS:
        raise ValueError(f"unsupported target.contig_naming in metadata: {contig_naming!r}")


def infer_contig_naming(labels: Sequence[str]) -> Optional[str]:
    seen = [label.strip() for label in labels if label and label.strip()]
    if not seen:
        return None
    for naming in ("ncbi", "ucsc", "plink", "plink_splitx"):
        if all(canonical_contig_from_label(label, naming) is not None for label in seen):
            return naming
    return None


def make_vtable_metadata(
    genome_build: str,
    contig_naming: Optional[str],
    provenance: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    metadata: Dict[str, object] = {
        "object_type": "variant_table",
        "genome_build": genome_build,
    }
    if contig_naming is not None:
        metadata["contig_naming"] = contig_naming
    if provenance:
        metadata.update(provenance)
    validate_vtable_metadata(metadata)
    return metadata


def make_vmap_metadata(
    target_meta: Dict[str, object],
    provenance: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    metadata: Dict[str, object] = {
        "object_type": "variant_map",
        "target": {
            "genome_build": target_meta["genome_build"],
        },
    }
    if target_meta.get("contig_naming") is not None:
        metadata["target"]["contig_naming"] = target_meta["contig_naming"]
    if provenance:
        metadata.update(provenance)
    validate_vmap_metadata(metadata)
    return metadata


def read_vtable(path: Path) -> List[VariantRow]:
    rows: List[VariantRow] = []
    with open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 5:
                raise ValueError(f"invalid vtable row in {path}: {line.strip()}")
            rows.append(VariantRow(*parts))
    validate_vtable_rows(rows, label=str(path))
    return rows


def write_vtable(path: Path, rows: Iterable[VariantRow]) -> None:
    rows_list = list(rows)
    validate_vtable_rows(rows_list, label=str(path))
    ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for row in rows_list:
            handle.write("\t".join([row.chrom, row.pos, row.id, row.a1, row.a2]) + "\n")


def read_vmap(path: Path) -> List[VMapRow]:
    rows: List[VMapRow] = []
    with open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 8:
                raise ValueError(f"invalid vmap row in {path}: {line.strip()}")
            try:
                source_index = int(parts[6])
            except ValueError as exc:
                raise ValueError(f"invalid vmap row in {path}: {line.strip()}") from exc
            rows.append(VMapRow(parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], source_index, parts[7]))
    validate_vmap_rows(rows)
    return rows


def write_vmap(path: Path, rows: Iterable[VMapRow]) -> None:
    rows_list = list(rows)
    validate_vmap_rows(rows_list)
    ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for row in rows_list:
            handle.write(
                "\t".join(
                    [
                        row.chrom,
                        row.pos,
                        row.id,
                        row.a1,
                        row.a2,
                        row.source_shard,
                        str(row.source_index),
                        row.allele_op,
                    ]
                )
                + "\n"
            )


def validate_vtable_rows(rows: Sequence[VariantRow], label: str = "vtable") -> None:
    for row in rows:
        if not row.chrom:
            raise ValueError(f"{label} row is missing chrom")
        try:
            pos = int(row.pos)
        except ValueError as exc:
            raise ValueError(f"{label} row has invalid pos: {row.pos!r}") from exc
        if pos <= 0:
            raise ValueError(f"{label} row has invalid pos: {row.pos!r}")
        if not row.id:
            raise ValueError(f"{label} row is missing id")
        validate_allele_value(row.a1, label=label)
        validate_allele_value(row.a2, label=label)


def validate_vmap_rows(rows: Sequence[VMapRow]) -> None:
    seen = set()
    for row in rows:
        validate_allele_value(row.a1, label="vmap")
        validate_allele_value(row.a2, label="vmap")
        try:
            pos = int(row.pos)
        except ValueError as exc:
            raise ValueError(f"vmap row has invalid pos: {row.pos!r}") from exc
        if pos <= 0:
            raise ValueError(f"vmap row has invalid pos: {row.pos!r}")
        key = target_row_key(row)
        if key in seen:
            raise ValueError("duplicate chrom:pos:a1:a2 in vmap target rows")
        seen.add(key)
        if row.source_index == -1:
            if row.source_shard != MISSING_SOURCE_SHARD:
                raise ValueError("vmap row with source_index=-1 must have source_shard='.'")
            if row.allele_op != "missing":
                raise ValueError("vmap row with source_index=-1 must have allele_op=missing")
            continue
        if row.source_index < 0:
            raise ValueError("vmap row source_index out of range")
        if not row.source_shard:
            raise ValueError("vmap row with source_index>=0 must define source_shard")
        if row.allele_op not in VALID_ALLELE_OPS - {"missing"}:
            raise ValueError(f"invalid allele_op in vmap row: {row.allele_op!r}")


def variant_rows_from_vmap_rows(rows: Sequence[VMapRow]) -> List[VariantRow]:
    return [VariantRow(row.chrom, row.pos, row.id, row.a1, row.a2) for row in rows]


def target_row_key(row: VariantRow | VMapRow) -> tuple[str, str, str, str]:
    return (row.chrom, row.pos, row.a1.upper(), row.a2.upper())


def duplicate_target_row_keys(rows: Sequence[VariantRow | VMapRow]) -> set[tuple[str, str, str, str]]:
    seen = set()
    duplicate_keys = set()
    for row in rows:
        key = target_row_key(row)
        if key in seen:
            duplicate_keys.add(key)
        else:
            seen.add(key)
    return duplicate_keys


def write_vmap_status_qc(path: Path, rows: Iterable[tuple[str, int, str, str]]) -> None:
    ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("source_shard\tsource_index\tsource_id\tstatus\n")
        for source_shard, source_index, source_id, status in rows:
            handle.write(f"{source_shard}\t{source_index}\t{source_id}\t{status}\n")


def load_variant_object(path: Path) -> LoadedVariantObject:
    metadata = load_metadata(path)
    if path.name.endswith(".vtable"):
        validate_vtable_metadata(metadata)
        return LoadedVariantObject(
            path=path,
            object_type="variant_table",
            target_rows=read_vtable(path),
            target_metadata={
                "genome_build": metadata["genome_build"],
                **({"contig_naming": metadata["contig_naming"]} if "contig_naming" in metadata else {}),
            },
            raw_metadata=metadata,
            base_vmap_rows=None,
        )
    if path.name.endswith(".vmap"):
        validate_vmap_metadata(metadata)
        vmap_rows = read_vmap(path)
        return LoadedVariantObject(
            path=path,
            object_type="variant_map",
            target_rows=variant_rows_from_vmap_rows(vmap_rows),
            target_metadata=dict(metadata["target"]),
            raw_metadata=metadata,
            base_vmap_rows=vmap_rows,
        )
    raise ValueError(f"unsupported variant object path: {path}")


def compose_allele_ops(first: str, second: str) -> str:
    if first == "missing" or second == "missing":
        return "missing"
    composed = ALLELE_OP_COMPOSITION.get((first, second))
    if composed is None:
        raise ValueError(f"cannot compose allele operations: {first!r}, {second!r}")
    return composed


def shard_local_provenance(shards: Sequence[str]) -> List[Tuple[str, int]]:
    counts: Dict[str, int] = {}
    out: List[Tuple[str, int]] = []
    for shard in shards:
        if not shard:
            raise ValueError("source shard is missing")
        index = counts.get(shard, 0)
        out.append((shard, index))
        counts[shard] = index + 1
    return out


def possible_fasta_contigs(label: str, contig_naming: str) -> List[str]:
    canonical = canonical_contig_from_label(label, contig_naming)
    if not canonical:
        return []
    names = [contig_label_for_naming(canonical, contig_naming)]
    ucsc_name = contig_label_for_naming(canonical, "ucsc")
    ncbi_name = contig_label_for_naming(canonical, "ncbi")
    if ucsc_name not in names:
        names.append(ucsc_name)
    if ncbi_name not in names:
        names.append(ncbi_name)
    return names


def validate_allele_value(value: str, *, label: str) -> None:
    token = value.strip().upper()
    if not token:
        raise ValueError(f"invalid allele code in {label}: {value!r}")
    if any(base not in COMPLEMENT for base in token):
        raise ValueError(f"invalid allele code in {label}: {value!r}")


def validate_snv_alleles(row: VariantRow, *, label: str) -> None:
    for allele_name, allele in (("a1", row.a1), ("a2", row.a2)):
        token = allele.strip().upper()
        if len(token) != 1 or token not in COMPLEMENT:
            raise ValueError(
                f"{label} currently supports only single-base A/C/G/T alleles; "
                f"found {allele_name}={allele!r} at {row.chrom}:{row.pos}"
            )


def require_contig_naming(metadata: Dict[str, object], *, label: str) -> str:
    contig_naming = metadata.get("contig_naming")
    if contig_naming not in SUPPORTED_CONTIG_NAMINGS:
        raise ValueError(
            f"{label} metadata must define contig_naming before this operation can run; "
            "run normalize_contigs.py first"
        )
    return str(contig_naming)


def canonical_contig_from_row(row: object, contig_naming: str, *, label: str) -> str:
    chrom = getattr(row, "chrom", None)
    if not isinstance(chrom, str):
        raise ValueError(f"{label} row is missing chrom")
    canonical = canonical_contig_from_label(chrom, contig_naming)
    if canonical is None:
        raise ValueError(
            f"{label} has contigs inconsistent with declared contig_naming={contig_naming!r}; "
            f"first invalid label: {chrom!r}. Run normalize_contigs.py first"
        )
    if canonical == UNKNOWN_CONTIG:
        raise ValueError(f"{label} contains chrom=unknown; run normalize_contigs.py first")
    return canonical


def declared_coordinate_sort_key(row: object, contig_naming: str, *, label: str) -> Tuple[int, int]:
    canonical = canonical_contig_from_row(row, contig_naming, label=label)
    try:
        pos = int(getattr(row, "pos"))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} row has invalid pos: {getattr(row, 'pos', None)!r}") from exc
    if pos <= 0:
        raise ValueError(f"{label} row has invalid pos: {getattr(row, 'pos', None)!r}")
    return (CANONICAL_CONTIG_RANK[canonical], pos)


def sort_target_rows_by_declared_coordinate(rows: Sequence[object], contig_naming: str, *, label: str) -> List[object]:
    require_rows_match_contig_naming(rows, contig_naming, label=label)
    return sorted(rows, key=lambda row: declared_coordinate_sort_key(row, contig_naming, label=label))


def require_rows_match_contig_naming(rows: Sequence[object], contig_naming: str, *, label: str) -> None:
    invalid: List[str] = []
    for row in rows:
        chrom = getattr(row, "chrom", None)
        if not isinstance(chrom, str):
            invalid.append("<missing>")
            continue
        canonical = canonical_contig_from_label(chrom, contig_naming)
        if canonical is None or canonical == UNKNOWN_CONTIG:
            invalid.append(chrom)
            if len(invalid) == 3:
                break
    if invalid:
        examples = ", ".join(repr(item) for item in invalid)
        raise ValueError(
            f"{label} has contigs inconsistent with declared contig_naming={contig_naming!r}: {examples}. "
            "Run normalize_contigs.py first"
        )


def parse_chr2use(raw: Optional[str]) -> Tuple[List[str], bool]:
    explicit = raw is not None and raw.strip() != ""
    if not explicit:
        return [], False
    spec = raw.strip()
    tokens = [token.strip() for token in spec.split(",") if token.strip()]
    if not tokens:
        raise ValueError("--chr2use is empty")
    labels: List[str] = []
    seen = set()
    for token in tokens:
        token_lower = token.lower()
        if token_lower.startswith("chr"):
            token_lower = token_lower[3:]
        if "-" in token_lower:
            start_raw, end_raw = token_lower.split("-", 1)
            if not start_raw.isdigit() or not end_raw.isdigit():
                raise ValueError(f"invalid --chr2use token: {token!r}")
            start = int(start_raw)
            end = int(end_raw)
            if start <= 0 or end <= 0 or start > end:
                raise ValueError(f"invalid --chr2use range: {token!r}")
            for value in range(start, end + 1):
                canonical = normalize_chrom_label(str(value))
                if not canonical:
                    raise ValueError(f"invalid --chr2use value: {value}")
                if canonical not in seen:
                    seen.add(canonical)
                    labels.append(canonical)
            continue
        canonical = normalize_chrom_label(token)
        if not canonical:
            raise ValueError(f"invalid --chr2use token: {token!r}")
        if canonical not in seen:
            seen.add(canonical)
            labels.append(canonical)
    return labels, explicit


def filter_variant_rows_by_chr(
    rows: Sequence[VariantRow],
    labels: Sequence[str],
    contig_naming: str,
) -> List[VariantRow]:
    allowed = set(labels)
    return [row for row in rows if canonical_contig_from_row(row, contig_naming, label="variant table") in allowed]


def complement_allele(value: str) -> str:
    token = value.upper()
    if any(base not in COMPLEMENT for base in token):
        raise ValueError(f"cannot complement non-ACGT allele: {value!r}")
    return "".join(COMPLEMENT[base] for base in reversed(token))


def classify_allele_operation(
    source_a1: str,
    source_a2: str,
    target_a1: str,
    target_a2: str,
    allow_strand_flips: bool,
) -> str:
    source = (source_a1.upper(), source_a2.upper())
    target = (target_a1.upper(), target_a2.upper())
    if source == target:
        return "identity"
    if (source[1], source[0]) == target:
        return "swap"
    if not allow_strand_flips:
        return "missing"
    try:
        flipped = (complement_allele(source[0]), complement_allele(source[1]))
    except ValueError:
        return "missing"
    if flipped == target:
        return "flip"
    if (flipped[1], flipped[0]) == target:
        return "flip_swap"
    return "missing"
