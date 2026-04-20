#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from .tabular_rows import VMapRow, VMapRowsTable, VariantRow, VariantRowsTable
from .contig_utils import (
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
POSITIVE_INT_TOKEN_PATTERN = r"\+?0*[1-9]\d*"
SIGNED_INT_TOKEN_PATTERN = r"[+-]?\d+"


@dataclass(frozen=True)
class LoadedVariantObject:
    path: Path
    object_type: str
    target_rows: List[VariantRow]
    target_metadata: Dict[str, object]
    raw_metadata: Dict[str, object]
    base_vmap_rows: Optional[List[VMapRow]]


@dataclass(frozen=True)
class LoadedVariantObjectTables:
    path: Path
    object_type: str
    target_rows_table: VariantRowsTable
    target_metadata: Dict[str, object]
    raw_metadata: Dict[str, object]
    base_vmap_table: Optional[VMapRowsTable]


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


WRITE_CHUNK_ROWS = 100_000


def normalize_allele_token(value: str) -> str:
    return value.strip().upper()


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
    return read_vtable_table(path).to_rows()


def read_vtable_table(path: Path) -> VariantRowsTable:
    # Assumes: path points to a tab-delimited vtable payload.
    # Performs: PN(allele canonicalization), PV(row-shape/type parseability), SV/CV(frame validation).
    # Guarantees: VariantRowsTable with canonical uppercase a1/a2 and validated vtable invariants.
    try:
        frame = pd.read_table(
            path,
            sep="\t",
            header=None,
            names=["chrom", "pos", "id", "a1", "a2"],
            dtype="object",
            keep_default_na=False,
            na_values=[],
            skip_blank_lines=True,
            compression="infer",
        )
    except Exception as exc:
        raise ValueError(f"invalid vtable row in {path}") from exc
    if frame.isna().any(axis=None):
        raise ValueError(f"invalid vtable row in {path}")
    table = VariantRowsTable.from_frame(frame, copy=False)
    frame = table.to_frame(copy=False)
    frame["a1"] = frame["a1"].str.strip().str.upper()
    frame["a2"] = frame["a2"].str.strip().str.upper()
    _validate_vtable_frame(frame, label=str(path), assume_normalized_alleles=True)
    return table


def write_vtable(path: Path, rows: Iterable[VariantRow]) -> None:
    rows_list = list(rows)
    validate_vtable_rows(rows_list, label=str(path))
    ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for row in rows_list:
            handle.write("\t".join([row.chrom, row.pos, row.id, row.a1, row.a2]) + "\n")


def write_vtable_table(path: Path, table: VariantRowsTable, *, assume_validated: bool = False) -> None:
    # Assumes: table columns conform to VariantRowsTable schema.
    # Performs: CV(vtable invariants before emission) unless assume_validated, file emission.
    # Guarantees: on-disk .vtable rows satisfy contract-level schema/invariants.
    frame = table.to_frame(copy=False)
    if not assume_validated:
        _validate_vtable_frame(frame, label=str(path))
    ensure_parent_dir(path)
    lines = (
        frame["chrom"].astype(str) + "\t"
        + frame["pos"].astype(str) + "\t"
        + frame["id"].astype(str) + "\t"
        + frame["a1"].astype(str) + "\t"
        + frame["a2"].astype(str) + "\n"
    )
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        _write_series_in_chunks(handle, lines)


def read_vmap(path: Path) -> List[VMapRow]:
    return read_vmap_table(path).to_rows()


def read_vmap_table(path: Path) -> VMapRowsTable:
    # Assumes: path points to a tab-delimited vmap payload.
    # Performs: PN(allele canonicalization), PV(row-shape/type parseability), SV/CV(frame validation).
    # Guarantees: VMapRowsTable with canonical uppercase a1/a2, int64 source_index, and validated vmap invariants.
    try:
        frame = pd.read_table(
            path,
            sep="\t",
            header=None,
            names=["chrom", "pos", "id", "a1", "a2", "source_shard", "source_index", "allele_op"],
            dtype="object",
            keep_default_na=False,
            na_values=[],
            skip_blank_lines=True,
            compression="infer",
        )
    except Exception as exc:
        raise ValueError(f"invalid vmap row in {path}") from exc
    if frame.isna().any(axis=None):
        raise ValueError(f"invalid vmap row in {path}")
    source_index = frame["source_index"].str.strip()
    source_index_numeric = pd.to_numeric(source_index, errors="coerce")
    valid_source_index_tokens = source_index.str.fullmatch(SIGNED_INT_TOKEN_PATTERN).fillna(False)
    invalid_source_index = ~valid_source_index_tokens | source_index_numeric.isna()
    if bool(invalid_source_index.any()):
        raise ValueError(f"invalid vmap row in {path}")
    frame["source_index"] = source_index_numeric.astype("int64")
    table = VMapRowsTable.from_frame(frame, copy=False)
    frame = table.to_frame(copy=False)
    frame["a1"] = frame["a1"].str.strip().str.upper()
    frame["a2"] = frame["a2"].str.strip().str.upper()
    _validate_vmap_frame(frame, assume_normalized_alleles=True)
    return table


def write_vmap(path: Path, rows: Iterable[VMapRow], *, assume_validated: bool = False) -> None:
    rows_list = list(rows)
    if not assume_validated:
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


def write_vmap_table(path: Path, table: VMapRowsTable, *, assume_validated: bool = False) -> None:
    # Assumes: table columns conform to VMapRowsTable schema.
    # Performs: CV(vmap invariants before emission) unless assume_validated, file emission.
    # Guarantees: on-disk .vmap rows satisfy contract-level schema/invariants.
    frame = table.to_frame(copy=False)
    if not assume_validated:
        _validate_vmap_frame(frame)
    ensure_parent_dir(path)
    frame.loc[:, ["chrom", "pos", "id", "a1", "a2", "source_shard", "source_index", "allele_op"]].to_csv(
        path,
        sep="\t",
        header=False,
        index=False,
        lineterminator="\n",
        compression="infer",
        chunksize=WRITE_CHUNK_ROWS,
    )


def _write_series_in_chunks(handle, lines: pd.Series, chunk_rows: int = WRITE_CHUNK_ROWS) -> None:
    # PERF: loop retained for bounded-memory file output; vectorizing full join increases peak memory.
    for start in range(0, len(lines), chunk_rows):
        stop = start + chunk_rows
        handle.write("".join(lines.iloc[start:stop].tolist()))


def _invalid_positive_int_token_mask(values: pd.Series) -> pd.Series:
    # Assumes: `values` is a scalar-like Series expected to encode positive integer coordinates.
    # Performs: vectorized token normalization plus numeric parsing and positive-integer validation.
    # Guarantees: True at indices where token is absent, non-integer, or <= 0.
    tokens = values.str.strip()
    numeric = pd.to_numeric(tokens, errors="coerce")
    valid_positive_tokens = tokens.str.fullmatch(POSITIVE_INT_TOKEN_PATTERN).fillna(False)
    return ~valid_positive_tokens | numeric.isna() | numeric.le(0)


def _validate_vtable_frame(
    frame: pd.DataFrame,
    label: str,
    *,
    assume_normalized_alleles: bool = False,
) -> None:
    if frame["chrom"].eq("").any():
        raise ValueError(f"{label} row is missing chrom")
    invalid_pos = _invalid_positive_int_token_mask(frame["pos"])
    if invalid_pos.any():
        first_idx = int(invalid_pos.idxmax())
        raise ValueError(f"{label} row has invalid pos: {frame.at[first_idx, 'pos']!r}")
    if frame["id"].eq("").any():
        raise ValueError(f"{label} row is missing id")
    validate_allele_values(
        pd.concat([frame["a1"], frame["a2"]], ignore_index=True),
        label=label,
        assume_normalized=assume_normalized_alleles,
    )


def _validate_vmap_frame(
    frame: pd.DataFrame,
    *,
    assume_normalized_alleles: bool = False,
) -> None:
    validate_allele_values(
        pd.concat([frame["a1"], frame["a2"]], ignore_index=True),
        label="vmap",
        assume_normalized=assume_normalized_alleles,
    )
    invalid_pos = _invalid_positive_int_token_mask(frame["pos"])
    if invalid_pos.any():
        first_idx = int(invalid_pos.idxmax())
        raise ValueError(f"vmap row has invalid pos: {frame.at[first_idx, 'pos']!r}")
    duplicate_mask = frame.duplicated(subset=["chrom", "pos", "a1", "a2"], keep="first")
    if duplicate_mask.any():
        raise ValueError("duplicate chrom:pos:a1:a2 in vmap target rows")
    missing_match = frame["source_index"].eq(-1)
    bad_missing_shard = missing_match & frame["source_shard"].ne(MISSING_SOURCE_SHARD)
    if bad_missing_shard.any():
        raise ValueError("vmap row with source_index=-1 must have source_shard='.'")
    bad_missing_op = missing_match & frame["source_shard"].eq(MISSING_SOURCE_SHARD) & frame["allele_op"].ne("missing")
    if bad_missing_op.any():
        raise ValueError("vmap row with source_index=-1 must have allele_op=missing")
    out_of_range = frame["source_index"].lt(0) & ~missing_match
    if out_of_range.any():
        raise ValueError("vmap row source_index out of range")
    required_shard = frame["source_index"].ge(0) & frame["source_shard"].eq("")
    if required_shard.any():
        raise ValueError("vmap row with source_index>=0 must define source_shard")
    valid_nonmissing_ops = VALID_ALLELE_OPS - {"missing"}
    invalid_op = frame["source_index"].ge(0) & frame["source_shard"].ne("") & ~frame["allele_op"].isin(valid_nonmissing_ops)
    if invalid_op.any():
        first_idx = int(invalid_op.idxmax())
        raise ValueError(f"invalid allele_op in vmap row: {frame.at[first_idx, 'allele_op']!r}")


def validate_vtable_rows(rows: Sequence[VariantRow], label: str = "vtable") -> None:
    _validate_vtable_frame(
        VariantRowsTable.from_rows(list(rows), keep_row_idx=False).to_frame(copy=False),
        label,
    )


def validate_vmap_rows(rows: Sequence[VMapRow]) -> None:
    _validate_vmap_frame(VMapRowsTable.from_rows(list(rows), keep_row_idx=False).to_frame(copy=False))


def variant_rows_from_vmap_rows(rows: Sequence[VMapRow]) -> List[VariantRow]:
    return [VariantRow(row.chrom, row.pos, row.id, row.a1, row.a2) for row in rows]



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


def load_variant_object_tables(path: Path) -> LoadedVariantObjectTables:
    # Assumes: metadata sidecar exists and matches object suffix.
    # Performs: PV(object-type/metadata boundary checks), PN/PV via table readers.
    # Guarantees: typed table-first loaded object with target metadata and optional provenance table.
    metadata = load_metadata(path)
    if path.name.endswith(".vtable"):
        validate_vtable_metadata(metadata)
        return LoadedVariantObjectTables(
            path=path,
            object_type="variant_table",
            target_rows_table=read_vtable_table(path),
            target_metadata={
                "genome_build": metadata["genome_build"],
                **({"contig_naming": metadata["contig_naming"]} if "contig_naming" in metadata else {}),
            },
            raw_metadata=metadata,
            base_vmap_table=None,
        )
    if path.name.endswith(".vmap"):
        validate_vmap_metadata(metadata)
        vmap_table = read_vmap_table(path)
        target_frame = vmap_table.to_frame(copy=False).loc[:, ["chrom", "pos", "id", "a1", "a2"]]
        return LoadedVariantObjectTables(
            path=path,
            object_type="variant_map",
            target_rows_table=VariantRowsTable.from_frame(target_frame, copy=False),
            target_metadata=dict(metadata["target"]),
            raw_metadata=metadata,
            base_vmap_table=vmap_table,
        )
    raise ValueError(f"unsupported variant object path: {path}")


def compose_allele_ops(first: str, second: str) -> str:
    if first == "missing" or second == "missing":
        return "missing"
    composed = ALLELE_OP_COMPOSITION.get((first, second))
    if composed is None:
        raise ValueError(f"cannot compose allele operations: {first!r}, {second!r}")
    return composed


def compose_allele_ops_series(first_ops: pd.Series, second_ops: pd.Series) -> pd.Series:
    # Assumes: both inputs are position-aligned and encode allele-op tokens.
    # Performs: vectorized pairwise allele-op composition via unique-pair remap.
    # Guarantees: one composed allele-op token per input row; output index matches first_ops.index.
    if len(first_ops) != len(second_ops):
        raise ValueError("compose_allele_ops_series requires aligned series lengths")
    original_index = first_ops.index
    first = first_ops.reset_index(drop=True).astype(str)
    second = second_ops.reset_index(drop=True).astype(str)
    pairs = pd.DataFrame({"first": first, "second": second}, dtype="object")
    unique_pairs = pairs.drop_duplicates(ignore_index=True)
    composed_values: List[str] = []
    # PERF: loop retained over tiny unique-op pair set (bounded by allele-op domain), not row count.
    for first_op, second_op in unique_pairs.itertuples(index=False, name=None):
        composed_values.append(compose_allele_ops(str(first_op), str(second_op)))
    unique_pairs["composed"] = pd.Series(composed_values, dtype="object")
    result = pairs.merge(unique_pairs, on=["first", "second"], how="left", sort=False)["composed"]
    return result.set_axis(original_index)


def duplicate_target_rows_mask_table(frame: pd.DataFrame) -> pd.Series:
    # Assumes: frame represents target-row columns for variant rows/maps.
    # Performs: vectorized duplicate detection on canonical target key.
    # Guarantees: boolean mask marking all duplicate target rows.
    required = {"chrom", "pos", "a1", "a2"}
    missing = [col for col in required if col not in frame.columns]
    if missing:
        missing_csv = ",".join(missing)
        raise ValueError(f"target frame is missing required columns: {missing_csv}")
    return frame.duplicated(subset=["chrom", "pos", "a1", "a2"], keep=False)


def sort_target_table_by_declared_coordinate(frame: pd.DataFrame, contig_naming: str, *, label: str) -> pd.DataFrame:
    # Assumes: target frame has at least chrom/pos columns in declared contig naming.
    # Performs: PV/SV checks for sortable keys and deterministic stable sort.
    # Guarantees: frame sorted by canonical contig rank then integer position.
    required = {"chrom", "pos"}
    missing = [col for col in required if col not in frame.columns]
    if missing:
        missing_csv = ",".join(missing)
        raise ValueError(f"{label} frame is missing required columns: {missing_csv}")
    out = frame.copy()
    canonical = out["chrom"].map(lambda chrom: canonical_contig_from_label(str(chrom), contig_naming))
    invalid_mask = canonical.isna() | canonical.eq(UNKNOWN_CONTIG)
    if invalid_mask.any():
        first_invalid = str(out.loc[invalid_mask, "chrom"].iloc[0])
        raise ValueError(
            f"{label} has contigs inconsistent with declared contig_naming={contig_naming!r}: {first_invalid!r}. "
            "Run normalize_contigs.py first"
        )
    pos_numeric = pd.to_numeric(out["pos"], errors="coerce")
    invalid_pos = pos_numeric.isna() | (pos_numeric <= 0)
    if invalid_pos.any():
        first_pos = out.loc[invalid_pos, "pos"].iloc[0]
        raise ValueError(f"{label} row has invalid pos: {first_pos!r}")
    out["__rank"] = canonical.map(lambda token: CANONICAL_CONTIG_RANK[str(token)]).astype("int64")
    out["__pos"] = pos_numeric.astype("int64")
    out.sort_values(by=["__rank", "__pos"], kind="mergesort", inplace=True)
    out.drop(columns=["__rank", "__pos"], inplace=True)
    out.reset_index(drop=True, inplace=True)
    return out


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
    token = normalize_allele_token(value)
    if not token:
        raise ValueError(f"invalid allele code in {label}: {value!r}")
    if any(base not in COMPLEMENT for base in token):
        raise ValueError(f"invalid allele code in {label}: {value!r}")


def validate_allele_values(
    values: Sequence[str] | pd.Series,
    *,
    label: str,
    assume_normalized: bool = False,
) -> None:
    if isinstance(values, pd.Series):
        raw = values.reset_index(drop=True)
    else:
        raw = pd.Series(list(values), dtype="object")
    if assume_normalized:
        normalized = raw.astype(str)
    else:
        normalized = raw.astype(str).str.strip().str.upper()
    empty_mask = normalized.isna() | normalized.eq("")
    if bool(empty_mask.any()):
        first_empty = int(empty_mask.idxmax())
        raise ValueError(f"invalid allele code in {label}: {raw.at[first_empty]!r}")
    valid_mask = normalized.str.fullmatch(r"[ACGT]+").fillna(False)
    invalid_mask = ~valid_mask
    if bool(invalid_mask.any()):
        first_invalid = int(invalid_mask.idxmax())
        raise ValueError(f"invalid allele code in {label}: {raw.at[first_invalid]!r}")


def validate_snv_alleles(row: VariantRow, *, label: str) -> None:
    for allele_name, allele in (("a1", row.a1), ("a2", row.a2)):
        token = normalize_allele_token(allele)
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


def require_table_matches_contig_naming(table: VariantRowsTable, contig_naming: str, *, label: str) -> None:
    frame = table.to_frame(copy=False)
    canonical = frame["chrom"].map(lambda chrom: canonical_contig_from_label(str(chrom), contig_naming))
    invalid_mask = canonical.isna() | canonical.eq(UNKNOWN_CONTIG)
    if invalid_mask.any():
        examples = frame.loc[invalid_mask, "chrom"].astype(str).head(3).tolist()
        examples_str = ", ".join(repr(c) for c in examples)
        raise ValueError(
            f"{label} has contigs inconsistent with declared contig_naming={contig_naming!r}: {examples_str}. "
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
    token = normalize_allele_token(value)
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
    source = (source_a1, source_a2)
    target = (target_a1, target_a2)
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
