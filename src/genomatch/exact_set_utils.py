#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pandas as pd

from .tabular_rows import VMapRowsTable, VariantRowsTable
from .vtable_utils import (
    iter_vmap_table_chunks,
    iter_vtable_table_chunks,
    load_metadata,
    read_vmap_table,
    read_vtable_table,
    require_contig_naming,
    validate_variants_count_if_present,
    validate_vmap_metadata,
    validate_vtable_metadata,
)

ExactVariantKeys = pd.MultiIndex
TARGET_KEY_COLUMNS = ["chrom", "pos", "a1", "a2"]


@dataclass(frozen=True)
class TargetObjectInfo:
    path: Path
    object_type: str
    target_metadata: dict
    raw_metadata: dict
    variants_count: int | None


def load_target_object_info(path: Path) -> TargetObjectInfo:
    metadata = load_metadata(path)
    if path.name.endswith(".vtable"):
        validate_vtable_metadata(metadata)
        return TargetObjectInfo(
            path=path,
            object_type="variant_table",
            target_metadata={
                "genome_build": metadata["genome_build"],
                **({"contig_naming": metadata["contig_naming"]} if "contig_naming" in metadata else {}),
            },
            raw_metadata=metadata,
            variants_count=metadata.get("variants_count"),
        )
    if path.name.endswith(".vmap"):
        validate_vmap_metadata(metadata)
        return TargetObjectInfo(
            path=path,
            object_type="variant_map",
            target_metadata=dict(metadata["target"]),
            raw_metadata=metadata,
            variants_count=metadata.get("variants_count"),
        )
    raise ValueError(f"unsupported variant object path: {path}")


def format_metadata_summary(infos: Sequence[TargetObjectInfo]) -> str:
    lines = ["input metadata:"]
    for info in infos:
        meta = info.target_metadata
        contig_naming = meta.get("contig_naming", "<missing>")
        lines.append(f"- {info.path}: genome_build={meta['genome_build']}, contig_naming={contig_naming}")
    return "\n".join(lines)


def require_shared_target_metadata(infos: Sequence[TargetObjectInfo]) -> None:
    builds = {str(info.target_metadata["genome_build"]) for info in infos}
    metadata_summary = format_metadata_summary(infos)
    contig_namings = set()
    for info in infos:
        try:
            contig_namings.add(require_contig_naming(info.target_metadata, label="variant object"))
        except ValueError as exc:
            raise ValueError(f"{exc}\n{metadata_summary}") from exc
    if len(builds) != 1:
        raise ValueError(f"all inputs must have the same genome_build\n{metadata_summary}")
    if len(contig_namings) != 1:
        raise ValueError(f"all inputs must have the same contig_naming\n{metadata_summary}")


def vtable_metadata_from_target_info(info: TargetObjectInfo) -> dict:
    if info.object_type == "variant_table":
        return dict(info.raw_metadata)
    return {"object_type": "variant_table", **dict(info.target_metadata)}


def read_target_table(path: Path) -> VariantRowsTable:
    if path.name.endswith(".vtable"):
        return read_vtable_table(path)
    if path.name.endswith(".vmap"):
        vmap_table = read_vmap_table(path)
        target_frame = vmap_table.to_frame(copy=False).loc[:, ["chrom", "pos", "id", "a1", "a2"]]
        return VariantRowsTable.from_frame(target_frame, copy=False)
    raise ValueError(f"unsupported variant object path: {path}")


def read_mapped_vmap_table(path: Path) -> VMapRowsTable:
    table = read_vmap_table(path)
    require_mapped_only_vmap_table(table, path=path)
    return table


def require_mapped_only_vmap_table(table: VMapRowsTable, *, path: Path) -> None:
    frame = table.to_frame(copy=False)
    if frame["source_index"].lt(0).any() or frame["allele_op"].eq("missing").any():
        raise ValueError(f"{path} is not a mapped-only .vmap")


def exact_key_index(frame: pd.DataFrame) -> ExactVariantKeys:
    missing = [column for column in TARGET_KEY_COLUMNS if column not in frame.columns]
    if missing:
        missing_csv = ",".join(missing)
        raise ValueError(f"target frame is missing required columns: {missing_csv}")
    return pd.MultiIndex.from_frame(frame.loc[:, TARGET_KEY_COLUMNS])


def drop_duplicate_target_identities(frame: pd.DataFrame) -> pd.DataFrame:
    duplicate_mask = exact_key_index(frame).duplicated(keep="first")
    if not duplicate_mask.any():
        return frame.reset_index(drop=True)
    return frame.loc[~duplicate_mask].reset_index(drop=True)


def assign_target_derived_ids(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["id"] = (
        out["chrom"].astype(str)
        + ":"
        + out["pos"].astype(str)
        + ":"
        + out["a1"].astype(str)
        + ":"
        + out["a2"].astype(str)
    )
    return out


def empty_key_index() -> ExactVariantKeys:
    return pd.MultiIndex.from_tuples([], names=TARGET_KEY_COLUMNS)


def restrict_frame_to_membership_chunks(driver_frame: pd.DataFrame, restriction_path: Path) -> tuple[pd.DataFrame, int]:
    driver_keys = exact_key_index(driver_frame)
    matched_keys = empty_key_index()
    seen_vmap_keys = empty_key_index()
    row_count = 0
    is_vmap = restriction_path.name.endswith(".vmap")
    if restriction_path.name.endswith(".vtable"):
        chunk_iter = iter_vtable_table_chunks(restriction_path, label="variant table row")
    elif is_vmap:
        chunk_iter = iter_vmap_table_chunks(restriction_path)
    else:
        raise ValueError(f"unsupported variant object path: {restriction_path}")

    for table in chunk_iter:
        frame = table.to_frame(copy=False)
        if is_vmap:
            frame = frame.loc[:, ["chrom", "pos", "id", "a1", "a2"]]
        row_count += len(frame)
        chunk_keys = exact_key_index(frame)
        if is_vmap:
            if chunk_keys.isin(seen_vmap_keys).any():
                raise ValueError("duplicate chrom:pos:a1:a2 in vmap target rows")
            seen_vmap_keys = seen_vmap_keys.union(chunk_keys)
        matched_keys = matched_keys.union(chunk_keys[chunk_keys.isin(driver_keys)])

    retain_mask = driver_keys.isin(matched_keys)
    return driver_frame.loc[retain_mask].reset_index(drop=True), row_count


def restrict_frame_to_membership_frame(driver_frame: pd.DataFrame, membership_frame: pd.DataFrame) -> pd.DataFrame:
    driver_keys = exact_key_index(driver_frame)
    membership_keys = exact_key_index(membership_frame)
    retain_mask = driver_keys.isin(membership_keys)
    return driver_frame.loc[retain_mask].reset_index(drop=True)


def validate_loaded_row_count(info: TargetObjectInfo, observed_count: int) -> None:
    validate_variants_count_if_present(info.path, info.raw_metadata, observed_count)
