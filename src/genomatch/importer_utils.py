from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

import pandas as pd

_ACGT_RE = re.compile(r"^[ACGT]+$")

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
    normalize_allele_token,
    VMapRowsTable,
    VariantRow,
    infer_contig_naming,
    make_vmap_metadata,
    parse_chr2use,
    write_metadata,
    write_vmap_table,
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


IMPORTED_ROWS_COLUMNS = ["chrom", "pos", "id", "a1", "a2", "source_shard", "source_index"]
IMPORT_QC_COLUMNS = ["source_shard", "source_index", "reason"]


def _require_columns(frame: pd.DataFrame, required: Sequence[str], *, label: str) -> None:
    missing = [name for name in required if name not in frame.columns]
    if missing:
        missing_csv = ",".join(missing)
        raise ValueError(f"{label} frame is missing required columns: {missing_csv}")


def _rows_frame_from_dataclass_rows(rows: Iterable[ImportedVariantRow]) -> pd.DataFrame:
    rows_list = list(rows)
    frame = pd.DataFrame(
        {
            "chrom": pd.array([item.row.chrom for item in rows_list], dtype=object),
            "pos": pd.array([item.row.pos for item in rows_list], dtype=object),
            "id": pd.array([item.row.id for item in rows_list], dtype=object),
            "a1": pd.array([item.row.a1 for item in rows_list], dtype=object),
            "a2": pd.array([item.row.a2 for item in rows_list], dtype=object),
            "source_shard": pd.array([item.source_shard for item in rows_list], dtype=object),
            "source_index": [item.source_index for item in rows_list],
        },
        columns=IMPORTED_ROWS_COLUMNS,
    )
    return frame


def _qc_frame_from_dataclass_rows(qc_rows: Sequence[ImportQcRow]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_shard": [item.source_shard for item in qc_rows],
            "source_index": [item.source_index for item in qc_rows],
            "reason": [item.reason for item in qc_rows],
        },
        columns=IMPORT_QC_COLUMNS,
    )


def _deduplicate_target_identity_first_wins(rows_view: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Assumes: rows_view already contains IMPORTED_ROWS_COLUMNS in importer discovery order.
    # Performs: deterministic first-occurrence deduplication on target identity (chrom,pos,a1,a2).
    # Guarantees: retained rows preserve first-seen order; dropped duplicates are reported as import QC rows.
    duplicate_mask = rows_view.duplicated(subset=["chrom", "pos", "a1", "a2"], keep="first")
    if not bool(duplicate_mask.any()):
        return rows_view.reset_index(drop=True), pd.DataFrame(columns=IMPORT_QC_COLUMNS)
    dropped = rows_view.loc[duplicate_mask, ["source_shard", "source_index"]].copy()
    dropped["reason"] = "duplicate_target"
    retained = rows_view.loc[~duplicate_mask].reset_index(drop=True)
    return retained, dropped.loc[:, IMPORT_QC_COLUMNS].reset_index(drop=True)


def finalize_imported_vmap_vectorized(
    *,
    output_path: Path,
    rows_frame: pd.DataFrame,
    genome_build: str,
    target_contig_naming: str | None = None,
    infer_target_contig_naming: bool = True,
    created_by: str,
    derived_from: Path,
    qc_rows_frame: pd.DataFrame | None = None,
) -> None:
    # Assumes: `rows_frame` was produced by trusted importer-stage PN/PV and contains
    # canonical row payload columns defined by `IMPORTED_ROWS_COLUMNS`.
    # Performs: CV(required-column presence, source_index numeric->int64 coercion,
    # vmap row-table schema at write boundary, metadata schema via writer helpers).
    # Guarantees: emits `.vmap` + metadata from validated row payload and, when present,
    # emits QC TSV with required `IMPORT_QC_COLUMNS` schema.
    _require_columns(rows_frame, IMPORTED_ROWS_COLUMNS, label="imported rows")
    rows_view = rows_frame.loc[:, IMPORTED_ROWS_COLUMNS].copy(deep=False)
    rows_view, duplicate_qc = _deduplicate_target_identity_first_wins(rows_view)
    chroms = rows_view["chrom"].astype(str).tolist()
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

    vmap_frame = rows_view.assign(allele_op="identity")
    vmap_frame["source_index"] = pd.to_numeric(vmap_frame["source_index"], errors="raise").astype("int64")
    write_vmap_table(
        output_path,
        VMapRowsTable.from_frame(vmap_frame.loc[:, [
            "chrom",
            "pos",
            "id",
            "a1",
            "a2",
            "source_shard",
            "source_index",
            "allele_op",
        ]], copy=False),
        assume_validated=True,
    )
    metadata = make_vmap_metadata(
        {"genome_build": genome_build, "contig_naming": contig_naming},
        provenance={"created_by": created_by, "derived_from": str(derived_from)},
    )
    write_metadata(output_path, metadata)

    if qc_rows_frame is None:
        qc_view = pd.DataFrame(columns=IMPORT_QC_COLUMNS)
    else:
        _require_columns(qc_rows_frame, IMPORT_QC_COLUMNS, label="import qc")
        qc_view = qc_rows_frame.loc[:, IMPORT_QC_COLUMNS]
    if not duplicate_qc.empty:
        qc_view = pd.concat([qc_view, duplicate_qc], ignore_index=True)
    if qc_view.empty:
        return
    write_import_qc(output_path.with_name(output_path.name + ".qc.tsv"), qc_view)


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
    finalize_imported_vmap_vectorized(
        output_path=output_path,
        rows_frame=_rows_frame_from_dataclass_rows(rows),
        genome_build=genome_build,
        target_contig_naming=target_contig_naming,
        infer_target_contig_naming=infer_target_contig_naming,
        created_by=created_by,
        derived_from=derived_from,
        qc_rows_frame=_qc_frame_from_dataclass_rows(qc_rows),
    )


def write_import_qc(path: Path, rows: Sequence[ImportQcRow] | pd.DataFrame) -> None:
    # Assumes: caller passes dataclass QC rows or DataFrame intended for import QC output.
    # Performs: CV(QC schema column presence when DataFrame input, deterministic row order
    # preservation, newline-delimited tabular serialization).
    # Guarantees: writes UTF-8 QC TSV with header `source_shard/source_index/reason`.
    ensure_parent_dir(path)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("source_shard\tsource_index\treason\n")
        if isinstance(rows, pd.DataFrame):
            _require_columns(rows, IMPORT_QC_COLUMNS, label="import qc")
            # PERF: row loop retained for streamed text emission; vectorization would still
            # require per-row string formatting and does not materially reduce I/O-bound cost.
            for source_shard, source_index, reason in rows.loc[:, IMPORT_QC_COLUMNS].itertuples(index=False, name=None):
                handle.write(f"{source_shard}\t{source_index}\t{reason}\n")
            return
        for row in rows:
            handle.write(f"{row.source_shard}\t{row.source_index}\t{row.reason}\n")


def resolve_import_input_paths(
    input_arg: str,
    *,
    kind_label: str,
    explicit_shards_csv: str | None = None,
) -> List[DiscoveredInputShard]:
    input_path = Path(input_arg)
    if "@" not in input_arg:
        if explicit_shards_csv is not None:
            raise ValueError("--shards requires --input to contain '@'")
        if not input_path.exists():
            raise ValueError(f"input {kind_label} not found: {input_path}")
        return [DiscoveredInputShard(path=input_path, source_shard=MISSING_SOURCE_SHARD)]
    if explicit_shards_csv is not None:
        tokens = explicit_shards_csv.split(",")
        if any(token == "" for token in tokens):
            raise ValueError("--shards contains an empty token")
        if len(set(tokens)) != len(tokens):
            raise ValueError("--shards contains duplicate tokens")
        explicit: List[DiscoveredInputShard] = []
        for token in tokens:
            candidate = Path(str(input_path).replace("@", token))
            if not candidate.is_file():
                raise ValueError(f"requested shard not found for token {token!r}: {candidate}")
            explicit.append(DiscoveredInputShard(path=candidate, source_shard=token))
        return explicit
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


def is_canonical_allele_token(token: str) -> bool:
    return bool(token) and bool(_ACGT_RE.match(token))


def is_canonical_import_allele(value: str) -> bool:
    return is_canonical_allele_token(normalize_allele_token(value))


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
