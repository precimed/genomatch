from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, List, Sequence

import pandas as pd

"""Typed row/table containers with schema checks only.

Performance-contract role:
- Assumes: caller provides columns aligned to required schemas.
- Performs: PV(required-column presence checks) for table wrappers.
- Guarantees: strongly typed row/table containers preserving caller values.

This module does not own PN canonicalization (for example allele trim/upper or dtype
canonicalization). That ownership belongs to parser/loader boundaries.
"""


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


def _require_columns(frame: pd.DataFrame, required: Sequence[str], *, label: str) -> None:
    # Assumes: caller intends to construct a typed table wrapper for `label`.
    # Performs: PV(required-column structural check).
    # Guarantees: all required columns exist, or raises a deterministic error.
    missing = [name for name in required if name not in frame.columns]
    if missing:
        missing_csv = ",".join(missing)
        raise ValueError(f"{label} frame is missing required columns: {missing_csv}")


@dataclass(frozen=True)
class VariantRowsTable:
    """Vectorized container for VariantRow records."""

    frame: pd.DataFrame

    REQUIRED_COLUMNS: ClassVar[tuple[str, ...]] = ("chrom", "pos", "id", "a1", "a2")
    OPTIONAL_COLUMNS: ClassVar[tuple[str, ...]] = ("row_idx",)

    def __post_init__(self) -> None:
        _require_columns(self.frame, self.REQUIRED_COLUMNS, label="variant")

    @classmethod
    def from_rows(cls, rows: Sequence[VariantRow], *, keep_row_idx: bool = True) -> "VariantRowsTable":
        # Assumes: rows already satisfy caller-owned normalization semantics.
        # Performs: schema-preserving row -> frame conversion.
        # Guarantees: VariantRowsTable with required columns (and optional row_idx), object dtype for all string columns.
        data = {
            "chrom": pd.array([row.chrom for row in rows], dtype="object"),
            "pos": pd.array([row.pos for row in rows], dtype="object"),
            "id": pd.array([row.id for row in rows], dtype="object"),
            "a1": pd.array([row.a1 for row in rows], dtype="object"),
            "a2": pd.array([row.a2 for row in rows], dtype="object"),
        }
        if keep_row_idx:
            data["row_idx"] = pd.array(list(range(len(rows))), dtype="int64")
        return cls(pd.DataFrame(data))

    @classmethod
    def from_frame(cls, frame: pd.DataFrame, *, copy: bool = True) -> "VariantRowsTable":
        # Assumes: input frame is caller-owned and already semantically canonical if needed.
        # Performs: PV(required-column check) via wrapper construction.
        # Guarantees: table wrapper around copied or original frame.
        return cls(frame.copy(deep=True) if copy else frame)

    def to_rows(self) -> List[VariantRow]:
        # Assumes: required schema columns exist.
        # Performs: typed projection to row dataclasses.
        # Guarantees: row list preserving field values (string/int casts as needed for row types).
        rows: List[VariantRow] = []
        for record in self.frame.loc[:, self.REQUIRED_COLUMNS].itertuples(index=False, name=None):
            chrom, pos, row_id, a1, a2 = record
            rows.append(VariantRow(str(chrom), str(pos), str(row_id), str(a1), str(a2)))
        return rows

    def to_frame(self, *, copy: bool = True) -> pd.DataFrame:
        return self.frame.copy(deep=True) if copy else self.frame

    def __len__(self) -> int:
        return len(self.frame)


@dataclass(frozen=True)
class VMapRowsTable:
    """Vectorized container for VMapRow records."""

    frame: pd.DataFrame

    REQUIRED_COLUMNS: ClassVar[tuple[str, ...]] = (
        "chrom",
        "pos",
        "id",
        "a1",
        "a2",
        "source_shard",
        "source_index",
        "allele_op",
    )
    OPTIONAL_COLUMNS: ClassVar[tuple[str, ...]] = ("row_idx",)

    def __post_init__(self) -> None:
        _require_columns(self.frame, self.REQUIRED_COLUMNS, label="vmap")

    @classmethod
    def from_rows(cls, rows: Sequence[VMapRow], *, keep_row_idx: bool = True) -> "VMapRowsTable":
        # Assumes: rows already satisfy caller-owned normalization semantics.
        # Performs: schema-preserving row -> frame conversion.
        # Guarantees: VMapRowsTable with required columns (and optional row_idx), object dtype for string columns, int64 for source_index.
        data = {
            "chrom": pd.array([row.chrom for row in rows], dtype="object"),
            "pos": pd.array([row.pos for row in rows], dtype="object"),
            "id": pd.array([row.id for row in rows], dtype="object"),
            "a1": pd.array([row.a1 for row in rows], dtype="object"),
            "a2": pd.array([row.a2 for row in rows], dtype="object"),
            "source_shard": pd.array([row.source_shard for row in rows], dtype="object"),
            "source_index": pd.array([row.source_index for row in rows], dtype="int64"),
            "allele_op": pd.array([row.allele_op for row in rows], dtype="object"),
        }
        if keep_row_idx:
            data["row_idx"] = pd.array(list(range(len(rows))), dtype="int64")
        return cls(pd.DataFrame(data))

    @classmethod
    def from_frame(cls, frame: pd.DataFrame, *, copy: bool = True) -> "VMapRowsTable":
        # Assumes: input frame is caller-owned and already semantically canonical if needed.
        # Performs: PV(required-column check) via wrapper construction.
        # Guarantees: table wrapper around copied or original frame.
        return cls(frame.copy(deep=True) if copy else frame)

    def to_rows(self) -> List[VMapRow]:
        # Assumes: required schema columns exist.
        # Performs: typed projection to row dataclasses.
        # Guarantees: row list preserving field values (string/int casts as needed for row types).
        rows: List[VMapRow] = []
        for record in self.frame.loc[:, self.REQUIRED_COLUMNS].itertuples(index=False, name=None):
            chrom, pos, row_id, a1, a2, source_shard, source_index, allele_op = record
            rows.append(
                VMapRow(
                    str(chrom),
                    str(pos),
                    str(row_id),
                    str(a1),
                    str(a2),
                    str(source_shard),
                    int(source_index),
                    str(allele_op),
                )
            )
        return rows

    def to_frame(self, *, copy: bool = True) -> pd.DataFrame:
        return self.frame.copy(deep=True) if copy else self.frame

    def __len__(self) -> int:
        return len(self.frame)
