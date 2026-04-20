from __future__ import annotations

from collections.abc import Callable, Sequence

import pandas as pd


def require_columns(frame: pd.DataFrame, required: Sequence[str], *, label: str) -> None:
    # Assumes: caller declares stage-required columns.
    # Performs: PV(required-column structural check).
    # Guarantees: all required columns exist, or deterministic missing-column error.
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} frame is missing required columns: {','.join(missing)}")


def map_unique_values(series: pd.Series, mapper: Callable[[str], str]) -> pd.Series:
    # Assumes: mapper is deterministic for equal string inputs.
    # Performs: SV(unique-value map/remap optimization for expensive per-label transforms).
    # Guarantees: mapped output preserving input index order.
    # PERF: loop is over unique values (small cardinality relative to row count).
    mapping: dict[str, str] = {}
    for value in pd.unique(series):
        key = str(value)
        mapping[key] = mapper(key)
    return series.map(lambda value: mapping[str(value)])


def strict_int_series(series: pd.Series) -> tuple[pd.Series, pd.Series]:
    # Assumes: caller expects integer-like textual tokens.
    # Performs: PN(text canonicalization) + PV(strict integer parseability mask).
    # Guarantees: Int64 output series plus validity mask aligned to input index.
    text = series.astype(str).str.strip()
    valid = text.str.fullmatch(r"[+-]?\d+")
    out = pd.Series(pd.NA, index=series.index, dtype="Int64")
    out.loc[valid] = text.loc[valid].astype("int64")
    return out, valid


def first_true_index(mask: pd.Series) -> int:
    # Assumes: caller checked that mask has at least one True entry.
    # Performs: deterministic first-failing-index extraction.
    # Guarantees: integer index of first True value under pandas idxmax semantics.
    return int(mask.idxmax())


