#!/usr/bin/env python3
from __future__ import annotations

import re
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, TextIO, Tuple

import pandas as pd

from .vtable_utils import open_text


def detect_delimiter(line: str) -> Optional[str]:
    if "\t" in line:
        return "\t"
    if "," in line and line.count(",") > line.count("\t"):
        return ","
    return None


def split_line(line: str, delimiter: Optional[str]) -> List[str]:
    line = line.rstrip("\n")
    return line.split(delimiter) if delimiter else line.split()


def join_line(columns: List[str], delimiter: Optional[str]) -> str:
    return (delimiter or " ").join(columns)


def load_metadata(path: Path) -> Dict[str, object]:
    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise ValueError("PyYAML is required to read summary-stat metadata YAML") from exc
    data = yaml.safe_load(path.read_text(encoding="utf-8", errors="replace"))
    if not isinstance(data, dict):
        raise ValueError("metadata YAML must be a mapping")
    return data


def metadata_key_index(data: Dict[str, object]) -> Dict[str, str]:
    by_lower: Dict[str, str] = {}
    for key in data:
        if not isinstance(key, str):
            continue
        key_lower = key.lower()
        existing = by_lower.get(key_lower)
        if existing is not None and existing != key:
            raise ValueError(
                f"summary-stat metadata has ambiguous top-level keys differing only by case: {existing!r}, {key!r}"
            )
        by_lower[key_lower] = key
    return by_lower


def find_metadata_value(data: Dict[str, object], key: str) -> object | None:
    if key in data:
        return data[key]
    actual_key = metadata_key_index(data).get(key.lower())
    if actual_key is None:
        return None
    return data[actual_key]


JOINED_FIELD_TOKENIZER = re.compile(r"[:_]")
JOINED_FIELD_WITH_DELIMITERS = re.compile(r"([:_])")


def split_joined_variant_value(raw: str) -> List[str]:
    return [token for token in JOINED_FIELD_TOKENIZER.split(raw.strip()) if token]


def extract_variant_field(raw: str, field_name: str) -> str:
    value = raw.strip()
    if not value:
        return value
    if ":" not in value and "_" not in value:
        return value
    parts = split_joined_variant_value(value)
    field_index = {
        "CHR": 0,
        "POS": 1,
        "EffectAllele": 2,
        "OtherAllele": 3,
    }.get(field_name)
    if field_index is None or len(parts) <= field_index:
        return value
    return parts[field_index]


def rewrite_variant_fields(raw: str, replacements: Dict[str, str]) -> str:
    if not replacements:
        return raw
    value = raw.strip()
    if ":" not in value and "_" not in value:
        if len(replacements) == 1:
            return next(iter(replacements.values()))
        raise ValueError("shared summary-stat variant column requires joined values when rewriting multiple fields")
    pieces = JOINED_FIELD_WITH_DELIMITERS.split(value)
    tokens = pieces[0::2]
    delimiters = pieces[1::2]
    field_order = {
        "CHR": 0,
        "POS": 1,
        "EffectAllele": 2,
        "OtherAllele": 3,
    }
    max_index = max(field_order[field_name] for field_name in replacements)
    if len(tokens) <= max_index:
        raise ValueError(f"joined summary-stat variant value is missing fields required for rewrite: {raw!r}")
    for field_name, new_value in replacements.items():
        tokens[field_order[field_name]] = new_value
    rebuilt = tokens[0]
    for delimiter, token in zip(delimiters, tokens[1:]):
        rebuilt += delimiter + token
    return rebuilt


@dataclass(frozen=True)
class VariantColumnMapping:
    chr: int
    pos: Optional[int]
    snp: Optional[int]
    effect_allele: int
    other_allele: int


@dataclass(frozen=True)
class EffectColumnMapping:
    signed: List[Optional[int]]
    invert: List[Optional[int]]
    frequency: List[Optional[int]]


MISSING_VALUE_TOKENS: frozenset = frozenset({"", "nan", "none"})


def is_missing_token_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(MISSING_VALUE_TOKENS)


@dataclass(frozen=True)
class SumstatsTable:
    path: Path
    header_line: str
    header: List[str]
    delimiter: Optional[str]
    frame: pd.DataFrame
    source_index: pd.Series


def build_sumstats_read_csv_kwargs(
    path: Path,
    delimiter: Optional[str],
    *,
    keep_default_na: bool,
    dtype: object | None = None,
) -> Dict[str, object]:
    kwargs: Dict[str, object] = {
        "filepath_or_buffer": path,
        "comment": "#",
        "skip_blank_lines": True,
        "keep_default_na": keep_default_na,
    }
    if dtype is not None:
        kwargs["dtype"] = dtype
    if delimiter is None:
        kwargs["sep"] = r"\s+"
        kwargs["engine"] = "python"
    else:
        kwargs["sep"] = delimiter
    return kwargs


def resolve_column(header: List[str], value: object, label: str, *, required: bool = False) -> Optional[int]:
    if value is None:
        if required:
            raise ValueError(f"missing required column mapping: {label}")
        return None
    if not isinstance(value, str):
        raise ValueError(f"column mapping for {label} must be a string column name")
    raw = value.strip()
    if not raw:
        if required:
            raise ValueError(f"missing required column mapping: {label}")
        return None
    by_name = {name: idx for idx, name in enumerate(header)}
    by_lower = {name.lower(): idx for idx, name in enumerate(header)}
    if raw in by_name:
        return by_name[raw]
    if raw.lower() in by_lower:
        return by_lower[raw.lower()]
    raise ValueError(f"column not found for {label}: {raw}")


def resolve_variant_columns(
    header: List[str],
    metadata: Dict[str, object],
    *,
    require_pos: bool,
) -> VariantColumnMapping:
    return VariantColumnMapping(
        chr=resolve_column(header, find_metadata_value(metadata, "col_CHR"), "col_CHR", required=True),
        pos=resolve_column(header, find_metadata_value(metadata, "col_POS"), "col_POS", required=require_pos),
        snp=resolve_column(header, find_metadata_value(metadata, "col_SNP"), "col_SNP"),
        effect_allele=resolve_column(
            header,
            find_metadata_value(metadata, "col_EffectAllele"),
            "col_EffectAllele",
            required=True,
        ),
        other_allele=resolve_column(
            header,
            find_metadata_value(metadata, "col_OtherAllele"),
            "col_OtherAllele",
            required=True,
        ),
    )


def resolve_effect_columns(header: List[str], metadata: Dict[str, object]) -> EffectColumnMapping:
    return EffectColumnMapping(
        signed=[
            resolve_column(header, find_metadata_value(metadata, "col_BETA"), "col_BETA"),
            resolve_column(header, find_metadata_value(metadata, "col_Z"), "col_Z"),
        ],
        invert=[
            resolve_column(header, find_metadata_value(metadata, "col_OR"), "col_OR"),
            resolve_column(header, find_metadata_value(metadata, "col_ORL95"), "col_ORL95"),
            resolve_column(header, find_metadata_value(metadata, "col_ORU95"), "col_ORU95"),
        ],
        frequency=[
            resolve_column(header, find_metadata_value(metadata, "col_EAF"), "col_EAF"),
            resolve_column(header, find_metadata_value(metadata, "col_OAF"), "col_OAF"),
            resolve_column(header, find_metadata_value(metadata, "col_CaseEAF"), "col_CaseEAF"),
            resolve_column(header, find_metadata_value(metadata, "col_CaseOAF"), "col_CaseOAF"),
            resolve_column(header, find_metadata_value(metadata, "col_ControlEAF"), "col_ControlEAF"),
            resolve_column(header, find_metadata_value(metadata, "col_ControlOAF"), "col_ControlOAF"),
        ],
    )


@contextmanager
def open_sumstats_data(path: Path) -> Iterator[Tuple[TextIO, str, List[str], Optional[str]]]:
    with open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            header_line = line.rstrip("\n")
            delimiter = detect_delimiter(header_line)
            yield handle, header_line, split_line(header_line, delimiter), delimiter
            return
    raise ValueError(f"sumstats file is missing a header line: {path}")


def read_sumstats_table(path: Path) -> SumstatsTable:
    """
    Assumes: path points to a single-file sumstats payload with one header row (after optional comment/blank preamble).
    Performs: PN(parse delimiter/header and tabular rows via pandas as string-safe payload values), PV(header presence and deterministic source_index assignment).
    Guarantees: table rows preserve raw data-row order with source_index=0..N-1 over rows only (excluding comments/blanks/header), while semantic required-column validation is deferred to importer/apply kernels.
    """
    with open_sumstats_data(path) as (_handle, header_line, header, delimiter):
        pass

    main_read_csv_kwargs = build_sumstats_read_csv_kwargs(
        path,
        delimiter,
        keep_default_na=False,
        dtype=object,
    )
    frame = pd.read_csv(**main_read_csv_kwargs)

    source_index = pd.Series(range(len(frame)), dtype="int64")
    return SumstatsTable(
        path=path,
        header_line=header_line,
        header=header,
        delimiter=delimiter,
        frame=frame,
        source_index=source_index,
    )
