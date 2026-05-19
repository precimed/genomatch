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


def refine_delimiter_from_data_row(
    header_line: str,
    header: List[str],
    delimiter: Optional[str],
    data_lines: List[str],
) -> Tuple[List[str], Optional[str]]:
    if delimiter != "\t" or not data_lines:
        return header, delimiter
    tab_field_counts = [len(split_line(data_line, delimiter)) for data_line in data_lines]
    if all(count == len(header) for count in tab_field_counts):
        return header, delimiter
    whitespace_header = split_line(header_line, None)
    whitespace_field_counts = [len(split_line(data_line, None)) for data_line in data_lines]
    if len(whitespace_header) == len(header) and all(count == len(header) for count in whitespace_field_counts):
        return whitespace_header, None
    return header, delimiter


def load_metadata(path: Path) -> Dict[str, object]:
    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise ValueError("PyYAML is required to read summary-stat metadata YAML") from exc
    data = yaml.safe_load(path.read_text(encoding="utf-8", errors="replace"))
    if not isinstance(data, dict):
        raise ValueError("metadata YAML must be a mapping")
    return data


def resolve_sumstats_input_path(
    input_path_raw: str | None,
    *,
    metadata_path: Path,
    metadata: Dict[str, object],
    consumer_label: str,
) -> Path:
    if input_path_raw:
        return Path(input_path_raw)
    value = find_metadata_value(metadata, "path_sumStats")
    if value is None:
        raise ValueError(f"{consumer_label} requires --input or metadata path_sumStats")
    if not isinstance(value, str):
        raise ValueError("metadata path_sumStats must be a string filename")
    filename = value.strip()
    if not filename:
        raise ValueError("metadata path_sumStats must be a non-empty filename")
    if "/" in filename or "\\" in filename:
        raise ValueError("metadata path_sumStats must be a filename without '/'")
    if "@" in filename:
        raise ValueError("metadata path_sumStats must not contain '@'")
    return metadata_path.parent / filename


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
    chr: Optional[int]
    pos: Optional[int]
    snp: Optional[int]
    effect_allele: int
    other_allele: int


@dataclass(frozen=True)
class EffectColumnMapping:
    signed: List[Optional[int]]
    invert: List[Optional[int]]
    frequency: List[Optional[int]]


MISSING_VALUE_TOKENS: frozenset = frozenset({"", ".", "na", "n/a", "nan", "none", "null"})
SUMSTATS_READ_CHUNK_ROWS = 5_000_000
SUMSTATS_DELIMITER_SAMPLE_ROWS = 10


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


def read_sumstats_header(path: Path) -> Tuple[str, List[str], Optional[str], int]:
    with open_sumstats_data(path) as (_handle, header_line, header, delimiter, header_line_number):
        return header_line, header, delimiter, header_line_number


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


def _lookup_header_index(header: List[str], raw: str, *, label: str) -> Optional[int]:
    by_name = {name: idx for idx, name in enumerate(header)}
    by_lower = {name.lower(): idx for idx, name in enumerate(header)}
    if raw in by_name:
        return by_name[raw]
    raw_lower = raw.lower()
    if raw_lower in by_lower:
        return by_lower[raw_lower]

    stripped_matches = [idx for idx, name in enumerate(header) if name.strip() == raw]
    if len(stripped_matches) == 1:
        return stripped_matches[0]
    if len(stripped_matches) > 1:
        raise ValueError(f"column mapping for {label} is ambiguous after trimming header whitespace: {raw}")

    stripped_lower_matches = [idx for idx, name in enumerate(header) if name.strip().lower() == raw_lower]
    if len(stripped_lower_matches) == 1:
        return stripped_lower_matches[0]
    if len(stripped_lower_matches) > 1:
        raise ValueError(f"column mapping for {label} is ambiguous after trimming header whitespace: {raw}")
    return None


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
    idx = _lookup_header_index(header, raw, label=label)
    if idx is not None:
        return idx
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
def open_sumstats_data(path: Path) -> Iterator[Tuple[TextIO, str, List[str], Optional[str], int]]:
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle):
            if not line.strip():
                continue
            # Preserve historical preamble behavior (`##...` and `# ...` are comments), but
            # accept single-`#` tabular headers such as `#chrom\tpos\t...`.
            if line.startswith("##") or line.startswith("# "):
                continue
            header_line = line.rstrip("\n")
            delimiter = detect_delimiter(header_line)
            header = split_line(header_line, delimiter)
            if header_line.startswith("#") and len(header) <= 1:
                continue
            data_lines: List[str] = []
            for data_line in handle:
                if not data_line.strip() or data_line.startswith("#"):
                    continue
                data_lines.append(data_line.rstrip("\n"))
                if len(data_lines) >= SUMSTATS_DELIMITER_SAMPLE_ROWS:
                    break
            header, delimiter = refine_delimiter_from_data_row(
                header_line,
                header,
                delimiter,
                data_lines,
            )
            yield handle, header_line, header, delimiter, line_number
            return
    raise ValueError(f"sumstats file is missing a header line: {path}")


def read_sumstats_table(path: Path, *, usecols: Optional[List[str]] = None) -> SumstatsTable:
    """
    Assumes: path points to a single-file sumstats payload with one header row (after optional comment/blank preamble).
    Performs: PN(parse delimiter/header and tabular rows via pandas as string-safe payload values), PV(header presence and deterministic source_index assignment).
    Guarantees: table rows preserve raw data-row order with source_index=0..N-1 over rows only (excluding comments/blanks/header), while semantic required-column validation is deferred to importer/apply kernels.
    """
    chunks = list(iter_sumstats_table_chunks(path, usecols=usecols))
    if not chunks:
        header_line, header, delimiter, _header_line_number = read_sumstats_header(path)
        frame = pd.DataFrame(columns=header if usecols is None else usecols, dtype=object)
        source_index = pd.Series([], dtype="int64")
        return SumstatsTable(
            path=path,
            header_line=header_line,
            header=header,
            delimiter=delimiter,
            frame=frame,
            source_index=source_index,
        )
    frame = pd.concat([chunk.frame for chunk in chunks], ignore_index=True)
    source_index = pd.concat([chunk.source_index for chunk in chunks], ignore_index=True)
    first = chunks[0]
    return SumstatsTable(
        path=path,
        header_line=first.header_line,
        header=first.header,
        delimiter=first.delimiter,
        frame=frame,
        source_index=source_index,
    )


def iter_sumstats_table_chunks(
    path: Path,
    *,
    usecols: Optional[List[str]] = None,
    chunk_rows: int = SUMSTATS_READ_CHUNK_ROWS,
) -> Iterator[SumstatsTable]:
    """
    Assumes: path points to a single-file sumstats payload with one header row (after optional comment/blank preamble).
    Performs: PN(parse delimiter/header and tabular rows via pandas as string-safe payload values), PV(header presence and deterministic source_index assignment).
    Guarantees: yields chunks that preserve raw data-row order with source_index=0..N-1 over rows only (excluding comments/blanks/header), while semantic required-column validation is deferred to importer/apply kernels.
    """
    if chunk_rows <= 0:
        raise ValueError("sumstats chunk_rows must be positive")
    header_line, header, delimiter, header_line_number = read_sumstats_header(path)
    main_read_csv_kwargs = build_sumstats_read_csv_kwargs(
        path,
        delimiter,
        keep_default_na=False,
        dtype=object,
    )
    # Ensure we keep a tabular `#...` header token (e.g. `#chrom`) as a real column name.
    main_read_csv_kwargs["header"] = None
    main_read_csv_kwargs["names"] = header
    main_read_csv_kwargs["skiprows"] = header_line_number + 1
    if usecols is not None:
        main_read_csv_kwargs["usecols"] = usecols
    main_read_csv_kwargs["chunksize"] = chunk_rows
    row_offset = 0
    for frame in pd.read_csv(**main_read_csv_kwargs):
        frame = frame.reset_index(drop=True)
        chunk_len = len(frame)
        source_index = pd.Series(range(row_offset, row_offset + chunk_len), dtype="int64")
        row_offset += chunk_len
        yield SumstatsTable(
            path=path,
            header_line=header_line,
            header=header,
            delimiter=delimiter,
            frame=frame,
            source_index=source_index,
        )
