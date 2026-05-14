#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from ._cli_utils import run_cli
from .apply_vmap_utils import output_variant_id
from .sumstats_clean import harmonize_clean_sumstats_frame, resolve_clean_metadata_columns
from .sumstats_utils import (
    find_metadata_value,
    join_line,
    load_metadata,
    MISSING_VALUE_TOKENS,
    open_text,
    read_sumstats_table,
    read_sumstats_header,
    resolve_sumstats_input_path,
    _lookup_header_index,
    resolve_column,
    resolve_effect_columns,
    rewrite_variant_fields,
    SumstatsTable,
)
from .vtable_utils import (
    load_metadata as load_variant_metadata,
    require_contig_naming,
    require_table_matches_contig_naming,
    read_vmap_table,
    validate_vmap_metadata,
)


NA_NUMERIC = "n/a"
JOINED_VARIANT_FIELD_ORDER = ["CHR", "POS", "EffectAllele", "OtherAllele"]
logger = logging.getLogger(__name__)


def warn_once(warning_keys: Set[Tuple[str, str]], key: Tuple[str, str], message: str) -> None:
    if key in warning_keys:
        return
    warning_keys.add(key)
    logger.warning("%s", message)


def parse_finite_float(raw: str) -> float:
    value = float(raw)
    if not math.isfinite(value):
        raise ValueError(f"non-finite numeric value: {raw!r}")
    return value


def maybe_negate(raw: str, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> str:
    try:
        return str(-parse_finite_float(raw))
    except Exception:
        warn_once(
            warning_keys,
            ("negate", column_name),
            f"could not negate non-numeric value in column {column_name!r}; writing {NA_NUMERIC}",
        )
        return NA_NUMERIC


def maybe_invert(raw: str, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> str:
    try:
        value = parse_finite_float(raw)
        if value == 0:
            raise ValueError("cannot invert zero")
        return str(1.0 / value)
    except Exception:
        warn_once(
            warning_keys,
            ("invert", column_name),
            f"could not invert zero or non-numeric value in column {column_name!r}; writing {NA_NUMERIC}",
        )
        return NA_NUMERIC


def maybe_complement(raw: str, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> str:
    try:
        value = parse_finite_float(raw)
        return str(1.0 - value)
    except Exception:
        warn_once(
            warning_keys,
            ("complement", column_name),
            f"could not complement non-numeric value in column {column_name!r}; writing {NA_NUMERIC}",
        )
        return NA_NUMERIC


def maybe_invert_interval(
    lower_raw: str,
    upper_raw: str,
    *,
    lower_column: str,
    upper_column: str,
    warning_keys: Set[Tuple[str, str]],
) -> Tuple[str, str]:
    try:
        lower = parse_finite_float(lower_raw)
        upper = parse_finite_float(upper_raw)
        if lower == 0 or upper == 0:
            raise ValueError("cannot invert zero")
        return str(1.0 / upper), str(1.0 / lower)
    except Exception:
        warn_once(
            warning_keys,
            ("invert_interval", f"{lower_column}|{upper_column}"),
            "could not invert zero or non-numeric odds-ratio interval in columns "
            f"{lower_column!r} and {upper_column!r}; writing {NA_NUMERIC}",
        )
        return NA_NUMERIC, NA_NUMERIC


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply a .vmap to summary statistics.")
    parser.add_argument("--input", help="Input summary statistics file (optional when metadata defines path_sumStats)")
    parser.add_argument("--sumstats-metadata", required=True, help="Cleansumstats-style metadata YAML")
    parser.add_argument("--vmap", required=True, help="Input .vmap")
    parser.add_argument("--output", required=True, help="Output summary statistics file")
    parser.add_argument("--clean", action="store_true", help="Emit canonical cleaned summary statistics")
    parser.add_argument("--fill-mode", choices=["column", "row"], default="column", help="Clean-mode fill behavior")
    parser.add_argument("--use-af-inference", action="store_true", help="Enable AF-based clean-mode derivation rules")
    parser.add_argument(
        "--only-mapped-target",
        action="store_true",
        help="Drop target rows with source_index=-1 instead of writing unmatched target rows with missing payload values",
    )
    parser.add_argument(
        "--retain-snp-id",
        action="store_true",
        help="Use retained target-side .vmap id values as output SNP IDs instead of generated chrom:pos:a1:a2 IDs",
    )
    return parser.parse_args()


def metadata_column_name(metadata: Dict[str, object], key: str) -> Optional[str]:
    value = find_metadata_value(metadata, key)
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"column mapping for {key} must be a string column name")
    raw = value.strip()
    return raw or None


def output_column_index(header: List[str], column_name: Optional[str]) -> Optional[int]:
    if column_name is None:
        return None
    return _lookup_header_index(header, column_name, label=column_name)


def load_sumstats_rows_single_file(
    sumstats_table: SumstatsTable,
    required_variant_columns: Dict[str, int],
    needed_source_indices: Set[int],
) -> Tuple[str, List[str], Optional[str], Dict[Tuple[str, int], List[str]]]:
    """
    Assumes: `sumstats_table` came from shared PN/PV parse boundary (`read_sumstats_table`) and source provenance has already been validated as single-file.
    Performs: SV required-variant-field presence validation and source-shard/index lookup materialization.
    Guarantees: returns rows keyed by (source_shard='.', source_index) using shared zero-based source_index semantics over parsed data rows.
    """
    frame = sumstats_table.frame
    if required_variant_columns:
        field_missing: Dict[str, bool] = {}
        for field_name, col_idx in required_variant_columns.items():
            if col_idx >= len(frame.columns):
                raise ValueError(
                    f"required variant column {field_name!r} index {col_idx} is out of range "
                    f"for sumstats frame with {len(frame.columns)} columns: {sumstats_table.path}"
                )
            # Allowed redundancy: normalization repeats PN already done by extract_variant_field in
            # import_sumstats, but apply_vmap operates on a different trust boundary (pre-vmap sumstats
            # re-read at apply time, not at import time).
            column_values = frame.iloc[:, col_idx].astype(str).str.strip().str.lower()
            field_missing[field_name] = bool(column_values.isin(MISSING_VALUE_TOKENS).any())
        if any(field_missing.values()):
            missing_fields = [name for name, is_missing in field_missing.items() if is_missing]
            raise ValueError(
                "sumstats input has missing values in required variant columns: "
                + ", ".join(missing_fields)
            )

    rows_by_provenance: Dict[Tuple[str, int], List[str]] = {}
    if needed_source_indices:
        if len(sumstats_table.source_index) != len(frame):
            raise ValueError(f"sumstats source_index length mismatch for parsed frame: {sumstats_table.path}")
        payload_rows_by_source_index = {
            int(source_index): [str(value) for value in row_values]
            for source_index, row_values in zip(
                sumstats_table.source_index.tolist(),
                frame.itertuples(index=False, name=None),
            )
        }
        for row_index in sorted(needed_source_indices):
            row_values = payload_rows_by_source_index.get(row_index)
            if row_values is not None:
                rows_by_provenance[(".", row_index)] = row_values
    return sumstats_table.header_line, sumstats_table.header, sumstats_table.delimiter, rows_by_provenance


def rewrite_or_synthesize_variant_value(raw: str, replacements: Dict[str, str]) -> str:
    value = raw.strip()
    if value and (":" in value or "_" in value):
        return rewrite_variant_fields(value, replacements)
    if len(replacements) == 1:
        return next(iter(replacements.values()))
    tokens = [NA_NUMERIC] * (max(JOINED_VARIANT_FIELD_ORDER.index(name) for name in replacements) + 1)
    for field_name, new_value in replacements.items():
        tokens[JOINED_VARIANT_FIELD_ORDER.index(field_name)] = new_value
    rebuilt = tokens[0]
    for idx, token in enumerate(tokens[1:], start=1):
        rebuilt += ("_" if idx == 1 else ":") + token
    return rebuilt


def build_missing_payload_row(
    output_width: int,
    output_header: List[str],
    vrow,
    idx_chr: int,
    idx_pos: Optional[int],
    idx_id: Optional[int],
    idx_a1: Optional[int],
    idx_a2: Optional[int],
    output_idx_pos: Optional[int],
    output_idx_id: Optional[int],
    output_idx_a1: Optional[int],
    output_idx_a2: Optional[int],
    retain_snp_id: bool,
) -> List[str]:
    cols = [NA_NUMERIC] * output_width
    out_id = output_variant_id(vrow, retain_snp_id=retain_snp_id)
    variant_rewrites: Dict[int, Dict[str, str]] = {}
    variant_rewrites.setdefault(idx_chr, {})["CHR"] = vrow.chrom
    if idx_pos is not None:
        variant_rewrites.setdefault(idx_pos, {})["POS"] = vrow.pos
    if idx_id is not None:
        cols[idx_id] = out_id
    if idx_a1 is not None:
        variant_rewrites.setdefault(idx_a1, {})["EffectAllele"] = vrow.a1
    if idx_a2 is not None:
        variant_rewrites.setdefault(idx_a2, {})["OtherAllele"] = vrow.a2
    for idx, replacements in variant_rewrites.items():
        cols[idx] = rewrite_or_synthesize_variant_value(cols[idx], replacements)
    if idx_pos is None and output_idx_pos is not None:
        cols[output_idx_pos] = vrow.pos
    if idx_id is None and output_idx_id is not None:
        cols[output_idx_id] = out_id
    if idx_a1 is None and output_idx_a1 is not None:
        cols[output_idx_a1] = vrow.a1
    if idx_a2 is None and output_idx_a2 is not None:
        cols[output_idx_a2] = vrow.a2
    return cols


def maybe_negate_clean_series(series: pd.Series, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce").astype(float)
    invalid_mask = ~np.isfinite(numeric)
    if invalid_mask.any():
        warn_once(
            warning_keys,
            ("negate_clean", column_name),
            f"could not negate non-numeric value in column {column_name!r}; writing missing value",
        )
    numeric = -numeric
    numeric[invalid_mask] = np.nan
    return numeric


def maybe_invert_clean_series(series: pd.Series, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce").astype(float)
    invalid_mask = ~np.isfinite(numeric) | (numeric == 0)
    if invalid_mask.any():
        warn_once(
            warning_keys,
            ("invert_clean", column_name),
            f"could not invert zero or non-numeric value in column {column_name!r}; writing missing value",
        )
    numeric = 1.0 / numeric
    numeric[invalid_mask] = np.nan
    return numeric


def maybe_complement_clean_series(series: pd.Series, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce").astype(float)
    invalid_mask = ~np.isfinite(numeric)
    if invalid_mask.any():
        warn_once(
            warning_keys,
            ("complement_clean", column_name),
            f"could not complement non-numeric value in column {column_name!r}; writing missing value",
        )
    numeric = 1.0 - numeric
    numeric[invalid_mask] = np.nan
    return numeric


def maybe_invert_interval_clean_series(
    lower_series: pd.Series,
    upper_series: pd.Series,
    *,
    lower_column: str,
    upper_column: str,
    warning_keys: Set[Tuple[str, str]],
) -> Tuple[pd.Series, pd.Series]:
    lower = pd.to_numeric(lower_series, errors="coerce").astype(float)
    upper = pd.to_numeric(upper_series, errors="coerce").astype(float)
    invalid_mask = ~np.isfinite(lower) | ~np.isfinite(upper) | (lower == 0) | (upper == 0)
    if invalid_mask.any():
        warn_once(
            warning_keys,
            ("invert_interval_clean", f"{lower_column}|{upper_column}"),
            "could not invert zero or non-numeric odds-ratio interval in columns "
            f"{lower_column!r} and {upper_column!r}; writing missing value",
        )
    lower_out = 1.0 / upper
    upper_out = 1.0 / lower
    lower_out[invalid_mask] = np.nan
    upper_out[invalid_mask] = np.nan
    return lower_out, upper_out


def _unique_names_in_header_order(header: Sequence[str], indices: Sequence[int]) -> List[str]:
    names: List[str] = []
    seen: Set[str] = set()
    for idx in sorted(set(indices)):
        name = header[idx]
        if name not in seen:
            names.append(name)
            seen.add(name)
    return names


def clean_sumstats_usecols(preview_header: Sequence[str], metadata: Dict[str, object]) -> List[str]:
    resolved = resolve_clean_metadata_columns(preview_header, metadata, include_variant_columns=True)
    return _unique_names_in_header_order(preview_header, list(resolved.values()))


def clean_variant_column_names(preview_header: Sequence[str], metadata: Dict[str, object]) -> List[str]:
    resolved = resolve_clean_metadata_columns(preview_header, metadata, include_variant_columns=True)
    variant_indices = [
        resolved[key]
        for key in ("col_CHR", "col_POS", "col_SNP", "col_EffectAllele", "col_OtherAllele")
        if key in resolved
    ]
    return _unique_names_in_header_order(preview_header, variant_indices)


def filtered_vmap_frame(vmap_frame: pd.DataFrame, *, only_mapped_target: bool) -> pd.DataFrame:
    if not only_mapped_target:
        return vmap_frame.reset_index(drop=True)
    return vmap_frame.loc[vmap_frame["source_index"].ne(-1)].reset_index(drop=True)


def require_single_file_sumstats_provenance(vmap_frame: pd.DataFrame) -> pd.Series:
    mapped_frame = vmap_frame.loc[vmap_frame["source_index"].ne(-1), ["source_shard", "source_index"]]
    if mapped_frame.empty:
        return mapped_frame["source_index"].astype("int64")
    invalid_shard_mask = mapped_frame["source_shard"].ne(".")
    if invalid_shard_mask.any():
        unsupported = sorted(mapped_frame.loc[invalid_shard_mask, "source_shard"].astype(str).unique().tolist())
        raise ValueError(
            "apply_vmap_to_sumstats.py supports single-file payload lookup only for source_shard='.'; "
            f"found {unsupported!r}"
        )
    return mapped_frame["source_index"].astype("int64")


def collect_single_file_source_indices(vmap_frame: pd.DataFrame) -> Set[int]:
    return set(require_single_file_sumstats_provenance(vmap_frame).tolist())


def validate_source_indices_in_range(source_indices: pd.Series, row_count: int) -> None:
    if source_indices.empty:
        return
    invalid_mask = source_indices.lt(0) | source_indices.ge(row_count)
    if invalid_mask.any():
        first_missing = int(source_indices.loc[invalid_mask].iloc[0])
        raise ValueError(
            f"vmap source provenance out of range for summary statistics input: "
            f"first missing source_shard='.', source_index={first_missing}"
        )


def validate_required_variant_columns(
    frame: pd.DataFrame,
    preview_header: Sequence[str],
    required_variant_columns: Dict[str, int],
    *,
    path: Path,
) -> None:
    field_missing: Dict[str, bool] = {}
    for field_name, col_idx in required_variant_columns.items():
        if col_idx >= len(preview_header):
            raise ValueError(
                f"required variant column {field_name!r} index {col_idx} is out of range "
                f"for sumstats header with {len(preview_header)} columns: {path}"
            )
        column_name = preview_header[col_idx]
        if column_name not in frame.columns:
            raise ValueError(f"required variant column {field_name!r} was not loaded from {path}: {column_name!r}")
        column_values = frame[column_name].astype(str).str.strip().str.lower()
        field_missing[field_name] = bool(column_values.isin(MISSING_VALUE_TOKENS).any())
    if any(field_missing.values()):
        missing_fields = [name for name, is_missing in field_missing.items() if is_missing]
        raise ValueError(
            "sumstats input has missing values in required variant columns: "
            + ", ".join(missing_fields)
        )


def collect_clean_frame(
    preview_header: Sequence[str],
    metadata: Dict[str, object],
    vmap_frame: pd.DataFrame,
    sumstats_table: SumstatsTable,
) -> Tuple[List[str], pd.DataFrame]:
    resolved = resolve_clean_metadata_columns(preview_header, metadata, include_variant_columns=True)
    variant_indices = {
        resolved[key]
        for key in ("col_CHR", "col_POS", "col_SNP", "col_EffectAllele", "col_OtherAllele")
        if key in resolved
    }
    payload_indices = [
        idx
        for idx in range(len(preview_header))
        if idx not in variant_indices and preview_header[idx] in sumstats_table.frame.columns
    ]
    payload_header = [preview_header[idx] for idx in payload_indices]
    if not payload_header:
        return payload_header, pd.DataFrame(index=range(len(vmap_frame)))

    source_indices = vmap_frame["source_index"].to_numpy(dtype=np.int64, copy=False)
    mapped_mask = source_indices != -1
    payload_frame = pd.DataFrame(index=range(len(vmap_frame)))
    identity_order = (
        mapped_mask.all()
        and len(vmap_frame) == len(sumstats_table.frame)
        and np.array_equal(source_indices, np.arange(len(source_indices)))
    )
    mapped_source_indices = source_indices[mapped_mask]
    for column in payload_header:
        if identity_order:
            payload_frame[column] = sumstats_table.frame[column].reset_index(drop=True)
        elif mapped_mask.all():
            payload_frame[column] = pd.Series(
                sumstats_table.frame[column].iloc[source_indices].to_numpy(copy=False),
                index=payload_frame.index,
                dtype=object,
            )
        else:
            values = pd.Series([None] * len(vmap_frame), index=payload_frame.index, dtype=object)
            values.loc[mapped_mask] = sumstats_table.frame[column].iloc[mapped_source_indices].to_numpy(copy=False)
            payload_frame[column] = values
    return payload_header, payload_frame


def clear_clean_unmatched_rows(payload_sumstats: pd.DataFrame, vmap_frame: pd.DataFrame) -> None:
    unmatched_mask = vmap_frame["source_index"].eq(-1).to_numpy()
    if unmatched_mask.any():
        payload_sumstats.loc[unmatched_mask, :] = float("nan")


def run_clean_apply(
    args: argparse.Namespace,
    *,
    metadata: Dict[str, object],
    preview_header: List[str],
    vmap_frame: pd.DataFrame,
    sumstats_table: SumstatsTable,
    output_path: Path,
) -> int:
    payload_header, payload_frame = collect_clean_frame(preview_header, metadata, vmap_frame, sumstats_table)
    payload_sumstats, _clean_metadata = harmonize_clean_sumstats_frame(
        payload_header,
        payload_frame,
        metadata,
        fill_mode=args.fill_mode,
        use_af_inference=args.use_af_inference,
        warn=lambda message: logger.warning("%s", message),
    )
    clear_clean_unmatched_rows(payload_sumstats, vmap_frame)
    warning_keys: Set[Tuple[str, str]] = set()
    swap_mask = vmap_frame["allele_op"].isin({"swap", "flip_swap"}).to_numpy()

    for column in ("BETA", "Z"):
        if column in payload_sumstats.columns:
            payload_sumstats.loc[swap_mask, column] = maybe_negate_clean_series(
                payload_sumstats.loc[swap_mask, column],
                column_name=column,
                warning_keys=warning_keys,
            )

    if "ORL95" in payload_sumstats.columns and "ORU95" in payload_sumstats.columns:
        lower_series, upper_series = maybe_invert_interval_clean_series(
            payload_sumstats.loc[swap_mask, "ORL95"],
            payload_sumstats.loc[swap_mask, "ORU95"],
            lower_column="ORL95",
            upper_column="ORU95",
            warning_keys=warning_keys,
        )
        payload_sumstats.loc[swap_mask, "ORL95"] = lower_series
        payload_sumstats.loc[swap_mask, "ORU95"] = upper_series
    else:
        for column in ("ORL95", "ORU95"):
            if column in payload_sumstats.columns:
                payload_sumstats.loc[swap_mask, column] = maybe_invert_clean_series(
                    payload_sumstats.loc[swap_mask, column],
                    column_name=column,
                    warning_keys=warning_keys,
                )

    if "OR" in payload_sumstats.columns:
        payload_sumstats.loc[swap_mask, "OR"] = maybe_invert_clean_series(
            payload_sumstats.loc[swap_mask, "OR"],
            column_name="OR",
            warning_keys=warning_keys,
        )

    for column in ("EAF", "CaseEAF", "ControlEAF"):
        if column in payload_sumstats.columns:
            payload_sumstats.loc[swap_mask, column] = maybe_complement_clean_series(
                payload_sumstats.loc[swap_mask, column],
                column_name=column,
                warning_keys=warning_keys,
            )

    retained_mask = np.ones(len(vmap_frame), dtype=bool)
    if args.only_mapped_target:
        if "P" not in payload_sumstats.columns:
            retained_mask[:] = False
        else:
            retained_mask &= ~payload_sumstats["P"].isna().to_numpy()
        if not retained_mask.any():
            raise ValueError("no retained target rows remain after applying --only-mapped-target")

    output_header = ["CHR", "POS", "SNP", "EffectAllele", "OtherAllele", *payload_sumstats.columns]
    chunk_size = 100_000
    with open_text(output_path, "wt") as handle:
        handle.write("\t".join(output_header) + "\n")
        for start in range(0, len(vmap_frame), chunk_size):
            stop = min(start + chunk_size, len(vmap_frame))
            chunk_offsets = np.flatnonzero(retained_mask[start:stop])
            if len(chunk_offsets) == 0:
                continue
            chunk_indices = chunk_offsets + start
            vmap_chunk = vmap_frame.iloc[chunk_indices]
            chunk_df = pd.DataFrame({
                "CHR": vmap_chunk["chrom"].to_numpy(copy=False),
                "POS": vmap_chunk["pos"].to_numpy(copy=False),
                "SNP": [
                    output_variant_id(row, retain_snp_id=args.retain_snp_id)
                    for row in vmap_chunk.itertuples(index=False)
                ],
                "EffectAllele": vmap_chunk["a1"].to_numpy(copy=False),
                "OtherAllele": vmap_chunk["a2"].to_numpy(copy=False),
            })
            for col in payload_sumstats.columns:
                chunk_df[col] = payload_sumstats[col].iloc[chunk_indices].to_numpy(copy=False)
            chunk_df.to_csv(
                handle,
                sep="\t",
                header=False,
                index=False,
                na_rep="",
                float_format=lambda x: format(x, ".15g") if math.isfinite(x) else "",
            )
    return 0


def run_legacy_apply(
    *,
    preview_header: List[str],
    preview_delimiter: Optional[str],
    metadata: Dict[str, object],
    vmap_frame: pd.DataFrame,
    rows_by_provenance: Dict[Tuple[str, int], List[str]],
    output_path: Path,
    retain_snp_id: bool,
) -> int:
    col_name_chr = metadata_column_name(metadata, "col_CHR")
    col_name_pos = metadata_column_name(metadata, "col_POS")
    col_name_id = metadata_column_name(metadata, "col_SNP")
    col_name_a1 = metadata_column_name(metadata, "col_EffectAllele")
    col_name_a2 = metadata_column_name(metadata, "col_OtherAllele")
    idx_chr = resolve_column(preview_header, col_name_chr, "col_CHR", required=True)
    effect_columns = resolve_effect_columns(preview_header, metadata)
    idx_pos = output_column_index(preview_header, col_name_pos)
    idx_id = output_column_index(preview_header, col_name_id)
    idx_a1 = output_column_index(preview_header, col_name_a1)
    idx_a2 = output_column_index(preview_header, col_name_a2)
    output_header = list(preview_header)
    output_idx_pos = idx_pos
    output_idx_id = idx_id
    output_idx_a1 = idx_a1
    output_idx_a2 = idx_a2
    missing_variant_columns: Dict[str, str] = {}
    for field_name, column_name, existing_idx in (
        ("POS", col_name_pos, idx_pos),
        ("SNP", col_name_id, idx_id),
        ("EffectAllele", col_name_a1, idx_a1),
        ("OtherAllele", col_name_a2, idx_a2),
    ):
        if column_name is None or existing_idx is not None:
            continue
        if column_name in missing_variant_columns:
            raise ValueError(
                f"apply_vmap_to_sumstats.py cannot infer a joined output representation for missing "
                f"variant columns mapped to {column_name!r}"
            )
        missing_variant_columns[column_name] = field_name
        output_header.append(column_name)
        new_idx = len(output_header) - 1
        if field_name == "POS":
            output_idx_pos = new_idx
        elif field_name == "SNP":
            output_idx_id = new_idx
        elif field_name == "EffectAllele":
            output_idx_a1 = new_idx
        elif field_name == "OtherAllele":
            output_idx_a2 = new_idx
    signed_cols = effect_columns.signed
    invert_cols = effect_columns.invert
    freq_cols = effect_columns.frequency
    warning_keys: Set[Tuple[str, str]] = set()
    with open_text(output_path, "wt") as handle:
        handle.write(join_line(output_header, preview_delimiter) + "\n")
        for vrow in vmap_frame.itertuples(index=False):
            if vrow.source_index == -1:
                cols = build_missing_payload_row(
                    len(output_header),
                    output_header,
                    vrow,
                    idx_chr,
                    idx_pos,
                    idx_id,
                    idx_a1,
                    idx_a2,
                    output_idx_pos,
                    output_idx_id,
                    output_idx_a1,
                    output_idx_a2,
                    retain_snp_id,
                )
                handle.write(join_line(cols, preview_delimiter) + "\n")
                continue
            cols = list(rows_by_provenance[(vrow.source_shard, vrow.source_index)])
            if len(cols) < len(output_header):
                cols.extend([""] * (len(output_header) - len(cols)))
            out_id = output_variant_id(vrow, retain_snp_id=retain_snp_id)
            variant_rewrites: Dict[int, Dict[str, str]] = {}
            variant_rewrites.setdefault(idx_chr, {})["CHR"] = vrow.chrom
            if idx_pos is not None:
                variant_rewrites.setdefault(idx_pos, {})["POS"] = vrow.pos
            if idx_id is not None:
                cols[idx_id] = out_id
            if idx_a1 is not None:
                variant_rewrites.setdefault(idx_a1, {})["EffectAllele"] = vrow.a1
            if idx_a2 is not None:
                variant_rewrites.setdefault(idx_a2, {})["OtherAllele"] = vrow.a2
            for idx, replacements in variant_rewrites.items():
                cols[idx] = rewrite_or_synthesize_variant_value(cols[idx], replacements)
            if idx_pos is None and output_idx_pos is not None:
                cols[output_idx_pos] = vrow.pos
            if idx_id is None and output_idx_id is not None:
                cols[output_idx_id] = out_id
            if idx_a1 is None and output_idx_a1 is not None:
                cols[output_idx_a1] = vrow.a1
            if idx_a2 is None and output_idx_a2 is not None:
                cols[output_idx_a2] = vrow.a2
            if vrow.allele_op in {"swap", "flip_swap"}:
                for idx in signed_cols:
                    if idx is not None:
                        cols[idx] = maybe_negate(
                            cols[idx],
                            column_name=preview_header[idx],
                            warning_keys=warning_keys,
                        )
                or_l95_idx, or_u95_idx = invert_cols[1], invert_cols[2]
                if or_l95_idx is not None and or_u95_idx is not None:
                    cols[or_l95_idx], cols[or_u95_idx] = maybe_invert_interval(
                        cols[or_l95_idx],
                        cols[or_u95_idx],
                        lower_column=preview_header[or_l95_idx],
                        upper_column=preview_header[or_u95_idx],
                        warning_keys=warning_keys,
                    )
                else:
                    for idx in invert_cols[1:]:
                        if idx is not None:
                            cols[idx] = maybe_invert(
                                cols[idx],
                                column_name=preview_header[idx],
                                warning_keys=warning_keys,
                            )
                or_idx = invert_cols[0]
                if or_idx is not None:
                    cols[or_idx] = maybe_invert(
                        cols[or_idx],
                        column_name=preview_header[or_idx],
                        warning_keys=warning_keys,
                    )
                for idx in freq_cols:
                    if idx is not None:
                        cols[idx] = maybe_complement(
                            cols[idx],
                            column_name=preview_header[idx],
                            warning_keys=warning_keys,
                        )
            handle.write(join_line(cols, preview_delimiter) + "\n")
    return 0


def main() -> int:
    args = parse_args()
    meta_path = Path(args.sumstats_metadata)
    vmap_path = Path(args.vmap)
    output_path = Path(args.output)
    logger.info("apply_vmap_to_sumstats.py: applying %s to %s -> %s", vmap_path, args.input or "<metadata path>", output_path)
    if "@" in args.output:
        raise ValueError("apply_vmap_to_sumstats.py does not accept '@' paths")
    if not meta_path.exists():
        raise ValueError(f"metadata file not found: {meta_path}")
    if not vmap_path.exists():
        raise ValueError(f"vmap file not found: {vmap_path}")

    metadata: Dict[str, object] = load_metadata(meta_path)
    input_path = resolve_sumstats_input_path(
        args.input,
        metadata_path=meta_path,
        metadata=metadata,
        consumer_label="apply_vmap_to_sumstats.py",
    )
    if "@" in str(input_path):
        raise ValueError("apply_vmap_to_sumstats.py does not accept '@' paths")
    if not input_path.exists():
        raise ValueError(f"sumstats file not found: {input_path}")
    vmap_meta = load_variant_metadata(vmap_path)
    validate_vmap_metadata(vmap_meta)
    all_vmap_table = read_vmap_table(vmap_path)
    if len(all_vmap_table) == 0:
        raise ValueError("empty vmap")
    require_table_matches_contig_naming(
        all_vmap_table,
        require_contig_naming(dict(vmap_meta["target"]), label="variant map target"),
        label="variant map target",
    )
    vmap_frame = filtered_vmap_frame(all_vmap_table.to_frame(copy=False), only_mapped_target=args.only_mapped_target)
    _preview_header_line, preview_header, preview_delimiter, _preview_header_line_number = read_sumstats_header(input_path)
    required_variant_columns: Dict[str, int]

    if args.clean:
        clean_resolved = resolve_clean_metadata_columns(preview_header, metadata, include_variant_columns=True)
        if not clean_resolved:
            raise ValueError("clean summary-stat application requires metadata-defined columns")
        required_variant_columns = {
            "CHR": clean_resolved["col_CHR"],
            "POS": clean_resolved["col_POS"],
            "EffectAllele": clean_resolved["col_EffectAllele"],
            "OtherAllele": clean_resolved["col_OtherAllele"],
        }
    else:
        col_name_chr = metadata_column_name(metadata, "col_CHR")
        idx_chr = resolve_column(preview_header, col_name_chr, "col_CHR", required=True)
        idx_a1 = output_column_index(preview_header, metadata_column_name(metadata, "col_EffectAllele"))
        idx_a2 = output_column_index(preview_header, metadata_column_name(metadata, "col_OtherAllele"))
        idx_pos_required = output_column_index(preview_header, metadata_column_name(metadata, "col_POS"))
        required_variant_columns = {
            "CHR": idx_chr,
        }
        if idx_a1 is not None:
            required_variant_columns["EffectAllele"] = idx_a1
        if idx_a2 is not None:
            required_variant_columns["OtherAllele"] = idx_a2
        if idx_pos_required is not None:
            required_variant_columns["POS"] = idx_pos_required
    _effect_columns = resolve_effect_columns(preview_header, metadata)

    if args.clean:
        source_indices = require_single_file_sumstats_provenance(vmap_frame)
        sumstats_table = read_sumstats_table(input_path, usecols=clean_sumstats_usecols(preview_header, metadata))
        validate_required_variant_columns(sumstats_table.frame, preview_header, required_variant_columns, path=input_path)
        sumstats_table.frame.drop(columns=clean_variant_column_names(preview_header, metadata), inplace=True)
        validate_source_indices_in_range(source_indices, len(sumstats_table.frame))
        rc = run_clean_apply(
            args,
            metadata=metadata,
            preview_header=preview_header,
            vmap_frame=vmap_frame,
            sumstats_table=sumstats_table,
            output_path=output_path,
        )
    else:
        needed_source_indices = collect_single_file_source_indices(vmap_frame)
        sumstats_table = read_sumstats_table(input_path)
        if needed_source_indices:
            _header_line, _header, _delimiter, rows_by_provenance = load_sumstats_rows_single_file(
                sumstats_table,
                required_variant_columns,
                needed_source_indices,
            )
        else:
            rows_by_provenance = {}
        missing_keys = sorted(
            source_index
            for source_index in needed_source_indices
            if (".", source_index) not in rows_by_provenance
        )
        if missing_keys:
            source_index = missing_keys[0]
            raise ValueError(
                f"vmap source provenance out of range for summary statistics input: "
                f"first missing source_shard='.', source_index={source_index}"
            )
        rc = run_legacy_apply(
            preview_header=preview_header,
            preview_delimiter=preview_delimiter,
            metadata=metadata,
            vmap_frame=vmap_frame,
            rows_by_provenance=rows_by_provenance,
            output_path=output_path,
            retain_snp_id=args.retain_snp_id,
        )
    logger.info("apply_vmap_to_sumstats.py: wrote %s with %s retained target rows", output_path, len(vmap_frame))
    return rc


if __name__ == "__main__":
    raise SystemExit(run_cli(main))
