#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from ._cli_utils import run_cli
from .apply_vmap_utils import build_needed_source_indices, filtered_vmap_rows
from .sumstats_clean import harmonize_clean_sumstats, resolve_clean_metadata_columns
from .sumstats_utils import (
    find_metadata_value,
    join_line,
    load_metadata,
    open_text,
    open_sumstats_data,
    resolve_column,
    resolve_effect_columns,
    rewrite_variant_fields,
    split_line,
)
from .vtable_utils import (
    load_metadata as load_variant_metadata,
    read_vmap,
    require_contig_naming,
    require_rows_match_contig_naming,
    validate_vmap_metadata,
)


NA_NUMERIC = "n/a"
JOINED_VARIANT_FIELD_ORDER = ["CHR", "POS", "EffectAllele", "OtherAllele"]


def warn_once(warning_keys: Set[Tuple[str, str]], key: Tuple[str, str], message: str) -> None:
    if key in warning_keys:
        return
    warning_keys.add(key)
    print(f"Warning: {message}", file=sys.stderr)


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
    parser.add_argument("--input", required=True, help="Input summary statistics file")
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
    return parser.parse_args()


def max_required_column(*indices: int | None) -> int:
    return max(idx for idx in indices if idx is not None)


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
    by_name = {name: idx for idx, name in enumerate(header)}
    by_lower = {name.lower(): idx for idx, name in enumerate(header)}
    if column_name in by_name:
        return by_name[column_name]
    return by_lower.get(column_name.lower())


def load_sumstats_rows_single_file(
    input_path: Path,
    min_columns: int,
    needed_by_shard: Dict[str, Set[int]],
) -> Tuple[str, List[str], Optional[str], Dict[Tuple[str, int], List[str]]]:
    if set(needed_by_shard) - {"."}:
        unsupported = sorted(set(needed_by_shard) - {"."})
        raise ValueError(
            "apply_vmap_to_sumstats.py supports single-file payload lookup only for source_shard='.'; "
            f"found {unsupported!r}"
        )
    with open_sumstats_data(input_path) as (handle, header_line, header, delimiter):
        rows_by_provenance: Dict[Tuple[str, int], List[str]] = {}
        row_index = 0
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = split_line(line, delimiter)
            if len(cols) < min_columns:
                raise ValueError(f"sumstats row has fewer columns than expected: {input_path}")
            if row_index in needed_by_shard.get(".", set()):
                rows_by_provenance[(".", row_index)] = cols
            row_index += 1
        return header_line, header, delimiter, rows_by_provenance


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
) -> List[str]:
    cols = [NA_NUMERIC] * output_width
    variant_rewrites: Dict[int, Dict[str, str]] = {}
    variant_rewrites.setdefault(idx_chr, {})["CHR"] = vrow.chrom
    if idx_pos is not None:
        variant_rewrites.setdefault(idx_pos, {})["POS"] = vrow.pos
    if idx_id is not None:
        cols[idx_id] = vrow.id
    if idx_a1 is not None:
        variant_rewrites.setdefault(idx_a1, {})["EffectAllele"] = vrow.a1
    if idx_a2 is not None:
        variant_rewrites.setdefault(idx_a2, {})["OtherAllele"] = vrow.a2
    for idx, replacements in variant_rewrites.items():
        cols[idx] = rewrite_or_synthesize_variant_value(cols[idx], replacements)
    if idx_pos is None and output_idx_pos is not None:
        cols[output_idx_pos] = vrow.pos
    if idx_id is None and output_idx_id is not None:
        cols[output_idx_id] = vrow.id
    if idx_a1 is None and output_idx_a1 is not None:
        cols[output_idx_a1] = vrow.a1
    if idx_a2 is None and output_idx_a2 is not None:
        cols[output_idx_a2] = vrow.a2
    return cols


def parse_clean_numeric(raw: object) -> float:
    if raw is None:
        raise ValueError("missing value")
    if isinstance(raw, str) and not raw.strip():
        raise ValueError("missing value")
    value = float(raw)
    if not math.isfinite(value):
        raise ValueError(f"non-finite numeric value: {raw!r}")
    return value


def maybe_negate_clean(raw: object, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> float:
    try:
        return -parse_clean_numeric(raw)
    except Exception:
        warn_once(
            warning_keys,
            ("negate_clean", column_name),
            f"could not negate non-numeric value in column {column_name!r}; writing missing value",
        )
        return float("nan")


def maybe_invert_clean(raw: object, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> float:
    try:
        value = parse_clean_numeric(raw)
        if value == 0:
            raise ValueError("cannot invert zero")
        return 1.0 / value
    except Exception:
        warn_once(
            warning_keys,
            ("invert_clean", column_name),
            f"could not invert zero or non-numeric value in column {column_name!r}; writing missing value",
        )
        return float("nan")


def maybe_complement_clean(raw: object, *, column_name: str, warning_keys: Set[Tuple[str, str]]) -> float:
    try:
        value = parse_clean_numeric(raw)
        return 1.0 - value
    except Exception:
        warn_once(
            warning_keys,
            ("complement_clean", column_name),
            f"could not complement non-numeric value in column {column_name!r}; writing missing value",
        )
        return float("nan")


def maybe_invert_interval_clean(
    lower_raw: object,
    upper_raw: object,
    *,
    lower_column: str,
    upper_column: str,
    warning_keys: Set[Tuple[str, str]],
) -> Tuple[float, float]:
    try:
        lower = parse_clean_numeric(lower_raw)
        upper = parse_clean_numeric(upper_raw)
        if lower == 0 or upper == 0:
            raise ValueError("cannot invert zero")
        return 1.0 / upper, 1.0 / lower
    except Exception:
        warn_once(
            warning_keys,
            ("invert_interval_clean", f"{lower_column}|{upper_column}"),
            "could not invert zero or non-numeric odds-ratio interval in columns "
            f"{lower_column!r} and {upper_column!r}; writing missing value",
        )
        return float("nan"), float("nan")


def collect_clean_rows(
    preview_header: Sequence[str],
    metadata: Dict[str, object],
    vmap_rows,
    rows_by_provenance: Dict[Tuple[str, int], List[str]],
) -> Tuple[List[str], List[List[object]]]:
    resolved = resolve_clean_metadata_columns(preview_header, metadata, include_variant_columns=True)
    variant_indices = {
        resolved[key]
        for key in ("col_CHR", "col_POS", "col_SNP", "col_EffectAllele", "col_OtherAllele")
        if key in resolved
    }
    payload_indices = [idx for idx in range(len(preview_header)) if idx not in variant_indices]
    payload_header = [preview_header[idx] for idx in payload_indices]
    payload_rows: List[List[object]] = []
    for vrow in vmap_rows:
        if vrow.source_index == -1:
            payload_rows.append([None] * len(payload_indices))
            continue
        cols = rows_by_provenance[(vrow.source_shard, vrow.source_index)]
        payload_rows.append([cols[idx] for idx in payload_indices])
    return payload_header, payload_rows


def clear_clean_unmatched_rows(payload_sumstats, vmap_rows) -> None:
    unmatched_indices = [idx for idx, vrow in enumerate(vmap_rows) if vrow.source_index == -1]
    if unmatched_indices:
        payload_sumstats.iloc[unmatched_indices, :] = float("nan")


def format_clean_value(value: object) -> str:
    if value is None:
        return ""
    try:
        if value != value:
            return ""
    except Exception:
        pass
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return format(value, ".15g")
    return str(value)


def run_clean_apply(
    args: argparse.Namespace,
    *,
    metadata: Dict[str, object],
    preview_header: List[str],
    vmap_rows,
    rows_by_provenance: Dict[Tuple[str, int], List[str]],
    output_path: Path,
) -> int:
    payload_header, payload_rows = collect_clean_rows(preview_header, metadata, vmap_rows, rows_by_provenance)
    payload_sumstats, _clean_metadata = harmonize_clean_sumstats(
        payload_header,
        payload_rows,
        metadata,
        fill_mode=args.fill_mode,
        use_af_inference=args.use_af_inference,
        warn=lambda message: print(f"Warning: {message}", file=sys.stderr),
    )
    clear_clean_unmatched_rows(payload_sumstats, vmap_rows)
    warning_keys: Set[Tuple[str, str]] = set()
    for idx, vrow in enumerate(vmap_rows):
        if vrow.allele_op not in {"swap", "flip_swap"}:
            continue
        for column in ("BETA", "Z"):
            if column in payload_sumstats.columns:
                payload_sumstats.at[idx, column] = maybe_negate_clean(
                    payload_sumstats.at[idx, column],
                    column_name=column,
                    warning_keys=warning_keys,
                )
        if "ORL95" in payload_sumstats.columns and "ORU95" in payload_sumstats.columns:
            payload_sumstats.at[idx, "ORL95"], payload_sumstats.at[idx, "ORU95"] = maybe_invert_interval_clean(
                payload_sumstats.at[idx, "ORL95"],
                payload_sumstats.at[idx, "ORU95"],
                lower_column="ORL95",
                upper_column="ORU95",
                warning_keys=warning_keys,
            )
        else:
            for column in ("ORL95", "ORU95"):
                if column in payload_sumstats.columns:
                    payload_sumstats.at[idx, column] = maybe_invert_clean(
                        payload_sumstats.at[idx, column],
                        column_name=column,
                        warning_keys=warning_keys,
                    )
        if "OR" in payload_sumstats.columns:
            payload_sumstats.at[idx, "OR"] = maybe_invert_clean(
                payload_sumstats.at[idx, "OR"],
                column_name="OR",
                warning_keys=warning_keys,
            )
        for column in ("EAF", "CaseEAF", "ControlEAF"):
            if column in payload_sumstats.columns:
                payload_sumstats.at[idx, column] = maybe_complement_clean(
                    payload_sumstats.at[idx, column],
                    column_name=column,
                    warning_keys=warning_keys,
                )

    retained_indices = list(range(len(vmap_rows)))
    if args.only_mapped_target:
        if "P" not in payload_sumstats.columns:
            retained_indices = []
        else:
            retained_indices = [idx for idx in retained_indices if format_clean_value(payload_sumstats.at[idx, "P"]) != ""]
        if not retained_indices:
            raise ValueError("no retained target rows remain after applying --only-mapped-target")

    with open_text(output_path, "wt") as handle:
        output_header = ["CHR", "POS", "SNP", "EffectAllele", "OtherAllele", *list(payload_sumstats.columns)]
        handle.write("\t".join(output_header) + "\n")
        for idx in retained_indices:
            vrow = vmap_rows[idx]
            row_values = [
                vrow.chrom,
                vrow.pos,
                vrow.id,
                vrow.a1,
                vrow.a2,
                *[format_clean_value(payload_sumstats.at[idx, column]) for column in payload_sumstats.columns],
            ]
            handle.write("\t".join(row_values) + "\n")
    return 0


def run_legacy_apply(
    *,
    preview_header: List[str],
    preview_delimiter: Optional[str],
    metadata: Dict[str, object],
    vmap_rows,
    rows_by_provenance: Dict[Tuple[str, int], List[str]],
    output_path: Path,
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
        for vrow in vmap_rows:
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
                )
                handle.write(join_line(cols, preview_delimiter) + "\n")
                continue
            cols = list(rows_by_provenance[(vrow.source_shard, vrow.source_index)])
            if len(cols) < len(output_header):
                cols.extend([""] * (len(output_header) - len(cols)))
            variant_rewrites: Dict[int, Dict[str, str]] = {}
            variant_rewrites.setdefault(idx_chr, {})["CHR"] = vrow.chrom
            if idx_pos is not None:
                variant_rewrites.setdefault(idx_pos, {})["POS"] = vrow.pos
            if idx_id is not None:
                cols[idx_id] = vrow.id
            if idx_a1 is not None:
                variant_rewrites.setdefault(idx_a1, {})["EffectAllele"] = vrow.a1
            if idx_a2 is not None:
                variant_rewrites.setdefault(idx_a2, {})["OtherAllele"] = vrow.a2
            for idx, replacements in variant_rewrites.items():
                cols[idx] = rewrite_or_synthesize_variant_value(cols[idx], replacements)
            if idx_pos is None and output_idx_pos is not None:
                cols[output_idx_pos] = vrow.pos
            if idx_id is None and output_idx_id is not None:
                cols[output_idx_id] = vrow.id
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
    input_path = Path(args.input)
    meta_path = Path(args.sumstats_metadata)
    vmap_path = Path(args.vmap)
    output_path = Path(args.output)
    if "@" in args.input or "@" in args.output:
        raise ValueError("apply_vmap_to_sumstats.py does not accept '@' paths")
    if not input_path.exists():
        raise ValueError(f"sumstats file not found: {input_path}")
    if not meta_path.exists():
        raise ValueError(f"metadata file not found: {meta_path}")
    if not vmap_path.exists():
        raise ValueError(f"vmap file not found: {vmap_path}")

    metadata: Dict[str, object] = load_metadata(meta_path)
    vmap_meta = load_variant_metadata(vmap_path)
    validate_vmap_metadata(vmap_meta)
    all_vmap_rows = read_vmap(vmap_path)
    if not all_vmap_rows:
        raise ValueError("empty vmap")
    require_rows_match_contig_naming(
        all_vmap_rows,
        require_contig_naming(dict(vmap_meta["target"]), label="variant map target"),
        label="variant map target",
    )
    vmap_rows = filtered_vmap_rows(all_vmap_rows, only_mapped_target=args.only_mapped_target)
    with open_sumstats_data(input_path) as (_preview_handle, _header_line, preview_header, preview_delimiter):
        pass

    if args.clean:
        clean_resolved = resolve_clean_metadata_columns(preview_header, metadata, include_variant_columns=True)
        if not clean_resolved:
            raise ValueError("clean summary-stat application requires metadata-defined columns")
        min_columns = max(clean_resolved.values()) + 1
    else:
        col_name_chr = metadata_column_name(metadata, "col_CHR")
        idx_chr = resolve_column(preview_header, col_name_chr, "col_CHR", required=True)
        effect_columns = resolve_effect_columns(preview_header, metadata)
        min_columns = max_required_column(idx_chr, *effect_columns.signed, *effect_columns.invert, *effect_columns.frequency) + 1

    needed_by_shard = build_needed_source_indices(vmap_rows)
    if needed_by_shard:
        _header_line, _header, _delimiter, rows_by_provenance = load_sumstats_rows_single_file(
            input_path,
            min_columns,
            needed_by_shard,
        )
    else:
        rows_by_provenance = {}
    missing_keys = sorted(
        (source_shard, source_index)
        for source_shard, indices in needed_by_shard.items()
        for source_index in indices
        if (source_shard, source_index) not in rows_by_provenance
    )
    if missing_keys:
        source_shard, source_index = missing_keys[0]
        raise ValueError(
            f"vmap source provenance out of range for summary statistics input: "
            f"first missing source_shard={source_shard!r}, source_index={source_index}"
        )

    if args.clean:
        return run_clean_apply(
            args,
            metadata=metadata,
            preview_header=preview_header,
            vmap_rows=vmap_rows,
            rows_by_provenance=rows_by_provenance,
            output_path=output_path,
        )
    return run_legacy_apply(
        preview_header=preview_header,
        preview_delimiter=preview_delimiter,
        metadata=metadata,
        vmap_rows=vmap_rows,
        rows_by_provenance=rows_by_provenance,
        output_path=output_path,
    )


if __name__ == "__main__":
    raise SystemExit(run_cli(main))
