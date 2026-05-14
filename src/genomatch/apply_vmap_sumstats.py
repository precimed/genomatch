#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from ._cli_utils import run_cli
from .sumstats_clean import (
    PAYLOAD_METADATA_TO_CANONICAL,
    harmonize_clean_sumstats_frame,
)
from .sumstats_metadata import (
    retained_payload_columns as resolve_retained_payload_columns,
)
from .sumstats_utils import (
    MISSING_VALUE_TOKENS,
    SumstatsTable,
    find_metadata_value,
    load_metadata,
    open_text,
    read_sumstats_header,
    read_sumstats_table,
    resolve_sumstats_input_path,
    split_line,
)
from .vtable_utils import (
    load_metadata as load_variant_metadata,
    require_contig_naming,
    require_table_matches_contig_naming,
    read_vmap_table,
    validate_vmap_metadata,
)


VARIANT_OUTPUT_COLUMNS = ("CHR", "POS", "SNP", "EffectAllele", "OtherAllele")
logger = logging.getLogger(__name__)


def warn_once(warning_keys: Set[Tuple[str, str]], key: Tuple[str, str], message: str) -> None:
    if key in warning_keys:
        return
    warning_keys.add(key)
    logger.warning("%s", message)


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
        "--retain-snp-id",
        action="store_true",
        help="Use retained target-side .vmap id values as output SNP IDs instead of generated chrom:pos:a1:a2 IDs",
    )
    return parser.parse_args()


def retained_payload_columns(
    preview_header: Sequence[str],
    metadata: Dict[str, object],
) -> Tuple[List[Tuple[str, int, str]], Dict[str, str]]:
    return resolve_retained_payload_columns(
        preview_header,
        metadata,
        payload_keys=list(PAYLOAD_METADATA_TO_CANONICAL),
        reserved_output_columns=VARIANT_OUTPUT_COLUMNS,
    )


def require_single_file_sumstats_provenance(vmap_frame: pd.DataFrame) -> pd.Series:
    provenance_frame = vmap_frame.loc[:, ["source_shard", "source_index"]]
    invalid_shard_mask = provenance_frame["source_shard"].ne(".")
    if invalid_shard_mask.any():
        unsupported = sorted(provenance_frame.loc[invalid_shard_mask, "source_shard"].astype(str).unique().tolist())
        raise ValueError(
            "apply_vmap_to_sumstats.py supports single-file payload lookup only for source_shard='.'; "
            f"found {unsupported!r}"
        )
    return provenance_frame["source_index"].astype("int64")


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


def validate_required_payload_row_widths(
    input_path: Path,
    *,
    delimiter: Optional[str],
    header_line_number: int,
    retained_columns: Sequence[Tuple[str, int, str]],
) -> None:
    if not retained_columns:
        return
    max_required_idx = max(idx for _key, idx, _output_name in retained_columns)
    with open_text(input_path, "rt") as handle:
        for line_number, line in enumerate(handle):
            if line_number <= header_line_number:
                continue
            if not line.strip() or line.startswith("#"):
                continue
            if len(split_line(line, delimiter)) <= max_required_idx:
                raise ValueError(
                    "sumstats input row has fewer fields than required by metadata-declared payload columns: "
                    f"{input_path}:{line_number + 1}"
                )


def input_missing_value(metadata: Dict[str, object]) -> Optional[str]:
    value = find_metadata_value(metadata, "missing_value")
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError("metadata missing_value must be a string")
    return value


def normalize_missing_values(series: pd.Series, *, missing_value: Optional[str]) -> pd.Series:
    out = series.copy()
    missing_mask = out.isna()
    text = out.astype(str).str.strip()
    missing_mask = missing_mask | text.str.lower().isin(MISSING_VALUE_TOKENS)
    if missing_value is not None:
        missing_mask = missing_mask | out.astype(str).eq(missing_value)
    out.loc[missing_mask] = np.nan
    return out


def collect_payload_frame(
    sumstats_table: SumstatsTable,
    vmap_frame: pd.DataFrame,
    retained_columns: Sequence[Tuple[str, int, str]],
    *,
    metadata: Dict[str, object],
) -> pd.DataFrame:
    payload_frame = pd.DataFrame(index=range(len(vmap_frame)))
    if not retained_columns:
        return payload_frame
    source_indices = vmap_frame["source_index"].to_numpy(dtype=np.int64, copy=False)
    missing_value = input_missing_value(metadata)
    for _key, idx, output_name in retained_columns:
        values = sumstats_table.frame.iloc[source_indices, idx].reset_index(drop=True)
        payload_frame[output_name] = normalize_missing_values(values, missing_value=missing_value)
    return payload_frame


def output_variant_ids(vmap_frame: pd.DataFrame, *, retain_snp_id: bool) -> pd.Series:
    if retain_snp_id:
        return vmap_frame["id"].astype(str).reset_index(drop=True)
    return (
        vmap_frame["chrom"].astype(str)
        + ":"
        + vmap_frame["pos"].astype(str)
        + ":"
        + vmap_frame["a1"].astype(str)
        + ":"
        + vmap_frame["a2"].astype(str)
    ).reset_index(drop=True)


def build_variant_output_frame(vmap_frame: pd.DataFrame, *, retain_snp_id: bool) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "CHR": vmap_frame["chrom"].to_numpy(copy=False),
            "POS": vmap_frame["pos"].to_numpy(copy=False),
            "SNP": output_variant_ids(vmap_frame, retain_snp_id=retain_snp_id),
            "EffectAllele": vmap_frame["a1"].to_numpy(copy=False),
            "OtherAllele": vmap_frame["a2"].to_numpy(copy=False),
        },
        index=range(len(vmap_frame)),
    )


def numeric_transform(
    series: pd.Series,
    *,
    column_name: str,
    warning_keys: Set[Tuple[str, str]],
    warning_key: str,
    invalid_message: str,
    transform,
    invalid_extra=None,
) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce").astype(float)
    invalid_mask = ~np.isfinite(numeric)
    if invalid_extra is not None:
        invalid_mask = invalid_mask | invalid_extra(numeric)
    if invalid_mask.any():
        warn_once(warning_keys, (warning_key, column_name), invalid_message.format(column=column_name))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        out = transform(numeric)
    out.loc[invalid_mask] = np.nan
    return out


def apply_swap_effect_transforms(
    payload_frame: pd.DataFrame,
    vmap_frame: pd.DataFrame,
    column_by_key: Dict[str, str],
) -> pd.DataFrame:
    if payload_frame.empty:
        return payload_frame
    out = payload_frame.copy()
    swap_mask = vmap_frame["allele_op"].isin({"swap", "flip_swap"}).to_numpy()
    if not swap_mask.any():
        return out
    warning_keys: Set[Tuple[str, str]] = set()

    for key in ("col_BETA", "col_Z"):
        column = column_by_key.get(key)
        if column in out.columns:
            out.loc[swap_mask, column] = numeric_transform(
                out.loc[swap_mask, column],
                column_name=column,
                warning_keys=warning_keys,
                warning_key="negate",
                invalid_message="could not negate non-numeric value in column {column!r}; writing missing value",
                transform=lambda numeric: -numeric,
            )

    lower_column = column_by_key.get("col_ORL95")
    upper_column = column_by_key.get("col_ORU95")
    if lower_column in out.columns and upper_column in out.columns:
        lower = pd.to_numeric(out.loc[swap_mask, lower_column], errors="coerce").astype(float)
        upper = pd.to_numeric(out.loc[swap_mask, upper_column], errors="coerce").astype(float)
        invalid_mask = ~np.isfinite(lower) | ~np.isfinite(upper) | (lower == 0) | (upper == 0)
        if invalid_mask.any():
            warn_once(
                warning_keys,
                ("invert_interval", f"{lower_column}|{upper_column}"),
                "could not invert zero or non-numeric odds-ratio interval in columns "
                f"{lower_column!r} and {upper_column!r}; writing missing value",
            )
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            lower_out = 1.0 / upper
            upper_out = 1.0 / lower
        lower_out.loc[invalid_mask] = np.nan
        upper_out.loc[invalid_mask] = np.nan
        out.loc[swap_mask, lower_column] = lower_out
        out.loc[swap_mask, upper_column] = upper_out
    else:
        for column in (lower_column, upper_column):
            if column in out.columns:
                out.loc[swap_mask, column] = numeric_transform(
                    out.loc[swap_mask, column],
                    column_name=column,
                    warning_keys=warning_keys,
                    warning_key="invert",
                    invalid_message="could not invert zero or non-numeric value in column {column!r}; writing missing value",
                    transform=lambda numeric: 1.0 / numeric,
                    invalid_extra=lambda numeric: numeric == 0,
                )

    or_column = column_by_key.get("col_OR")
    if or_column in out.columns:
        out.loc[swap_mask, or_column] = numeric_transform(
            out.loc[swap_mask, or_column],
            column_name=or_column,
            warning_keys=warning_keys,
            warning_key="invert",
            invalid_message="could not invert zero or non-numeric value in column {column!r}; writing missing value",
            transform=lambda numeric: 1.0 / numeric,
            invalid_extra=lambda numeric: numeric == 0,
        )

    for key in ("col_EAF", "col_CaseEAF", "col_ControlEAF"):
        column = column_by_key.get(key)
        if column in out.columns:
            out.loc[swap_mask, column] = numeric_transform(
                out.loc[swap_mask, column],
                column_name=column,
                warning_keys=warning_keys,
                warning_key="complement",
                invalid_message="could not complement non-numeric value in column {column!r}; writing missing value",
                transform=lambda numeric: 1.0 - numeric,
            )
    return out


def write_sumstats_output(output_path: Path, output_frame: pd.DataFrame) -> None:
    with open_text(output_path, "wt") as handle:
        output_frame.to_csv(
            handle,
            sep="\t",
            header=True,
            index=False,
            na_rep="",
            float_format=lambda x: format(x, ".15g") if math.isfinite(x) else "",
        )


def output_sidecar_path(output_path: Path) -> Path:
    return output_path.with_name(output_path.name + ".meta.yaml")


def output_stats_model(metadata: Dict[str, object]) -> str:
    value = find_metadata_value(metadata, "stats_Model")
    if isinstance(value, str) and value in {
        "linear",
        "logistic",
        "ordinal",
        "linear mixed-model",
        "logistic mixed-model",
        "cox",
        "other",
    }:
        return value
    return "other"


def write_output_sidecar(
    output_path: Path,
    *,
    vmap_meta: Dict[str, object],
    source_metadata: Dict[str, object],
    payload_mappings: Dict[str, str],
    clean: bool,
    fill_mode: str,
    use_af_inference: bool,
) -> None:
    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise ValueError("PyYAML is required to write summary-stat metadata YAML") from exc

    target_meta = vmap_meta["target"]
    sidecar: Dict[str, object] = {
        "cleansumstats_metafile_kind": "minimal",
        "path_sumStats": output_path.name,
        "delimiter": "\t",
        "missing_value": "",
        "genome_build": target_meta["genome_build"],
        "contig_naming": target_meta["contig_naming"],
        "stats_Model": output_stats_model(source_metadata),
        "col_CHR": "CHR",
        "col_POS": "POS",
        "col_SNP": "SNP",
        "col_EffectAllele": "EffectAllele",
        "col_OtherAllele": "OtherAllele",
        "clean": clean,
    }
    for key in PAYLOAD_METADATA_TO_CANONICAL:
        if key in payload_mappings:
            sidecar[key] = payload_mappings[key]
    if clean:
        sidecar["fill_mode"] = fill_mode
        sidecar["use_af_inference"] = use_af_inference

    output_sidecar_path(output_path).write_text(yaml.safe_dump(sidecar, sort_keys=False), encoding="utf-8")


def run_sumstats_apply(
    args: argparse.Namespace,
    *,
    metadata: Dict[str, object],
    retained_columns: Sequence[Tuple[str, int, str]],
    projection_payload_mappings: Dict[str, str],
    vmap_frame: pd.DataFrame,
    vmap_meta: Dict[str, object],
    sumstats_table: SumstatsTable,
    output_path: Path,
) -> int:
    if not retained_columns:
        raise ValueError("apply_vmap_to_sumstats.py requires at least one metadata-declared payload/stat column")
    source_indices = require_single_file_sumstats_provenance(vmap_frame)
    validate_source_indices_in_range(source_indices, len(sumstats_table.frame))

    payload_frame = collect_payload_frame(
        sumstats_table,
        vmap_frame,
        retained_columns,
        metadata=metadata,
    )
    if args.clean:
        payload_sumstats, clean_metadata = harmonize_clean_sumstats_frame(
            list(payload_frame.columns),
            payload_frame,
            metadata,
            fill_mode=args.fill_mode,
            use_af_inference=args.use_af_inference,
            warn=lambda message: logger.warning("%s", message),
        )
        clean_payload_mappings = {
            key: canonical
            for key, canonical in PAYLOAD_METADATA_TO_CANONICAL.items()
            if canonical in payload_sumstats.columns and find_metadata_value(clean_metadata, key) is not None
        }
        payload_sumstats = apply_swap_effect_transforms(payload_sumstats, vmap_frame, clean_payload_mappings)
        payload_mappings = clean_payload_mappings
    else:
        payload_sumstats = apply_swap_effect_transforms(payload_frame, vmap_frame, projection_payload_mappings)
        payload_mappings = projection_payload_mappings

    output_frame = pd.concat(
        [build_variant_output_frame(vmap_frame, retain_snp_id=args.retain_snp_id), payload_sumstats.reset_index(drop=True)],
        axis=1,
    )
    write_sumstats_output(output_path, output_frame)
    write_output_sidecar(
        output_path,
        vmap_meta=vmap_meta,
        source_metadata=metadata,
        payload_mappings=payload_mappings,
        clean=args.clean,
        fill_mode=args.fill_mode,
        use_af_inference=args.use_af_inference,
    )
    return 0


def main() -> int:
    args = parse_args()
    meta_path = Path(args.sumstats_metadata)
    vmap_path = Path(args.vmap)
    output_path = Path(args.output)
    logger.info(
        "apply_vmap_to_sumstats.py: applying %s to %s -> %s",
        vmap_path,
        args.input or "<metadata path>",
        output_path,
    )
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
    vmap_frame = all_vmap_table.to_frame(copy=False).reset_index(drop=True)
    _preview_header_line, preview_header, preview_delimiter, preview_header_line_number = read_sumstats_header(input_path)
    # Resolve metadata-declared payload columns before reading rows so duplicate
    # mappings and output-name collisions fail at the metadata/header boundary.
    retained_columns, projection_payload_mappings = retained_payload_columns(preview_header, metadata)
    if not retained_columns:
        raise ValueError("apply_vmap_to_sumstats.py requires at least one metadata-declared payload/stat column")
    validate_required_payload_row_widths(
        input_path,
        delimiter=preview_delimiter,
        header_line_number=preview_header_line_number,
        retained_columns=retained_columns,
    )
    sumstats_table = read_sumstats_table(input_path)
    rc = run_sumstats_apply(
        args,
        metadata=metadata,
        retained_columns=retained_columns,
        projection_payload_mappings=projection_payload_mappings,
        vmap_frame=vmap_frame,
        vmap_meta=vmap_meta,
        sumstats_table=sumstats_table,
        output_path=output_path,
    )
    logger.info("apply_vmap_to_sumstats.py: wrote %s with %s retained target rows", output_path, len(vmap_frame))
    return rc


if __name__ == "__main__":
    raise SystemExit(run_cli(main))
