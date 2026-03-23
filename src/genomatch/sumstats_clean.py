from __future__ import annotations

import math
import re
from copy import deepcopy
from statistics import NormalDist
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

from .sumstats_utils import find_metadata_value


VARIANT_METADATA_KEYS = (
    "col_CHR",
    "col_POS",
    "col_SNP",
    "col_EffectAllele",
    "col_OtherAllele",
)

PAYLOAD_METADATA_TO_CANONICAL = {
    "col_BETA": "BETA",
    "col_SE": "SE",
    "col_Z": "Z",
    "col_P": "P",
    "col_OR": "OR",
    "col_ORL95": "ORL95",
    "col_ORU95": "ORU95",
    "col_N": "N",
    "col_CaseN": "CaseN",
    "col_ControlN": "ControlN",
    "col_StudyN": "StudyN",
    "col_INFO": "INFO",
    "col_Direction": "Direction",
    "col_EAF": "EAF",
    "col_CaseEAF": "CaseEAF",
    "col_ControlEAF": "ControlEAF",
    "col_OAF": "OAF",
    "col_CaseOAF": "CaseOAF",
    "col_ControlOAF": "ControlOAF",
}

ALL_METADATA_TO_CANONICAL = {
    **{
        "col_CHR": "CHR",
        "col_POS": "POS",
        "col_SNP": "SNP",
        "col_EffectAllele": "EffectAllele",
        "col_OtherAllele": "OtherAllele",
    },
    **PAYLOAD_METADATA_TO_CANONICAL,
}

PAYLOAD_OUTPUT_ORDER = [
    "P",
    "Z",
    "N",
    "BETA",
    "OR",
    "SE",
    "ORL95",
    "ORU95",
    "CaseN",
    "ControlN",
    "StudyN",
    "INFO",
    "Direction",
    "EAF",
    "CaseEAF",
    "ControlEAF",
]

FREQUENCY_COLUMNS = ("EAF", "CaseEAF", "ControlEAF", "OAF", "CaseOAF", "ControlOAF")
NORMAL_DIST = NormalDist()
NON_ALNUM = re.compile(r"[^A-Za-z0-9]")


def normalize_header_token(value: str) -> str:
    return NON_ALNUM.sub("", value).lower().strip()


def is_missing_value(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    try:
        return bool(pd.isna(value))
    except Exception:
        return False


def normalize_cell_value(value: object) -> object:
    return None if is_missing_value(value) else value


def resolve_clean_metadata_columns(
    header: Sequence[str],
    metadata: Dict[str, object],
    *,
    include_variant_columns: bool,
) -> Dict[str, int]:
    keys = list(PAYLOAD_METADATA_TO_CANONICAL)
    if include_variant_columns:
        keys = [*VARIANT_METADATA_KEYS, *keys]
    normalized_header = [normalize_header_token(str(column)) for column in header]
    resolved: Dict[str, int] = {}
    seen_canonical: Dict[str, str] = {}
    for key in keys:
        raw_value = find_metadata_value(metadata, key)
        if raw_value is None:
            continue
        if not isinstance(raw_value, str):
            raise ValueError(f"column mapping for {key} must be a string column name")
        normalized_value = normalize_header_token(raw_value)
        matches = [idx for idx, token in enumerate(normalized_header) if token == normalized_value]
        if len(matches) != 1:
            raise ValueError(f"column mapping for {key} could not be identified uniquely after header normalization")
        canonical = ALL_METADATA_TO_CANONICAL[key]
        existing = seen_canonical.get(canonical)
        if existing is not None and existing != key:
            raise ValueError(
                f"metadata-declared input columns collapse to the same canonical output name {canonical!r}"
            )
        seen_canonical[canonical] = key
        resolved[key] = matches[0]
    return resolved


def metadata_numeric_value(metadata: Dict[str, object], key: str) -> float | None:
    value = find_metadata_value(metadata, key)
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except Exception:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def warn_message(warn: Callable[[str], None] | None, message: str) -> None:
    if warn is not None:
        warn(message)


def normalize_model(metadata: Dict[str, object], warn: Callable[[str], None] | None) -> str | None:
    raw_model = find_metadata_value(metadata, "stats_Model")
    if raw_model is None:
        return None
    if not isinstance(raw_model, str):
        raise ValueError("stats_Model must be a string")
    model = raw_model.strip().lower()
    if model == "linear mixed-model":
        model = "linear"
    elif model == "logistic mixed-model":
        model = "logistic"
    if model in {"ordinal", "cox", "other"}:
        raise ValueError("not implemented")
    if model not in {"linear", "logistic"}:
        return model
    metadata["stats_Model"] = model
    if model == "logistic" and find_metadata_value(metadata, "stats_TotalN") is not None:
        warn_message(warn, "ignoring stats_TotalN for logistic model")
    return model


def drop_all_missing_columns(sumstats: pd.DataFrame, metadata: Dict[str, object]) -> pd.DataFrame:
    for key, canonical in PAYLOAD_METADATA_TO_CANONICAL.items():
        if canonical not in sumstats.columns:
            metadata.pop(key, None)
            continue
        if sumstats[canonical].map(is_missing_value).all():
            sumstats = sumstats.drop(columns=[canonical])
            metadata.pop(key, None)
    return sumstats


def _set_missing_on_mask(series: pd.Series, mask: pd.Series) -> pd.Series:
    series = series.copy()
    series.loc[mask] = np.nan
    return series


def validate_numeric_ranges(sumstats: pd.DataFrame) -> pd.DataFrame:
    for column in list(sumstats.columns):
        if column == "Direction":
            continue
        numeric = pd.to_numeric(sumstats[column], errors="coerce")
        finite_mask = np.isfinite(numeric)
        numeric = _set_missing_on_mask(numeric, ~finite_mask)
        if column == "SE":
            numeric = _set_missing_on_mask(numeric, numeric <= 0)
        elif column == "P":
            numeric = _set_missing_on_mask(numeric, (numeric < 0) | (numeric > 1))
        elif column in {"OR", "ORL95", "ORU95"}:
            numeric = _set_missing_on_mask(numeric, numeric <= 0)
        elif column in {"N", "CaseN", "ControlN"}:
            numeric = _set_missing_on_mask(numeric, numeric < 0)
        elif column == "INFO":
            numeric = _set_missing_on_mask(numeric, numeric < 0)
        elif column in FREQUENCY_COLUMNS:
            numeric = _set_missing_on_mask(numeric, (numeric < 0) | (numeric > 1))
        sumstats[column] = numeric
    return sumstats


def safe_formula(func: Callable[..., float], *args: object) -> float:
    clean_args: List[float] = []
    for arg in args:
        if is_missing_value(arg):
            return np.nan
        try:
            value = float(arg)
        except Exception:
            return np.nan
        if not math.isfinite(value):
            return np.nan
        clean_args.append(value)
    try:
        value = func(*clean_args)
    except Exception:
        return np.nan
    return value if math.isfinite(value) else np.nan


def sign(value: float) -> float:
    if value > 0:
        return 1.0
    if value < 0:
        return -1.0
    return 0.0


def qnorm(value: float) -> float:
    if not 0 < value < 1:
        raise ValueError("qnorm input out of range")
    return NORMAL_DIST.inv_cdf(value)


def pnorm(value: float) -> float:
    return NORMAL_DIST.cdf(value)


def pval_from_neglog10p(value: float) -> float:
    return 10 ** (-value)


def pval_from_log10p(value: float) -> float:
    return 10**value


def zscore_from_beta_se(beta: float, se: float) -> float:
    return beta / se


def zscore_from_pval_beta(p_value: float, beta: float) -> float:
    return sign(beta) * abs(qnorm(p_value / 2.0))


def pval_from_zscore(z_value: float) -> float:
    return 2.0 * pnorm(-abs(z_value))


def beta_from_zscore_se(z_value: float, se: float) -> float:
    return z_value * se


def beta_from_zscore_N_af(z_value: float, n_value: float, af_value: float) -> float:
    return z_value / math.sqrt(2.0 * af_value * (1.0 - af_value) * (n_value + z_value**2))


def se_from_zscore_beta(z_value: float, beta: float) -> float:
    return beta / z_value


def se_from_zscore_N_af(z_value: float, n_value: float, af_value: float) -> float:
    return 1.0 / math.sqrt(2.0 * af_value * (1.0 - af_value) * (n_value + z_value**2))


def n_from_zscore_beta_af(z_value: float, beta: float, af_value: float) -> float:
    return z_value**2 / (2.0 * af_value * (1.0 - af_value) * beta**2) - z_value**2


def beta_from_oddsratio(odds_ratio: float) -> float:
    return math.log(odds_ratio)


def se_from_ORu95_ORl95(or_u95: float, or_l95: float) -> float:
    return (math.log(or_u95) - math.log(or_l95)) / (2.0 * qnorm(0.975))


def neff_from_nca_nco(n_case: float, n_control: float) -> float:
    return 4.0 / ((1.0 / n_case) + (1.0 / n_control))


def af_from_case_control(case_af: float, control_af: float, n_case: float, n_control: float) -> float:
    return (case_af * n_case + control_af * n_control) / (n_case + n_control)


def apply_constant_rule(sumstats: pd.DataFrame, output_column: str, value: float | None, *, fill_mode: str) -> pd.DataFrame:
    if value is None:
        return sumstats
    if fill_mode == "column":
        if output_column in sumstats.columns:
            return sumstats
        sumstats[output_column] = pd.Series([value] * len(sumstats), index=sumstats.index, dtype=float)
        return sumstats
    if output_column not in sumstats.columns:
        sumstats[output_column] = pd.Series([np.nan] * len(sumstats), index=sumstats.index, dtype=float)
    mask = sumstats[output_column].isna()
    if mask.any():
        sumstats.loc[mask, output_column] = value
    return sumstats


def apply_formula_rule(
    sumstats: pd.DataFrame,
    output_column: str,
    input_columns: Sequence[str],
    formula: Callable[..., float],
    *,
    fill_mode: str,
) -> pd.DataFrame:
    if fill_mode == "column":
        if output_column in sumstats.columns:
            return sumstats
        if any(column not in sumstats.columns for column in input_columns):
            return sumstats
        sumstats[output_column] = pd.Series(
            [
                safe_formula(formula, *(sumstats[column].iloc[idx] for column in input_columns))
                for idx in range(len(sumstats))
            ],
            index=sumstats.index,
            dtype=float,
        )
        return sumstats
    if any(column not in sumstats.columns for column in input_columns):
        return sumstats
    if output_column not in sumstats.columns:
        sumstats[output_column] = pd.Series([np.nan] * len(sumstats), index=sumstats.index, dtype=float)
    mask = sumstats[output_column].isna()
    if not mask.any():
        return sumstats
    sumstats.loc[mask, output_column] = [
        safe_formula(formula, *(sumstats[column].iloc[idx] for column in input_columns))
        for idx in sumstats.index[mask]
    ]
    return sumstats


def maybe_add_from_oaf(sumstats: pd.DataFrame, metadata: Dict[str, object], output_column: str, input_column: str) -> pd.DataFrame:
    if output_column not in sumstats.columns and input_column in sumstats.columns:
        sumstats[output_column] = 1.0 - sumstats[input_column]
        metadata[f"col_{output_column}"] = output_column
    return sumstats


def harmonize_clean_sumstats(
    header: Sequence[str],
    rows: Sequence[Sequence[object]],
    metadata: Dict[str, object],
    *,
    fill_mode: str,
    use_af_inference: bool,
    warn: Callable[[str], None] | None = None,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    if fill_mode not in {"column", "row"}:
        raise ValueError(f"unsupported fill-mode: {fill_mode!r}")

    metadata_out = deepcopy(metadata)
    model = normalize_model(metadata_out, warn)
    resolved = resolve_clean_metadata_columns(header, metadata_out, include_variant_columns=False)

    data: Dict[str, List[object]] = {}
    seen_canonical: Dict[str, str] = {}
    for key, idx in resolved.items():
        canonical = PAYLOAD_METADATA_TO_CANONICAL[key]
        existing_key = seen_canonical.get(canonical)
        if existing_key is not None and existing_key != key:
            raise ValueError(f"metadata-declared input columns collapse to the same canonical output name {canonical!r}")
        seen_canonical[canonical] = key
        data[canonical] = [normalize_cell_value(row[idx]) for row in rows]
        metadata_out[key] = canonical

    sumstats = pd.DataFrame(data)
    for key in PAYLOAD_METADATA_TO_CANONICAL:
        if key not in resolved:
            metadata_out.pop(key, None)

    sumstats = drop_all_missing_columns(sumstats, metadata_out)

    neglog10p = bool(find_metadata_value(metadata_out, "stats_neglog10P"))
    log10p = bool(find_metadata_value(metadata_out, "stats_log10P"))
    if neglog10p and log10p:
        raise ValueError("stats_neglog10P and stats_log10P cannot both be true")
    if (neglog10p or log10p) and "P" not in sumstats.columns:
        raise ValueError("col_P is required when stats_neglog10P or stats_log10P is true")
    if "P" in sumstats.columns:
        p_values = pd.to_numeric(sumstats["P"], errors="coerce")
        if neglog10p:
            p_values = np.power(10.0, -p_values)
        elif log10p:
            p_values = np.power(10.0, p_values)
        sumstats["P"] = p_values

    sumstats = validate_numeric_ranges(sumstats)

    sumstats = maybe_add_from_oaf(sumstats, metadata_out, "EAF", "OAF")
    sumstats = maybe_add_from_oaf(sumstats, metadata_out, "CaseEAF", "CaseOAF")
    sumstats = maybe_add_from_oaf(sumstats, metadata_out, "ControlEAF", "ControlOAF")
    for old_key, old_column in (("col_OAF", "OAF"), ("col_CaseOAF", "CaseOAF"), ("col_ControlOAF", "ControlOAF")):
        if old_column in sumstats.columns:
            sumstats = sumstats.drop(columns=[old_column])
        metadata_out.pop(old_key, None)

    stats_total_n = metadata_numeric_value(metadata_out, "stats_TotalN") if model == "linear" else None
    stats_case_n = metadata_numeric_value(metadata_out, "stats_CaseN") if model == "logistic" else None
    stats_control_n = metadata_numeric_value(metadata_out, "stats_ControlN") if model == "logistic" else None

    if model == "linear":
        sumstats = apply_constant_rule(sumstats, "N", stats_total_n, fill_mode=fill_mode)
    if model == "logistic":
        sumstats = apply_constant_rule(sumstats, "CaseN", stats_case_n, fill_mode=fill_mode)
        sumstats = apply_constant_rule(sumstats, "ControlN", stats_control_n, fill_mode=fill_mode)
        sumstats = apply_formula_rule(sumstats, "N", ("CaseN", "ControlN"), neff_from_nca_nco, fill_mode=fill_mode)
        sumstats = apply_formula_rule(
            sumstats,
            "EAF",
            ("CaseEAF", "ControlEAF", "CaseN", "ControlN"),
            af_from_case_control,
            fill_mode=fill_mode,
        )
        sumstats = apply_formula_rule(sumstats, "BETA", ("OR",), beta_from_oddsratio, fill_mode=fill_mode)
        sumstats = apply_formula_rule(sumstats, "SE", ("ORU95", "ORL95"), se_from_ORu95_ORl95, fill_mode=fill_mode)

    sumstats = apply_formula_rule(sumstats, "Z", ("P", "BETA"), zscore_from_pval_beta, fill_mode=fill_mode)
    sumstats = apply_formula_rule(sumstats, "Z", ("BETA", "SE"), zscore_from_beta_se, fill_mode=fill_mode)
    sumstats = apply_formula_rule(sumstats, "P", ("Z",), pval_from_zscore, fill_mode=fill_mode)
    sumstats = apply_formula_rule(sumstats, "SE", ("Z", "BETA"), se_from_zscore_beta, fill_mode=fill_mode)
    sumstats = apply_formula_rule(sumstats, "BETA", ("Z", "SE"), beta_from_zscore_se, fill_mode=fill_mode)

    if use_af_inference:
        sumstats = apply_formula_rule(sumstats, "N", ("Z", "BETA", "EAF"), n_from_zscore_beta_af, fill_mode=fill_mode)
        sumstats = apply_formula_rule(sumstats, "BETA", ("Z", "N", "EAF"), beta_from_zscore_N_af, fill_mode=fill_mode)
        sumstats = apply_formula_rule(sumstats, "SE", ("Z", "N", "EAF"), se_from_zscore_N_af, fill_mode=fill_mode)

    sumstats = validate_numeric_ranges(sumstats)

    for key, canonical in PAYLOAD_METADATA_TO_CANONICAL.items():
        if canonical in sumstats.columns:
            metadata_out[key] = canonical
        else:
            metadata_out.pop(key, None)

    ordered_columns = [column for column in PAYLOAD_OUTPUT_ORDER if column in sumstats.columns]
    return sumstats.loc[:, ordered_columns], metadata_out
