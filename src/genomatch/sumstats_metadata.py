from __future__ import annotations

import re
from typing import Dict, List, Optional, Sequence, Tuple

from .sumstats_utils import find_metadata_value


NON_ALNUM = re.compile(r"[^A-Za-z0-9]")


def normalize_header_token(value: str) -> str:
    return NON_ALNUM.sub("", value).lower().strip()


def metadata_column_name(metadata: Dict[str, object], key: str) -> Optional[str]:
    value = find_metadata_value(metadata, key)
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"column mapping for {key} must be a string column name")
    raw = value.strip()
    return raw or None


def resolve_metadata_column_indices(
    header: Sequence[str],
    metadata: Dict[str, object],
    keys: Sequence[str],
) -> Dict[str, int]:
    normalized_header = [normalize_header_token(str(column)) for column in header]
    resolved: Dict[str, int] = {}
    for key in keys:
        raw_value = metadata_column_name(metadata, key)
        if raw_value is None:
            continue
        normalized_value = normalize_header_token(raw_value)
        matches = [idx for idx, token in enumerate(normalized_header) if token == normalized_value]
        if not matches:
            raise ValueError(f"column not found for {key}: {raw_value}")
        if len(matches) > 1:
            raise ValueError(
                f"column mapping for {key} could not be identified uniquely after header normalization"
            )
        resolved[key] = matches[0]
    return resolved


def retained_payload_columns(
    preview_header: Sequence[str],
    metadata: Dict[str, object],
    *,
    payload_keys: Sequence[str],
    reserved_output_columns: Sequence[str],
) -> Tuple[List[Tuple[str, int, str]], Dict[str, str]]:
    resolved = resolve_metadata_column_indices(preview_header, metadata, payload_keys)
    retained: List[Tuple[str, int, str]] = []
    output_by_key: Dict[str, str] = {}
    seen_output: Dict[str, str] = {}
    for key, idx in sorted(resolved.items(), key=lambda item: item[1]):
        output_name = metadata_column_name(metadata, key)
        if output_name is None:
            continue
        if output_name in reserved_output_columns:
            raise ValueError(
                f"metadata-declared payload column {key} would collide with canonical output column {output_name!r}"
            )
        existing = seen_output.get(output_name)
        if existing is not None and existing != key:
            raise ValueError(f"metadata-declared payload columns produce duplicate output name {output_name!r}")
        seen_output[output_name] = key
        retained.append((key, idx, output_name))
        output_by_key[key] = output_name
    return retained, output_by_key
