from __future__ import annotations

import os
import shutil
from functools import lru_cache
from pathlib import Path
from typing import Dict, Sequence, Tuple

from .vtable_utils import SUPPORTED_GENOME_BUILDS

DEFAULT_MATCH_CONFIG = Path("/ref/config.yaml")
REFERENCE_ACCESS_MODES = {"BULK", "LEGACY"}
_REFERENCE_GENOME_BUILDS = SUPPORTED_GENOME_BUILDS - {"unknown"}
PRIMARY_UCSC_CONTIGS = tuple([f"chr{idx}" for idx in range(1, 23)] + ["chrX", "chrY", "chrM"])


def _load_yaml(path: Path) -> Dict[str, object]:
    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise ValueError("PyYAML is required to read config files") from exc
    data = yaml.safe_load(path.read_text(encoding="utf-8", errors="replace"))
    if not isinstance(data, dict):
        raise ValueError(f"YAML document must be a mapping: {path}")
    return data


def resolve_config_path() -> Path:
    override = os.environ.get("MATCH_CONFIG")
    if override:
        return Path(override)
    if DEFAULT_MATCH_CONFIG.exists():
        return DEFAULT_MATCH_CONFIG
    raise ValueError(
        "MATCH_CONFIG is required unless /ref/config.yaml exists; point it to a user-provided config.yaml"
    )


def load_match_config() -> Dict[str, object]:
    path = resolve_config_path()
    return _load_match_config_cached(str(path))


@lru_cache(maxsize=None)
def _load_match_config_cached(path_raw: str) -> Dict[str, object]:
    path = Path(path_raw)
    if not path.exists():
        raise ValueError(f"match config not found: {path}")
    return _load_yaml(path)


def _resolve_config_path_value(raw_value: object, *, label: str) -> Path:
    if not isinstance(raw_value, str) or not raw_value:
        raise ValueError(f"config entry {label} must be a non-empty path string")
    path = Path(raw_value)
    if not path.is_absolute():
        path = resolve_config_path().parent / path
    if not path.exists():
        raise ValueError(f"{label} not found: {path}")
    return path


def _reject_legacy_config_shape(config: Dict[str, object]) -> None:
    if "references" in config or "chain" in config:
        raise ValueError(
            "legacy match config shape is no longer supported; update config.yaml to use "
            "'builds.<build>.ucsc_fasta' and list-valued 'liftover' entries"
        )


def _require_supported_reference_build(genome_build: str, *, label: str) -> None:
    if genome_build not in _REFERENCE_GENOME_BUILDS:
        accepted = ", ".join(sorted(_REFERENCE_GENOME_BUILDS))
        raise ValueError(f"unsupported {label}: {genome_build!r}; expected one of: {accepted}")


def _resolve_ucsc_fasta_entry(raw_entry: object, *, genome_build: str) -> Path:
    if not isinstance(raw_entry, dict):
        raise ValueError(f"config builds.{genome_build} must be a mapping with ucsc_fasta")
    return _resolve_config_path_value(
        raw_entry.get("ucsc_fasta"),
        label=f"builds.{genome_build}.ucsc_fasta",
    )


def validate_primary_ucsc_fasta(fasta_path: Path, *, genome_build: str) -> None:
    references = _indexed_fasta_references(str(fasta_path))
    missing = [contig for contig in PRIMARY_UCSC_CONTIGS if contig not in references]
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "" if len(missing) <= 5 else f", ... ({len(missing)} missing total)"
        raise ValueError(
            f"configured UCSC FASTA for genome build {genome_build!r} is missing required primary contigs: "
            f"{preview}{suffix}"
        )


def resolve_internal_reference_fasta(genome_build: str) -> Path:
    _require_supported_reference_build(genome_build, label="genome build")
    config = load_match_config()
    _reject_legacy_config_shape(config)
    builds = config.get("builds")
    if not isinstance(builds, dict):
        raise ValueError("config is missing 'builds' section")
    build_entry = builds.get(genome_build)
    if build_entry is None:
        raise ValueError(f"config has no builds.{genome_build}.ucsc_fasta entry")
    fasta_path = _resolve_ucsc_fasta_entry(build_entry, genome_build=genome_build)
    validate_primary_ucsc_fasta(fasta_path, genome_build=genome_build)
    return fasta_path


def _configured_liftover_edges(config: Dict[str, object]) -> Dict[Tuple[str, str], object]:
    liftover = config.get("liftover")
    if not isinstance(liftover, list):
        raise ValueError("config is missing list-valued 'liftover' section")
    edges: Dict[Tuple[str, str], object] = {}
    duplicate_edges: set[Tuple[str, str]] = set()
    for index, raw_edge in enumerate(liftover):
        if not isinstance(raw_edge, dict):
            raise ValueError(f"config liftover[{index}] must be a mapping")
        source = raw_edge.get("source")
        target = raw_edge.get("target")
        if not isinstance(source, str) or not source:
            raise ValueError(f"config liftover[{index}].source must be a non-empty string")
        if not isinstance(target, str) or not target:
            raise ValueError(f"config liftover[{index}].target must be a non-empty string")
        chain = raw_edge.get("chain")
        if not isinstance(chain, str) or not chain:
            raise ValueError(f"config liftover[{index}].chain must be a non-empty path string")
        _require_supported_reference_build(source, label=f"liftover[{index}].source")
        _require_supported_reference_build(target, label=f"liftover[{index}].target")
        key = (source, target)
        if key in edges:
            duplicate_edges.add(key)
            continue
        edges[key] = chain
    if duplicate_edges:
        rendered = ", ".join(f"{source}->{target}" for source, target in sorted(duplicate_edges))
        raise ValueError(f"config has duplicate liftover edges: {rendered}")
    return edges


def resolve_liftover_chain(source_build: str, target_build: str) -> Path:
    _require_supported_reference_build(source_build, label="source genome build")
    _require_supported_reference_build(target_build, label="target genome build")
    config = load_match_config()
    _reject_legacy_config_shape(config)
    edges = _configured_liftover_edges(config)
    chain_value = edges.get((source_build, target_build))
    if chain_value is None:
        raise ValueError(f"unsupported liftover direction: {source_build!r} -> {target_build!r}")
    return _resolve_config_path_value(
        chain_value,
        label=f"liftover edge {source_build!r}->{target_build!r} chain",
    )


def resolve_liftover_assets(source_build: str, target_build: str) -> Tuple[Path, Path, Path]:
    source_fasta = resolve_internal_reference_fasta(source_build)
    target_fasta = resolve_internal_reference_fasta(target_build)
    chain_path = resolve_liftover_chain(source_build, target_build)
    return source_fasta, target_fasta, chain_path


def resolve_bcftools_binary() -> str:
    override = os.environ.get("MATCH_BCFTOOLS")
    if override:
        if shutil.which(override) or Path(override).exists():
            return override
        raise ValueError(f"bcftools not found: {override}")
    found = shutil.which("bcftools")
    if found:
        return found
    raise ValueError("bcftools not found; install it or set MATCH_BCFTOOLS")


def resolve_reference_access_mode() -> str:
    raw_mode = os.environ.get("MATCH_REFERENCE_ACCESS_MODE")
    if raw_mode is None or raw_mode.strip() == "":
        return "BULK"
    mode = raw_mode.strip().upper()
    if mode not in REFERENCE_ACCESS_MODES:
        accepted = ", ".join(sorted(REFERENCE_ACCESS_MODES))
        raise ValueError(
            f"unsupported MATCH_REFERENCE_ACCESS_MODE {raw_mode!r}; expected one of: {accepted}"
        )
    return mode


@lru_cache(maxsize=None)
def open_indexed_fasta(path: str):
    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise ValueError("pysam is required for indexed FASTA access") from exc
    fasta_path = Path(path)
    fai_path = Path(str(fasta_path) + ".fai")
    if not fai_path.exists():
        raise ValueError(f"FASTA index not found: {fai_path}")
    try:
        return pysam.FastaFile(str(fasta_path))
    except Exception as exc:
        raise ValueError(f"failed to open indexed FASTA {fasta_path}: {exc}") from exc


@lru_cache(maxsize=None)
def _indexed_fasta_references(path: str) -> frozenset[str]:
    fasta = open_indexed_fasta(path)
    fasta_path = Path(path)
    try:
        return frozenset(fasta.references)
    except Exception as exc:
        raise ValueError(f"failed to inspect indexed FASTA {fasta_path}: {exc}") from exc


@lru_cache(maxsize=None)
def _fetch_contig_sequence(path: str, contig: str) -> str:
    fasta = open_indexed_fasta(path)
    fasta_path = Path(path)
    references = _indexed_fasta_references(path)
    if contig not in references:
        return ""
    try:
        return fasta.fetch(contig).upper()
    except Exception as exc:
        raise ValueError(
            f"failed to fetch reference contig from indexed FASTA {fasta_path} at {contig}: {exc}"
        ) from exc


def fetch_reference_base_legacy(fasta_path: Path, contig: str, pos: int) -> str:
    if pos <= 0:
        return ""
    path_str = str(fasta_path)
    fasta = open_indexed_fasta(path_str)
    references = _indexed_fasta_references(path_str)
    if contig not in references:
        return ""
    try:
        contig_length = fasta.get_reference_length(contig)
    except Exception as exc:
        raise ValueError(f"failed to inspect contig {contig!r} in indexed FASTA {fasta_path}: {exc}") from exc
    if pos > contig_length:
        return ""
    try:
        seq = fasta.fetch(contig, pos - 1, pos)
    except Exception as exc:
        raise ValueError(
            f"failed to fetch reference base from indexed FASTA {fasta_path} at {contig}:{pos}: {exc}"
        ) from exc
    if len(seq) != 1:
        return ""
    return seq.upper()


def fetch_reference_bases_bulk(
    fasta_path: Path,
    queries: Sequence[Tuple[str, int]],
) -> Dict[Tuple[str, int], str]:
    path_str = str(fasta_path)
    out: Dict[Tuple[str, int], str] = {}
    positions_by_contig: Dict[str, list[int]] = {}
    for contig, pos in queries:
        if pos <= 0:
            out[(contig, pos)] = ""
            continue
        positions_by_contig.setdefault(contig, []).append(pos)
    for contig, positions in positions_by_contig.items():
        sequence = _fetch_contig_sequence(path_str, contig)
        if not sequence:
            for pos in positions:
                out[(contig, pos)] = ""
            continue
        contig_length = len(sequence)
        for pos in positions:
            if pos > contig_length:
                out[(contig, pos)] = ""
                continue
            out[(contig, pos)] = sequence[pos - 1]
    return out


def fetch_reference_bases(
    fasta_path: Path,
    queries: Sequence[Tuple[str, int]],
) -> Dict[Tuple[str, int], str]:
    if resolve_reference_access_mode() == "LEGACY":
        return {
            (contig, pos): fetch_reference_base_legacy(fasta_path, contig, pos)
            for contig, pos in queries
        }
    return fetch_reference_bases_bulk(fasta_path, queries)


def fetch_reference_base(fasta_path: Path, contig: str, pos: int) -> str:
    mode = resolve_reference_access_mode()
    if mode == "LEGACY":
        return fetch_reference_base_legacy(fasta_path, contig, pos)
    sequence = _fetch_contig_sequence(str(fasta_path), contig)
    if pos <= 0:
        return ""
    if not sequence:
        return ""
    if pos > len(sequence):
        return ""
    return sequence[pos - 1]
