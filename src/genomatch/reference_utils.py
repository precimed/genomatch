from __future__ import annotations

import os
import shutil
from functools import lru_cache
from pathlib import Path
from typing import Dict, Tuple

DEFAULT_MATCH_CONFIG = Path("/ref/config.yaml")


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


def _resolve_fasta_entry(raw_entry: object, *, label: str) -> Path:
    if isinstance(raw_entry, dict):
        return _resolve_config_path_value(raw_entry.get("fasta"), label=label)
    return _resolve_config_path_value(raw_entry, label=label)


def resolve_internal_reference_fasta(genome_build: str) -> Path:
    config = load_match_config()
    references = config.get("references")
    if not isinstance(references, dict):
        raise ValueError("config is missing 'references' section")
    ucsc = references.get("ucsc")
    if not isinstance(ucsc, dict):
        raise ValueError("config is missing 'references.ucsc' section")
    build_entry = ucsc.get(genome_build)
    if build_entry is None:
        raise ValueError(f"config has no internal UCSC reference for genome build {genome_build!r}")
    return _resolve_fasta_entry(
        build_entry,
        label=f"internal UCSC reference FASTA for genome build {genome_build!r}",
    )


def resolve_liftover_chain(source_build: str, target_build: str) -> Path:
    config = load_match_config()
    chain_key_map = {
        ("GRCh37", "GRCh38"): "hg19ToHg38",
        ("GRCh38", "GRCh37"): "hg38ToHg19",
    }
    chain_name = chain_key_map.get((source_build, target_build))
    if chain_name is None:
        raise ValueError(f"unsupported liftover direction: {source_build!r} -> {target_build!r}")
    chain_block = config.get("chain")
    if not isinstance(chain_block, dict):
        raise ValueError("config is missing 'chain' section")
    if chain_name not in chain_block:
        raise ValueError(f"config has no chain.{chain_name} entry")
    return _resolve_config_path_value(
        chain_block.get(chain_name),
        label=f"chain.{chain_name}",
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


def fetch_reference_base(fasta_path: Path, contig: str, pos: int) -> str:
    if pos <= 0:
        return ""
    fasta = open_indexed_fasta(str(fasta_path))
    try:
        references = set(fasta.references)
    except Exception as exc:
        raise ValueError(f"failed to inspect indexed FASTA {fasta_path}: {exc}") from exc
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
