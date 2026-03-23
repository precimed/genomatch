from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
from typing import Iterable, List

import pytest

pysam = pytest.importorskip("pysam")

from utils import REPO_ROOT, run_cmd, write_match_config


REF_DIR = REPO_ROOT / "ref"
GRCH37_UCSC_SOURCE_FASTA = REF_DIR / "ucsc" / "GRCh37" / "hg19.p13.plusMT.no_alt_analysis_set.fa"
GRCH38_UCSC_SOURCE_FASTA = REF_DIR / "ucsc" / "GRCh38" / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CHAIN_37_TO_38 = REF_DIR / "chain" / "hg19ToHg38.over.chain.gz"
CHAIN_38_TO_37 = REF_DIR / "chain" / "hg38ToHg19.over.chain.gz"
MAPPING_PATH = REPO_ROOT / "tests" / "data" / "generated_reference" / "dbsnp_cleansumstat_reference_GRCh38_GRCh37.txt"

ALLOWED_BASES = {"A", "C", "G", "T"}


@dataclass(frozen=True)
class MappingEntry:
    chrom: str
    bp38: str
    bp37: str
    rsid: str
    a1: str
    a2: str


def ensure_real_liftover_assets() -> None:
    required = [
        GRCH37_UCSC_SOURCE_FASTA,
        GRCH37_UCSC_SOURCE_FASTA.with_suffix(".fa.fai"),
        GRCH38_UCSC_SOURCE_FASTA,
        GRCH38_UCSC_SOURCE_FASTA.with_suffix(".fna.fai"),
        CHAIN_37_TO_38,
        CHAIN_38_TO_37,
        MAPPING_PATH,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        pytest.skip(f"missing real liftover assets: {', '.join(missing)}")


def resolve_bcftools_with_liftover() -> str:
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools not found in PATH")
    probe = run_cmd([bcftools, "+liftover"])
    text = (probe.stdout + probe.stderr).lower()
    if 'plugin "liftover" was not found' in text or 'could not load "liftover"' in text:
        pytest.skip("bcftools +liftover plugin not available in PATH")
    return bcftools


def fetch_ref_base(fasta: "pysam.FastaFile", contig: str, pos: int) -> str:
    try:
        base = fasta.fetch(contig, pos - 1, pos).upper()
    except (KeyError, OSError, ValueError):
        return ""
    return base if len(base) == 1 else ""


def canonical_alleles_for_build(entry: MappingEntry, *, build_name: str) -> tuple[str, str]:
    fasta_path = GRCH37_UCSC_SOURCE_FASTA if build_name == "GRCh37" else GRCH38_UCSC_SOURCE_FASTA
    chrom = f"chr{entry.chrom}"
    pos = int(entry.bp37 if build_name == "GRCh37" else entry.bp38)
    fasta = pysam.FastaFile(str(fasta_path))
    try:
        ref_base = fetch_ref_base(fasta, chrom, pos)
    finally:
        fasta.close()
    if ref_base == entry.a1 and ref_base != entry.a2:
        return (entry.a2, entry.a1)
    if ref_base == entry.a2 and ref_base != entry.a1:
        return (entry.a1, entry.a2)
    raise ValueError(f"mapping entry is not reference-anchored for {build_name}: {entry}")


def select_mapping_entries(*, chrom: str = "10", max_rows: int = 6) -> List[MappingEntry]:
    ensure_real_liftover_assets()
    selected: List[MappingEntry] = []
    ref37 = pysam.FastaFile(str(GRCH37_UCSC_SOURCE_FASTA))
    ref38 = pysam.FastaFile(str(GRCH38_UCSC_SOURCE_FASTA))
    try:
        with open(MAPPING_PATH, "r", encoding="utf-8", newline="\n") as handle:
            for line in handle:
                if not line.strip():
                    continue
                coord38, coord37, rsid, a1, a2 = line.rstrip("\n").split()[:5]
                chrom38, bp38 = coord38.split(":")
                chrom37, bp37 = coord37.split(":")
                if chrom38 != chrom or chrom37 != chrom:
                    continue
                if a1 not in ALLOWED_BASES or a2 not in ALLOWED_BASES:
                    continue
                base37 = fetch_ref_base(ref37, f"chr{chrom37}", int(bp37))
                base38 = fetch_ref_base(ref38, f"chr{chrom38}", int(bp38))
                if base37 not in {a1, a2} or base38 not in {a1, a2}:
                    continue
                selected.append(MappingEntry(chrom, bp38, bp37, rsid, a1, a2))
                if len(selected) >= max_rows:
                    break
    finally:
        ref37.close()
        ref38.close()
    if len(selected) < max_rows:
        pytest.skip(f"not enough mapping entries for chr{chrom}: {len(selected)} found")
    return selected


def _write_fasta_record(handle, name: str, sequence: str) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), 60):
        handle.write(sequence[start : start + 60] + "\n")


def write_ucsc_subset_fasta(
    output_path: Path,
    *,
    source_fasta: Path,
    chroms: Iterable[str],
    source_labels_are_ucsc: bool,
) -> Path:
    fasta = pysam.FastaFile(str(source_fasta))
    try:
        chrom_list = sorted(set(chroms), key=int)
        with open(output_path, "w", encoding="utf-8", newline="\n") as handle:
            for chrom in chrom_list:
                source_contig = f"chr{chrom}" if source_labels_are_ucsc else chrom
                target_contig = f"chr{chrom}"
                _write_fasta_record(handle, target_contig, fasta.fetch(source_contig))
    finally:
        fasta.close()
    pysam.faidx(str(output_path))
    return output_path


def write_real_match_config(tmp_path: Path, *, chroms: Iterable[str]) -> Path:
    ensure_real_liftover_assets()
    config = tmp_path / "config.yaml"
    grch37_ucsc = write_ucsc_subset_fasta(
        tmp_path / "GRCh37.ucsc.fa",
        source_fasta=GRCH37_UCSC_SOURCE_FASTA,
        chroms=chroms,
        source_labels_are_ucsc=True,
    )
    grch38_ucsc = write_ucsc_subset_fasta(
        tmp_path / "GRCh38.ucsc.fa",
        source_fasta=GRCH38_UCSC_SOURCE_FASTA,
        chroms=chroms,
        source_labels_are_ucsc=True,
    )
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=CHAIN_37_TO_38,
        chain38to37=CHAIN_38_TO_37,
    )
    return config
