from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
MATCH_SCRIPT_DIR = REPO_ROOT / "src" / "genomatch"
MATCH_TEST_DIR = REPO_ROOT / "tests"
SYNTHETIC_DIR = MATCH_TEST_DIR / "synthetic"
BASE_PREFIX = SYNTHETIC_DIR / "base"
BASE_VCF = SYNTHETIC_DIR / "base.vcf"

if str(MATCH_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(MATCH_SCRIPT_DIR))

from bfile_utils import read_bed_selected, write_bed_matrix

PLINK_BIN = os.environ.get("PLINK_BIN", "plink")
PLINK2_BIN = os.environ.get("PLINK2_BIN", "plink2")


@dataclass(frozen=True)
class BimRow:
    chrom: str
    snp: str
    cm: str
    bp: str
    a1: str
    a2: str


@dataclass
class BaseData:
    rows: List[BimRow]
    genos: List[List[int]]
    n_samples: int

def resolve_plink() -> List[str]:
    if shutil.which(PLINK_BIN) is None:
        raise RuntimeError(f"plink not found: {PLINK_BIN}")
    return [PLINK_BIN]


def resolve_plink2() -> List[str]:
    if shutil.which(PLINK2_BIN) is None:
        raise RuntimeError(f"plink2 not found: {PLINK2_BIN}")
    return [PLINK2_BIN]


def run_cmd(cmd: Iterable[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(c) for c in cmd],
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=False,
    )


def run_py(script_name: str, *args: str) -> subprocess.CompletedProcess:
    return run_cmd([sys.executable, str(MATCH_SCRIPT_DIR / script_name), *map(str, args)])


def run_py_with_env(script_name: str, env: dict[str, str], *args: str) -> subprocess.CompletedProcess:
    full_env = os.environ.copy()
    full_env.update(env)
    return subprocess.run(
        [sys.executable, str(MATCH_SCRIPT_DIR / script_name), *map(str, args)],
        capture_output=True,
        text=True,
        check=False,
        env=full_env,
    )


def ensure_base_bfile() -> None:
    SYNTHETIC_DIR.mkdir(parents=True, exist_ok=True)
    bed = BASE_PREFIX.with_suffix(".bed")
    bim = BASE_PREFIX.with_suffix(".bim")
    fam = BASE_PREFIX.with_suffix(".fam")
    if bed.exists() and bim.exists() and fam.exists():
        return

    generator = MATCH_TEST_DIR / "generate_synthetic_data.py"
    result = run_cmd([sys.executable, str(generator), "--out-dir", str(SYNTHETIC_DIR), "--vcf-name", BASE_VCF.name])
    if result.returncode != 0:
        raise RuntimeError(f"Synthetic data generation failed: {result.stderr}\n{result.stdout}")

    result = run_cmd(resolve_plink() + ["--vcf", str(BASE_VCF), "--double-id", "--make-bed", "--out", str(BASE_PREFIX)])
    if result.returncode != 0:
        raise RuntimeError(f"plink failed to create base bfile: {result.stderr}\n{result.stdout}")


def copy_bfile(src_prefix: Path, dst_prefix: Path) -> None:
    for ext in (".bed", ".bim", ".fam"):
        shutil.copy2(src_prefix.with_suffix(ext), dst_prefix.with_suffix(ext))


def read_bim(path: Path) -> List[BimRow]:
    rows: List[BimRow] = []
    with open(path, "r", newline="\n") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 6:
                raise ValueError(f"Invalid BIM row: {line}")
            rows.append(BimRow(*parts[:6]))
    return rows


def write_bim(path: Path, rows: List[BimRow]) -> None:
    with open(path, "w", newline="\n") as handle:
        for row in rows:
            handle.write("\t".join([row.chrom, row.snp, row.cm, row.bp, row.a1, row.a2]) + "\n")


def read_fam_sample_count(path: Path) -> int:
    count = 0
    with open(path, "r", newline="\n") as handle:
        for line in handle:
            if line.strip():
                count += 1
    return count


def read_bed_matrix(bed_path: Path, n_samples: int, n_snps: int) -> List[List[int]]:
    with open(bed_path, "rb") as handle:
        header = handle.read(3)
        if header != b"\x6c\x1b\x01":
            raise ValueError("Unsupported .bed format (expected SNP-major)")
        bytes_per_snp = (n_samples + 3) // 4
        data = handle.read()

    expected_bytes = bytes_per_snp * n_snps
    if len(data) < expected_bytes:
        raise ValueError(".bed file is shorter than expected")

    matrix: List[List[int]] = []
    offset = 0
    for _ in range(n_snps):
        chunk = data[offset : offset + bytes_per_snp]
        offset += bytes_per_snp
        genos: List[int] = []
        for byte in chunk:
            for shift in range(0, 8, 2):
                if len(genos) == n_samples:
                    break
                genos.append((byte >> shift) & 0b11)
        matrix.append(genos)
    return matrix


def swap_genotypes(genos: List[int]) -> List[int]:
    swapped = []
    for geno in genos:
        if geno == 0:
            swapped.append(3)
        elif geno == 3:
            swapped.append(0)
        else:
            swapped.append(geno)
    return swapped


def write_empty_bfile(prefix: Path, fam_source: Path) -> None:
    shutil.copy2(fam_source, prefix.with_suffix(".fam"))
    with open(prefix.with_suffix(".bim"), "w", newline="\n") as handle:
        handle.write("")
    with open(prefix.with_suffix(".bed"), "wb") as handle:
        handle.write(b"\x6c\x1b\x01")


def apply_bim_updates(rows: List[BimRow], updates: List[Tuple[int, dict]]) -> List[BimRow]:
    updated = list(rows)
    for idx, changes in updates:
        row = updated[idx]
        data = {
            "chrom": row.chrom,
            "snp": row.snp,
            "cm": row.cm,
            "bp": row.bp,
            "a1": row.a1,
            "a2": row.a2,
        }
        data.update(changes)
        updated[idx] = BimRow(
            data["chrom"],
            data["snp"],
            data["cm"],
            data["bp"],
            data["a1"],
            data["a2"],
        )
    return updated


def max_bp(rows: List[BimRow]) -> int:
    return max(int(row.bp) for row in rows)


def assert_warning_emitted(result: subprocess.CompletedProcess) -> None:
    text = (result.stdout + result.stderr).lower()
    if "warn" not in text and "missing" not in text:
        raise AssertionError("Expected warning was not emitted")


def assert_error_message(result: subprocess.CompletedProcess, keywords: Iterable[str]) -> None:
    text = (result.stdout + result.stderr).lower()
    if not text.strip():
        raise AssertionError("Expected an error message, but none was emitted")
    if keywords and not any(keyword in text for keyword in keywords):
        raise AssertionError(f"Expected error message to include one of: {keywords}")


def write_lines(path: Path, lines: List[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_tsv(path: Path) -> List[List[str]]:
    return [line.rstrip("\n").split("\t") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    lines = []
    for name, seq in sequences.items():
        lines.append(f">{name}")
        lines.append(seq)
    write_lines(path, lines)
    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise RuntimeError("pysam is required for indexed FASTA test fixtures") from exc
    pysam.faidx(str(path))


def write_match_config(
    path: Path,
    *,
    grch37_fasta: Path | None = None,
    grch38_fasta: Path | None = None,
    grch37_ucsc_fasta: Path | None = None,
    grch38_ucsc_fasta: Path | None = None,
    chain37to38: Path | None = None,
    chain38to37: Path | None = None,
) -> None:
    grch37_ucsc_fasta = grch37_ucsc_fasta or grch37_fasta
    grch38_ucsc_fasta = grch38_ucsc_fasta or grch38_fasta
    if grch37_ucsc_fasta is None or grch38_ucsc_fasta is None:
        raise ValueError("UCSC FASTA paths must be provided for GRCh37 and GRCh38")
    lines = [
        "references:",
        "  ucsc:",
        "    GRCh37:",
        f"      fasta: {grch37_ucsc_fasta}",
        "    GRCh38:",
        f"      fasta: {grch38_ucsc_fasta}",
    ]
    if grch37_fasta is not None or grch38_fasta is not None:
        lines.extend(
            [
                "  ncbi:",
                "    GRCh37:",
                f"      fasta: {grch37_fasta or ''}",
                "    GRCh38:",
                f"      fasta: {grch38_fasta or ''}",
            ]
        )
    lines.extend(
        [
            "chain:",
            f"  hg19ToHg38: {chain37to38 or ''}",
            f"  hg38ToHg19: {chain38to37 or ''}",
        ]
    )
    write_lines(path, lines)


def write_bfile(prefix: Path, bim_lines: List[str], fam_lines: List[str], genotypes: List[List[int]]) -> None:
    write_lines(prefix.with_suffix(".bim"), bim_lines)
    write_lines(prefix.with_suffix(".fam"), fam_lines)
    write_bed_matrix(prefix.with_suffix(".bed"), genotypes, len(fam_lines))
