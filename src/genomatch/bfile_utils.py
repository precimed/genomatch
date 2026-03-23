#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from .haploid_utils import HAPLOID_SCHEMA_PATH as HAPLOID_SCHEMA
from .haploid_utils import HAPLOID_CHROMS, expected_ploidy_pair, is_sex_dependent_ploidy


ALLOWED_BASES = {"A", "C", "G", "T"}


@dataclass(frozen=True)
class BimRow:
    chrom: str
    snp: str
    cm: str
    bp: str
    a1: str
    a2: str


def read_bim(path: Path) -> List[BimRow]:
    rows: List[BimRow] = []
    with open(path, "r", newline="\n") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 6:
                raise ValueError(f"invalid BIM row: {line.strip()}")
            rows.append(BimRow(*parts[:6]))
    return rows


def write_bim(path: Path, rows: List[BimRow]) -> None:
    with open(path, "w", newline="\n") as handle:
        for row in rows:
            handle.write("\t".join([row.chrom, row.snp, row.cm, row.bp, row.a1, row.a2]) + "\n")


def read_fam_rows(path: Path) -> List[Tuple[str, str, str, str, str, str]]:
    rows: List[Tuple[str, str, str, str, str, str]] = []
    with open(path, "r", newline="\n") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 6:
                raise ValueError("invalid .fam row")
            rows.append(tuple(parts[:6]))
    return rows


def validate_alleles(rows: List[BimRow], label: str) -> None:
    for row in rows:
        for allele in (row.a1, row.a2):
            allele_u = allele.upper()
            if not allele_u:
                raise ValueError(f"invalid allele code in {label}: {allele}")
            if any(base not in ALLOWED_BASES for base in allele_u):
                raise ValueError(f"invalid allele code in {label}: {allele}")


def read_fam_samples(path: Path) -> Tuple[List[int], List[str]]:
    sexes: List[int] = []
    sample_ids: List[str] = []
    for fid, iid, _father, _mother, sex_raw, _pheno in read_fam_rows(path):
        try:
            sex = int(sex_raw)
        except ValueError:
            sex = 0
        if sex not in (0, 1, 2):
            sex = 0
        sexes.append(sex)
        sample_ids.append(f"{fid}:{iid}")
    return sexes, sample_ids


def ploidy_for_variant(chrom: str, bp_raw: str, *, genome_build: str) -> Tuple[int, int]:
    return expected_ploidy_pair(chrom, bp_raw, genome_build=genome_build)


def count_target_ploidy_genotype_issues(
    genos: Sequence[int],
    sexes: List[int],
    ploidy_male: int,
    ploidy_female: int,
    row: BimRow,
) -> Tuple[int, int, int]:
    haploid_het_incompatible = 0
    absent_nonmissing_incompatible = 0
    unknown_sex_unvalidated = 0
    ploidy_pair = (ploidy_male, ploidy_female)
    for idx, g in enumerate(genos):
        sex = sexes[idx]
        if sex == 0 and is_sex_dependent_ploidy(ploidy_pair):
            unknown_sex_unvalidated += 1
            continue
        ploidy = ploidy_male if sex == 1 else ploidy_female
        if ploidy == 2:
            continue
        if ploidy == 1:
            if g == 2:
                haploid_het_incompatible += 1
        elif ploidy == 0:
            if g != 1:
                absent_nonmissing_incompatible += 1
        else:
            raise ValueError(f"unsupported ploidy value {ploidy} at {row.chrom}:{row.bp}")
    return haploid_het_incompatible, absent_nonmissing_incompatible, unknown_sex_unvalidated


def decode_bed_chunk(chunk: bytes, n_samples: int) -> bytes:
    genos = bytearray(n_samples)
    idx = 0
    for byte in chunk:
        for shift in range(0, 8, 2):
            if idx == n_samples:
                break
            genos[idx] = (byte >> shift) & 0b11
            idx += 1
    return bytes(genos)


def read_bed_selected(
    bed_path: Path, n_samples: int, n_snps: int, selected_indices: Iterable[int]
) -> Dict[int, bytes]:
    indices = sorted(set(selected_indices))
    if not indices:
        return {}
    bytes_per_snp = (n_samples + 3) // 4
    with open(bed_path, "rb") as handle:
        header = handle.read(3)
        if header != b"\x6c\x1b\x01":
            raise ValueError("unsupported .bed format (expected SNP-major)")
        handle.seek(0, 2)
        file_size = handle.tell()
        expected_size = 3 + bytes_per_snp * n_snps
        if file_size < expected_size:
            raise ValueError(".bed file is shorter than expected")
        data: Dict[int, bytes] = {}
        for idx in indices:
            offset = 3 + idx * bytes_per_snp
            handle.seek(offset)
            chunk = handle.read(bytes_per_snp)
            if len(chunk) != bytes_per_snp:
                raise ValueError(".bed file is shorter than expected")
            data[idx] = decode_bed_chunk(chunk, n_samples)
    return data


def write_bed_matrix(path: Path, matrix: Iterable[Sequence[int]], n_samples: int) -> None:
    bytes_per_snp = (n_samples + 3) // 4
    with open(path, "wb") as handle:
        handle.write(b"\x6c\x1b\x01")
        for genos in matrix:
            if len(genos) != n_samples:
                raise ValueError("genotype row length does not match sample count")
            out = bytearray()
            for i in range(0, n_samples, 4):
                byte = 0
                for j in range(4):
                    idx = i + j
                    geno = 0 if idx >= n_samples else genos[idx]
                    if geno < 0 or geno > 3:
                        raise ValueError("genotype code out of range")
                    byte |= (geno & 0b11) << (2 * j)
                out.append(byte)
            if len(out) != bytes_per_snp:
                raise ValueError("unexpected .bed encoding length")
            handle.write(out)


def swap_genotypes(genos: Sequence[int]) -> bytes:
    swapped = bytearray(len(genos))
    for idx, g in enumerate(genos):
        if g == 0:
            swapped[idx] = 3
        elif g == 3:
            swapped[idx] = 0
        else:
            swapped[idx] = g
    return bytes(swapped)
