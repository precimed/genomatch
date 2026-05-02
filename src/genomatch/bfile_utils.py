#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from .haploid_utils import HAPLOID_SCHEMA_PATH as HAPLOID_SCHEMA
from .haploid_utils import HAPLOID_CHROMS, expected_ploidy_pair, is_sex_dependent_ploidy


ALLOWED_BASES = {"A", "C", "G", "T"}
BED_HEADER = b"\x6c\x1b\x01"


@dataclass(frozen=True)
class BimRow:
    chrom: str
    snp: str
    cm: str
    bp: str
    a1: str
    a2: str


@dataclass(frozen=True)
class PackedBedRemapPlan:
    output_sample_count: int
    output_bytes_per_snp: int
    base_missing_chunk: np.ndarray
    output_byte_indices: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    source_byte_indices: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    source_shifts: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


@dataclass(frozen=True)
class PackedPloidyValidationPlan:
    unknown_sex_unvalidated: int
    haploid_output_byte_indices: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    absent_output_byte_indices: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


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


def _count_genotype_code_by_slot(code: int) -> np.ndarray:
    out = np.zeros((4, 256), dtype=np.uint8)
    for slot in range(4):
        for value in range(256):
            geno = (value >> (2 * slot)) & 0b11
            out[slot, value] = 1 if geno == code else 0
    return out


_HET_BY_SLOT = _count_genotype_code_by_slot(2)
_MISSING_BY_SLOT = _count_genotype_code_by_slot(1)
_NONMISSING_BY_SLOT = np.uint8(1) - _MISSING_BY_SLOT


def build_packed_ploidy_validation_plan(
    sexes: Sequence[int],
    ploidy_male: int,
    ploidy_female: int,
) -> PackedPloidyValidationPlan:
    unknown_sex_unvalidated = 0
    ploidy_pair = (ploidy_male, ploidy_female)
    sex_dependent = is_sex_dependent_ploidy(ploidy_pair)

    haploid_per_slot: list[list[int]] = [[], [], [], []]
    absent_per_slot: list[list[int]] = [[], [], [], []]
    for sample_idx, sex in enumerate(sexes):
        if sex == 0 and sex_dependent:
            unknown_sex_unvalidated += 1
            continue
        ploidy = ploidy_male if sex == 1 else ploidy_female
        if ploidy == 2:
            continue
        slot = sample_idx % 4
        byte_idx = sample_idx // 4
        if ploidy == 1:
            haploid_per_slot[slot].append(byte_idx)
        elif ploidy == 0:
            absent_per_slot[slot].append(byte_idx)
        else:
            raise ValueError(f"unsupported ploidy value {ploidy}")

    return PackedPloidyValidationPlan(
        unknown_sex_unvalidated=unknown_sex_unvalidated,
        haploid_output_byte_indices=tuple(np.asarray(indices, dtype=np.intp) for indices in haploid_per_slot),  # type: ignore[arg-type]
        absent_output_byte_indices=tuple(np.asarray(indices, dtype=np.intp) for indices in absent_per_slot),  # type: ignore[arg-type]
    )


def count_target_ploidy_genotype_issues_packed(
    chunk: bytes,
    plan: PackedPloidyValidationPlan,
) -> Tuple[int, int, int]:
    packed = np.frombuffer(chunk, dtype=np.uint8)
    haploid_het_incompatible = 0
    absent_nonmissing_incompatible = 0

    for slot in range(4):
        haploid_indices = plan.haploid_output_byte_indices[slot]
        if haploid_indices.size:
            haploid_het_incompatible += int(_HET_BY_SLOT[slot, packed[haploid_indices]].sum(dtype=np.int64))
        absent_indices = plan.absent_output_byte_indices[slot]
        if absent_indices.size:
            absent_nonmissing_incompatible += int(_NONMISSING_BY_SLOT[slot, packed[absent_indices]].sum(dtype=np.int64))

    return haploid_het_incompatible, absent_nonmissing_incompatible, plan.unknown_sex_unvalidated


def bytes_per_bed_row(n_samples: int) -> int:
    return (n_samples + 3) // 4


def _validate_bed_header_and_size(handle, *, n_samples: int, n_snps: int) -> int:
    bytes_per_snp = bytes_per_bed_row(n_samples)
    header = handle.read(len(BED_HEADER))
    if header != BED_HEADER:
        raise ValueError("unsupported .bed format (expected SNP-major)")
    handle.seek(0, 2)
    file_size = handle.tell()
    expected_size = len(BED_HEADER) + bytes_per_snp * n_snps
    if file_size < expected_size:
        raise ValueError(".bed file is shorter than expected")
    return bytes_per_snp


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


def _swap_packed_bed_byte(value: int) -> int:
    swapped = 0
    for shift in range(0, 8, 2):
        geno = (value >> shift) & 0b11
        if geno == 0:
            new_geno = 3
        elif geno == 3:
            new_geno = 0
        else:
            new_geno = geno
        swapped |= new_geno << shift
    return swapped


SWAP_BED_CHUNK_TABLE = bytes(_swap_packed_bed_byte(value) for value in range(256))


def swap_bed_chunk(chunk: bytes) -> bytes:
    return chunk.translate(SWAP_BED_CHUNK_TABLE)


def encode_bed_row(genos: Sequence[int], n_samples: int) -> bytes:
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
    if len(out) != bytes_per_bed_row(n_samples):
        raise ValueError("unexpected .bed encoding length")
    return bytes(out)


def missing_bed_row(n_samples: int) -> bytes:
    full_bytes, remainder = divmod(n_samples, 4)
    out = bytearray([0b01010101]) * full_bytes
    if remainder:
        byte = 0
        for shift in range(0, remainder * 2, 2):
            byte |= 0b01 << shift
        out.append(byte)
    return bytes(out)


def build_packed_bed_remap_plan(local_to_output: Sequence[int], output_sample_count: int) -> PackedBedRemapPlan:
    output_to_local = [-1] * output_sample_count
    for local_idx, output_idx in enumerate(local_to_output):
        if output_idx == -1:
            continue
        if output_idx < 0 or output_idx >= output_sample_count:
            raise ValueError(
                f"sample-axis mapping output index out of range: {output_idx} for output_sample_count={output_sample_count}"
            )
        if output_to_local[output_idx] != -1:
            raise ValueError(f"duplicate output sample index in sample-axis mapping: {output_idx}")
        output_to_local[output_idx] = local_idx

    slot_output_byte_indices: list[np.ndarray] = []
    slot_source_byte_indices: list[np.ndarray] = []
    slot_source_shifts: list[np.ndarray] = []
    for slot in range(4):
        output_byte_indices: list[int] = []
        source_byte_indices: list[int] = []
        source_shifts: list[int] = []
        for output_idx, local_idx in enumerate(output_to_local):
            if output_idx % 4 != slot or local_idx == -1:
                continue
            output_byte_indices.append(output_idx // 4)
            source_byte_indices.append(local_idx // 4)
            source_shifts.append((local_idx % 4) * 2)
        slot_output_byte_indices.append(np.asarray(output_byte_indices, dtype=np.intp))
        slot_source_byte_indices.append(np.asarray(source_byte_indices, dtype=np.intp))
        slot_source_shifts.append(np.asarray(source_shifts, dtype=np.uint8))

    return PackedBedRemapPlan(
        output_sample_count=output_sample_count,
        output_bytes_per_snp=bytes_per_bed_row(output_sample_count),
        base_missing_chunk=np.frombuffer(missing_bed_row(output_sample_count), dtype=np.uint8).copy(),
        output_byte_indices=tuple(slot_output_byte_indices),  # type: ignore[arg-type]
        source_byte_indices=tuple(slot_source_byte_indices),  # type: ignore[arg-type]
        source_shifts=tuple(slot_source_shifts),  # type: ignore[arg-type]
    )


def remap_bed_chunk(chunk: bytes, plan: PackedBedRemapPlan) -> bytes:
    source = np.frombuffer(chunk, dtype=np.uint8)
    out = plan.base_missing_chunk.copy()
    for slot in range(4):
        output_byte_indices = plan.output_byte_indices[slot]
        if output_byte_indices.size == 0:
            continue
        source_values = source[plan.source_byte_indices[slot]]
        source_values = (source_values >> plan.source_shifts[slot]) & np.uint8(0b11)
        output_shift = np.uint8(slot * 2)
        output_mask = np.uint8(0b11 << (slot * 2))
        current = out[output_byte_indices]
        out[output_byte_indices] = (current & np.uint8(~int(output_mask) & 0xFF)) | (source_values << output_shift)
    return out.tobytes()


def read_bed_selected_chunks(
    bed_path: Path, n_samples: int, n_snps: int, selected_indices: Iterable[int]
) -> Dict[int, bytes]:
    indices = sorted(set(selected_indices))
    if not indices:
        return {}
    with open(bed_path, "rb") as handle:
        bytes_per_snp = _validate_bed_header_and_size(handle, n_samples=n_samples, n_snps=n_snps)
        data: Dict[int, bytes] = {}
        for idx in indices:
            offset = len(BED_HEADER) + idx * bytes_per_snp
            handle.seek(offset)
            chunk = handle.read(bytes_per_snp)
            if len(chunk) != bytes_per_snp:
                raise ValueError(".bed file is shorter than expected")
            data[idx] = chunk
    return data


def read_bed_selected(
    bed_path: Path, n_samples: int, n_snps: int, selected_indices: Iterable[int]
) -> Dict[int, bytes]:
    indices = sorted(set(selected_indices))
    if not indices:
        return {}
    bytes_per_snp = (n_samples + 3) // 4
    with open(bed_path, "rb") as handle:
        bytes_per_snp = _validate_bed_header_and_size(handle, n_samples=n_samples, n_snps=n_snps)
        data: Dict[int, bytes] = {}
        for idx in indices:
            offset = len(BED_HEADER) + idx * bytes_per_snp
            handle.seek(offset)
            chunk = handle.read(bytes_per_snp)
            if len(chunk) != bytes_per_snp:
                raise ValueError(".bed file is shorter than expected")
            data[idx] = decode_bed_chunk(chunk, n_samples)
    return data


def write_bed_chunks(path: Path, chunks: Iterable[bytes], *, bytes_per_snp: int) -> None:
    with open(path, "wb") as handle:
        handle.write(BED_HEADER)
        for chunk in chunks:
            if len(chunk) != bytes_per_snp:
                raise ValueError("unexpected .bed chunk length")
            handle.write(chunk)


def write_bed_matrix(path: Path, matrix: Iterable[Sequence[int]], n_samples: int) -> None:
    write_bed_chunks(
        path,
        (encode_bed_row(genos, n_samples) for genos in matrix),
        bytes_per_snp=bytes_per_bed_row(n_samples),
    )


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
