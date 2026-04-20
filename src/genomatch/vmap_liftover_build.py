#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from .reference_utils import fetch_reference_bases, resolve_bcftools_binary, resolve_liftover_assets
from .haploid_utils import expected_ploidy_pair
from .vtable_utils import (
    VMapRow,
    VariantRow,
    compose_allele_ops,
    convert_contig_label,
    declared_coordinate_sort_key,
    duplicate_target_row_keys,
    load_variant_object,
    normalize_chrom_label,
    require_contig_naming,
    require_rows_match_contig_naming,
    shard_local_provenance,
    normalize_allele_token,
    target_row_key,
    validate_allele_value,
    variant_rows_from_vmap_rows,
    write_metadata,
    write_vmap_status_qc,
    write_vmap,
    write_vtable,
)


@dataclass(frozen=True)
class LiftoverSourceRecord:
    row: VariantRow
    source_shard: str
    source_index: int
    upstream_allele_op: str
    input_allele_op: str


@dataclass(frozen=True)
class PreparedLiftoverRow:
    row_id: str
    source_record: LiftoverSourceRecord
    status: str
    ref: str | None = None
    alt: str | None = None


@dataclass(frozen=True)
class LiftoverPreparation:
    status: str
    ref: str | None = None
    alt: str | None = None
    input_allele_op: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Liftover a .vtable/.vmap across genome builds, preserving the input object type. "
            "Input contigs may be UCSC or numeric/NCBI-style, but the internal reference path is UCSC-only."
        )
    )
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--output", required=True, help="Output .vtable or .vmap")
    parser.add_argument("--target-build", required=True, help="Target genome build")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse retained bcftools +liftover output if it already exists, then rerun parsing/QC/output steps",
    )
    return parser.parse_args()


def convert_rows_to_ucsc(rows: Sequence[VariantRow], contig_naming: str) -> List[VariantRow]:
    out_rows: List[VariantRow] = []
    for row in rows:
        try:
            chrom = convert_contig_label(row.chrom, contig_naming, "ucsc")
        except ValueError as exc:
            raise ValueError(
                f"unable to normalize input contig {row.chrom!r} from declared contig_naming={contig_naming!r} "
                "to internal UCSC naming for liftover"
            ) from exc
        out_rows.append(VariantRow(chrom, row.pos, row.id, row.a1, row.a2))
    return out_rows


def parse_source_position(row: VariantRow) -> int:
    try:
        return int(row.pos)
    except ValueError as exc:
        raise ValueError(f"invalid source position for liftover input: {row.chrom}:{row.pos}") from exc


def resolve_ref_alt(row: VariantRow, ref_base: str) -> LiftoverPreparation:
    validate_allele_value(row.a1, label="liftover_build.py input")
    validate_allele_value(row.a2, label="liftover_build.py input")
    if ref_base not in {"A", "C", "G", "T"}:
        raise ValueError(f"source reference base not found for liftover input: {row.chrom}:{row.pos}")
    a1 = row.a1
    a2 = row.a2
    is_snv = len(a1) == 1 and len(a2) == 1
    is_supported_non_snv = (len(a1) == 1) != (len(a2) == 1)
    if a2 == ref_base and a1 != ref_base:
        return LiftoverPreparation("ready", ref_base, a1, "swap")
    if a1 == ref_base and a2 != ref_base:
        return LiftoverPreparation("ready", ref_base, a2, "identity")
    if is_snv:
        raise ValueError(
            f"source alleles do not match reference for liftover input: {row.chrom}:{row.pos} "
            f"({row.a1}/{row.a2}, ref={ref_base})"
        )
    if is_supported_non_snv:
        return LiftoverPreparation("unsupported_non_snv")
    return LiftoverPreparation("unsupported_non_snv")


def prepare_liftover_rows(
    rows: Sequence[VariantRow],
    source_provenance: Sequence[Tuple[str, int]],
    upstream_allele_ops: Sequence[str],
    source_fasta: Path,
) -> List[PreparedLiftoverRow]:
    positions: List[int] = []
    queries: List[Tuple[str, int]] = []
    for row in rows:
        pos = parse_source_position(row)
        positions.append(pos)
        queries.append((row.chrom, pos))
    reference_bases = fetch_reference_bases(source_fasta, queries)

    prepared_rows: List[PreparedLiftoverRow] = []
    for idx, (row, (source_shard, source_index), upstream_allele_op) in enumerate(
        zip(rows, source_provenance, upstream_allele_ops)
    ):
        row_id = f"row{idx}"
        ref_base = reference_bases.get((row.chrom, positions[idx]), "")
        preparation = resolve_ref_alt(row, ref_base)
        source_record = LiftoverSourceRecord(
            row,
            source_shard,
            source_index,
            upstream_allele_op,
            preparation.input_allele_op or "identity",
        )
        prepared_rows.append(
            PreparedLiftoverRow(
                row_id=row_id,
                source_record=source_record,
                status=preparation.status,
                ref=preparation.ref,
                alt=preparation.alt,
            )
        )
    return prepared_rows


def write_temp_vcf(
    path: Path,
    prepared_rows: Sequence[PreparedLiftoverRow],
) -> Dict[str, LiftoverSourceRecord]:
    row_lookup: Dict[str, LiftoverSourceRecord] = {}
    ready_rows = [prepared for prepared in prepared_rows if prepared.status == "ready"]
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##INFO=<ID=SRC_IDX,Number=1,Type=Integer,Description="Source row index">\n')
        for chrom in sorted({prepared.source_record.row.chrom for prepared in ready_rows}, key=chrom_sort_key):
            handle.write(f"##contig=<ID={chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for idx, prepared in enumerate(prepared_rows):
            if prepared.status != "ready":
                continue
            row = prepared.source_record.row
            row_lookup[prepared.row_id] = prepared.source_record
            handle.write(
                "\t".join(
                    [
                        row.chrom,
                        row.pos,
                        prepared.row_id,
                        prepared.ref or "",
                        prepared.alt or "",
                        ".",
                        "PASS",
                        f"SRC_IDX={idx}",
                    ]
                )
                + "\n"
            )
    return row_lookup


def run_bcftools_liftover(
    input_vcf: Path,
    output_vcf: Path,
    source_fasta: Path,
    target_fasta: Path,
    chain_path: Path,
) -> None:
    bcftools = resolve_bcftools_binary()
    cmd = [
        bcftools,
        "+liftover",
        str(input_vcf),
        "-o",
        str(output_vcf),
        "--",
        "-s",
        str(source_fasta),
        "-f",
        str(target_fasta),
        "-c",
        str(chain_path),
        "--flip-tag",
        "FLIP",
        "--swap-tag",
        "SWAP",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        text = (result.stderr or result.stdout).strip()
        raise ValueError(f"bcftools +liftover failed: {text}")
def chrom_sort_key(label: str) -> Tuple[int, str]:
    canonical = normalize_chrom_label(label)
    if canonical.isdigit():
        return (int(canonical), "")
    order = {"X": 23, "Y": 24, "MT": 25}
    if canonical in order:
        return (order[canonical], "")
    return (999, label)


def parse_lifted_vcf(
    path: Path,
    row_lookup: Dict[str, LiftoverSourceRecord],
    final_contig_naming: str,
    source_build: str,
    target_build: str,
) -> Tuple[List[VMapRow], Dict[str, str]]:
    parsed_rows: List[Tuple[str, VMapRow]] = []
    status_by_row_id: Dict[str, str] = {}

    vcf_data: List[List[str]] = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 8:
                raise ValueError(f"invalid lifted VCF row: {raw_line.strip()}")
            vcf_data.append(parts[:8])

    acgt_set = {"A", "C", "G", "T"}
    for parts in vcf_data:
        chrom, pos, row_id, ref, alt, _qual, filt, info_raw = parts
        if row_id not in row_lookup:
            raise ValueError(f"unexpected liftover row ID: {row_id}")
        if filt not in {".", "PASS"}:
            continue
        record = row_lookup[row_id]

        if "," in alt:
            status_by_row_id[row_id] = "unsupported_non_snv"
            continue

        ref_token = normalize_allele_token(ref)
        alt_token = normalize_allele_token(alt)
        if (
            not ref_token
            or not alt_token
            or not set(ref_token) <= acgt_set
            or not set(alt_token) <= acgt_set
        ):
            status_by_row_id[row_id] = "unsupported_non_snv"
            continue

        if len(ref_token) > 1 and len(alt_token) > 1:
            status_by_row_id[row_id] = "unsupported_non_snv"
            continue

        try:
            final_chrom = convert_contig_label(chrom, "ucsc", final_contig_naming)
        except ValueError:
            status_by_row_id[row_id] = "unsupported_target_contig"
            continue

        source_ploidy_pair = expected_ploidy_pair(
            record.row.chrom,
            record.row.pos,
            genome_build=source_build,
        )
        target_ploidy_pair = expected_ploidy_pair(
            final_chrom,
            pos,
            genome_build=target_build,
        )
        if target_ploidy_pair != source_ploidy_pair:
            status_by_row_id[row_id] = "ploidy_class_changed"
            continue

        liftover_op = allele_op_from_liftover_info(info_raw)
        # Compose the full allele-orientation history:
        # final_op = upstream ∘ input_to_vcf ∘ bcftools_liftover ∘ vcf_to_canonical_output
        #
        # vcf_to_canonical_output is always "swap" by construction because lifted VCF rows
        # are parsed as REF/ALT, while the emitted canonical row is written as a1=ALT, a2=REF.
        #
        # input_to_vcf is the source-side normalization required to express the input row as
        # a valid source-build VCF REF/ALT record. When liftover_build.py is used after the
        # recommended restrict_build_compatible.py step, this is expected to be "swap"
        # because that tool canonicalizes retained rows to a1=non-reference, a2=reference.
        output_allele_op = "swap"
        parsed_rows.append(
            (
                row_id,
                VMapRow(
                    final_chrom,
                    pos,
                    record.row.id,
                    alt_token,
                    ref_token,
                    record.source_shard,
                    record.source_index,
                    compose_allele_ops(
                        record.upstream_allele_op,
                        compose_allele_ops(
                            record.input_allele_op,
                            compose_allele_ops(liftover_op, output_allele_op),
                        ),
                    ),
                ),
            )
        )
        status_by_row_id[row_id] = "lifted"

    parsed_rows = sorted(
        parsed_rows,
        key=lambda item: declared_coordinate_sort_key(item[1], final_contig_naming, label="liftover_build.py output"),
    )
    duplicate_keys = duplicate_target_row_keys([row for _row_id, row in parsed_rows])

    out_rows: List[VMapRow] = []
    for row_id, row in parsed_rows:
        if target_row_key(row) in duplicate_keys:
            status_by_row_id[row_id] = "duplicate_target"
            continue
        out_rows.append(row)
    return out_rows, status_by_row_id


def allele_op_from_liftover_info(info_raw: str) -> str:
    tags = info_raw.split(";") if info_raw and info_raw != "." else []
    has_flip = "FLIP" in tags
    swap_value: str | None = None
    for tag in tags:
        if tag.startswith("SWAP="):
            swap_value = tag.split("=", 1)[1]
            break
    if swap_value not in {None, "1"}:
        raise ValueError(
            f"liftover output contains unsupported SWAP annotation {swap_value!r}; "
            "cannot represent this variant in canonical v1 allele_op"
        )
    if has_flip and swap_value == "1":
        return "flip_swap"
    if has_flip:
        return "flip"
    if swap_value == "1":
        return "swap"
    return "identity"


def liftover_debug_vcf_path(output_path: Path, *, kind: str) -> Path:
    return output_path.with_name(output_path.name + f".bcftools_{kind}.vcf")


def print_resume_skip(output_vcf: Path) -> None:
    print(
        f"liftover_build.py: skipping bcftools +liftover; retained output exists at {output_vcf}",
        file=sys.stderr,
    )


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.exists():
        raise ValueError(f"input not found: {input_path}")

    loaded = load_variant_object(input_path)
    source_rows = loaded.target_rows
    source_meta = dict(loaded.target_metadata)
    source_build = str(source_meta.get("genome_build"))
    source_contig_naming = require_contig_naming(source_meta, label="variant object")
    if source_contig_naming == "plink_splitx":
        raise ValueError(
            "liftover_build.py does not support contig_naming=plink_splitx; "
            "normalize to build-independent ncbi/ucsc/plink first, liftover, then normalize back if needed"
        )
    target_build = str(args.target_build)
    if source_build == "unknown":
        raise ValueError(
            "source metadata genome_build must be known for liftover "
            f"(declared genome_build={source_build!r}, declared contig_naming={source_contig_naming!r}, "
            "internal_reference_naming='ucsc')"
        )
    if source_build == target_build:
        raise ValueError("source and target builds are identical; liftover is not required")
    require_rows_match_contig_naming(source_rows, source_contig_naming, label="variant object")

    source_fasta, target_fasta, chain_path = resolve_liftover_assets(source_build, target_build)
    ucsc_rows = convert_rows_to_ucsc(source_rows, source_contig_naming)
    if loaded.base_vmap_rows is not None:
        source_provenance = [(row.source_shard, row.source_index) for row in loaded.base_vmap_rows]
        upstream_allele_ops = [row.allele_op for row in loaded.base_vmap_rows]
    else:
        source_provenance = shard_local_provenance(["."] * len(source_rows))
        upstream_allele_ops = ["identity"] * len(source_rows)

    input_vcf = liftover_debug_vcf_path(output_path, kind="input")
    output_vcf = liftover_debug_vcf_path(output_path, kind="output")
    input_vcf.parent.mkdir(parents=True, exist_ok=True)
    prepared_rows = prepare_liftover_rows(
        ucsc_rows,
        source_provenance,
        upstream_allele_ops,
        source_fasta,
    )
    row_lookup = write_temp_vcf(input_vcf, prepared_rows)
    if args.resume and output_vcf.exists():
        print_resume_skip(output_vcf)
    else:
        run_bcftools_liftover(input_vcf, output_vcf, source_fasta, target_fasta, chain_path)
    raw_rows, parse_status_by_row_id = parse_lifted_vcf(
        output_vcf,
        row_lookup,
        source_contig_naming,
        source_build,
        target_build,
    )

    qc_rows: List[Tuple[str, int, str, str]] = []
    for prepared in prepared_rows:
        row = prepared.source_record.row
        source_shard = prepared.source_record.source_shard
        source_index = prepared.source_record.source_index
        if prepared.status != "ready":
            status = prepared.status
        else:
            status = parse_status_by_row_id.get(prepared.row_id, "unmapped")
        qc_rows.append((source_shard, source_index, row.id, status))

    if loaded.base_vmap_rows is not None:
        out_rows = raw_rows
        out_meta = dict(loaded.raw_metadata)
        out_meta["target"] = dict(out_meta["target"])
        out_meta["target"]["genome_build"] = target_build
        write_vmap(output_path, out_rows)
    else:
        out_rows = variant_rows_from_vmap_rows(raw_rows)
        out_meta = dict(loaded.raw_metadata)
        out_meta["genome_build"] = target_build
        write_vtable(output_path, out_rows)
    write_metadata(output_path, out_meta)
    write_vmap_status_qc(output_path.with_name(output_path.name + ".qc.tsv"), qc_rows)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
