#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import numpy as np
import pgenlib

from ._cli_utils import run_cli
from .apply_vmap_utils import build_needed_source_indices, filtered_vmap_rows, output_variant_id
from .contig_utils import supported_exact_contig_tokens
from .haploid_utils import expected_ploidy_pair, is_sex_dependent_ploidy
from .sample_axis_utils import (
    SAMPLE_ID_MODE_CHOICES,
    SAMPLE_ID_MODE_FID_IID,
    build_sample_axis_plan,
    compute_reconciliation_missingness_summary,
    parse_psam_table,
    require_identical_sample_signatures,
    require_psam_fid_presence_consistent,
)
from .vtable_utils import (
    MISSING_SOURCE_SHARD,
    VariantRow,
    load_metadata,
    read_vmap,
    require_contig_naming,
    require_rows_match_contig_naming,
    validate_vmap_metadata,
)


SUPPORTED_CHANNEL_HARDCALL = "hardcall"
SUPPORTED_CHANNEL_PHASE = "phase"
SUPPORTED_CHANNEL_DOSAGE = "dosage"
logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SourceShard:
    prefix: Path


@dataclass(frozen=True)
class SourcePVarRow:
    chrom: str
    pos: str
    variant_id: str
    ref: str
    alt: str
    allele_ct: int


@dataclass(frozen=True)
class PreparedSourceShard:
    pgen_path: Path
    pvar_rows: Sequence[SourcePVarRow]
    hardcall_phase_present: bool
    n_samples: int


@dataclass(frozen=True)
class SourceRowPlan:
    channel: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply a .vmap to a PLINK 2 PFILE.")
    parser.add_argument("--source-prefix", required=True, help="Source PLINK 2 prefix")
    parser.add_argument("--vmap", required=True, help="Input .vmap")
    parser.add_argument("--output-prefix", required=True, help="Output PLINK 2 prefix")
    parser.add_argument("--target-psam", help="Explicit target .psam defining the output sample axis")
    parser.add_argument(
        "--sample-id-mode",
        choices=SAMPLE_ID_MODE_CHOICES,
        default="fid_iid",
        help="Subject-key mode for explicit target-sample reconciliation",
    )
    parser.add_argument(
        "--only-mapped-target",
        action="store_true",
        help="Drop target rows with source_index=-1 instead of writing unmatched target rows with missing payload data",
    )
    parser.add_argument(
        "--retain-snp-id",
        action="store_true",
        help="Use retained target-side .vmap id values as output IDs instead of generated chrom:pos:a1:a2 IDs",
    )
    return parser.parse_args()


def pfile_component(prefix: Path, suffix: str) -> Path:
    return Path(str(prefix) + suffix)


def discover_source_shards(source_prefix_arg: str) -> Dict[str, SourceShard]:
    source_prefix = Path(source_prefix_arg)
    name_prefix, name_suffix = source_prefix.name.split("@", 1)
    discovered: Dict[str, SourceShard] = {}
    for token in supported_exact_contig_tokens():
        shard_prefix = source_prefix.parent / f"{name_prefix}{token}{name_suffix}"
        components = {suffix: pfile_component(shard_prefix, suffix) for suffix in (".pgen", ".pvar", ".psam")}
        present_suffixes = [suffix for suffix, component in components.items() if component.exists()]
        if not present_suffixes:
            continue
        if len(present_suffixes) != len(components):
            missing = next(component for suffix, component in components.items() if suffix not in set(present_suffixes))
            raise ValueError(f"incomplete source shard for token={token!r}: missing {missing}")
        discovered[token] = SourceShard(prefix=shard_prefix)
    if not discovered:
        raise ValueError(f"no source shards found for template: {source_prefix}")
    return discovered


def shard_for_source_shard(discovered: Dict[str, SourceShard], source_shard: str) -> SourceShard:
    shard = discovered.get(source_shard)
    if shard is None:
        raise ValueError(f"missing required source shard for source_shard={source_shard!r}")
    return shard


def grouped_output_indices(vmap_rows, output_prefix_arg: str) -> List[Tuple[Path, List[int]]]:
    if "@" not in output_prefix_arg:
        return [(Path(output_prefix_arg), list(range(len(vmap_rows))))]
    grouped: Dict[str, List[int]] = {}
    order: List[str] = []
    for idx, row in enumerate(vmap_rows):
        if row.chrom not in grouped:
            grouped[row.chrom] = []
            order.append(row.chrom)
        grouped[row.chrom].append(idx)
    return [(Path(output_prefix_arg.replace("@", chrom_label)), grouped[chrom_label]) for chrom_label in order]


def resolve_chunk_size() -> int:
    raw = os.environ.get("MATCH_PFILE_APPLY_CHUNK_SIZE", "1024").strip()
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"invalid MATCH_PFILE_APPLY_CHUNK_SIZE: {raw!r}") from exc
    if value <= 0:
        raise ValueError(f"MATCH_PFILE_APPLY_CHUNK_SIZE must be positive; found {raw!r}")
    return value


def chunked_indices(indices: Sequence[int], chunk_size: int) -> Iterable[Sequence[int]]:
    for start in range(0, len(indices), chunk_size):
        yield indices[start : start + chunk_size]


def parse_pvar_rows(path: Path) -> List[SourcePVarRow]:
    rows: List[SourcePVarRow] = []
    with open(path, "r", encoding="utf-8", newline="\n") as handle:
        header = None
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
        if header is None:
            raise ValueError(".pvar header not found")
        columns = {name: idx for idx, name in enumerate(header)}
        required = ["CHROM", "POS", "ID", "REF", "ALT"]
        missing = [name for name in required if name not in columns]
        if missing:
            raise ValueError(f".pvar is missing required columns: {', '.join(missing)}")
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            max_required_idx = max(columns[name] for name in required)
            if len(parts) <= max_required_idx:
                raise ValueError(f"malformed .pvar row in {path}")
            alt = parts[columns["ALT"]]
            allele_ct = 1 + len(alt.split(",")) if alt else 1
            rows.append(
                SourcePVarRow(
                    chrom=parts[columns["CHROM"]],
                    pos=parts[columns["POS"]],
                    variant_id=parts[columns["ID"]] or ".",
                    ref=parts[columns["REF"]],
                    alt=alt,
                    allele_ct=allele_ct,
                )
            )
    if not rows:
        raise ValueError(f"empty source .pvar: {path}")
    return rows


def prepare_source_payloads(args: argparse.Namespace, needed_by_shard: Dict[str, Set[int]]):
    prepared: Dict[str, PreparedSourceShard] = {}
    source_tables: Dict[str, object] = {}
    if "@" in args.source_prefix:
        discovered_shards = discover_source_shards(args.source_prefix)
        for source_shard in sorted(needed_by_shard):
            shard = shard_for_source_shard(discovered_shards, source_shard)
            shard_pgen = pfile_component(shard.prefix, ".pgen")
            shard_pvar = pfile_component(shard.prefix, ".pvar")
            shard_psam = pfile_component(shard.prefix, ".psam")
            pvar_rows = parse_pvar_rows(shard_pvar)
            shard_table = parse_psam_table(
                shard_psam,
                sample_id_mode=args.sample_id_mode,
                label="source .psam",
            )
            source_tables[source_shard] = shard_table
            reader = pgenlib.PgenReader(os.fsencode(str(shard_pgen)), shard_table.sample_count)
            try:
                phase_present = reader.hardcall_phase_present()
            finally:
                reader.close()
            prepared[source_shard] = PreparedSourceShard(
                pgen_path=shard_pgen,
                pvar_rows=pvar_rows,
                hardcall_phase_present=phase_present,
                n_samples=shard_table.sample_count,
            )
        return prepared, source_tables

    if set(needed_by_shard) - {MISSING_SOURCE_SHARD}:
        unexpected = sorted(set(needed_by_shard) - {MISSING_SOURCE_SHARD})
        raise ValueError(
            "apply_vmap_to_pfile.py single-file source prefixes support only source_shard='.' provenance; "
            f"found {unexpected!r}"
        )

    source_prefix = Path(args.source_prefix)
    source_pgen = pfile_component(source_prefix, ".pgen")
    source_pvar = pfile_component(source_prefix, ".pvar")
    source_psam = pfile_component(source_prefix, ".psam")
    pvar_rows = parse_pvar_rows(source_pvar)
    source_table = parse_psam_table(
        source_psam,
        sample_id_mode=args.sample_id_mode,
        label="source .psam",
    )
    source_tables[MISSING_SOURCE_SHARD] = source_table
    reader = pgenlib.PgenReader(os.fsencode(str(source_pgen)), source_table.sample_count)
    try:
        phase_present = reader.hardcall_phase_present()
    finally:
        reader.close()
    prepared[MISSING_SOURCE_SHARD] = PreparedSourceShard(
        pgen_path=source_pgen,
        pvar_rows=pvar_rows,
        hardcall_phase_present=phase_present,
        n_samples=source_table.sample_count,
    )
    return prepared, source_tables


def load_retained_vmap_rows(vmap_path: Path, *, only_mapped_target: bool):
    metadata = load_metadata(vmap_path)
    validate_vmap_metadata(metadata)
    target_build = str(metadata["target"]["genome_build"])
    target_contig_naming = require_contig_naming(dict(metadata["target"]), label="variant map target")
    all_vmap_rows = read_vmap(vmap_path)
    if not all_vmap_rows:
        raise ValueError("empty vmap")
    require_rows_match_contig_naming(all_vmap_rows, target_contig_naming, label="variant map target")
    vmap_rows = filtered_vmap_rows(all_vmap_rows, only_mapped_target=only_mapped_target)
    if not vmap_rows:
        raise ValueError("no retained target rows remain after applying --only-mapped-target")
    if all(row.source_index == -1 for row in vmap_rows):
        raise ValueError(
            "vmap contains only unmatched target rows; apply_vmap_to_pfile.py will not emit "
            "an all-missing PLINK 2 payload"
        )
    return target_build, vmap_rows


def allele_array_for_row(reader: pgenlib.PgenReader, source_index: int, n_samples: int) -> np.ndarray:
    alleles = np.empty(n_samples * 2, dtype=np.int32)
    reader.read_alleles(source_index, alleles)
    return alleles


def allele_and_phase_for_row(reader: pgenlib.PgenReader, source_index: int, n_samples: int) -> Tuple[np.ndarray, np.ndarray]:
    alleles = np.empty(n_samples * 2, dtype=np.int32)
    phasepresent = np.zeros(n_samples, dtype=np.uint8)
    reader.read_alleles_and_phasepresent(source_index, alleles, phasepresent)
    return alleles, phasepresent


def dosages_for_row(reader: pgenlib.PgenReader, source_index: int, n_samples: int) -> np.ndarray:
    dosages = np.empty(n_samples, dtype=np.float32)
    reader.read_dosages(source_index, dosages)
    return dosages


def hardcall_dosages(alleles: np.ndarray) -> np.ndarray:
    n_samples = len(alleles) // 2
    out = np.empty(n_samples, dtype=np.float32)
    for sample_idx in range(n_samples):
        a1 = int(alleles[2 * sample_idx])
        a2 = int(alleles[2 * sample_idx + 1])
        if a1 == -9 and a2 == -9:
            out[sample_idx] = -9.0
            continue
        if a1 == -9 or a2 == -9:
            if a1 == -9:
                out[sample_idx] = float(a2)
            else:
                out[sample_idx] = float(a1)
            continue
        out[sample_idx] = float(a1 + a2)
    return out


def alleles_have_calls(alleles: np.ndarray) -> bool:
    return any(int(value) != -9 for value in alleles)


def dosages_require_preservation(dosages: np.ndarray, hardcall_dos: np.ndarray) -> bool:
    for dosage, hardcall_dosage in zip(dosages, hardcall_dos):
        dosage_value = float(dosage)
        if dosage_value == -9.0 and hardcall_dosage == -9.0:
            continue
        if dosage_value == -9.0 or hardcall_dosage == -9.0:
            return True
        if abs(dosage_value - hardcall_dosage) > 1e-6:
            return True
    return False


def dosage_only_representation(dosages: np.ndarray, hardcall_dos: np.ndarray) -> bool:
    for dosage, hardcall_dosage in zip(dosages, hardcall_dos):
        dosage_value = float(dosage)
        if dosage_value == -9.0:
            continue
        if hardcall_dosage == -9.0:
            return True
        if abs(dosage_value - round(dosage_value)) > 1e-6:
            return True
    return False


def plan_source_rows(vmap_rows, prepared_sources: Dict[str, PreparedSourceShard]):
    readers: Dict[str, pgenlib.PgenReader] = {}
    plans: Dict[Tuple[str, int], SourceRowPlan] = {}
    channels_seen: Set[str] = set()
    try:
        for row in vmap_rows:
            if row.source_index == -1:
                continue
            key = (row.source_shard, row.source_index)
            if key in plans:
                continue
            prepared = prepared_sources.get(row.source_shard)
            if prepared is None:
                raise ValueError(f"missing required source shard for source_shard={row.source_shard!r}")
            if row.source_index >= len(prepared.pvar_rows):
                if row.source_shard == MISSING_SOURCE_SHARD:
                    raise ValueError(
                        f"vmap source provenance out of range for source .pvar: "
                        f"source_shard={MISSING_SOURCE_SHARD!r}, source_index={row.source_index}"
                    )
                raise ValueError(
                    f"vmap source provenance out of range for source_shard={row.source_shard!r}: "
                    f"first missing source_index={row.source_index}"
                )
            if prepared.pvar_rows[row.source_index].allele_ct != 2:
                raise ValueError(
                    "retained mapped source row is non-biallelic in source .pvar: "
                    f"source_shard={row.source_shard!r}, source_index={row.source_index}"
                )
            reader = readers.get(row.source_shard)
            if reader is None:
                reader = pgenlib.PgenReader(os.fsencode(str(prepared.pgen_path)), prepared.n_samples)
                readers[row.source_shard] = reader

            if prepared.hardcall_phase_present:
                alleles, phasepresent = allele_and_phase_for_row(reader, row.source_index, prepared.n_samples)
            else:
                alleles = allele_array_for_row(reader, row.source_index, prepared.n_samples)
                phasepresent = np.zeros(prepared.n_samples, dtype=np.uint8)
            dosages = dosages_for_row(reader, row.source_index, prepared.n_samples)
            hardcall_dos = hardcall_dosages(alleles)
            has_calls = alleles_have_calls(alleles)
            dosage_needed = dosages_require_preservation(dosages, hardcall_dos)
            if has_calls and dosage_needed and not dosage_only_representation(dosages, hardcall_dos):
                raise ValueError(
                    "retained mapped source row cannot be preserved through the supported pgenlib path: "
                    f"source_shard={row.source_shard!r}, source_index={row.source_index}"
                )
            if dosage_needed:
                channel = SUPPORTED_CHANNEL_DOSAGE
            elif prepared.hardcall_phase_present and any(int(flag) != 0 for flag in phasepresent):
                channel = SUPPORTED_CHANNEL_PHASE
            else:
                channel = SUPPORTED_CHANNEL_HARDCALL
            plans[key] = SourceRowPlan(channel=channel)
            channels_seen.add(channel)
    finally:
        for reader in readers.values():
            reader.close()

    if len(channels_seen) > 1:
        rendered = ", ".join(sorted(channels_seen))
        logger.warning(
            "retained mapped source rows have incoherent supported PFILE channel availability; "
            "preserving rows with channels: %s.",
            rendered,
        )
    return plans


def swap_alleles(alleles: np.ndarray) -> np.ndarray:
    swapped = alleles.copy()
    for idx, value in enumerate(swapped):
        if value == 0:
            swapped[idx] = 1
        elif value == 1:
            swapped[idx] = 0
    return swapped


def swap_dosages(dosages: np.ndarray) -> np.ndarray:
    swapped = dosages.copy()
    mask = swapped != -9.0
    swapped[mask] = 2.0 - swapped[mask]
    return swapped


def write_output_pvar(path: Path, rows: Sequence[VariantRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for row in rows:
            handle.write(f"{row.chrom}\t{row.pos}\t{row.id}\t{row.a2}\t{row.a1}\n")


def target_variant_rows(vmap_rows, *, retain_snp_id: bool) -> List[VariantRow]:
    return [
        VariantRow(row.chrom, row.pos, output_variant_id(row, retain_snp_id=retain_snp_id), row.a1, row.a2)
        for row in vmap_rows
    ]


def resolve_target_ploidy_rows(rows: Sequence[VariantRow], *, target_build: str) -> List[Tuple[int, int]]:
    return [expected_ploidy_pair(row.chrom, row.pos, genome_build=target_build) for row in rows]


def classify_hardcall_allele_pair(a1: int, a2: int) -> str:
    if a1 == -9 and a2 == -9:
        return "fully_missing"
    if (a1 in {0, 1} and a2 == -9) or (a2 in {0, 1} and a1 == -9):
        return "haploid_like"
    if a1 == 0 and a2 == 0:
        return "diploid_hom_ref"
    if (a1, a2) in {(0, 1), (1, 0)}:
        return "diploid_het"
    if a1 == 1 and a2 == 1:
        return "diploid_hom_alt"
    return "malformed_partial"


def resolved_sample_ploidy(sex: int, ploidy_pair: Tuple[int, int]) -> int | None:
    if sex == 0 and is_sex_dependent_ploidy(ploidy_pair):
        return None
    if sex == 1:
        return ploidy_pair[0]
    return ploidy_pair[1]


def validate_hardcall_ploidy(
    alleles: np.ndarray,
    sexes: Sequence[int],
    ploidy_pair: Tuple[int, int],
) -> Tuple[int, int, int, int]:
    diploid_het_in_haploid = 0
    malformed_partial = 0
    nonmissing_hardcall_in_absent = 0
    unknown_sex_unvalidated = 0
    for sample_idx in range(len(alleles) // 2):
        ploidy = resolved_sample_ploidy(sexes[sample_idx], ploidy_pair)
        if ploidy is None:
            unknown_sex_unvalidated += 1
            continue
        state = classify_hardcall_allele_pair(
            int(alleles[2 * sample_idx]),
            int(alleles[2 * sample_idx + 1]),
        )
        if ploidy == 2:
            if state == "malformed_partial":
                malformed_partial += 1
            continue
        if ploidy == 1:
            if state == "diploid_het":
                diploid_het_in_haploid += 1
            elif state == "malformed_partial":
                malformed_partial += 1
            continue
        if ploidy == 0:
            if state != "fully_missing":
                nonmissing_hardcall_in_absent += 1
            continue
        raise ValueError(f"unsupported ploidy value {ploidy}")
    return (
        diploid_het_in_haploid,
        malformed_partial,
        nonmissing_hardcall_in_absent,
        unknown_sex_unvalidated,
    )


def validate_dosage_ploidy(
    dosages: np.ndarray,
    sexes: Sequence[int],
    ploidy_pair: Tuple[int, int],
) -> Tuple[int, int]:
    nonmissing_dosage_in_absent = 0
    unknown_sex_unvalidated = 0
    for sample_idx, dosage in enumerate(dosages):
        ploidy = resolved_sample_ploidy(sexes[sample_idx], ploidy_pair)
        if ploidy is None:
            unknown_sex_unvalidated += 1
            continue
        if ploidy == 0 and float(dosage) != -9.0:
            nonmissing_dosage_in_absent += 1
    return nonmissing_dosage_in_absent, unknown_sex_unvalidated


def scatter_alleles(local_alleles: np.ndarray, local_to_output: Sequence[int], output_sample_count: int) -> np.ndarray:
    out = np.full(output_sample_count * 2, -9, dtype=np.int32)
    for local_idx, output_idx in enumerate(local_to_output):
        if output_idx == -1:
            continue
        out[2 * output_idx] = int(local_alleles[2 * local_idx])
        out[2 * output_idx + 1] = int(local_alleles[2 * local_idx + 1])
    return out


def scatter_phasepresent(local_phase: np.ndarray, local_to_output: Sequence[int], output_sample_count: int) -> np.ndarray:
    out = np.zeros(output_sample_count, dtype=np.uint8)
    for local_idx, output_idx in enumerate(local_to_output):
        if output_idx != -1:
            out[output_idx] = int(local_phase[local_idx])
    return out


def scatter_dosages(local_dosages: np.ndarray, local_to_output: Sequence[int], output_sample_count: int) -> np.ndarray:
    out = np.full(output_sample_count, -9.0, dtype=np.float32)
    for local_idx, output_idx in enumerate(local_to_output):
        if output_idx != -1:
            out[output_idx] = float(local_dosages[local_idx])
    return out


def main() -> int:
    args = parse_args()
    vmap_path = Path(args.vmap)
    source_prefix = Path(args.source_prefix)
    output_prefix = Path(args.output_prefix)
    logger.info("apply_vmap_to_pfile.py: applying %s to %s -> %s", vmap_path, source_prefix, output_prefix)
    source_pgen = pfile_component(source_prefix, ".pgen")
    source_pvar = pfile_component(source_prefix, ".pvar")
    source_psam = pfile_component(source_prefix, ".psam")

    if "@" in args.vmap:
        raise ValueError("apply_vmap_to_pfile.py requires a single-file .vmap input")
    if "@" not in args.source_prefix:
        for path, label in ((source_pgen, ".pgen"), (source_pvar, ".pvar"), (source_psam, ".psam"), (vmap_path, ".vmap")):
            if not path.exists():
                raise ValueError(f"required input not found ({label}): {path}")
    elif not vmap_path.exists():
        raise ValueError(f"required input not found (.vmap): {vmap_path}")

    target_build, vmap_rows = load_retained_vmap_rows(vmap_path, only_mapped_target=args.only_mapped_target)
    needed_by_shard = build_needed_source_indices(vmap_rows)
    prepared_sources, source_tables = prepare_source_payloads(args, needed_by_shard)
    if args.sample_id_mode == SAMPLE_ID_MODE_FID_IID:
        require_psam_fid_presence_consistent(source_tables.values(), label="referenced source .psam files")
    target_table = None
    if args.target_psam:
        target_table = parse_psam_table(
            Path(args.target_psam),
            sample_id_mode=args.sample_id_mode,
            label="target .psam",
        )
        if args.sample_id_mode == SAMPLE_ID_MODE_FID_IID:
            require_psam_fid_presence_consistent(
                [*source_tables.values(), target_table],
                label="referenced source .psam files and --target-psam",
            )
        output_psam_source = Path(args.target_psam)
    else:
        output_psam_source = require_identical_sample_signatures(source_tables, descriptor=".psam")
    sample_axis_plan = build_sample_axis_plan(
        source_tables,
        output_sample_path=output_psam_source,
        explicit_target_table=target_table,
    )
    n_samples = sample_axis_plan.output_sample_count
    source_row_plans = plan_source_rows(vmap_rows, prepared_sources)

    missing_count = sum(1 for row in vmap_rows if row.source_index == -1)
    hardcall_phase_present = any(plan.channel == SUPPORTED_CHANNEL_PHASE for plan in source_row_plans.values())
    dosage_present = any(plan.channel == SUPPORTED_CHANNEL_DOSAGE for plan in source_row_plans.values())
    chunk_size = resolve_chunk_size()
    target_rows = target_variant_rows(vmap_rows, retain_snp_id=args.retain_snp_id)
    ploidy_rows = resolve_target_ploidy_rows(target_rows, target_build=target_build)
    diploid_het_in_haploid_count = 0
    malformed_partial_count = 0
    nonmissing_hardcall_in_absent_count = 0
    nonmissing_dosage_in_absent_count = 0
    unknown_sex_unvalidated_count = 0

    readers: Dict[str, pgenlib.PgenReader] = {}
    try:
        for shard_out_prefix, indices in grouped_output_indices(vmap_rows, args.output_prefix):
            out_pgen = pfile_component(shard_out_prefix, ".pgen")
            out_pvar = pfile_component(shard_out_prefix, ".pvar")
            out_psam = pfile_component(shard_out_prefix, ".psam")
            shard_target_rows = [target_rows[idx] for idx in indices]
            write_output_pvar(out_pvar, shard_target_rows)
            shutil.copyfile(sample_axis_plan.output_sample_path, out_psam)
            writer = pgenlib.PgenWriter(
                os.fsencode(str(out_pgen)),
                n_samples,
                variant_ct=len(indices),
                nonref_flags=False,
                allele_ct_limit=2,
                hardcall_phase_present=hardcall_phase_present,
                dosage_present=dosage_present,
            )
            try:
                for chunk in chunked_indices(indices, chunk_size):
                    alleles_cache: Dict[Tuple[str, int], np.ndarray] = {}
                    phase_cache: Dict[Tuple[str, int], np.ndarray] = {}
                    dosages_cache: Dict[Tuple[str, int], np.ndarray] = {}
                    for idx in chunk:
                        row = vmap_rows[idx]
                        if row.source_index == -1:
                            if dosage_present:
                                writer.append_dosages(np.full(n_samples, -9.0, dtype=np.float32))
                            else:
                                writer.append_alleles(np.full(n_samples * 2, -9, dtype=np.int32))
                            continue

                        key = (row.source_shard, row.source_index)
                        prepared = prepared_sources[row.source_shard]
                        reader = readers.get(row.source_shard)
                        if reader is None:
                            reader = pgenlib.PgenReader(os.fsencode(str(prepared.pgen_path)), prepared.n_samples)
                            readers[row.source_shard] = reader
                        plan = source_row_plans[key]
                        local_to_output = sample_axis_plan.source_local_to_output[row.source_shard]

                        if plan.channel == SUPPORTED_CHANNEL_DOSAGE:
                            dosages = dosages_cache.get(key)
                            if dosages is None:
                                dosages = dosages_for_row(reader, row.source_index, prepared.n_samples)
                                dosages_cache[key] = dosages
                            out_dosages = dosages if row.allele_op in {"identity", "flip"} else swap_dosages(dosages)
                            out_dosages = scatter_dosages(out_dosages, local_to_output, n_samples)
                            dosage_issues, unknown_sex_issues = validate_dosage_ploidy(
                                out_dosages,
                                sample_axis_plan.output_sexes,
                                ploidy_rows[idx],
                            )
                            nonmissing_dosage_in_absent_count += dosage_issues
                            unknown_sex_unvalidated_count += unknown_sex_issues
                            writer.append_dosages(out_dosages)
                            continue

                        alleles = alleles_cache.get(key)
                        if alleles is None:
                            alleles = allele_array_for_row(reader, row.source_index, prepared.n_samples)
                            alleles_cache[key] = alleles
                        out_alleles = alleles if row.allele_op in {"identity", "flip"} else swap_alleles(alleles)
                        out_alleles = scatter_alleles(out_alleles, local_to_output, n_samples)
                        hardcall_issues = validate_hardcall_ploidy(
                            out_alleles,
                            sample_axis_plan.output_sexes,
                            ploidy_rows[idx],
                        )
                        diploid_het_in_haploid_count += hardcall_issues[0]
                        malformed_partial_count += hardcall_issues[1]
                        nonmissing_hardcall_in_absent_count += hardcall_issues[2]
                        unknown_sex_unvalidated_count += hardcall_issues[3]

                        if plan.channel == SUPPORTED_CHANNEL_PHASE:
                            phasepresent = phase_cache.get(key)
                            if phasepresent is None:
                                _, phasepresent = allele_and_phase_for_row(reader, row.source_index, prepared.n_samples)
                                phase_cache[key] = phasepresent
                            writer.append_partially_phased(
                                out_alleles,
                                scatter_phasepresent(phasepresent, local_to_output, n_samples),
                            )
                        else:
                            writer.append_alleles(out_alleles)
            finally:
                writer.close()
        reconciliation_summary = compute_reconciliation_missingness_summary(sample_axis_plan, vmap_rows)
        if reconciliation_summary is not None:
            logger.info(
                "Sample-axis reconciliation summary: introduced "
                f"{reconciliation_summary.total_missing_cells} missing row/sample cells across "
                f"{reconciliation_summary.mapped_variant_count} retained mapped variants and "
                f"{reconciliation_summary.output_sample_count} output subjects."
            )
            if reconciliation_summary.subjects_over_threshold:
                logger.warning(
                    "sample-axis reconciliation caused >50%% added missingness for %s output subjects.",
                    reconciliation_summary.subjects_over_threshold,
                )
            if reconciliation_summary.variants_over_threshold:
                logger.warning(
                    "sample-axis reconciliation caused >50%% added missingness for %s retained mapped variants.",
                    reconciliation_summary.variants_over_threshold,
                )
        if missing_count:
            logger.warning("%s variants missing from source; filled with missing payload rows.", missing_count)
        ploidy_warning_parts: List[str] = []
        if diploid_het_in_haploid_count:
            ploidy_warning_parts.append(
                f"found {diploid_het_in_haploid_count} incompatible diploid-het hardcalls in haploid target ploidy"
            )
        if malformed_partial_count:
            ploidy_warning_parts.append(
                f"found {malformed_partial_count} incompatible malformed partial hardcalls under target ploidy validation"
            )
        if nonmissing_hardcall_in_absent_count:
            ploidy_warning_parts.append(
                f"found {nonmissing_hardcall_in_absent_count} incompatible nonmissing hardcalls in absent target regions"
            )
        if nonmissing_dosage_in_absent_count:
            ploidy_warning_parts.append(
                f"found {nonmissing_dosage_in_absent_count} incompatible nonmissing dosages in absent target regions"
            )
        if ploidy_warning_parts:
            logger.warning("%s.", "; ".join(ploidy_warning_parts))
        if unknown_sex_unvalidated_count:
            logger.warning(
                "skipped ploidy validation for %s sex-dependent row/sample cells with unknown output sex.",
                unknown_sex_unvalidated_count,
            )
        logger.info("apply_vmap_to_pfile.py: completed projection for %s retained target rows", len(vmap_rows))
        return 0
    finally:
        for reader in readers.values():
            reader.close()


if __name__ == "__main__":
    raise SystemExit(run_cli(main))
