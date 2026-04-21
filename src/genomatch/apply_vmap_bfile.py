#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

from ._cli_utils import run_cli
from .apply_vmap_utils import build_needed_source_indices, filtered_vmap_rows
from .bfile_utils import (
    HAPLOID_SCHEMA,
    BimRow,
    PackedBedRemapPlan,
    build_packed_bed_remap_plan,
    bytes_per_bed_row,
    count_target_ploidy_genotype_issues,
    decode_bed_chunk,
    missing_bed_row,
    read_bim,
    read_bed_selected_chunks,
    remap_bed_chunk,
    swap_bed_chunk,
    validate_alleles,
    write_bim,
    write_bed_chunks,
)
from .haploid_utils import expected_ploidy_pair, has_non_diploid_ploidy
from .contig_utils import supported_exact_contig_tokens
from .sample_axis_utils import (
    SAMPLE_ID_MODE_CHOICES,
    build_sample_axis_plan,
    compute_reconciliation_missingness_summary,
    parse_fam_table,
    require_identical_sample_signatures,
)
from .vtable_utils import load_metadata, MISSING_SOURCE_SHARD, read_vmap, require_contig_naming, require_rows_match_contig_naming, validate_vmap_metadata

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SourceShard:
    prefix: Path


@dataclass(frozen=True)
class PreparedSourceShard:
    bed_path: Path
    n_snps: int
    n_samples: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply a .vmap to a PLINK bfile.")
    parser.add_argument("--source-prefix", required=True, help="Source PLINK prefix")
    parser.add_argument("--vmap", required=True, help="Input .vmap")
    parser.add_argument("--output-prefix", required=True, help="Output PLINK prefix")
    parser.add_argument("--target-fam", help="Explicit target .fam defining the output sample axis")
    parser.add_argument(
        "--sample-id-mode",
        choices=SAMPLE_ID_MODE_CHOICES,
        default="fid_iid",
        help="Subject-key mode for explicit target-sample reconciliation",
    )
    parser.add_argument(
        "--only-mapped-target",
        action="store_true",
        help="Drop target rows with source_index=-1 instead of writing unmatched target rows with missing genotypes",
    )
    return parser.parse_args()


def vmap_rows_to_bim_rows(vmap_rows) -> List[BimRow]:
    return [BimRow(row.chrom, row.id, "0", row.pos, row.a1, row.a2) for row in vmap_rows]


def bfile_component(prefix: Path, suffix: str) -> Path:
    return Path(str(prefix) + suffix)


def discover_source_shards(source_prefix_arg: str) -> Dict[str, SourceShard]:
    source_prefix = Path(source_prefix_arg)
    name_prefix, name_suffix = source_prefix.name.split("@", 1)
    discovered: Dict[str, SourceShard] = {}
    for token in supported_exact_contig_tokens():
        shard_prefix = source_prefix.parent / f"{name_prefix}{token}{name_suffix}"
        components = {
            suffix: bfile_component(shard_prefix, suffix)
            for suffix in (".bim", ".bed", ".fam")
        }
        present_suffixes = [suffix for suffix, component in components.items() if component.exists()]
        if not present_suffixes:
            continue
        if len(present_suffixes) != len(components):
            missing = next(
                component for suffix, component in components.items() if suffix not in set(present_suffixes)
            )
            raise ValueError(
                f"incomplete source shard for token={token!r}: missing {missing}"
            )
        discovered[token] = SourceShard(
            prefix=shard_prefix,
        )
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
    return [
        (Path(output_prefix_arg.replace("@", chrom_label)), grouped[chrom_label])
        for chrom_label in order
    ]

def resolve_chunk_size() -> int:
    raw = os.environ.get("MATCH_BFILE_APPLY_CHUNK_SIZE", "4096").strip()
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"invalid MATCH_BFILE_APPLY_CHUNK_SIZE: {raw!r}") from exc
    if value <= 0:
        raise ValueError(f"MATCH_BFILE_APPLY_CHUNK_SIZE must be positive; found {raw!r}")
    return value


def chunked_indices(indices: Sequence[int], chunk_size: int) -> Iterable[Sequence[int]]:
    for start in range(0, len(indices), chunk_size):
        yield indices[start : start + chunk_size]


def chunk_count(total_items: int, chunk_size: int) -> int:
    return (total_items + chunk_size - 1) // chunk_size


def has_identity_sample_axis_remap(local_to_output: Sequence[int], output_sample_count: int) -> bool:
    return len(local_to_output) == output_sample_count and all(
        output_idx == local_idx for local_idx, output_idx in enumerate(local_to_output)
    )


def prepare_source_payloads(
    args: argparse.Namespace,
    needed_by_shard: Dict[str, Set[int]],
    source_bim: Path,
    source_bed: Path,
    source_fam: Path,
) -> Tuple[Dict[str, PreparedSourceShard], Dict[str, object]]:
    prepared: Dict[str, PreparedSourceShard] = {}
    source_tables: Dict[str, object] = {}
    if "@" in args.source_prefix:
        discovered_shards = discover_source_shards(args.source_prefix)
        for source_shard in sorted(needed_by_shard):
            shard = shard_for_source_shard(discovered_shards, source_shard)
            shard_bim = bfile_component(shard.prefix, ".bim")
            shard_bed = bfile_component(shard.prefix, ".bed")
            shard_fam = bfile_component(shard.prefix, ".fam")
            shard_rows = read_bim(shard_bim)
            if not shard_rows:
                raise ValueError(f"empty source .bim for source_shard={source_shard!r}")
            validate_alleles(shard_rows, f"source shard {source_shard}")
            shard_table = parse_fam_table(
                shard_fam,
                sample_id_mode=args.sample_id_mode,
                label="source .fam",
            )
            source_tables[source_shard] = shard_table
            prepared[source_shard] = PreparedSourceShard(shard_bed, len(shard_rows), shard_table.sample_count)
        return prepared, source_tables

    if set(needed_by_shard) - {MISSING_SOURCE_SHARD}:
        unexpected = sorted(set(needed_by_shard) - {MISSING_SOURCE_SHARD})
        raise ValueError(
            "apply_vmap_to_bfile.py single-file source prefixes support only source_shard='.' provenance; "
            f"found {unexpected!r}"
        )
    source_rows = read_bim(source_bim)
    if not source_rows:
        raise ValueError("empty source .bim")
    validate_alleles(source_rows, "source")
    source_table = parse_fam_table(
        source_fam,
        sample_id_mode=args.sample_id_mode,
        label="source .fam",
    )
    source_tables[MISSING_SOURCE_SHARD] = source_table
    prepared[MISSING_SOURCE_SHARD] = PreparedSourceShard(source_bed, len(source_rows), source_table.sample_count)
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
            "vmap contains only unmatched target rows; apply_vmap_to_bfile.py will not emit "
            "an all-missing PLINK payload"
        )
    return target_build, vmap_rows


def resolve_target_ploidy_rows(target_rows: Sequence[BimRow], *, target_build: str) -> List[Tuple[int, int]]:
    if not HAPLOID_SCHEMA.exists():
        raise ValueError(f"haploid schema not found: {HAPLOID_SCHEMA}")
    return [
        expected_ploidy_pair(row.chrom, row.bp, genome_build=target_build)
        for row in target_rows
    ]


def load_chunk_source_bed_chunks(
    chunk_indices: Sequence[int],
    vmap_rows,
    prepared_sources: Dict[str, PreparedSourceShard],
) -> Dict[Tuple[str, int], bytes]:
    needed: Dict[str, set[int]] = {}
    for idx in chunk_indices:
        row = vmap_rows[idx]
        if row.source_index == -1:
            continue
        needed.setdefault(row.source_shard, set()).add(row.source_index)
    source_chunks: Dict[Tuple[str, int], bytes] = {}
    for source_shard, selected_indices in sorted(needed.items()):
        prepared = prepared_sources.get(source_shard)
        if prepared is None:
            raise ValueError(f"missing required source shard for source_shard={source_shard!r}")
        missing_indices = sorted(idx for idx in selected_indices if idx >= prepared.n_snps)
        if missing_indices:
            if source_shard == MISSING_SOURCE_SHARD:
                raise ValueError(
                    f"vmap source provenance out of range for source .bim: "
                    f"source_shard={MISSING_SOURCE_SHARD!r}, source_index={missing_indices[0]}"
                )
            raise ValueError(
                f"vmap source provenance out of range for source_shard={source_shard!r}: "
                f"first missing source_index={missing_indices[0]}"
            )
        shard_chunks = read_bed_selected_chunks(prepared.bed_path, prepared.n_samples, prepared.n_snps, selected_indices)
        for source_index, chunk in shard_chunks.items():
            source_chunks[(source_shard, source_index)] = chunk
    return source_chunks


def main() -> int:
    args = parse_args()
    source_prefix = Path(args.source_prefix)
    vmap_path = Path(args.vmap)
    out_prefix = Path(args.output_prefix)
    logger.info("apply_vmap_to_bfile.py: applying %s to %s -> %s", vmap_path, source_prefix, out_prefix)

    source_bim = bfile_component(source_prefix, ".bim")
    source_bed = bfile_component(source_prefix, ".bed")
    source_fam = bfile_component(source_prefix, ".fam")
    if "@" in args.vmap:
        raise ValueError("apply_vmap_to_bfile.py requires a single-file .vmap input")
    if "@" not in args.source_prefix:
        for path, label in ((source_bim, ".bim"), (source_bed, ".bed"), (source_fam, ".fam"), (vmap_path, ".vmap")):
            if not path.exists():
                raise ValueError(f"required input not found ({label}): {path}")
    elif not vmap_path.exists():
        raise ValueError(f"required input not found (.vmap): {vmap_path}")

    target_build, vmap_rows = load_retained_vmap_rows(vmap_path, only_mapped_target=args.only_mapped_target)

    target_rows = vmap_rows_to_bim_rows(vmap_rows)
    validate_alleles(target_rows, "target")
    needed_by_shard = build_needed_source_indices(vmap_rows)
    prepared_sources, source_tables = prepare_source_payloads(
        args,
        needed_by_shard,
        source_bim,
        source_bed,
        source_fam,
    )
    target_table = None
    if args.target_fam:
        target_table = parse_fam_table(
            Path(args.target_fam),
            sample_id_mode=args.sample_id_mode,
            label="target .fam",
        )
        output_fam_source = Path(args.target_fam)
    else:
        output_fam_source = require_identical_sample_signatures(source_tables, descriptor=".fam")
    sample_axis_plan = build_sample_axis_plan(
        source_tables,
        output_sample_path=output_fam_source,
        explicit_target_table=target_table,
    )
    n_samples = sample_axis_plan.output_sample_count
    missing_count = sum(1 for row in vmap_rows if row.source_index == -1)
    packed_missing_row = missing_bed_row(n_samples)
    output_bytes_per_snp = bytes_per_bed_row(n_samples)
    chunk_size = resolve_chunk_size()
    identity_sample_axis = not sample_axis_plan.reconciliation_active
    identity_remap_by_shard: Dict[str, bool] = {}
    packed_sample_axis_plans: Dict[str, PackedBedRemapPlan] = {}
    if not identity_sample_axis:
        identity_remap_by_shard = {
            shard: has_identity_sample_axis_remap(local_to_output, n_samples)
            for shard, local_to_output in sample_axis_plan.source_local_to_output.items()
        }
        packed_sample_axis_plans = {
            shard: build_packed_bed_remap_plan(local_to_output, n_samples)
            for shard, local_to_output in sample_axis_plan.source_local_to_output.items()
            if not identity_remap_by_shard[shard]
        }

    ploidy_rows = resolve_target_ploidy_rows(
        target_rows,
        target_build=target_build,
    )
    haploid_het_incompatible_count = 0
    absent_nonmissing_incompatible_count = 0
    unknown_sex_unvalidated_count = 0

    def row_bed_chunk(
        idx: int,
        source_chunks: Dict[Tuple[str, int], bytes],
        swapped_packed_cache: Dict[Tuple[str, int], bytes],
        remapped_packed_cache: Dict[Tuple[str, int, bool], bytes],
    ) -> bytes:
        nonlocal haploid_het_incompatible_count, absent_nonmissing_incompatible_count, unknown_sex_unvalidated_count
        row = vmap_rows[idx]
        ploidy_pair = ploidy_rows[idx]
        needs_ploidy_validation = has_non_diploid_ploidy(ploidy_pair)
        needs_swap = row.allele_op in {"swap", "flip_swap"}
        if row.source_index == -1:
            packed_chunk = packed_missing_row
        else:
            key = (row.source_shard, row.source_index)
            packed_chunk = source_chunks.get(key)
            if packed_chunk is None:
                raise ValueError(
                    f"missing source genotypes for requested SNP: "
                    f"source_shard={row.source_shard!r}, source_index={row.source_index}"
                )
            if needs_swap:
                swapped_packed = swapped_packed_cache.get(key)
                if swapped_packed is None:
                    swapped_packed = swap_bed_chunk(packed_chunk)
                    swapped_packed_cache[key] = swapped_packed
                packed_chunk = swapped_packed
            if not identity_sample_axis and not identity_remap_by_shard[row.source_shard]:
                remapped_key = (row.source_shard, row.source_index, needs_swap)
                remapped_chunk = remapped_packed_cache.get(remapped_key)
                if remapped_chunk is None:
                    remapped_chunk = remap_bed_chunk(packed_chunk, packed_sample_axis_plans[row.source_shard])
                    remapped_packed_cache[remapped_key] = remapped_chunk
                packed_chunk = remapped_chunk
        if not needs_ploidy_validation:
            return packed_chunk

        genos = decode_bed_chunk(packed_chunk, n_samples)
        haploid_issues, absent_issues, unknown_sex_issues = count_target_ploidy_genotype_issues(
            genos,
            sample_axis_plan.output_sexes,
            ploidy_pair[0],
            ploidy_pair[1],
            target_rows[idx],
        )
        haploid_het_incompatible_count += haploid_issues
        absent_nonmissing_incompatible_count += absent_issues
        unknown_sex_unvalidated_count += unknown_sex_issues
        return packed_chunk

    def iter_output_rows(indices: Sequence[int], output_bed_path: Path) -> Iterable[bytes]:
        total_chunks = chunk_count(len(indices), chunk_size)
        for done_chunks, chunk in enumerate(chunked_indices(indices, chunk_size), start=1):
            source_chunks = load_chunk_source_bed_chunks(chunk, vmap_rows, prepared_sources)
            swapped_packed_cache: Dict[Tuple[str, int], bytes] = {}
            remapped_packed_cache: Dict[Tuple[str, int, bool], bytes] = {}
            for idx in chunk:
                yield row_bed_chunk(
                    idx,
                    source_chunks,
                    swapped_packed_cache,
                    remapped_packed_cache,
                )
            logger.info(
                "apply_vmap_to_bfile.py progress: %s/%s chunks -> %s",
                done_chunks,
                total_chunks,
                output_bed_path,
            )

    for shard_out_prefix, indices in grouped_output_indices(vmap_rows, args.output_prefix):
        out_bim = bfile_component(shard_out_prefix, ".bim")
        out_bed = bfile_component(shard_out_prefix, ".bed")
        out_fam = bfile_component(shard_out_prefix, ".fam")
        out_ploidy = bfile_component(shard_out_prefix, ".ploidy")
        shard_target_rows = [target_rows[idx] for idx in indices]
        write_bim(out_bim, shard_target_rows)
        shutil.copyfile(sample_axis_plan.output_sample_path, out_fam)
        shard_ploidy_rows = [ploidy_rows[idx] for idx in indices]
        if any(has_non_diploid_ploidy(ploidy_pair) for ploidy_pair in shard_ploidy_rows):
            with open(out_ploidy, "w", encoding="utf-8", newline="\n") as handle:
                for male, female in shard_ploidy_rows:
                    handle.write(f"{male}\t{female}\n")
        write_bed_chunks(out_bed, iter_output_rows(indices, out_bed), bytes_per_snp=output_bytes_per_snp)
    if haploid_het_incompatible_count or absent_nonmissing_incompatible_count:
        parts: List[str] = []
        if haploid_het_incompatible_count:
            parts.append(
                f"found {haploid_het_incompatible_count} incompatible heterozygous haploid-target genotypes"
            )
        if absent_nonmissing_incompatible_count:
            parts.append(
                f"found {absent_nonmissing_incompatible_count} incompatible absent-region nonmissing genotypes"
            )
        logger.warning("%s.", "; ".join(parts))
    if unknown_sex_unvalidated_count:
        logger.warning(
            "skipped ploidy validation for %s sex-dependent row/sample cells with unknown output sex.",
            unknown_sex_unvalidated_count,
        )
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
        logger.warning("%s variants missing from source; filled with missing genotypes.", missing_count)
    logger.info("apply_vmap_to_bfile.py: completed projection for %s retained target rows", len(vmap_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli(main))
