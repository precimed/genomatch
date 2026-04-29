#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ._cli_utils import run_cli
from .contig_utils import SUPPORTED_CONTIG_NAMINGS, contig_label_for_naming, normalize_chrom_label
from .importer_utils import DiscoveredInputShard
from .vmap_prepare_variants import all_retained_stage_outputs
from .vtable_utils import CANONICAL_CONTIG_RANK, load_metadata, metadata_path_for, write_metadata
from .workflow_wrapper_utils import (
    delete_variant_object,
    planned_existing_variant_outputs,
    read_target_metadata,
    run_command,
    tool_command,
    variant_object_exists,
    variant_object_path,
)

WRAPPER_NAME = "prepare_variants_sharded.py"
INPUT_FORMATS = ("bim", "pvar", "vcf")

logger = logging.getLogger(__name__)


def _discovery_tokens() -> list[str]:
    autosomes = [str(idx) for idx in range(1, 23)]
    chr_autosomes = [f"chr{idx}" for idx in range(1, 23)]
    extras = [
        "X",
        "chrX",
        "XY",
        "chrXY",
        "PAR",
        "par",
        "NONPAR",
        "nonpar",
        "23",
        "25",
        "Y",
        "chrY",
        "24",
        "M",
        "MT",
        "chrM",
        "chrMT",
        "26",
    ]
    return autosomes + chr_autosomes + extras


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare chromosome-sharded raw input variants by running prepare_variants.py per target-contig group, "
            "then concatenating per-group final .vmap outputs in declared contig order."
        )
    )
    parser.add_argument("--input", required=True, help="Sharded raw importer input template containing @")
    parser.add_argument("--input-format", required=True, choices=INPUT_FORMATS, help="Importer selection")
    parser.add_argument("--output", required=True, help="Non-sharded output stem; final artifact is <output>.vmap")
    parser.add_argument("--prefix", required=True, help="Sharded retained-intermediate prefix containing @")
    parser.add_argument("--dst-build", default="GRCh38", help="Destination genome build (default: GRCh38)")
    parser.add_argument(
        "--dst-contig-naming",
        default="ncbi",
        choices=sorted(SUPPORTED_CONTIG_NAMINGS),
        help="Destination contig naming for prepared artifacts (default: ncbi)",
    )
    parser.add_argument(
        "--allow-strand-flips",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Pass through to prepare_variants.py (default: enabled)",
    )
    parser.add_argument(
        "--norm-indels",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Pass through to prepare_variants.py (default: enabled)",
    )
    parser.add_argument("--drop-strand-ambiguous", action="store_true", help="Pass through to prepare_variants.py")
    parser.add_argument("--max-allele-length", type=int, default=150, help="Pass through to prepare_variants.py")
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--resume", action="store_true", help="Resume by skipping completed per-group finals")
    mode_group.add_argument("--force", action="store_true", help="Delete wrapper-managed outputs first, then rerun")
    return parser.parse_args()


def discover_input_shards(input_template: Path) -> list[DiscoveredInputShard]:
    discovered: list[DiscoveredInputShard] = []
    for token in _discovery_tokens():
        candidate = Path(str(input_template).replace("@", token))
        if candidate.is_file():
            discovered.append(DiscoveredInputShard(path=candidate, source_shard=token))
    if not discovered:
        raise ValueError(f"no input shards found for template: {input_template}")
    return sorted(discovered, key=lambda shard: str(shard.path))


def group_shards_by_canonical(shards: list[DiscoveredInputShard]) -> dict[str, list[DiscoveredInputShard]]:
    grouped: dict[str, list[DiscoveredInputShard]] = {}
    for shard in shards:
        canonical = normalize_chrom_label(shard.source_shard)
        if canonical not in CANONICAL_CONTIG_RANK:
            raise ValueError(f"{WRAPPER_NAME}: unsupported shard token for grouping: {shard.source_shard!r}")
        grouped.setdefault(canonical, []).append(shard)
    return grouped


def ordered_group_labels(grouped: dict[str, list[DiscoveredInputShard]]) -> list[str]:
    return sorted(grouped, key=lambda canonical: CANONICAL_CONTIG_RANK[canonical])


def group_prefix(prefix_template: Path, canonical_group: str) -> Path:
    return Path(str(prefix_template).replace("@", canonical_group))


def prepare_command(
    args: argparse.Namespace,
    *,
    input_template: Path,
    output_prefix: Path,
    canonical_group: str,
    group_shards: list[DiscoveredInputShard],
) -> list[str]:
    target_contig = contig_label_for_naming(canonical_group, args.dst_contig_naming)
    shard_tokens = ",".join(shard.source_shard for shard in group_shards)
    cmd = tool_command(
        "prepare_variants.py",
        "--input",
        str(input_template),
        "--input-format",
        args.input_format,
        "--shards",
        shard_tokens,
        "--contigs",
        target_contig,
        "--prefix",
        str(output_prefix),
        "--output",
        str(output_prefix),
        "--dst-build",
        args.dst_build,
        "--dst-contig-naming",
        args.dst_contig_naming,
        "--max-allele-length",
        str(args.max_allele_length),
    )
    cmd.append("--allow-strand-flips" if args.allow_strand_flips else "--no-allow-strand-flips")
    cmd.append("--norm-indels" if args.norm_indels else "--no-norm-indels")
    if args.drop_strand_ambiguous:
        cmd.append("--drop-strand-ambiguous")
    if args.resume:
        cmd.append("--resume")
    if args.force:
        cmd.append("--force")
    return cmd


def per_group_managed_outputs(prefix_template: Path, group_labels: list[str]) -> list[Path]:
    outputs: list[Path] = []
    for label in group_labels:
        prefix = group_prefix(prefix_template, label)
        outputs.extend(all_retained_stage_outputs(prefix))
        outputs.append(variant_object_path(prefix))
    return outputs


def require_outputs_absent(paths: list[Path]) -> None:
    existing = planned_existing_variant_outputs(paths)
    if not existing:
        return
    rendered = "\n".join(f"- {path}" for path in sorted(existing))
    raise ValueError(
        "planned output or intermediate files already exist; rerun with --resume to reuse them "
        "or --force to delete them first and rerun cleanly:\n"
        f"{rendered}"
    )


def delete_outputs(paths: list[Path]) -> None:
    for path in paths:
        delete_variant_object(path)


def validate_group_metadata(final_paths: list[Path], *, input_template_raw: str) -> dict:
    if not final_paths:
        raise ValueError("no prepared group outputs to merge")
    for path in final_paths:
        if not variant_object_exists(path):
            raise ValueError(f"missing required per-group final .vmap: {path}")

    first_metadata = load_metadata(final_paths[0])
    first_target = read_target_metadata(final_paths[0])
    first_derived_from = first_metadata.get("derived_from")
    if first_derived_from != input_template_raw:
        raise ValueError("per-group final .vmap metadata derived_from must equal the original --input template")

    for path in final_paths[1:]:
        target = read_target_metadata(path)
        metadata = load_metadata(path)
        if target.get("genome_build") != first_target.get("genome_build"):
            raise ValueError("per-group final .vmap files have mismatched genome_build")
        if target.get("contig_naming") != first_target.get("contig_naming"):
            raise ValueError("per-group final .vmap files have mismatched contig_naming")
        if metadata.get("derived_from") != first_derived_from:
            raise ValueError("per-group final .vmap files have mismatched derived_from")

    return first_metadata


def write_concatenated_vmap(path: Path, final_paths: list[Path], metadata: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_output = path.with_name(path.stem + ".tmp.vmap")
    temp_meta = metadata_path_for(temp_output)
    final_meta = metadata_path_for(path)
    temp_qc = temp_output.with_name(temp_output.name + ".qc.tsv")
    final_qc = path.with_name(path.name + ".qc.tsv")
    if temp_output.exists():
        temp_output.unlink()
    if temp_meta.exists():
        temp_meta.unlink()
    if temp_qc.exists():
        temp_qc.unlink()
    try:
        with open(temp_output, "w", encoding="utf-8", newline="\n") as output:
            for final_path in final_paths:
                with open(final_path, "r", encoding="utf-8", newline="") as handle:
                    for line in handle:
                        output.write(line)
        write_metadata(temp_output, metadata)
        temp_output.replace(path)
        temp_meta.replace(final_meta)
        if temp_qc.exists():
            temp_qc.unlink()
    except Exception:
        if temp_output.exists():
            temp_output.unlink()
        if temp_meta.exists():
            temp_meta.unlink()
        if temp_qc.exists():
            temp_qc.unlink()
        raise
    if final_qc.exists():
        final_qc.unlink()


def main() -> int:
    args = parse_args()
    if args.dst_contig_naming == "plink_splitx":
        raise ValueError("prepare_variants_sharded.py does not support --dst-contig-naming=plink_splitx")
    if "@" not in args.input:
        raise ValueError("prepare_variants_sharded.py requires --input to contain '@'")
    if "@" not in args.prefix:
        raise ValueError("prepare_variants_sharded.py requires --prefix to contain '@'")
    if "@" in args.output:
        raise ValueError("prepare_variants_sharded.py requires non-sharded --output without '@'")

    input_template = Path(args.input)
    prefix_template = Path(args.prefix)
    final_output = variant_object_path(Path(args.output))

    shards = discover_input_shards(input_template)
    grouped = group_shards_by_canonical(shards)
    group_labels = ordered_group_labels(grouped)

    managed_outputs = [*per_group_managed_outputs(prefix_template, group_labels), final_output]
    if args.force:
        delete_outputs(managed_outputs)
    elif not args.resume:
        require_outputs_absent(managed_outputs)

    if args.resume and variant_object_exists(final_output):
        logger.info("%s: wrote %s", WRAPPER_NAME, final_output)
        return 0

    final_paths: list[Path] = []
    for group_label in group_labels:
        output_prefix = group_prefix(prefix_template, group_label)
        group_final = variant_object_path(output_prefix)
        final_paths.append(group_final)
        if args.resume and variant_object_exists(group_final):
            logger.info("%s: skipping group %s; output exists at %s", WRAPPER_NAME, group_label, group_final)
            continue
        run_command(
            prepare_command(
                args,
                input_template=input_template,
                output_prefix=output_prefix,
                canonical_group=group_label,
                group_shards=grouped[group_label],
            )
        )

    metadata = validate_group_metadata(final_paths, input_template_raw=args.input)
    write_concatenated_vmap(final_output, final_paths, metadata)
    logger.info("%s: wrote %s", WRAPPER_NAME, final_output)
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
