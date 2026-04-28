#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
import tempfile
from pathlib import Path

from ._cli_utils import run_cli
from .importer_utils import DiscoveredInputShard, resolve_import_input_paths
from .vmap_prepare_variants import all_retained_stage_outputs
from .contig_utils import SUPPORTED_CONTIG_NAMINGS
from .vtable_utils import load_metadata, write_metadata
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
KIND_LABEL_BY_FORMAT = {
    "bim": ".bim",
    "pvar": ".pvar",
    "vcf": "VCF",
}

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare sharded raw input variants by running prepare_variants.py per discovered shard, "
            "then concatenating and sorting the per-shard final .vmap outputs."
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
    parser.add_argument("--chr2use", "--contigs", dest="chr2use", help="Pass through to prepare_variants.py")
    parser.add_argument("--max-allele-length", type=int, default=150, help="Pass through to prepare_variants.py")
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--resume", action="store_true", help="Resume by skipping completed per-shard finals")
    mode_group.add_argument("--force", action="store_true", help="Delete wrapper-managed outputs first, then rerun")
    return parser.parse_args()


def shard_prefix(prefix_template: Path, source_shard: str) -> Path:
    return Path(str(prefix_template).replace("@", source_shard))


def sort_scratch_prefix(prefix_template: Path) -> Path:
    return Path(str(prefix_template).replace("@", "all_targets") + ".sort_tmp")


def one_shard_template(input_template: Path, shard: DiscoveredInputShard, workspace: Path) -> Path:
    name_prefix, name_suffix = input_template.name.split("@", 1)
    template = workspace / input_template.name
    link_path = workspace / f"{name_prefix}{shard.source_shard}{name_suffix}"
    link_path.symlink_to(shard.path.resolve())
    return template


def prepare_command(args: argparse.Namespace, *, input_template: Path, shard_output_prefix: Path) -> list[str]:
    cmd = tool_command(
        "prepare_variants.py",
        "--input",
        str(input_template),
        "--input-format",
        args.input_format,
        "--prefix",
        str(shard_output_prefix),
        "--output",
        str(shard_output_prefix),
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
    if args.chr2use:
        cmd.extend(["--chr2use", args.chr2use])
    if args.resume:
        cmd.append("--resume")
    if args.force:
        cmd.append("--force")
    return cmd


def per_shard_managed_outputs(prefix_template: Path, shards: list[DiscoveredInputShard]) -> list[Path]:
    outputs: list[Path] = []
    for shard in shards:
        prefix = shard_prefix(prefix_template, shard.source_shard)
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


def validate_shard_metadata(final_paths: list[Path]) -> dict:
    if not final_paths:
        raise ValueError("no prepared shard outputs to merge")
    first_metadata = load_metadata(final_paths[0])
    first_target = read_target_metadata(final_paths[0])
    for path in final_paths[1:]:
        target = read_target_metadata(path)
        if target.get("genome_build") != first_target.get("genome_build"):
            raise ValueError("per-shard final .vmap files have mismatched genome_build")
        if target.get("contig_naming") != first_target.get("contig_naming"):
            raise ValueError("per-shard final .vmap files have mismatched contig_naming")
    return first_metadata


def write_concatenated_vmap(path: Path, final_paths: list[Path], metadata: dict, derived_from: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as output:
        for shard_path in final_paths:
            with open(shard_path, "r", encoding="utf-8", newline="") as handle:
                for line in handle:
                    output.write(line)
    out_metadata = json.loads(json.dumps(metadata))
    out_metadata["derived_from"] = derived_from
    write_metadata(path, out_metadata)


def sort_final_vmap(args: argparse.Namespace, *, concatenated_path: Path, final_output: Path) -> None:
    run_command(
        tool_command(
            "sort_variants.py",
            "--input",
            str(concatenated_path),
            "--output",
            str(final_output),
            "--drop-duplicates",
            "--prefix",
            str(sort_scratch_prefix(Path(args.prefix))),
        )
    )


def main() -> int:
    args = parse_args()
    if "@" not in args.input:
        raise ValueError("prepare_variants_sharded.py requires --input to contain '@'")
    if "@" not in args.prefix:
        raise ValueError("prepare_variants_sharded.py requires --prefix to contain '@'")
    if "@" in args.output:
        raise ValueError("prepare_variants_sharded.py requires non-sharded --output without '@'")

    input_template = Path(args.input)
    prefix_template = Path(args.prefix)
    final_output = variant_object_path(Path(args.output))
    shards = resolve_import_input_paths(args.input, kind_label=KIND_LABEL_BY_FORMAT[args.input_format])
    managed_outputs = [*per_shard_managed_outputs(prefix_template, shards), final_output]

    if args.force:
        delete_outputs(managed_outputs)
    elif not args.resume:
        require_outputs_absent(managed_outputs)

    if args.resume and variant_object_exists(final_output):
        logger.info("%s: wrote %s", WRAPPER_NAME, final_output)
        return 0

    processed_any = False
    final_paths: list[Path] = []
    with tempfile.TemporaryDirectory(prefix="genomatch-prepare-sharded.") as workspace_raw:
        workspace = Path(workspace_raw)
        for shard in shards:
            output_prefix = shard_prefix(prefix_template, shard.source_shard)
            shard_final = variant_object_path(output_prefix)
            final_paths.append(shard_final)
            if args.resume and variant_object_exists(shard_final):
                logger.info("%s: skipping shard %s; output exists at %s", WRAPPER_NAME, shard.source_shard, shard_final)
                continue
            single_input_template = one_shard_template(input_template, shard, workspace)
            run_command(prepare_command(args, input_template=single_input_template, shard_output_prefix=output_prefix))
            processed_any = True

    metadata = validate_shard_metadata(final_paths)
    if args.resume and not processed_any and variant_object_exists(final_output):
        logger.info("%s: wrote %s", WRAPPER_NAME, final_output)
        return 0

    with tempfile.TemporaryDirectory(prefix="genomatch-prepare-sharded-merge.") as merge_workspace_raw:
        concatenated_path = Path(merge_workspace_raw) / "concatenated.vmap"
        write_concatenated_vmap(concatenated_path, final_paths, metadata, args.input)
        sort_final_vmap(args, concatenated_path=concatenated_path, final_output=final_output)

    logger.info("%s: wrote %s", WRAPPER_NAME, final_output)
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
