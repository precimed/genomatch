#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ._cli_utils import run_cli
from .vtable_utils import SUPPORTED_CONTIG_NAMINGS
from .workflow_wrapper_utils import (
    copy_variant_object,
    planned_existing_variant_outputs,
    print_skip_resolved_build,
    read_target_build,
    read_target_metadata,
    require_resume_no_hole,
    run_command,
    run_command_if_needed,
    tool_command,
    variant_object_exists,
    variant_object_path,
    delete_variant_object,
)

WRAPPER_NAME = "prepare_variants.py"
logger = logging.getLogger(__name__)


IMPORTER_BY_FORMAT = {
    "bim": "import_bim.py",
    "pvar": "import_pvar.py",
    "sumstats": "import_sumstats.py",
    "vcf": "import_vcf.py",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare raw input variants through the canonical import -> optional contig normalization -> "
            "metadata update -> restrict -> optional liftover -> optional strand/contig restriction "
            "chain, retaining stage outputs and copying the last retained stage to <output>.vmap."
        )
    )
    parser.add_argument("--input", help="Raw importer input (optional for --input-format=sumstats)")
    parser.add_argument("--input-format", required=True, choices=sorted(IMPORTER_BY_FORMAT), help="Importer selection")
    parser.add_argument("--sumstats-metadata", help="Summary-stat metadata YAML required for --input-format=sumstats")
    parser.add_argument("--id-vtable", help="Optional .vtable for --input-format=sumstats ID-based coordinate enrichment")
    parser.add_argument("--output", required=True, help="Output stem; the final prepared artifact is <output>.vmap")
    parser.add_argument(
        "--prefix",
        help="Optional base prefix for retained intermediates; defaults to --output",
    )
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
        help="Pass through to restrict_build_compatible.py (default: enabled)",
    )
    parser.add_argument(
        "--norm-indels",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Pass through to restrict_build_compatible.py (default: enabled)",
    )
    parser.add_argument("--drop-strand-ambiguous", action="store_true", help="Run drop_strand_ambiguous.py")
    parser.add_argument("--chr2use", "--contigs", dest="chr2use", help="If set, run restrict_contigs.py with this value")
    parser.add_argument("--max-allele-length", type=int, default=150, help="Maximum allele length (default: 150); passed through to selected importer")
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--resume", action="store_true", help="Resume by skipping already completed retained stages")
    mode_group.add_argument(
        "--force",
        action="store_true",
        help="Delete wrapper-managed outputs first, then rerun cleanly from scratch",
    )
    return parser.parse_args()


def importer_command(args: argparse.Namespace, output_path: Path) -> list[str]:
    importer = IMPORTER_BY_FORMAT[args.input_format]
    cmd = tool_command(importer, "--output", str(output_path))
    if args.input_format == "sumstats":
        if not args.sumstats_metadata:
            raise ValueError("--sumstats-metadata is required for --input-format=sumstats")
        if args.input:
            cmd.extend(["--input", args.input])
        cmd.extend(["--sumstats-metadata", args.sumstats_metadata])
        if args.id_vtable:
            cmd.extend(["--id-vtable", args.id_vtable])
    else:
        if not args.input:
            raise ValueError(f"--input is required for --input-format={args.input_format}")
        cmd.extend(["--input", args.input])
        if args.sumstats_metadata:
            raise ValueError("--sumstats-metadata is supported only for --input-format=sumstats")
        if args.id_vtable:
            raise ValueError("--id-vtable is supported only for --input-format=sumstats")
    cmd.extend(["--max-allele-length", str(args.max_allele_length)])
    return cmd


def transform_command(script_name: str, input_path: Path, output_path: Path, extra_args: list[str] | None = None) -> list[str]:
    input_flag = "--source" if script_name == "restrict_build_compatible.py" else "--input"
    cmd = tool_command(script_name, input_flag, str(input_path), "--output", str(output_path))
    if extra_args:
        cmd.extend(extra_args)
    return cmd


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


def delete_wrapper_managed_outputs(paths: list[Path]) -> None:
    for path in paths:
        delete_variant_object(path)


def all_retained_stage_outputs(prefix: Path) -> list[Path]:
    return [
        variant_object_path(prefix, ".imported"),
        variant_object_path(prefix, ".normalized"),
        variant_object_path(prefix, ".build_compatible"),
        variant_object_path(prefix, ".lifted"),
        variant_object_path(prefix, ".splitx"),
        variant_object_path(prefix, ".strand_filtered"),
        variant_object_path(prefix, ".contigs"),
    ]


def ensure_stage_ready_for_resume(output_path: Path, later_outputs: list[Path], *, resume: bool) -> None:
    if not resume or variant_object_exists(output_path):
        return
    require_resume_no_hole(WRAPPER_NAME, output_path, later_outputs)


def should_defer_plink_splitx(current_meta: dict[str, object], dst_contig_naming: str) -> bool:
    return dst_contig_naming == "plink_splitx" and str(current_meta.get("genome_build")) == "unknown"


def main() -> int:
    args = parse_args()
    prefix = Path(args.prefix) if args.prefix is not None else Path(args.output)
    final_output = variant_object_path(Path(args.output))

    imported_path = variant_object_path(prefix, ".imported")
    normalized_path = variant_object_path(prefix, ".normalized")
    restricted_path = variant_object_path(prefix, ".build_compatible")
    lifted_path = variant_object_path(prefix, ".lifted")
    splitx_path = variant_object_path(prefix, ".splitx")
    strand_path = variant_object_path(prefix, ".strand_filtered")
    contig_path = variant_object_path(prefix, ".contigs")

    planned_variant_outputs = [
        imported_path,
        normalized_path,
        restricted_path,
        lifted_path,
    ]
    if args.dst_contig_naming == "plink_splitx":
        planned_variant_outputs.append(splitx_path)
    if args.drop_strand_ambiguous:
        planned_variant_outputs.append(strand_path)
    if args.chr2use:
        planned_variant_outputs.append(contig_path)
    planned_variant_outputs.append(final_output)
    force_variant_outputs = [*all_retained_stage_outputs(prefix), final_output]

    if args.force:
        delete_wrapper_managed_outputs(force_variant_outputs)
    elif not args.resume:
        require_outputs_absent(planned_variant_outputs)

    if args.resume and variant_object_exists(final_output):
        logger.info("%s: wrote %s", WRAPPER_NAME, final_output)
        return 0

    for path in planned_variant_outputs:
        path.parent.mkdir(parents=True, exist_ok=True)

    run_command_if_needed(importer_command(args, imported_path), imported_path, resume=args.resume, wrapper_name=WRAPPER_NAME)
    current_path = imported_path
    current_meta = read_target_metadata(current_path)
    defer_splitx = should_defer_plink_splitx(current_meta, args.dst_contig_naming)
    initial_normalize_to = None
    if current_meta.get("contig_naming") != args.dst_contig_naming:
        if defer_splitx:
            if current_meta.get("contig_naming") is None:
                initial_normalize_to = "plink"
        else:
            initial_normalize_to = args.dst_contig_naming

    normalize_later_outputs = [
        restricted_path,
        lifted_path,
        *( [splitx_path] if defer_splitx else [] ),
        *( [strand_path] if args.drop_strand_ambiguous else [] ),
        *( [contig_path] if args.chr2use else [] ),
        final_output,
    ]
    if initial_normalize_to is not None:
        ensure_stage_ready_for_resume(normalized_path, normalize_later_outputs, resume=args.resume)
        run_command_if_needed(
            transform_command("normalize_contigs.py", current_path, normalized_path, ["--to", initial_normalize_to]),
            normalized_path,
            resume=args.resume,
            wrapper_name=WRAPPER_NAME,
        )
        current_path = normalized_path

    if args.resume and read_target_build(current_path) != "unknown":
        print_skip_resolved_build(WRAPPER_NAME, current_path)
    else:
        run_command(tool_command("guess_build.py", "--input", str(current_path), "--write"))

    restricted_later_outputs = [
        lifted_path,
        *( [splitx_path] if defer_splitx else [] ),
        *( [strand_path] if args.drop_strand_ambiguous else [] ),
        *( [contig_path] if args.chr2use else [] ),
        final_output,
    ]
    ensure_stage_ready_for_resume(restricted_path, restricted_later_outputs, resume=args.resume)
    current_build = read_target_build(current_path)
    restrict_extra_args: list[str] = []
    if args.allow_strand_flips:
        restrict_extra_args.append("--allow-strand-flips")
    if args.norm_indels:
        restrict_extra_args.append("--norm-indels")
    if current_build == args.dst_build:
        restrict_extra_args.append("--sort")
    run_command_if_needed(
        transform_command("restrict_build_compatible.py", current_path, restricted_path, restrict_extra_args),
        restricted_path,
        resume=args.resume,
        wrapper_name=WRAPPER_NAME,
    )
    current_path = restricted_path

    stage_after_build_outputs = [
        *( [splitx_path] if defer_splitx else [] ),
        *( [strand_path] if args.drop_strand_ambiguous else [] ),
        *( [contig_path] if args.chr2use else [] ),
        final_output,
    ]
    if current_build != args.dst_build:
        ensure_stage_ready_for_resume(lifted_path, stage_after_build_outputs, resume=args.resume)
        liftover_cmd = tool_command(
            "liftover_build.py",
            "--input",
            str(current_path),
            "--output",
            str(lifted_path),
            "--target-build",
            args.dst_build,
        )
        if args.resume:
            liftover_cmd.append("--resume")
        run_command_if_needed(liftover_cmd, lifted_path, resume=args.resume, wrapper_name=WRAPPER_NAME)
        current_path = lifted_path

    if defer_splitx:
        splitx_later_outputs = [
            *( [strand_path] if args.drop_strand_ambiguous else [] ),
            *( [contig_path] if args.chr2use else [] ),
            final_output,
        ]
        ensure_stage_ready_for_resume(splitx_path, splitx_later_outputs, resume=args.resume)
        run_command_if_needed(
            transform_command("normalize_contigs.py", current_path, splitx_path, ["--to", "plink_splitx"]),
            splitx_path,
            resume=args.resume,
            wrapper_name=WRAPPER_NAME,
        )
        current_path = splitx_path

    if args.drop_strand_ambiguous:
        strand_later_outputs = [*( [contig_path] if args.chr2use else [] ), final_output]
        ensure_stage_ready_for_resume(strand_path, strand_later_outputs, resume=args.resume)
        run_command_if_needed(
            transform_command("drop_strand_ambiguous.py", current_path, strand_path),
            strand_path,
            resume=args.resume,
            wrapper_name=WRAPPER_NAME,
        )
        current_path = strand_path

    if args.chr2use:
        ensure_stage_ready_for_resume(contig_path, [final_output], resume=args.resume)
        run_command_if_needed(
            transform_command("restrict_contigs.py", current_path, contig_path, ["--chr2use", args.chr2use]),
            contig_path,
            resume=args.resume,
            wrapper_name=WRAPPER_NAME,
        )
        current_path = contig_path

    copy_variant_object(current_path, final_output)
    logger.info("%s: wrote %s", WRAPPER_NAME, final_output)
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
