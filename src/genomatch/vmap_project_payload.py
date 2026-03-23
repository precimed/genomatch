#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ._cli_utils import run_cli
from .apply_vmap_utils import build_needed_source_indices, filtered_vmap_rows
from .contig_utils import supported_exact_contig_tokens
from .sample_axis_utils import (
    SAMPLE_ID_MODE_CHOICES,
    SAMPLE_ID_MODE_FID_IID,
    parse_fam_table,
    parse_psam_table,
    require_psam_fid_presence_consistent,
    sex_to_label,
)
from .vtable_utils import read_vmap
from .workflow_wrapper_utils import (
    delete_bfile_outputs,
    delete_pfile_outputs,
    delete_variant_object,
    existing_bfile_artifacts,
    existing_pfile_artifacts,
    existing_variant_object_artifacts,
    run_command,
    sidecar_output_path,
    tool_command,
    variant_object_path,
)

WRAPPER_NAME = "project_payload.py"
AUTO_PREFIX_TOKEN = "all_targets"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Project a raw payload into a prepared target by retaining the matched mapping as <prefix>.vmap "
            "and then applying it back to the original payload."
        )
    )
    parser.add_argument("--input", required=True, help="Raw payload input")
    parser.add_argument(
        "--input-format",
        required=True,
        choices=["bfile", "pfile", "sumstats", "sumstats-clean"],
        help="Payload input format",
    )
    parser.add_argument(
        "--sumstats-metadata",
        help="Summary-stat metadata YAML required for --input-format=sumstats or sumstats-clean",
    )
    parser.add_argument("--source-vmap", help="Prepared source .vmap")
    parser.add_argument("--target", required=True, help="Target .vtable or .vmap")
    parser.add_argument(
        "--output",
        required=True,
        help="Rewritten payload path (sumstats) or PLINK output prefix (bfile/pfile)",
    )
    parser.add_argument(
        "--fill-mode",
        choices=["column", "row"],
        default="column",
        help="Clean-sumstats fill behavior for --input-format=sumstats-clean",
    )
    parser.add_argument(
        "--use-af-inference",
        action="store_true",
        help="Enable AF-based clean-sumstats derivation rules for --input-format=sumstats-clean",
    )
    parser.add_argument("--prefix", help="Optional prefix for the retained matched mapping; defaults to --output")
    parser.add_argument(
        "--full-target",
        action="store_true",
        help="Keep unmatched target rows in the apply step instead of the wrapper default mapped-only projection",
    )
    parser.add_argument("--target-fam", help="Explicit target .fam for bfile input")
    parser.add_argument("--target-psam", help="Explicit target .psam for pfile input")
    parser.add_argument(
        "--sample-id-mode",
        choices=SAMPLE_ID_MODE_CHOICES,
        default="fid_iid",
        help="Subject-key mode for genotype-payload sample reconciliation",
    )
    parser.add_argument(
        "--sample-axis",
        choices=["union"],
        help="Wrapper-only sample-axis convenience for sharded genotype payloads",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Delete wrapper-managed outputs first, including the retained matched .vmap, then rerun cleanly",
    )
    return parser.parse_args()


def require_supported_args(args: argparse.Namespace) -> None:
    if not str(args.target).endswith((".vtable", ".vmap")):
        raise ValueError("--target must point to a .vtable or .vmap")
    if args.source_vmap is not None and not str(args.source_vmap).endswith(".vmap"):
        raise ValueError("--source-vmap must point to a .vmap")
    if Path(args.target).suffix == ".vtable" and args.source_vmap is None:
        raise ValueError("--source-vmap is required when --target points to a .vtable")
    if args.prefix is not None and "@" in str(args.prefix):
        raise ValueError("--prefix must not contain @")
    if args.input_format in {"sumstats", "sumstats-clean"}:
        if args.target_fam:
            raise ValueError("--target-fam is supported only for --input-format=bfile")
        if args.target_psam:
            raise ValueError("--target-psam is supported only for --input-format=pfile")
        if args.sample_axis is not None:
            raise ValueError("--sample-axis union is supported only for --input-format=bfile or pfile")
        if args.sample_id_mode != "fid_iid":
            raise ValueError("--sample-id-mode is supported only for --input-format=bfile or pfile")
        if not args.sumstats_metadata:
            raise ValueError(f"--sumstats-metadata is required for --input-format={args.input_format}")
        if "@" in str(args.output):
            raise ValueError(f"--output must not contain @ for --input-format={args.input_format}")
        if args.input_format != "sumstats-clean":
            if args.fill_mode != "column":
                raise ValueError("--fill-mode is supported only for --input-format=sumstats-clean")
            if args.use_af_inference:
                raise ValueError("--use-af-inference is supported only for --input-format=sumstats-clean")
    elif args.sumstats_metadata:
        raise ValueError("--sumstats-metadata is supported only for --input-format=sumstats or sumstats-clean")
    else:
        if args.fill_mode != "column":
            raise ValueError("--fill-mode is supported only for --input-format=sumstats-clean")
        if args.use_af_inference:
            raise ValueError("--use-af-inference is supported only for --input-format=sumstats-clean")
        if args.input_format == "bfile":
            if args.target_psam:
                raise ValueError("--target-psam is supported only for --input-format=pfile")
        if args.input_format == "pfile":
            if args.target_fam:
                raise ValueError("--target-fam is supported only for --input-format=bfile")
        if args.sample_axis is not None and (args.target_fam or args.target_psam):
            raise ValueError("--sample-axis union cannot be combined with --target-fam or --target-psam")


def bfile_source_prefix(input_path: str) -> str:
    path = Path(input_path)
    if path.suffix != ".bim":
        raise ValueError("project_payload.py bfile input requires --input to end with .bim")
    return str(path.with_suffix(""))


def pfile_source_prefix(input_path: str) -> str:
    path = Path(input_path)
    if path.suffix != ".pvar":
        raise ValueError("project_payload.py pfile input requires --input to end with .pvar")
    return str(path.with_suffix(""))


def resolved_prefix(args: argparse.Namespace) -> Path:
    if args.prefix is not None:
        return Path(args.prefix)
    if args.input_format in {"bfile", "pfile"} and "@" in str(args.output):
        return Path(str(args.output).replace("@", AUTO_PREFIX_TOKEN))
    return Path(args.output)


def use_target_vmap_directly(args: argparse.Namespace) -> bool:
    return args.source_vmap is None and Path(args.target).suffix == ".vmap"


def target_match_command(source_path: str, target_path: str, output_path: Path) -> list[str]:
    return tool_command(
        "match_vmap_to_target.py",
        "--source",
        source_path,
        "--target",
        target_path,
        "--output",
        str(output_path),
    )


def apply_command(args: argparse.Namespace, matched_path: Path) -> list[str]:
    if args.input_format in {"sumstats", "sumstats-clean"}:
        cmd = tool_command(
            "apply_vmap_to_sumstats.py",
            "--input",
            args.input,
            "--sumstats-metadata",
            args.sumstats_metadata,
            "--vmap",
            str(matched_path),
            "--output",
            args.output,
        )
        if args.input_format == "sumstats-clean":
            cmd.extend(["--clean", "--fill-mode", args.fill_mode])
            if args.use_af_inference:
                cmd.append("--use-af-inference")
    elif args.input_format == "bfile":
        cmd = tool_command(
            "apply_vmap_to_bfile.py",
            "--source-prefix",
            bfile_source_prefix(args.input),
            "--vmap",
            str(matched_path),
            "--output-prefix",
            args.output,
        )
        cmd.extend(["--sample-id-mode", args.sample_id_mode])
        if args.target_fam:
            cmd.extend(["--target-fam", args.target_fam])
    else:
        cmd = tool_command(
            "apply_vmap_to_pfile.py",
            "--source-prefix",
            pfile_source_prefix(args.input),
            "--vmap",
            str(matched_path),
            "--output-prefix",
            args.output,
        )
        cmd.extend(["--sample-id-mode", args.sample_id_mode])
        if args.target_psam:
            cmd.extend(["--target-psam", args.target_psam])
    if not args.full_target:
        cmd.append("--only-mapped-target")
    return cmd


def synthesized_target_sample_path(args: argparse.Namespace, prefix: Path) -> Path | None:
    if args.sample_axis != "union" or args.input_format not in {"bfile", "pfile"}:
        return None
    suffix = ".target_samples.fam" if args.input_format == "bfile" else ".target_samples.psam"
    return sidecar_output_path(prefix, suffix)


def source_prefix_for_args(args: argparse.Namespace) -> str:
    if args.input_format == "bfile":
        return bfile_source_prefix(args.input)
    if args.input_format == "pfile":
        return pfile_source_prefix(args.input)
    raise ValueError("source prefix is defined only for genotype payload inputs")


def discover_source_shards(source_prefix_arg: str, suffixes: tuple[str, ...]) -> dict[str, Path]:
    source_prefix = Path(source_prefix_arg)
    name_prefix, name_suffix = source_prefix.name.split("@", 1)
    discovered: dict[str, Path] = {}
    for token in supported_exact_contig_tokens():
        shard_prefix = source_prefix.parent / f"{name_prefix}{token}{name_suffix}"
        components = [Path(str(shard_prefix) + suffix) for suffix in suffixes]
        present = [component.exists() for component in components]
        if not any(present):
            continue
        if not all(present):
            missing = next(component for component in components if not component.exists())
            raise ValueError(f"incomplete source shard for token={token!r}: missing {missing}")
        discovered[token] = shard_prefix
    if not discovered:
        raise ValueError(f"no source shards found for template: {source_prefix}")
    return discovered


def build_union_target_sample(args: argparse.Namespace, prefix: Path, vmap_path: Path, *, allow_overwrite: bool) -> Path | None:
    if args.sample_axis != "union":
        return None
    source_prefix = source_prefix_for_args(args)
    if "@" not in source_prefix:
        print(
            f"{WRAPPER_NAME}: --sample-axis union is a no-op for non-sharded source input",
            file=sys.stderr,
        )
        return None
    vmap_rows = filtered_vmap_rows(read_vmap(vmap_path), only_mapped_target=not args.full_target)
    referenced_shards = sorted(build_needed_source_indices(vmap_rows))
    if len(referenced_shards) <= 1:
        print(
            f"{WRAPPER_NAME}: --sample-axis union is a no-op when only one referenced source shard remains",
            file=sys.stderr,
        )
        return None
    out_path = synthesized_target_sample_path(args, prefix)
    if out_path is None:
        return None

    if args.input_format == "bfile":
        discovered = discover_source_shards(source_prefix, (".bed", ".bim", ".fam"))
        tables = {}
        for shard in referenced_shards:
            shard_prefix = discovered.get(shard)
            if shard_prefix is None:
                raise ValueError(f"missing required source shard for source_shard={shard!r}")
            tables[shard] = parse_fam_table(
                Path(str(shard_prefix) + ".fam"),
                sample_id_mode=args.sample_id_mode,
                label="source .fam",
            )
    else:
        discovered = discover_source_shards(source_prefix, (".pgen", ".pvar", ".psam"))
        tables = {}
        for shard in referenced_shards:
            shard_prefix = discovered.get(shard)
            if shard_prefix is None:
                raise ValueError(f"missing required source shard for source_shard={shard!r}")
            tables[shard] = parse_psam_table(
                Path(str(shard_prefix) + ".psam"),
                sample_id_mode=args.sample_id_mode,
                label="source .psam",
            )
        if args.sample_id_mode == SAMPLE_ID_MODE_FID_IID:
            require_psam_fid_presence_consistent(tables.values(), label="participating source .psam files")
    if out_path.exists() and not allow_overwrite:
        raise ValueError(
            "planned output files already exist; rerun with --force to delete wrapper-managed outputs first:\n"
            f"- {out_path}"
        )

    ordered_keys: list[tuple[str, ...]] = []
    first_rows: dict[tuple[str, ...], object] = {}
    resolved_sex: dict[tuple[str, ...], int] = {}
    missing_plus_known_count = 0
    for shard in referenced_shards:
        table = tables[shard]
        for key, row in zip(table.keys, table.rows):
            if key not in first_rows:
                ordered_keys.append(key)
                first_rows[key] = row
                resolved_sex[key] = row.sex
                continue
            prev_sex = resolved_sex[key]
            if prev_sex in (1, 2) and row.sex in (1, 2) and prev_sex != row.sex:
                raise ValueError("conflicting known sex codes while synthesizing --sample-axis union")
            if prev_sex != row.sex and 0 in (prev_sex, row.sex):
                missing_plus_known_count += 1
                if prev_sex == 0 and row.sex in (1, 2):
                    resolved_sex[key] = row.sex
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if args.input_format == "bfile":
        lines = []
        for key in ordered_keys:
            row = first_rows[key]
            lines.append(f"{row.fid} {row.iid} 0 0 {sex_to_label(resolved_sex[key])} -9")
    else:
        include_sex = any(resolved_sex[key] in (1, 2) for key in ordered_keys)
        header = []
        if args.sample_id_mode == SAMPLE_ID_MODE_FID_IID:
            header.append("#FID")
        else:
            header.append("#IID")
        if args.sample_id_mode == SAMPLE_ID_MODE_FID_IID:
            header.append("IID")
        if include_sex:
            header.append("SEX")
        lines = ["\t".join(header)]
        for key in ordered_keys:
            row = first_rows[key]
            parts = []
            if args.sample_id_mode == SAMPLE_ID_MODE_FID_IID:
                parts.extend([row.fid, row.iid])
            else:
                parts.append(row.iid)
            if include_sex:
                parts.append(sex_to_label(resolved_sex[key]))
            lines.append("\t".join(parts))
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if missing_plus_known_count:
        print(
            f"Warning: {WRAPPER_NAME} --sample-axis union kept known sex when merging with missing sex for "
            f"{missing_plus_known_count} subject occurrences.",
            file=sys.stderr,
        )
    return out_path


def existing_wrapper_outputs(args: argparse.Namespace, matched_path: Path, *, retain_matched_vmap: bool) -> list[Path]:
    existing = existing_variant_object_artifacts(matched_path) if retain_matched_vmap else []
    if args.input_format in {"sumstats", "sumstats-clean"}:
        output_path = Path(args.output)
        if output_path.exists():
            existing.append(output_path)
        return existing
    if args.input_format == "bfile":
        return existing + existing_bfile_artifacts(Path(args.output), Path(args.target))
    return existing + existing_pfile_artifacts(Path(args.output), Path(args.target))


def delete_wrapper_outputs(args: argparse.Namespace, matched_path: Path, *, retain_matched_vmap: bool) -> None:
    if retain_matched_vmap:
        delete_variant_object(matched_path)
    synthesized_path = synthesized_target_sample_path(args, resolved_prefix(args))
    if synthesized_path is not None and synthesized_path.exists():
        synthesized_path.unlink()
    if args.input_format in {"sumstats", "sumstats-clean"}:
        output_path = Path(args.output)
        if output_path.exists():
            output_path.unlink()
        return
    if args.input_format == "bfile":
        delete_bfile_outputs(Path(args.output), Path(args.target))
        return
    delete_pfile_outputs(Path(args.output), Path(args.target))


def main() -> int:
    args = parse_args()
    require_supported_args(args)

    prefix = resolved_prefix(args)
    matched_path = variant_object_path(prefix)
    retain_matched_vmap = not use_target_vmap_directly(args)

    if args.force:
        delete_wrapper_outputs(args, matched_path, retain_matched_vmap=retain_matched_vmap)
    else:
        existing = existing_wrapper_outputs(args, matched_path, retain_matched_vmap=retain_matched_vmap)
        if existing:
            rendered = "\n".join(f"- {path}" for path in sorted(existing))
            raise ValueError(
                "planned output files already exist; rerun with --force to delete wrapper-managed outputs first:\n"
                f"{rendered}"
            )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    matched_created = False
    if retain_matched_vmap:
        matched_path.parent.mkdir(parents=True, exist_ok=True)
        assert args.source_vmap is not None
        run_command(target_match_command(args.source_vmap, args.target, matched_path))
        apply_vmap_path = matched_path
        matched_created = True
    else:
        print(
            f"{WRAPPER_NAME}: --source-vmap omitted and --target is a .vmap; skipping match_vmap_to_target.py",
            file=sys.stderr,
        )
        apply_vmap_path = Path(args.target)

    try:
        synthesized_path = build_union_target_sample(args, prefix, apply_vmap_path, allow_overwrite=args.force)
    except Exception:
        if matched_created and retain_matched_vmap:
            delete_variant_object(matched_path)
        raise
    if synthesized_path is not None:
        if args.input_format == "bfile":
            args.target_fam = str(synthesized_path)
        else:
            args.target_psam = str(synthesized_path)

    run_command(apply_command(args, apply_vmap_path))

    print(f"{WRAPPER_NAME}: wrote {args.output}", file=sys.stderr)
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
