from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

from vtable_utils import load_variant_object, metadata_path_for


def variant_object_path(prefix: Path, suffix: str = "") -> Path:
    base = str(prefix)
    if base.endswith(".vmap"):
        stem = base[:-5]
    else:
        stem = base
    return Path(f"{stem}{suffix}.vmap")


def sidecar_output_path(prefix: Path, suffix: str) -> Path:
    base = str(prefix)
    if base.endswith(".vmap"):
        stem = base[:-5]
    else:
        stem = base
    return Path(f"{stem}{suffix}")


def print_command(cmd: list[str]) -> None:
    rendered = " ".join(str(part) for part in cmd)
    print(f"+ {rendered}", file=sys.stderr)


def run_command(cmd: list[str]) -> None:
    print_command(cmd)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.stdout:
        print(result.stdout, end="", file=sys.stderr)
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    if result.returncode != 0:
        raise ValueError(f"{Path(cmd[1]).name} failed with exit code {result.returncode}")


def variant_object_exists(path: Path) -> bool:
    return path.exists() and metadata_path_for(path).exists()


def variant_object_artifacts(path: Path) -> list[Path]:
    return [
        path,
        metadata_path_for(path),
        path.with_name(path.name + ".qc.tsv"),
    ]


def existing_variant_object_artifacts(path: Path) -> list[Path]:
    return [artifact for artifact in variant_object_artifacts(path) if artifact.exists()]


def delete_variant_object(path: Path) -> None:
    for artifact in existing_variant_object_artifacts(path):
        artifact.unlink()


def copy_variant_object(source_path: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    delete_variant_object(output_path)
    shutil.copy2(source_path, output_path)
    shutil.copy2(metadata_path_for(source_path), metadata_path_for(output_path))
    qc_source = source_path.with_name(source_path.name + ".qc.tsv")
    qc_output = output_path.with_name(output_path.name + ".qc.tsv")
    if qc_source.exists():
        shutil.copy2(qc_source, qc_output)
    elif qc_output.exists():
        qc_output.unlink()


def read_target_metadata(path: Path) -> dict[str, object]:
    metadata = json.loads(metadata_path_for(path).read_text(encoding="utf-8"))
    if metadata.get("object_type") == "variant_map":
        return dict(metadata["target"])
    return dict(metadata)


def read_target_build(path: Path) -> str:
    return str(read_target_metadata(path)["genome_build"])


def print_skip(wrapper_name: str, step_label: str, output_path: Path) -> None:
    print(f"{wrapper_name}: skipping {step_label}; output exists at {output_path}", file=sys.stderr)


def print_skip_resolved_build(wrapper_name: str, path: Path) -> None:
    print(
        f"{wrapper_name}: skipping guess_build.py; genome_build is already resolved in {path}",
        file=sys.stderr,
    )


def run_command_if_needed(cmd: list[str], output_path: Path, *, resume: bool, wrapper_name: str) -> None:
    if resume and variant_object_exists(output_path):
        print_skip(wrapper_name, Path(cmd[1]).name, output_path)
        return
    run_command(cmd)


def run_command_if_needed_plain_output(cmd: list[str], output_path: Path, *, force: bool, wrapper_name: str) -> None:
    if not force and output_path.exists():
        raise ValueError(f"{wrapper_name}: planned output already exists: {output_path}")
    run_command(cmd)


BFILE_REQUIRED_SUFFIXES = (".bed", ".bim", ".fam")
BFILE_OPTIONAL_SUFFIXES = (".ploidy",)
PFILE_REQUIRED_SUFFIXES = (".pgen", ".pvar", ".psam")


def expected_bfile_output_prefixes(output_prefix: Path, target_path: Path) -> list[Path]:
    if "@" not in str(output_prefix):
        return [output_prefix]
    chrom_order: list[str] = []
    for row in load_variant_object(target_path).target_rows:
        if row.chrom not in chrom_order:
            chrom_order.append(row.chrom)
    return [Path(str(output_prefix).replace("@", chrom_label)) for chrom_label in chrom_order]


def planned_bfile_artifacts(output_prefix: Path, target_path: Path | None = None) -> list[Path]:
    if "@" in str(output_prefix):
        if target_path is None:
            raise ValueError("target_path is required to plan sharded PLINK outputs")
        prefixes = expected_bfile_output_prefixes(output_prefix, target_path)
    else:
        prefixes = [output_prefix]
    artifacts: list[Path] = []
    for prefix in prefixes:
        for suffix in (*BFILE_REQUIRED_SUFFIXES, *BFILE_OPTIONAL_SUFFIXES):
            artifacts.append(sidecar_output_path(prefix, suffix))
    return artifacts


def existing_bfile_artifacts(output_prefix: Path, target_path: Path | None = None) -> list[Path]:
    return [artifact for artifact in planned_bfile_artifacts(output_prefix, target_path) if artifact.exists()]


def delete_bfile_outputs(output_prefix: Path, target_path: Path | None = None) -> None:
    for artifact in existing_bfile_artifacts(output_prefix, target_path):
        artifact.unlink()


def bfile_output_exists(output_prefix: Path, target_path: Path) -> bool:
    expected_prefixes = expected_bfile_output_prefixes(output_prefix, target_path)
    if not expected_prefixes:
        return False
    return all(
        all(sidecar_output_path(shard_prefix, suffix).exists() for suffix in BFILE_REQUIRED_SUFFIXES)
        for shard_prefix in expected_prefixes
    )


def planned_pfile_artifacts(output_prefix: Path, target_path: Path | None = None) -> list[Path]:
    if "@" in str(output_prefix):
        if target_path is None:
            raise ValueError("target_path is required to plan sharded PLINK outputs")
        prefixes = expected_bfile_output_prefixes(output_prefix, target_path)
    else:
        prefixes = [output_prefix]
    artifacts: list[Path] = []
    for prefix in prefixes:
        for suffix in PFILE_REQUIRED_SUFFIXES:
            artifacts.append(sidecar_output_path(prefix, suffix))
    return artifacts


def existing_pfile_artifacts(output_prefix: Path, target_path: Path | None = None) -> list[Path]:
    return [artifact for artifact in planned_pfile_artifacts(output_prefix, target_path) if artifact.exists()]


def delete_pfile_outputs(output_prefix: Path, target_path: Path | None = None) -> None:
    for artifact in existing_pfile_artifacts(output_prefix, target_path):
        artifact.unlink()


def planned_existing_variant_outputs(paths: list[Path]) -> list[Path]:
    existing: list[Path] = []
    for path in paths:
        existing.extend(existing_variant_object_artifacts(path))
    return existing


def require_resume_no_hole(wrapper_name: str, missing_output: Path, later_outputs: list[Path]) -> None:
    existing: list[Path] = []
    for output_path in later_outputs:
        existing.extend(existing_variant_object_artifacts(output_path))
    if not existing:
        return
    rendered = "\n".join(f"- {path}" for path in sorted(existing))
    raise ValueError(
        f"{wrapper_name} --resume found a hole in retained stage outputs; missing required stage {missing_output} "
        "while later planned outputs already exist:\n"
        f"{rendered}"
    )
