#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import pysam
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
FIXTURE_DIR = REPO_ROOT / "examples" / "tutorial-1"
SCRATCH_DIR = Path(os.environ.get("TMPDIR", "/tmp")) / "genomatch-tutorial-1-generator"
SEED = 20260514
CHROMS = ("21", "22", "X")
ALLOWED_BASES = ("A", "C", "G", "T")
COMPLEMENT = str.maketrans("ACGT", "TGCA")

COUNTS_BY_CHROM = {
    "21": {"all": 24, "ref_study": 5, "ref_cohort": 5, "ref_only": 21, "study_cohort": 4, "study_only": 9, "cohort_only": 9},
    "22": {"all": 31, "ref_study": 7, "ref_cohort": 8, "ref_only": 29, "study_cohort": 6, "study_only": 12, "cohort_only": 11},
    "X": {"all": 29, "ref_study": 6, "ref_cohort": 7, "ref_only": 28, "study_cohort": 6, "study_only": 11, "cohort_only": 10},
}

ISSUE_LEGEND = {
    "allele_swap": "source alleles are reversed relative to the reference-oriented target row",
    "strand_flip": "source alleles are reverse-complemented relative to the reference-oriented target row",
    "strand_ambiguous": "SNP has an A/T or C/G allele pair",
    "indel": "row is an insertion/deletion represented with A/C/G/T alleles",
    "multiallelic": "VCF row has multiple ALT alleles and is dropped by the importer",
    "mixed_contig_label": "study rows intentionally mix accepted contig-label spellings",
    "duplicate_rsid": "same rsID appears at more than one coordinate",
    "placeholder_id": "input ID is '.'",
    "duplicate_target": "duplicate chrom/position/allele identity within one input",
    "liftover_collision": "two study rows collapse to the same target identity after preparation",
    "invalid_allele": "row contains a non-A/C/G/T allele and is dropped during import",
    "lowercase_allele": "lowercase alleles are accepted and normalized",
    "long_allele": "allele exceeds the tutorial max allele length and is dropped during import",
    "invalid_stat": "invalid clean summary-statistic values become missing in projection",
    "missing_genotype": "PLINK genotype row includes one or more missing calls",
    "x_par": "chrX PAR locus",
    "x_nonpar": "chrX non-PAR locus",
    "unsorted": "input rows are intentionally not sorted by coordinate",
    "reference_only": "present only in the target-universe VCF",
    "study_only": "present only in the study after preparation",
    "cohort_only": "present only in the cohort after preparation",
}


@dataclass(frozen=True)
class Candidate:
    chrom: str
    pos38: int
    ref38: str
    alt38: str
    pos37: int
    ref37: str
    alt37: str


@dataclass(frozen=True)
class Scenario:
    name: str
    category: str
    candidate: Candidate
    issue_tags: tuple[str, ...]
    reference: bool
    study: bool
    cohort: bool
    variant_kind: str = "snv"
    reference_id: str | None = None
    study_id: str | None = None
    cohort_id: str | None = None
    allele_mode_study: str = "identity"
    allele_mode_cohort: str = "identity"
    invalid_stat: bool = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate the checked-in Tutorial 1 fixture.")
    parser.add_argument("--skip-validation", action="store_true", help="Write fixture bytes without running tutorial CLI validation")
    return parser.parse_args()


def require_match_config() -> Path:
    raw = os.environ.get("MATCH_CONFIG")
    if not raw:
        raise SystemExit("MATCH_CONFIG must point to a config.yaml with real GRCh37/GRCh38 FASTA and chain assets")
    path = Path(raw)
    if not path.exists():
        raise SystemExit(f"MATCH_CONFIG does not exist: {path}")
    return path


def resolve_config_path(config_path: Path, raw: object, label: str) -> Path:
    if not isinstance(raw, str) or not raw:
        raise SystemExit(f"config entry {label} must be a non-empty path string")
    path = Path(raw)
    if not path.is_absolute():
        path = config_path.parent / path
    if not path.exists():
        raise SystemExit(f"required reference asset is missing: {label} -> {path}")
    return path


def load_assets(config_path: Path) -> dict[str, Path]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise SystemExit(f"config must be a mapping: {config_path}")
    builds = config.get("builds")
    liftover = config.get("liftover")
    if not isinstance(builds, dict) or not isinstance(liftover, list):
        raise SystemExit("config must use the documented builds/liftover shape")
    assets = {
        "grch37_fasta": resolve_config_path(config_path, builds.get("GRCh37", {}).get("ucsc_fasta"), "builds.GRCh37.ucsc_fasta"),
        "grch38_fasta": resolve_config_path(config_path, builds.get("GRCh38", {}).get("ucsc_fasta"), "builds.GRCh38.ucsc_fasta"),
    }
    for key, source, target in (("chain37to38", "GRCh37", "GRCh38"), ("chain38to37", "GRCh38", "GRCh37")):
        chain_value = None
        for edge in liftover:
            if isinstance(edge, dict) and edge.get("source") == source and edge.get("target") == target:
                chain_value = edge.get("chain")
                break
        if chain_value is None:
            raise SystemExit(f"config is missing liftover edge {source}->{target}")
        assets[key] = resolve_config_path(config_path, chain_value, f"liftover {source}->{target}")
    for fasta_key in ("grch37_fasta", "grch38_fasta"):
        if not Path(str(assets[fasta_key]) + ".fai").exists():
            raise SystemExit(f"FASTA index is missing: {assets[fasta_key]}.fai")
    return assets


def complement(allele: str) -> str:
    return allele.translate(COMPLEMENT)[::-1]


def alt_for(ref: str, index: int) -> str:
    if index % 17 == 0:
        return {"A": "T", "T": "A", "C": "G", "G": "C"}[ref]
    options = [base for base in ALLOWED_BASES if base != ref]
    return options[index % len(options)]


def fetch_base(fasta: pysam.FastaFile, chrom: str, pos: int) -> str:
    try:
        value = fasta.fetch(f"chr{chrom}", pos - 1, pos).upper()
    except Exception:
        return ""
    return value if value in ALLOWED_BASES else ""


def write_liftover_input(path: Path, rows: Iterable[tuple[str, int, str, str, str]]) -> None:
    rows = list(rows)
    contigs = sorted({row[0] for row in rows}, key=lambda value: (23 if value == "X" else int(value)))
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom in contigs:
            handle.write(f"##contig=<ID=chr{chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, row_id, ref, alt in rows:
            handle.write(f"chr{chrom}\t{pos}\t{row_id}\t{ref}\t{alt}\t.\tPASS\t.\n")


def run_bcftools_liftover(input_vcf: Path, output_vcf: Path, source_fasta: Path, target_fasta: Path, chain: Path) -> None:
    cmd = [
        "bcftools",
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
        str(chain),
        "--flip-tag",
        "FLIP",
        "--swap-tag",
        "SWAP",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"bcftools +liftover failed:\n{result.stderr or result.stdout}")


def parse_lifted(path: Path) -> dict[str, tuple[str, int, str, str]]:
    records: dict[str, tuple[str, int, str, str]] = {}
    duplicate_ids: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, pos, row_id, ref, alt = line.rstrip("\n").split("\t")[:5]
            chrom = chrom[3:] if chrom.startswith("chr") else chrom
            if row_id in records:
                duplicate_ids.add(row_id)
                continue
            records[row_id] = (chrom, int(pos), ref.upper(), alt.upper())
    for row_id in duplicate_ids:
        records.pop(row_id, None)
    return records


def liftover_candidate_batch(
    rows: list[tuple[str, int, str, str, str]],
    *,
    source_fasta: Path,
    target_fasta: Path,
    chain: Path,
    tmp_dir: Path,
    label: str,
) -> dict[str, tuple[str, int, str, str]]:
    input_vcf = tmp_dir / f"{label}.input.vcf"
    output_vcf = tmp_dir / f"{label}.output.vcf"
    write_liftover_input(input_vcf, rows)
    run_bcftools_liftover(input_vcf, output_vcf, source_fasta, target_fasta, chain)
    return parse_lifted(output_vcf)


def reverse_validated_candidates(
    chrom: str,
    desired: int,
    *,
    assets: dict[str, Path],
    tmp_dir: Path,
) -> list[Candidate]:
    fasta38 = pysam.FastaFile(str(assets["grch38_fasta"]))
    raw_rows: list[tuple[str, int, str, str, str]] = []
    try:
        if chrom == "X":
            starts = [50000, 5_000_000]
            strides = [7919, 10037]
            quotas = [max(24, desired // 4), desired * 12]
        else:
            chrom_int = int(chrom)
            starts = [100_000 + (SEED % 9973) + chrom_int * 173]
            strides = [10007 + chrom_int * 10]
            quotas = [desired * 12]
        row_index = 0
        target_rows = 0
        for start, stride, quota in zip(starts, strides, quotas):
            target_rows += quota
            pos = start
            while len(raw_rows) < target_rows:
                ref = fetch_base(fasta38, chrom, pos)
                if ref:
                    raw_rows.append((chrom, pos, f"{chrom}_{row_index}", ref, alt_for(ref, row_index)))
                    row_index += 1
                pos += stride
    finally:
        fasta38.close()

    lifted37 = liftover_candidate_batch(
        raw_rows,
        source_fasta=assets["grch38_fasta"],
        target_fasta=assets["grch37_fasta"],
        chain=assets["chain38to37"],
        tmp_dir=tmp_dir,
        label=f"{chrom}.38to37",
    )
    reverse_rows: list[tuple[str, int, str, str, str]] = []
    original_by_id = {row_id: (chrom_raw, pos, ref, alt) for chrom_raw, pos, row_id, ref, alt in raw_rows}
    for row_id, (chrom37, pos37, ref37, alt37) in lifted37.items():
        if chrom37 != chrom or len(ref37) != 1 or len(alt37) != 1 or ref37 not in ALLOWED_BASES or alt37 not in ALLOWED_BASES:
            continue
        reverse_rows.append((chrom37, pos37, row_id, ref37, alt37))
    lifted38 = liftover_candidate_batch(
        reverse_rows,
        source_fasta=assets["grch37_fasta"],
        target_fasta=assets["grch38_fasta"],
        chain=assets["chain37to38"],
        tmp_dir=tmp_dir,
        label=f"{chrom}.37to38",
    )
    candidates: list[Candidate] = []
    for row_id in sorted(lifted38, key=lambda value: int(value.rsplit("_", 1)[1])):
        chrom38, pos38, ref38, alt38 = lifted38[row_id]
        original = original_by_id[row_id]
        if (chrom38, pos38, ref38, alt38) != original:
            continue
        chrom37, pos37, ref37, alt37 = lifted37[row_id]
        candidates.append(Candidate(chrom, pos38, ref38, alt38, pos37, ref37, alt37))
        if len(candidates) == desired:
            return candidates
    raise RuntimeError(f"not enough uniquely round-tripping candidates on chr{chrom}: needed {desired}, got {len(candidates)}")


def build_scenarios(assets: dict[str, Path]) -> list[Scenario]:
    scenarios: list[Scenario] = []
    with tempfile.TemporaryDirectory(prefix="tutorial-1-liftover.") as tmp:
        tmp_dir = Path(tmp)
        for chrom in CHROMS:
            required = sum(COUNTS_BY_CHROM[chrom].values()) + 8
            candidates = reverse_validated_candidates(chrom, required, assets=assets, tmp_dir=tmp_dir)
            cursor = 0
            for category, count in COUNTS_BY_CHROM[chrom].items():
                for _ in range(count):
                    candidate = candidates[cursor]
                    cursor += 1
                    reference = category in {"all", "ref_study", "ref_cohort", "ref_only"}
                    study = category in {"all", "ref_study", "study_cohort", "study_only"}
                    cohort = category in {"all", "ref_cohort", "study_cohort", "cohort_only"}
                    category_tag_by_name = {
                        "ref_only": "reference_only",
                        "study_only": "study_only",
                        "cohort_only": "cohort_only",
                    }
                    tags = [category_tag_by_name[category]] if category in category_tag_by_name else []
                    if chrom == "X":
                        tags.append("x_par" if candidate.pos38 <= 2_781_479 else "x_nonpar")
                    if {candidate.ref38, candidate.alt38} in ({"A", "T"}, {"C", "G"}):
                        tags.append("strand_ambiguous")
                    idx = len(scenarios)
                    study_mode = ("identity", "swap", "flip", "flip_swap")[idx % 4]
                    cohort_mode = ("identity", "swap", "flip")[idx % 3]
                    if study and study_mode in {"swap", "flip_swap"}:
                        tags.append("allele_swap")
                    if study and study_mode in {"flip", "flip_swap"}:
                        tags.append("strand_flip")
                    if cohort and cohort_mode == "swap":
                        tags.append("allele_swap")
                    if cohort and cohort_mode == "flip":
                        tags.append("strand_flip")
                    scenarios.append(
                        Scenario(
                            name=f"t1_{chrom}_{idx:03d}",
                            category=category,
                            candidate=candidate,
                            issue_tags=tuple(dict.fromkeys(tags)),
                            reference=reference,
                            study=study,
                            cohort=cohort,
                            reference_id=f"rsTUT{100000 + idx}",
                            study_id=f"rsTUT{100000 + idx}",
                            cohort_id=f"rsTUT{100000 + idx}",
                            allele_mode_study=study_mode,
                            allele_mode_cohort=cohort_mode,
                        )
                    )
    scenarios = add_special_flags_and_indels(scenarios)
    return scenarios


def add_special_flags_and_indels(scenarios: list[Scenario]) -> list[Scenario]:
    out = list(scenarios)
    duplicate_rsid_indices = [3, 97]
    for idx in duplicate_rsid_indices:
        out[idx] = replace(out[idx], study_id="rsDUPLICATE", reference_id="rsDUPLICATE", issue_tags=tuple(dict.fromkeys((*out[idx].issue_tags, "duplicate_rsid"))))
    placeholder_idx = 11
    out[placeholder_idx] = replace(out[placeholder_idx], reference_id=".", cohort_id=".", issue_tags=tuple(dict.fromkeys((*out[placeholder_idx].issue_tags, "placeholder_id"))))
    invalid_stat_idx = next(i for i, row in enumerate(out) if row.study and row.reference)
    out[invalid_stat_idx] = replace(out[invalid_stat_idx], invalid_stat=True, issue_tags=tuple(dict.fromkeys((*out[invalid_stat_idx].issue_tags, "invalid_stat"))))

    indel_slots = [i for i, row in enumerate(out) if row.reference and row.cohort and not row.study][:4]
    for offset, idx in enumerate(indel_slots):
        row = out[idx]
        ref = row.candidate.ref38
        alt = ref + ("A" if ref != "A" else "C")
        candidate = replace(row.candidate, alt38=alt)
        out[idx] = replace(out[idx], candidate=candidate, variant_kind="insertion", issue_tags=tuple(dict.fromkeys((*row.issue_tags, "indel"))))
    return out


def source_alleles(ref: str, alt: str, mode: str) -> tuple[str, str]:
    if mode == "identity":
        return alt, ref
    if mode == "swap":
        return ref, alt
    if mode == "flip":
        return complement(alt), complement(ref)
    if mode == "flip_swap":
        return complement(ref), complement(alt)
    raise ValueError(mode)


def study_contig_label(chrom: str, source_index: int) -> str:
    if chrom == "21":
        return "chr21" if source_index % 3 else "21"
    if chrom == "22":
        return "22" if source_index % 2 else "chr22"
    return "chrX"


def write_reference_vcf(scenarios: list[Scenario], output: Path) -> None:
    rows = [row for row in scenarios if row.reference]
    extras = [
        ("chr21", 1_000_001, "t1_multiallelic_drop", "A", "C,G", ("multiallelic",)),
        ("chr22", 1_000_003, "t1_invalid_ref_drop", "N", "A", ("invalid_allele",)),
        ("chrX", 1_000_005, "t1_long_allele_drop", "A", "A" * 151, ("long_allele",)),
    ]
    rng = random.Random(SEED)
    rng.shuffle(rows)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom in CHROMS:
            handle.write(f"##contig=<ID=chr{chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            candidate = row.candidate
            handle.write(f"chr{candidate.chrom}\t{candidate.pos38}\t{row.reference_id or row.name}\t{candidate.ref38}\t{candidate.alt38}\t.\tPASS\t.\n")
        for chrom, pos, row_id, ref, alt, _tags in extras:
            handle.write(f"{chrom}\t{pos}\t{row_id}\t{ref}\t{alt}\t.\tPASS\t.\n")


def write_sumstats(scenarios: list[Scenario], output: Path) -> None:
    rows = [row for row in scenarios if row.study]
    collision_source, lowercase_source = study_duplicate_source_rows(scenarios)
    rows.append(replace(collision_source, name="t1_liftover_collision_partner", study_id="rsTUT_COLLISION", allele_mode_study="swap", issue_tags=("liftover_collision", "duplicate_target")))
    rows.append(replace(lowercase_source, name="t1_lowercase_alleles", study_id="rsTUT_LOWER", allele_mode_study="identity", issue_tags=("lowercase_allele", "duplicate_target")))
    rng = random.Random(SEED + 1)
    rng.shuffle(rows)
    header = [
        "marker",
        "chromosome",
        "base_pair",
        "effect_allele_raw",
        "other_allele_raw",
        "odds_ratio",
        "p_value",
        "case_count",
        "control_count",
        "case_other_af",
        "control_other_af",
        "info_score",
        "direction",
    ]
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for source_index, row in enumerate(rows):
            candidate = row.candidate
            ea, oa = source_alleles(candidate.ref37, candidate.alt37, row.allele_mode_study)
            if "lowercase_allele" in row.issue_tags:
                ea, oa = ea.lower(), oa.lower()
            p_value = "1.4" if row.invalid_stat else f"{0.0001 + (source_index % 97) / 100000:.7f}"
            case_oaf = "1.2" if row.invalid_stat else f"{0.15 + (source_index % 40) / 100:.3f}"
            writer.writerow(
                [
                    row.study_id or row.name,
                    study_contig_label(candidate.chrom, source_index),
                    candidate.pos37,
                    ea,
                    oa,
                    f"{1.02 + (source_index % 20) / 100:.3f}",
                    p_value,
                    420 + (source_index % 17),
                    690 + (source_index % 19),
                    case_oaf,
                    f"{0.20 + (source_index % 35) / 100:.3f}",
                    f"{0.82 + (source_index % 12) / 100:.3f}",
                    "+" if source_index % 2 == 0 else "-",
                ]
            )
        writer.writerow(["t1_invalid_allele_drop", "chr21", 101001, "I", "A", "1.1", "0.01", "400", "600", "0.2", "0.3", "0.9", "+"])
        writer.writerow(["t1_long_allele_drop", "chr22", 101003, "A" * 151, "C", "1.1", "0.01", "400", "600", "0.2", "0.3", "0.9", "+"])


def write_sumstats_metadata(path: Path) -> None:
    text = """cleansumstats_metafile_kind: minimal
path_sumStats: study.grch37.tsv
delimiter: "\\t"
missing_value: ""
genome_build: GRCh37
contig_naming: ucsc
stats_Model: logistic
stats_TraitType: case-control
stats_neglog10P: false
stats_log10P: false
col_SNP: marker
col_CHR: chromosome
col_POS: base_pair
col_EffectAllele: effect_allele_raw
col_OtherAllele: other_allele_raw
col_OR: odds_ratio
col_P: p_value
col_CaseN: case_count
col_ControlN: control_count
col_CaseOAF: case_other_af
col_ControlOAF: control_other_af
col_INFO: info_score
col_Direction: direction
"""
    path.write_text(text, encoding="utf-8")


def fam_lines() -> list[str]:
    return [
        "F1 I1 0 0 1 -9",
        "F2 I2 0 0 2 -9",
        "F3 I3 0 0 0 -9",
        "F4 I4 0 0 1 -9",
        "F5 I5 0 0 2 -9",
        "F6 I6 0 0 1 -9",
    ]


def genotype_row(row_index: int) -> list[int]:
    patterns = [
        [0, 2, 3, 1, 0, 3],
        [3, 0, 2, 1, 3, 0],
        [0, 1, 0, 3, 2, 1],
        [1, 3, 0, 2, 0, 3],
    ]
    return patterns[row_index % len(patterns)]


def write_bed_matrix(path: Path, matrix: list[list[int]], n_samples: int) -> None:
    from genomatch.bfile_utils import write_bed_matrix as write_matrix

    write_matrix(path, matrix, n_samples)


def plink_path(prefix: Path, suffix: str) -> Path:
    return Path(str(prefix) + suffix)


def write_cohort_bfiles(scenarios: list[Scenario], output_dir: Path) -> None:
    rows_by_chrom: dict[str, list[Scenario]] = {chrom: [] for chrom in CHROMS}
    for row in scenarios:
        if row.cohort:
            rows_by_chrom[row.candidate.chrom].append(row)
    duplicate_source = next(row for row in scenarios if row.cohort and row.reference)
    rows_by_chrom[duplicate_source.candidate.chrom].append(
        replace(duplicate_source, name="t1_duplicate_bim_key", cohort_id="rsTUT_DUPKEY", issue_tags=("duplicate_target",))
    )
    rng = random.Random(SEED + 2)
    for chrom in CHROMS:
        rows = rows_by_chrom[chrom]
        rng.shuffle(rows)
        bim_lines: list[str] = []
        matrix: list[list[int]] = []
        for idx, row in enumerate(rows):
            candidate = row.candidate
            a1, a2 = source_alleles(candidate.ref38, candidate.alt38, row.allele_mode_cohort)
            tags = set(row.issue_tags)
            if idx % 23 == 0:
                tags.add("missing_genotype")
            bim_lines.append("\t".join([chrom, row.cohort_id or row.name, "0", str(candidate.pos38), a1, a2]))
            matrix.append(genotype_row(idx))
        prefix = output_dir / f"cohort.{chrom}"
        plink_path(prefix, ".bim").write_text("\n".join(bim_lines) + "\n", encoding="utf-8")
        plink_path(prefix, ".fam").write_text("\n".join(fam_lines()) + "\n", encoding="utf-8")
        write_bed_matrix(plink_path(prefix, ".bed"), matrix, 6)


def study_duplicate_source_rows(scenarios: list[Scenario]) -> tuple[Scenario, Scenario]:
    rows = [row for row in scenarios if row.study]
    collision_source = next(row for row in rows if row.reference and row.study and row.variant_kind == "snv")
    lowercase_source = rows[7]
    return collision_source, lowercase_source


def dropped_study_row_names(scenarios: list[Scenario]) -> set[str]:
    duplicate_names = {row.name for row in study_duplicate_source_rows(scenarios)}
    strand_ambiguous_names = {
        row.name
        for row in scenarios
        if row.study and "strand_ambiguous" in row.issue_tags
    }
    return duplicate_names | strand_ambiguous_names


def effective_membership(row: Scenario, dropped_study_row_names: set[str]) -> tuple[bool, bool, bool]:
    return row.reference, row.study and row.name not in dropped_study_row_names, row.cohort


def venn_counts(scenarios: list[Scenario]) -> dict[str, int]:
    dropped_names = dropped_study_row_names(scenarios)
    effective = [effective_membership(row, dropped_names) for row in scenarios]
    return {
        "reference": sum(reference for reference, _study, _cohort in effective),
        "study": sum(study for _reference, study, _cohort in effective),
        "cohort": sum(cohort for _reference, _study, cohort in effective),
        "3-way": sum(reference and study and cohort for reference, study, cohort in effective),
        "ref+study": sum(reference and study and not cohort for reference, study, cohort in effective),
        "ref+cohort": sum(reference and cohort and not study for reference, study, cohort in effective),
        "study+cohort": sum(study and cohort and not reference for reference, study, cohort in effective),
        "reference_only": sum(reference and not study and not cohort for reference, study, cohort in effective),
        "study_only": sum(study and not reference and not cohort for reference, study, cohort in effective),
        "cohort_only": sum(cohort and not reference and not study for reference, study, cohort in effective),
    }


def write_manifest(scenarios: list[Scenario], output: Path) -> None:
    dropped_names = dropped_study_row_names(scenarios)
    duplicate_drop_names = {row.name for row in study_duplicate_source_rows(scenarios)}
    fieldnames = [
        "row_id",
        "membership",
        "chrom38",
        "pos38",
        "ref38",
        "alt38",
        "chrom37",
        "pos37",
        "ref37",
        "alt37",
        "reference_id",
        "study_id",
        "cohort_id",
        "issue_tags",
        "notes",
    ]
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in scenarios:
            candidate = row.candidate
            reference, study, cohort = effective_membership(row, dropped_names)
            membership = "+".join(name for name, enabled in (("reference", reference), ("study", study), ("cohort", cohort)) if enabled)
            issue_tags = tuple(row.issue_tags)
            notes = "generated from real local FASTA and UCSC chain liftover"
            note_parts: list[str] = []
            if row.name in duplicate_drop_names:
                issue_tags = tuple(dict.fromkeys((*issue_tags, "duplicate_target")))
                note_parts.append(
                    "raw study row exists, but the intentional duplicate-target study example "
                    "causes build-compatible preparation to drop it from effective study membership"
                )
            if row.study and "strand_ambiguous" in row.issue_tags:
                note_parts.append(
                    "raw study row exists, but Tutorial 1 uses --drop-strand-ambiguous so it is excluded "
                    "from effective study membership"
                )
            if note_parts:
                notes = "; ".join(note_parts)
            writer.writerow(
                {
                    "row_id": row.name,
                    "membership": membership,
                    "chrom38": candidate.chrom,
                    "pos38": candidate.pos38,
                    "ref38": candidate.ref38,
                    "alt38": candidate.alt38,
                    "chrom37": f"chr{candidate.chrom}",
                    "pos37": candidate.pos37,
                    "ref37": candidate.ref37,
                    "alt37": candidate.alt37,
                    "reference_id": row.reference_id or "",
                    "study_id": row.study_id or "",
                    "cohort_id": row.cohort_id or "",
                    "issue_tags": ",".join(issue_tags),
                    "notes": notes,
                }
            )
        extra_rows = [
            ("t1_multiallelic_drop", "reference", "multiallelic", "VCF importer reports multiallelic"),
            ("t1_invalid_allele_drop", "study/reference", "invalid_allele", "importers report non_actg_allele"),
            ("t1_long_allele_drop", "study/reference", "long_allele", "importers report allele_too_long"),
            ("t1_liftover_collision_partner", "study", "liftover_collision,duplicate_target", "reverse allele row collapses during preparation"),
            ("t1_lowercase_alleles", "study", "lowercase_allele,duplicate_target", "lowercase study alleles are accepted before duplicate-target QC"),
            ("t1_duplicate_bim_key", "cohort", "duplicate_target", "duplicate BIM target identity is reported during cohort preparation"),
            ("fixture_unsorted_inputs", "reference/study/cohort", "unsorted", "raw input rows are shuffled deterministically before writing"),
            ("fixture_mixed_study_contigs", "study", "mixed_contig_label", "study chromosome column intentionally mixes chr-prefixed and plain labels"),
            ("fixture_missing_genotypes", "cohort", "missing_genotype", "deterministic PLINK BED rows include missing genotype calls"),
        ]
        for row_id, membership, tags, notes in extra_rows:
            writer.writerow({key: "" for key in fieldnames} | {"row_id": row_id, "membership": membership, "issue_tags": tags, "notes": notes})


def write_fixture_readme(counts: dict[str, int], output: Path) -> None:
    lines = [
        "# Tutorial 1 Fixture",
        "",
        "This directory contains a small, deterministic fixture for the Tutorial 1 workflow.",
        "",
        "- `reference.grch38.vcf` is the GRCh38 target universe and uses UCSC contig labels.",
        "- `study.grch37.tsv` is a GRCh37 summary-statistics payload with non-canonical column names.",
        "- `cohort.{21,22,X}.{bed,bim,fam}` is a sharded GRCh38 PLINK 1 BFILE payload with PLINK contig labels.",
        "",
        "The fixture is intentionally messy. Warnings during preparation and projection are expected; inspect `manifest.tsv` for the reason each special row exists.",
        "",
        "## Effective Membership Counts",
        "",
        "These counts are after preparation QC and the tutorial's study-side strand-ambiguous filter, so intentionally duplicated study rows and A/T or C/G study SNPs are excluded from study membership.",
        "",
        f"`reference={counts['reference']}  study={counts['study']}  cohort={counts['cohort']}`",
        "",
        "| subset | count |",
        "| --- | ---: |",
        f"| reference + study + cohort | {counts['3-way']} |",
        f"| reference + study only | {counts['ref+study']} |",
        f"| reference + cohort only | {counts['ref+cohort']} |",
        f"| study + cohort only | {counts['study+cohort']} |",
        f"| reference only | {counts['reference_only']} |",
        f"| study only | {counts['study_only']} |",
        f"| cohort only | {counts['cohort_only']} |",
        "",
        "```text",
        "                  study",
        "             ref+study       study+cohort",
        "reference        3-way             cohort",
        "             ref+cohort",
        "```",
        "",
        "## issue_tags",
        "",
    ]
    for tag, description in sorted(ISSUE_LEGEND.items()):
        lines.append(f"- `{tag}`: {description}")
    lines.extend(["", "See ../../docs/tutorial-1.md for the end-to-end commands."])
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_static_files(fixture_dir: Path) -> None:
    for chrom in CHROMS:
        prefix = fixture_dir / f"cohort.{chrom}"
        if plink_path(prefix, ".bed").read_bytes()[:3] != b"\x6c\x1b\x01":
            raise RuntimeError(f"invalid PLINK BED header: {plink_path(prefix, '.bed')}")
        for line in plink_path(prefix, ".bim").read_text(encoding="utf-8").splitlines():
            if line.strip() and len(line.split()) < 6:
                raise RuntimeError(f"invalid BIM row: {line}")
    text = (fixture_dir / "reference.grch38.vcf").read_text(encoding="utf-8")
    if "#CHROM\tPOS\tID\tREF\tALT" not in text:
        raise RuntimeError("reference VCF header is missing #CHROM")
    metadata = yaml.safe_load((fixture_dir / "study.grch37.meta.yaml").read_text(encoding="utf-8"))
    for key in ("col_SNP", "col_CHR", "col_POS", "col_EffectAllele", "col_OtherAllele"):
        if key not in metadata:
            raise RuntimeError(f"sumstats metadata missing {key}")


def run_cli(cmd: list[str], *, cwd: Path, env: dict[str, str]) -> None:
    result = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        rendered = " ".join(cmd)
        raise RuntimeError(f"validation command failed: {rendered}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")


def validate_roundtrip(fixture_dir: Path) -> None:
    if SCRATCH_DIR.exists():
        shutil.rmtree(SCRATCH_DIR)
    try:
        (SCRATCH_DIR / "work").mkdir(parents=True)
        (SCRATCH_DIR / "out").mkdir(parents=True)
        env = os.environ.copy()
        env["PYTHONPATH"] = str(SRC_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
        commands = [
            [sys.executable, "-m", "genomatch.prepare_variants", "--input", str(fixture_dir / "reference.grch38.vcf"), "--input-format", "vcf", "--dst-build", "GRCh38", "--dst-contig-naming", "ncbi", "--output", str(SCRATCH_DIR / "work/reference")],
            [sys.executable, "-m", "genomatch.prepare_variants", "--input-format", "sumstats", "--sumstats-metadata", str(fixture_dir / "study.grch37.meta.yaml"), "--dst-build", "GRCh38", "--dst-contig-naming", "ncbi", "--drop-strand-ambiguous", "--output", str(SCRATCH_DIR / "work/study")],
            [sys.executable, "-m", "genomatch.prepare_variants_sharded", "--input", str(fixture_dir / "cohort.@.bim"), "--input-format", "bim", "--prefix", str(SCRATCH_DIR / "work/cohort.@"), "--output", str(SCRATCH_DIR / "work/cohort"), "--dst-build", "GRCh38", "--dst-contig-naming", "ncbi", "--shards", "21,22,X"],
            [sys.executable, "-m", "genomatch.restrict_vmap", str(SCRATCH_DIR / "work/study.vmap"), str(SCRATCH_DIR / "work/reference.vmap"), "--output", str(SCRATCH_DIR / "work/study.reference.vmap")],
            [sys.executable, "-m", "genomatch.restrict_vmap", str(SCRATCH_DIR / "work/cohort.vmap"), str(SCRATCH_DIR / "work/reference.vmap"), "--output", str(SCRATCH_DIR / "work/cohort.reference.vmap")],
            [sys.executable, "-m", "genomatch.project_payload", "--input-format", "sumstats-clean", "--sumstats-metadata", str(fixture_dir / "study.grch37.meta.yaml"), "--vmap", str(SCRATCH_DIR / "work/study.reference.vmap"), "--output", str(SCRATCH_DIR / "out/study.reference.tsv"), "--use-af-inference"],
            [sys.executable, "-m", "genomatch.project_payload", "--input", str(fixture_dir / "cohort.@.bim"), "--input-format", "bfile", "--vmap", str(SCRATCH_DIR / "work/cohort.reference.vmap"), "--output", str(SCRATCH_DIR / "out/cohort.@"), "--skip-ploidy-check"],
        ]
        for cmd in commands:
            run_cli(cmd, cwd=REPO_ROOT, env=env)
    finally:
        if SCRATCH_DIR.exists():
            shutil.rmtree(SCRATCH_DIR)


def main() -> int:
    args = parse_args()
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))
    config_path = require_match_config()
    assets = load_assets(config_path)
    scenarios = build_scenarios(assets)
    counts = venn_counts(scenarios)

    if FIXTURE_DIR.exists():
        shutil.rmtree(FIXTURE_DIR)
    FIXTURE_DIR.mkdir(parents=True)
    write_reference_vcf(scenarios, FIXTURE_DIR / "reference.grch38.vcf")
    write_sumstats(scenarios, FIXTURE_DIR / "study.grch37.tsv")
    write_sumstats_metadata(FIXTURE_DIR / "study.grch37.meta.yaml")
    write_cohort_bfiles(scenarios, FIXTURE_DIR)
    write_manifest(scenarios, FIXTURE_DIR / "manifest.tsv")
    write_fixture_readme(counts, FIXTURE_DIR / "README.md")
    validate_static_files(FIXTURE_DIR)
    if not args.skip_validation:
        validate_roundtrip(FIXTURE_DIR)
    print(f"reference={counts['reference']}  study={counts['study']}  cohort={counts['cohort']}")
    print(f"3-way={counts['3-way']}  ref+study={counts['ref+study']}  ref+cohort={counts['ref+cohort']}  study+cohort={counts['study+cohort']}")
    print(f"reference_only={counts['reference_only']}  study_only={counts['study_only']}  cohort_only={counts['cohort_only']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
