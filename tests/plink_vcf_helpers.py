from __future__ import annotations

from dataclasses import dataclass
import json
import random
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from utils import REPO_ROOT, read_bim, read_tsv, run_cmd, run_py, write_json, write_lines


BASES = ("A", "C", "G", "T")
PLINK_SPECIAL = {"23": "X", "24": "Y", "25": "X", "26": "MT"}


@dataclass(frozen=True)
class VcfVariant:
    chrom: str
    pos: int
    vid: str
    ref: str
    alt: str
    genotypes: List[str]


@dataclass(frozen=True)
class HaploidVariant:
    chrom: str
    snp: str
    pos: int
    a1: str
    a2: str
    region: str


def assert_plink_ok(
    result,
    allow_warnings: bool = False,
    allowed_warning_substrings: Optional[List[str]] = None,
) -> None:
    if result.returncode != 0:
        raise AssertionError(f"plink failed: {result.stderr}\n{result.stdout}")
    text = result.stdout + result.stderr
    if "error" in text.lower():
        raise AssertionError(f"plink emitted warnings:\n{text}")
    if "warning" in text.lower():
        if not allow_warnings:
            raise AssertionError(f"plink emitted warnings:\n{text}")
        if allowed_warning_substrings:
            for line in text.splitlines():
                if "warning" not in line.lower():
                    continue
                if not any(sub.lower() in line.lower() for sub in allowed_warning_substrings):
                    raise AssertionError(f"unexpected plink warning:\n{line}\nFull output:\n{text}")


def assert_py_ok(result) -> None:
    if result.returncode != 0:
        raise AssertionError(f"script failed: {result.stderr}\n{result.stdout}")


def normalize_test_chrom(label: str, *, to_naming: str) -> str:
    token = label.strip()
    if to_naming == "plink":
        if token.isdigit():
            try:
                value = int(token)
            except ValueError:
                value = -1
            if 1 <= value <= 24 or token == "26":
                return str(value)
        token_upper = token.upper()
        if token_upper in {"X", "Y", "MT"}:
            return {"X": "23", "Y": "24", "MT": "26"}[token_upper]
    if to_naming == "plink_splitx":
        if token.isdigit():
            try:
                value = int(token)
            except ValueError:
                value = -1
            if 1 <= value <= 26:
                return str(value)
        token_upper = token.upper()
        if token_upper in {"X", "Y", "MT"}:
            return {"X": "23", "Y": "24", "MT": "26"}[token_upper]
    token_lower = token.lower()
    if token_lower.startswith("chr"):
        token_lower = token_lower[3:]
        if token_lower == "m":
            canonical = "MT"
        elif token_lower in {"x", "y"}:
            canonical = token_lower.upper()
        elif token_lower.isdigit() and 1 <= int(token_lower) <= 22:
            canonical = str(int(token_lower))
        else:
            canonical = token
    elif token in {"X", "Y", "MT"}:
        canonical = token
    elif token.isdigit():
        canonical = PLINK_SPECIAL.get(token, str(int(token)) if 1 <= int(token) <= 22 else token)
    else:
        canonical = token
    if to_naming == "ncbi":
        return canonical
    if to_naming == "ucsc":
        return "chrM" if canonical == "MT" else f"chr{canonical}"
    return canonical


def write_random_vcf(
    path: Path,
    *,
    n_samples: int,
    n_snps: int,
    seed: int,
    missing_rate: float,
) -> List[str]:
    rng = random.Random(seed)
    samples = [f"ID{idx:04d}" for idx in range(1, n_samples + 1)]

    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)
        handle.write(header + "\n")

        for snp_idx in range(n_snps):
            pos = 100000 + snp_idx * 10
            ref = rng.choice(BASES)
            alt = rng.choice([base for base in BASES if base != ref])
            maf = rng.uniform(0.05, 0.5)
            row = ["1", str(pos), f"rs{snp_idx + 1}", ref, alt, ".", "PASS", ".", "GT"]
            for _ in range(n_samples):
                if rng.random() < missing_rate:
                    gt = "./."
                else:
                    dose = 0
                    if rng.random() < maf:
                        dose += 1
                    if rng.random() < maf:
                        dose += 1
                    gt = "0/0" if dose == 0 else "0/1" if dose == 1 else "1/1"
                row.append(gt)
            handle.write("\t".join(row) + "\n")

    return samples


def parse_vcf(path: Path) -> Tuple[List[str], List[VcfVariant]]:
    samples: List[str] = []
    variants: List[VcfVariant] = []
    with open(path, "r", encoding="utf-8", newline="\n") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            variants.append(
                VcfVariant(
                    chrom=parts[0],
                    pos=int(parts[1]),
                    vid=parts[2],
                    ref=parts[3],
                    alt=parts[4],
                    genotypes=[field.split(":", 1)[0] for field in parts[9:]],
                )
            )
    return samples, variants


def parse_vcf_genotypes(path: Path) -> Tuple[List[str], List[Tuple[str, str]], List[List[str]]]:
    samples: List[str] = []
    alleles: List[Tuple[str, str]] = []
    genotypes: List[List[str]] = []
    with open(path, "r", encoding="utf-8", newline="\n") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            alleles.append((parts[3], parts[4]))
            genotypes.append([field.split(":", 1)[0] for field in parts[9:]])
    return samples, alleles, genotypes


def gt_to_dosage(gt: str) -> Optional[int]:
    if gt in ("./.", ".|.", "."):
        return None
    sep = "/" if "/" in gt else "|"
    alleles = gt.split(sep)
    if len(alleles) != 2 or "." in alleles:
        return None
    return int(alleles[0]) + int(alleles[1])


def write_fam_sexes(fam_path: Path, sexes: List[int]) -> None:
    lines = []
    with open(fam_path, "r", encoding="utf-8", newline="\n") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            parts[4] = str(sexes[len(lines)])
            lines.append("\t".join(parts))
    write_lines(fam_path, lines)


def write_vcf(path: Path, samples: List[str], variants: List[HaploidVariant], genotypes: List[List[str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)
        handle.write(header + "\n")
        for variant, gts in zip(variants, genotypes):
            row = [
                variant.chrom,
                str(variant.pos),
                variant.snp,
                variant.a2,
                variant.a1,
                ".",
                "PASS",
                ".",
                "GT",
            ]
            row.extend(gts)
            handle.write("\t".join(row) + "\n")


def load_regions(schema_path: Path, build_name: str) -> Dict[str, dict]:
    with open(schema_path, "r", encoding="utf-8", newline="\n") as handle:
        data = json.load(handle)
    build = data["builds"][build_name]
    return {region["name"]: region for region in build["regions"]}


def in_interval(pos: int, interval: dict) -> bool:
    end = interval["end"] if interval["end"] is not None else 10**12
    return interval["start"] <= pos <= end


def classify_region(chrom: str, pos: int, regions: Dict[str, dict]) -> str:
    if chrom == "X":
        par1 = regions["PAR1"]["chromosomes"]["X"]
        par2 = regions["PAR2"]["chromosomes"]["X"]
        if in_interval(pos, par1):
            return "X_par1"
        if in_interval(pos, par2):
            return "X_par2"
        return "X_nonpar"
    if chrom == "Y":
        return "Y_msy"
    if chrom == "MT":
        return "MT"
    return "22"


def ploidy_for_region(region: str, sex: int) -> int:
    if region in ("X_par1", "X_par2"):
        return 2
    if region == "X_nonpar":
        return 1 if sex == 1 else 2
    if region == "Y_msy":
        return 1 if sex == 1 else 0
    if region == "MT":
        return 1
    return 2


def region_intervals(region: dict, chrom: str) -> List[Tuple[int, int]]:
    raw = region["chromosomes"][chrom]
    interval_list = raw if isinstance(raw, list) else [raw]
    out = []
    for interval in interval_list:
        start = int(interval["start"])
        end = int(interval["end"]) if interval["end"] is not None else start + 1_000_000
        out.append((start, end))
    return out


def generate_region_positions(region: dict, chrom: str, count: int) -> List[int]:
    positions: List[int] = []
    for start, end in region_intervals(region, chrom):
        pos = start
        while pos <= end and len(positions) < count:
            positions.append(pos)
            pos += 1
        if len(positions) >= count:
            break
    if len(positions) < count:
        raise RuntimeError(f"insufficient interval space for {chrom}: wanted {count}, found {len(positions)}")
    return positions


def load_haploid_variants(bim_path: Path, *, build_name: str, counts: Dict[str, int]) -> List[HaploidVariant]:
    schema_path = REPO_ROOT / "match" / "schemas" / "human_haploid_regions_grch37_grch38.json"
    regions = load_regions(schema_path, build_name)
    rows_by_chrom: Dict[str, List[object]] = {"22": [], "X": [], "Y": [], "MT": []}
    for row in read_bim(bim_path):
        if row.chrom in rows_by_chrom:
            rows_by_chrom[row.chrom].append(row)
    region_specs = {
        "22": ("22", None),
        "X_par1": ("X", "PAR1"),
        "X_par2": ("X", "PAR2"),
        "X_nonpar": ("X", "X_nonPAR"),
        "Y_msy": ("Y", "Y_MSY"),
        "MT": ("MT", "MT"),
    }
    row_offsets: Dict[str, int] = {chrom: 0 for chrom in rows_by_chrom}
    selected: List[HaploidVariant] = []
    for region_name, want in counts.items():
        chrom, schema_region_name = region_specs[region_name]
        chrom_rows = rows_by_chrom[chrom]
        start = row_offsets[chrom]
        stop = start + want
        if len(chrom_rows) < stop:
            raise RuntimeError(f"insufficient variants for {region_name}: {len(chrom_rows) - start} < {want}")
        row_offsets[chrom] = stop
        if schema_region_name is None:
            positions = list(range(1_000_001, 1_000_001 + want))
        else:
            positions = generate_region_positions(regions[schema_region_name], chrom, want)
        for row, pos in zip(chrom_rows[start:stop], positions):
            selected.append(HaploidVariant(chrom, row.snp, pos, row.a1, row.a2, region_name))
    return selected


def ensure_bim_import_contig_naming(vmap_path: Path, *, source_prefix: Path) -> None:
    meta_path = vmap_path.with_name(vmap_path.name + ".meta.json")
    metadata = json.loads(meta_path.read_text(encoding="utf-8"))
    target_meta = metadata["target"]
    if target_meta.get("contig_naming") is None:
        bim_rows = read_bim(source_prefix.with_suffix(".bim"))
        target_meta["contig_naming"] = "plink_splitx" if any(row.chrom == "25" for row in bim_rows) else "plink"
        meta_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_variant_object_target_metadata(path: Path) -> Dict[str, object]:
    meta_path = path.with_name(path.name + ".meta.json")
    metadata = json.loads(meta_path.read_text(encoding="utf-8"))
    if metadata["object_type"] == "variant_map":
        return dict(metadata["target"])
    return dict(metadata)


def build_source_vmap(source_prefix: Path, *, genome_build: str, output_path: Path) -> None:
    result = run_py(
        "import_bim.py",
        "--input",
        source_prefix.with_suffix(".bim"),
        "--output",
        output_path,
        "--genome-build",
        genome_build,
    )
    assert_py_ok(result)
    ensure_bim_import_contig_naming(output_path, source_prefix=source_prefix)


def build_source_vtable(source_prefix: Path, *, genome_build: str, output_path: Path) -> None:
    imported_vmap = output_path.with_name(output_path.stem + ".imported.vmap")
    build_source_vmap(source_prefix, genome_build=genome_build, output_path=imported_vmap)
    result = run_py(
        "convert_vmap_to_target.py",
        "--source",
        imported_vmap,
        "--output",
        output_path,
    )
    assert_py_ok(result)


def write_target_vtable(path: Path, rows: List[List[str]], *, genome_build: str, contig_naming: str = "ncbi") -> None:
    normalized_rows = []
    for row in rows:
        normalized_rows.append([normalize_test_chrom(row[0], to_naming=contig_naming), *row[1:]])
    write_lines(path, ["\t".join(row) for row in normalized_rows])
    write_json(
        path.with_name(path.name + ".meta.json"),
        {
            "object_type": "variant_table",
            "genome_build": genome_build,
            "contig_naming": contig_naming,
        },
    )


def build_vmap(source_vmap: Path, target_vtable: Path, output_path: Path) -> None:
    result = run_py(
        "match_vmap_to_target.py",
        "--source",
        source_vmap,
        "--target",
        target_vtable,
        "--output",
        output_path,
    )
    assert_py_ok(result)


def apply_vmap_to_bfile(source_prefix: Path, vmap_path: Path, output_prefix: Path):
    return run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap_path,
        "--output-prefix",
        output_prefix,
    )


def recode_bfile_to_vcf(plink_cmd: List[str], prefix: Path):
    result = run_cmd(
        plink_cmd
        + ["--bfile", str(prefix), "--recode", "vcf", "--keep-allele-order", "--out", str(prefix)]
    )
    assert_plink_ok(
        result,
        allow_warnings=True,
        allowed_warning_substrings=["het. haploid genotypes present"],
    )
    return result


def read_ploidy(path: Path) -> List[Tuple[int, int]]:
    return [(int(row[0]), int(row[1])) for row in read_tsv(path)]
