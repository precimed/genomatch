from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Tuple

HAPLOID_SCHEMA_PATH = Path(__file__).resolve().parent / "schemas" / "human_haploid_regions.json"
PLOIDY_MAP = {"diploid": 2, "haploid": 1, "absent": 0}
HAPLOID_CHROMS = {"X", "Y", "MT"}


@dataclass(frozen=True)
class RegionDef:
    name: str
    chrom: str
    intervals: List[Tuple[int, int | None]]
    ploidy_male: int
    ploidy_female: int


@lru_cache(maxsize=None)
def load_haploid_schema() -> Dict[str, object]:
    with open(HAPLOID_SCHEMA_PATH, "r", newline="\n") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"haploid schema must be a JSON object: {HAPLOID_SCHEMA_PATH}")
    return data


def x_par_intervals_for_build(genome_build: str) -> List[Tuple[int, int]]:
    schema = load_haploid_schema()
    builds = schema.get("builds")
    if not isinstance(builds, dict):
        raise ValueError(f"haploid schema is missing builds: {HAPLOID_SCHEMA_PATH}")
    build = builds.get(genome_build)
    if not isinstance(build, dict):
        raise ValueError(f"haploid schema is missing build {genome_build!r}: {HAPLOID_SCHEMA_PATH}")
    regions = build.get("regions")
    if not isinstance(regions, list):
        raise ValueError(f"haploid schema build {genome_build!r} is missing regions: {HAPLOID_SCHEMA_PATH}")
    intervals: List[Tuple[int, int]] = []
    for region in regions:
        if not isinstance(region, dict):
            continue
        if region.get("name") not in {"PAR1", "PAR2"}:
            continue
        chroms = region.get("chromosomes")
        if not isinstance(chroms, dict):
            continue
        x_region = chroms.get("X")
        if not isinstance(x_region, dict):
            continue
        start = x_region.get("start")
        end = x_region.get("end")
        if isinstance(start, int) and isinstance(end, int):
            intervals.append((start, end))
    if not intervals:
        raise ValueError(f"haploid schema build {genome_build!r} has no X PAR intervals: {HAPLOID_SCHEMA_PATH}")
    return intervals


def load_haploid_regions(schema_path: Path = HAPLOID_SCHEMA_PATH, build_name: str = "GRCh38") -> Dict[str, List[RegionDef]]:
    if schema_path != HAPLOID_SCHEMA_PATH:
        with open(schema_path, "r", newline="\n") as handle:
            data = json.load(handle)
        if not isinstance(data, dict):
            raise ValueError(f"haploid schema must be a JSON object: {schema_path}")
    else:
        data = load_haploid_schema()
    builds = data.get("builds")
    if not isinstance(builds, dict):
        raise ValueError(f"haploid schema is missing builds: {schema_path}")
    build = builds.get(build_name)
    if not isinstance(build, dict):
        raise ValueError(f"haploid schema missing {build_name} build definition")
    regions: Dict[str, List[RegionDef]] = {}
    raw_regions = build.get("regions")
    if not isinstance(raw_regions, list):
        raise ValueError(f"haploid schema build {build_name!r} is missing regions: {schema_path}")
    for region in raw_regions:
        if not isinstance(region, dict):
            continue
        ploidy = region.get("ploidy")
        if not isinstance(ploidy, dict):
            continue
        male = PLOIDY_MAP[str(ploidy["male"])]
        female = PLOIDY_MAP[str(ploidy["female"])]
        chroms = region.get("chromosomes")
        if not isinstance(chroms, dict):
            continue
        for chrom, intervals in chroms.items():
            interval_list = [intervals] if isinstance(intervals, dict) else intervals
            region_def = RegionDef(
                name=str(region["name"]),
                chrom=str(chrom),
                intervals=[(item["start"], item["end"]) for item in interval_list],
                ploidy_male=male,
                ploidy_female=female,
            )
            regions.setdefault(str(chrom), []).append(region_def)
    return regions


@lru_cache(maxsize=None)
def haploid_regions_for_build(build_name: str) -> Dict[str, List[RegionDef]]:
    return load_haploid_regions(HAPLOID_SCHEMA_PATH, build_name=build_name)


def expected_ploidy_pair(chrom: str, pos_raw: str, *, genome_build: str) -> Tuple[int, int]:
    from .vtable_utils import normalize_chrom_label

    canonical = normalize_chrom_label(chrom)
    if canonical not in HAPLOID_CHROMS:
        return (2, 2)
    try:
        pos = int(pos_raw)
    except ValueError as exc:
        raise ValueError(f"invalid position for ploidy model: {chrom}:{pos_raw}") from exc
    if pos <= 0:
        raise ValueError(f"invalid position for ploidy model: {chrom}:{pos_raw}")
    for region in haploid_regions_for_build(genome_build).get(canonical, []):
        for start, end in region.intervals:
            end_val = end if end is not None else 1_000_000_000_000
            if start <= pos <= end_val:
                return (region.ploidy_male, region.ploidy_female)
    raise ValueError(f"haploid region not found for {chrom}:{pos_raw}")


def has_non_diploid_ploidy(ploidy_pair: Tuple[int, int]) -> bool:
    return ploidy_pair != (2, 2)


def is_sex_dependent_ploidy(ploidy_pair: Tuple[int, int]) -> bool:
    return ploidy_pair[0] != ploidy_pair[1]
