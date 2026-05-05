from __future__ import annotations

from typing import Dict, List, Optional

from .haploid_utils import load_haploid_schema, x_par_intervals_for_build


SUPPORTED_CONTIG_NAMINGS = {"ncbi", "ucsc", "plink", "plink_splitx"}
NUMERIC_CONTIG_NAMINGS = {"ncbi", "plink", "plink_splitx"}
UNKNOWN_CONTIG = "unknown"

CHROM_SYNONYMS = {
    "x": "X",
    "xy": "X",
    "y": "Y",
    "m": "MT",
    "mt": "MT",
    "par": "X",
    "par1": "X",
    "par2": "X",
    "x_par1": "X",
    "x_par2": "X",
    "x_nonpar": "X",
    "nonpar": "X",
}
PLINK_NUMERIC_MAP = {
    23: "X",
    24: "Y",
    25: "X",
    26: "MT",
}
CANONICAL_CONTIG_ORDER = [str(idx) for idx in range(1, 23)] + ["X", "Y", "MT"]


def supported_exact_contig_tokens() -> Dict[str, str]:
    tokens: Dict[str, str] = {}
    for canonical in CANONICAL_CONTIG_ORDER:
        tokens[canonical] = canonical
        if canonical == "MT":
            tokens["chrM"] = canonical
            tokens["26"] = canonical
            continue
        tokens[f"chr{canonical}"] = canonical
        if canonical == "X":
            tokens["23"] = canonical
            tokens["25"] = canonical
            tokens["XY"] = canonical
            tokens["chrXY"] = canonical
        elif canonical == "Y":
            tokens["24"] = canonical
    return tokens


def canonical_contig_from_label(
    label: str,
    naming: str,
    *,
    allow_unknown: bool = False,
) -> Optional[str]:
    token = label.strip()
    if not token:
        return None
    if token == UNKNOWN_CONTIG:
        return UNKNOWN_CONTIG if allow_unknown else None
    if naming == "ncbi":
        if token in {"X", "Y", "MT"}:
            return token
        if token.isdigit():
            value = int(token)
            if 1 <= value <= 22:
                return str(value)
        return None
    if naming == "ucsc":
        if token == "chrM":
            return "MT"
        if token.startswith("chr"):
            raw = token[3:]
            if raw in {"X", "Y"}:
                return raw
            if raw.isdigit():
                value = int(raw)
                if 1 <= value <= 22:
                    return str(value)
        return None
    if naming == "plink":
        if not token.isdigit():
            return None
        value = int(token)
        if 1 <= value <= 22:
            return str(value)
        if value == 23:
            return "X"
        if value == 24:
            return "Y"
        if value == 26:
            return "MT"
        return None
    if naming == "plink_splitx":
        if not token.isdigit():
            return None
        value = int(token)
        if 1 <= value <= 22:
            return str(value)
        if value in {23, 25}:
            return "X"
        if value == 24:
            return "Y"
        if value == 26:
            return "MT"
        return None
    raise ValueError(f"unsupported contig naming: {naming}")


def normalize_chrom_label(raw: Optional[str]) -> str:
    if raw is None:
        return ""
    tmp = raw.strip()
    if not tmp:
        return ""
    tmp_lower = tmp.lower()
    if tmp_lower.startswith("chr"):
        tmp_lower = tmp_lower[3:]
    if tmp_lower in CHROM_SYNONYMS:
        return CHROM_SYNONYMS[tmp_lower]
    if tmp_lower.isdigit():
        value = int(tmp_lower)
        if 1 <= value <= 22:
            return str(value)
        return PLINK_NUMERIC_MAP.get(value, "")
    if tmp in {str(idx) for idx in range(1, 23)} | {"X", "Y", "MT"}:
        return tmp
    return ""


def canonical_contig_from_any_supported_label(label: str) -> Optional[str]:
    canonical_values = {
        canonical
        for naming in SUPPORTED_CONTIG_NAMINGS
        for canonical in [canonical_contig_from_label(label, naming)]
        if canonical not in {None, UNKNOWN_CONTIG}
    }
    if len(canonical_values) == 1:
        return next(iter(canonical_values))
    if not canonical_values:
        return None
    raise ValueError(f"ambiguous contig label across supported namings: {label!r}")


def contig_label_for_naming(canonical: str, naming: str) -> str:
    if naming == "ncbi":
        return canonical
    if naming == "ucsc":
        if canonical == "MT":
            return "chrM"
        return f"chr{canonical}"
    if naming == "plink":
        reverse = {"X": "23", "Y": "24", "MT": "26"}
        return reverse.get(canonical, canonical)
    if naming == "plink_splitx":
        reverse = {"X": "23", "Y": "24", "MT": "26"}
        return reverse.get(canonical, canonical)
    raise ValueError(f"unsupported contig naming: {naming}")


def convert_contig_label(label: str, from_naming: str, to_naming: str) -> str:
    if from_naming == to_naming:
        return label
    canonical = canonical_contig_from_label(label, from_naming)
    if not canonical:
        raise ValueError(f"unable to normalize contig label: {label!r}")
    return contig_label_for_naming(canonical, to_naming)


def resolve_plink_splitx_build(genome_build: str | None) -> str:
    if genome_build in {"GRCh37", "GRCh38", "T2T-CHM13v2.0"}:
        return str(genome_build)
    if genome_build in {None, "", "unknown"}:
        raise ValueError(
            "normalize_contigs.py --to plink_splitx requires metadata genome_build=GRCh37, GRCh38, or T2T-CHM13v2.0"
        )
    raise ValueError(
        f"normalize_contigs.py --to plink_splitx does not support genome_build={genome_build!r}; "
        "expected GRCh37, GRCh38, or T2T-CHM13v2.0"
    )


def contig_label_for_plink_splitx(canonical: str, pos: str, genome_build: str | None) -> str:
    if canonical != "X":
        return contig_label_for_naming(canonical, "plink_splitx")
    build_name = resolve_plink_splitx_build(genome_build)
    try:
        pos_value = int(pos)
    except ValueError as exc:
        raise ValueError(f"normalize_contigs.py --to plink_splitx requires integer pos for X rows; found {pos!r}") from exc
    if pos_value <= 0:
        raise ValueError(f"normalize_contigs.py --to plink_splitx requires positive pos for X rows; found {pos!r}")
    xy_par_code = str(load_haploid_schema()["plink_chromosome_codes_plink1"]["XY_PAR"])
    for start, end in x_par_intervals_for_build(build_name):
        if start <= pos_value <= end:
            return xy_par_code
    return "23"


def repair_contig_label(label: str, to_naming: str, *, pos: str | None = None, genome_build: str | None = None) -> str:
    canonical = normalize_chrom_label(label)
    if not canonical:
        return UNKNOWN_CONTIG
    if to_naming == "plink_splitx":
        if pos is None:
            raise ValueError("normalize_contigs.py --to plink_splitx requires row positions")
        return contig_label_for_plink_splitx(canonical, pos, genome_build)
    return contig_label_for_naming(canonical, to_naming)


def normalize_contig_for_reference(label: str, declared_naming: str, reference_naming: str = "ucsc") -> str:
    canonical = canonical_contig_from_label(label, declared_naming)
    if not canonical:
        raise ValueError(
            f"unable to normalize contig label {label!r} from declared contig_naming={declared_naming!r} "
            f"to internal reference naming {reference_naming!r}"
        )
    return contig_label_for_naming(canonical, reference_naming)
