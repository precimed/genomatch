#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from ._cli_utils import run_cli
from .reference_utils import fetch_reference_bases, resolve_internal_reference_fasta
from .vtable_utils import (
    VariantRow,
    load_variant_object,
    normalize_contig_for_reference,
    require_contig_naming,
    require_rows_match_contig_naming,
    write_metadata,
)

MIN_GUESS_COMPATIBILITY_RATE = 0.60
MIN_GUESS_MARGIN = 0.10
HIGH_CONFIDENCE_COMPATIBILITY_RATE = 0.85


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Guess target-side genome build metadata for a .vtable or .vmap using UCSC-style internal reference FASTA assets."
    )
    parser.add_argument("--input", required=True, help="Input .vtable or .vmap")
    parser.add_argument("--write", action="store_true", help="Write updated metadata sidecar")
    parser.add_argument("--force", action="store_true", help="Overwrite existing metadata fields")
    return parser.parse_args()


def compatibility_rate(rows: List[VariantRow], fasta_path: Path, contig_naming: str) -> Tuple[int, int]:
    candidate_rows: List[Tuple[VariantRow, str, int]] = []
    queries: List[Tuple[str, int]] = []
    for row in rows:
        try:
            pos = int(row.pos)
        except ValueError:
            continue
        try:
            ucsc_contig = normalize_contig_for_reference(row.chrom, contig_naming, "ucsc")
        except ValueError:
            continue
        candidate_rows.append((row, ucsc_contig, pos))
        queries.append((ucsc_contig, pos))
    reference_bases = fetch_reference_bases(fasta_path, queries)

    compatible = 0
    considered = 0
    for row, ucsc_contig, pos in candidate_rows:
        ref_base = reference_bases.get((ucsc_contig, pos), "")
        if ref_base not in {"A", "C", "G", "T"}:
            continue
        considered += 1
        alleles = {row.a1.upper(), row.a2.upper()}
        if ref_base in alleles:
            compatible += 1
    return compatible, considered


def guess_build(rows: List[VariantRow], contig_naming: str) -> Dict[str, object]:
    evidence = []
    build_names = ["GRCh37", "GRCh38"]
    if not build_names:
        raise ValueError("config is missing UCSC reference entries for GRCh37/GRCh38")
    for build in build_names:
        fasta_path = resolve_internal_reference_fasta(build)
        compatible, considered = compatibility_rate(rows, fasta_path, contig_naming)
        rate = compatible / considered if considered else 0.0
        evidence.append(
            {
                "genome_build": build,
                "fasta": str(fasta_path),
                "internal_reference_naming": "ucsc",
                "normalization": "none" if contig_naming == "ucsc" else f"{contig_naming}->ucsc",
                "compatible_rows": compatible,
                "considered_rows": considered,
                "compatibility_rate": rate,
            }
        )
    evidence.sort(key=lambda item: float(item["compatibility_rate"]), reverse=True)
    best = evidence[0] if evidence else None
    second_rate = float(evidence[1]["compatibility_rate"]) if len(evidence) > 1 else 0.0
    guessed = "unknown"
    confidence = "low"
    if best is not None:
        best_rate = float(best["compatibility_rate"])
        if best_rate >= MIN_GUESS_COMPATIBILITY_RATE and (best_rate - second_rate) >= MIN_GUESS_MARGIN:
            guessed = str(best["genome_build"])
            confidence = "high" if best_rate >= HIGH_CONFIDENCE_COMPATIBILITY_RATE else "medium"
    return {
        "genome_build": guessed,
        "confidence": confidence,
        "evidence": evidence,
    }


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.exists():
        raise ValueError(f"input variant object not found: {input_path}")

    loaded = load_variant_object(input_path)
    current_meta = dict(loaded.raw_metadata)
    contig_naming = require_contig_naming(loaded.target_metadata, label="variant object")
    rows = loaded.target_rows
    if not rows:
        raise ValueError("empty variant object")
    require_rows_match_contig_naming(rows, contig_naming, label="variant object")
    build_guess = guess_build(rows, contig_naming)

    output_meta = dict(current_meta)
    if loaded.object_type == "variant_map":
        target_meta = dict(output_meta["target"])
        metadata_can_change = args.force or target_meta.get("genome_build") == "unknown"
    else:
        metadata_can_change = args.force or current_meta.get("genome_build") == "unknown"
    if metadata_can_change:
        if loaded.object_type == "variant_map":
            output_meta["target"] = dict(output_meta["target"])
            output_meta["target"]["genome_build"] = build_guess["genome_build"]
        else:
            output_meta["genome_build"] = build_guess["genome_build"]

    summary = {
        "input": str(input_path),
        "object_type": loaded.object_type,
        "contig_naming": contig_naming,
        "internal_reference_naming": "ucsc",
        "normalization": "none" if contig_naming == "ucsc" else f"{contig_naming}->ucsc",
        "genome_build": build_guess["genome_build"],
        "confidence": build_guess["confidence"],
        "write_requested": bool(args.write),
        "metadata_updated": False,
        "evidence": build_guess["evidence"],
    }
    if args.write and metadata_can_change:
        write_metadata(input_path, output_meta)
        summary["metadata_updated"] = True

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
