#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

from ._cli_utils import run_cli
from .importer_utils import (
    filter_import_rows_by_chr2use,
    finalize_imported_vmap,
    ImportQcRow,
    ImportedVariantRow,
    is_canonical_import_allele,
    is_valid_import_position,
    resolve_import_input_paths,
)
from .vtable_utils import normalize_allele_token, VariantRow


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert a PLINK .bim file to .vmap.")
    parser.add_argument("--input", required=True, help="Input .bim file")
    parser.add_argument("--output", required=True, help="Output .vmap file")
    parser.add_argument("--genome-build", default="unknown", help="Genome build for metadata")
    parser.add_argument("--chr2use", "--contigs", dest="chr2use", help="Comma-separated chromosome list or ranges")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_path = Path(args.output)
    rows: List[ImportedVariantRow] = []
    qc_rows: List[ImportQcRow] = []
    for shard in resolve_import_input_paths(args.input, kind_label=".bim"):
        source_index = 0
        with open(shard.path, "r", encoding="utf-8", newline="\n") as handle:
            for line in handle:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split()
                if len(parts) < 6:
                    qc_rows.append(ImportQcRow(shard.source_shard, source_index, "malformed_row"))
                    source_index += 1
                    continue
                chrom, snp, _cm, bp, a1, a2 = parts[:6]
                a1 = normalize_allele_token(a1)
                a2 = normalize_allele_token(a2)
                if not chrom or not is_valid_import_position(bp):
                    qc_rows.append(ImportQcRow(shard.source_shard, source_index, "malformed_row"))
                    source_index += 1
                    continue
                if not is_canonical_import_allele(a1) or not is_canonical_import_allele(a2):
                    qc_rows.append(ImportQcRow(shard.source_shard, source_index, "non_actg_allele"))
                    source_index += 1
                    continue
                rows.append(
                    ImportedVariantRow(
                        VariantRow(chrom, bp, snp or ".", a1, a2),
                        shard.source_shard,
                        source_index,
                    )
                )
                source_index += 1
    rows, chr_qc_rows = filter_import_rows_by_chr2use(rows, args.chr2use)
    qc_rows.extend(chr_qc_rows)
    finalize_imported_vmap(
        output_path=output_path,
        rows=rows,
        genome_build=args.genome_build,
        created_by="import_bim.py",
        derived_from=Path(args.input),
        qc_rows=qc_rows,
    )
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
