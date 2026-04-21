#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
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
from .vtable_utils import normalize_allele_token, VariantRow, open_text

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert a VCF text file to .vmap.")
    parser.add_argument("--input", required=True, help="Input .vcf or .vcf.gz file")
    parser.add_argument("--output", required=True, help="Output .vmap file")
    parser.add_argument("--genome-build", default="unknown", help="Genome build for metadata")
    parser.add_argument("--chr2use", "--contigs", dest="chr2use", help="Comma-separated chromosome list or ranges")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_path = Path(args.output)
    logger.info("import_vcf.py: importing %s -> %s", args.input, output_path)
    rows: List[ImportedVariantRow] = []
    qc_rows: List[ImportQcRow] = []
    for shard in resolve_import_input_paths(args.input, kind_label="VCF"):
        source_index = 0
        with open_text(shard.path, "rt") as handle:
            found_header = False
            for line in handle:
                if not line.strip():
                    continue
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    found_header = True
                    continue
                if not found_header:
                    raise ValueError("VCF header line #CHROM not found")
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    qc_rows.append(ImportQcRow(shard.source_shard, source_index, "malformed_row"))
                    source_index += 1
                    continue
                chrom, pos, vid, ref, alt = parts[:5]
                if not chrom or not is_valid_import_position(pos):
                    qc_rows.append(ImportQcRow(shard.source_shard, source_index, "malformed_row"))
                    source_index += 1
                    continue
                if "," in alt:
                    qc_rows.append(ImportQcRow(shard.source_shard, source_index, "multiallelic"))
                    source_index += 1
                    continue
                ref = normalize_allele_token(ref)
                alt = normalize_allele_token(alt)
                if not is_canonical_import_allele(alt) or not is_canonical_import_allele(ref):
                    qc_rows.append(ImportQcRow(shard.source_shard, source_index, "non_actg_allele"))
                    source_index += 1
                    continue
                rows.append(
                    ImportedVariantRow(
                        VariantRow(chrom, pos, vid or ".", alt, ref),
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
        created_by="import_vcf.py",
        derived_from=Path(args.input),
        qc_rows=qc_rows,
    )
    logger.info("import_vcf.py: wrote %s with %s retained rows and %s QC rows", output_path, len(rows), len(qc_rows))
    return 0


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
