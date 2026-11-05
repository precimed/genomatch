#!/usr/bin/env python3
"""Generate schema 1.0 objects with the initial schema 1.0 writer.

Run from the repository root in the match-liftover environment:
    python tests/fixtures/schema/generate_1_0.py

Use the writer revision recorded in README.md. These are compatibility anchors,
not routine regenerated test inputs. This generator refuses a different current
writer schema; it must not be adapted to relabel a newer format as schema 1.0.
"""
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / 'src'))

from genomatch.vtable_utils import (
    CURRENT_SCHEMA_VERSION, VMapRow, VariantRow, make_vmap_metadata,
    make_vtable_metadata, metadata_with_variants_count, write_metadata,
    write_vmap, write_vtable,
)


def main():
    if CURRENT_SCHEMA_VERSION != '1.0':
        raise RuntimeError('Use the original schema 1.0 writer checkout; see README.md')
    directory = Path(__file__).resolve().parent
    rows = [
        VMapRow('2', '20', 'rs_snv', 'T', 'C', 'chr2', 7, 'identity'),
        VMapRow('1', '100', '1:100:AC:A', 'A', 'AC', 'chr1', 2, 'swap'),
        VMapRow('1', '100', 'rs_multibase', 'GT', 'ACG', 'chr1', 8, 'flip_swap'),
    ]
    provenance = {'created_by': 'generate_1_0.py', 'derived_from': 'schema-fixture'}
    vmap = directory / 'vmap_1_0.vmap'
    write_vmap(vmap, rows)
    write_metadata(vmap, metadata_with_variants_count(make_vmap_metadata(
        {'genome_build': 'GRCh37', 'contig_naming': 'ncbi'}, provenance,
    ), len(rows)))
    vtable = directory / 'vtable_1_0.vtable'
    write_vtable(vtable, [VariantRow(row.chrom, row.pos, row.id, row.a1, row.a2) for row in rows])
    write_metadata(vtable, metadata_with_variants_count(make_vtable_metadata(
        'GRCh37', 'ncbi', provenance,
    ), len(rows)))


if __name__ == '__main__':
    main()
