import pandas as pd
import pytest

from genomatch.apply_vmap_sumstats import output_variant_ids
from genomatch.apply_vmap_utils import output_variant_id
from genomatch.exact_set_utils import format_target_derived_ids
from genomatch.vmap_liftover_build import parse_lifted_vcf
from genomatch.vtable_utils import (
    VariantRow,
    format_variant_id,
    format_variant_ids,
    sort_target_rows_by_declared_coordinate,
    sort_target_table_by_declared_coordinate,
)
from utils import read_tsv, run_py, write_json, write_lines


def test_generated_ids_agree_across_scalar_vector_and_payload_renderers():
    # No metadata/reference lookup is needed, including for unverified rows.
    frame = pd.DataFrame([
        ('1', '100', 'keep-snv', 'T', 'G'),
        ('1', '100', 'keep-del', 'G', 'GA'),
        ('2', '9', 'keep-mnv', 'AC', 'TG'),
    ], columns=['chrom', 'pos', 'id', 'a1', 'a2'], index=[8, 3, 5])
    expected = ['1:100:G:T', '1:100:GA:G', '2:9:TG:AC']
    assert [format_variant_id(row) for row in frame.itertuples()] == expected
    assert format_variant_ids(frame).tolist() == expected
    assert format_target_derived_ids(frame).tolist() == expected
    assert output_variant_ids(frame, retain_snp_id=False).tolist() == expected
    assert [output_variant_id(row, retain_snp_id=False) for row in frame.itertuples()] == expected
    assert output_variant_ids(frame, retain_snp_id=True).tolist() == frame['id'].tolist()
    assert output_variant_ids(frame, retain_snp_id=False).index.tolist() == [0, 1, 2]
    assert [output_variant_id(row, retain_snp_id=True) for row in frame.itertuples()] == frame['id'].tolist()


def test_coordinate_sort_helpers_use_numeric_positions_and_ref_alt_ties():
    rows = [VariantRow(*row) for row in [
        ('10', '1', 'chr10', 'A', 'C'),
        ('2', '1', 'chr2', 'A', 'C'),
        ('1', '100', 'deletion', 'G', 'GA'),
        ('1', '100', 'snv-t', 'T', 'G'),
        ('1', '100', 'snv-c', 'C', 'G'),
        ('1', '20', 'early', 'A', 'C'),
        ('1', '100', 'same-key', 'T', 'G'),
    ]]
    expected = ['early', 'snv-c', 'snv-t', 'same-key', 'deletion', 'chr2', 'chr10']
    assert [row.id for row in sort_target_rows_by_declared_coordinate(rows, 'ncbi', label='test')] == expected
    frame = pd.DataFrame(rows)
    assert sort_target_table_by_declared_coordinate(frame, 'ncbi', label='test')['id'].tolist() == expected


@pytest.mark.parametrize('tool', ['sort_variants.py', 'intersect_variants.py', 'union_variants.py'])
def test_coordinate_order_for_colocated_snv_and_deletion(tmp_path, tool):
    source = tmp_path / 'source.vtable'
    out = tmp_path / 'out.vtable'
    write_lines(source, ['1\t100\tdel\tG\tGA', '1\t100\tsnv\tT\tG'])
    write_json(source.with_name(source.name + '.meta.json'), {
        'object_type': 'variant_table', 'genome_build': 'GRCh37',
        'contig_naming': 'ncbi', 'variants_count': 2,
    })
    args = ['--input', source] if tool == 'sort_variants.py' else [source, source]
    result = run_py(tool, *args, '--output', out)
    assert result.returncode == 0, result.stderr
    rows = read_tsv(out)
    assert [row[3:] for row in rows] == [['T', 'G'], ['G', 'GA']]
    assert [row[2] for row in rows] == (
        ['snv', 'del'] if tool == 'sort_variants.py' else ['1:100:G:T', '1:100:GA:G']
    )


def test_liftover_coordinate_tie_break_preserves_alleles_and_ids(tmp_path):
    lifted = tmp_path / 'lifted.vcf'
    write_lines(lifted, [
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
        'chr1\t100\t0\tGA\tG\t.\tPASS\t.',
        'chr1\t100\t1\tG\tT\t.\tPASS\t.',
    ])
    lookup = pd.DataFrame([
        ('0', '1', '100', 'del', 'identity', 0),
        ('1', '1', '100', 'snv', 'identity', 1),
    ], columns=['row_id', 'chrom', 'pos', 'id', 'upstream_allele_op', 'row_idx']).set_index('row_id', drop=False)
    table, statuses = parse_lifted_vcf(lifted, lookup, 'ncbi', 'GRCh37', 'GRCh38', is_vmap=False)
    assert table.to_frame(copy=False)[['id', 'a1', 'a2']].values.tolist() == [
        ['snv', 'T', 'G'], ['del', 'G', 'GA'],
    ]
    assert set(statuses.values()) == {'lifted'}
