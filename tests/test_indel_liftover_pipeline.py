import json
import random
from pathlib import Path

import pandas as pd
import pytest

from genomatch import vmap_liftover_build as liftover
from test_liftover_build import write_fake_bcftools
from test_indel_reference_assignment import write_object
from utils import (read_bed_matrix, read_tsv, run_py_with_env, write_bfile, write_json,
                   write_lines, write_match_config, write_primary_ucsc_fasta)


def pipeline_env(tmp_path, swap=False):
    source_fasta, target_fasta = tmp_path/'source.fa', tmp_path/'target.fa'
    write_primary_ucsc_fasta(source_fasta, {'chr1': 'AGCAGCAGC'})
    write_primary_ucsc_fasta(target_fasta, {'chr1': 'AGCAGCAGC'})
    chain = tmp_path/'chain'
    write_lines(chain, [f'chr1\t1\tchr1\t4' + ('\tswap=1' if swap else '')])
    config = tmp_path/'config.yaml'
    write_match_config(config, grch37_ucsc_fasta=source_fasta, grch38_ucsc_fasta=target_fasta, chain37to38=chain)
    binary = tmp_path/'bcftools'
    write_fake_bcftools(binary)
    return {'MATCH_CONFIG': str(config), 'MATCH_BCFTOOLS': str(binary)}


@pytest.mark.parametrize('swap', [False, True])
def test_deletion_restriction_liftover_and_genotype_projection(tmp_path, swap):
    env = pipeline_env(tmp_path, swap)
    source, restricted, lifted = [tmp_path/name for name in ('source.vmap', 'restricted.vmap', 'lifted.vmap')]
    # Already prepared deletion, with a prior swap from the raw payload order.
    write_object(source, ['1\t1\tdeletion\tA\tAG'], 'vmap')
    # write_object uses a nontrivial raw locator; make this single-file payload row 0.
    write_lines(source, ['1\t1\tdeletion\tA\tAG\t.\t0\tswap'])
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', source, '--output', restricted, '--include-indels')
    assert result.returncode == 0, result.stderr
    assert read_tsv(restricted) == [['1', '1', 'deletion', 'A', 'AG', '.', '0', 'swap']]
    result = run_py_with_env('liftover_build.py', env, '--input', restricted, '--output', lifted, '--target-build', 'GRCh38')
    assert result.returncode == 0, result.stderr
    input_records = read_tsv(lifted.with_name(lifted.name+'.bcftools_input.vcf'))
    record = next(r for r in input_records if not r[0].startswith('#'))
    assert record[3:5] == ['AG', 'A']
    expected_op = 'identity' if swap else 'swap'
    assert read_tsv(lifted) == [['1', '4', 'deletion', *(['AG', 'A'] if swap else ['A', 'AG']), '.', '0', expected_op]]
    raw, projected = tmp_path/'raw', tmp_path/'projected'
    write_bfile(raw, ['1\tdeletion\t0\t1\tAG\tA'], ['F\tS1\t0\t0\t1\t-9', 'F\tS2\t0\t0\t2\t-9'], [[0, 3]])
    result = run_py_with_env('apply_vmap_to_bfile.py', env, '--vmap', lifted, '--source-prefix', raw, '--output-prefix', projected)
    assert result.returncode == 0, result.stderr
    assert read_bed_matrix(projected.with_suffix('.bed'), 2, 1) == ([[0, 3]] if swap else [[3, 0]])
    study, metadata, effects = tmp_path/'study.tsv', tmp_path/'study.yaml', tmp_path/'effects.tsv'
    write_lines(study, ['CHR\tPOS\tSNP\tEA\tOA\tBETA', '1\t1\tdeletion\tAG\tA\t0.5'])
    write_lines(metadata, ['col_CHR: CHR', 'col_POS: POS', 'col_SNP: SNP', 'col_EffectAllele: EA', 'col_OtherAllele: OA', 'col_BETA: BETA'])
    result = run_py_with_env('apply_vmap_to_sumstats.py', env, '--vmap', lifted, '--input', study, '--sumstats-metadata', metadata, '--output', effects)
    assert result.returncode == 0, result.stderr
    result_table = pd.read_table(effects)
    assert result_table['BETA'].tolist() == ([0.5] if swap else [-0.5])


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
@pytest.mark.parametrize('ref,alt,pos,info,status', [
    ('AG', 'CC', '2', '.', None), ('AG', 'A', '2', 'FLIP', None),
    ('A', 'C,G', '2', '.', 'unsupported_representation'),
    ('A', '<DEL>', '2', '.', 'unsupported_representation'),
    ('A', 'A', '2', '.', 'unsupported_representation'),
    ('A', 'C', '0', '.', 'unsupported_representation'),
    ('A', 'C', 'bad', '.', 'unsupported_representation'),
    ('A', 'C', '2', 'SWAP=-1', 'unsupported_representation'),
    ('A', 'C', '2', 'SWAP=2', 'unsupported_representation'),
    ('A', 'C', '2', 'FLIP=1', 'unsupported_representation'),
])
def test_lifted_output_acceptance(tmp_path, suffix, ref, alt, pos, info, status):
    is_vmap = suffix == 'vmap'
    lookup = pd.DataFrame({'row_id': ['row0'], 'chrom': ['chr1'], 'pos': ['1'], 'id': ['variant'],
                           'upstream_allele_op': ['swap'],
                           **({'source_shard': ['original'], 'source_index': [91]} if is_vmap else {'row_idx': [7]})}).set_index('row_id', drop=False)
    vcf = tmp_path/'lifted.vcf'
    write_lines(vcf, [f'chr1\t{pos}\trow0\t{ref}\t{alt}\t.\tPASS\t{info}'])
    table, statuses = liftover.parse_lifted_vcf(vcf, lookup, 'ncbi', 'GRCh37', 'GRCh38', is_vmap=is_vmap)
    assert statuses == {'row0': status or 'lifted'}
    if status:
        assert table.frame.empty
    else:
        assert table.frame[['a1', 'a2']].values.tolist() == [[alt, ref]]
        assert table.frame['source_index' if is_vmap else 'row_idx'].tolist() == [91 if is_vmap else 7]
        if is_vmap:
            assert table.frame['allele_op'].tolist() == ['flip_swap' if info == 'FLIP' else 'swap']


def test_multiple_results_and_invalid_candidate_do_not_poison_duplicate_detection(tmp_path):
    lookup = pd.DataFrame({'row_id': ['row0', 'row1', 'row2'], 'chrom': ['chr1']*3, 'pos': ['1', '2', '3'], 'id': ['multiple', 'bad', 'good'], 'upstream_allele_op': ['identity']*3, 'row_idx': [8, 2, 9]}).set_index('row_id', drop=False)
    path = tmp_path/'output.vcf'
    write_lines(path, ['chr1\t5\trow0\tA\tC\t.\tPASS\t.', 'chr1\t6\trow0\tA\tG\t.\tPASS\t.', 'chr1\t5\trow1\tA\tC\t.\tPASS\tSWAP=2', 'chr1\t5\trow2\tA\tC\t.\tPASS\t.'])
    table, statuses = liftover.parse_lifted_vcf(path, lookup, 'ncbi', 'GRCh37', 'GRCh38', is_vmap=False)
    assert table.frame['row_idx'].tolist() == [9]
    assert statuses == {'row0': 'unsupported_representation', 'row1': 'unsupported_representation', 'row2': 'lifted'}


@pytest.mark.parametrize('reverse', [False, True])
def test_real_plugin_multibase_and_reverse_strand(tmp_path, reverse):
    from real_liftover_helpers import resolve_bcftools_with_liftover
    binary = resolve_bcftools_with_liftover()
    rng = random.Random(71)
    sequence = ''.join(rng.choice('ACGT') for _ in range(600))
    source_fasta, target_fasta = tmp_path/'source.fa', tmp_path/'target.fa'
    target_sequence = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1] if reverse else sequence
    write_primary_ucsc_fasta(source_fasta, {'chr1': sequence})
    write_primary_ucsc_fasta(target_fasta, {'chr1': target_sequence})
    chain = tmp_path/'chain'
    write_lines(chain, [f'chain 1000 chr1 600 + 0 600 chr1 600 {"-" if reverse else "+"} 0 600 1', '600', ''])
    config = tmp_path/'config.yaml'
    write_match_config(config, grch37_ucsc_fasta=source_fasta, grch38_ucsc_fasta=target_fasta, chain37to38=chain)
    source, output = tmp_path/'source.vmap', tmp_path/'lifted.vmap'
    # Deletion, MNV, delins; entire source a2 is reference-compatible.
    refs = [sequence[p-1:p+2] for p in (101, 201, 301)]
    alts = [refs[0][0], refs[1].translate(str.maketrans('ACGT', 'CGTA')), refs[2][0]+'CCCCC']
    write_object(source, [f'1\t{p}\tvariant{i}\t{alt}\t{ref}' for i, (p, alt, ref) in enumerate(zip((101, 201, 301), alts, refs))], 'vmap')
    result = run_py_with_env('liftover_build.py', {'MATCH_CONFIG': str(config), 'MATCH_BCFTOOLS': binary}, '--input', source, '--output', output, '--target-build', 'GRCh38')
    assert result.returncode == 0, result.stderr
    rows = read_tsv(output)
    assert len(rows) == 3
    assert any(len(row[3]) > 1 and len(row[4]) > 1 for row in rows)
    assert {r[-1] for r in rows} == ({'flip_swap'} if reverse else {'swap'})
    emitted_input = [r for r in read_tsv(output.with_name(output.name+'.bcftools_input.vcf')) if not r[0].startswith('#')]
    assert [r[3:5] for r in emitted_input] == [[ref, alt] for ref, alt in zip(refs, alts)]
    for row in rows:
        pos = int(row[1])
        assert target_sequence[pos-1:pos-1+len(row[4])] == row[4]


@pytest.mark.parametrize('source_chrom,source_pos,target_chrom,info,expected', [
    ('chrY', '1', 'chr1', '.', 'ploidy_class_changed'),
    ('chrMT', '1', 'chr1', '.', 'ploidy_class_changed'),
    ('chrX', '60001', 'chr1', '.', 'lifted'),  # source PAR and target autosome are both (2,2)
    ('chrY', '1', 'chr1', 'SWAP=2', 'unsupported_representation'),
    ('chr1', '1', 'chrUn', 'SWAP=2', 'unsupported_representation'),
])
def test_ploidy_and_independent_annotation_exclusions(tmp_path, source_chrom, source_pos, target_chrom, info, expected):
    lookup = pd.DataFrame({'row_id': ['row0'], 'chrom': [source_chrom], 'pos': [source_pos], 'id': ['variant'], 'upstream_allele_op': ['identity'], 'row_idx': [11]}).set_index('row_id', drop=False)
    path = tmp_path/'lifted.vcf'
    write_lines(path, [f'{target_chrom}\t5\trow0\tAG\tCC\t.\tPASS\t{info}'])
    table, statuses = liftover.parse_lifted_vcf(path, lookup, 'ncbi', 'GRCh37', 'GRCh38', is_vmap=False)
    assert statuses == {'row0': expected}
    assert len(table.frame) == (1 if expected == 'lifted' else 0)


@pytest.mark.parametrize('a1,a2', [('A', 'AT'), ('CC', 'AG'), ('AG', 'CC')])
def test_source_vcf_preparation_never_reselects_alleles(a1, a2):
    frame = pd.DataFrame({'row_id': ['stable'], 'chrom': ['chr1'], 'pos': ['1'], 'id': ['variant'], 'a1': [a1], 'a2': [a2], 'upstream_allele_op': ['swap'], 'row_idx': [17]})
    prepared = liftover.prepare_liftover_table(frame, is_vmap=False)
    assert prepared[['ref', 'alt', 'status', 'row_idx']].values.tolist() == [[a2, a1, 'ready', 17]]
