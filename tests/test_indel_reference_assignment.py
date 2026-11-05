"""Ordered reference assignment and normalization boundary regressions."""
import itertools
import json
from pathlib import Path

import pandas as pd
import pytest

from genomatch import vmap_restrict_build_compatible as restriction
from genomatch.tabular_rows import VariantRowsTable
from test_restrict_build_compatible import restrict_env_with_fake_norm
from utils import read_tsv, run_py_with_env, write_object_metadata, write_lines


def write_object(path, rows, suffix):
    write_lines(path, [f'{row}\toriginal\t{100 + i}\tswap' if suffix == 'vmap' else row for i, row in enumerate(rows)])
    target = {'genome_build': 'GRCh37', 'contig_naming': 'ncbi'}
    write_object_metadata(path.with_name(path.name + '.meta.json'), {'object_type': 'variant_map', 'target': target} if suffix == 'vmap' else {'object_type': 'variant_table', **target})


@pytest.mark.parametrize('include,swaps,flips,norm', itertools.product([False, True], repeat=4))
def test_all_input_classes_and_independent_permissions(tmp_path, monkeypatch, include, swaps, flips, norm):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={}, grch37_sequence='AGCAAG')
    fasta = tmp_path / 'GRCh37.ucsc.fa'
    # Each class includes supplied REF, only-a1, partial REF, and neither matches.
    rows = [
        ('1', '1', 'snv', 'C', 'A'), ('1', '1', 'snv_swap', 'A', 'C'),
        ('1', '1', 'snv_flip', 'G', 'T'), ('1', '1', 'snv_flip_swap', 'T', 'G'),
        ('1', '1', 'deletion', 'A', 'AG'), ('1', '1', 'partial', 'A', 'AT'),
        ('1', '1', 'neither', 'T', 'CC'), ('1', '1', 'both', 'AGC', 'AG'),
        ('1', '1', 'both_swap', 'AG', 'AC'), ('1', '1', 'both_neither', 'CC', 'TT'),
        ('1', '4', 'no_multibase_flip', 'A', 'TT'),
        ('1', '1', 'identical_snv', 'A', 'A'), ('1', '1', 'identical_multi', 'AG', 'AG'),
    ]
    calls = []
    def normalize(candidates, **kwargs):
        calls.extend(candidates)
        return {c.candidate_id: restriction.NormalizationOutcome(c.row, c.local_op) for c in candidates}
    monkeypatch.setattr(restriction, 'normalize_candidates', normalize)
    table = VariantRowsTable.from_frame(pd.DataFrame(rows, columns=['chrom', 'pos', 'id', 'a1', 'a2']))
    outcome = restriction.build_restriction_outcomes_table(table, fasta_path=fasta, contig_naming='ncbi', allow_strand_flips=flips, norm_indels=norm, genome_build='GRCh37', output_path=tmp_path/'out', include_indels=include, allow_indel_swaps=swaps)
    expected = [True, True, flips, flips, include, include and swaps, False, include, include and swaps, False, include and swaps, False, False]
    assert outcome['keep'].tolist() == expected
    assert outcome.loc[:3, 'local_op'].tolist() == ['identity', 'swap', 'flip' if flips else 'identity', 'flip_swap' if flips else 'identity']
    for i in range(4, 11):
        assert outcome.at[i, 'local_op'] == ('swap' if expected[i] and i in (5, 8, 10) else 'identity')
    assert outcome.loc[11:, 'status'].tolist() == ['identical_alleles'] * 2
    if not include:
        assert outcome.loc[4:10, 'status'].tolist() == ['non_snv_input'] * 7
    assert len(calls) == (sum(expected[4:11]) if norm else 0)
    assert all(c.row.a2 in ('A', 'AG') for c in calls)


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
@pytest.mark.parametrize('mode', ['BULK', 'LEGACY'])
@pytest.mark.parametrize('norm', [False, True])
def test_deletion_and_qc_keep_original_locators(tmp_path, suffix, mode, norm):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={}, grch37_sequence='AGCATC')
    env['MATCH_REFERENCE_ACCESS_MODE'] = mode
    source, output = tmp_path/f'in.{suffix}', tmp_path/f'out.{suffix}'
    write_object(source, ['1\t1\tdeletion\tA\tAG', '1\t4\tpartial\tA\tAG', '1\t1\tidentical\tag\tAG'], suffix)
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', source, '--output', output, '--include-indels', *(['--norm-indels'] if norm else []))
    assert result.returncode == 0, result.stderr
    assert read_tsv(output)[0][:5] == ['1', '1', 'deletion', 'A', 'AG']
    if suffix == 'vmap':
        assert read_tsv(output)[0][5:] == ['original', '100', 'swap']
    qc = read_tsv(output.with_name(output.name+'.qc.tsv'))
    assert [r[-1] for r in qc[1:]] == ['no_reference_match', 'identical_alleles']
    assert [r[:2] if suffix == 'vmap' else r[:1] for r in qc[1:]] == ([['original', '101'], ['original', '102']] if suffix == 'vmap' else [['1'], ['2']])


@pytest.mark.parametrize('norm', [False, True])
def test_cleanup_and_noop_warning_without_inclusion(tmp_path, norm):
    env, log = restrict_env_with_fake_norm(tmp_path, behavior={})
    source, output = tmp_path/'in.vtable', tmp_path/'out.vtable'
    write_object(source, ['1\t1\tkeep\tC\tA', '1\t1\texcluded\tA\tAA'], 'vtable')
    for suffix in ('_input_e.vcf', '_output_x.log', '_batch99.input.vcf'):
        output.with_name(output.name+'.bcftools_norm'+suffix).write_text('stale')
    unrelated = output.with_name(output.name+'.other')
    unrelated.write_text('keep')
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', source, '--output', output, '--allow-indel-swaps', *(['--norm-indels'] if norm else []))
    assert result.returncode == 0, result.stderr
    assert 'no-ops' in result.stderr
    assert not list(tmp_path.glob('out.vtable.bcftools_norm*'))
    assert unrelated.read_text() == 'keep'
    assert not log.exists()
    assert read_tsv(output.with_name(output.name+'.qc.tsv'))[-1] == ['1', 'excluded', 'non_snv_input']


@pytest.mark.parametrize('failure', ['missing', 'unknown', 'malformed', 'ref_mismatch', 'multiallelic', 'identical', 'position', 'contig', 'multiple', 'symbolic'])
@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
def test_normalization_validation_and_no_assignment_fallback(tmp_path, failure, suffix):
    record = {'pos': 1, 'ref': 'AG', 'alt': 'A'}
    behavior = {'records': [record]}
    expected = None
    if failure == 'missing':
        behavior, expected = {'action': 'drop'}, 'norm_ref_mismatch'
    elif failure == 'unknown':
        record['id'] = 'unknown'
    elif failure == 'malformed':
        record['chrom'] = 'chr1\nbroken'
    elif failure == 'ref_mismatch':
        record['ref'] = 'AT'
    elif failure == 'multiallelic':
        record['alt'], expected = 'A,C', 'norm_multiallelic'
    elif failure == 'identical':
        record['alt'], expected = 'AG', 'norm_identical_ref_alt_alleles'
    elif failure == 'position':
        record['pos'], expected = 0, 'norm_invalid_position'
    elif failure == 'contig':
        record['chrom'], expected = 'chrUn', 'unsupported_target_contig'
    elif failure == 'multiple':
        behavior, expected = {'records': [record, record]}, 'norm_multiple_output_records'
    elif failure == 'symbolic':
        record['alt'], expected = '<DEL>', 'norm_not_atcg_alleles'
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={'row0': behavior}, grch37_sequence='AGC')
    source, output = tmp_path/f'in.{suffix}', tmp_path/f'out.{suffix}'
    # Both a2 and a1 match; the selected identity must never fall back to swap.
    write_object(source, ['1\t1\tboth\tAGC\tAG'], suffix)
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', source, '--output', output, '--include-indels', '--allow-indel-swaps', '--norm-indels')
    if expected is None:
        assert result.returncode != 0
    else:
        assert result.returncode == 0, result.stderr
        assert read_tsv(output) == []
        assert read_tsv(output.with_name(output.name+'.qc.tsv'))[-1][-1] == expected
    input_vcf = next(tmp_path.glob('out.*.bcftools_norm*.input.vcf')).read_text()
    assert '\tAG\tAGC\t' in input_vcf
    assert input_vcf.count('\trow0\t') == 1


def test_batch_round_trips_preserve_ids_and_cleanup_fewer_batches(tmp_path, monkeypatch):
    import sys
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={'row2': {'action': 'drop'}}, grch37_sequence='AGCAGCAGC')
    source, output = tmp_path/'in.vtable', tmp_path/'out.vtable'
    write_object(source, ['1\t1\tclass2\tA\tAG', '1\t4\tclass3\tAGC\tAG', '1\t7\tmissing\tA\tAG'], 'vtable')
    monkeypatch.setattr(restriction, 'NORMALIZATION_BATCH_ROWS', 1)
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    args = ['restrict', '--source', str(source), '--output', str(output), '--include-indels', '--norm-indels']
    monkeypatch.setattr(sys, 'argv', args)
    assert restriction.main() == 0
    first_batches = list(tmp_path.glob('out.vtable.bcftools_norm*.input.vcf'))
    assert len(first_batches) > 1
    assert all(f"\trow{i}\t" in "".join(p.read_text() for p in first_batches) for i in range(3))
    assert read_tsv(output.with_name(output.name+'.qc.tsv'))[-1] == ['2', 'missing', 'norm_ref_mismatch']
    monkeypatch.setattr(restriction, 'NORMALIZATION_BATCH_ROWS', 10)
    assert restriction.main() == 0
    assert len(list(tmp_path.glob('out.vtable.bcftools_norm*.input.vcf'))) < len(first_batches)


def test_pre_normalization_ref_mismatch_is_fatal(tmp_path, monkeypatch):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={}, grch37_sequence='AGC')
    candidate = restriction.NormalizationCandidate('row0', 0, restriction.VariantRow('1', '1', 'x', 'A', 'AT'), 'identity')
    with pytest.raises(ValueError, match='unexpected REF mismatch'):
        restriction.normalize_candidates([candidate], fasta_path=tmp_path/'GRCh37.ucsc.fa', contig_naming='ncbi', output_path=tmp_path/'out.vtable')


@pytest.mark.parametrize('a1,a2', [('CC', 'AG'), ('AGC', 'AG'), ('A', 'AG')])
def test_normalization_accepts_both_multibase_output(tmp_path, a1, a2):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={'row0': {'records': [{'pos': 1, 'ref': 'AG', 'alt': 'CC'}]}}, grch37_sequence='AGC')
    source, output = tmp_path/'in.vtable', tmp_path/'out.vtable'
    write_object(source, [f'1\t1\tvariant\t{a1}\t{a2}'], 'vtable')
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', source, '--output', output, '--include-indels', '--norm-indels')
    assert result.returncode == 0, result.stderr
    assert read_tsv(output) == [['1', '1', 'variant', 'CC', 'AG']]


def test_real_normalization_trimming_collisions_and_mnv(tmp_path):
    from real_liftover_helpers import resolve_bcftools_with_liftover
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={}, grch37_sequence='AGCAGC')
    env['MATCH_BCFTOOLS'] = resolve_bcftools_with_liftover()
    source, output = tmp_path/'in.vtable', tmp_path/'out.vtable'
    write_object(source, ['1\t1\tpadded\tATC\tAGC', '1\t2\tsnv\tT\tG', '1\t4\tmnv\tCCT\tAGC'], 'vtable')
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', source, '--output', output, '--include-indels', '--norm-indels')
    assert result.returncode == 0, result.stderr
    assert read_tsv(output) == [['1', '4', 'mnv', 'CCT', 'AGC']]
    assert read_tsv(output.with_name(output.name+'.qc.tsv'))[1:] == [['0', 'padded', 'duplicate_target'], ['1', 'snv', 'duplicate_target']]


def test_bulk_normalization_results_preserve_interleaved_rows_and_long_alleles(tmp_path, monkeypatch):
    """Filtered middle candidates must not shift later normalized values to SNVs."""
    restrict_env_with_fake_norm(tmp_path, behavior={}, grch37_sequence='AGC' * 10)
    long_alt = 'C' * 1000
    source = pd.DataFrame([
        ('1', '1', 'snv', 'C', 'A'),
        ('1', '4', 'normalized', 'A', 'AG'),
        ('1', '7', 'snv2', 'G', 'A'),
        ('1', '10', 'missing', 'A', 'AG'),
        ('1', '13', 'long', 'CC', 'AG'),
    ], columns=['chrom', 'pos', 'id', 'a1', 'a2'], dtype=object)
    def normalize(candidates, **kwargs):
        return {
            candidates[0].candidate_id: restriction.NormalizationOutcome(
                restriction.VariantRow('1', '5', 'normalized', 'T', 'G'), 'identity'),
            candidates[1].candidate_id: restriction.NormalizationOutcome(None, 'identity', 'norm_ref_mismatch'),
            candidates[2].candidate_id: restriction.NormalizationOutcome(
                restriction.VariantRow('1', '14', 'long', long_alt, 'G'), 'identity'),
        }
    monkeypatch.setattr(restriction, 'normalize_candidates', normalize)
    outcomes = restriction.build_restriction_outcomes_table(
        VariantRowsTable.from_frame(source), fasta_path=tmp_path/'GRCh37.ucsc.fa',
        contig_naming='ncbi', allow_strand_flips=False, norm_indels=True,
        genome_build='GRCh37', output_path=tmp_path/'out.vtable', include_indels=True,
    )
    assert outcomes['keep'].tolist() == [True, True, True, False, True]
    assert outcomes['out_pos'].tolist() == ['1', '5', '7', '10', '14']
    assert outcomes['out_a1'].tolist() == ['C', 'T', 'G', 'A', long_alt]
    assert outcomes['status'].tolist() == [None, None, None, 'norm_ref_mismatch', None]
    assert outcomes['out_a1'].dtype == object
