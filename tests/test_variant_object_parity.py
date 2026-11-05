"""Shared target behavior and QC identity across variant object types."""
import json

import pandas as pd
import pytest

from genomatch import vtable_utils
from genomatch.vtable_utils import read_vtable_table, iter_vtable_table_chunks
from test_liftover_build import write_fake_bcftools
from test_restrict_build_compatible import restrict_env_with_fake_norm
from utils import read_tsv, run_py_with_env, write_object_metadata, write_lines, write_match_config, write_primary_ucsc_fasta


TOOLS = [
    'guess_build', 'normalize_contigs', 'restrict_contigs',
    'drop_strand_ambiguous', 'restrict_build_compatible', 'liftover_build', 'sort_variants',
]


def write_input(tmp_path, suffix, rows, *, build='GRCh37', naming='ncbi'):
    source = tmp_path / f'input.{suffix}'
    target = {'genome_build': build}
    if naming is not None:
        target['contig_naming'] = naming
    if suffix == 'vmap':
        # Non-contiguous locators and repeated shard-local indexes expose reindexing.
        rows = [f'{row}\tshard{i % 2}\t{20 + i // 2}\tflip' for i, row in enumerate(rows)]
        metadata = {'object_type': 'variant_map', 'target': target}
    else:
        metadata = {'object_type': 'variant_table', **target}
    write_lines(source, rows)
    write_object_metadata(source.with_name(source.name + '.meta.json'), metadata)
    return source


def qc_locator(suffix, ordinal):
    return [f'shard{ordinal % 2}', str(20 + ordinal // 2)] if suffix == 'vmap' else [str(ordinal)]


@pytest.fixture
def reference_env(tmp_path):
    ref37 = tmp_path / 'hg19.fa'
    ref38 = tmp_path / 'hg38.fa'
    chain = tmp_path / 'lift.chain'
    config = tmp_path / 'config.yaml'
    bcftools = tmp_path / 'bcftools'
    write_primary_ucsc_fasta(ref37, {'chr1': 'AAAAAA'})
    write_primary_ucsc_fasta(ref38, {'chr1': 'CCCCCC'})
    # Two source keys collide, one is unmapped, and the retained rows reorder.
    write_lines(chain, ['chr1\t1\tchr1\t5', 'chr1\t2\tchr1\t5',
                        'chr1\t4\tchr1\t3', 'chr1\t5\tchr1\t1'])
    write_fake_bcftools(bcftools)
    write_match_config(config, grch37_ucsc_fasta=ref37, grch38_ucsc_fasta=ref38, chain37to38=chain)
    return {'MATCH_CONFIG': str(config), 'MATCH_BCFTOOLS': str(bcftools)}


def tool_args(tool, source, output):
    args = ['--source' if tool == 'restrict_build_compatible' else '--input', source]
    if tool == 'guess_build':
        return [*args, '--write']
    args += ['--output', output]
    if tool == 'normalize_contigs':
        args += ['--to', 'ucsc']
    elif tool == 'restrict_contigs':
        args += ['--contigs', '1']
    elif tool == 'liftover_build':
        args += ['--target-build', 'GRCh38']
    return args


@pytest.mark.parametrize('tool', TOOLS)
@pytest.mark.parametrize('empty', [False, True])
def test_target_behavior_parity(tmp_path, reference_env, tool, empty):
    rows = [] if empty else ['1\t5\tshared_id\tc\ta', '1\t4\tshared_id\tt\ta', '1\t3\tother\tg\ta']
    results = []
    for suffix in ['vtable', 'vmap']:
        source = write_input(tmp_path, suffix, rows, build='unknown' if tool == 'guess_build' else 'GRCh37')
        output = tmp_path / f'output.{suffix}'
        qc = output.with_name(output.name + '.qc.tsv')
        write_lines(qc, ['stale'])
        result = run_py_with_env(f'{tool}.py', reference_env, *tool_args(tool, source, output))
        if empty and tool == 'guess_build':
            assert result.returncode != 0
            assert 'empty variant object' in result.stderr
            continue
        assert result.returncode == 0, result.stderr
        artifact = source if tool == 'guess_build' else output
        target_rows = [row[:5] for row in read_tsv(artifact)]
        metadata = json.loads(artifact.with_name(artifact.name + '.meta.json').read_text())
        target_metadata = metadata['target'] if suffix == 'vmap' else {key: metadata[key] for key in ['genome_build', 'contig_naming']}
        if tool != 'guess_build':
            assert metadata['variants_count'] == len(target_rows)
            assert all(len(row) == (8 if suffix == 'vmap' else 5) for row in read_tsv(artifact))
            if qc.exists():
                qc_rows = read_tsv(qc)
                assert qc_rows[0] == (['source_shard', 'source_index', 'source_id', 'status']
                                      if suffix == 'vmap' else ['row_idx', 'source_id', 'status'])
                # Only the third input row is unmapped in this fixture.
                assert qc_rows[1][:-2] == qc_locator(suffix, 2)
                statuses = [row[-2:] for row in qc_rows[1:]]
            else:
                statuses = []
        else:
            statuses = []
        if suffix == 'vmap':
            expected_locators = {row.split('\t')[2]: [] for row in rows}
            for i, row in enumerate(rows):
                expected_locators[row.split('\t')[2]].append(qc_locator(suffix, i))
            for row in read_tsv(artifact):
                assert row[5:7] in expected_locators[row[2]]
        results.append((target_rows, target_metadata, statuses))
    if results:
        assert results[0] == results[1]


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
@pytest.mark.parametrize('tool', TOOLS)
def test_missing_contig_metadata_fails_for_both_types(tmp_path, reference_env, suffix, tool):
    source = write_input(tmp_path, suffix, ['1\t1\tid\tC\tA'], naming=None)
    output = tmp_path / f'out.{suffix}'
    result = run_py_with_env(f'{tool}.py', reference_env, *tool_args(tool, source, output))
    if tool == 'normalize_contigs':
        assert result.returncode == 0, result.stderr  # The explicit repair path.
    else:
        assert result.returncode != 0
        assert 'contig_naming' in result.stderr


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
@pytest.mark.parametrize('tool', [tool for tool in TOOLS if tool not in {'guess_build', 'normalize_contigs'}])
@pytest.mark.parametrize('alt', ['C', 'AC'])
def test_duplicate_input_fails_with_explicit_cleanup_route(tmp_path, reference_env, suffix, tool, alt):
    source = write_input(tmp_path, suffix, [
        f'1\t1\tfirst\t{alt}\tA', '1\t2\tother\tG\tA', f'1\t1\tsecond\t{alt}\tA',
    ])
    output = tmp_path / f'out.{suffix}'
    metadata = source.with_name(source.name + '.meta.json')
    original = metadata.read_bytes()
    result = run_py_with_env(f'{tool}.py', reference_env, *tool_args(tool, source, output))
    assert result.returncode != 0
    assert 'duplicate chrom:pos:a1:a2' in result.stderr
    assert 'sort_variants.py --drop-duplicates' in result.stderr
    assert metadata.read_bytes() == original


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
@pytest.mark.parametrize('mixed', [False, True])
def test_contig_repair_accepts_exact_repeats_and_audits_all_output_duplicates(tmp_path, suffix, mixed):
    source = write_input(tmp_path, suffix, [
        'chr1\t1\tfirst\tC\tA', 'chr1\t1\texact_repeat\tC\tA',
        f'{"1" if mixed else "chr1"}\t1\tcollision\tC\tA',
        'chr1\t2\tretained\tG\tA',
    ], naming=None if mixed else 'ucsc')
    output = tmp_path / f'repaired.{suffix}'
    result = run_py_with_env('normalize_contigs.py', {}, '--input', source, '--output', output,
                             '--to', 'ncbi' if mixed else 'ucsc')
    assert result.returncode == 0, result.stderr
    assert [row[:5] for row in read_tsv(output)] == [[('1' if mixed else 'chr1'), '2', 'retained', 'G', 'A']]
    qc = read_tsv(output.with_name(output.name + '.qc.tsv'))
    assert qc[0] == (['source_shard', 'source_index', 'source_id', 'status'] if suffix == 'vmap'
                     else ['row_idx', 'source_id', 'status'])
    assert qc[1:] == [[*qc_locator(suffix, i), row_id, 'duplicate_target']
                      for i, row_id in enumerate(['first', 'exact_repeat', 'collision'])]
    metadata = json.loads(output.with_name(output.name + '.meta.json').read_text())
    target = metadata['target'] if suffix == 'vmap' else metadata
    assert target['contig_naming'] == ('ncbi' if mixed else 'ucsc')
    assert metadata['variants_count'] == 1
    if suffix == 'vmap':
        assert read_tsv(output)[0][5:] == [*qc_locator(suffix, 3), 'flip']


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
@pytest.mark.parametrize('sample_rows', [0, 3])
def test_build_guess_preserves_repeated_row_weighting_and_only_updates_metadata(tmp_path, suffix, sample_rows):
    source = write_input(tmp_path, suffix, [
        '1\t1\tfirst\tC\tA', '1\t1\trepeat1\tC\tA', '1\t1\trepeat2\tC\tA', '1\t2\tother\tT\tG',
    ], build='unknown')
    fasta37, fasta38 = tmp_path / '37.fa', tmp_path / '38.fa'
    write_primary_ucsc_fasta(fasta37, {'chr1': 'AA'})
    write_primary_ucsc_fasta(fasta38, {'chr1': 'GG'})
    config = tmp_path / 'config.yaml'
    write_match_config(config, grch37_ucsc_fasta=fasta37, grch38_ucsc_fasta=fasta38)
    original = source.read_bytes()
    result = run_py_with_env('guess_build.py', {'MATCH_CONFIG': str(config)}, '--input', source,
                             '--sample-rows', str(sample_rows), '--write')
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary['genome_build'] == 'GRCh37'  # Unique-key weighting would tie the builds.
    assert summary['sample_rows_used'] == (sample_rows or 4)
    evidence = next(item for item in summary['evidence'] if item['genome_build'] == 'GRCh37')
    assert evidence['considered_rows'] == (sample_rows or 4)
    assert evidence['compatible_rows'] >= 2
    metadata = json.loads(source.with_name(source.name + '.meta.json').read_text())
    assert metadata['variants_count'] == 4
    assert (metadata['target'] if suffix == 'vmap' else metadata)['genome_build'] == 'GRCh37'
    assert source.read_bytes() == original


@pytest.mark.parametrize('tool', [
    'restrict_vmap', 'assign_vmap_ids', 'convert_vmap_to_target', 'project_payload',
    'apply_vmap_to_bfile', 'apply_vmap_to_pfile', 'apply_vmap_to_sumstats',
])
def test_source_mapping_readers_reject_duplicate_keys(tmp_path, tool):
    source = write_input(tmp_path, 'vmap', ['1\t1\tfirst\tC\tA', '1\t1\tsecond\tC\tA'])
    membership = tmp_path / 'membership.vtable'
    write_lines(membership, ['1\t1\tlookup\tC\tA'])
    write_object_metadata(membership.with_name(membership.name + '.meta.json'),
               {'object_type': 'variant_table', 'genome_build': 'GRCh37', 'contig_naming': 'ncbi'})
    output = tmp_path / 'out.vmap'
    if tool == 'restrict_vmap':
        args = [source, membership, '--output', output]
    elif tool == 'assign_vmap_ids':
        args = ['--vmap', source, '--id-source', membership, '--output', output]
    elif tool == 'convert_vmap_to_target':
        args = ['--source', source, '--output', tmp_path / 'out.vtable']
    elif tool == 'project_payload':
        args = ['--vmap', source, '--input-format', 'bfile', '--input', tmp_path / 'raw.@.bim', '--output', tmp_path / 'out.@']
    elif tool in {'apply_vmap_to_bfile', 'apply_vmap_to_pfile'}:
        args = ['--vmap', source, '--source-prefix', tmp_path / 'raw.@', '--output-prefix', tmp_path / 'out.@']
    else:
        raw, metadata = tmp_path / 'raw.tsv', tmp_path / 'raw.yaml'
        write_lines(raw, ['CHR\tPOS\tEA\tOA', '1\t1\tC\tA'])
        write_lines(metadata, ['col_CHR: CHR', 'col_POS: POS', 'col_EffectAllele: EA', 'col_OtherAllele: OA'])
        args = ['--vmap', source, '--input', raw, '--sumstats-metadata', metadata, '--output', tmp_path / 'out.tsv']
    result = run_py_with_env(f'{tool}.py', {}, *args)
    assert result.returncode != 0
    assert 'duplicate chrom:pos:a1:a2' in result.stderr
    assert 'sort_variants.py --drop-duplicates' in result.stderr


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
def test_contig_collisions_use_original_qc_locators(tmp_path, suffix):
    source = write_input(tmp_path, suffix, [
        'chrUn\t1\tunresolved\tA\tG', 'chr1\t2\tfirst\tA\tG',
        '1\t2\tsecond\tA\tG', 'chr1\t3\tretained\tC\tA',
    ], naming=None)
    output = tmp_path / f'out.{suffix}'
    result = run_py_with_env('normalize_contigs.py', {}, '--input', source, '--output', output, '--to', 'ncbi')
    assert result.returncode == 0, result.stderr
    assert [row[:5] for row in read_tsv(output)] == [['1', '3', 'retained', 'C', 'A']]
    assert read_tsv(output.with_name(output.name + '.qc.tsv'))[1:] == [
        [*qc_locator(suffix, 1), 'first', 'duplicate_target'],
        [*qc_locator(suffix, 2), 'second', 'duplicate_target'],
    ]


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
def test_reference_collisions_and_exclusions_use_original_locators(tmp_path, reference_env, suffix):
    source = write_input(tmp_path, suffix, [
        '1\t5\tretained\tA\tT', '1\t1\tfirst\tA\tG',
        '1\t3\tmismatch\tC\tG', '1\t1\tsecond\tG\tA',
        '1\t2\talso_retained\tC\tA',
    ])
    output = tmp_path / f'out.{suffix}'
    result = run_py_with_env('restrict_build_compatible.py', reference_env, '--source', source, '--output', output)
    assert result.returncode == 0, result.stderr
    expected = [['1', '5', 'retained', 'T', 'A'], ['1', '2', 'also_retained', 'C', 'A']]
    assert [row[:5] for row in read_tsv(output)] == expected
    assert read_tsv(output.with_name(output.name + '.qc.tsv'))[1:] == [
        [*qc_locator(suffix, 2), 'mismatch', 'no_reference_match'],
        [*qc_locator(suffix, 1), 'first', 'duplicate_target'],
        [*qc_locator(suffix, 3), 'second', 'duplicate_target'],
    ]
    summary = json.loads(result.stdout)
    assert (summary['output_rows'], summary['dropped_rows'], summary['duplicate_target_rows']) == (2, 3, 2)


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
def test_liftover_collisions_and_reordered_results_preserve_qc_locators(tmp_path, reference_env, suffix):
    source = write_input(tmp_path, suffix, [f'1\t{i}\tid{i}\tC\tA' for i in range(1, 6)])
    output = tmp_path / f'out.{suffix}'
    result = run_py_with_env('liftover_build.py', reference_env, '--input', source, '--output', output, '--target-build', 'GRCh38')
    assert result.returncode == 0, result.stderr
    assert [row[:5] for row in read_tsv(output)] == [['1', '1', 'id5', 'C', 'A'], ['1', '3', 'id4', 'C', 'A']]
    assert read_tsv(output.with_name(output.name + '.qc.tsv'))[1:] == [
        [*qc_locator(suffix, 0), 'id1', 'duplicate_target'],
        [*qc_locator(suffix, 1), 'id2', 'duplicate_target'],
        [*qc_locator(suffix, 2), 'id3', 'unmapped'],
    ]


def test_vtable_loader_ordinals_continue_across_chunks_and_survive_reindexing(tmp_path, monkeypatch):
    source = write_input(tmp_path, 'vtable', [f'1\t{i + 1}\tid{i}\tC\tA' for i in range(5)])
    monkeypatch.setattr(vtable_utils, 'WRITE_CHUNK_ROWS', 2)
    chunks = list(iter_vtable_table_chunks(source))
    frame = pd.concat([chunk.frame for chunk in chunks], ignore_index=True)
    assert frame['row_idx'].tolist() == [0, 1, 2, 3, 4]
    assert read_vtable_table(source).frame['row_idx'].tolist() == frame['row_idx'].tolist()
    transformed = vtable_utils.sort_target_table_by_declared_coordinate(
        frame.iloc[[4, 1, 3]].reset_index(drop=True), 'ncbi', label='test',
    ).reset_index(drop=True)
    assert transformed['row_idx'].tolist() == [1, 3, 4]


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
def test_normalization_round_trip_preserves_qc_locators(tmp_path, suffix):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={
        'row1': {'records': [{'pos': 1, 'ref': 'AA', 'alt': 'A'}]},
        'row2': {'records': [{'pos': 1, 'ref': 'AA', 'alt': 'A'}]},
    })
    source = write_input(tmp_path, suffix, [
        '1\t6\tretained\tC\tA', '1\t2\tfirst\tA\tAT', '1\t4\tsecond\tA\tAC',
    ])
    output = tmp_path / f'out.{suffix}'
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', source, '--output', output, '--norm-indels', '--include-indels', '--allow-indel-swaps')
    assert result.returncode == 0, result.stderr
    assert [row[:5] for row in read_tsv(output)] == [['1', '6', 'retained', 'C', 'A']]
    assert read_tsv(output.with_name(output.name + '.qc.tsv'))[1:] == [
        [*qc_locator(suffix, 1), 'first', 'duplicate_target'],
        [*qc_locator(suffix, 2), 'second', 'duplicate_target'],
    ]


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
def test_sort_cleanup_keeps_first_exact_key_without_qc(tmp_path, suffix):
    source = write_input(tmp_path, suffix, [
        '1\t2\tlater\tC\tA', '1\t1\tfirst\tG\tA', '1\t1\tsecond\tG\tA',
        '1\t1\treversed\tA\tG',
    ])
    output = tmp_path / f'out.{suffix}'
    qc = output.with_name(output.name + '.qc.tsv')
    write_lines(qc, ['stale'])
    result = run_py_with_env('sort_variants.py', {}, '--input', source, '--output', output, '--drop-duplicates')
    assert result.returncode == 0, result.stderr
    assert [row[:5] for row in read_tsv(output)] == [
        ['1', '1', 'first', 'G', 'A'], ['1', '1', 'reversed', 'A', 'G'], ['1', '2', 'later', 'C', 'A'],
    ]
    assert not qc.exists()
    if suffix == 'vmap':
        assert read_tsv(output)[0][5:] == [*qc_locator(suffix, 1), 'flip']


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
def test_normalized_row_colliding_with_unchanged_row_drops_both(tmp_path, suffix):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={
        'row1': {'records': [{'pos': 1, 'ref': 'A', 'alt': 'C'}]},
    })
    source = write_input(tmp_path, suffix, ['1\t1\tunchanged\tC\tA', '1\t4\tnormalized\tA\tAC'])
    output = tmp_path / f'out.{suffix}'
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', source, '--output', output, '--norm-indels', '--include-indels', '--allow-indel-swaps')
    assert result.returncode == 0, result.stderr
    assert read_tsv(output) == []
    assert read_tsv(output.with_name(output.name + '.qc.tsv'))[1:] == [
        [*qc_locator(suffix, 0), 'unchanged', 'duplicate_target'],
        [*qc_locator(suffix, 1), 'normalized', 'duplicate_target'],
    ]


@pytest.mark.parametrize('field', ['source_shard', 'source_index'])
def test_vtable_qc_rejects_provenance_fields(tmp_path, field):
    frame = pd.DataFrame({'row_idx': [7], 'id': ['id'], 'status': ['duplicate_target'], field: [0]})
    with pytest.raises(ValueError, match='must not contain .vmap provenance fields'):
        vtable_utils.write_variant_status_qc(tmp_path / 'out.qc.tsv', frame, is_vmap=False)
    assert not (tmp_path / 'out.qc.tsv').exists()


def test_liftover_vtable_keeps_row_idx_through_all_stages(tmp_path, reference_env, monkeypatch):
    import sys
    from genomatch import vmap_liftover_build as liftover

    source = write_input(tmp_path, 'vtable', ['1\t3\tunmapped\tC\tA', '1\t1\tretained\tG\tA'])
    output = tmp_path / 'out.vtable'
    for key, value in reference_env.items():
        monkeypatch.setenv(key, value)
    original_prepare = liftover.prepare_liftover_table
    original_write = liftover.write_temp_vcf
    original_parse = liftover.parse_lifted_vcf

    def check_identifiers(frame):
        assert 'row_idx' in frame.columns
        assert 'source_index' not in frame.columns
        assert 'source_shard' not in frame.columns

    def prepare(frame, *, is_vmap):
        assert not is_vmap
        check_identifiers(frame)
        result = original_prepare(frame, is_vmap=is_vmap)
        check_identifiers(result)
        return result

    def write(path, prepared, *, is_vmap):
        result = original_write(path, prepared, is_vmap=is_vmap)
        check_identifiers(result)
        return result

    def parse(path, lookup, *args, is_vmap):
        check_identifiers(lookup)
        result, statuses = original_parse(path, lookup, *args, is_vmap=is_vmap)
        check_identifiers(result.frame)
        assert result.frame['row_idx'].tolist() == [1]
        return result, statuses

    monkeypatch.setattr(liftover, 'prepare_liftover_table', prepare)
    monkeypatch.setattr(liftover, 'write_temp_vcf', write)
    monkeypatch.setattr(liftover, 'parse_lifted_vcf', parse)
    monkeypatch.setattr(sys, 'argv', ['liftover_build', '--input', str(source), '--output', str(output), '--target-build', 'GRCh38'])
    assert liftover.main() == 0
    assert read_tsv(output.with_name(output.name + '.qc.tsv')) == [
        ['row_idx', 'source_id', 'status'], ['0', 'unmapped', 'unmapped'],
    ]


@pytest.mark.parametrize('suffix', ['vtable', 'vmap'])
def test_set_membership_reads_exclude_qc_identifiers(tmp_path, suffix):
    from genomatch.exact_set_utils import read_target_table

    source = write_input(tmp_path, suffix, ['1\t1\tid\tC\tA'])
    assert list(read_target_table(source).frame.columns) == ['chrom', 'pos', 'id', 'a1', 'a2']
