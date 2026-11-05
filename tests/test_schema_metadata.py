"""Schema boundaries, historical bytes, writer selection, and resume validation."""
import importlib
import json
from pathlib import Path
import shutil
import sys

import pytest

from genomatch import vtable_utils as objects
from genomatch import workflow_wrapper_utils as workflow
from utils import REPO_ROOT, read_tsv, run_py, write_json, write_lines, write_object_metadata

FIXTURES = REPO_ROOT / 'tests/fixtures/schema'
MISSING = object()
BAD_VERSIONS = [
    pytest.param(MISSING, id='missing'), pytest.param(None, id='null'),
    pytest.param(1.0, id='number'), pytest.param(True, id='boolean'),
    pytest.param('', id='empty'), pytest.param(' 1.0', id='leading-space'),
    pytest.param('1.0 ', id='trailing-space'), pytest.param('2.0', id='unsupported'),
    pytest.param(['1.0'], id='list'), pytest.param({'version': '1.0'}, id='object'),
]
EXPECTED_ROWS = [
    ['2', '20', 'rs_snv', 'T', 'C'],
    ['1', '100', '1:100:AC:A', 'A', 'AC'],
    ['1', '100', 'rs_multibase', 'GT', 'ACG'],
]


def copy_object(tmp_path, suffix, name='input', version='1.0', malformed=False):
    source = FIXTURES / f'{suffix}_1_0.{suffix}'
    output = tmp_path / f'{name}.{suffix}'
    shutil.copyfile(source, output)
    metadata = json.loads(objects.metadata_path_for(source).read_text())
    if version is MISSING:
        metadata.pop('schema_version')
    else:
        metadata['schema_version'] = version
    write_json(objects.metadata_path_for(output), metadata)
    if malformed:
        write_lines(output, ['malformed row'])
    return output


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
def test_checked_in_schema_1_0_fixtures_load_with_original_meaning(suffix):
    path = FIXTURES / f'{suffix}_1_0.{suffix}'
    assert (FIXTURES / 'generate_1_0.py').is_file()
    assert (FIXTURES / 'README.md').is_file()
    before = (path.read_bytes(), objects.metadata_path_for(path).read_bytes())
    for loader in (objects.load_variant_object, objects.load_variant_object_tables):
        loaded = loader(path)
        rows = loaded.target_rows if hasattr(loaded, 'target_rows') else loaded.target_rows_table.to_rows()
        assert [[row.chrom, row.pos, row.id, row.a1, row.a2] for row in rows] == EXPECTED_ROWS
        assert loaded.target_metadata == {'genome_build': 'GRCh37', 'contig_naming': 'ncbi'}
        expected_meta = {
            'object_type': 'variant_map' if suffix == 'vmap' else 'variant_table',
            'schema_version': '1.0', 'variants_count': 3,
            'created_by': 'generate_1_0.py', 'derived_from': 'schema-fixture',
        }
        if suffix == 'vmap':
            expected_meta['target'] = loaded.target_metadata
            mapped = loaded.base_vmap_rows if hasattr(loaded, 'base_vmap_rows') else loaded.base_vmap_table.to_rows()
            assert [(row.source_shard, row.source_index, row.allele_op) for row in mapped] == [
                ('chr2', 7, 'identity'), ('chr1', 2, 'swap'), ('chr1', 8, 'flip_swap'),
            ]
        else:
            expected_meta.update(loaded.target_metadata)
        assert loaded.raw_metadata == expected_meta
    assert (path.read_bytes(), objects.metadata_path_for(path).read_bytes()) == before


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
@pytest.mark.parametrize('version', BAD_VERSIONS)
def test_reads_and_metadata_updates_reject_invalid_versions(tmp_path, suffix, version):
    path = copy_object(tmp_path, suffix, version=version, malformed=True)
    original = objects.metadata_path_for(path).read_bytes()
    metadata = json.loads(original)
    for operation in (lambda: objects.load_variant_object_tables(path),
                      lambda: objects.update_metadata_in_place(path, metadata),
                      lambda: workflow.read_target_metadata(path)):
        with pytest.raises(ValueError, match='schema_version') as error:
            operation()
        assert str(path) in str(error.value)
        assert 'regenerate from the original raw inputs' in str(error.value)
    assert objects.metadata_path_for(path).read_bytes() == original


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
def test_current_writer_roundtrip_and_optional_row_count(tmp_path, suffix):
    source = objects.load_variant_object(FIXTURES / f'{suffix}_1_0.{suffix}')
    output = tmp_path / f'written.{suffix}'
    if suffix == 'vmap':
        objects.write_vmap(output, source.base_vmap_rows)
        metadata = objects.make_vmap_metadata(source.target_metadata)
    else:
        objects.write_vtable(output, source.target_rows)
        metadata = objects.make_vtable_metadata('GRCh37', 'ncbi')
    assert metadata['schema_version'] == '1.0'
    assert 'variants_count' not in metadata
    objects.write_metadata(output, metadata)
    assert objects.load_variant_object(output).target_rows == source.target_rows


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
@pytest.mark.parametrize('version', BAD_VERSIONS)
def test_metadata_writer_sets_current_version_without_mutating_input(tmp_path, suffix, version):
    source = copy_object(tmp_path, suffix)
    metadata = objects.load_metadata(source)
    if version is MISSING:
        metadata.pop('schema_version')
    else:
        metadata['schema_version'] = version
    original = json.loads(json.dumps(metadata))
    output = tmp_path / f'output.{suffix}'
    objects.write_metadata(output, metadata)
    assert objects.load_metadata(output) == {**original, 'schema_version': '1.0'}
    assert metadata == original


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
def test_metadata_writer_validates_other_fields_before_overwriting(tmp_path, suffix):
    path = copy_object(tmp_path, suffix)
    metadata = objects.load_metadata(path)
    before = objects.metadata_path_for(path).read_bytes()
    metadata['variants_count'] = -1
    with pytest.raises(ValueError, match='variants_count'):
        objects.write_metadata(path, metadata)
    assert objects.metadata_path_for(path).read_bytes() == before


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
def test_metadata_only_update_rejects_a_different_supported_version(tmp_path, monkeypatch, suffix):
    monkeypatch.setattr(objects, 'SUPPORTED_OBJECT_SCHEMAS', {
        kind: frozenset({'1.0', 'test-old'}) for kind in ('variant_map', 'variant_table')
    })
    path = copy_object(tmp_path, suffix, version='test-old')
    metadata = objects.load_metadata(path)
    metadata['schema_version'] = '1.0'
    before = (path.read_bytes(), objects.metadata_path_for(path).read_bytes())
    with pytest.raises(ValueError, match='must preserve schema_version'):
        objects.update_metadata_in_place(path, metadata)
    assert (path.read_bytes(), objects.metadata_path_for(path).read_bytes()) == before


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
def test_metadata_only_update_validates_replacement_before_overwriting(tmp_path, suffix):
    path = copy_object(tmp_path, suffix)
    metadata = objects.load_metadata(path)
    metadata.pop('schema_version')
    before = objects.metadata_path_for(path).read_bytes()
    with pytest.raises(ValueError, match='schema_version'):
        objects.update_metadata_in_place(path, metadata)
    assert objects.metadata_path_for(path).read_bytes() == before


def test_sidecar_helper_defaults_version_and_preserves_explicit_values(tmp_path):
    path = tmp_path / 'sidecar.json'
    payload = {'object_type': 'variant_table'}
    write_object_metadata(path, payload)
    assert json.loads(path.read_text()) == {**payload, 'schema_version': objects.CURRENT_SCHEMA_VERSION}
    assert payload == {'object_type': 'variant_table'}
    write_object_metadata(path, {**payload, 'schema_version': None})
    assert json.loads(path.read_text())['schema_version'] is None
    write_json(path, payload)
    assert 'schema_version' not in json.loads(path.read_text())


@pytest.mark.parametrize('old_index', [0, 1])
@pytest.mark.parametrize('include_count', [False, True])
@pytest.mark.parametrize('empty', [False, True])
def test_raw_concatenation_rejects_any_supported_noncurrent_input(
        tmp_path, monkeypatch, old_index, include_count, empty):
    from genomatch.vmap_prepare_variants_sharded import write_concatenated_vmap

    monkeypatch.setattr(objects, 'SUPPORTED_OBJECT_SCHEMAS', {
        kind: frozenset({'1.0', 'test-old'}) for kind in ('variant_map', 'variant_table')
    })
    inputs = [copy_object(tmp_path, 'vmap', name=f'group{i}') for i in range(2)]
    old = inputs[old_index]
    metadata = objects.load_metadata(old)
    metadata['schema_version'] = 'test-old'
    if empty:
        old.write_bytes(b'')
        metadata['variants_count'] = 0
    if not include_count:
        metadata.pop('variants_count')
    write_json(objects.metadata_path_for(old), metadata)
    output = copy_object(tmp_path, 'vmap', name='output')
    original = (output.read_bytes(), objects.metadata_path_for(output).read_bytes())
    with pytest.raises(ValueError, match='raw concatenation requires schema_version') as error:
        write_concatenated_vmap(output, inputs, objects.load_metadata(inputs[0]))
    assert str(old) in str(error.value)
    assert 'test-old' in str(error.value)
    assert 'regenerate from the original raw inputs' in str(error.value)
    assert (output.read_bytes(), objects.metadata_path_for(output).read_bytes()) == original
    assert not (tmp_path / 'output.tmp.vmap').exists()


@pytest.mark.parametrize('include_count', [False, True])
def test_raw_concatenation_preserves_current_schema_row_bytes(tmp_path, include_count):
    from genomatch.vmap_prepare_variants_sharded import write_concatenated_vmap

    inputs = [copy_object(tmp_path, 'vmap', name=f'group{i}') for i in range(2)]
    if not include_count:
        for path in inputs:
            metadata = objects.load_metadata(path)
            metadata.pop('variants_count')
            write_json(objects.metadata_path_for(path), metadata)
    output = tmp_path / 'output.vmap'
    write_concatenated_vmap(output, inputs, objects.load_metadata(inputs[0]))
    assert output.read_bytes() == b''.join(path.read_bytes() for path in inputs)
    metadata = objects.load_metadata(output)
    assert metadata['schema_version'] == objects.CURRENT_SCHEMA_VERSION
    assert metadata['variants_count'] == 6


# Every command-line object-reading path exercises the same schema boundary.
# Calling main with argv keeps the malformed-input matrix independent of shell
# startup and external FASTA / PLINK assets.
READ_PATHS = [
    (tool, suffix) for tool in ('normalize_contigs', 'restrict_contigs',
        'drop_strand_ambiguous', 'restrict_build_compatible', 'liftover_build',
        'sort_variants', 'guess_build', 'intersect_variants', 'union_variants')
    for suffix in ('vmap', 'vtable')
] + [
    ('convert_vmap_to_target', 'vmap'),
    ('restrict_vmap-source', 'vmap'), ('restrict_vmap-membership', 'vmap'),
    ('restrict_vmap-membership', 'vtable'), ('assign_vmap_ids-source', 'vmap'),
    ('assign_vmap_ids-membership', 'vmap'), ('assign_vmap_ids-membership', 'vtable'),
    ('import_sumstats', 'vmap'), ('import_sumstats', 'vtable'),
    ('apply_vmap_to_bfile', 'vmap'), ('apply_vmap_to_pfile', 'vmap'),
    ('apply_vmap_to_sumstats', 'vmap'), ('project_payload', 'vmap'),
]


def command_args(tool, source, partner, output, tmp_path):
    raw = tmp_path / 'raw.tsv'
    yaml = tmp_path / 'raw.yaml'
    write_lines(raw, ['SNP\tEA\tOA\tBETA', 'rs_snv\tT\tC\t0.5'])
    write_lines(yaml, ['col_SNP: SNP', 'col_EffectAllele: EA',
                       'col_OtherAllele: OA', 'col_BETA: BETA'])
    if tool == 'import_sumstats':
        return ['--input', raw, '--sumstats-metadata', yaml, '--id-lookup', source, '--output', output]
    if tool in ('apply_vmap_to_sumstats', 'project_payload'):
        return ['--input', raw, '--sumstats-metadata', yaml, '--vmap', source,
                '--output', output, *(['--input-format', 'sumstats'] if tool == 'project_payload' else [])]
    if tool in ('apply_vmap_to_bfile', 'apply_vmap_to_pfile'):
        prefix = tmp_path / 'raw'
        for ext in ('.bed', '.bim', '.fam', '.pgen', '.pvar', '.psam'):
            prefix.with_suffix(ext).touch()
        return ['--source-prefix', prefix, '--vmap', source, '--output-prefix', output]
    if tool.startswith(('restrict_vmap', 'assign_vmap_ids')):
        first, second = (source, partner) if tool.endswith('source') else (partner, source)
        if tool.startswith('assign'):
            return ['--vmap', first, '--id-source', second, '--output', output]
        return [first, second, '--output', output]
    if tool in ('intersect_variants', 'union_variants'):
        return [partner, source, '--output', output]
    flag = '--source' if tool in ('restrict_build_compatible', 'convert_vmap_to_target') else '--input'
    args = [flag, source]
    if tool != 'guess_build':
        args += ['--output', output]
    extras = {'normalize_contigs': ['--to', 'ncbi'], 'restrict_contigs': ['--contigs', '1'],
              'liftover_build': ['--target-build', 'GRCh38'], 'guess_build': ['--write']}
    return args + extras.get(tool, [])


@pytest.mark.parametrize('tool,suffix', READ_PATHS)
@pytest.mark.parametrize('version', BAD_VERSIONS)
def test_cli_read_paths_check_schema_before_malformed_rows(tmp_path, monkeypatch, tool, suffix, version):
    source = copy_object(tmp_path, suffix, version=version, malformed=True)
    partner = copy_object(tmp_path, 'vmap', name='partner')
    output_suffix = 'vtable' if tool.split('-')[0] in ('intersect_variants', 'union_variants', 'convert_vmap_to_target') else suffix
    if tool.split('-')[0] in ('restrict_vmap', 'assign_vmap_ids', 'import_sumstats'):
        output_suffix = 'vmap'
    output = tmp_path / f'output.{output_suffix}'
    args = command_args(tool, source, partner, output, tmp_path)
    module = importlib.import_module(f"genomatch.{tool.split('-')[0]}")
    monkeypatch.setattr(sys, 'argv', [tool, *map(str, args)])
    with pytest.raises(ValueError, match='schema_version') as error:
        module.main()
    assert str(source) in str(error.value)
    assert 'regenerate from the original raw inputs' in str(error.value)
    assert not output.exists()


@pytest.mark.parametrize('wrapper,stage', [
    ('prepare_variants', 'final'), ('prepare_variants', 'imported'),
    ('prepare_variants_sharded', 'final'), ('prepare_variants_sharded', 'group'),
])
@pytest.mark.parametrize('version', BAD_VERSIONS)
def test_resume_validates_final_and_skipped_objects(tmp_path, monkeypatch, wrapper, stage, version):
    sharded = wrapper.endswith('_sharded')
    raw = tmp_path / ('raw.1.bim' if sharded else 'raw.bim')
    write_lines(raw, ['1\trs1\t0\t1\tG\tA'])
    name = {'final': 'output', 'imported': 'output.imported', 'group': 'group.1'}[stage]
    existing = copy_object(tmp_path, 'vmap', name=name, version=version, malformed=True)
    original = (existing.read_bytes(), objects.metadata_path_for(existing).read_bytes())
    args = ['--input', str(tmp_path / 'raw.@.bim') if sharded else str(raw),
            '--input-format', 'bim', '--output', str(tmp_path / 'output'), '--resume']
    if sharded:
        args += ['--prefix', str(tmp_path / 'group.@')]
    monkeypatch.setattr(sys, 'argv', [wrapper, *args])
    module = importlib.import_module(f'genomatch.{wrapper}')
    with pytest.raises(ValueError, match='schema_version') as error:
        module.main()
    assert str(existing) in str(error.value)
    assert 'regenerate from the original raw inputs' in str(error.value)
    assert (existing.read_bytes(), objects.metadata_path_for(existing).read_bytes()) == original
    if stage != 'final':
        assert not (tmp_path / 'output.vmap').exists()


@pytest.mark.parametrize('wrapper', ['prepare_variants', 'prepare_variants_sharded'])
def test_resume_supported_final_succeeds_without_rewriting(tmp_path, monkeypatch, wrapper):
    sharded = wrapper.endswith('_sharded')
    raw = tmp_path / ('raw.1.bim' if sharded else 'raw.bim')
    write_lines(raw, ['1\trs1\t0\t1\tG\tA'])
    existing = copy_object(tmp_path, 'vmap', name='output')
    original = (existing.read_bytes(), objects.metadata_path_for(existing).read_bytes())
    args = ['--input', str(tmp_path / 'raw.@.bim') if sharded else str(raw),
            '--input-format', 'bim', '--output', str(tmp_path / 'output'), '--resume']
    if sharded:
        args += ['--prefix', str(tmp_path / 'group.@')]
    monkeypatch.setattr(sys, 'argv', [wrapper, *args])
    assert importlib.import_module(f'genomatch.{wrapper}').main() == 0
    assert (existing.read_bytes(), objects.metadata_path_for(existing).read_bytes()) == original


def test_cli_reports_schema_error_with_remedy(tmp_path):
    source = copy_object(tmp_path, 'vmap', version=MISSING, malformed=True)
    result = run_py('guess_build.py', '--input', source, '--write')
    assert result.returncode != 0
    assert str(source) in result.stderr
    assert 'schema_version' in result.stderr
    assert 'regenerate from the original raw inputs' in result.stderr


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
@pytest.mark.parametrize('version', [pytest.param(MISSING, id='missing'), 1.0, '2.0'])
def test_preparation_id_lookup_rejects_schema_before_rows(tmp_path, suffix, version):
    source = copy_object(tmp_path, suffix, version=version, malformed=True)
    raw, yaml = tmp_path / 'raw.tsv', tmp_path / 'raw.yaml'
    write_lines(raw, ['SNP\tEA\tOA\tBETA', 'rs_snv\tT\tC\t0.5'])
    write_lines(yaml, ['col_SNP: SNP', 'col_EffectAllele: EA', 'col_OtherAllele: OA', 'col_BETA: BETA'])
    result = run_py('prepare_variants.py', '--input', raw, '--input-format', 'sumstats',
                    '--sumstats-metadata', yaml, '--id-lookup', source,
                    '--output', tmp_path / 'prepared', '--src-build', 'GRCh37')
    assert result.returncode != 0
    assert str(source) in result.stderr
    assert 'schema_version' in result.stderr
    assert 'regenerate from the original raw inputs' in result.stderr
    assert not (tmp_path / 'prepared.vmap').exists()


PRODUCERS = [
    (tool, suffix) for tool in ('sort_variants', 'normalize_contigs', 'restrict_contigs',
        'drop_strand_ambiguous', 'restrict_build_compatible', 'liftover_build',
        'intersect_variants', 'union_variants') for suffix in ('vmap', 'vtable')
] + [('convert_vmap_to_target', 'vmap'), ('restrict_vmap-source', 'vmap'),
     ('assign_vmap_ids-source', 'vmap')]


@pytest.mark.parametrize('tool,suffix', PRODUCERS)
def test_serializing_tools_choose_writer_schema_instead_of_inheriting(tmp_path, monkeypatch, tool, suffix):
    from test_prepare_variants import base_env
    # Simulate a reader supporting another version with the same row layout.
    # This is test-only support, not a second public schema or a relabelled
    # historical compatibility fixture. It exposes blind metadata copying.
    monkeypatch.setattr(objects, 'SUPPORTED_OBJECT_SCHEMAS', {
        kind: frozenset({'1.0', 'test-old'}) for kind in ('variant_map', 'variant_table')
    })
    source = tmp_path / f'source.{suffix}'
    target = {'genome_build': 'GRCh37', 'contig_naming': 'ncbi'}
    meta = {'object_type': 'variant_map', 'schema_version': 'test-old', 'target': target} if suffix == 'vmap' else {
        'object_type': 'variant_table', 'schema_version': 'test-old', **target,
    }
    rows = ['1\t1\trs1\tG\tA', '1\t2\trs2\tC\tA']
    if suffix == 'vmap':
        rows = [f'{row}\t.\t{idx}\tidentity' for idx, row in enumerate(rows)]
    write_lines(source, rows)
    write_json(objects.metadata_path_for(source), meta)
    partner = copy_object(tmp_path, 'vmap', name='partner')
    # Use matching partner keys so set tools and ID assignment emit rows.
    write_lines(partner, ['1\t1\tassigned1\tG\tA\t.\t0\tidentity',
                          '1\t2\tassigned2\tC\tA\t.\t1\tidentity'])
    partner_meta = objects.load_metadata(partner)
    partner_meta['variants_count'] = 2
    write_json(objects.metadata_path_for(partner), partner_meta)
    output_suffix = 'vtable' if tool in ('convert_vmap_to_target', 'intersect_variants', 'union_variants') else suffix
    output = tmp_path / f'output.{output_suffix}'
    env = base_env(tmp_path, grch37_sequences={'chr1': 'AAAAAA'}, grch38_sequences={'chr1': 'AAAAAA'})
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    args = command_args(tool, source, partner, output, tmp_path)
    # Sets must choose a possibly old-schema driver too.
    if tool in ('intersect_variants', 'union_variants'):
        args = [source, partner, '--output', output]
    monkeypatch.setattr(sys, 'argv', [tool, *map(str, args)])
    assert importlib.import_module(f"genomatch.{tool.split('-')[0]}").main() == 0
    assert read_tsv(output)
    assert objects.load_metadata(output)['schema_version'] == '1.0'


@pytest.mark.parametrize('tool', ['import_bim', 'import_pvar', 'import_vcf', 'import_sumstats',
                                  'prepare_variants', 'prepare_variants_sharded'])
def test_importers_and_wrappers_emit_current_schema(tmp_path, tool):
    from test_prepare_variants import base_env
    from utils import run_py_with_env
    source = tmp_path / 'raw.bim'
    write_lines(source, ['1\trs1\t0\t1\tG\tA'])
    output = tmp_path / 'output.vmap'
    args = ['--input', source, '--output', output]
    if tool in ('import_pvar', 'import_vcf'):
        source = tmp_path / ('raw.pvar' if tool == 'import_pvar' else 'raw.vcf')
        write_lines(source, ['##fileformat=VCFv4.2', '#CHROM\tPOS\tID\tREF\tALT', '1\t1\trs1\tA\tG'])
        args = ['--input', source, '--output', output]
    elif tool == 'import_sumstats':
        source, yaml = tmp_path / 'raw.tsv', tmp_path / 'raw.yaml'
        write_lines(source, ['CHR\tPOS\tEA\tOA', '1\t1\tG\tA'])
        write_lines(yaml, ['col_CHR: CHR', 'col_POS: POS', 'col_EffectAllele: EA', 'col_OtherAllele: OA'])
        args = ['--input', source, '--sumstats-metadata', yaml, '--output', output]
    elif tool.startswith('prepare_variants'):
        args += ['--input-format', 'bim', '--src-build', 'GRCh37', '--dst-build', 'GRCh37']
        if tool.endswith('_sharded'):
            source.rename(tmp_path / 'raw.1.bim')
            args[1] = tmp_path / 'raw.@.bim'
            args += ['--prefix', tmp_path / 'group.@']
    env = base_env(tmp_path, grch37_sequences={'chr1': 'AAAAAA'}, grch38_sequences={'chr1': 'CCCCCC'})
    result = run_py_with_env(f'{tool}.py', env, *args)
    assert result.returncode == 0, result.stderr
    assert read_tsv(output)
    for sidecar in tmp_path.glob('*.vmap.meta.json'):
        assert json.loads(sidecar.read_text())['schema_version'] == '1.0', sidecar
    assert objects.load_metadata(output)['schema_version'] == '1.0'


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
def test_guess_build_metadata_edit_preserves_supported_input_schema(tmp_path, monkeypatch, suffix):
    module = importlib.import_module('genomatch.guess_build')
    monkeypatch.setattr(objects, 'SUPPORTED_OBJECT_SCHEMAS', {
        kind: frozenset({'1.0', 'test-old'}) for kind in ('variant_map', 'variant_table')
    })
    source = copy_object(tmp_path, suffix, version='test-old')
    original = source.read_bytes()
    monkeypatch.setattr(module, 'guess_build', lambda rows, naming: {
        'genome_build': 'GRCh38', 'confidence': 'high', 'evidence': [],
    })
    monkeypatch.setattr(sys, 'argv', ['guess_build', '--input', str(source), '--write', '--force'])
    assert module.main() == 0
    assert objects.load_metadata(source)['schema_version'] == 'test-old'
    assert workflow.read_target_build(source) == 'GRCh38'
    assert source.read_bytes() == original
