"""Unit tests for `.vmap` and `.vtable` schema versions."""
import importlib
import json
import shutil
import sys

import pytest

from genomatch import vtable_utils as objects
from utils import REPO_ROOT, write_json, write_lines

FIXTURES = REPO_ROOT / 'tests/fixtures/schema'
MISSING = object()
EXPECTED_ROWS = [
    ['2', '20', 'rs_snv', 'T', 'C'],
    ['1', '100', '1:100:AC:A', 'A', 'AC'],
    ['1', '100', 'rs_multibase', 'GT', 'ACG'],
]


def fixture_path(suffix):
    return FIXTURES / f'{suffix}_1_0.{suffix}'


def copy_object(tmp_path, suffix, name='input', version='1.0', malformed=False):
    source = fixture_path(suffix)
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


def object_bytes(path):
    return path.read_bytes(), objects.metadata_path_for(path).read_bytes()


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
def test_checked_in_fixture_loads_with_original_meaning(suffix):
    loaded = objects.load_variant_object(fixture_path(suffix))
    assert [[row.chrom, row.pos, row.id, row.a1, row.a2] for row in loaded.target_rows] == EXPECTED_ROWS
    assert loaded.target_metadata == {'genome_build': 'GRCh37', 'contig_naming': 'ncbi'}
    assert loaded.raw_metadata['schema_version'] == '1.0'
    assert loaded.raw_metadata['variants_count'] == 3
    if suffix == 'vmap':
        assert [(row.source_shard, row.source_index, row.allele_op) for row in loaded.base_vmap_rows] == [
            ('chr2', 7, 'identity'), ('chr1', 2, 'swap'), ('chr1', 8, 'flip_swap'),
        ]


@pytest.mark.parametrize('suffix', ['vmap', 'vtable'])
def test_write_then_read_sets_current_version(tmp_path, suffix):
    source = objects.load_variant_object(fixture_path(suffix))
    output = tmp_path / f'written.{suffix}'
    if suffix == 'vmap':
        objects.write_vmap(output, source.base_vmap_rows)
    else:
        objects.write_vtable(output, source.target_rows)
    # An input's version must not carry over to newly written rows.
    metadata = {**source.raw_metadata, 'schema_version': 'older'}
    objects.write_metadata(output, metadata)
    assert metadata['schema_version'] == 'older'
    written = objects.load_variant_object(output)
    assert written.raw_metadata['schema_version'] == objects.CURRENT_SCHEMA_VERSION
    assert written.target_rows == source.target_rows


@pytest.mark.parametrize('suffix,version', [
    pytest.param('vmap', MISSING, id='vmap-missing'),
    pytest.param('vtable', '2.0', id='vtable-unsupported'),
])
def test_read_refuses_invalid_version_before_parsing_rows(tmp_path, suffix, version):
    path = copy_object(tmp_path, suffix, version=version, malformed=True)
    with pytest.raises(ValueError, match='schema_version') as error:
        objects.load_variant_object_tables(path)
    assert str(path) in str(error.value)
    assert 'regenerate from the original raw inputs' in str(error.value)


def test_resume_refuses_legacy_final_output(tmp_path, monkeypatch):
    raw = tmp_path / 'raw.bim'
    write_lines(raw, ['1\trs1\t0\t1\tG\tA'])
    final = copy_object(tmp_path, 'vmap', name='output', version=MISSING)
    before = object_bytes(final)
    monkeypatch.setattr(sys, 'argv', [
        'prepare_variants', '--input', str(raw), '--input-format', 'bim',
        '--output', str(tmp_path / 'output'), '--resume',
    ])
    module = importlib.import_module('genomatch.prepare_variants')
    with pytest.raises(ValueError, match='schema_version') as error:
        module.main()
    assert str(final) in str(error.value)
    assert object_bytes(final) == before


def test_metadata_only_update_refuses_version_change(tmp_path, monkeypatch):
    monkeypatch.setattr(objects, 'SUPPORTED_OBJECT_SCHEMAS', {
        kind: frozenset({'1.0', 'older'}) for kind in ('variant_map', 'variant_table')
    })
    path = copy_object(tmp_path, 'vmap', version='older')
    metadata = {**objects.load_metadata(path), 'schema_version': '1.0'}
    before = object_bytes(path)
    with pytest.raises(ValueError, match='must preserve schema_version'):
        objects.update_metadata_in_place(path, metadata)
    assert object_bytes(path) == before
