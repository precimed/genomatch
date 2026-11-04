import json
from pathlib import Path

import pytest

from test_restrict_build_compatible import restrict_env_with_fake_norm
from utils import read_tsv, run_py_with_env, write_lines


def wrapper_args(tmp_path, sharded, rows):
    if sharded:
        write_lines(tmp_path/'raw.1.bim', rows)
        return ['prepare_variants_sharded.py', '--input', tmp_path/'raw.@.bim', '--prefix', tmp_path/'work.@', '--output', tmp_path/'final', '--input-format', 'bim', '--src-build', 'GRCh37', '--dst-build', 'GRCh37'], tmp_path/'work.1'
    write_lines(tmp_path/'raw.bim', rows)
    return ['prepare_variants.py', '--input', tmp_path/'raw.bim', '--prefix', tmp_path/'work', '--output', tmp_path/'final', '--input-format', 'bim', '--src-build', 'GRCh37', '--dst-build', 'GRCh37'], tmp_path/'work'


@pytest.mark.parametrize('sharded', [False, True])
@pytest.mark.parametrize('options,expected', [
    ([], {'--allow-strand-flips'}),
    (['--include-indels'], {'--include-indels', '--allow-strand-flips', '--allow-indel-swaps', '--norm-indels'}),
    (['--include-indels', '--no-allow-indel-swaps'], {'--include-indels', '--allow-strand-flips', '--norm-indels'}),
    (['--include-indels', '--no-norm-indels'], {'--include-indels', '--allow-strand-flips', '--allow-indel-swaps'}),
    (['--include-indels', '--no-norm-indels', '--no-allow-indel-swaps', '--no-allow-strand-flips'], {'--include-indels'}),
    (['--allow-indel-swaps', '--norm-indels'], {'--allow-strand-flips', '--allow-indel-swaps', '--norm-indels'}),
])
def test_actual_wrapper_forwarding_and_input_filter_placement(tmp_path, sharded, options, expected):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={}, grch37_sequence='AGCAGC')
    args, prefix = wrapper_args(tmp_path, sharded, ['1\tsnv\t0\t4\tC\tA', '1\tdeletion\t0\t1\tA\tAG'])
    result = run_py_with_env(args[0], env, *args[1:], *options)
    assert result.returncode == 0, result.stderr
    command = next(line for line in result.stderr.splitlines() if '+ restrict_build_compatible.py ' in line)
    tokens = set(command.split())
    controls = {'--include-indels', '--allow-strand-flips', '--allow-indel-swaps', '--norm-indels'}
    assert tokens & controls == expected
    assert len(read_tsv(Path(str(prefix)+'.imported.vmap'))) == 2
    assert len(read_tsv(tmp_path/'final.vmap')) == (2 if '--include-indels' in expected else 1)
    if '--include-indels' not in expected:
        assert read_tsv(Path(str(prefix)+'.build_compatible.vmap.qc.tsv'))[-1][-1] == 'non_snv_input'
        if '--norm-indels' in expected:
            assert 'no-ops' in result.stderr
    if sharded:
        command = next(line for line in result.stderr.splitlines() if '+ prepare_variants.py ' in line)
        assert ('--norm-indels' in command) == ('--norm-indels' in options)
        assert ('--no-norm-indels' in command) == ('--no-norm-indels' in options)


@pytest.mark.parametrize('sharded', [False, True])
def test_force_without_indels_cleans_all_normalization_diagnostics(tmp_path, sharded):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={}, grch37_sequence='AGC')
    args, prefix = wrapper_args(tmp_path, sharded, ['1\tdeletion\t0\t1\tA\tAG'])
    first = run_py_with_env(args[0], env, *args[1:], '--include-indels')
    assert first.returncode == 0, first.stderr
    namespace = Path(str(prefix)+'.build_compatible.vmap')
    assert list(tmp_path.glob(namespace.name+'.bcftools_norm*'))
    second = run_py_with_env(args[0], env, *args[1:], '--force')
    assert second.returncode == 0, second.stderr
    assert not list(tmp_path.glob(namespace.name+'.bcftools_norm*'))
    assert read_tsv(tmp_path/'final.vmap') == []


@pytest.mark.parametrize('sharded', [False, True])
def test_all_filtered_input_succeeds_across_builds(tmp_path, sharded):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={'__fail__': 'must not invoke'}, grch37_sequence='AGC')
    args, prefix = wrapper_args(tmp_path, sharded, ['1\tdeletion\t0\t1\tA\tAG'])
    args[-1] = 'GRCh38'
    # Resolve a valid edge even though no records require a plugin call.
    config = Path(env['MATCH_CONFIG'])
    chain = tmp_path/'chain'
    chain.write_text('unused')
    config.write_text(config.read_text().replace('liftover: []', f'liftover:\n  - source: GRCh37\n    target: GRCh38\n    chain: {chain}'))
    result = run_py_with_env(args[0], env, *args[1:])
    assert result.returncode == 0, result.stderr
    for stage in ('build_compatible', 'lifted', 'sorted'):
        assert read_tsv(Path(str(prefix)+f'.{stage}.vmap')) == []
    metadata = json.loads((tmp_path/'final.vmap.meta.json').read_text())
    assert metadata['target'] == {'genome_build': 'GRCh38', 'contig_naming': 'ncbi'}
    assert metadata['variants_count'] == 0
    assert metadata['derived_from'] == str(tmp_path/('raw.@.bim' if sharded else 'raw.bim'))


@pytest.mark.parametrize('sharded', [False, True])
def test_empty_unknown_source_requests_explicit_build(tmp_path, sharded):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={})
    args, _ = wrapper_args(tmp_path, sharded, [])
    i = args.index('--src-build')
    del args[i:i+2]
    result = run_py_with_env(args[0], env, *args[1:])
    assert result.returncode != 0
    assert '--src-build' in result.stderr


@pytest.mark.parametrize('all_empty', [False, True])
def test_sharded_empty_groups_and_metadata(tmp_path, all_empty):
    env, _ = restrict_env_with_fake_norm(tmp_path, behavior={}, grch37_sequence='AGC')
    args, prefix = wrapper_args(tmp_path, True, ['1\tdropped\t0\t1\tA\tAG'])
    write_lines(tmp_path/'raw.2.bim', ['2\tsecond\t0\t1\t' + ('A\tAG' if all_empty else 'C\tA')])
    result = run_py_with_env(args[0], env, *args[1:])
    assert result.returncode == 0, result.stderr
    rows = read_tsv(tmp_path/'final.vmap')
    assert len(rows) == (0 if all_empty else 1)
    metadata = json.loads((tmp_path/'final.vmap.meta.json').read_text())
    assert metadata['variants_count'] == len(rows)
    assert metadata['derived_from'] == str(tmp_path/'raw.@.bim')
