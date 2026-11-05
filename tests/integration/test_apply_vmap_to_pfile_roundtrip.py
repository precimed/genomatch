"""Read projected PGEN hardcalls and dosages through PLINK 2, independently of pgenlib."""
from pathlib import Path

import pytest

from utils import resolve_plink2, run_cmd, run_py, write_object_metadata, write_lines


def export_pfile(plink2: list[str], prefix: Path, output: Path):
    result = run_cmd(plink2 + [
        '--pfile', str(prefix), '--export', 'vcf', 'vcf-dosage=DS-force', '--out', str(output),
    ])
    assert result.returncode == 0, result.stdout + result.stderr
    records = []
    samples = []
    for line in output.with_suffix('.vcf').read_text().splitlines():
        if line.startswith('##'):
            continue
        fields = line.split('\t')
        if line.startswith('#CHROM'):
            samples = fields[9:]
            continue
        sample_fields = [dict(zip(fields[8].split(':'), value.split(':'))) for value in fields[9:]]
        records.append((fields[:5], sample_fields))
    return samples, records


def test_pfile_roundtrip_preserves_hardcalls_dosages_order_and_swaps(tmp_path):
    plink2 = resolve_plink2()
    vcf = tmp_path / 'input.vcf'
    source, output = tmp_path / 'source', tmp_path / 'output'
    vmap = tmp_path / 'target.vmap'
    write_lines(vcf, [
        '##fileformat=VCFv4.2',
        '##contig=<ID=1,length=1000>',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\ts5',
        '1\t100\tsnv\tA\tG\t.\tPASS\t.\tGT:DS\t0/0:0\t0/1:1.0625\t1/1:2\t./.:.\t./.:0.5',
        '1\t200\tdeletion\tAC\tA\t.\tPASS\t.\tGT:DS\t1/1:2\t0/0:0\t0/1:0.9375\t./.:.\t./.:1.5',
        '1\t300\tsnv2\tC\tT\t.\tPASS\t.\tGT:DS\t0/1:1\t1/1:2\t0/0:0\t./.:.\t./.:0.5',
    ])
    result = run_cmd(plink2 + [
        '--vcf', str(vcf), 'dosage=DS', '--make-pgen', '--out', str(source),
    ])
    assert result.returncode == 0, result.stdout + result.stderr
    # Reorder, swap both an SNV and an indel, and reuse the indel source row.
    # Expected REF/ALT are fixed here, independently of what the writer emits.
    target_rows = [
        ['1', '300', 'target2', 'C', 'T', '.', '2', 'swap'],
        ['1', '100', 'target0', 'G', 'A', '.', '0', 'identity'],
        ['1', '200', 'target1', 'AC', 'A', '.', '1', 'swap'],
        ['1', '250', 'reused', 'A', 'AC', '.', '1', 'identity'],
    ]
    write_lines(vmap, ['\t'.join(row) for row in target_rows])
    write_object_metadata(vmap.with_name(vmap.name + '.meta.json'), {
        'object_type': 'variant_map',
        'target': {'genome_build': 'GRCh37', 'contig_naming': 'ncbi'},
        'variants_count': 4,
    })
    result = run_py('apply_vmap_to_pfile.py', '--source-prefix', source,
                    '--vmap', vmap, '--output-prefix', output)
    assert result.returncode == 0, result.stderr
    expected_variants = [
        ['1', '300', '1:300:T:C', 'T', 'C'],
        ['1', '100', '1:100:A:G', 'A', 'G'],
        ['1', '200', '1:200:A:AC', 'A', 'AC'],
        ['1', '250', '1:250:AC:A', 'AC', 'A'],
    ]
    pvar_rows = [line.split('\t') for line in output.with_suffix('.pvar').read_text().splitlines()
                 if line and not line.startswith('#')]
    assert pvar_rows == expected_variants
    source_samples, source_records = export_pfile(plink2, source, tmp_path / 'source_export')
    output_samples, output_records = export_pfile(plink2, output, tmp_path / 'output_export')
    assert source_samples == output_samples == ['s1', 's2', 's3', 's4', 's5']
    assert len(source_records) == 3
    assert len(output_records) == 4
    # Prove that the source fixture retained fractional dosages and missing calls.
    assert float(source_records[0][1][1]['DS']) == pytest.approx(1.0625)
    assert source_records[0][1][3]['GT'] == './.'
    assert source_records[0][1][3]['DS'] == '.'
    assert float(source_records[0][1][4]['DS']) == pytest.approx(0.5)
    for target, (variant, samples), expected_variant in zip(target_rows, output_records, expected_variants):
        assert variant == expected_variant
        source_fields = source_records[int(target[6])][1]
        swapped = target[7] == 'swap'
        for actual, original in zip(samples, source_fields):
            expected_gt = original['GT'].translate(str.maketrans('01', '10')) if swapped else original['GT']
            # Unphased heterozygotes may be emitted in either allele order.
            assert sorted(actual['GT'].split('/')) == sorted(expected_gt.split('/'))
            if original['DS'] == '.':
                assert actual['DS'] == '.'
            else:
                dose = float(original['DS'])
                assert float(actual['DS']) == pytest.approx(2 - dose if swapped else dose, abs=1e-4)
