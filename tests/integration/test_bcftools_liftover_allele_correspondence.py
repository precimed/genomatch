"""Check real-plugin annotations against genotypes and controlled haplotypes."""
import random

import pandas as pd
import pytest

from genomatch import vmap_liftover_build as liftover
from real_liftover_helpers import resolve_bcftools_with_liftover
from utils import write_fasta, write_lines


SENTINELS = ['0/0', '1/1', '0|1', '1|0', './.']
SWAPPED_SENTINELS = ['1/1', '0/0', '1|0', '0|1', './.']
UPSTREAM_OPS = ['identity', 'swap', 'flip', 'flip_swap']
COMPOSED_OPS = {
    'identity': ['identity', 'swap', 'flip', 'flip_swap'],
    'swap': ['swap', 'identity', 'flip_swap', 'flip'],
    'flip': ['flip', 'flip_swap', 'identity', 'swap'],
    'flip_swap': ['flip_swap', 'flip', 'swap', 'identity'],
}


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def haplotype(sequence, pos, ref, allele):
    start = pos - 1
    assert sequence[start:start + len(ref)] == ref
    return sequence[:start] + allele + sequence[start + len(ref):]


def run_fixture(tmp_path, monkeypatch, *, kind, reverse, swap=False, new_ref=False):
    tmp_path.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv('MATCH_BCFTOOLS', resolve_bcftools_with_liftover())
    rng = random.Random(71)
    source = ''.join(rng.choice('ACGT') for _ in range(600))
    pos = 201
    if kind.startswith('palindromic-'):
        source = source[:pos - 1] + ('A' if kind == 'palindromic-at' else 'C') + source[pos:]
    elif kind == 'repeat-deletion':
        # ACA/A at 205 must left-align across the CA repeat to GCA/G at
        # 201. The terminal T also fixes the reverse-strand repeat boundary.
        source = source[:200] + 'GCACACACT' + source[209:]
        pos = 205
    ref_length = 1 if kind in ('snv', 'insertion', 'palindromic-at', 'palindromic-cg') else 3
    ref = source[pos - 1:pos - 1 + ref_length]
    alt = {
        'snv': ref.translate(str.maketrans('ACGT', 'CGTA')),
        'mnv': ref.translate(str.maketrans('ACGT', 'CGTA')),
        'deletion': ref[0],
        'repeat-deletion': ref[0],
        'insertion': ref + 'AC',
        'delins': ref[0] + 'CCCCC',
        'palindromic-at': 'T',
        'palindromic-cg': 'G',
    }[kind]
    target_allele = alt if swap else ref
    if new_ref:
        assert kind == 'snv'
        target_allele = next(base for base in 'ACGT' if base not in (ref, alt))
    target = haplotype(source, pos, ref, target_allele)
    if reverse:
        target = reverse_complement(target)

    source_fasta, target_fasta = tmp_path / 'source.fa', tmp_path / 'target.fa'
    for fasta, sequence in ((source_fasta, source), (target_fasta, target)):
        write_fasta(fasta, {'chr1': sequence})
    chain = tmp_path / 'mapping.chain'
    # For unequal lengths, map the shared first base and downstream sequence;
    # the replaced suffix lies in a chain gap on either or both references.
    blocks = ['600'] if len(target_allele) == len(ref) else [
        f'{pos} {len(ref) - 1} {len(target_allele) - 1}',
        str(len(source) - (pos - 1 + len(ref))),
    ]
    write_lines(chain, [
        f'chain 1000 chr1 {len(source)} + 0 {len(source)} '
        f'chr1 {len(target)} {"-" if reverse else "+"} 0 {len(target)} 1',
        *blocks, '',
    ])

    prepared = liftover.prepare_liftover_table(pd.DataFrame({
        'row_id': ['row0'], 'chrom': ['chr1'], 'pos': [str(pos)],
        'id': ['variant'], 'a1': [alt], 'a2': [ref],
        'upstream_allele_op': ['identity'], 'source_shard': ['original'],
        'source_index': [91],
    }), is_vmap=True)
    input_vcf, output_vcf = tmp_path / 'input.vcf', tmp_path / 'output.vcf'
    lookup = liftover.write_temp_vcf(input_vcf, prepared, is_vmap=True)
    # Add test-only samples to the production-emitted REF/ALT assignments.
    lines = []
    for line in input_vcf.read_text().splitlines():
        if line.startswith('#CHROM'):
            lines.extend([
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                line + '\tFORMAT\tref_hom\talt_hom\tref_alt_phased\talt_ref_phased\tmissing',
            ])
        elif line.startswith('#'):
            lines.append(line)
        else:
            lines.append(line + '\tGT\t' + '\t'.join(SENTINELS))
    write_lines(input_vcf, lines)
    liftover.run_bcftools_liftover(input_vcf, output_vcf, source_fasta, target_fasta, chain)
    records = [line.split('\t') for line in output_vcf.read_text().splitlines()
               if line and not line.startswith('#')]
    assert len(records) == 1
    record = records[0]
    assert record[0] == 'chr1'
    assert record[2] == 'row0'
    assert record[8] == 'GT'
    return source, target, pos, ref, alt, output_vcf, lookup, record


@pytest.mark.parametrize('kind', [
    'snv', 'mnv', 'deletion', 'insertion', 'delins', 'repeat-deletion',
    'palindromic-at', 'palindromic-cg',
])
@pytest.mark.parametrize('reverse', [False, True], ids=['forward', 'reverse'])
@pytest.mark.parametrize('swap', [False, True], ids=['same-ref', 'exchanged-ref'])
def test_liftover_tags_match_sentinels_haplotypes_and_composed_ops(
    tmp_path, monkeypatch, kind, reverse, swap,
):
    source, target, pos, ref, alt, output, lookup, record = run_fixture(
        tmp_path, monkeypatch, kind=kind, reverse=reverse, swap=swap,
    )
    # Expected annotations come from the constructed references and chain,
    # independently of the toolkit's INFO decoder.
    info = set(record[7].split(';'))
    assert ('FLIP' in info) == reverse
    assert {token for token in info if token.startswith('SWAP')} == ({'SWAP=1'} if swap else set())
    assert record[9:] == (SWAPPED_SENTINELS if swap else SENTINELS)
    target_pos, target_ref, target_alt = int(record[1]), record[3], record[4]
    if kind == 'repeat-deletion':
        # Assert a concrete shift as well as haplotype equivalence, so a
        # future fixture edit cannot accidentally remove left-alignment coverage.
        expected_pos, expected_ref, expected_alt = (392, 'AGT', 'A') if reverse else (201, 'GCA', 'G')
        assert target_pos == expected_pos
        assert (target_ref, target_alt) == ((expected_alt, expected_ref) if swap else (expected_ref, expected_alt))
    elif kind == 'insertion' and reverse and not swap:
        # The anchor must precede the mirrored source base; local repeats may
        # move it farther left than the mandatory one-base reanchoring.
        assert target_pos < len(source) - pos + 1
    target_alleles = (target_ref, target_alt)
    for index, allele in enumerate((ref, alt)):
        expected = haplotype(source, pos, ref, allele)
        if reverse:
            expected = reverse_complement(expected)
        corresponding = target_alleles[1 - index if swap else index]
        assert haplotype(target, target_pos, target_ref, corresponding) == expected

    local_op = ('flip_swap' if swap else 'flip') if reverse else ('swap' if swap else 'identity')
    for upstream, expected_op in zip(UPSTREAM_OPS, COMPOSED_OPS[local_op]):
        table, statuses = liftover.parse_lifted_vcf(
            output, lookup.assign(upstream_allele_op=upstream),
            'ncbi', 'GRCh37', 'GRCh38', is_vmap=True,
        )
        assert statuses == {'row0': 'lifted'}
        assert table.frame[['a1', 'a2']].values.tolist() == [[target_alt, target_ref]]
        assert table.frame['allele_op'].tolist() == [expected_op]
        assert table.frame[['source_shard', 'source_index']].values.tolist() == [['original', 91]]


@pytest.mark.parametrize('kind', ['palindromic-at', 'palindromic-cg'])
def test_palindromic_strings_do_not_determine_genotype_exchange(tmp_path, monkeypatch, kind):
    *_, flipped = run_fixture(
        tmp_path / 'flip', monkeypatch, kind=kind, reverse=True,
    )
    *_, swapped = run_fixture(
        tmp_path / 'swap', monkeypatch, kind=kind, reverse=False, swap=True,
    )
    assert flipped[3:5] == swapped[3:5]
    assert set(flipped[7].split(';')) == {'FLIP'}
    assert set(swapped[7].split(';')) == {'SWAP=1'}
    assert flipped[9:] == SENTINELS
    assert swapped[9:] == SWAPPED_SENTINELS


@pytest.mark.parametrize('reverse', [False, True], ids=['forward', 'reverse'])
def test_liftover_new_reference_is_not_coerced_into_biallelic_op(tmp_path, monkeypatch, reverse):
    _, _, _, _, _, output, lookup, record = run_fixture(
        tmp_path, monkeypatch, kind='snv', reverse=reverse, new_ref=True,
    )
    info = set(record[7].split(';'))
    assert 'SWAP=-1' in info
    assert ('FLIP' in info) == reverse
    assert ',' in record[4]
    # A third allele was introduced: the source REF/ALT now have indices 1/2.
    assert record[9:] == ['1/1', '2/2', '1|2', '2|1', './.']
    table, statuses = liftover.parse_lifted_vcf(
        output, lookup, 'ncbi', 'GRCh37', 'GRCh38', is_vmap=True,
    )
    assert table.frame.empty
    assert statuses == {'row0': 'unsupported_representation'}
