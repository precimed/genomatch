"""Audit the real normalizer's allele correspondence with synthetic genotypes."""
import shutil

import pytest

from genomatch.vmap_restrict_build_compatible import run_bcftools_norm
from utils import run_cmd, write_fasta, write_lines


SENTINELS = ['0/0', '1/1', '0|1', '1|0', './.']


@pytest.fixture
def bcftools(monkeypatch):
    binary = shutil.which('bcftools')
    if binary is None:
        pytest.skip('bcftools not found in PATH')
    monkeypatch.setenv('MATCH_BCFTOOLS', binary)
    return binary


def write_input(tmp_path, sequence, pos, ref, alt):
    fasta, vcf = tmp_path / 'reference.fa', tmp_path / 'input.vcf'
    write_fasta(fasta, {'chr1': sequence})
    write_lines(vcf, [
        '##fileformat=VCFv4.2',
        f'##contig=<ID=chr1,length={len(sequence)}>',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        '\tref_hom\talt_hom\tref_alt_phased\talt_ref_phased\tmissing',
        f'chr1\t{pos}\tvariant\t{ref}\t{alt}\t.\tPASS\t.\tGT\t'
        + '\t'.join(SENTINELS),
    ])
    return fasta, vcf


def read_record(path):
    records = [line.split('\t') for line in path.read_text().splitlines()
               if line and not line.startswith('#')]
    assert len(records) == 1
    record = records[0]
    assert record[0] == 'chr1'
    assert record[2] == 'variant'
    assert record[8] == 'GT'
    return record


def allele_haplotype(sequence, pos, ref, allele):
    start = pos - 1
    assert sequence[start:start + len(ref)] == ref
    return sequence[:start] + allele + sequence[start + len(ref):]


@pytest.mark.parametrize('sequence,pos,ref,alt,expected', [
    pytest.param('CAAAAAGT', 4, 'AA', 'A', (1, 'CA', 'C'), id='homopolymer-deletion'),
    pytest.param('CAAAAAGT', 4, 'A', 'AA', (1, 'C', 'CA'), id='homopolymer-insertion'),
    pytest.param('CATATATG', 2, 'ATAT', 'AT', (1, 'CAT', 'C'), id='repeat-deletion'),
    pytest.param('CATATATG', 4, 'AT', 'ATAT', (1, 'C', 'CAT'), id='repeat-insertion'),
    pytest.param('AGCAGC', 1, 'AGC', 'ATC', (2, 'G', 'T'), id='padded-snv'),
    pytest.param('AGCAGC', 1, 'AGC', 'CCT', (1, 'AGC', 'CCT'), id='mnv'),
    pytest.param('AGCAGC', 1, 'AGC', 'ACCC', (2, 'G', 'CC'), id='delins'),
    pytest.param('AGCAGC', 1, 'AG', 'AGC', (2, 'G', 'GC'), id='both-alleles-match-reference'),
    pytest.param('AAAAC', 1, 'AA', 'A', (1, 'AA', 'A'), id='contig-start'),
])
def test_norm_preserves_sentinels_and_ordered_haplotypes(
    tmp_path, bcftools, sequence, pos, ref, alt, expected,
):
    fasta, vcf = write_input(tmp_path, sequence, pos, ref, alt)
    output = tmp_path / 'normalized.vcf'
    # Exercise the production command builder, rather than restating its flags.
    run_bcftools_norm(vcf, output, fasta, log_path=tmp_path / 'norm.log')
    record = read_record(output)
    actual_pos, actual_ref, actual_alt = int(record[1]), record[3], record[4]
    assert (actual_pos, actual_ref, actual_alt) == expected
    assert record[9:] == SENTINELS
    # Independent of GT handling: each allele must still encode its original
    # full-contig haplotype, even after trimming or shifting through repeats.
    for original, normalized in zip((ref, alt), (actual_ref, actual_alt)):
        assert allele_haplotype(sequence, pos, ref, original) == allele_haplotype(
            sequence, actual_pos, actual_ref, normalized,
        )


@pytest.mark.parametrize('ref,alt,expected', [
    pytest.param('G', 'A', (1, 'A', 'G'), id='snv'),
    pytest.param('AT', 'AGC', (2, 'GC', 'T'), id='unequal-length-alleles'),
])
def test_norm_sentinels_detect_ref_alt_repair(tmp_path, bcftools, ref, alt, expected):
    # Positive control only: production prohibits -c s. The deliberately wrong
    # REF forces repair, while ALT matches the reference sequence in full.
    fasta, vcf = write_input(tmp_path, 'AGCAGC', 1, ref, alt)
    output = tmp_path / 'repaired.vcf'
    result = run_cmd([bcftools, 'norm', '-f', str(fasta), '-c', 's',
                      '-o', str(output), '-Ov', str(vcf)])
    assert result.returncode == 0, result.stderr
    record = read_record(output)
    assert (int(record[1]), record[3], record[4]) == expected
    assert record[9:] == ['1/1', '0/0', '1|0', '0|1', './.']


def test_norm_ref_replacement_does_not_exchange_sentinels(tmp_path, bcftools):
    # Neither submitted allele matches FASTA. Unlike a swap, replacing REF
    # leaves GT unchanged, illustrating why sentinels alone are insufficient.
    fasta, vcf = write_input(tmp_path, 'AGCAGC', 1, 'G', 'C')
    output = tmp_path / 'repaired.vcf'
    result = run_cmd([bcftools, 'norm', '-f', str(fasta), '-c', 's',
                      '-o', str(output), '-Ov', str(vcf)])
    assert result.returncode == 0, result.stderr
    record = read_record(output)
    assert (int(record[1]), record[3], record[4]) == (1, 'A', 'C')
    assert record[9:] == SENTINELS


@pytest.mark.parametrize('ref,alt', [('G', 'A'), ('AT', 'AGC'), ('G', 'C')])
def test_norm_mismatched_ref_fails_instead_of_swapping(tmp_path, bcftools, ref, alt):
    fasta, vcf = write_input(tmp_path, 'AGCAGC', 1, ref, alt)
    with pytest.raises(ValueError, match='Reference allele mismatch'):
        run_bcftools_norm(vcf, tmp_path / 'output.vcf', fasta,
                          log_path=tmp_path / 'norm.log')
