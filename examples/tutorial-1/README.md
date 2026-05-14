# Tutorial 1 Fixture

This directory contains a small, deterministic fixture for the Tutorial 1 workflow.

- `reference.grch38.vcf` is the GRCh38 target universe and uses UCSC contig labels.
- `study.grch37.tsv` is a GRCh37 summary-statistics payload with non-canonical column names.
- `cohort.{21,22,X}.{bed,bim,fam}` is a sharded GRCh38 PLINK 1 BFILE payload with PLINK contig labels.

The fixture is intentionally messy. Warnings during preparation and projection are expected; inspect `manifest.tsv` for the reason each special row exists.

## Effective Membership Counts

These counts are after preparation QC and the tutorial's study-side strand-ambiguous filter, so intentionally duplicated study rows and A/T or C/G study SNPs are excluded from study membership.

`reference=200  study=94  cohort=150`

| subset | count |
| --- | ---: |
| reference + study + cohort | 52 |
| reference + study only | 13 |
| reference + cohort only | 52 |
| study + cohort only | 11 |
| reference only | 83 |
| study only | 18 |
| cohort only | 35 |

```text
                  study
             ref+study       study+cohort
reference        3-way             cohort
             ref+cohort
```

## issue_tags

- `allele_swap`: source alleles are reversed relative to the reference-oriented target row
- `cohort_only`: present only in the cohort after preparation
- `duplicate_rsid`: same rsID appears at more than one coordinate
- `duplicate_target`: duplicate chrom/position/allele identity within one input
- `indel`: row is an insertion/deletion represented with A/C/G/T alleles
- `invalid_allele`: row contains a non-A/C/G/T allele and is dropped during import
- `invalid_stat`: invalid clean summary-statistic values become missing in projection
- `liftover_collision`: two study rows collapse to the same target identity after preparation
- `long_allele`: allele exceeds the tutorial max allele length and is dropped during import
- `lowercase_allele`: lowercase alleles are accepted and normalized
- `missing_genotype`: PLINK genotype row includes one or more missing calls
- `mixed_contig_label`: study rows intentionally mix accepted contig-label spellings
- `multiallelic`: VCF row has multiple ALT alleles and is dropped by the importer
- `placeholder_id`: input ID is '.'
- `reference_only`: present only in the target-universe VCF
- `strand_ambiguous`: SNP has an A/T or C/G allele pair
- `strand_flip`: source alleles are reverse-complemented relative to the reference-oriented target row
- `study_only`: present only in the study after preparation
- `unsorted`: input rows are intentionally not sorted by coordinate
- `x_nonpar`: chrX non-PAR locus
- `x_par`: chrX PAR locus

See ../../docs/tutorial-1.md for the end-to-end commands.
