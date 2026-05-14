# Tutorial 1: Reference-Universe Projection

This tutorial uses the checked-in fixture under `examples/tutorial-1/` to show the normal prepare -> restrict -> project workflow on small but realistic data.

The reference VCF defines the target variant universe. The study summary statistics and the cohort BFILE are payloads: each keeps its own source-row provenance and is projected independently through its own restricted `.vmap`.

## Inputs

```text
examples/tutorial-1/
  reference.grch38.vcf
  study.grch37.tsv
  study.grch37.meta.yaml
  cohort.21.{bed,bim,fam}
  cohort.22.{bed,bim,fam}
  cohort.X.{bed,bim,fam}
```

Warnings are expected. The fixture includes allele swaps, strand flips, ambiguous SNPs, indels, duplicate IDs, invalid alleles, long alleles, a multiallelic VCF row, liftover collisions, invalid summary-statistic values, and missing genotype calls. See `examples/tutorial-1/manifest.tsv` for row-level intent.

## Setup

Reference-aware steps require `MATCH_CONFIG` to point at local GRCh37/GRCh38 FASTA and chain assets. See `docs/install.md` and `docs/downloads.md`.

```bash
export MATCH_CONFIG=/path/to/ref/config.yaml
mkdir -p work/tutorial-1 out/tutorial-1
```

## 1. Prepare the Reference Universe

```bash
prepare_variants.py \
  --input examples/tutorial-1/reference.grch38.vcf \
  --input-format vcf \
  --dst-build GRCh38 \
  --dst-contig-naming ncbi \
  --output work/tutorial-1/reference
```

This writes `work/tutorial-1/reference.vmap`. It is the target universe used for membership.

## 2. Prepare the Study Payload

```bash
prepare_variants.py \
  --input-format sumstats \
  --sumstats-metadata examples/tutorial-1/study.grch37.meta.yaml \
  --dst-build GRCh38 \
  --dst-contig-naming ncbi \
  --drop-strand-ambiguous \
  --output work/tutorial-1/study
```

The study starts on GRCh37 and is lifted to GRCh38 during preparation. The metadata uses non-canonical column names so the clean projection can derive fields such as `BETA`, `Z`, effective `N`, and allele frequencies.

This fixture deliberately includes many allele flips and several A/T or C/G SNPs. In such summary-statistics data, strand orientation for those ambiguous SNPs cannot be resolved confidently, so the tutorial drops them from the study payload. In real datasets, if alleles are known to be reported on the positive strand for a known reference genome, users usually do not need `--drop-strand-ambiguous`.

## 3. Prepare the Sharded Cohort Payload

```bash
prepare_variants_sharded.py \
  --input examples/tutorial-1/cohort.@.bim \
  --input-format bim \
  --prefix work/tutorial-1/cohort.@ \
  --output work/tutorial-1/cohort \
  --dst-build GRCh38 \
  --dst-contig-naming ncbi \
  --shards 21,22,X
```

This imports the PLINK 1 BIM shards, prepares each chromosome group, and writes one concatenated `work/tutorial-1/cohort.vmap`.

## 4. Restrict Each Payload to the Reference Universe

```bash
restrict_vmap.py \
  work/tutorial-1/study.vmap \
  work/tutorial-1/reference.vmap \
  --output work/tutorial-1/study.reference.vmap

restrict_vmap.py \
  work/tutorial-1/cohort.vmap \
  work/tutorial-1/reference.vmap \
  --output work/tutorial-1/cohort.reference.vmap
```

The restricted `.vmap` files preserve payload-specific row provenance. The reference controls membership; it is not applied as a payload here.

## 5. Project the Study

```bash
project_payload.py \
  --input-format sumstats-clean \
  --sumstats-metadata examples/tutorial-1/study.grch37.meta.yaml \
  --vmap work/tutorial-1/study.reference.vmap \
  --output out/tutorial-1/study.reference.tsv \
  --use-af-inference
```

The output is clean, canonical summary statistics in the GRCh38 reference universe. Invalid numeric values in retained payload rows are written as missing values where the clean pipeline permits that.

## 6. Project the Cohort

```bash
project_payload.py \
  --input examples/tutorial-1/cohort.@.bim \
  --input-format bfile \
  --vmap work/tutorial-1/cohort.reference.vmap \
  --output out/tutorial-1/cohort.@ \
  --skip-ploidy-check
```

This writes sharded PLINK 1 BFILE outputs under `out/tutorial-1/cohort.<contig>`. The tutorial uses `--skip-ploidy-check` because the fixture focuses on variant mapping and includes compact toy chrX genotypes.

## Regenerating the Fixture

The checked-in fixture is generated from local reference assets:

```bash
MATCH_CONFIG=/path/to/ref/config.yaml \
python scripts/generate_tutorial_1_fixture.py
```

The generator rewrites every file under `examples/tutorial-1/`, including `README.md` and `manifest.tsv`, and validates the tutorial command path before completing.
