# Reference Downloads and Config

This document covers the reference assets required by the toolkit and the `config.yaml` / `MATCH_CONFIG` setup that points tools at those assets.

## Current reference model

Current reference-aware tools use UCSC-style internal reference assets:

- `guess_build.py`
- `restrict_build_compatible.py`
- `liftover_build.py`

These tools resolve UCSC-style internal FASTA assets from the user-provided file pointed to by `MATCH_CONFIG`. `liftover_build.py` also resolves directed UCSC chain-file edges from the same config.

User inputs may still declare either UCSC or NCBI-style contig naming as per [spec/contigs-and-metadata.md](spec/contigs-and-metadata.md).

Optional FASTA assets with NCBI-style contig naming may still be kept locally under `ref/ncbi/`, but they are not used by the current implementation.

## Config setup

Reference assets are intentionally not tracked in git. Use
[`config.example.yaml`](config.example.yaml) as the source of truth for the
expected `config.yaml` shape and relative reference layout.

Reference-aware tools do not accept CLI overrides for FASTA or chain assets.
They resolve everything from config plus object metadata. Save an adapted copy
of `config.example.yaml` as `ref/config.yaml`, then export:

```bash
export MATCH_CONFIG=/path/to/ref/config.yaml
```

If `MATCH_CONFIG` is unset, reference-aware tools also fall back to
`/ref/config.yaml`, which is convenient for container runs that bind the host
reference tree to `/ref`. Paths in the config may be absolute or relative to
the config file.

## Automatic download helper

Fetch the required UCSC assets automatically with:

```bash
bash scripts/download_reference.sh
```

`scripts/download_reference.sh` downloads the GRCh37/GRCh38 UCSC FASTA and chain files, skips files that already exist, and runs `samtools faidx` only when the FASTA index is missing. T2T-CHM13v2.0 assets are optional unless a T2T workflow is requested; download them manually from the sources below and add them to `config.yaml` when needed.

## Manual UCSC downloads

These are the required FASTA files for the current implementation.

### GRCh37 UCSC-style FASTA

The recommended GRCh37 asset is UCSC's `hg19.p13.plusMT.no_alt_analysis_set.fa.gz`. This is the closest GRCh37 analogue to the GRCh38 `no_alt_analysis_set` file and is the correct asset for the current toolkit's internal GRCh37 reference path.

```bash
wget -O- \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
  | gzip -d > ref/ucsc/GRCh37/hg19.p13.plusMT.no_alt_analysis_set.fa

samtools faidx ref/ucsc/GRCh37/hg19.p13.plusMT.no_alt_analysis_set.fa
```

### GRCh38 UCSC-style FASTA

```bash
wget -O- \
  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
  | gzip -d > ref/ucsc/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

samtools faidx ref/ucsc/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

## Manual UCSC chain downloads

These are required only for `liftover_build.py`. Same-build reference-compatible restriction with `restrict_build_compatible.py` does not require chain files.

```bash
wget -O ref/chain/hg19ToHg38.over.chain.gz \
  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

wget -O ref/chain/hg38ToHg19.over.chain.gz \
  http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

## Optional T2T-CHM13v2.0 downloads

T2T-CHM13v2.0 support uses the UCSC `hs1`/CHM13v2.0 primary-contig surface (`chr1` through `chr22`, `chrX`, `chrY`, `chrM`). The full `chm13v2.0.fa` includes `chrY` sequence from HG002/NA24385; the toolkit uses that full primary-contig FASTA because the primary-contig contract includes `Y`.

Download starting points:

- FASTA: `https://42basepairs.com/browse/s3/human-pangenomics/T2T/CHM13/assemblies?file=chm13v2.0.fa`
- Chain files: `https://42basepairs.com/browse/s3/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo`

Upstream chain sources may provide plain `.chain` files or differently named `.over.chain.gz` files. Save or gzip them under names that match your `liftover` config entries, then inspect chain headers to confirm source and target contigs are UCSC-style.

## Optional NCBI-style FASTA assets

These are not used by the current reference-aware implementation, but they may still be useful for inspection, fixture generation, or future utilities.

### Optional GRCh37 FASTA with NCBI-style contig naming

`human_g1k_v37.fasta` uses NCBI-style chromosome labels (`1`, `2`, ..., `X`, `Y`, `MT`). It is not interchangeable with the UCSC-style GRCh37 FASTA above.

```bash
wget -O- \
  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz \
  | gzip -d > ref/ncbi/GRCh37/human_g1k_v37.fasta

samtools faidx ref/ncbi/GRCh37/human_g1k_v37.fasta
```

### Optional GRCh38 FASTA with NCBI-style contig naming

```bash
wget -O- \
  https://ftp.ensembl.org/pub/grch38/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  | gzip -d > ref/ncbi/GRCh38/GRCh38.primary.fa

samtools faidx ref/ncbi/GRCh38/GRCh38.primary.fa
```
