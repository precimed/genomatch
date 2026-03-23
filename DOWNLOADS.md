# Reference Downloads and Config

This document covers the reference assets required by the toolkit and the `config.yaml` / `MATCH_CONFIG` setup that points tools at those assets.

## Current reference model

Current reference-aware tools use UCSC-style internal reference assets:

- `guess_build.py`
- `restrict_build_compatible.py`
- `liftover_build.py`

These tools resolve UCSC-style internal FASTA assets from the user-provided file pointed to by `MATCH_CONFIG`. `liftover_build.py` also resolves UCSC chain files from the same config.

User inputs may still declare either UCSC or NCBI-style contig naming as per [spec/contigs-and-metadata.md](spec/contigs-and-metadata.md).

Optional FASTA assets with NCBI-style contig naming may still be kept locally under `ref/ncbi/`, but they are not used by the current implementation.

## Reference layout

Reference assets are intentionally not tracked in git. The recommended local layout is:

```text
ref/
  ucsc/
    GRCh37/
      hg19.p13.plusMT.no_alt_analysis_set.fa
      hg19.p13.plusMT.no_alt_analysis_set.fa.fai
    GRCh38/
      GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
      GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
  ncbi/
    GRCh37/
      human_g1k_v37.fasta
      human_g1k_v37.fasta.fai
    GRCh38/
      GRCh38.primary.fa
      GRCh38.primary.fa.fai
  chain/
    hg19ToHg38.over.chain.gz
    hg38ToHg19.over.chain.gz
```

`ref/ucsc/` is the active reference set for the current toolkit. `ref/ncbi/` is optional and currently unused by `guess_build.py`, `restrict_build_compatible.py`, and `liftover_build.py`.

Create directories:

```bash
mkdir -p ref/ucsc/GRCh37 ref/ucsc/GRCh38
mkdir -p ref/ncbi/GRCh37 ref/ncbi/GRCh38
mkdir -p ref/chain
```

## Automatic download helper

Fetch the required UCSC assets automatically with:

```bash
bash scripts/download_reference.sh
```

`scripts/download_reference.sh` downloads only the required UCSC FASTA and chain files, skips files that already exist, and runs `samtools faidx` only when the FASTA index is missing.

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

## config.yaml shape

Reference-aware tools do not accept CLI overrides for FASTA or chain assets. They resolve everything from config plus object metadata.

The current config shape is:

```yaml
references:
  ucsc:
    GRCh37:
      fasta: ucsc/GRCh37/hg19.p13.plusMT.no_alt_analysis_set.fa
    GRCh38:
      fasta: ucsc/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
  ncbi:
    GRCh37:
      fasta: ncbi/GRCh37/human_g1k_v37.fasta
    GRCh38:
      fasta: ncbi/GRCh38/GRCh38.primary.fa
chain:
  hg19ToHg38: chain/hg19ToHg38.over.chain.gz
  hg38ToHg19: chain/hg38ToHg19.over.chain.gz
```

Notes:

- Paths may be absolute or relative to the config file.
- The recommended layout is to place `config.yaml` inside the root of the `ref/` tree, so relative paths like `ucsc/...` and `chain/...` resolve naturally.
- `references.ucsc.*` is the active internal reference set for current reference-aware tools.
- `references.ncbi.*` is optional and currently unused by the implementation.
- `chain.*` is the active chain-file configuration block and matches the on-disk `ref/chain/` layout.
- Current reference-aware tools do not treat UCSC and NCBI-style FASTA assets as interchangeable.

Save this as `ref/config.yaml`, then export:

```bash
export MATCH_CONFIG=/path/to/ref/config.yaml
```

If `MATCH_CONFIG` is unset, reference-aware tools also fall back to `/ref/config.yaml`, which is convenient for container runs that bind the host reference tree to `/ref`.
