# Installation

This document covers supported end-user installation paths. Development/source-checkout workflows are documented in [DEVELOPMENT.md](DEVELOPMENT.md). Reference assets, `config.yaml`, and `MATCH_CONFIG` are documented in [DOWNLOADS.md](DOWNLOADS.md).

## Option 1: `pip install genomatch`

Create a user runtime environment with the external toolchain first:

```bash
conda create -n genomatch -c conda-forge -c bioconda \
  python=3.12 bcftools bcftools-liftover-plugin samtools pip

conda activate genomatch
```

This path assumes the following runtime dependencies come from that environment or another system package manager:

- `bcftools`
- `bcftools-liftover-plugin`
- `samtools`
- Python 3.12

Then install the Python package:

```bash
python -m pip install genomatch
```

Validate the CLI tools:

```bash
prepare_variants.py --help
project_payload.py --help
```

To update an existing environment to the newest published package version:

```bash
conda activate genomatch
python -m pip install --upgrade genomatch
```

For reference-aware commands, also configure your reference tree and `MATCH_CONFIG` as described in [DOWNLOADS.md](DOWNLOADS.md).

## Option 2: prebuilt Singularity / Apptainer image from GHCR

Pull the published GHCR container image with Apptainer via its `docker://` reference, then run the toolkit inside that image with your `config.yaml` and reference tree mounted in.

Typical use looks like:

```bash
apptainer pull genomatch.sif docker://ghcr.io/precimed/genomatch:<tag>

apptainer exec \
  --home "$PWD":/work \
  --pwd /work \
  --bind /path/to/ref:/ref \
  genomatch.sif \
  prepare_variants.py --help
```

This example assumes:

- your current directory contains the input and output files you want to work with
- the host reference tree lives at `/path/to/ref`
- the reference tree includes `/path/to/ref/config.yaml`
- that config uses paths relative to itself, such as `ucsc/...` and `chain/...`

With that layout, the container finds `/ref/config.yaml` automatically and resolves the rest of the reference assets under `/ref/...`. If you use a different layout, set `MATCH_CONFIG` explicitly and make sure every path referenced from that config exists inside the container at the same path written in the file. See [DOWNLOADS.md](DOWNLOADS.md) for the required reference layout and config shape.
