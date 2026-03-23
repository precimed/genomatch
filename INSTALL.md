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

## Option 2: prebuilt Singularity / Apptainer image via ORAS

Pull the published image from GHCR via its release-specific ORAS reference, then run the toolkit inside that image with your `config.yaml` and reference tree mounted in.

Typical use looks like:

```bash
apptainer pull genomatch.sif oras://ghcr.io/<publisher>/<package>:<tag>

apptainer exec \
  --bind "$PWD":/work \
  --env MATCH_CONFIG=/work/config.yaml \
  genomatch.sif \
  prepare_variants.py --help
```

Any paths referenced from `MATCH_CONFIG` must also be visible inside the container. See [DOWNLOADS.md](DOWNLOADS.md) for the required reference layout and config shape.
