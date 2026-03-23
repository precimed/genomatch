# Installation

This document covers supported end-user installation paths. Development/source-checkout workflows are documented in [DEVELOPMENT.md](DEVELOPMENT.md). Reference assets, `config.yaml`, and `MATCH_CONFIG` are documented in [DOWNLOADS.md](DOWNLOADS.md).

## Option 1: `pip install genomatch`

This path assumes the non-Python runtime dependencies are already available through your environment management channel:

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
