# Installation

This document covers software/runtime setup for the toolkit. Reference assets, `config.yaml`, and `MATCH_CONFIG` are documented in [DOWNLOADS.md](DOWNLOADS.md).

## Local Conda environment

Recommended local setup is a Conda environment with `bcftools`, the `bcftools-liftover-plugin`, and `samtools`.

```bash
conda create -n match-liftover -c conda-forge -c bioconda \
  python=3.12 bcftools bcftools-liftover-plugin samtools numpy pandas pysam pytest pyyaml pgenlib

conda activate match-liftover
```

Check that the required tools are available:

```bash
bcftools --version
samtools --version
pytest --version
python -c "import yaml; print('PyYAML OK')"
```

## Docker image

The repository root contains a `Dockerfile` that packages the toolkit together with the current runtime dependencies from this document, including Python 3.12, `bcftools`, the `bcftools-liftover-plugin`, and `samtools`.

Build the image from a clean checkout:

```bash
docker build -t genomatch:latest .
```

The image does not include reference FASTA files, chain files, or your `config.yaml`. Mount those from the host and set `MATCH_CONFIG` to the in-container config path:

```bash
docker run --rm -it \
  -v "$PWD":/work \
  -e MATCH_CONFIG=/work/config.yaml \
  genomatch:latest \
  prepare_variants.py --help
```

Any paths referenced from `MATCH_CONFIG` must also be visible inside the container, typically through bind mounts under `/work`. See [DOWNLOADS.md](DOWNLOADS.md) for the required reference layout and config shape.

## `MATCH_BCFTOOLS`

If `bcftools` is not on `PATH`, you may set `MATCH_BCFTOOLS=/path/to/bcftools`.

## External `bcftools +liftover` setups

`liftover_build.py` resolves the `bcftools` executable from `PATH` by default, or from `MATCH_BCFTOOLS` when that environment variable is set.

Current behavior:

- `MATCH_BCFTOOLS` must point to one executable path or executable name
- it is not parsed as a shell command string
- there is currently no separate match-specific override for the `+liftover` plugin path

This means that values like:

```bash
MATCH_BCFTOOLS="singularity exec image.sif bcftools"
```

do not work, because the implementation expects one executable, not a composite command.

If `bcftools +liftover` is available only inside a Singularity container or another external runtime, the recommended pattern is to point `MATCH_BCFTOOLS` at a small wrapper script that executes `bcftools` in that runtime and exports any plugin-path variables it needs, such as `BCFTOOLS_PLUGINS`.

Example wrapper:

```bash
#!/usr/bin/env bash
export BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
exec singularity exec /path/to/container.sif bcftools "$@"
```

Then use:

```bash
export MATCH_BCFTOOLS=/absolute/path/to/bcftools-wrapper.sh
python3 src/genomatch/liftover_build.py \
  --input input.vmap \
  --output lifted.vmap \
  --target-build GRCh38
```

If you use this pattern, make sure the wrapper's runtime can also see the FASTA and chain paths resolved from `MATCH_CONFIG`.

## Typical setup flow

```bash
conda create -n match-liftover -c conda-forge -c bioconda \
  python=3.12 bcftools bcftools-liftover-plugin samtools pytest pyyaml

conda activate match-liftover
bash scripts/download_reference.sh
# create your own config.yaml based on config.example.yaml
# set MATCH_CONFIG as documented in DOWNLOADS.md
pytest tests/test_guess_build.py tests/test_restrict_build_compatible.py tests/test_liftover_build.py
```
