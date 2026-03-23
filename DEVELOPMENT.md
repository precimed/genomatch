# Development

This document covers local development, test execution, and maintainer-facing packaging notes. End-user installation is documented in [INSTALL.md](INSTALL.md). Reference assets and `MATCH_CONFIG` setup are documented in [DOWNLOADS.md](DOWNLOADS.md).

## Local development environment

The canonical developer workflow is:

1. clone this repository locally
2. create or activate a local `match-liftover` Conda environment
3. install the local checkout into that environment with `pip install -e .`
4. run the CLI tools from that environment, including on real data

The expected runtime environment is Python 3.12 plus external `bcftools`, the `bcftools-liftover-plugin`, and `samtools`.

Create and activate the environment:

```bash
conda create -n match-liftover -c conda-forge -c bioconda \
  python=3.12 bcftools bcftools-liftover-plugin samtools numpy pandas pysam pytest pyyaml pgenlib

conda activate match-liftover
python -m pip install -e .
```

Validate the core toolchain:

```bash
bcftools --version
samtools --version
pytest --version
python -c "import yaml; print('PyYAML OK')"
prepare_variants.py --help
```

## Running local workflows

After `pip install -e .`, invoke the tools directly:

```bash
prepare_variants.py --help
intersect_variants.py --help
project_payload.py --help
```

This is the preferred way to exercise real-data workflows during development. Keep source-path execution like `python src/genomatch/...` out of the normal workflow.

## Running tests

The canonical test environment is the local `match-liftover` Conda environment above with the repository installed editable.

Run the narrowest relevant tests first:

```bash
pytest tests/test_guess_build.py
pytest tests/test_restrict_build_compatible.py
pytest tests/test_liftover_build.py
pytest tests/test_package_install.py
```

Run the full suite before merging code changes:

```bash
pytest tests
```

See [TESTS.md](TESTS.md) for the authoritative test contract.

## Building runtime images

Build the local Docker image from the repository root:

```bash
docker build -t genomatch:latest .
```

If you need a local Singularity/Apptainer image derived from that Docker image, build it outside the repository workflow using your site-standard tooling. A typical local pattern is:

```bash
apptainer build genomatch.sif docker-daemon://genomatch:latest
```

## Maintainer notes

This file is the place to document maintainer workflows for:

- building and publishing container artifacts
- preparing and publishing PyPI releases

Keep end-user install instructions in [INSTALL.md](INSTALL.md). Keep reference setup in [DOWNLOADS.md](DOWNLOADS.md). Record maintainer-only Docker, Singularity/Apptainer, ORAS, and PyPI release commands here as those workflows are finalized.
