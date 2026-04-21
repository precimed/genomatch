# Development

This document covers local development, test execution, and maintainer-facing packaging notes. End-user installation is documented in [INSTALL.md](INSTALL.md). Reference assets and `MATCH_CONFIG` setup are documented in [DOWNLOADS.md](DOWNLOADS.md).

## Table of contents

- [Local development environment](#local-development-environment)
- [Running local workflows](#running-local-workflows)
- [Running tests](#running-tests)
- [Building runtime images](#building-runtime-images)
- [Pre-release validation](#pre-release-validation)
- [Actual release](#actual-release)

## Local development environment

The canonical developer workflow is:

1. clone this repository locally
2. create or activate a local `match-liftover` Conda environment
3. install the local checkout into that environment with `pip install -e .`
4. run the CLI tools from that environment, including on real data

The expected runtime environment is Python 3.12 plus external `bcftools`, the `bcftools-liftover-plugin`, and `samtools`.
Development adds `plink`/`plink2` for payload/apply test coverage, plus Python packaging/test tools.

Create and activate the environment:

```bash
conda create -n match-liftover -c conda-forge -c bioconda \
  python=3.12 bcftools bcftools-liftover-plugin samtools plink plink2 pip

conda activate match-liftover
python -m pip install -e .
python -m pip install pytest build twine
```

Validate the core toolchain:

```bash
bcftools --version
samtools --version
plink --version
plink2 --version
pytest --version
prepare_variants.py --help
python -m build --version
python -m twine --version
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
apptainer build genomatch.sif docker-daemon:genomatch:latest
```

## Pre-release validation

The supported release path is the GitHub Actions workflow in `.github/workflows/release.yaml`.
Tag pushes matching `v*` publish:

- `genomatch==<version>` to PyPI
- `ghcr.io/precimed/genomatch:<version>` to GHCR
- `ghcr.io/precimed/genomatch:latest` to GHCR

Manual `workflow_dispatch` runs are used to test the release path before tagging. They can publish:

- `genomatch==<testpypi_release_label>` to TestPyPI
- `ghcr.io/precimed/genomatch:<custom-tag>` to GHCR, for example `dev`

### Release readiness

Treat a release as ready when all of the following are true:

- the repository version is correct and consistent in `pyproject.toml` and `src/genomatch/_version.py`
- `CHANGELOG.md` follows the release convention:
  ongoing work stays under `## [Unreleased]` on the development branch, and immediately before a release that heading is replaced by `## [<version>] - YYYY-MM-DD` for the release being tagged
- repository settings are in place:
  - PyPI trusted publisher for `precimed/genomatch`, workflow `release.yaml`, environment `pypi`
  - TestPyPI trusted publisher for `precimed/genomatch`, workflow `release.yaml`, environment `testpypi`
  - GitHub Actions enabled for the repository
  - GHCR publishing enabled for the repository
- a manual `workflow_dispatch` publish test succeeds for the release surfaces you care about
- the published manual workflow outputs are validated locally

An optional dry run with empty `testpypi_release_label` and empty `ghcr_tag` is a useful first check, but it is not required once the publish test path succeeds.

### Test the release workflow

Use `Actions -> Release -> Run workflow`.

Manual dispatch always runs tests, builds Python artifacts, and builds the local container image. Publishing is opt-in:

- leave `testpypi_release_label` empty to skip TestPyPI publishing
- set `testpypi_release_label` to a non-empty PEP 440 version to publish to TestPyPI
- leave `ghcr_tag` empty to skip manual GHCR publishing
- set `ghcr_tag` to a non-empty tag (for example `dev`) to publish that tag to GHCR

To test the Python package path:

1. set `testpypi_release_label` to a PEP 440 version (for example `0.1.0.dev1`)
2. leave `ghcr_tag` empty unless you also want a container publish

To test the container path:

1. set `ghcr_tag` to a custom tag (for example `dev`)
2. leave `testpypi_release_label` empty unless you also want a TestPyPI publish

Manual GHCR dispatch runs push only the custom tag you request. They do not push `latest`.

If you expect to upload multiple test package builds, use a PEP 440 development version such as `0.1.0.dev1`, `0.1.0.dev2`, and so on. TestPyPI, like PyPI, does not let you overwrite an existing version.

After a successful manual publish, validate the outputs locally.

For TestPyPI:

```bash
conda create -n genomatch-testpypi -c conda-forge -c bioconda \
  python=3.12 bcftools bcftools-liftover-plugin samtools pip

conda activate genomatch-testpypi
python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple \
  genomatch==<testpypi_release_label>

prepare_variants.py --help
```

If you published a development version such as `0.1.0.dev1`, install that exact version string.

For GHCR:

```bash
apptainer pull genomatch-dev.sif docker://ghcr.io/precimed/genomatch:dev
apptainer exec genomatch-dev.sif prepare_variants.py --help
```

## Actual release

After release readiness is established:

1. update the version in `pyproject.toml` and `src/genomatch/_version.py`
2. update `CHANGELOG.md` for the release:
   replace `## [Unreleased]` with `## [<version>] - YYYY-MM-DD` and keep the pending changes under that new release heading
3. push that release commit to GitHub
4. create and push the version tag:

```bash
git tag -a "v<version>" -m "genomatch <version>"
git push origin "v<version>"
```

5. watch the `Release` workflow in GitHub Actions publish all artifacts
6. verify the published release locally

### Verify the published release

After publication, verify both release surfaces:

```bash
python -m pip install --upgrade genomatch==<version>
prepare_variants.py --help
```

```bash
apptainer pull genomatch.sif docker://ghcr.io/precimed/genomatch:<version>
apptainer exec genomatch.sif prepare_variants.py --help
```

### Pre-release checklist

Run this before pushing the release tag:

- confirm the intended version in `pyproject.toml` and `src/genomatch/_version.py`
- confirm `CHANGELOG.md` has been converted from `## [Unreleased]` to `## [<version>] - YYYY-MM-DD` for the release commit
- run a successful manual `workflow_dispatch` publish test for TestPyPI and, if relevant, the GHCR `dev` tag
- validate the published manual workflow outputs locally
