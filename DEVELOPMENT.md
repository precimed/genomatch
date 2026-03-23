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
  python=3.12 bcftools bcftools-liftover-plugin samtools numpy pandas pysam pytest pyyaml pgenlib python-build twine

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
apptainer build genomatch.sif docker-daemon://genomatch:latest
```

## Maintainer notes

This file is the place to document maintainer workflows for:

- building and publishing container artifacts
- preparing and publishing PyPI releases

Keep end-user install instructions in [INSTALL.md](INSTALL.md). Keep reference setup in [DOWNLOADS.md](DOWNLOADS.md).

## Releasing a version

The primary release path is the tag-driven GitHub Actions workflow in `.github/workflows/release.yaml`.

The current release targets are:

- PyPI package: `genomatch==<version>`
- GHCR container: `ghcr.io/precimed/genomatch:<version>`
- GHCR moving tag: `ghcr.io/precimed/genomatch:latest`

### Release readiness for `0.1.0`

`0.1.0` is ready to release from the repository side:

- package metadata and console scripts are present in `pyproject.toml`
- local package build succeeds
- `twine check` passes on the built wheel and sdist
- the Docker image builds successfully
- the test suite passes in `match-liftover`

What still depends on maintainer credentials outside the repo:

- PyPI upload permission for the `genomatch` project
- GHCR push permission for `ghcr.io/precimed/genomatch`
- GitHub repository settings configured for the release workflow

### Primary release path: GitHub Actions workflow

The workflow filename is `release.yaml`. It triggers on pushed tags matching `v*`.
It also supports manual `workflow_dispatch` runs.

Before the workflow can publish `0.1.0`, configure:

- a PyPI trusted publisher for the `genomatch` project pointing at `precimed/genomatch`
- workflow name set to `release.yaml`
- environment name set to `pypi`
- GitHub Actions enabled for this repository
- GHCR package publishing enabled for this repository

The workflow publishes:

1. `genomatch==<version>` to PyPI
2. `ghcr.io/precimed/genomatch:<version>` to GHCR
3. `ghcr.io/precimed/genomatch:latest` to GHCR

The workflow will:

1. verify that the tag version matches `pyproject.toml` and `src/genomatch/_version.py`
2. run the full test suite in a `match-liftover`-like micromamba environment
3. build and validate the PyPI wheel and sdist
4. publish the package to PyPI
5. build and push `ghcr.io/precimed/genomatch:<version>` and `ghcr.io/precimed/genomatch:latest`

### Manual workflow dry run

You can test the workflow without publishing anything from the GitHub Actions UI:

1. open the `Release` workflow
2. choose `Run workflow`
3. leave `publish_pypi` unchecked
4. leave `publish_ghcr` unchecked
5. optionally leave `release_version` empty to validate the version already declared in the repository

In this mode the workflow will still:

- verify the package version metadata
- run the full test suite
- build and validate the wheel and sdist
- build the Docker image locally in CI

but it will skip:

- PyPI publication
- GHCR pushes

### Preferred release flow

1. Make sure the release commit is already on GitHub and includes `.github/workflows/release.yaml`.
2. Create and push the version tag:

```bash
git tag -a "v<version>" -m "genomatch <version>"
git push origin "v<version>"
```

3. Watch the `Release` workflow in GitHub Actions.
4. After success, verify the published PyPI package and GHCR image.

### Local fallback only

The remaining sections describe the equivalent local fallback workflow for maintainers. Use it only when you need a manual release path outside GitHub Actions.

### 1. Prepare the release commit

If you are preparing a new release version, update both version declarations in the same commit:

```bash
$EDITOR pyproject.toml
$EDITOR src/genomatch/_version.py
git commit -am "Release genomatch <version>"
git tag -a "v<version>" -m "genomatch <version>"
```

For `0.1.0`, both files should already say `0.1.0`.

### 2. Build the PyPI artifacts

Run the build from the repository root with `python -P` so the current working directory is not prepended to `sys.path`. This avoids local artifact directories such as `build/` shadowing the PyPA `build` package.

Also use `--no-isolation`, because the default isolated build mode tries to download build dependencies from the network.

```bash
rm -rf build dist src/genomatch.egg-info
python -P -m build --no-isolation .
```

Validate the artifacts:

```bash
python -m twine check dist/*
```

### 3. Upload the PyPI release

Upload only after the test suite and artifact checks pass:

```bash
python -m twine upload dist/*
```

If you want to test the publication path first, upload to TestPyPI instead:

```bash
python -m twine upload --repository testpypi dist/*
```

### 4. Build and publish the GHCR container

Build and tag the OCI image with both the release tag and `latest`:

```bash
docker build \
  -t ghcr.io/precimed/genomatch:<version> \
  -t ghcr.io/precimed/genomatch:latest \
  .
```

Log in to GHCR and push both tags:

```bash
echo "$GITHUB_TOKEN" | docker login ghcr.io -u "$GITHUB_USER" --password-stdin
docker push ghcr.io/precimed/genomatch:<version>
docker push ghcr.io/precimed/genomatch:latest
```

For `0.1.0`, that becomes:

```bash
docker build \
  -t ghcr.io/precimed/genomatch:0.1.0 \
  -t ghcr.io/precimed/genomatch:latest \
  .

docker push ghcr.io/precimed/genomatch:0.1.0
docker push ghcr.io/precimed/genomatch:latest
```

The published image is the artifact Apptainer users will consume via:

```bash
apptainer pull genomatch.sif oras://ghcr.io/precimed/genomatch:0.1.0
```

or:

```bash
apptainer pull genomatch.sif oras://ghcr.io/precimed/genomatch:latest
```

### 5. Verify the published release

After publication, verify both release surfaces:

```bash
python -m pip install --upgrade genomatch==<version>
prepare_variants.py --help
```

```bash
apptainer pull genomatch.sif oras://ghcr.io/precimed/genomatch:<version>
apptainer exec genomatch.sif prepare_variants.py --help
```

### Pre-release checklist

Run this before uploading or pushing release tags:

```bash
conda activate match-liftover
pytest tests
docker build -t genomatch:test .
rm -rf build dist src/genomatch.egg-info
python -P -m build --no-isolation .
python -m twine check dist/*
```
