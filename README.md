# Variant Table Toolkit

This toolkit harmonizes genetic variant data across common research formats and reference assemblies. It supports `GRCh37`, `GRCh38`, and `T2T-CHM13v2.0`, chromosomes `1-22`, `X`, `Y`, and `MT`, common contig naming modes (`ncbi`, `ucsc`, `plink`), and biallelic variants, including SNPs and supported indels.

The common workflow rewrites `chr` / `bp` / `a1` / `a2` / `snp` fields in BFILE, PFILE, VCF, or summary-statistics inputs into a standardized variant key while adjusting the attached data, such as genotypes or summary-statistic columns, accordingly. Users request the target build, contig naming, and optional filtering or normalization flags; the pipeline handles build guessing in the source data, liftover between builds (if needed), allele swaps, reference-anchored allele ordering, sorting, and duplicate removal. The workflow is split into preparation and projection phases, so users can save and reuse a prepared variant set, or project only a user-defined subset of source variants.

## Start here

1. Install the runtime using one of the supported paths in [docs/install.md](docs/install.md).
2. Download the reference FASTA/chain assets and configure `config.yaml` as described in [docs/downloads.md](docs/downloads.md).
3. Run through the worked example in [docs/tutorial.md](docs/tutorial.md).

## Documentation

- [Concepts](docs/concepts.md): `.vmap`, `.vtable`, payloads, source-row mapping, object metadata, and allele ordering.
- [Workflows](docs/workflows.md): the common prepare, combine, restrict, and project workflow.
- [Primitive tools](docs/primitives.md): lower-level import, normalization, liftover, mapping, and direct payload-application tools.

## Specifications

For exact schema and edge-case rules, see [SPEC.md](SPEC.md) and the detailed specs in [spec/](spec/). Wrapper behavior for `prepare_variants.py`, `prepare_variants_sharded.py`, and `project_payload.py` is defined in [spec/workflow.md](spec/workflow.md). Payload-application semantics are defined in [spec/payload-application.md](spec/payload-application.md).
