# SPEC.md

## Scope

This file is the technical contract for the `match/` toolkit after the provenance-first `.vmap` / `.vtable` redesign.

## Test environment

- Canonical test runs for `match/` should use the `match-liftover` conda environment.
- This applies to both focused `match-test/` runs and full-repository `pytest` runs when validating `match/` changes.
- The `match-liftover` environment is the expected place for the real `bcftools +liftover` integration tests and the migrated default test contract.

## Normative files

The normative files for the `match/` contract are:

- [README.md](README.md)
- [INSTALL.md](INSTALL.md)
- [TESTS.md](TESTS.md)
- [SPEC.md](SPEC.md)
- [core-objects.md](spec/core-objects.md)
- [contigs-and-metadata.md](spec/contigs-and-metadata.md)
- [shard-discovery.md](spec/shard-discovery.md)
- [importers.md](spec/importers.md)
- [mapping.md](spec/mapping.md)
- [variant-transforms.md](spec/variant-transforms.md)
- [ploidy-model.md](spec/ploidy-model.md)
- [payload-application.md](spec/payload-application.md)
- [sumstats-harmonization.md](spec/sumstats-harmonization.md)
- [workflow.md](spec/workflow.md)

Within [README.md](README.md), the following sections are intentionally normative and may serve as source of truth:

- the top-level object and metadata summary
- `Canonical tools`
- `Using .vmap vs .vtable as input across tools`

## Object model

- `import_*` emits `.vmap`, not `.vtable`
- importer-emitted `.vmap` is the normal starting artifact for provenance-preserving workflows
- `.vmap` provenance is `(source_shard, source_index)`
- `source_shard` is the exact discovered shard label at import time, or `.` for a single-file source
- `.vtable` remains in the system as the provenance-free materialized target-side object, produced primarily by `convert_vmap_to_target.py`, `intersect_variants.py`, and `union_variants.py`
- `.vmap` metadata remains target-side only
- importer-emitted rows preserve source-format allele meaning at import time, while reference-aware target rows emitted by `restrict_build_compatible.py` and `liftover_build.py` use ordered alleles `a1=non-reference` and `a2=reference`
- only `import_*` tools originate new provenance in v1
- for tools that accept both `.vtable` and `.vmap`, the default contract is to emit the same type as the input and preserve provenance only for `.vmap` input
- for tools that accept `.vmap` as a table-like input, `.vmap` is interpreted as its target side unless a tool-specific spec says otherwise
- `match_vmap_to_target.py` may accept a target `.vtable` or a target `.vmap`; if the target is `.vmap`, only its target side participates in matching and its provenance is ignored
- target-side transforms on `.vmap` preserve provenance unless the tool explicitly materializes a provenance-free `.vtable` or otherwise documents a provenance break

## Canonical tools surface

- `import_bim.py`
- `import_pvar.py`
- `import_vcf.py`
- `import_sumstats.py`
- `guess_build.py`
- `normalize_contigs.py`
- `restrict_contigs.py`
- `drop_strand_ambiguous.py`
- `restrict_build_compatible.py`
- `liftover_build.py`
- `sort_variants.py`
- `match_vmap_to_target.py`
- `convert_vmap_to_target.py`
- `intersect_variants.py`
- `union_variants.py`
- `apply_vmap_to_sumstats.py`
- `apply_vmap_to_bfile.py`
- `apply_vmap_to_pfile.py`

The summary-stat metadata contract is `match/schemas/raw-sumstats-metadata.yaml`.

## Workflow tools surface

- `prepare_variants.py`
- `project_payload.py`

Wrapper behavior for `prepare_variants.py` and `project_payload.py` is defined in [workflow.md](spec/workflow.md). Those wrappers are orchestration only and reuse the canonical-tool semantics defined in the topic-specific spec files.

## Reference model

The reference/config model for the `match/` toolkit is defined in [INSTALL.md](INSTALL.md).

`INSTALL.md` is the source of truth for:

- current UCSC-internal reference behavior
- required reference assets
- reference layout
- config shape
- asset resolution
- override policy
