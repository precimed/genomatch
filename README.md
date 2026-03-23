# Variant Table Toolkit

This toolkit harmonizes genetic variant data across common research formats and reference assemblies with explicit provenance and auditable row drops. It supports `GRCh37` and `GRCh38`, chromosomes `1-22`, `X`, `Y`, and `MT`, common contig naming modes (`ncbi`, `ucsc`, `plink`, `plink_splitx`), and biallelic `A/C/G/T` variants. The toolkit is designed for primary reference contigs only; alternate loci, decoy contigs, and other non-primary reference sequences are out of scope.

## Table of contents

This README is organized from concepts to workflows to reference details. It starts with the core mental model, then shows the common `prepare_variants.py` / `intersect_variants.py` / `project_payload.py` workflow with a concrete usage example. After that it covers core objects, metadata, and canonical tools in more detail, and finishes with short wrapper overviews that point to [spec/workflow.md](spec/workflow.md) for authoritative wrapper semantics.

- [Core mental model](#core-mental-model)
- [Canonical tools vs common workflow](#canonical-tools-vs-common-workflow)
- [Workflow tools](#workflow-tools)
- [Core objects and metadata](#core-objects-and-metadata)
- [Canonical tools](#canonical-tools)
- [Usage example for canonical tools](#usage-example-for-canonical-tools)

## Getting started

Before running the reference-aware tools and wrappers:

1. Ensure `bcftools` with the `+liftover` plugin is available.
2. Download the required reference data and place your own `config.yaml` next to that reference tree.
3. Set `MATCH_CONFIG=/path/to/config.yaml`.

`MATCH_CONFIG` is required. Reference paths in the config are resolved relative to the config file location. See [INSTALL.md](INSTALL.md) for the reference layout and full setup details.

## Core mental model

Conceptually, the toolkit treats raw input as three pieces that can be separated and recombined: the variant table itself (`.vtable`), the attached payload, and the provenance linking variants back to the raw source:

- `raw data = .vtable + payload + provenance`
- `.vtable` is the ordered variant table, one row per biallelic variant
- `payload` is the external format-specific data attached to those variants, such as summary-stat columns, PLINK 1 BED/BIM/FAM genotype matrices, or PLINK 2 PGEN/PVAR/PSAM genotype matrices; payloads can be split across multiple files (shards)
- `provenance` is the link back from variant to its raw payload, represented as `(source_shard, source_index)` plus `allele_op` to capture `a1 <-> a2` allele swaps and strand flips
- `.vmap = .vtable + provenance`; `.vmap` is therefore the object emitted by import operations, and used for most downstream transforms: it keeps transformed table of variants together with enough provenance to rewrite the original payload later, but it does not contain the payload itself.

The toolkit provides a set of tools for transforming tables of variants: importing raw inputs, normalizing representation, restricting rows, liftover between builds, performing intersections and unions of variant tables. These transformations can carry provenance linking each transformed variant row back to the original payload. That is what `.vmap` stores: a transformed variant table together with provenance, but without the payload itself. We refer to that transformed variant table as the *target side* of `.vmap`, while its *source side* is the original payload that the stored provenance links back to. Target side of `.vmap` can be seen as a `.vtable` object that only stores the transformed variant table, without provenance. `apply_vmap_*` reconnects the transformed variant rows in `.vmap` to the original payload and writes a rewritten payload in that transformed row order.

- In the following diagram `.vmap` is shown on the right, and its provenance points to source payload on the left:
  ```text
                         +---- .vmap -----+
       payload           |   .vtable +    |
       (source) <--------|---provenance   |
                         +----------------+
  ```
- Transformations will update `.vmap` while still keeping its provenance pointing to the same source payload

Intersection and matching of variants is performed on `chr:bp:a1:a2`, ignoring variant `id`. A `prepare_variants.py` pipeline combines operations that standardize `chr:bp:a1:a2`, making it into a unique variant identifier.
Intersection and matching support biallelic SNPs and non-SNPs, while reference-aware restriction and liftover are narrower and are specified in the spec docs.

For exact schema and edge-case rules, see [SPEC.md](SPEC.md), [core-objects.md](spec/core-objects.md), [contigs-and-metadata.md](spec/contigs-and-metadata.md), [shard-discovery.md](spec/shard-discovery.md), [importers.md](spec/importers.md), [variant-transforms.md](spec/variant-transforms.md), [mapping.md](spec/mapping.md), [payload-application.md](spec/payload-application.md), and [workflow.md](spec/workflow.md). Reference-aware setup and assets are described in [INSTALL.md](INSTALL.md).

## Canonical tools vs workflow tools

You can use this toolkit at two levels:

- use the canonical tools when you need explicit control over import, normalization, liftover, matching, intersection, and payload application
- use `prepare_variants.py`, `intersect_variants.py`, `union_variants.py`, and `project_payload.py` workflow when you want to
  - prepare several inputs (e.g. cohort, reference, and summary-statistics payloads) into the same build and contig naming
  - intersect a subset of inputs to get one shared variant table, or union prepared inputs when that is the desired target set
  - project all original payloads into that shared target table

For individual-level data, `prepare_variants.py` and `project_payload.py` can be executed
within trusted research environments hosting sensitive data.
`prepare_variants.py` produces only lists of variants which can be copied to a joint location for
`intersect_variants.py`.

## Workflow tools

Use the tables below as a quick wrapper reference. The authoritative canonical-tool and wrapper contracts remain in the spec.

`prepare_variants.py` is the convenience wrapper for preparing one raw input into a final `<output>.vmap` that can later be intersected or projected back onto the original payload. It orchestrates import, contig normalization, metadata resolution, same-build restriction or liftover, and optional strand-ambiguous / contig filtering, while retaining the intermediate `.vmap` stages. For the exact stage graph and `--resume` / `--force` behavior, see [spec/workflow.md](spec/workflow.md).

| `prepare_variants.py` argument | Meaning |
| --- | --- |
| `--input-format` | Required importer format: `bim`, `pvar`, `vcf`, or `sumstats`. See [spec/importers.md](spec/importers.md) |
| `--input`, `--output` | Required raw input payload and output stem. Genotype payload inputs may be sharded via `@`; see [spec/shard-discovery.md](spec/shard-discovery.md). The final prepared output is always a single `<output>.vmap` |
| `--sumstats-metadata`, `--id-vtable` | Sumstats-only controls: `--sumstats-metadata` is required for `--input-format sumstats` and follows the raw metadata schema at [schemas/raw-sumstats-metadata.yaml](schemas/raw-sumstats-metadata.yaml), which is the same specification used by [BioPsyk/cleansumstats](https://github.com/BioPsyk/cleansumstats). `--id-vtable` optionally fills missing `chr`/`pos` by matching on variant ID, and becomes required when those fields are absent from the summary-stat metadata |
| `--dst-contig-naming` | Target contig naming; defaults to `ncbi`. Supported values are `ncbi`, `ucsc`, `plink`, and `plink_splitx`; `plink_splitx` follows PLINK `--split-x` style X/XY_PAR labeling. See [spec/contigs-and-metadata.md](spec/contigs-and-metadata.md) |
| `--dst-build` | Target genome build; defaults to `GRCh38`. If the input build differs, liftover rewrites both coordinates and alleles via `bcftools +liftover`, then emits canonical target rows with `a2=reference`. See [spec/variant-transforms.md](spec/variant-transforms.md) |
| `--[no-]allow-strand-flips` | Control strand-flip allowance during same-build reference-aware restriction; enabled by default |
| `--[no-]norm-indels` | Control indel normalization during same-build reference-aware restriction; `bcftools norm` is used internally when enabled |
| `--drop-strand-ambiguous` | Optional removal of strand-ambiguous variants after build compatibility / liftover |
| `--chr2use` / `--contigs` | Optional restriction of final prepared output to selected contigs. `--contigs` is an alias for `--chr2use`. For non-autosomes and ploidy-related behavior, see [spec/ploidy-model.md](spec/ploidy-model.md) |
| `--prefix` | Optional retained-intermediate stem; defaults to `--output` |
| `--resume`, `--force` | Wrapper execution controls; mutually exclusive |

`project_payload.py` is the convenience wrapper for projecting an original payload into an already prepared target `.vtable` or `.vmap`. `--target` always defines the target row set. If `--source-vmap` is supplied, the wrapper first matches source provenance onto that target and retains `<prefix>.vmap`; if `--source-vmap` is omitted, `--target` itself must be a `.vmap` and is applied directly. For genotype payloads, the wrapper can also restrict, reorder, or reconcile subject axes via explicit target sample files or wrapper-level union synthesis. For the exact projection flow, output retention rules, and target-sample reconciliation semantics, see [spec/workflow.md](spec/workflow.md).

| `project_payload.py` argument | Meaning |
| --- | --- |
| `--input-format` | Required payload format: `bfile`, `pfile`, `sumstats`, or `sumstats-clean` (`bfile` maps to `bim` from `prepare_variants.py`; `pfile` maps to prepared `pvar`). See [spec/payload-application.md](spec/payload-application.md) |
| `--input`, `--output` | Required raw payload to project and rewritten payload destination. Genotype payload `--input` may be sharded via `@`; for `bfile` / `pfile`, `--output` is a PLINK output prefix and may also be sharded via `@` |
| `--target` | Required target `.vtable` or `.vmap`; always defines the target row set. See [spec/mapping.md](spec/mapping.md) |
| `--source-vmap` | Optional prepared source `.vmap`; required when `--target` is a `.vtable`. If provided, it is the source of provenance even when `--target` is a `.vmap` |
| `--full-target` | Optional wrapper opt-out from the default mapped-only projection. By default, `project_payload.py` keeps only mapped target rows; `--full-target` retains unmatched target rows in the projected payload |
| `--sumstats-metadata`, `--fill-mode`, `--use-af-inference` | Summary-stat controls: `--sumstats-metadata` is required for `sumstats` and `sumstats-clean`, and `--fill-mode` / `--use-af-inference` apply to `sumstats-clean`. See [spec/sumstats-harmonization.md](spec/sumstats-harmonization.md) |
| `--target-fam` / `--target-psam` / `--sample-axis union` / `--sample-id-mode {fid_iid,iid}` | Optional genotype payload subject-axis controls for explicit target sample files, wrapper-synthesized union targets, and subject keying. See [Sample-axis reconciliation for genotype payloads](spec/payload-application.md#sample-axis-reconciliation-for-genotype-payloads) |
| `--prefix` | Optional retained matched-mapping prefix. It defaults to `--output`, except that sharded `@` PLINK outputs derive an `all_targets` prefix so the retained single-file `<prefix>.vmap` has a concrete non-sharded name |
| `--force` | Delete wrapper-managed outputs first, then rerun cleanly |

### Common workflow example

#### Step 1: prepare each input once

Prepare a cohort `.bim` payload:

```bash
python3 match/prepare_variants.py \
  --input cohort.bim \
  --input-format bim \
  --output work/cohort.prepared
```

Prepare a sharded reference `.bim` payload, split across multiple inputs:

```bash
python3 match/prepare_variants.py \
  --input reference.@.bim \
  --input-format bim \
  --output work/reference.prepared
```

Prepare a summary-statistics payload:

```bash
python3 match/prepare_variants.py \
  --input study.tsv.gz \
  --input-format sumstats \
  --sumstats-metadata study.yaml \
  --output work/study.prepared
```

Here `--output` is a stem, so these commands retain stage-specific `.vmap` outputs under prefixes that default to the chosen `--output`, and also copy the final prepared objects to:

- `work/cohort.prepared.vmap`
- `work/reference.prepared.vmap`
- `work/study.prepared.vmap`

#### Step 2: intersect the prepared outputs

Once the inputs are prepared, intersect the prepared cohort and reference `.vmap` outputs:

```bash
python3 match/intersect_variants.py \
  --inputs work/cohort.prepared.vmap work/reference.prepared.vmap \
  --output work/shared.vtable
```

The output `work/shared.vtable` is the shared target set.
If you want a merged target set instead of an exact intersection, use `union_variants.py` on the prepared inputs instead.

#### Step 3: project payloads into the shared target set

Project the summary-statistics payload:

```bash
python3 match/project_payload.py \
  --input study.tsv.gz \
  --input-format sumstats-clean \
  --sumstats-metadata study.yaml \
  --source-vmap work/study.prepared.vmap \
  --target work/shared.vtable \
  --output study.shared.tsv.gz
```
As summary-statistics `.vmap` was not part of the intersection
it's likely that there are variants in `--target` that are missing from summary statistics.
`project_payload.py` now does that mapped-only projection by default. Use `--full-target`
if you want to retain unmatched target rows as explicit missing rows. With
`--input-format sumstats-clean`, the rewritten output is harmonized TSV.
This clean harmonization step is a missing-value and missing-column completion pipeline. It may infer missing fields from other available fields, but it does not try to verify semantic consistency between overlapping values that are already present.

Project the PLINK payload:

```bash
python3 match/project_payload.py \
  --input cohort.bim \
  --input-format bfile \
  --source-vmap work/cohort.prepared.vmap \
  --target work/shared.vtable \
  --output cohort.shared
```

Project the sharded reference payload:

```bash
python3 match/project_payload.py \
  --input reference.@.bim \
  --input-format bfile \
  --source-vmap work/reference.prepared.vmap \
  --target work/shared.vtable \
  --output reference.shared.@
```

This writes the sharded payload outputs under `reference.shared.<contig>` and retains the matched mapping as `reference.shared.all_targets.vmap` unless `--prefix` is set explicitly.

## Core objects and metadata

Each `.vtable` or `.vmap` object is a single file with JSON sidecars at `<object>.meta.json`.

| Object | Purpose | Required metadata |
| --- | --- | --- |
| `.vmap` | Provenance-carrying object for import, cleanup, reference-compatible restriction, liftover, matching, and application | `object_type=variant_map`, `target.genome_build`, `target.contig_naming` |
| `.vtable` | Provenance-free table, produced mainly by `convert_vmap_to_target.py`, `intersect_variants.py`, and `union_variants.py` | `object_type=variant_table`, `genome_build`, `contig_naming` |

Downstream transforms other than `normalize_contigs.py` require declared `contig_naming`. If an imported object is unresolved or mixed, clean it with `normalize_contigs.py` before continuing. See [core-objects.md](spec/core-objects.md) for the on-disk column contracts.

`a1` and `a2` are ordered alleles, not an unordered allele set. Importers preserve source-format allele meaning at import time, while reference-aware tools such as `restrict_build_compatible.py` and `liftover_build.py` emit reference-aligned variants with `a1=non-reference` and `a2=reference`. This keeps canonical reference-aware rows aligned with ALT/effect/counting-style semantics while still making the reference allele explicit in `a2`.

## Canonical tools

| Task | Tool | Typical input | Emits | Notes |
| --- | --- | --- | --- | --- |
| Import | `import_bim.py` | raw `.bim` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_pvar.py` | raw `.pvar` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_vcf.py` | raw `.vcf` / `.vcf.gz` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_sumstats.py` | raw summary statistics plus metadata YAML | `.vmap` | Creates provenance; supports `--chr2use` / `--contigs` and optional `--id-vtable`; single-file only |
| Normalize / metadata | `guess_build.py` | `.vmap` or `.vtable` | same as input | Metadata-only update; on `.vmap`, updates target-side metadata only |
| Normalize / metadata | `normalize_contigs.py` | `.vmap` or `.vtable` | same as input | Repairs or standardizes target contigs; target naming is selected with `--to` |
| Normalize / metadata | `restrict_contigs.py` | `.vmap` or `.vtable` | same as input | Target-side contig filter via `--chr2use` / `--contigs` |
| Normalize / metadata | `drop_strand_ambiguous.py` | `.vmap` or `.vtable` | same as input | Drops strand-ambiguous target rows |
| Validate / liftover | `restrict_build_compatible.py` | `.vmap` or `.vtable` | same as input | Same-build, reference-aware filtering; key flags are `--allow-strand-flips`, `--norm-indels`, and `--sort` |
| Validate / liftover | `liftover_build.py` | `.vmap` or `.vtable` | same as input | Explicit build conversion to `--target-build`; preserves `.vmap` provenance and re-sorts into declared coordinate order |
| Order  | `sort_variants.py` | `.vmap` or `.vtable` | same as input | Explicit standalone target-row sorting by declared contig order, then numeric position |
| Match / intersect / materialize | `match_vmap_to_target.py` | source `.vmap` plus target `.vtable` or `.vmap` | `.vmap` | Exact same-build mapping onto target rows; target `.vmap` provenance is ignored with a warning |
| Match / intersect / materialize | `intersect_variants.py` | two or more `.vmap` / `.vtable` inputs | `.vtable` | Exact intersection; output order follows the first input |
| Match / intersect / materialize | `union_variants.py` | two or more `.vmap` / `.vtable` inputs | `.vtable` | Exact union with first-occurrence deduplication; output is re-sorted into declared coordinate order |
| Match / intersect / materialize | `convert_vmap_to_target.py` | `.vmap` | `.vtable` | Materializes the `.vmap` target side and intentionally drops provenance |
| Apply to payloads | `apply_vmap_to_sumstats.py` | `.vmap` plus summary-stat payload | rewritten payload | Single-file payload only; keeps full target-row order by default, supports `--only-mapped-target`, and exposes `--clean` for canonical cleaned output |
| Apply to payloads | `apply_vmap_to_bfile.py` | `.vmap` plus PLINK 1 BED/BIM/FAM payload | rewritten payload | Shared genotype-payload `apply_vmap_*` contract; supports `--only-mapped-target`, payload `@` sharding, bounded-chunk processing, explicit `--target-fam`, and `--sample-id-mode {fid_iid,iid}` |
| Apply to payloads | `apply_vmap_to_pfile.py` | `.vmap` plus PLINK 2 PGEN/PVAR/PSAM payload | rewritten payload | Shared genotype-payload `apply_vmap_*` contract; supports `--only-mapped-target`, payload `@` sharding, bounded-chunk processing, explicit `--target-psam`, `--sample-id-mode {fid_iid,iid}`, preserves supported biallelic hardcalls/phase/unphased dosages, and emits `.pgen/.pvar/.psam` |

### Using .vmap vs .vtable as input across tools

- `import_*` tools are the only tools that originate new provenance, so they emit `.vmap`.
- For tools that accept both `.vmap` and `.vtable`, the default rule is to act on target rows and emit the same artifact type as the input.
- `match_vmap_to_target.py` is stricter: source must be `.vmap`, target may be `.vtable` or `.vmap`, and the output is `.vmap`; if the target is `.vmap`, matching uses only its target side, ignores its provenance entirely, and emits a warning. Use `convert_vmap_to_target.py` first if you want a provenance-free target and no warning.
- `apply_vmap_*` consumes `.vmap` as a mapping object, not as a table. The payload must still match the `.vmap` source provenance exactly.

### Important usage notes

- Import tools always preserve source provenance for retained rows, also when it drops unrepresentable rows or rows filtered by `--chr2use` / `--contigs`
- import tools do not repair mixed or invalid contig labels. Use `normalize_contigs.py` for cleanup.
- `restrict_build_compatible.py` is the same-build, reference-aware canonicalization step: it filters to reference-compatible rows, can optionally re-sort into declared coordinate order, and can optionally normalize supported indels via retained `bcftools norm` intermediates while still emitting canonical biallelic `a1=non-reference`, `a2=reference` rows; it does not liftover.
- `liftover_build.py` supports SNVs and one-side-multibase reference-anchored rows, but not lifted outputs where both final alleles remain multi-base.
- `intersect_variants.py` preserves first-input order. `union_variants.py` first deduplicates by first occurrence across the full input stream and then uses declared coordinate order. `match_vmap_to_target.py` follows target-side order. `liftover_build.py` and `sort_variants.py` use declared coordinate order: declared contig order, then numeric position, with stable input order for ties.
- coordinate-changing transforms and genotype-payload application follow the shared expected-ploidy contract in [spec/ploidy-model.md](spec/ploidy-model.md).
- genotype-payload `apply_vmap_*` tools apply `.vmap` provenance back to the original payload, support `@` source discovery and `@` target-contig output sharding, and follow the shared ploidy contract.
- both genotype-payload `apply_vmap_*` tools support explicit target sample files plus `--sample-id-mode {fid_iid,iid}` for subject restriction, reordering, and reconciliation. Explicit target sample files define output subject order exactly and drive reconciliation-added missingness reporting.
- `apply_vmap_to_bfile.py` uses `--target-fam`; `apply_vmap_to_pfile.py` uses `--target-psam`. Without explicit target sample files, sharded source payloads must still have identical sample metadata across referenced shards.
- `project_payload.py --sample-axis union` is a wrapper-only convenience for sharded genotype payloads. It unions subject keys across referenced retained mapped shards only, synthesizes and retains `<prefix>.target_samples.fam` or `<prefix>.target_samples.psam`, then calls the canonical apply tool with that explicit target sample file. For non-sharded input, or when only one referenced shard remains, it warns and is a no-op.
- `@` sharding is supported for heavy raw inputs and genotype payloads; summary-stat payloads, `.vmap`, and `.vtable` objects are single-file.

### Usage example for canonical tools

Use case: liftover a summary-statistics file to `GRCh38` before downstream matching or application.

```bash
python3 match/import_sumstats.py \
  --input study.tsv.gz \
  --sumstats-metadata study_meta.yaml \
  --output study.raw.vmap

python3 match/normalize_contigs.py --to ncbi \
  --input study.raw.vmap \
  --output study.clean.vmap

python3 match/guess_build.py \
  --input study.clean.vmap \
  --write

python3 match/restrict_build_compatible.py \
  --source study.clean.vmap \
  --allow-strand-flips \
  --norm-indels \
  --output study.validated.vmap

python3 match/liftover_build.py --target-build GRCh38 \
  --input study.validated.vmap \
  --output study.grch38.vmap

python3 match/apply_vmap_to_sumstats.py \
  --input study.tsv.gz \
  --sumstats-metadata study_meta.yaml \
  --vmap study.grch38.vmap \
  --output study.grch38.tsv.gz
```

Use case: build one shared `GRCh38` variant universe for chromosome-sharded cohort genotypes, a target reference panel such as 1000 Genomes, and the lifted summary statistics from the first workflow, then emit a constrained cohort PLINK dataset plus aligned summary statistics in that same universe.

```bash
python3 match/import_vcf.py \
  --input kgp.vcf.gz \
  --output kgp.raw.vmap \
  --genome-build GRCh37

python3 match/normalize_contigs.py --to ncbi \
  --input kgp.raw.vmap \
  --output kgp.clean.vmap

python3 match/restrict_build_compatible.py \
  --source kgp.clean.vmap \
  --allow-strand-flips \
  --norm-indels \
  --output kgp.validated.vmap

python3 match/liftover_build.py --target-build GRCh38 \
  --input kgp.validated.vmap \
  --output kgp.grch38.vmap

python3 match/convert_vmap_to_target.py \
  --source kgp.grch38.vmap \
  --output kgp.grch38.vtable

python3 match/import_bim.py \
  --input cohort.@.bim \
  --output cohort.raw.vmap \
  --genome-build GRCh38

python3 match/normalize_contigs.py --to ncbi \
  --input cohort.raw.vmap \
  --output cohort.clean.vmap

python3 match/intersect_variants.py \
  --inputs cohort.clean.vmap kgp.grch38.vtable \
  --output shared.grch38.vtable

python3 match/match_vmap_to_target.py \
  --source cohort.clean.vmap \
  --target shared.grch38.vtable \
  --output cohort.shared.vmap

python3 match/apply_vmap_to_bfile.py \
  --source-prefix cohort.@ \
  --vmap cohort.shared.vmap \
  --output-prefix cohort.shared.@

python3 match/match_vmap_to_target.py \
  --source study.grch38.vmap \
  --target shared.grch38.vtable \
  --output study.shared.vmap

python3 match/apply_vmap_to_sumstats.py \
  --input study.tsv.gz \
  --sumstats-metadata study_meta.yaml \
  --vmap study.shared.vmap \
  --output study.shared.tsv.gz
```

For a PLINK 2 payload, replace the final genotype-application step with `apply_vmap_to_pfile.py`; it follows the same shared row-selection, provenance, and sharding contract, but writes `.pgen/.pvar/.psam` and preserves supported biallelic hardcalls, hardcall phase, and unphased dosages.
