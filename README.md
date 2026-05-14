# Variant Table Toolkit

This toolkit harmonizes genetic variant data across common research formats and reference assemblies with explicit provenance and auditable row drops. It supports `GRCh37`, `GRCh38`, and `T2T-CHM13v2.0`, chromosomes `1-22`, `X`, `Y`, and `MT`, common contig naming modes (`ncbi`, `ucsc`, `plink`, `plink_splitx`), and biallelic `A/C/G/T` variants. The toolkit is designed for primary reference contigs only; alternate loci, decoy contigs, and other non-primary reference sequences are out of scope.

## Table of contents

This README starts with end-user setup pointers, then introduces the core mental model, the common `prepare_variants.py` / `prepare_variants_sharded.py` / `intersect_variants.py` / `restrict_vmap.py` / `project_payload.py` workflow, and the canonical tool surface. It finishes with the detailed object and tool summaries that point to [spec/workflow.md](spec/workflow.md) and the rest of the spec for authoritative semantics.

- [Getting started](#getting-started)
- [Core mental model](#core-mental-model)
- [Canonical tools vs common workflow](#canonical-tools-vs-common-workflow)
- [Workflow tools](#workflow-tools)
- [Core objects and metadata](#core-objects-and-metadata)
- [Canonical tools](#canonical-tools)
- [Usage example for canonical tools](#usage-example-for-canonical-tools)

## Getting started

For end-user setup:

1. Follow [INSTALL.md](INSTALL.md) to install the runtime, either with `pip install genomatch` on top of local `bcftools` plus `+liftover`, or via the published Apptainer image.
2. Follow [DOWNLOADS.md](DOWNLOADS.md) to download the required reference data, place `config.yaml`, and set `MATCH_CONFIG`.

`MATCH_CONFIG` is required for reference-aware commands. Reference paths in the config are resolved relative to the config file location. The published container image packages the runtime only; users must still supply reference assets and `config.yaml` separately.

## Core mental model

Conceptually, the toolkit treats raw input as three pieces that can be separated and recombined: the variant table itself (`.vtable`), the attached payload, and the provenance linking variants back to the raw source:

- `raw data = .vtable + payload + provenance`
- `.vtable` is the ordered variant table, one row per biallelic variant
- `payload` is the external format-specific data attached to those variants, such as summary-stat columns, PLINK 1 BED/BIM/FAM genotype matrices, or PLINK 2 PGEN/PVAR/PSAM genotype matrices; payloads can be split across multiple files (shards)
- `provenance` is the link back from variant to its raw payload, represented as `(source_shard, source_index)` plus `allele_op` to capture `a1 <-> a2` allele swaps and strand flips
- `.vmap = .vtable + provenance`; `.vmap` is therefore the object emitted by import operations, and used for most downstream transforms: it keeps transformed table of mapped variants together with enough provenance to rewrite the original payload later, but it does not contain the payload itself. A canonical `.vmap` is mapped-only; rows missing from a source payload are absent rather than represented by sentinel rows.

The toolkit provides a set of tools for transforming tables of variants: importing raw inputs, normalizing representation, restricting rows, liftover between builds, performing intersections and unions of variant tables. These transformations can carry provenance linking each transformed variant row back to the original payload. That is what `.vmap` stores: a transformed variant table together with provenance, but without the payload itself. We refer to that transformed variant table as the *target side* of `.vmap`, while its *source side* is the original payload that the stored provenance links back to. Target side of `.vmap` can be seen as a `.vtable` object that only stores the transformed variant table, without provenance. `apply_vmap_*` reconnects the transformed variant rows in `.vmap` to the original payload and writes a rewritten payload in that transformed row order.

- In the following diagram `.vmap` is shown on the right, and its provenance points to source payload on the left:
  ```text
                         +---- .vmap -----+
       payload           |   .vtable +    |
       (source) <--------|---provenance   |
                         +----------------+
  ```
- Transformations will update `.vmap` while still keeping its provenance pointing to the same source payload

Intersection of variants is performed on exact `chr:bp:a1:a2`, ignoring variant `id`. A `prepare_variants.py` pipeline combines operations that standardize `chr:bp:a1:a2`, making it into a unique variant identifier.
Intersection supports biallelic SNPs and non-SNPs, while reference-aware restriction and liftover are narrower and are specified in the spec docs.

For exact schema and edge-case rules, see [SPEC.md](SPEC.md), [core-objects.md](spec/core-objects.md), [contigs-and-metadata.md](spec/contigs-and-metadata.md), [shard-discovery.md](spec/shard-discovery.md), [importers.md](spec/importers.md), [variant-transforms.md](spec/variant-transforms.md), [mapping.md](spec/mapping.md), [payload-application.md](spec/payload-application.md), and [workflow.md](spec/workflow.md). End-user installation is described in [INSTALL.md](INSTALL.md), and reference-aware assets/config in [DOWNLOADS.md](DOWNLOADS.md).

## Canonical tools vs workflow tools

You can use this toolkit at two levels:

- use the canonical tools when you need explicit control over import, normalization, liftover, intersection, and payload application
- use `prepare_variants.py` / `prepare_variants_sharded.py`, `intersect_variants.py`, `union_variants.py`, and `project_payload.py` workflow when you want to
  - prepare several inputs (e.g. cohort, reference, and summary-statistics payloads) into the same build and contig naming
  - intersect a subset of inputs to get a shared membership `.vtable`, or union prepared inputs when that is the desired membership set
  - apply original payloads through explicit payload-specific `.vmap` restrictions

For individual-level data, `prepare_variants.py`, `prepare_variants_sharded.py`, and `project_payload.py` can be executed
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
| `--dst-build` | Target genome build (`GRCh37`, `GRCh38`, or `T2T-CHM13v2.0`); defaults to `GRCh38`. If the input build differs, liftover rewrites both coordinates and alleles via `bcftools +liftover`, then emits canonical target rows with `a2=reference`. See [spec/variant-transforms.md](spec/variant-transforms.md) |
| `--[no-]allow-strand-flips` | Control strand-flip allowance during same-build reference-aware restriction; enabled by default |
| `--[no-]norm-indels` | Control indel normalization during same-build reference-aware restriction; `bcftools norm` is used internally when enabled |
| `--drop-strand-ambiguous` | Optional removal of strand-ambiguous variants after build compatibility / liftover |
| `--chr2use` / `--contigs` | Optional restriction of final prepared output to selected contigs. `--contigs` is an alias for `--chr2use`. For non-autosomes and ploidy-related behavior, see [spec/ploidy-model.md](spec/ploidy-model.md) |
| `--prefix` | Optional retained-intermediate stem; defaults to `--output` |
| `--resume`, `--force` | Wrapper execution controls; mutually exclusive |

`prepare_variants_sharded.py` is the memory-bounded sharded-input variant of `prepare_variants.py`. It accepts `bim`, `pvar`, and `vcf` inputs whose `--input` contains `@`, groups selected input shards by unambiguous target contig, runs the full `prepare_variants.py` stage graph once per target-contig group, and concatenates the per-contig-group final `.vmap` files into one final `<output>.vmap`. Its stage semantics are inherited from `prepare_variants.py`; see [spec/workflow.md](spec/workflow.md).

| `prepare_variants_sharded.py` argument | Meaning |
| --- | --- |
| `--input`, `--input-format` | Required sharded raw input. `--input` must contain `@`; `--input-format` may be `bim`, `pvar`, or `vcf` |
| `--output` | Required non-sharded output stem. The final prepared output is `<output>.vmap` |
| `--prefix` | Required retained-intermediate stem and must contain `@`; per-contig-group retained outputs replace `@` with the target contig group token |
| Other preparation flags | Same meaning as `prepare_variants.py` and passed through to each per-contig-group invocation |
| `--resume`, `--force` | Wrapper execution controls; mutually exclusive |

`project_payload.py` is the convenience wrapper for applying an original payload through an explicit prepared `.vmap`. It requires `--vmap`, does not accept a target `.vtable`, and does not run `restrict_vmap.py`; users who want to apply a payload to a shared membership table must first run `restrict_vmap.py` themselves. Payload row order is always `.vmap` row order. For genotype payloads, the wrapper can also reconcile subject axes via explicit target sample files or wrapper-level union synthesis. For the exact payload-application flow, output rules, and target-sample reconciliation semantics, see [spec/workflow.md](spec/workflow.md).

| `project_payload.py` argument | Meaning |
| --- | --- |
| `--input-format` | Required payload format: `bfile`, `pfile`, `sumstats`, or `sumstats-clean` (`bfile` maps to `bim` from `prepare_variants.py`; `pfile` maps to prepared `pvar`). See [spec/payload-application.md](spec/payload-application.md) |
| `--input`, `--output` | Required raw payload to project and rewritten payload destination. Genotype payload `--input` may be sharded via `@`; for `bfile` / `pfile`, `--output` is a PLINK output prefix and may also be sharded via `@` |
| `--vmap` | Required mapped-only `.vmap` to apply; it defines output variant rows, row order, IDs, source provenance, and `allele_op` |
| `--retain-snp-id` | Optional payload projection control. By default, `project_payload.py` writes corrected output IDs as `chrom:pos:a1:a2` from `.vmap` target rows; `--retain-snp-id` writes target-side `.vmap` IDs instead |
| `--sumstats-metadata`, `--fill-mode`, `--use-af-inference` | Summary-stat controls: `--sumstats-metadata` is required for `sumstats` and `sumstats-clean`, and `--fill-mode` / `--use-af-inference` apply to `sumstats-clean`. See [spec/sumstats-harmonization.md](spec/sumstats-harmonization.md) |
| `--target-fam` / `--target-psam` / `--sample-axis {union,native}` / `--sample-id-mode {fid_iid,iid}` / `--skip-ploidy-check` | Optional genotype payload subject-axis controls for explicit target sample files, wrapper-synthesized union targets, native per-output-shard sample axes, subject keying, and ploidy-validation skipping. See [Sample-axis reconciliation for genotype payloads](spec/payload-application.md#sample-axis-reconciliation-for-genotype-payloads) |
| `--force` | Delete wrapper-managed outputs first, then rerun cleanly |

### Common workflow example

#### Step 1: prepare each input once

Prepare a cohort `.bim` payload:

```bash
prepare_variants.py \
  --input cohort.bim \
  --input-format bim \
  --output work/cohort.prepared
```

Prepare a sharded reference `.bim` payload, split across multiple inputs:

```bash
prepare_variants_sharded.py \
  --input reference.@.bim \
  --input-format bim \
  --prefix work/reference.@.prepared \
  --output work/reference.prepared
```

Prepare a summary-statistics payload:

```bash
prepare_variants.py \
  --input study.tsv.gz \
  --input-format sumstats \
  --sumstats-metadata study.yaml \
  --output work/study.prepared
```

Here `--output` is a stem, so these commands write the final prepared objects to:

- `work/cohort.prepared.vmap`
- `work/reference.prepared.vmap`
- `work/study.prepared.vmap`

`prepare_variants.py` retains stage-specific `.vmap` outputs under a prefix that defaults to `--output`. `prepare_variants_sharded.py` requires an explicit sharded `--prefix` and retains per-contig-group stage-specific `.vmap` outputs by replacing `@` with each target contig group token.

#### Step 2: intersect the prepared outputs

Once the inputs are prepared, intersect the prepared cohort and reference `.vmap` outputs:

```bash
intersect_variants.py \
  work/cohort.prepared.vmap work/reference.prepared.vmap \
  --output work/shared.vtable
```

The output `work/shared.vtable` is the shared membership table. Use `restrict_vmap.py` to make one payload-specific `.vmap` per payload that should be applied.
If you want a merged target set instead of an exact intersection, use `union_variants.py` on the prepared inputs instead.

#### Step 3: restrict each `.vmap` to the shared membership table

```bash
restrict_vmap.py \
  work/study.prepared.vmap work/shared.vtable \
  --output work/study.shared.vmap

restrict_vmap.py \
  work/cohort.prepared.vmap work/shared.vtable \
  --output work/cohort.shared.vmap

restrict_vmap.py \
  work/reference.prepared.vmap work/shared.vtable \
  --output work/reference.shared.vmap
```

Each restricted `.vmap` preserves source-specific row order, IDs, provenance, and `allele_op`.

#### Step 4: project payloads through explicit `.vmap` inputs

Project the summary-statistics payload:

```bash
project_payload.py \
  --input study.tsv.gz \
  --input-format sumstats-clean \
  --sumstats-metadata study.yaml \
  --vmap work/study.shared.vmap \
  --output study.shared.tsv.gz
```
With `--input-format sumstats-clean`, the rewritten output is harmonized TSV.
This clean harmonization step is a missing-value and missing-column completion pipeline. It may infer missing fields from other available fields, but it does not try to verify semantic consistency between overlapping values that are already present.

Project the PLINK payload:

```bash
project_payload.py \
  --input cohort.bim \
  --input-format bfile \
  --vmap work/cohort.shared.vmap \
  --output cohort.shared
```

Project the sharded reference payload:

```bash
project_payload.py \
  --input reference.@.bim \
  --input-format bfile \
  --vmap work/reference.shared.vmap \
  --output reference.shared.@
```

This writes the sharded payload outputs under `reference.shared.<contig>`.

## Core objects and metadata

Each `.vtable` or `.vmap` object is a single file with JSON sidecars at `<object>.meta.json`.

| Object | Purpose | Required metadata |
| --- | --- | --- |
| `.vmap` | Mapped-only provenance-carrying object for import, cleanup, reference-compatible restriction, liftover, `restrict_vmap.py`, and application | `object_type=variant_map`, `target.genome_build`, `target.contig_naming`, `variants_count` |
| `.vtable` | Provenance-free table, produced mainly by `convert_vmap_to_target.py`, `intersect_variants.py`, and `union_variants.py` | `object_type=variant_table`, `genome_build`, `contig_naming`, `variants_count` |

Downstream transforms other than `normalize_contigs.py` require declared `contig_naming`. If an imported object is unresolved or mixed, clean it with `normalize_contigs.py` before continuing. See [core-objects.md](spec/core-objects.md) for the on-disk column contracts.

Newly written `.vmap` and `.vtable` sidecars include top-level `variants_count` row-count metadata. Readers accept missing `variants_count` for backward compatibility. Set-operation tools may use it to choose memory-efficient scan order, but not to change documented output order, ID policy, provenance, or emitted object type.

`a1` and `a2` are ordered alleles, not an unordered allele set. Importers preserve source-format allele meaning at import time, while reference-aware tools such as `restrict_build_compatible.py` and `liftover_build.py` emit reference-aligned variants with `a1=non-reference` and `a2=reference`. This keeps canonical reference-aware rows aligned with ALT/effect/counting-style semantics while still making the reference allele explicit in `a2`.

## Canonical tools

| Task | Tool | Typical input | Emits | Notes |
| --- | --- | --- | --- | --- |
| Import | `import_bim.py` | raw `.bim` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_pvar.py` | raw `.pvar` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_vcf.py` | raw `.vcf` / `.vcf.gz` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_sumstats.py` | raw summary statistics plus metadata YAML | `.vmap` | Creates provenance; supports `--chr2use` / `--contigs` and optional `--id-vtable`; single-file only |
| Normalize / metadata | `guess_build.py` | `.vmap` or `.vtable` | same as input | Metadata-only update; on `.vmap`, updates target-side metadata only; build evidence uses a default random downsample of up to `10,000` rows (`--sample-rows`, use `0` for all rows) |
| Normalize / metadata | `normalize_contigs.py` | `.vmap` or `.vtable` | same as input | Repairs or standardizes target contigs; target naming is selected with `--to` |
| Normalize / metadata | `restrict_contigs.py` | `.vmap` or `.vtable` | same as input | Target-side contig filter via `--chr2use` / `--contigs` |
| Normalize / metadata | `drop_strand_ambiguous.py` | `.vmap` or `.vtable` | same as input | Drops strand-ambiguous target rows |
| Validate / liftover | `restrict_build_compatible.py` | `.vmap` or `.vtable` | same as input | Same-build, reference-aware filtering; key flags are `--allow-strand-flips`, `--norm-indels`, and `--sort` |
| Validate / liftover | `liftover_build.py` | `.vmap` or `.vtable` | same as input | Explicit build conversion to `--target-build`; preserves `.vmap` provenance and re-sorts into declared coordinate order |
| Order  | `sort_variants.py` | `.vmap` or `.vtable` | same as input | Explicit standalone target-row sorting by declared contig order, then numeric position, then `a1`/`a2` lexicographically; optional `--drop-duplicates` drops duplicate target identities after sorting |
| Restrict / intersect / materialize | `intersect_variants.py` | two or more `.vmap` / `.vtable` inputs | `.vtable` | Exact symmetric intersection; output order and IDs come from the first input, provenance is dropped |
| Restrict / intersect / materialize | `restrict_vmap.py` | source `.vmap` plus one or more `.vmap` / `.vtable` restrictions | `.vmap` | Exact source-order restriction by intersection of all restriction inputs; output order, IDs, provenance, and `allele_op` come from the source `.vmap` |
| Restrict / intersect / materialize | `union_variants.py` | two or more `.vmap` / `.vtable` inputs | `.vtable` | Exact source-agnostic target-key union; output is re-sorted into declared coordinate order |
| Restrict / intersect / materialize | `convert_vmap_to_target.py` | `.vmap` | `.vtable` | Materializes the `.vmap` target side and intentionally drops provenance |
| Apply to payloads | `apply_vmap_to_sumstats.py` | `.vmap` plus summary-stat payload | rewritten payload | Single-file payload only; emits rows in `.vmap` order, writes corrected target-based IDs by default with `--retain-snp-id` as a legacy opt-out, writes tab-delimited output plus `<output>.meta.yaml`, and exposes `--clean` for canonical cleaned output |
| Apply to payloads | `apply_vmap_to_bfile.py` | `.vmap` plus PLINK 1 BED/BIM/FAM payload | rewritten payload | Shared genotype-payload `apply_vmap_*` contract; emits rows in `.vmap` order, supports payload `@` sharding, bounded-chunk processing, corrected target-based IDs by default, explicit `--target-fam`, and `--sample-id-mode {fid_iid,iid}` |
| Apply to payloads | `apply_vmap_to_pfile.py` | `.vmap` plus PLINK 2 PGEN/PVAR/PSAM payload | rewritten payload | Shared genotype-payload `apply_vmap_*` contract; emits rows in `.vmap` order, supports payload `@` sharding, bounded-chunk processing, corrected target-based IDs by default, explicit `--target-psam`, `--sample-id-mode {fid_iid,iid}`, preserves supported biallelic hardcalls/phase/unphased dosages, and emits `.pgen/.pvar/.psam` |

### Using .vmap vs .vtable as input across tools

- `import_*` tools are the only tools that originate new provenance, so they emit `.vmap`.
- For tools that accept both `.vmap` and `.vtable`, the default rule is to act on target rows and emit the same artifact type as the input, unless the tool explicitly defines a provenance-free `.vtable` output.
- `intersect_variants.py` always emits provenance-free `.vtable`.
- `restrict_vmap.py` is the provenance-preserving exact restriction tool for `.vmap`; restriction inputs are membership filters only and multiple restriction inputs are combined by intersection.
- `apply_vmap_*` consumes `.vmap` as a mapping object, not as a table. The payload must still match the `.vmap` source provenance exactly.

### Important usage notes

- Import tools always preserve source provenance for retained rows, also when it drops unrepresentable rows or rows filtered by `--chr2use` / `--contigs`
- import tools do not repair mixed or invalid contig labels. Use `normalize_contigs.py` for cleanup.
- `restrict_build_compatible.py` is the same-build, reference-aware canonicalization step: it filters to reference-compatible rows, can optionally re-sort into declared coordinate order, and can optionally normalize supported indels via retained `bcftools norm` intermediates while still emitting canonical biallelic `a1=non-reference`, `a2=reference` rows; it does not liftover.
- `liftover_build.py` supports SNVs and one-side-multibase reference-anchored rows, but not lifted outputs where both final alleles remain multi-base.
- `intersect_variants.py` preserves first-input order and emits `.vtable`. `restrict_vmap.py` preserves source `.vmap` order and emits `.vmap`. If you need sorted output, sort the desired first input or source `.vmap` before intersecting or restricting. `union_variants.py` emits a provenance-free target-key union with target-derived IDs (`chrom:pos:a1:a2`) in declared coordinate order. `liftover_build.py` and `sort_variants.py` use declared coordinate order: declared contig order, then numeric position, then `a1`/`a2` lexicographically, with stable input order for exact `(chrom, pos, a1, a2)` ties.
- `intersect_variants.py`, `restrict_vmap.py`, and `union_variants.py` may use `variants_count` to start internal set operations from smaller objects. This is a memory strategy only; the output contracts above remain unchanged.
- If any input to those set-operation tools lacks `variants_count`, the tool warns and uses the user-declared input order for its internal scan plan.
- coordinate-changing transforms and genotype-payload application follow the shared expected-ploidy contract in [spec/ploidy-model.md](spec/ploidy-model.md).
- `apply_vmap_*` and `project_payload.py` write corrected final payload IDs as `chrom:pos:a1:a2` from `.vmap` target rows by default. Use `--retain-snp-id` only when the projected payload must preserve target-side `.vmap` IDs.
- genotype-payload `apply_vmap_*` tools apply `.vmap` provenance back to the original payload, support `@` source discovery and `@` target-contig output sharding, and follow the shared ploidy contract.
- both genotype-payload `apply_vmap_*` tools support explicit target sample files plus `--sample-id-mode {fid_iid,iid}` for subject restriction, reordering, and reconciliation. Explicit target sample files define output subject order exactly and drive reconciliation-added missingness reporting.
- both genotype-payload `apply_vmap_*` tools support `--sample-axis native` to keep each emitted output shard on its native source sample axis, provided every emitted shard has one unambiguous source sample axis.
- both genotype-payload `apply_vmap_*` tools support `--skip-ploidy-check` to skip secondary genotype-content ploidy validation; BFILE `.ploidy` output is still emitted according to the ploidy model.
- `apply_vmap_to_bfile.py` uses `--target-fam`; `apply_vmap_to_pfile.py` uses `--target-psam`. Without explicit target sample files or `--sample-axis native`, sharded source payloads must still have identical sample metadata across referenced shards.
- `project_payload.py --sample-axis union` is a wrapper-only convenience for sharded genotype payloads. It unions subject keys across exactly the source shards referenced by the input `.vmap`, synthesizes and retains `<output-prefix>.target_samples.fam` or `<output-prefix>.target_samples.psam`, then calls the canonical apply tool with that explicit target sample file. For non-sharded input, or when only one referenced shard remains, it warns and is a no-op.
- `@` sharding is supported for heavy raw inputs and genotype payloads; summary-stat payloads, `.vmap`, and `.vtable` objects are single-file.
- for filename-based `@` discovery, shard processing order follows deterministic lexical path order; this order is authoritative for first-wins behaviors in workflows that deduplicate duplicates across shards (for example, lexical labels like `chr10` sort before `chr2`).

### Usage example for canonical tools

Use case: liftover a summary-statistics file to `GRCh38` before downstream intersection or application.

```bash
import_sumstats.py \
  --input study.tsv.gz \
  --sumstats-metadata study_meta.yaml \
  --output study.raw.vmap

normalize_contigs.py --to ncbi \
  --input study.raw.vmap \
  --output study.clean.vmap

guess_build.py \
  --input study.clean.vmap \
  --write

restrict_build_compatible.py \
  --source study.clean.vmap \
  --allow-strand-flips \
  --norm-indels \
  --output study.validated.vmap

liftover_build.py --target-build GRCh38 \
  --input study.validated.vmap \
  --output study.grch38.vmap

apply_vmap_to_sumstats.py \
  --input study.tsv.gz \
  --sumstats-metadata study_meta.yaml \
  --vmap study.grch38.vmap \
  --output study.grch38.tsv.gz
```

Use case: build one shared `GRCh38` variant universe for chromosome-sharded cohort genotypes, a target reference panel such as 1000 Genomes, and the lifted summary statistics from the first workflow, then emit a constrained cohort PLINK dataset plus aligned summary statistics in that same universe.

```bash
import_vcf.py \
  --input kgp.vcf.gz \
  --output kgp.raw.vmap \
  --genome-build GRCh37

normalize_contigs.py --to ncbi \
  --input kgp.raw.vmap \
  --output kgp.clean.vmap

restrict_build_compatible.py \
  --source kgp.clean.vmap \
  --allow-strand-flips \
  --norm-indels \
  --output kgp.validated.vmap

liftover_build.py --target-build GRCh38 \
  --input kgp.validated.vmap \
  --output kgp.grch38.vmap

convert_vmap_to_target.py \
  --source kgp.grch38.vmap \
  --output kgp.grch38.vtable

import_bim.py \
  --input cohort.@.bim \
  --output cohort.raw.vmap \
  --genome-build GRCh38

normalize_contigs.py --to ncbi \
  --input cohort.raw.vmap \
  --output cohort.clean.vmap

intersect_variants.py \
  cohort.clean.vmap kgp.grch38.vtable \
  --output shared.grch38.vtable

restrict_vmap.py \
  cohort.clean.vmap shared.grch38.vtable \
  --output cohort.shared.vmap

apply_vmap_to_bfile.py \
  --source-prefix cohort.@ \
  --vmap cohort.shared.vmap \
  --output-prefix cohort.shared.@

restrict_vmap.py \
  study.grch38.vmap shared.grch38.vtable \
  --output study.shared.vmap

apply_vmap_to_sumstats.py \
  --input study.tsv.gz \
  --sumstats-metadata study_meta.yaml \
  --vmap study.shared.vmap \
  --output study.shared.tsv.gz
```

For a PLINK 2 payload, replace the final genotype-application step with `apply_vmap_to_pfile.py`; it follows the same shared row-selection, provenance, and sharding contract, but writes `.pgen/.pvar/.psam` and preserves supported biallelic hardcalls, hardcall phase, and unphased dosages.
