# Primitives

The primitive tools provide explicit control over import, normalization, liftover, intersection, and payload application. Most users should start with [workflow.md](workflow.md). For summary-statistics metadata, SNP-only imports, and clean projection examples, see [sumstats.md](sumstats.md).

## Object model reference

### Core mental model

Conceptually, the toolkit treats raw input as three pieces that can be separated and recombined: the variant table itself (`.vtable`), the attached payload, and the provenance linking variants back to the raw source:

- `raw data = .vtable + payload + provenance`
- `.vtable` is the ordered variant table, one row per biallelic variant
- `payload` is the external format-specific data attached to those variants, such as summary-stat columns, PLINK 1 BED/BIM/FAM genotype matrices, or PLINK 2 PGEN/PVAR/PSAM genotype matrices; payloads can be split across multiple files (shards)
- `provenance` is the link back from variant to its raw payload, represented as `(source_shard, source_index)` plus `allele_op` to capture `a1 <-> a2` allele swaps and strand flips
- `.vmap = .vtable + provenance`; `.vmap` is therefore the object emitted by import operations, and used for most downstream transforms: it keeps transformed table of mapped variants together with enough provenance to rewrite the original payload later, but it does not contain the payload itself. A canonical `.vmap` is mapped-only; rows missing from a source payload are absent rather than represented by sentinel rows.

The toolkit provides a set of tools for transforming tables of variants: importing raw inputs, normalizing representation, restricting rows, liftover between builds, and performing intersections and unions of variant tables. These transformations can carry provenance linking each transformed variant row back to the original payload. That is what `.vmap` stores: a transformed variant table together with provenance, but without the payload itself. A `.vtable` stores only the transformed variant table, without provenance. `apply_vmap_*` reconnects the transformed variant rows in `.vmap` to the original payload and writes a rewritten payload in that transformed row order.

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

For exact schema and edge-case rules, see [SPEC.md](../SPEC.md), [core-objects.md](../spec/core-objects.md), [contigs-and-metadata.md](../spec/contigs-and-metadata.md), [shard-discovery.md](../spec/shard-discovery.md), [importers.md](../spec/importers.md), [variant-transforms.md](../spec/variant-transforms.md), [mapping.md](../spec/mapping.md), [payload-application.md](../spec/payload-application.md), and [workflow.md](../spec/workflow.md). End-user installation is described in [install.md](install.md), and reference-aware assets/config in [downloads.md](downloads.md).

### Core objects and metadata

Each `.vtable` or `.vmap` object is a single file with JSON sidecars at `<object>.meta.json`.

| Object | Purpose | Required metadata |
| --- | --- | --- |
| `.vmap` | Mapped-only provenance-carrying object for import, cleanup, reference-compatible restriction, liftover, `restrict_vmap.py`, and application | `object_type=variant_map`, `target.genome_build`, `target.contig_naming`, `variants_count` |
| `.vtable` | Provenance-free table, produced mainly by `convert_vmap_to_target.py`, `intersect_variants.py`, and `union_variants.py` | `object_type=variant_table`, `genome_build`, `contig_naming`, `variants_count` |

Downstream transforms other than `normalize_contigs.py` require declared `contig_naming`. If an imported object is unresolved or mixed, clean it with `normalize_contigs.py` before continuing. See [core-objects.md](../spec/core-objects.md) for the on-disk column contracts.

Newly written `.vmap` and `.vtable` sidecars include top-level `variants_count` row-count metadata. Readers accept missing `variants_count` for backward compatibility. Set-operation tools may use it to choose memory-efficient scan order, but not to change documented output order, ID policy, provenance, or emitted object type.

`a1` and `a2` are ordered alleles, not an unordered allele set. Importers preserve source-format allele meaning at import time, while reference-aware tools such as `restrict_build_compatible.py` and `liftover_build.py` emit reference-aligned variants with `a1=non-reference` and `a2=reference`. This keeps canonical reference-aware rows aligned with ALT/effect/counting-style semantics while still making the reference allele explicit in `a2`.

### Using .vmap vs .vtable as input across tools

- `import_*` tools are the only tools that originate new provenance, so they emit `.vmap`.
- For tools that accept both `.vmap` and `.vtable`, the default rule is to act on target rows and emit the same artifact type as the input, unless the tool explicitly defines a provenance-free `.vtable` output.
- `intersect_variants.py` always emits provenance-free `.vtable`.
- `restrict_vmap.py` is the provenance-preserving exact restriction tool for `.vmap`; restriction inputs are membership filters only and multiple restriction inputs are combined by intersection.
- `apply_vmap_*` consumes `.vmap` as a mapping object, not as a table. The payload must still match the `.vmap` source provenance exactly.

## Primitive tools

| Task | Tool | Typical input | Emits | Notes |
| --- | --- | --- | --- | --- |
| Import | `import_bim.py` | raw `.bim` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_pvar.py` | raw `.pvar` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_vcf.py` | raw `.vcf` / `.vcf.gz` | `.vmap` | Creates provenance; supports raw-input `@` discovery and `--chr2use` / `--contigs` |
| Import | `import_sumstats.py` | raw summary statistics plus metadata YAML | `.vmap` | Creates provenance; supports `--chr2use` / `--contigs` and optional `--id-lookup`; single-file only |
| Normalize / metadata | `guess_build.py` | `.vmap` or `.vtable` | same as input | Metadata-only update; on `.vmap`, updates variant-table metadata only; build evidence uses a default random downsample of up to `10,000` rows (`--sample-rows`, use `0` for all rows) |
| Normalize / metadata | `normalize_contigs.py` | `.vmap` or `.vtable` | same as input | Repairs or standardizes target contigs; target naming is selected with `--to` |
| Normalize / metadata | `restrict_contigs.py` | `.vmap` or `.vtable` | same as input | Contig filter via `--chr2use` / `--contigs` |
| Normalize / metadata | `drop_strand_ambiguous.py` | `.vmap` or `.vtable` | same as input | Drops strand-ambiguous target rows |
| Validate / liftover | `restrict_build_compatible.py` | `.vmap` or `.vtable` | same as input | Same-build, reference-aware filtering; key flags are `--allow-strand-flips`, `--norm-indels`, and `--sort` |
| Validate / liftover | `liftover_build.py` | `.vmap` or `.vtable` | same as input | Explicit build conversion to `--target-build`; preserves `.vmap` provenance and re-sorts into declared coordinate order |
| Order  | `sort_variants.py` | `.vmap` or `.vtable` | same as input | Explicit standalone target-row sorting by declared contig order, then numeric position, then `a1`/`a2` lexicographically; optional `--drop-duplicates` drops duplicate target identities after sorting |
| Assign IDs | `assign_vmap_ids.py` | source `.vmap` plus `.vmap` / `.vtable` ID source | `.vmap` | Copies IDs from the ID source by exact `chrom:pos:a1:a2` match, preserves source `.vmap` provenance, and by default drops rows without a non-missing assigned ID |
| Restrict / intersect / materialize | `intersect_variants.py` | two or more `.vmap` / `.vtable` inputs | `.vtable` | Exact symmetric source-agnostic target-key intersection; output has target-derived IDs in declared coordinate order |
| Restrict / intersect / materialize | `restrict_vmap.py` | source `.vmap` plus one or more `.vmap` / `.vtable` restrictions | `.vmap` | Exact source-order restriction by intersection of all restriction inputs; output order, IDs, provenance, and `allele_op` come from the source `.vmap` |
| Restrict / intersect / materialize | `union_variants.py` | two or more `.vmap` / `.vtable` inputs | `.vtable` | Exact source-agnostic target-key union; output is re-sorted into declared coordinate order |
| Restrict / intersect / materialize | `convert_vmap_to_target.py` | `.vmap` | `.vtable` | Materializes the `.vmap` rows as a `.vtable` and intentionally drops provenance |
| Apply to payloads | `apply_vmap_to_sumstats.py` | `.vmap` plus summary-stat payload | rewritten payload | Single-file payload only; emits rows in `.vmap` order, writes corrected target-based IDs by default with `--retain-snp-id` as a legacy opt-out, writes tab-delimited output plus `<output>.meta.yaml`, and exposes `--clean` for canonical cleaned output |
| Apply to payloads | `apply_vmap_to_bfile.py` | `.vmap` plus PLINK 1 BED/BIM/FAM payload | rewritten payload | Shared genotype-payload `apply_vmap_*` contract; emits rows in `.vmap` order, supports payload `@` sharding, bounded-chunk processing, corrected target-based IDs by default, explicit `--target-fam`, and `--sample-id-mode {fid_iid,iid}` |
| Apply to payloads | `apply_vmap_to_pfile.py` | `.vmap` plus PLINK 2 PGEN/PVAR/PSAM payload | rewritten payload | Shared genotype-payload `apply_vmap_*` contract; emits rows in `.vmap` order, supports payload `@` sharding, bounded-chunk processing, corrected target-based IDs by default, explicit `--target-psam`, `--sample-id-mode {fid_iid,iid}`, preserves supported biallelic hardcalls/phase/unphased dosages, and emits `.pgen/.pvar/.psam` |

## Important usage notes

- Import tools always preserve source provenance for retained rows, also when it drops unrepresentable rows or rows filtered by `--chr2use` / `--contigs`
- import tools do not repair mixed or invalid contig labels. Use `normalize_contigs.py` for cleanup.
- `restrict_build_compatible.py` is the same-build, reference-aware canonicalization step: it filters to reference-compatible rows, can optionally re-sort into declared coordinate order, and can optionally normalize supported indels via retained `bcftools norm` intermediates while still emitting canonical biallelic `a1=non-reference`, `a2=reference` rows; it does not liftover.
- `liftover_build.py` supports SNVs and one-side-multibase reference-anchored rows, but not lifted outputs where both final alleles remain multi-base.
- `intersect_variants.py` and `union_variants.py` emit provenance-free target-key sets with target-derived IDs (`chrom:pos:a1:a2`) in declared coordinate order. `restrict_vmap.py` preserves source `.vmap` order and emits `.vmap`; if you need sorted restricted output, sort the source `.vmap` before restricting. `liftover_build.py` and `sort_variants.py` use declared coordinate order: declared contig order, then numeric position, then `a1`/`a2` lexicographically, with stable input order for exact `(chrom, pos, a1, a2)` ties.
- `intersect_variants.py`, `restrict_vmap.py`, and `union_variants.py` may use `variants_count` to start internal set operations from smaller objects. This is a memory strategy only; the output contracts above remain unchanged.
- `assign_vmap_ids.py` is useful after `restrict_vmap.py` and before projection when a final payload-specific `.vmap` should carry rsIDs or another preferred ID namespace. It matches prepared coordinates to IDs, not IDs to coordinates; use `import_sumstats.py --id-lookup` for preparation-time coordinate recovery from raw IDs. By default it drops rows without an ID-source match; use `--unmatched-id-policy variant-key` or `--unmatched-id-policy missing` to retain them with generated `chrom:pos:a1:a2` IDs or `.` IDs. Use `project_payload.py --retain-snp-id` to project assigned IDs into the payload.
- If any input to `intersect_variants.py` or `union_variants.py` lacks `variants_count`, the tool warns and uses the user-declared input order for its internal scan plan. `restrict_vmap.py` may use the same fallback silently because missing counts do not affect observable restriction semantics.
- coordinate-changing transforms and genotype-payload application follow the shared expected-ploidy contract in [spec/ploidy-model.md](../spec/ploidy-model.md).
- `apply_vmap_*` and `project_payload.py` write corrected final payload IDs as `chrom:pos:a1:a2` from `.vmap` rows by default. Use `--retain-snp-id` only when the projected payload must preserve `.vmap` IDs.
- genotype-payload `apply_vmap_*` tools apply `.vmap` provenance back to the original payload, support `@` source discovery and `@` target-contig output sharding, and follow the shared ploidy contract.
- both genotype-payload `apply_vmap_*` tools support explicit target sample files plus `--sample-id-mode {fid_iid,iid}` for subject restriction, reordering, and reconciliation. Explicit target sample files define output subject order exactly and drive reconciliation-added missingness reporting.
- both genotype-payload `apply_vmap_*` tools support `--sample-axis native` to keep each emitted output shard on its native source sample axis, provided every emitted shard has one unambiguous source sample axis.
- both genotype-payload `apply_vmap_*` tools support `--skip-ploidy-check` to skip secondary genotype-content ploidy validation; BFILE `.ploidy` output is still emitted according to the ploidy model.
- `apply_vmap_to_bfile.py` uses `--target-fam`; `apply_vmap_to_pfile.py` uses `--target-psam`. Without explicit target sample files or `--sample-axis native`, sharded source payloads must still have identical sample metadata across referenced shards.
- `project_payload.py --sample-axis union` is a wrapper-only convenience for sharded genotype payloads. It unions subject keys across exactly the source shards referenced by the input `.vmap`, synthesizes and retains `<output-prefix>.target_samples.fam` or `<output-prefix>.target_samples.psam`, then calls the canonical apply tool with that explicit target sample file. For non-sharded input, or when only one referenced shard remains, it warns and is a no-op.
- `@` sharding is supported for heavy raw inputs and genotype payloads; summary-stat payloads, `.vmap`, and `.vtable` objects are single-file.
- for filename-based `@` discovery, shard processing order follows deterministic lexical path order; this order is authoritative for first-wins behaviors in workflows that deduplicate duplicates across shards (for example, lexical labels like `chr10` sort before `chr2`).

## Usage example for primitive tools

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
