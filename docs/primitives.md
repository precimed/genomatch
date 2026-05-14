# Primitives

The primitive tools provide explicit control over import, normalization, liftover, intersection, and payload application. Most users should start with [workflows.md](workflows.md).

For the `.vtable` / `.vmap` object model and source-row mapping details behind these tools, see [concepts.md](concepts.md).

## Primitive tools

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

## Important usage notes

- Import tools always preserve source provenance for retained rows, also when it drops unrepresentable rows or rows filtered by `--chr2use` / `--contigs`
- import tools do not repair mixed or invalid contig labels. Use `normalize_contigs.py` for cleanup.
- `restrict_build_compatible.py` is the same-build, reference-aware canonicalization step: it filters to reference-compatible rows, can optionally re-sort into declared coordinate order, and can optionally normalize supported indels via retained `bcftools norm` intermediates while still emitting canonical biallelic `a1=non-reference`, `a2=reference` rows; it does not liftover.
- `liftover_build.py` supports SNVs and one-side-multibase reference-anchored rows, but not lifted outputs where both final alleles remain multi-base.
- `intersect_variants.py` preserves first-input order and emits `.vtable`. `restrict_vmap.py` preserves source `.vmap` order and emits `.vmap`. If you need sorted output, sort the desired first input or source `.vmap` before intersecting or restricting. `union_variants.py` emits a provenance-free target-key union with target-derived IDs (`chrom:pos:a1:a2`) in declared coordinate order. `liftover_build.py` and `sort_variants.py` use declared coordinate order: declared contig order, then numeric position, then `a1`/`a2` lexicographically, with stable input order for exact `(chrom, pos, a1, a2)` ties.
- `intersect_variants.py`, `restrict_vmap.py`, and `union_variants.py` may use `variants_count` to start internal set operations from smaller objects. This is a memory strategy only; the output contracts above remain unchanged.
- If any input to those set-operation tools lacks `variants_count`, the tool warns and uses the user-declared input order for its internal scan plan.
- coordinate-changing transforms and genotype-payload application follow the shared expected-ploidy contract in [spec/ploidy-model.md](../spec/ploidy-model.md).
- `apply_vmap_*` and `project_payload.py` write corrected final payload IDs as `chrom:pos:a1:a2` from `.vmap` target rows by default. Use `--retain-snp-id` only when the projected payload must preserve target-side `.vmap` IDs.
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
