# Workflow

Most users should use the workflow-level tools:

- [`prepare_variants.py`](#prepare-variants) prepares one raw input into a final `.vmap`.
- [`prepare_variants_sharded.py`](#prepare-sharded-variants) prepares chromosome-sharded genotype or VCF inputs into one final `.vmap`.
- [`intersect_variants.py`](#intersect-variants) derives a shared provenance-free membership `.vtable` from prepared inputs.
- [`union_variants.py`](#union-variants) derives a source-agnostic union membership `.vtable` from prepared inputs.
- [`restrict_vmap.py`](#restrict-vmap) filters one payload-specific `.vmap` to a membership set.
- [`project_payload.py`](#project-payloads) rewrites an original payload through an explicit mapped-only `.vmap`.

For a worked example, start with [tutorial-1.md](tutorial-1.md). Use the tables below as a quick wrapper reference. The detailed tool and wrapper contracts remain in the spec.

## Minimum mental model

- `prepare_variants.py` turns each raw input into one prepared `.vmap`.
- A `.vmap` is the prepared variant table plus the source provenance needed to project the same original payload later. It does not contain the payload itself.
- `intersect_variants.py` and `union_variants.py` combine prepared variant sets and emit a `.vtable`.
- A `.vtable` is only a shared membership table. It has no source-payload provenance, so it cannot be projected directly.
- Use `restrict_vmap.py` to filter one payload-specific `.vmap` to a shared membership table, producing a restricted `.vmap` for that payload.
- `project_payload.py` rewrites the original payload through its own restricted `.vmap`.
- Variant matching uses prepared `chr:bp:a1:a2` target rows, not variant IDs.
- Choose the target genome build and contig naming up front; reference-aware preparation needs `MATCH_CONFIG`.
- `@` in input and output paths denotes sharded genotype or VCF paths.
- Summary-statistics inputs require metadata YAML; `sumstats-clean` projection writes harmonized canonical summary statistics.

For deeper object-model details, see [primitives.md](primitives.md#object-model-reference).

## Prepare variants

`prepare_variants.py` is the convenience wrapper for preparing one raw input into a final `<output>.vmap` that can later be intersected or projected back onto the original payload. It orchestrates import, contig normalization, metadata resolution, same-build restriction or liftover, and optional strand-ambiguous / contig filtering, while retaining the intermediate `.vmap` stages. For the exact stage graph and `--resume` / `--force` behavior, see [spec/workflow.md](../spec/workflow.md).

Illustrative stage flow:

| Stage | Description | Retained output |
| --- | --- | --- |
| Import | Always runs `import_<format>.py`; passes importer controls such as `--max-allele-length`, `--shards`, and sumstats-only `--id-vtable` when supplied | `<prefix>.imported.vmap` |
| Contig normalization | Runs `normalize_contigs.py --to <dst-contig-naming>` when the current object omits `contig_naming` or differs from `--dst-contig-naming`; final `plink_splitx` can be deferred until build is known | `<prefix>.normalized.vmap` |
| Build resolution | Runs `guess_build.py` in place, unless `--resume` can skip because build is already resolved | current retained `.vmap` |
| Same-build reference filtering | Always runs before any optional liftover; invokes `restrict_build_compatible.py`, passing `--allow-strand-flips` unless disabled, `--norm-indels` unless disabled, and `--sort --drop-duplicates` when current build already matches `--dst-build` | `<prefix>.build_compatible.vmap` |
| Liftover | Runs `liftover_build.py --target-build <dst-build>` when the current build differs from `--dst-build` after same-build filtering | `<prefix>.lifted.vmap` |
| Deferred split-X normalization | Runs `normalize_contigs.py --to plink_splitx` when `--dst-contig-naming=plink_splitx` and final split-X normalization was deferred until build resolution | `<prefix>.splitx.vmap` |
| Strand filter | Runs `drop_strand_ambiguous.py` when `--drop-strand-ambiguous` is supplied | `<prefix>.strand_filtered.vmap` |
| Contig filter | Runs `restrict_contigs.py --chr2use/--contigs <value>` when `--chr2use` / `--contigs` is supplied | `<prefix>.contigs.vmap` |
| Final copy | Always copies the last retained stage | `<output>.vmap` |

| `prepare_variants.py` argument | Meaning |
| --- | --- |
| `--input-format` | Required importer format: `bim`, `pvar`, `vcf`, or `sumstats`. See [spec/importers.md](../spec/importers.md) |
| `--input`, `--output` | Required raw input payload and output stem. Genotype payload inputs may be sharded via `@`; see [spec/shard-discovery.md](../spec/shard-discovery.md). The final prepared output is always a single `<output>.vmap` |
| `--sumstats-metadata`, `--id-vtable` | Sumstats-only controls: `--sumstats-metadata` is required for `--input-format sumstats` and follows the raw metadata schema at [schemas/raw-sumstats-metadata.yaml](../src/genomatch/schemas/raw-sumstats-metadata.yaml), which is the same specification used by [BioPsyk/cleansumstats](https://github.com/BioPsyk/cleansumstats). `--id-vtable` optionally fills missing `chr`/`pos` by matching on variant ID, and becomes required when those fields are absent from the summary-stat metadata |
| `--dst-contig-naming` | Target contig naming; defaults to `ncbi`. Supported values are `ncbi`, `ucsc`, `plink`, and `plink_splitx`; `plink_splitx` follows PLINK `--split-x` style X/XY_PAR labeling. See [spec/contigs-and-metadata.md](../spec/contigs-and-metadata.md) |
| `--dst-build` | Target genome build (`GRCh37`, `GRCh38`, or `T2T-CHM13v2.0`); defaults to `GRCh38`. If the input build differs, liftover rewrites both coordinates and alleles via `bcftools +liftover`, then emits canonical target rows with `a2=reference`. See [spec/variant-transforms.md](../spec/variant-transforms.md) |
| `--[no-]allow-strand-flips` | Control strand-flip allowance during same-build reference-aware restriction; enabled by default |
| `--[no-]norm-indels` | Control indel normalization during same-build reference-aware restriction; `bcftools norm` is used internally when enabled |
| `--drop-strand-ambiguous` | Optional removal of strand-ambiguous variants after build compatibility / liftover |
| `--chr2use` / `--contigs` | Optional restriction of final prepared output to selected contigs. `--contigs` is an alias for `--chr2use`. For non-autosomes and ploidy-related behavior, see [spec/ploidy-model.md](../spec/ploidy-model.md) |
| `--prefix` | Optional retained-intermediate stem; defaults to `--output` |
| `--resume`, `--force` | Wrapper execution controls; mutually exclusive |

## Prepare sharded variants

`prepare_variants_sharded.py` is the memory-bounded sharded-input variant of `prepare_variants.py`. It accepts `bim`, `pvar`, and `vcf` inputs whose `--input` contains `@`, groups selected input shards by unambiguous target contig, runs the full `prepare_variants.py` stage graph once per target-contig group, and concatenates the per-contig-group final `.vmap` files into one final `<output>.vmap`. Its stage semantics are inherited from `prepare_variants.py`; see [spec/workflow.md](../spec/workflow.md).

| `prepare_variants_sharded.py` argument | Meaning |
| --- | --- |
| `--input`, `--input-format` | Required sharded raw input. `--input` must contain `@`; `--input-format` may be `bim`, `pvar`, or `vcf` |
| `--output` | Required non-sharded output stem. The final prepared output is `<output>.vmap` |
| `--prefix` | Required retained-intermediate stem and must contain `@`; per-contig-group retained outputs replace `@` with the target contig group token |
| Other preparation flags | Same meaning as `prepare_variants.py` and passed through to each per-contig-group invocation |
| `--resume`, `--force` | Wrapper execution controls; mutually exclusive |

## Intersect variants

`intersect_variants.py` accepts two or more prepared `.vmap` / `.vtable` inputs and emits a provenance-free `.vtable` containing exact `chrom:pos:a1:a2` keys present in every input. Output IDs are target-derived `chrom:pos:a1:a2`, and output rows are in declared coordinate order.

## Union variants

`union_variants.py` accepts two or more prepared `.vmap` / `.vtable` inputs and emits a provenance-free `.vtable` containing the source-agnostic exact-key union. Output IDs are target-derived `chrom:pos:a1:a2`, and output rows are in declared coordinate order.

## Restrict vmap

`restrict_vmap.py` accepts one source `.vmap` plus one or more `.vmap` / `.vtable` membership inputs and emits a restricted `.vmap`. Output order, IDs, provenance, and `allele_op` come from the source `.vmap`; restriction inputs affect membership only.

Membership inputs can be either `.vtable` or `.vmap` files. When a `.vmap` is used as membership, only its target-side `chrom:pos:a1:a2` rows are considered; its source provenance and `allele_op` are ignored.

## Project payloads

`project_payload.py` is the convenience wrapper for applying an original payload through an explicit prepared `.vmap`. It requires `--vmap`, does not accept a target `.vtable`, and does not run `restrict_vmap.py`; users who want to apply a payload to a shared membership table must first run `restrict_vmap.py` themselves. Payload row order is always `.vmap` row order. For genotype payloads, the wrapper can also reconcile subject axes via explicit target sample files or wrapper-level union synthesis. For the exact payload-application contract, output rules, and target-sample reconciliation semantics, see [spec/workflow.md](../spec/workflow.md).

| `project_payload.py` argument | Meaning |
| --- | --- |
| `--input-format` | Required payload format: `bfile`, `pfile`, `sumstats`, or `sumstats-clean` (`bfile` maps to `bim` from `prepare_variants.py`; `pfile` maps to prepared `pvar`). See [spec/payload-application.md](../spec/payload-application.md) |
| `--input`, `--output` | Required raw payload to project and rewritten payload destination. Genotype payload `--input` may be sharded via `@`; for `bfile` / `pfile`, `--output` is a PLINK output prefix and may also be sharded via `@` |
| `--vmap` | Required mapped-only `.vmap` to apply; it defines output variant rows, row order, IDs, source provenance, and `allele_op` |
| `--retain-snp-id` | Optional payload projection control. By default, `project_payload.py` writes corrected output IDs as `chrom:pos:a1:a2` from `.vmap` target rows; `--retain-snp-id` writes target-side `.vmap` IDs instead |
| `--sumstats-metadata`, `--fill-mode`, `--use-af-inference` | Summary-stat controls: `--sumstats-metadata` is required for `sumstats` and `sumstats-clean`, and `--fill-mode` / `--use-af-inference` apply to `sumstats-clean`. See [spec/sumstats-harmonization.md](../spec/sumstats-harmonization.md) |
| `--target-fam` / `--target-psam` / `--sample-axis {union,native}` / `--sample-id-mode {fid_iid,iid}` / `--skip-ploidy-check` | Optional genotype payload subject-axis controls for explicit target sample files, wrapper-synthesized union targets, native per-output-shard sample axes, subject keying, and ploidy-validation skipping. See [Sample-axis reconciliation for genotype payloads](../spec/payload-application.md#sample-axis-reconciliation-for-genotype-payloads) |
| `--force` | Delete wrapper-managed outputs first, then rerun cleanly |
