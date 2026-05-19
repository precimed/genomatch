# Summary Statistics

Summary-statistics inputs use the same prepare, restrict, and project model as genotype payloads, with one extra requirement: a metadata YAML file describes the raw columns. The metadata follows the raw summary-statistics schema at [`src/genomatch/schemas/raw-sumstats-metadata.yaml`](../src/genomatch/schemas/raw-sumstats-metadata.yaml), which is compatible with the cleansumstats-style column naming contract used by this toolkit.

Most users should use the workflow wrappers:

- `prepare_variants.py --input-format sumstats` imports and prepares the variant rows into a `.vmap`.
- `restrict_vmap.py` filters that `.vmap` to a target membership set when needed.
- `project_payload.py --input-format sumstats` or `sumstats-clean` rewrites the original summary-statistics payload through the prepared `.vmap`.

For exact edge-case rules, see [`spec/importers.md`](../spec/importers.md), [`spec/workflow.md`](../spec/workflow.md), [`spec/payload-application.md`](../spec/payload-application.md), and [`spec/sumstats-harmonization.md`](../spec/sumstats-harmonization.md).

## Metadata Basics

The importer and projector resolve raw columns through metadata keys such as:

```yaml
cleansumstats_metafile_kind: minimal
path_sumStats: study.tsv.gz
delimiter: "\t"
missing_value: ""
genome_build: GRCh37
contig_naming: ncbi
stats_Model: logistic

col_CHR: CHR
col_POS: POS
col_SNP: SNP
col_EffectAllele: EA
col_OtherAllele: OA
col_OR: OR
col_P: P
col_CaseN: N_cases
col_ControlN: N_controls
```

`col_EffectAllele` and `col_OtherAllele` are required. Variant coordinates can be provided either as `col_CHR` plus `col_POS`, or through `col_SNP` with `--id-lookup` as described below. At least one payload/stat column, such as `col_BETA`, `col_OR`, `col_SE`, `col_Z`, `col_P`, sample-size columns, frequency columns, `col_INFO`, or `col_Direction`, must be present for projection.

If `--input` is omitted for `prepare_variants.py` or `project_payload.py`, the raw summary-statistics path is read from `path_sumStats` and resolved relative to the metadata file directory. `path_sumStats` is a filename-only field, so keep the metadata and raw sumstats file in the same directory for that mode.

Summary-statistics payloads are single-file inputs. They do not support `@` sharding.

## Prepare With CHR and POS

Use the normal path when the raw file has coordinate columns, or joined coordinate columns that the metadata points to as `col_CHR` and `col_POS`.

```bash
prepare_variants.py \
  --input study.tsv.gz \
  --input-format sumstats \
  --sumstats-metadata study.meta.yaml \
  --dst-build GRCh38 \
  --dst-contig-naming ncbi \
  --output work/study
```

This writes `work/study.vmap`. During preparation the wrapper imports the raw coordinates, normalizes contigs as needed, resolves or validates the build, applies reference-aware filtering, and lifts to the requested destination build when needed.

You may omit `--input` when `path_sumStats` in the metadata names the raw file:

```bash
prepare_variants.py \
  --input-format sumstats \
  --sumstats-metadata study.meta.yaml \
  --dst-build GRCh38 \
  --dst-contig-naming ncbi \
  --output work/study
```

## Import and Preparation Controls

`prepare_variants.py --input-format sumstats` supports these summary-statistics controls:

- `--id-lookup` enables the SNP-only coordinate lookup mode described above.
- `--src-build` declares the source genome build when it is already known and skips build guessing after import. In `--id-lookup` mode, the lookup object's target metadata still defines the imported build and contig naming.
- `--max-allele-length` is passed to the importer and drops rows whose effect or other allele exceeds the configured length. The default is `150`.
- `--chr2use` / `--contigs` restrict the final prepared `.vmap` to selected contigs.

## Prepare SNP-Only Files With ID Lookup

Some summary-statistics files have a variant ID column, often `SNP` or `rsid`, but no chromosome or base-pair columns. In that case, provide an auxiliary prepared or restricted `.vmap` with coordinates and IDs, and pass it with `--id-lookup`.

The metadata for this mode must define `col_SNP`, `col_EffectAllele`, and `col_OtherAllele`, and must omit `col_CHR` and `col_POS`:

```yaml
cleansumstats_metafile_kind: minimal
path_sumStats: study.snp_only.tsv.gz
delimiter: "\t"
missing_value: ""
stats_Model: logistic

col_SNP: rsid
col_EffectAllele: effect_allele
col_OtherAllele: other_allele
col_OR: odds_ratio
col_P: p_value
```

The lookup object should usually be a `.vmap` whose `id` column contains the raw summary-statistic IDs. Its metadata defines the starting build and contig naming before `prepare_variants.py` applies the requested `--dst-build` and `--dst-contig-naming`.

The usual source for an `--id-lookup` object is a `.vmap` that was already produced by `prepare_variants.py` from another dataset whose raw variant ID column contains dbSNP rsIDs or the same ID namespace used by the SNP-only summary statistics. Importers preserve those raw IDs in the `.vmap` `id` column, while coordinate and allele transforms update `chrom`, `pos`, `a1`, and `a2` as needed. `restrict_vmap.py` also preserves the source `.vmap` `id` values when filtering by target membership. That makes a prepared or restricted `.vmap` reusable as an rsID-to-coordinate lookup for later SNP-only summary-statistics imports.

The lookup object must also have its normal metadata sidecar. A `.vmap` sidecar lives at `<lookup>.vmap.meta.json` and stores the build and contig naming under `target`.

Warning: do not use `.vtable` files produced by `union_variants.py` or `intersect_variants.py` as rsID lookup inputs. Those tools intentionally rewrite the `id` column to target-derived `chrom:pos:a1:a2` IDs, so their IDs will not match SNP-only summary-statistic rsIDs. A `.vtable` is suitable for `--id-lookup` only if its `id` column is explicitly in the same ID namespace as the raw summary-statistic `SNP` values.

When a `.vmap` is supplied, only `chrom`, `pos`, and `id` are used. Its `source_shard`, `source_index`, `allele_op`, and source provenance are ignored, because this step is only recovering coordinates for raw summary-statistic IDs, not applying the lookup object's provenance.

Then prepare the SNP-only summary statistics:

```bash
prepare_variants.py \
  --input study.snp_only.tsv.gz \
  --input-format sumstats \
  --sumstats-metadata study.snp_only.meta.yaml \
  --id-lookup reference.vmap \
  --dst-build GRCh38 \
  --dst-contig-naming ncbi \
  --output work/study
```

In this mode, raw `SNP` values are matched against lookup IDs. Imported `chrom` and `pos` come from the lookup object, while imported `a1` and `a2` still come from raw `EffectAllele` and `OtherAllele`. The imported build and contig naming are inherited from the lookup object's metadata.

Rows with missing IDs or IDs absent from the lookup are dropped and recorded in `<output>.imported.vmap.qc.tsv` during preparation. The QC reasons are `invalid_id` and `id_not_found`. If a raw summary-statistic ID resolves to more than one distinct lookup `chrom` / `pos` pair, preparation fails clearly because the lookup object is ambiguous for a needed ID; duplicate lookup IDs that resolve to the same `chrom` / `pos` are allowed, regardless of lookup-object alleles.

`--id-lookup` is only a preparation-time coordinate-enrichment step. Downstream matching still uses exact prepared `chr:bp:a1:a2` rows and ignores variant IDs. It is independent of `--retain-snp-id`: `--retain-snp-id` is only a projection-time option for `project_payload.py` that decides whether the output `SNP` column uses `.vmap` IDs or generated `chrom:pos:a1:a2` IDs. It does not affect whether prepared `.vmap` files retain raw IDs for later use as `--id-lookup` objects.

For the opposite direction, use `assign_vmap_ids.py` after preparation or restriction. It matches already-prepared `.vmap` rows to an ID source by exact `chrom:pos:a1:a2` and copies IDs onto the `.vmap`; it does not fill missing coordinates by matching raw summary-statistic IDs. See [Assign vmap IDs](workflow.md#assign-vmap-ids) for policy details.

## Project Summary Statistics

Projection reconnects a prepared `.vmap` to the original raw summary-statistics file and writes rows in `.vmap` order. Use the same metadata that was used for preparation.

There are two projection modes:

- `--input-format sumstats` is the default projection mode. It writes the mapped target variant columns, keeps only payload/stat columns declared in the metadata, applies allele-orientation transforms where needed, and does not derive or clean statistical fields.
- `--input-format sumstats-clean` adds canonical summary-statistics harmonization. It normalizes recognized payload/stat columns to canonical names, validates numeric ranges, converts log-scaled p-values when metadata says they are log-scaled, drops all-missing payload columns, and derives missing fields where supported.

For the default projection mode:

```bash
project_payload.py \
  --input study.tsv.gz \
  --input-format sumstats \
  --sumstats-metadata study.meta.yaml \
  --vmap work/study.reference.vmap \
  --output out/study.reference.tsv.gz
```

The output always writes explicit canonical variant columns first:

```text
CHR  POS  SNP  EffectAllele  OtherAllele
```

By default, output `SNP` is generated from the `.vmap` row as `chrom:pos:a1:a2`. Use `--retain-snp-id` only when the output must preserve `.vmap` IDs instead.

The projected output also gets a metadata sidecar at `<output>.meta.yaml`. The sidecar describes the emitted tab-delimited file, target build and contig naming, canonical variant columns, retained payload columns, and whether clean harmonization was used.

## Clean Projection

For the clean projection mode:

```bash
project_payload.py \
  --input study.tsv.gz \
  --input-format sumstats-clean \
  --sumstats-metadata study.meta.yaml \
  --vmap work/study.reference.vmap \
  --output out/study.reference.clean.tsv.gz \
  --use-af-inference
```

Clean projection normalizes recognized payload/stat columns to canonical names, drops unrecognized payload columns, validates numeric ranges, converts `stats_neglog10P` or `stats_log10P` p-values to plain `P`, derives common fields where possible, and applies allele swap or flip/swap transforms to effect-size columns according to the `.vmap` provenance.

`--fill-mode` controls whether derived values are filled at column or row granularity:

- `--fill-mode column` is the default. A derivation creates a column only when the output column is absent.
- `--fill-mode row` fills missing cells in an existing output column when the required inputs are available for that row.

`--use-af-inference` enables additional allele-frequency-based derivations for `N`, `BETA`, and `SE` when the needed `Z`, `BETA`, `N`, and `EAF` inputs are available.

For logistic models, clean projection can derive `BETA` from `OR`, `SE` from odds-ratio confidence intervals, effective `N` from case/control counts, and combined `EAF` from case/control allele frequencies. For linear models, `stats_TotalN` can populate `N`. Other derivations include common relationships among `BETA`, `SE`, `Z`, and `P`.

The wrapper commands above call lower-level primitive tools internally. Users who need explicit control over individual stages can use `import_sumstats.py` and `apply_vmap_to_sumstats.py`; see [primitives.md](primitives.md) and the specs for those contracts.
