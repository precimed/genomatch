# Importers

## Tool rules

- `import_bim.py`, `import_pvar.py`, `import_vcf.py`, and `import_sumstats.py` emit `.vmap`, not `.vtable`.
- Import is the first provenance-aware mapping step in the canonical workflow.
- The importer-emitted `.vmap` target side contains retained canonical imported rows.
- The importer-emitted `.vmap` source side stores retained canonical target rows together with raw source provenance `(source_shard, source_index)`.
- Importers may infer `contig_naming`, but they must not normalize, repair, or rewrite input contig labels.
- Importers infer the first supported naming convention compatible with all retained rows in this fixed priority order: `ncbi`, `ucsc`, `plink`, `plink_splitx`.
- If importer contigs are mixed or invalid, the importer preserves retained rows as-is, omits `contig_naming` from metadata, emits a warning, and instructs the user to run `normalize_contigs.py`.
- All importers support `--chr2use` as an import-time row filter on retained target rows. `--contigs` is an alias.
- All importers support optional `--max-allele-length` (positive integer), defaulting to `150`.
- `import_bim.py`, `import_pvar.py`, and `import_vcf.py` may accept `@` in the raw input path as a shard-discovery convention.
- `import_sumstats.py` does not accept `@`; summary statistics are expected to come from a single file.

## Import-time allele ordering

Importer-emitted `.vmap` rows preserve source-format allele meaning at import time. They do not reorder alleles to the reference-second canonical target form used by downstream reference-aware tools.

Format-specific rules:

- `import_bim.py` writes PLINK `.bim` column 5 to `a1` and column 6 to `a2`
- `import_pvar.py` writes `ALT` to `a1` and `REF` to `a2`
- `import_vcf.py` writes `ALT` to `a1` and `REF` to `a2`
- `import_sumstats.py` writes `EffectAllele` to `a1` and `OtherAllele` to `a2`

Reference-aware canonicalization belongs to downstream tools such as `restrict_build_compatible.py` and `liftover_build.py`, which emit target rows in reference-second order with `a1=non-reference` and `a2=reference`.

## Retention, not repair

Import canonicalizes by retention, not by repair:

- rows representable in canonical target-side form may be retained
- rows not representable in canonical target-side form may be dropped
- importers must not silently mutate malformed or noncanonical source alleles into canonical ones

Examples of rows that may be dropped at import time include:

- `multiallelic`
- `non_actg_allele`
- `allele_too_long`
- `malformed_row`
- `filtered_by_chr2use`

Import-time allele length cap contract:

- apply the cap only at importer raw-input boundaries, using `--max-allele-length` (default `150`)
- if either imported allele exceeds the active cap, drop the source row with auditable QC reason `allele_too_long`
- this cap is importer-local and must not be re-validated by downstream transform tools, including `restrict_build_compatible.py` and `liftover_build.py` output parsing
- this cap is not re-validated after downstream allele normalization steps
- this cap is not re-validated when loading user-provided `.vtable` / `.vmap` objects

The preserved contract is:

- every retained target row has exact source provenance
- every dropped source row is auditable in import QC

## Source provenance

For importer-emitted `.vmap`:

- `source_shard` is the exact discovered shard label that replaced `@`
- for single-file imports, `source_shard` is `.`
- `source_index` is the zero-based row index within that `source_shard`
- `source_index` counts raw variant data rows only, not file preamble or schema lines
- for `.bim`, indexing starts at the first non-empty data row
- for `.pvar`, `##` meta lines and the single `#` header line are skipped and do not count toward `source_index`
- for `.vcf`, `##` meta lines and the `#CHROM` header line are skipped and do not count toward `source_index`

`source_shard` is a locator string, not a biology field. Importers must store it verbatim and must not interpret it as chromosome identity.

## `@` shard discovery for raw importer inputs

When `@` is present in the input path to `import_bim.py`, `import_pvar.py`, or `import_vcf.py`, the tool treats the input as a set of raw shards and applies the shared filename discovery contract from [shard-discovery.md](shard-discovery.md).

Import order rules:

- concatenate discovered shards in deterministic lexical path order
- preserve row order within each shard
- store the exact discovered replacement token as that shard's `source_shard`
- assign `source_index` from the row order within each discovered shard

For both single-file and sharded importer inputs, `--chr2use` / `--contigs` filters retained target rows by canonical contig while preserving retained row order. Rows dropped by that filter must be auditable in import QC.

## Import QC

If an importer drops any source rows, it must emit a QC sidecar at:

- `<output>.qc.tsv`

QC rows are keyed by original source provenance:

- `source_shard`
- `source_index`

and must include:

- `reason`

Additional raw source fields or identifiers may be included when useful, but the provenance fields and reason code are required.

## Importers

- `import_bim.py` writes `.vmap` and metadata
- `import_pvar.py` writes `.vmap` and metadata
- `import_vcf.py` writes `.vmap` and metadata
- `import_sumstats.py` writes `.vmap` and metadata
- importers preserve input contig labels exactly; they do not normalize or repair rows
- importers omit `contig_naming` and emit a warning when retained rows are mixed or invalid
- importers infer the first compatible supported naming in priority order `ncbi`, `ucsc`, `plink`, `plink_splitx`
- `import_sumstats.py` honors the cleansumstats metadata contract
- for `import_sumstats.py`, `--input` may be omitted; in that case, use `path_sumStats` from `--sumstats-metadata`
- when `import_sumstats.py --input` is omitted, resolve `path_sumStats` as `<directory of --sumstats-metadata>/<path_sumStats>`
- `path_sumStats` is a filename-only field (no `/`) and therefore always resolves within the metadata directory

The summary-stat metadata contract is `match/schemas/raw-sumstats-metadata.yaml`.

### `import_sumstats.py`: optional ID-based coordinate enrichment

- `import_sumstats.py` normally imports target-side variant rows from summary-stat columns described by the metadata contract.
- `import_sumstats.py` may accept optional `--id-vtable <path>` as a fallback mode for recovering missing target-side `chrom` / `pos` information from an auxiliary `.vtable`.
- `--id-vtable` is allowed only for `import_sumstats.py`; it is not a general matching mode for the toolkit.
- `--id-vtable` must name a `.vtable` whose `id` column is used as a lookup key.
- `--id-vtable` mode is valid only when the summary-stat metadata defines `SNP`, `EffectAllele`, and `OtherAllele`, and does not define `CHR` or `POS`
- in `--id-vtable` mode, `import_sumstats.py` must look up each raw summary-stat `SNP` value against `.vtable.id`
- in `--id-vtable` mode, imported target-side `chrom` and `pos` come from the matched `.vtable` row rather than from the raw summary-stat payload
- in `--id-vtable` mode, imported metadata `genome_build` and `contig_naming` are inherited from the matched `.vtable` object metadata rather than inferred from the raw summary-stat payload
- `a1` and `a2` still come from the raw summary-stat payload as `EffectAllele` and `OtherAllele`; `--id-vtable` does not change import-time allele ordering
- rows of `--id-vtable` whose `id` is missing, empty, or `.` are ignored by the lookup table and must emit a warning
- if a raw summary-stat `SNP` value is missing, empty, or `.`, the source row must be dropped with auditable import QC reason `invalid_id`
- if a raw summary-stat `SNP` value has no match in the filtered `.vtable.id` lookup set, the source row must be dropped with auditable import QC reason `id_not_found`
- if a raw summary-stat `SNP` value matches multiple rows in the filtered `.vtable.id` lookup set, the source row must be dropped with auditable import QC reason `ambiguous_id_match`
- if `--id-vtable` is not supplied and the importer can not determine target-side `chrom` and `pos` from the raw summary-stat input, the importer must fail clearly
- `--id-vtable` is an import-time coordinate-enrichment step only; downstream exact matching semantics remain `chr:bp:a1:a2`, ignoring `id`
