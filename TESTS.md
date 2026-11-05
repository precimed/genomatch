# TESTS.md

This file is the authoritative migration-aware test specification for `match/`.

## Test environment

- Canonical test execution for this toolkit should run under the `match-liftover` conda environment.
- Use that environment for both `pytest tests` and full `pytest` validation when reviewing `match/` changes.
- This keeps the real liftover integration coverage active instead of relying on skipped-test behavior in a lighter environment.
- `MATCH_REFERENCE_ACCESS_MODE` defaults to `BULK`, matching production behavior. End-to-end integration tests that rely on full real FASTA sequence assets must set `MATCH_REFERENCE_ACCESS_MODE=LEGACY`; preloading whole production FASTAs is multiple-GB behavior and is not appropriate for CI-sized tests.

Active migrated modules currently include:

- `tests/test_contig_cleanup_utils.py`
- `tests/test_importers.py`
- `tests/test_normalize_contigs.py`
- `tests/test_restrict_contigs.py`
- `tests/test_drop_strand_ambiguous.py`
- `tests/test_guess_build.py`
- `tests/test_restrict_build_compatible.py`
- `tests/test_liftover_build.py`
- `tests/test_sort_variants.py`
- `tests/test_variant_rendering.py`
- `tests/test_schema_metadata.py`
- `tests/test_convert_vmap_to_target.py`
- `tests/test_intersect_variants.py`
- `tests/test_restrict_vmap.py`
- `tests/test_assign_vmap_ids.py`
- `tests/test_union_variants.py`
- `tests/test_apply_vmap_to_sumstats.py`
- `tests/test_apply_vmap_to_bfile.py`
- `tests/test_apply_vmap_to_pfile.py`
- `tests/integration/test_apply_vmap_to_pfile_roundtrip.py`
- `tests/test_prepare_variants.py`
- `tests/test_prepare_variants_sharded.py`
- `tests/test_project_payload.py`
- `tests/test_plink_vcf_helpers.py`
- `tests/integration/test_liftover_build_bcftools_real.py`
- `tests/integration/test_end_to_end_sumstats_bcftools_real.py`
- `tests/integration/test_apply_vmap_to_bfile_vcf_roundtrip.py`
- `tests/integration/test_apply_vmap_to_bfile_haploid_roundtrip.py`
- `tests/integration/test_apply_vmap_to_sumstats_regressions.py`

Historical legacy coverage is preserved in the migration changelog under `match/changelog/TEST_2026_03_10.md`, but the active canonical suite is now `.vtable` / `.vmap` only.

## Preserved behavior

The migrated suite must keep coverage for:

- exact source-agnostic intersection with declared-coordinate output
- allele swap handling
- mapped-only `.vmap` rows with non-missing source provenance
- importer first-occurrence cleanup, explicit sort/dedup cleanup, and handling distinct input variant keys that become equal after transformation
- rsID ignored for exact set operations
- importer `--chr2use` / `--contigs` filtering and `restrict_contigs.py` target-row filtering
- summary-stat effect inversion for swapped alleles
- target-derived IDs for provenance-free set outputs

## Schema compatibility fixtures

Keep small checked-in fixtures with original metadata sidecars for every supported
[object/schema pair](spec/core-objects.md#schema-versions). Store them under
`tests/fixtures/schema/` with flat, version-labelled names such as
`vmap_1_0.vmap` and `vtable_1_0.vtable`.
Fixtures must include a multibase row; the `.vmap` fixture must also include a
non-`identity` operation and a source shard token other than `.`.

- Compatibility tests must load these committed files through the normal object
  readers and assert expected logical contents: row order, IDs, ordered alleles,
  metadata, and `.vmap` provenance and `allele_op`. Successful loading alone is
  insufficient; documented migration defaults, if any, must also be checked.
- Each fixture must have a schema-specific generator or reproduction instructions
  beside it recording the producing release/commit. Historical fixtures must come
  from the corresponding writer, not a current-writer output with its version
  label changed. Do not regenerate them during tests or routine fixture updates.
- Add fixtures when introducing a schema; retain existing fixtures for every
  supported older schema. Any intentional support removal must be documented
  and tested as a clear rejection rather than silently changing expectations.
- Keep current-writer save/load round trips separate from historical compatibility
  tests. Using temporary fixture copies, verify that missing, malformed, and
  unsupported versions are rejected clearly for both `.vmap` and `.vtable`.
- Exercise rejection through every CLI path that consumes canonical objects,
  including `--id-lookup` and wrapper metadata reads. Use malformed variant rows
  to prove the version check runs before row parsing, and assert the object/version
  diagnostic and regeneration remedy. Cover missing, numeric, null, empty,
  whitespace-altered, and unsupported string versions.
- Verify every canonical-object producer emits its explicit supported version,
  including transforms, conversions, and sharded wrapper outputs. `write_metadata`
  must set the current writer version without mutating the supplied metadata.
  Metadata-only updates must validate both the existing and replacement metadata
  and reject version changes; `guess_build.py --write` must preserve the validated
  input version. Use `write_object_metadata` in `tests/utils.py` for ordinary test
  sidecars, raw JSON for missing/invalid-version tests, and original fixture bytes
  for historical compatibility tests. A supported
  version without `variants_count` must remain readable under the existing rules.
- Raw sharded concatenation must reject a supported non-current schema on any
  input, including empty inputs and inputs without `variants_count`, before
  overwriting output. Current-schema inputs must preserve row bytes and order.
- For both preparation wrappers, cover an invalid final object under `--resume`
  when no stage would otherwise run, invalid skipped intermediates or per-group
  finals, and normal reuse of supported objects. Schema-invalid outputs must cause
  failure without overwriting them or returning success.

## Required migration tests

### Importers

- `import_bim.py`, `import_pvar.py`, `import_vcf.py`, and `import_sumstats.py` write `.vmap` and metadata
- importer-emitted `.vmap` target rows are the retained canonical imported rows
- importer exact-identity duplicates keep the first otherwise-retained occurrence in resolved shard and row order, ignoring `id`; subsequent duplicates have original source locators and reason `duplicate_target` in import QC
- importer-emitted `.vmap` source provenance is exact `(source_shard, source_index)`
- importer-emitted `.vmap` uses `source_shard=.` for single-file imports
- `import_bim.py`, `import_pvar.py`, and `import_vcf.py` support raw-input `@` discovery
- for `@` imports, `source_shard` stores the exact discovered filename substitution token verbatim
- raw-input `@` discovery accepts PLINK split-X PAR shard filename aliases `XY` and `chrXY`
- importers support `--chr2use` / `--contigs` as an import-time target-row filter
- importers preserve input contig labels exactly on retained rows; they do not normalize or repair rows
- importers may drop unsupported source rows instead of failing the whole import
- importers infer the first compatible supported naming in priority order `ncbi`, `ucsc`, `plink`, `plink_splitx`
- importers emit a warning only when retained contigs are mixed or invalid
- dropped importer rows must be auditable in `<output>.qc.tsv`
- importer QC rows are keyed by `source_shard` and `source_index` and include a reason code
- importer QC must cover at least representative dropped-row reasons such as `multiallelic`, `non_actg_allele`, `malformed_row`, and `filtered_by_chr2use`
- `import_sumstats.py` honors the cleansumstats metadata contract
- `import_sumstats.py` rejects `@` input paths
- `import_sumstats.py` supports optional `--max-allele-length` with default `150`
- in the summary-stats import path, source rows whose `EffectAllele` or `OtherAllele` exceeds the active cap are dropped with QC reason `allele_too_long`
- `import_sumstats.py --id-lookup` is valid only when summary-stat metadata defines `SNP`, `EffectAllele`, and `OtherAllele`, and does not define `CHR` or `POS`
- `import_sumstats.py --id-lookup` fills imported target-side `chrom` and `pos` from matched `.vtable` or `.vmap` target-side `id` rows while keeping imported `a1/a2` from raw `EffectAllele` / `OtherAllele`
- `import_sumstats.py --id-lookup` inherits imported `genome_build` and `contig_naming` from the lookup object's target-side metadata
- `import_sumstats.py --id-lookup` accepts `.vmap` lookup objects while ignoring source provenance
- `import_sumstats.py --id-lookup` ignores lookup-table rows whose `id` is missing, empty, or `.`, and emits a warning
- `import_sumstats.py --id-lookup` allows duplicate lookup IDs when they resolve to the same `chrom` and `pos`, ignoring lookup-object alleles
- `import_sumstats.py --id-lookup` fails clearly when a raw summary-stat `SNP` value resolves to more than one distinct lookup `chrom` / `pos` pair
- `import_sumstats.py --id-lookup` drops source rows whose `SNP` is missing, empty, or `.` with QC reason `invalid_id`
- `import_sumstats.py --id-lookup` drops unmatched source rows with QC reason `id_not_found`

### Metadata and contigs

- metadata sidecars are required by transformation tools
- metadata is authoritative once present; valid target rows under declared `contig_naming` are an object invariant, but downstream tools are not required to re-validate every row-level contig label at every boundary
- `guess_build.py` accepts `.vtable` and `.vmap`
- for `.vmap` input, `guess_build.py` updates target-side metadata only and does not change source provenance
- `guess_build.py` requires declared `contig_naming`
- `guess_build.py` writes inferred metadata only when `genome_build` is `unknown` or `--force` is used
- `guess_build.py` guesses build from configured reference assets and reports evidence
- `guess_build.py` defaults to random downsampling of at most 10,000 target rows when computing build evidence
- `guess_build.py --sample-rows 0` disables downsampling and evaluates all target rows
- `guess_build.py` rejects negative `--sample-rows`
- `guess_build.py` accepts UCSC-named input without extra normalization
- `guess_build.py` returns `genome_build=unknown` when evidence is ambiguous, below threshold, or unusable
- malformed metadata sidecars fail clearly during `guess_build.py`
- `guess_build.py` explicitly covers declared `plink` and `plink_splitx` inputs in addition to `ncbi` and `ucsc`
- `restrict_build_compatible.py` explicitly covers declared `plink` and `plink_splitx` inputs in addition to `ncbi` and `ucsc`, except that `--norm-indels` rejects `plink_splitx`
- `liftover_build.py` explicitly covers declared `plink` input in addition to `ncbi` and `ucsc`, and rejects `plink_splitx` clearly
- `normalize_contigs.py` supports `ncbi <-> ucsc`
- `normalize_contigs.py` supports `plink -> ncbi`
- `normalize_contigs.py` supports `plink_splitx -> ncbi`
- `normalize_contigs.py` may repair mixed or inconsistent target rows and must preserve row order
- `normalize_contigs.py` must not rewrite `source_shard`
- `normalize_contigs.py` drops unresolved target rows and reports normalized vs unresolved-drop counts
- `normalize_contigs.py` drops all rows whose distinct input variant keys become equal after contig normalization and audits them with status `duplicate_target` under the row-level QC rules
- `restrict_contigs.py` filters retained target rows and preserves retained target order
- reference-aware tools must accept UCSC-named input rows
- reference-aware tools must accept NCBI-style (`ncbi`) input rows
- reference-aware tools must accept PLINK-style (`plink`) input rows
- `guess_build.py` must accept PLINK split-X (`plink_splitx`) input rows
- `restrict_build_compatible.py` must accept PLINK split-X (`plink_splitx`) input rows when `--norm-indels` is off, and reject them clearly when `--norm-indels` is on
- `liftover_build.py` must reject PLINK split-X (`plink_splitx`) input rows clearly
- tests for `guess_build.py` and `restrict_build_compatible.py` must verify internal normalization to UCSC reference naming for each supported non-UCSC declared input naming mode (`ncbi`, `plink`, and `plink_splitx`) within their supported modes
- tests for `liftover_build.py` must verify internal normalization to UCSC reference naming for supported non-UCSC declared input naming modes (`ncbi` and `plink`)

### Reference access

- Toolkit FASTA lookup consumers (`guess_build.py` and `restrict_build_compatible.py`) must honor `MATCH_REFERENCE_ACCESS_MODE` with accepted values `BULK` and `LEGACY`, parsed case-insensitively
- when `MATCH_REFERENCE_ACCESS_MODE` is unset, these lookup consumers must run in default `BULK` mode
- for identical inputs and reference assets, `BULK` and `LEGACY` modes must produce consistent outputs and equivalent QC/error semantics for `guess_build.py` and `restrict_build_compatible.py`
- cover single-base and variable-length span queries under the [reference access policy](spec/variant-transforms.md#reference-access-policy), including distinct lengths at the same contig and position, first-base-only matches, contig boundaries, missing contigs, reference characters outside `A/C/G/T`, lowercase FASTA sequence, and empty query batches
- verify `BULK` avoids repeated per-variant indexed FASTA fetches and keeps results isolated by reference FASTA path, without prescribing cache layout, helper calls, or vectorization; preserve existing single-base lookup and build-guess behavior

### Tools accepting `.vmap` or `.vtable`

- for every primitive listed in [tools accepting `.vmap` or `.vtable`](spec/variant-transforms.md#tools-accepting-vmap-or-vtable), run equivalent ordered canonical target rows and metadata as `.vmap` and `.vtable` with the same options and reference assets; compare target rows, IDs, order, target metadata, retained/dropped counts, exclusion statuses, warnings, and failure behavior; also compare documented cleanup modes
- under the [duplicate-input rule](spec/variant-transforms.md#duplicate-containing-input), verify the five rejecting primitives and the specified source-mapping readers fail on repeated variant keys and name the cleanup route; cover both object types where supported; verify explicit keep-first sorting cleanup, contig repair accepting exact repeats and dropping all duplicate output keys to original-locator QC (including mixed labels without `contig_naming` and unchanged labels), and build guessing preserving repeated-row weighting and total row counts without rewriting input rows; set-membership and ID-lookup inputs keep their own repeated-key rules
- verify output type and metadata layout follow the input type; `.vmap` preserves original provenance and composes allele operations, while `.vtable` contains no provenance or `allele_op` columns
- exercise missing and malformed target metadata, reference mismatches, all-filtered input, duplicate cleanup, and distinct input variant keys becoming equal after transformation for both types within each primitive's supported behavior
- verify [variant keys](spec/core-objects.md#variant-keys) use the exact stored values before and after each transform and ignore IDs and provenance; include distinct contig spellings and allele orders that become equal only after the relevant transform
- for operations requiring row-level QC, check the same sidecar path and statuses for both types using the [shared QC row identifiers](spec/variant-transforms.md#row-level-qc); include filtering, sorting, and external-tool round trips to detect accidental reindexing of QC row identifiers
- verify the type-specific QC column layouts in that contract, keeping `.vtable` row ordinals distinct from `.vmap` provenance throughout processing; reject `.vtable` QC frames carrying provenance fields
- on successful runs with no rows dropped, including supported empty input, verify that no QC sidecar is emitted and any stale sidecar is removed for both types
- verify all affected rows are dropped when distinct input variant keys become equal after contig normalization, reference canonicalization, or liftover for both types, including normalized rows matching unchanged rows; separately verify keep-first cleanup in `sort_variants.py --drop-duplicates`
- verify explicit no-sidecar policies, including `sort_variants.py --drop-duplicates`, apply equally to both types
- preserve rejection of `.vtable` as the source mapping for `restrict_vmap.py`, `assign_vmap_ids.py`, and payload-application tools

### Mapping

- `intersect_variants.py` always emits provenance-free `.vtable`
- `intersect_variants.py` emits target-derived IDs in declared coordinate order
- `intersect_variants.py` supports both SNP and non-SNP biallelic exact intersection
- `intersect_variants.py` accepts two or more positional `.vmap` / `.vtable` inputs and named `--output`
- `restrict_vmap.py` emits `.vmap` preserving source `.vmap` row order, IDs, provenance, and existing `allele_op`
- `restrict_vmap.py` accepts positional inputs with source `.vmap` first, then one or more `.vmap` / `.vtable` restrictions, plus named `--output`
- `restrict_vmap.py` treats restriction inputs as exact `(chrom, pos, a1, a2)` membership filters only
- `restrict_vmap.py` accepts one or more `.vtable` / `.vmap` restriction inputs
- with multiple restriction inputs, `restrict_vmap.py` retains only source rows present in every restriction input
- `restrict_vmap.py` ignores restriction-input `.vmap` provenance
- `restrict_vmap.py` supports both SNP and non-SNP biallelic exact restriction
- for multi-base alleles, any `flip` / `flip_swap` semantics exercised by reference-aware tools use reverse-complement per allele string rather than basewise complementation without reversing order
- `intersect_variants.py` rejects `@` template paths and has no `--chr2use` surface
- `intersect_variants.py` fails cleanly when required metadata omits `contig_naming`; row-level contig-label validation is not part of its exact-intersection contract
- `restrict_vmap.py` rejects `@` template paths and has no `--chr2use` surface
- `restrict_vmap.py` fails cleanly when required metadata omits `contig_naming`; row-level contig-label validation is not part of its exact-restriction contract
- repeated membership keys produce the same intersection, union, and restriction-filter results as unique keys for either object type, including repetitions within and across chunk boundaries
- `drop_strand_ambiguous.py` drops target-side strand-ambiguous biallelic rows, defined by `reverse_complement(a1) == a2`
- `drop_strand_ambiguous.py` applies that rule to both SNP and non-SNP biallelic alleles
- `sort_variants.py` sorts target rows into declared coordinate order
- exact matching tools (`intersect_variants.py`, `union_variants.py`, `restrict_vmap.py`, `assign_vmap_ids.py`) reject any `genome_build: unknown` before reading rows or emitting output, including all-unknown inputs and an unknown input in either role, for both `.vmap` and `.vtable` membership inputs
- declared coordinate order follows the [row-order contract](spec/core-objects.md#row-order-contract)
- exercise the allele tie-break with a co-located SNV and deletion whose `a1`-first and `a2`-first orders differ, covering `sort_variants.py`, `intersect_variants.py`, `union_variants.py`, and `liftover_build.py`; verify numeric position and declared contig ordering independently of rendered ID strings
- generated-ID coverage must include every generation surface in the [generated-ID contract](spec/core-objects.md#generated-variant-ids), multibase alleles, and unverified rows whose IDs remain slot-derived without reference validation; copied and explicitly retained IDs remain unchanged
- `.vmap` rows carry `source_shard + source_index` provenance
- `.vmap` preserves original `source_shard` and shard-local `source_index`
- `.vmap` rows must have non-missing source provenance and non-negative `source_index`
- `.vmap` metadata carries target-side metadata only
- `convert_vmap_to_target.py` materializes target rows and breaks provenance intentionally
- `liftover_build.py` resolves chain and FASTA assets from config + metadata only
- `liftover_build.py` accepts declared `ncbi`, `ucsc`, and `plink` input contig naming while using UCSC internal reference assets
- `liftover_build.py` rejects declared `plink_splitx` input clearly
- `liftover_build.py` supports `--resume` for reparsing retained bcftools liftover output without rerunning bcftools
- `liftover_build.py` fixtures must cover prepared deletions, MNVs, and delins with multibase source `a2`, and inspect source VCF records for `REF=a2, ALT=a1` under the [required source reference compatibility](spec/variant-transforms.md#required-source-reference-compatibility)
- with a stubbed liftover plugin, verify structurally valid source rows are forwarded without source-FASTA sequence checks or source allele reselection, including `a1`-only and partial-REF matches; do not require biological correctness or mismatch detection for inputs violating the source reference compatibility requirement
- `liftover_build.py` audits unrepresentable lifted outputs under the [liftover acceptance contract](spec/variant-transforms.md#liftover), covering multiallelic, missing/symbolic/non-sequence, identical-allele, invalid-position, and multiple-result cases, plus unsupported allele-operation annotations such as `SWAP=-1` and `SWAP=2`; row-level exclusions do not fail the invocation
- `liftover_build.py` drops lifted rows on unsupported target contigs with QC status `unsupported_target_contig`
- `liftover_build.py` must drop rows with QC status `ploidy_class_changed` when lifted `chrom` or `pos` changes the expected `(male_ploidy, female_ploidy)` pair under the shared ploidy model
- for `ploidy_class_changed`, the comparison is against the final candidate row that would otherwise be emitted, immediately before final output/QC accounting
- `liftover_build.py` drops all lifted rows whose distinct input variant keys become equal on target `(chrom, pos, a1, a2)`, auditing them with status `duplicate_target`
- `liftover_build.py` audits non-retained input rows under the shared row-level QC contract
- `liftover_build.py` accepts UCSC-named input without changing declared output naming
- `liftover_build.py` fails clearly on missing chain config, missing bcftools, and missing FASTA indexes; trusting source reference compatibility does not waive schema, metadata, contig, or reference-asset validation
- `liftover_build.py` fails clearly on unsupported input contig labels, external-tool failure, unparseable output, and output that cannot be associated unambiguously with input rows
- cover lifted records with an unsupported allele-operation annotation and an independent contig or ploidy exclusion; require an applicable QC status without invocation failure, without pinning an internal check order
- verify lifted output uses the plugin's representation and annotation-derived allele operation without independent allele reassignment, normalization, or strand-ambiguity filtering
- `liftover_build.py` covers the degenerate all-dropped unrepresentable-output case with QC output
- tests for `liftover_build.py` must cover lifted one-side-multibase rows being retained as canonical `a1=non-reference`, `a2=reference`
- use a small real-plugin integration fixture set for both-multibase output and a reverse-strand multibase mapping, checking source VCF `REF=a2, ALT=a1` and composed allele operations
- `restrict_build_compatible.py` resolves reference FASTA from config + metadata only
- `restrict_build_compatible.py` accepts declared `ncbi`, `ucsc`, `plink`, and `plink_splitx` input contig naming while using UCSC internal reference FASTA
- `restrict_build_compatible.py --norm-indels` rejects declared `plink_splitx` input clearly
- regression coverage for [same-build reference compatibility](spec/variant-transforms.md#same-build-reference-compatibility) must include both-allele matches, `a1`-only matches, neither-allele matches, and partial-REF matches where only the first base agrees
- verify the specified exclusion status and original-provenance QC for no permitted reference match, including a match requiring a disabled operation
- cover all three input classes over inclusion enabled/disabled, primitive option defaults, and every combination of swap, strand-flip, and normalization permission in that contract
- cover primitive input exclusions and counts, checking the specified QC path/status and shared QC row identifiers
- cover identical single-base and multibase inputs, including mixed-case equivalents, with inclusion and normalization independently enabled/disabled; verify the specified identical-input status before length filtering or reference-allele assignment, using the shared row-level QC contract
- check warnings for either or both multibase-operation flags without input inclusion, including class 1-only input; the flags must not bypass exclusion or suppress class 1 reference-allele assignment
- include `a1=A, a2=AG` on reference `AGC`: the deletion must remain unchanged with normalization enabled and disabled
- include the same allele pair on `ATC` to exercise a mismatching supplied REF with and without swap permission
- cover operation composition with an existing `swap`, and effect/genotype projection, verifying that reference-allele assignment and normalization do not introduce an extra payload recode
- strand fixtures must exercise class 1 `flip` and `flip_swap`, and verify that enabling strand flips does not introduce either operation for classes 2 and 3, including equal-length multibase alleles
- include `a1=A, a2=TT` on reference `AAG` with strand flips enabled: without indel-swap permission it is dropped; with permission the forward swap is selected
- `restrict_build_compatible.py` keeps only same-build reference-compatible rows
- `restrict_build_compatible.py` canonicalizes retained rows to reference-second target ordering with `a1=non-reference` and `a2=reference`
- a tiny-FASTA PVAR import/restriction fixture with one reversed row among normal rows must yield `swap` only for that row and `identity` for the rest, preserving row order and original provenance
- `restrict_build_compatible.py` rejects removed flags `--sort` and `--drop-duplicates`
- cover class 3 inputs for which both REF/ALT assignments would pass REF validation: valid `a2` must be selected rather than dropping the row as ambiguous
- include padded SNVs, padded deletions, MNVs, and delins in class 3 fixtures, with input inclusion enabled and normalization enabled/disabled
- cover output validation failures of the selected REF/ALT allele assignment with an alternative assignment that would otherwise pass, to detect fallback after choosing an allele assignment
- verify that normalization receives only the selected REF/ALT allele assignment and that class 1 input bypasses `bcftools norm`
- with `--norm-indels`, rows whose normalized `chrom` or `pos` changes the expected `(male_ploidy, female_ploidy)` pair under the shared ploidy model are dropped with QC status `ploidy_class_changed`
- with `--norm-indels`, verify retained diagnostic VCFs and logs in the [normalization diagnostic namespace](spec/variant-transforms.md#normalization-diagnostics), including record coverage across actual calls or batches and selected `REF=a2` / `ALT=a1` assignments; require the namespace but not fixed suffixes, a fixed file count, or aggregation
- prepopulate `<output>.bcftools_norm*` with stale diagnostics, including legacy `_e`/`_x` files and higher-numbered batch files; verify complete cleanup on every restriction run with normalization enabled or disabled, including reruns producing fewer batches, while files outside the namespace remain untouched
- when no normalization runs, including when `--norm-indels` is omitted or no records require normalization, verify that no normalization subprocess runs and `<output>.bcftools_norm*` is empty after the restriction run
- with `--norm-indels`, `restrict_build_compatible.py` uses `bcftools norm` only for normalization against the configured same-build reference FASTA
- with `--norm-indels`, `restrict_build_compatible.py` does not use `bcftools norm` multiallelic split or join modes
- with `--norm-indels`, `restrict_build_compatible.py` does not use `bcftools norm` atomization
- with `--norm-indels`, `restrict_build_compatible.py` invokes `bcftools norm` with `-c e` and does not use `-c x`, `-c w`, or `-c s`
- with `--norm-indels`, normalized representation changes introduce no additional allele operation; reference-allele assignment for multibase input contributes only local `identity` or `swap`, composed with the existing `allele_op`
- cover every failure in the [normalization output-validation table](spec/variant-transforms.md#normalization), verifying exclusion with the specified status and shared QC row identifiers; retained output must remain one canonical biallelic row per retained input row
- cover an input submitted for normalization that yields no output record alongside a successfully normalized row; verify `norm_ref_mismatch` QC, retention of the successful row, and successful invocation; also cover all records submitted for normalization missing from the output
- verify malformed normalization output and output records with unknown input identifiers fail the invocation
- with `--norm-indels`, inject unexpected REF mismatches at the normalization boundary and in normalized output to verify fatal invariant failures after successful reference-allele assignment
- with `--norm-indels`, duplicate-target detection still applies after normalization on final canonical `(chrom, pos, a1, a2)`
- with `--norm-indels`, tests must cover indel left-normalization or trimming that changes target `pos` and allele strings while preserving canonical final `a1=non-reference`, `a2=reference`
- with `--norm-indels`, tests must cover the case where no normalization change is needed and the retained output remains otherwise identical to the non-`--norm-indels` path
- `restrict_build_compatible.py` drops all rows whose distinct input variant keys become equal after canonicalization and audits them with status `duplicate_target` under the row-level QC rules
- include a padded SNV that normalizes to the same target as an input SNV, verifying that both rows are dropped with `duplicate_target`
- for restriction and liftover, cover an excluded result sharing an output variant key with an otherwise valid candidate; the excluded result must not cause the valid candidate to be dropped as `duplicate_target`
- `restrict_build_compatible.py` emits JSON telemetry to stdout with input/output paths, declared metadata, normalization path, and retained/dropped row counts
- `restrict_build_compatible.py` accepts UCSC-named input without extra normalization
- `restrict_build_compatible.py` rejects unknown target contigs and rows inconsistent with declared `contig_naming`
- `restrict_build_compatible.py` fails clearly on missing reference config sections and missing FASTA indexes
- `restrict_build_compatible.py` fails clearly on bad relative reference paths and on explicit non-ACGT alleles
- restriction output preserves retained input order
- with `--norm-indels`, tests should verify clear failure on missing `bcftools`, failed `bcftools norm`, and missing FASTA indexes

### Intersections

- `intersect_variants.py` intersects exact `(chrom, pos, a1, a2)`
- `intersect_variants.py` accepts positional inputs and named `--output`
- `intersect_variants.py` requires the same `genome_build` and the same `contig_naming` across inputs
- `intersect_variants.py` fails cleanly when input metadata omits `contig_naming`
- `intersect_variants.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- output IDs follow the [generated-ID contract](spec/core-objects.md#generated-variant-ids)
- output order is declared coordinate order using the same ordering contract as `sort_variants.py`
- `intersect_variants.py` emits `.vtable`

### `.vmap` restrictions

- `restrict_vmap.py` restricts exact `(chrom, pos, a1, a2)`
- `restrict_vmap.py` accepts positional inputs and named `--output`
- `restrict_vmap.py` requires the same `genome_build` and the same `contig_naming` across inputs
- `restrict_vmap.py` fails cleanly when input metadata omits `contig_naming`
- `restrict_vmap.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- output IDs, order, provenance, and `allele_op` come from the source `.vmap`
- restriction inputs affect membership only
- multiple restriction inputs are combined by intersection, not union
- `restrict_vmap.py` emits `.vmap`

### `.vmap` ID assignment

- `assign_vmap_ids.py` accepts one source `.vmap`, one `.vmap` / `.vtable` ID source, and named `--output`
- `assign_vmap_ids.py` requires the same `genome_build` and the same `contig_naming` across inputs
- `assign_vmap_ids.py` performs no implicit normalization
- matching is by exact `(chrom, pos, a1, a2)`
- duplicate `(chrom, pos, a1, a2)` keys in the ID source are ignored when absent from the source `.vmap`
- duplicate `(chrom, pos, a1, a2)` keys in the ID source fail clearly when present in the source `.vmap`
- output IDs come from the ID source
- output row order, coordinates, provenance, and `allele_op` come from the source `.vmap`
- by default, rows without an ID-source match are dropped and audited in `<output>.qc.tsv` with `id_not_found`
- `--unmatched-id-policy variant-key` retains rows without an ID-source match with IDs following the [generated-ID contract](spec/core-objects.md#generated-variant-ids) and no QC row
- `--unmatched-id-policy missing` retains rows without an ID-source match with `id=.` and no QC row
- rows matched to a missing ID-source ID are dropped and audited with `missing_id`
- by default, duplicate retained non-missing output IDs fail clearly
- `--duplicate-id-policy allow` retains duplicate non-missing output IDs
- `--duplicate-id-policy drop-all` drops all rows carrying duplicate non-missing output IDs and audits them with `duplicate_id`
- duplicate-ID checks ignore missing IDs (`.` or empty)
- `assign_vmap_ids.py` emits `.vmap`

### Unions

- `union_variants.py` unions exact `(chrom, pos, a1, a2)`
- `union_variants.py` accepts positional inputs and named `--output`
- `union_variants.py` requires at least two inputs
- `union_variants.py` requires the same `genome_build` and the same `contig_naming` across inputs
- `union_variants.py` fails cleanly when input metadata omits `contig_naming`
- `union_variants.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- duplicate exact rows are deduplicated by exact target key, with output IDs following the [generated-ID contract](spec/core-objects.md#generated-variant-ids)
- `union_variants.py` re-sorts the deduplicated output into declared coordinate order using the same ordering contract as `sort_variants.py`
- ties under declared coordinate order preserve first-occurrence order
- `union_variants.py` emits `.vtable`

### Payload application

- `apply_vmap_to_sumstats.py` preserves `.vmap` row order
- `apply_vmap_to_sumstats.py` resolves source rows by exact `source_shard + source_index`
- for imported single-file sumstats payloads, `source_shard` is `.`
- `apply_vmap_to_sumstats.py` supports single-file payloads only
- `apply_vmap_to_sumstats.py` fails cleanly on out-of-range shard-local provenance
- `apply_vmap_to_sumstats.py` defines output variant columns, including `CHR`, from target `.vmap` rows
- `apply_vmap_to_sumstats.py` writes output `SNP` according to the [generated-ID contract](spec/core-objects.md#generated-variant-ids) from the `.vmap` target row by default
- `apply_vmap_to_sumstats.py --retain-snp-id` writes target-side `.vmap` `id` values as output `SNP`
- `apply_vmap_to_sumstats.py` does not require source-side POS / SNP / effect-allele / other-allele payload columns
- `apply_vmap_to_sumstats.py` writes tab-delimited output and empty missing fields in both projection and `--clean` modes
- `apply_vmap_to_sumstats.py` writes `<output>.meta.yaml` describing delimiter, missing-value encoding, output variant columns, retained payload mappings, target metadata, and clean-mode options when applicable
- `apply_vmap_to_sumstats.py` writes `<output>.meta.yaml` that validates against `src/genomatch/schemas/raw-sumstats-metadata.yaml`
- `apply_vmap_to_sumstats.py` records retained output payload mappings in `<output>.meta.yaml` as top-level `col_*` keys rather than a nested custom mapping structure
- `apply_vmap_to_sumstats.py` writes gzip-compressed output when `--output` ends with `.gz`
- `apply_vmap_to_sumstats.py` non-clean projection mode resolves metadata-declared input columns when input header casing differs from the metadata value
- `apply_vmap_to_sumstats.py` non-clean projection mode fails before row processing when two input headers match the same metadata-declared column after trimming whitespace and case-insensitive comparison
- `apply_vmap_to_sumstats.py` non-clean projection mode drops unrecognized input columns that are not declared in metadata
- `apply_vmap_to_sumstats.py` non-clean projection mode fails before row processing when retained metadata payload columns have duplicate trimmed output names
- `apply_vmap_to_sumstats.py` non-clean projection mode fails before row processing when a retained metadata payload column would collide with canonical output variant columns `CHR`, `POS`, `SNP`, `EffectAllele`, or `OtherAllele`
- `apply_vmap_to_sumstats.py` supports joined source variant fields as input metadata but emits explicit canonical output variant columns rather than joined target fields
- `apply_vmap_to_sumstats.py` supports inputs where the only source-side variant fields are joined fields such as `CHR:POS_OtherAllele_EffectAllele`
- `import_sumstats.py` reads combined variant columns as `CHR`, `POS`, `OtherAllele`, `EffectAllele` and logs the assumed order when both alleles come from one column; for `1:15791:C:T`, verify `EffectAllele=T` and correct effect sign through `sumstats-clean`
- `apply_vmap_to_sumstats.py` accepts tab-containing headers with whitespace-delimited data rows when the field count is unambiguous
- `apply_vmap_to_sumstats.py` swaps odds-ratio confidence bounds when inverting a swapped interval
- `apply_vmap_to_sumstats.py` applies swap-style numeric transforms for both `swap` and `flip_swap`
- `apply_vmap_to_sumstats.py` rejects `@` template paths in both `--input` and `--output`
- `apply_vmap_to_sumstats.py` relies on the shared pandas reader for ragged rows; missing retained payload fields are serialized as missing values when pandas can parse the row, and parser failures are surfaced from pandas
- `apply_vmap_to_sumstats.py` warns and writes a missing value when swapped numeric effect transforms cannot be applied because the payload value is non-numeric, non-finite, or non-invertible
- swapped alleles negate signed effects
- swapped alleles invert odds ratios
- swapped alleles complement effect frequencies
- `apply_vmap_to_sumstats.py --clean` normalizes metadata-declared payload headers with the cleaned-header rules before canonical stat selection
- `apply_vmap_to_sumstats.py --clean` renames recognized payload columns to cleaned canonical names, drops unrecognized payload columns, and drops all-missing payload columns
- `apply_vmap_to_sumstats.py --clean` converts `stats_neglog10P` and `stats_log10P` payloads to plain `P` and then enforces the numeric and range validation rules
- `apply_vmap_to_sumstats.py --clean` supports both `--fill-mode column` and `--fill-mode row` with the spec-defined derivation order
- `apply_vmap_to_sumstats.py --clean` derives `EAF` from `OAF` when needed and drops `OAF`-family columns afterward
- `apply_vmap_to_sumstats.py --clean` applies the second validation pass after derivation so invalid derived values are cleared back to missing
- `apply_vmap_to_sumstats.py --clean` applies swap / flip_swap transforms to cleaned numeric effect columns while preserving `Direction` unchanged
- `apply_vmap_to_sumstats.py --clean` writes tab-delimited output with canonical cleaned payload column ordering and serializes missing values as empty fields
- `apply_vmap_to_sumstats.py --clean` covers canonical fields such as `StudyN` only as explicitly defined by [sumstats-harmonization.md](spec/sumstats-harmonization.md)
- genotype-payload `apply_vmap_*` tools (`apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py`) preserve `.vmap` row order
- genotype-payload `apply_vmap_*` tools resolve source rows by exact `source_shard + source_index`
- for genotype-payload `apply_vmap_*` tools, `source_shard` lookup is exact stored provenance and not chromosome interpretation
- genotype-payload `apply_vmap_*` tools support single-file and `@` source-prefix discovery plus optional target-contig `@` output sharding
- genotype-payload `apply_vmap_*` tools replace output `@` with the retained target-side contig label exactly as it appears in `.vmap chrom`
- genotype-payload `apply_vmap_*` tools do not emit empty `@` output shards for contigs with no retained rows
- genotype-payload `apply_vmap_*` tools reject missing required source shards and out-of-range shard-local provenance
- genotype-payload `apply_vmap_*` tools fail cleanly on empty `.vmap` input rather than emitting an empty payload
- genotype-payload `apply_vmap_*` tools do not apply smart label resolution at lookup time
- genotype-payload `apply_vmap_*` tools use the shared abstract allele-op contract: `identity` and `flip` keep the mapped row in the same genotype orientation, while `swap` and `flip_swap` perform allele-order genotype rewrites
- for `@`-sharded PLINK payloads, lookup is by exact imported `source_shard`, not by chromosome interpretation
- `apply_vmap_to_bfile.py` `@` shard discovery accepts PLINK split-X PAR shard filename aliases `XY` and `chrXY`
- `apply_vmap_to_bfile.py` defines output `.bim` rows entirely from target `.vmap` rows with genetic position / cM fixed to `0`
- `apply_vmap_to_bfile.py` writes output `.bim` SNP IDs according to the [generated-ID contract](spec/core-objects.md#generated-variant-ids) from `.vmap` target rows by default
- `apply_vmap_to_bfile.py --retain-snp-id` writes retained target-side `.vmap` `id` values as output `.bim` SNP IDs
- `apply_vmap_to_bfile.py` propagates the payload `.fam` to every emitted output when `--target-fam` is not supplied
- `apply_vmap_to_bfile.py` accepts optional `--target-fam`
- with `--target-fam`, `apply_vmap_to_bfile.py` copies that target `.fam` exactly to every emitted output and uses it as the output sample axis
- with `--target-fam`, output sample order follows the target `.fam` exactly and subjects absent from a referenced source shard are emitted as missing for rows sourced from that shard
- without `--target-fam`, for sharded source payloads `apply_vmap_to_bfile.py` requires identical `.fam` contents across all referenced source shards in one invocation
- `apply_vmap_to_bfile.py` accepts `--sample-id-mode {fid_iid,iid}` and defaults to `fid_iid`
- under `--sample-id-mode=fid_iid`, BFILE subject keys are `(FID, IID)`; under `--sample-id-mode=iid`, BFILE subject keys are `IID`
- `apply_vmap_to_bfile.py` fails on duplicate subject keys within one source shard
- `apply_vmap_to_bfile.py` fails on duplicate subject keys in `--target-fam`
- with an explicit `--target-fam`, `apply_vmap_to_bfile.py` always summarizes reconciliation-added missingness on `.vmap` rows
- with an explicit `--target-fam`, `apply_vmap_to_bfile.py` warns when any retained output subject or retained output variant exceeds 50% reconciliation-added missingness
- `apply_vmap_to_bfile.py` treats `identity` and `flip` as unchanged genotype encoding
- `apply_vmap_to_bfile.py` treats `swap` and `flip_swap` as genotype-swapping operations
- `apply_vmap_to_bfile.py` follows the shared ploidy model for expected target-side ploidy and validation
- with `--target-fam`, BFILE ploidy validation uses the sex column from the target `.fam`
- `apply_vmap_to_bfile.py` validates ploidy subject by subject when sex is known, without disabling validation for all subjects when some output subjects have unknown sex
- for BFILE, target diploid accepts observed `hom-ref (0)`, `het (1)`, `hom-alt (2)`, and missing
- for BFILE, target haploid accepts observed `hom-ref (0)`, `hom-alt (2)`, and missing, while observed `het (1)` is incompatible
- for BFILE, target absent accepts only missing and treats any nonmissing observed state as incompatible
- `apply_vmap_to_bfile.py` does not rewrite offending genotype content to missing as part of ploidy validation
- `apply_vmap_to_bfile.py` emits aggregate auditable ploidy-validation warnings rather than one warning per offending genotype
- `apply_vmap_to_bfile.py` emits `.ploidy` whenever `.vmap` rows include any row that is non-diploid in at least one sex; `.ploidy` stores one `(male_ploidy, female_ploidy)` pair per emitted row and is not a per-sample validation report
- for `@`-sharded BFILE output, `.ploidy` emission is decided per emitted shard; shards with at least one non-diploid row emit `.ploidy`, while shards containing only `(2, 2)` rows do not
- `apply_vmap_to_bfile.py` supports non-sharded output prefixes that already contain dots
- `apply_vmap_to_bfile.py` covers shard merge, shard split, and arbitrary cross-contig relocation between `@`-sharded source input and `@`-sharded output, validating that the emitted genotype calls land in the correct target shard and row
- `apply_vmap_to_bfile.py` reads and writes genotype data in bounded chunks rather than loading the full genotype matrix for all variants into memory
- bounded-chunk behavior must hold for both single-file and `@`-sharded PLINK payloads without changing output semantics, including cases where `.vmap` order interleaves rows from multiple source shards
- a real `plink2` PFILE round trip must read source and projected `.pgen` output independently of pgenlib, compare hardcalls and fractional dosages under mixed `identity` / `swap` operations (including an indel), preserve missing values, sample order, target row order, and repeated source-row use, and assert emitted `.pvar` REF/ALT against fixed target expectations
- `apply_vmap_to_pfile.py` applies the same shared genotype-payload `apply_vmap_*` contract to PLINK 2 `.pgen/.pvar/.psam`
- `apply_vmap_to_pfile.py` defines output `.pvar` rows and allele columns from retained target `.vmap` rows, not from source `.pvar`
- `apply_vmap_to_pfile.py` writes output `.pvar` IDs according to the [generated-ID contract](spec/core-objects.md#generated-variant-ids) from `.vmap` target rows by default
- assert each default-generated `.pvar` ID equals `CHROM:POS:REF:ALT` assembled independently from that emitted row
- assert emitted `.pvar` `REF=a2`, `ALT=a1` and `.bim` column 5=`a1`, column 6=`a2`; for prepared reference-aligned fixtures, assert reference allele placement directly against the reference and sumstats `EffectAllele` is ALT with correctly oriented effects, rather than deriving expected orientation from emitted payloads
- `apply_vmap_to_pfile.py --retain-snp-id` writes target-side `.vmap` `id` values as output `.pvar` IDs
- `apply_vmap_to_pfile.py` propagates the payload `.psam` to every emitted output when `--target-psam` is not supplied
- `apply_vmap_to_pfile.py` accepts optional `--target-psam`
- with `--target-psam`, `apply_vmap_to_pfile.py` copies that target `.psam` exactly to every emitted output and uses it as the output sample axis
- with `--target-psam`, output sample order follows the target `.psam` exactly and subjects absent from a referenced source shard are emitted as missing for rows sourced from that shard
- without `--target-psam`, `apply_vmap_to_pfile.py` requires identical `.psam` contents across all referenced source shards
- `apply_vmap_to_pfile.py` accepts `--sample-id-mode {fid_iid,iid}` and defaults to `fid_iid`
- under `--sample-id-mode=iid`, PFILE subject keys are `IID`
- under `--sample-id-mode=fid_iid`, PFILE subject keys are `(FID, IID)`
- under `--sample-id-mode=fid_iid`, `apply_vmap_to_pfile.py` fails clearly if referenced source `.psam` files or `--target-psam` disagree on whether the header contains `FID`
- under `--sample-id-mode=fid_iid`, empty `FID` values do not trigger fallback to `IID`
- `apply_vmap_to_pfile.py` fails on duplicate subject keys within one source shard
- `apply_vmap_to_pfile.py` fails on duplicate subject keys in `--target-psam`
- with an explicit `--target-psam`, `apply_vmap_to_pfile.py` always summarizes reconciliation-added missingness on `.vmap` rows
- with an explicit `--target-psam`, `apply_vmap_to_pfile.py` warns when any retained output subject or retained output variant exceeds 50% reconciliation-added missingness
- with `--target-psam`, `apply_vmap_to_pfile.py` copies the explicit target `.psam` exactly and does not heuristically merge arbitrary extra source `.psam` metadata columns
- `apply_vmap_to_pfile.py` supports mapped biallelic hardcalls, hardcall phase information, unphased dosages, and haploid hardcalls or haploid dosages under PLINK/PGEN allele-count conventions
- `apply_vmap_to_pfile.py` fails on unsupported mapped inputs, including non-biallelic retained rows, phased dosage preservation, and mapped content not representable through the chosen pgenlib read or write path
- `apply_vmap_to_pfile.py` copies hardcalls, hardcall phase flags, and dosages unchanged for `identity` and `flip`
- `apply_vmap_to_pfile.py` swaps biallelic hardcall allele codes and rewrites every nonmissing unphased dosage as `2 - dosage` for `swap` and `flip_swap`, while preserving missing values and phase flags
- `apply_vmap_to_pfile.py` follows the shared ploidy model for expected target-side ploidy and validation
- `apply_vmap_to_pfile.py` validates ploidy subject by subject when sex is known, without disabling validation for all subjects when some output subjects have unknown sex
- for PFILE hardcalls, ploidy validation is based on emitted allele observations after sample-axis reconciliation and allele-op rewrite; hardcall phase flags do not affect ploidy compatibility
- for PFILE hardcalls, emitted allele pairs are classified as: fully missing `(-9,-9)`, haploid-like (exactly one allele is `0` or `1` and the other is `-9`), diploid hom-ref `(0,0)`, diploid het `(0,1)` or `(1,0)`, diploid hom-alt `(1,1)`, or malformed partial state for any other partial/malformed pair
- for PFILE hardcalls, target diploid accepts fully missing, haploid-like, diploid hom-ref, diploid het, and diploid hom-alt, while malformed partial states are incompatible
- for PFILE hardcalls, target haploid accepts fully missing, haploid-like, diploid hom-ref, and diploid hom-alt, while diploid het and malformed partial states are incompatible
- for PFILE hardcalls, target absent accepts only fully missing and treats all other states as incompatible
- for PFILE dosages, unphased dosage representation is treated as `0..2` allele-count scale even on haploid chromosomes
- for `restrict_build_compatible.py --norm-indels`, the `ploidy_class_changed` comparison is against the final normalized candidate row that would otherwise be emitted, immediately before final output/QC accounting
- `apply_vmap_to_pfile.py` does not rewrite offending genotype content to missing as part of ploidy validation
- `apply_vmap_to_pfile.py` emits aggregate auditable ploidy-validation warnings rather than one warning per offending genotype
- `apply_vmap_to_pfile.py` writes `.pgen` using one coherent file-wide output channel configuration that preserves hardcalls, hardcall phase information, and unphased dosages across all emitted rows
- `apply_vmap_to_pfile.py` reads and writes genotype data in bounded chunks rather than loading the full genotype matrix for all variants into memory

### Convenience wrapper

- `prepare_variants.py` accepts optional `--prefix`, defaulting it to `--output`
- `prepare_variants.py` requires one final `--output`
- `prepare_variants.py` dispatches to the appropriate `import_*` tool for the selected input format
- `prepare_variants.py` accepts optional `--id-lookup` only for `--input-format=sumstats` and passes it through to `import_sumstats.py`
- `prepare_variants.py` accepts optional `--max-allele-length` (default `150`) and passes it through unchanged to the selected importer
- `prepare_variants.py` defaults `--dst-build` to `GRCh38`
- `prepare_variants.py` defaults `--dst-contig-naming` to `ncbi`
- cover the [preparation input/option policy](spec/workflow.md#input-variant-filtering-and-indel-options) for both wrappers, including default invocations, `--include-indels`, independent positive/negative overrides, and actual subcommand forwarding
- use mixed-length input fixtures to check QC and provenance for input filtering, plus multibase lifted output originating from class 1 input to detect unintended output filtering
- check input-filter placement against retained imported/normalized stages and build-guess input, and warning behavior for explicit operation flags without input opt-in
- verify that both wrappers forward input inclusion to the restriction primitive, whose exclusions use the build-compatible stage QC and existing resume/force lifecycle; sharded forwarding must preserve omitted versus explicit indel-operation controls
- cover all-filtered input in both wrappers with resolved source metadata, including an empty group alongside a non-empty group, all groups empty, and a destination build different from the source; check empty retained stages, final metadata, counts, and exclusion audits
- cover a clear source-build error for empty or inconclusive input before filtering without a known build, without borrowing metadata from another shard group
- `prepare_variants.py` skips `normalize_contigs.py` only when the current object already declares the requested `--dst-contig-naming`
- `prepare_variants.py` runs `normalize_contigs.py` when importer-declared `contig_naming` differs from the requested `--dst-contig-naming`
- if `prepare_variants.py --dst-contig-naming=plink_splitx` encounters `genome_build=unknown`, it defers final `plink_splitx` normalization until build is known
- in that deferred `plink_splitx` case, if the current retained stage also omits `contig_naming`, `prepare_variants.py` first normalizes to interim build-independent `plink`
- `prepare_variants.py` skips `liftover_build.py` when the source already matches `--dst-build`
- `prepare_variants.py` runs `restrict_build_compatible.py` before any optional liftover
- `prepare_variants.py` never passes the removed sorting/duplicate-cleanup flags to restriction; it always invokes `sort_variants.py` as the final transform, after optional liftover, deferred split-X normalization, and filtering, then copies the sorted stage to final output
- in the deferred `plink_splitx` case, `prepare_variants.py` runs one final `normalize_contigs.py --to plink_splitx` stage only after build is known and after any required liftover
- in the deferred `plink_splitx` case, that final normalization happens before any optional `drop_strand_ambiguous.py` or `restrict_contigs.py`
- cover retained stage naming across all input-filter and operation controls
- under `--norm-indels`, `prepare_variants.py` retains one build-compatible stage
- `prepare_variants.py` retains `<prefix>.splitx.vmap` only for the deferred final `plink_splitx` normalization case
- `prepare_variants.py` copies the last retained `current_path` to `<output>.vmap`
- without `--resume` or `--force`, `prepare_variants.py` fails if any planned output or retained intermediate already exists
- `prepare_variants.py` treats `--output` as a stem and writes the final prepared artifact exactly to `<output>.vmap`
- `prepare_variants.py` defaults omitted `--prefix` to `--output`
- `prepare_variants.py` runs `drop_strand_ambiguous.py` only when explicitly requested
- `prepare_variants.py` runs `restrict_contigs.py` only when `--chr2use` / `--contigs` is supplied
- `prepare_variants.py` prints the invoked subcommands
- `prepare_variants.py` retains intermediate outputs by default
- `prepare_variants.py` rejects `--resume` together with `--force`
- with `--resume`, `prepare_variants.py` skips a planned stage only when that stage output already exists and never overwrites outputs
- with `--resume`, `prepare_variants.py` fails on a retained stage-graph hole where a later planned stage output exists but an earlier required stage output is missing
- with `--resume`, `prepare_variants.py` does nothing if final `<output>.vmap` already exists
- with `--force`, `prepare_variants.py` deletes the full retained `<prefix>.*.vmap` stage namespace plus final `<output>.vmap` before rerunning cleanly from scratch
- for both preparation wrappers, run with `--include-indels`, then rerun with `--force` without indel inclusion or normalization; verify the build-compatible stage's `<stage-output>.bcftools_norm*` namespace is empty for every processed group
- with `--resume`, `prepare_variants.py` propagates liftover resume so a missing final lifted object can be rebuilt from retained bcftools liftover output without rerunning bcftools
- retained intermediate filenames are deterministic and derived from the chosen prefix, where omitted `--prefix` means `--output`, plus clear step-specific suffixes, including final `<prefix>.sorted.vmap`
- test that deferred `plink_splitx` normalization first uses interim `plink` when both contig naming and build are unresolved
- test that deferred `plink_splitx` normalization runs after build resolution and after any required liftover
- test that deferred `plink_splitx` normalization retains `<prefix>.splitx.vmap`
- `prepare_variants.py` presents the final prepared output at `<output>.vmap`
- test that `--output` is treated as a stem and that the final wrapper artifact is exactly `<output>.vmap`
- test that omitted `--prefix` defaults to `--output`
- test that the wrapper determines the full planned stage graph before execution and fails by default if any planned wrapper-managed output already exists
- test that `--resume` and `--force` are mutually exclusive
- test that `--resume` skips an existing stage output only when all earlier required planned stage outputs also exist
- test that `--resume` fails on a hole in the retained stage chain, e.g. planned `a -> b -> c` with `a` and `c` present but `b` missing
- test that `--resume` never overwrites outputs
- test that if final `<output>.vmap` already exists in `--resume` mode, the wrapper performs no further work
- test that `--force` deletes wrapper-managed outputs first and then reruns cleanly from scratch
- test that `prepare_variants.py --input-format sumstats --id-lookup` passes `--id-lookup` through to `import_sumstats.py`
- test that `prepare_variants.py --input-format sumstats --max-allele-length` passes `--max-allele-length` through to `import_sumstats.py`
- test that the final sort stage uses no `--drop-duplicates` and participates in preflight, resume-hole detection, `--resume`, and `--force`, including same-build, lifted, and empty preparation results
- test that `prepare_variants.py --norm-indels` completes with a single retained build-compatible stage when all retained input rows have single-base alleles

- `project_payload.py` requires raw `--input`
- `project_payload.py` requires `--input-format`
- `project_payload.py` requires `--vmap`
- `project_payload.py --vmap` must name a prepared mapped-only `.vmap`
- `project_payload.py --vmap` defines final payload variant rows, row order, IDs, source provenance, and `allele_op`
- `project_payload.py` never invokes `restrict_vmap.py`
- `project_payload.py` requires final `--output`
- for `bfile` input, `project_payload.py` accepts optional `--target-fam` and rejects `--target-psam`
- for `pfile` input, `project_payload.py` accepts optional `--target-psam` and rejects `--target-fam`
- for `bfile` and `pfile` input, `project_payload.py` accepts optional `--sample-id-mode {fid_iid,iid}` and passes it through unchanged to the canonical apply tool
- for `bfile` and `pfile` input, `project_payload.py` accepts optional `--sample-axis union`
- `project_payload.py` rejects `--sample-axis union` for `sumstats` and `sumstats-clean`
- `project_payload.py` rejects `--sample-axis union` together with an explicit `--target-fam` or `--target-psam`
- for `project_payload.py --input-format sumstats-clean`, `--fill-mode` and `--use-af-inference` are accepted and passed through unchanged to `apply_vmap_to_sumstats.py --clean`
- `project_payload.py` requires `--sumstats-metadata` for `sumstats` and `sumstats-clean` input and rejects it for `bfile` and `pfile`
- `project_payload.py` supports `sumstats`, `sumstats-clean`, `bfile`, and `pfile` input
- if `project_payload.py --sample-axis union` synthesizes a target sample file, it retains that file exactly as `<output-prefix>.target_samples.fam` for `bfile` input or `<output-prefix>.target_samples.psam` for `pfile` input
- for output prefixes containing `@`, the synthesized target-sample `<output-prefix>` is derived by replacing each `@` with `all_targets`
- `project_payload.py --sample-axis union` unions subject keys across exactly the source shards referenced by the input `.vmap`
- `project_payload.py --sample-axis union` excludes shards discovered on disk that are not referenced by input `.vmap` rows
- `project_payload.py --sample-axis union` warns and is a no-op for non-`@` input or when only one referenced shard remains
- `project_payload.py --sample-axis union` synthesizes deterministic first-occurrence subject order across participating shards in deterministic shard order
- for `bfile` input, `project_payload.py --sample-axis union` synthesizes `.fam` rows as `FID IID 0 0 sex -9`
- for `pfile` input, `project_payload.py --sample-axis union` synthesizes `.psam` with `IID`
- for `pfile` input under `--sample-id-mode=fid_iid`, `project_payload.py --sample-axis union` also synthesizes `FID` and fails clearly if participating source `.psam` files disagree on whether the header contains `FID`
- for `pfile` input, synthesized `.psam` may include `SEX`, but it must omit arbitrary other source `.psam` columns rather than heuristically merging them
- in `project_payload.py --sample-axis union`, missing sex plus known sex for the same subject key warns and keeps the known sex
- in `project_payload.py --sample-axis union`, conflicting known male and known female codes for the same subject key fail
- in `project_payload.py --sample-axis union`, duplicate subject keys within one source shard fail
- `project_payload.py` then runs the format-appropriate `apply_vmap_*` tool
- `project_payload.py` dispatches `sumstats-clean` input to `apply_vmap_to_sumstats.py --clean`
- `project_payload.py` rejects `--fill-mode` and `--use-af-inference` outside `--input-format=sumstats-clean`
- `project_payload.py` does not support `--resume`
- `project_payload.py` supports `--force`
- by default, `project_payload.py` fails if any wrapper-managed output already exists
- for `bfile` and `pfile` output prefixes containing `@`, `project_payload.py` determines wrapper-managed PLINK outputs by replacing `@` with target contig labels present in the input `.vmap`, without glob-based inference
- output planning and `--force` cleanup reuse the validated map's contigs; the standalone prefix-planning fallback streams only the chromosome column, preserving exact labels and first-occurrence order, including empty input
- `project_payload.py --target-fam`, `--target-psam`, and `--sample-id-mode` are passed through unchanged to the format-appropriate canonical apply tool
- `project_payload.py --retain-snp-id` is passed through unchanged to the format-appropriate canonical apply tool; without it, projected payload IDs use the apply-time corrected target ID default
- `project_payload.py` prints the invoked subcommands
- for `sumstats` and `sumstats-clean` input, `project_payload.py --output` is the exact rewritten payload path and must not contain `@`
- for `bfile` and `pfile` input, `project_payload.py --output` is the PLINK output prefix, may contain `@`, and is passed through to the canonical apply tool
- test that `project_payload.py` does not retain or synthesize a restricted `.vmap`
- test that for `sumstats` and `sumstats-clean` input, `--output` is treated as the exact rewritten payload path and `@` is rejected
- test that for `bfile` and `pfile` input, `--output` is treated as the PLINK output prefix and `@` is allowed
- test that `project_payload.py` does not accept `--resume`
- test that by default `project_payload.py` fails if a retained synthesized target sample file or final projected payload output already exists
- test that `--force` deletes wrapper-managed outputs first, including retained synthesized target sample files, and then reruns cleanly from scratch
- test that `project_payload.py` prints the invoked subcommands

## Fixture strategy

- Tutorial 1 regeneration must reproduce the checked-in fixture bytes from local real references. Its end-to-end regression checks documented command arguments, independently predicted prepared keys, full REF compatibility, original provenance and allele operations, QC, `.vmap`/intersection `.vtable` membership equivalence, and projected genotypes and clean summary-statistic values. Corrupted outputs must fail validation rather than passing on command exit codes alone; see `tests/integration/test_tutorial_1.py`.
- Prefer small hand-written `.vtable`, `.vmap`, and sumstats fixtures.
- Keep external-tool integration out of the default unit suite unless the canonical tool requires it.
- Canonical tests must not depend on `.match`, source-side `.vmap` metadata, dual-build BIM, or hidden cross-build behavior.
- Follow the canonical cleanup workflow and operation controls in [workflow.md](spec/workflow.md) and [variant-transforms.md](spec/variant-transforms.md).
