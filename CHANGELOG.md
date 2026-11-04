# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project is intended to follow Semantic Versioning.

## [Unreleased]

### Changed

- Repeated variant keys are now accepted in `.vmap` inputs to `intersect_variants.py` and `union_variants.py`, and in restriction inputs to `restrict_vmap.py`. Previously, these inputs could be rejected for containing duplicates.

### Pending implementation

- Remove `--sort` and `--drop-duplicates` from `restrict_build_compatible.py`; use `sort_variants.py` for explicit ordering and optional keep-first cleanup. Preparation always retains a final `<prefix>.sorted.vmap` stage, covered by preflight and resume/force handling, including when liftover runs.
- Tools accepting `.vmap` or `.vtable` follow the [common transformation and QC contract](spec/variant-transforms.md#tools-accepting-vmap-or-vtable). Operations requiring row-level QC also write sidecars for `.vtable` inputs; QC row identifiers identify input rows without adding provenance to output tables.
- Restriction now drops all `.vtable` rows whose distinct input variant keys become equal after transformation and audits them with `duplicate_target`, replacing the previous retention of all rows with duplicate output variant keys or keep-first behavior with `--drop-duplicates`. Affected outputs will contain fewer rows; this also affects SNV-only inputs.
- Revised the [reference-compatibility contract](spec/variant-transforms.md#same-build-reference-compatibility) and [preparation input/option policy](spec/workflow.md#input-variant-filtering-and-indel-options); these changes are not yet implemented.
- `restrict_build_compatible.py` and both preparation wrappers now default to SNV-only input. Commands that previously retained multibase input must add `--include-indels`, including commands already passing `--norm-indels`.
- Direct `restrict_build_compatible.py` calls no longer swap indels automatically. To permit that behavior, supply both `--include-indels` and `--allow-indel-swaps`; preparation wrappers enable indel swaps when input inclusion is requested unless explicitly disabled.
- `liftover_build.py` now trusts callers to supply source-REF-compatible `a2` and writes `REF=a2, ALT=a1` without rechecking source sequence or reselecting source allele order. If source compatibility is not guaranteed, first run `restrict_build_compatible.py` with the required input/operation options; see the [required source reference compatibility](spec/variant-transforms.md#required-source-reference-compatibility).
- Same-build `--allow-strand-flips` now applies only to input SNVs. Review multibase inputs that previously relied on strand repair; chain-directed strand changes during liftover remain supported.
- Restriction with `--include-indels` and liftover now support both-multibase inputs and outputs, subject to their reference and QC contracts. Restriction can retain these representations without `--norm-indels`; normalization no longer requires an output allele to become single-base.
- Normalization retains diagnostic VCFs and logs in the `<output>.bcftools_norm*` namespace. Suffixes and partitioning within that namespace are implementation-defined, and no aggregation is required. Every restriction run removes stale files in the namespace, including when normalization is disabled.
- QC consumers must stop expecting `norm_ambiguous_orientation` and `norm_unsupported_complex_indel`. `norm_ref_mismatch` remains a row-level exclusion when an input submitted for normalization yields no output record; unexpected observed normalization REF mismatches now fail the invocation. Liftover replaces `unsupported_non_snv` with `unsupported_representation` for unrepresentable lifted outputs, including `FLIP`/`SWAP` annotations that do not map to a supported `allele_op`; these are row-level QC exclusions rather than invocation failures.
- When upgrading to the corrected indel reference-allele assignment and normalization algorithm, delete all outputs from previous runs, including intermediate maps, retained bcftools files, membership tables, restricted maps, and projected payloads, and regenerate from the original raw inputs. Do not reuse those outputs through `--resume`. Existing maps cannot be reliably repaired by re-normalization, and the unchanged file format does not identify maps produced by the old algorithm.

## [v0.4.4] - 2026-06-08

### Fixed
- Fixed a performance regression in `intersect_variants.py` and `restrict_vmap.py` introduced in `v0.4.0`, which made large inputs much slower and more memory-intensive than necessary.

## [v0.4.3] - 2026-05-19

### Added
- Added `assign_vmap_ids.py` for replacing IDs in an existing `.vmap` after preparation. IDs are copied from a prepared `.vmap` or `.vtable` lookup by exact `chrom:pos:a1:a2` match; row order and provenance are preserved. Rows without a non-missing assigned ID are dropped by default, with `--unmatched-id-policy` available to retain unmatched rows as generated variant-key IDs or `.` IDs. Duplicate retained non-missing output IDs fail by default; use `--duplicate-id-policy` to allow them or drop all duplicated-ID rows to QC.

## [v0.4.2] - 2026-05-19

### Added
- Added `prepare_variants.py` and `prepare_variants_sharded.py --src-build` to declare a known source genome build and skip build guessing during preparation.

## [v0.4.1] - 2026-05-16

### Added
- Added summary-statistics documentation covering metadata, SNP-only imports, projection, and clean mode.
- Added `--id-lookup` for SNP-only summary-statistics preparation. It can fill missing chromosome and position values by matching summary-statistic IDs against a compatible `.vmap` or `.vtable`, then continue through the normal preparation workflow.

### Fixed
- Fixed summary-statistics parsing for files whose header contains tabs but whose data rows use whitespace-separated trailing fields, such as ILAE-style summary-statistics tables.

## [v0.4.0] - 2026-05-14

### Added
- Added expanded user documentation under `docs/`, including Tutorial 1, a worked reference-universe projection example with downloadable release fixtures; the top-level README is now a concise entry point.
- `scripts/normalize_bed_contigs_ncbi.sh` for normalizing BED contig labels to canonical NCBI-style primary contigs and reporting dropped noncanonical rows.

### Changed
- `restrict_vmap.py` replaces `match_vmap_to_target.py` as the way to trim a prepared `.vmap` to a target set.
- `project_payload.py` now projects exactly the variants in `--vmap`. To project a subset or shared set, build that `.vmap` first with `restrict_vmap.py`; use `intersect_variants.py` or `union_variants.py` first when you need to define the shared set.
- `intersect_variants.py` and `union_variants.py` now write source-independent `.vtable` outputs with `chrom:pos:a1:a2` IDs in declared coordinate order; `intersect_variants.py` inputs are now positional.

### Removed
- Removed `--only-mapped-target` from all `apply_vmap_*` tools. This is now the only apply behavior: applying a `.vmap` can no longer add genotype or summary-statistics rows that have no corresponding source row.
- Removed missing-target sentinel rows from `.vmap` objects. Unmatched variants are now left out of the `.vmap` instead of written with `source_index=-1` and `allele_op=missing`.

### Fixed
- Fixed summary-statistics metadata handling for padded column headers by matching trimmed/case-insensitive names only when the result is unambiguous; also recognizes a wider set of missing-value tokens and allows `.gz` in `path_sumStats`.

### Performance
- Reduced peak memory use in `apply_vmap_to_sumstats.py` by reading payload columns in chunks, scattering retained rows into `.vmap` order, and writing output in row blocks. Sumstats projection now emits explicit canonical variant columns, writes tab-delimited output with empty missing fields, emits `<output>.meta.yaml`, and supports repeated source provenance.
- Exact set operations can now choose smaller driver inputs using `variants_count` metadata while preserving documented output order and ID semantics.

## [v0.3.2] - 2026-05-04

### Added
- T2T-CHM13v2.0 is now supported as a declared genome build by `restrict_build_compatible.py` and `liftover_build.py`. See `docs/downloads.md` for the updated download requirements.

### Changed
- The reference configuration format has changed. Users must update existing config files to match the new structure in `config.example.yaml`.

### Known limitations
- `guess_build.py` does not infer T2T-CHM13v2.0. If your inputs are in this build, pass `--genome-build T2T-CHM13v2.0` explicitly to the `import_*` tools.

## [v0.3.1] - 2026-05-04

### Added
- `apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py` now support `--sample-axis native` to keep each emitted output shard on its native source sample axis when that shard has one unambiguous sample signature.
- `apply_vmap_to_bfile.py`, `apply_vmap_to_pfile.py`, and `project_payload.py` now support `--skip-ploidy-check` for performance-sensitive genotype projection runs where secondary genotype-content ploidy validation is not needed.

## [v0.3.0] - 2026-05-03

### Added
- `restrict_build_compatible.py --drop-duplicates` (requires `--sort`) to drop duplicate target identities after sorting.
- Projected PLINK payloads now write a `<output>.qc.tsv` audit file when ploidy-incompatible genotypes are found — for example, heterozygous calls at haploid-target positions (chrX non-PAR in males, chrY, MT) or non-missing calls at absent-target positions (chrY in females). Each row identifies the affected variant and how many samples are involved. The file is omitted when no issues are detected, and removed on re-runs that produce no issues.

### Changed
- Declared coordinate order now sorts by declared contig order, numeric position, then `a1`/`a2` lexicographically. This affects tools that emit declared coordinate order (for example `sort_variants.py`, `union_variants.py`, and liftover output ordering).
- Ploidy validation for sex-chromosome and mitochondrial variants in projected PLINK payloads is substantially faster for large cohorts, particularly on chrX.

## [v0.2.2] - 2026-04-29

### Added
- `prepare_variants_sharded.py` for memory-bounded preparation of sharded BIM/PVAR/VCF inputs.
- `sort_variants.py --drop-duplicates` to explicitly retain the first row for duplicate target identities after sorting.
- `apply_vmap_*` and `project_payload.py --retain-snp-id` to opt out of generated projected payload IDs.

### Changed
- Projected payloads now write corrected output variant IDs as `chrom:pos:a1:a2` by default.
- Reduced peak memory use in `intersect_variants.py`; the effect is largest when the first input is small.

## [v0.2.1] - 2026-04-24

### Added
- `--max-allele-length` flag (default 150) across all importers and `prepare_variants.py` wrapper. Rows where either allele exceeds the cap are dropped with QC reason `allele_too_long`.

## [v0.2.0] - 2026-04-21

### Added
- Sumstats tools can now read `path_sumStats` from metadata when `--input` is omitted.
- Import/apply sumstats flows now support headers like `#chrom`.
- FASTA access modes are now configurable (`BULK` by default, `LEGACY` via `MATCH_REFERENCE_ACCESS_MODE`).
- CLI tools now emit INFO telemetry with timestamps on stderr.

### Changed
- Improved performance in high-volume paths, including packed BED remapping for `apply_vmap_to_bfile.py` and vectorized sumstats handling `import_sumstats.py` and `apply_vmap_to_sumstats.py`.

### Fixed
- Fixed importer behavior for duplicate target identities (`chrom,pos,a1,a2`) by keeping first-seen rows and auditing dropped duplicates in QC output.

## [0.1.0] - 2026-03-23

### Added
- Initial repository layout
