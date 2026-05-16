# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project is intended to follow Semantic Versioning.

## [Unreleased]

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
