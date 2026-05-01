# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project is intended to follow Semantic Versioning.

## [Unreleased]

### Added
- `restrict_build_compatible.py --drop-duplicates` (requires `--sort`) to drop duplicate target identities after sorting.

### Changed
- Declared coordinate order now sorts by declared contig order, numeric position, then `a1`/`a2` lexicographically. This affects tools that emit declared coordinate order (for example `sort_variants.py`, `union_variants.py`, and liftover output ordering).

## [v0.2.2] - 2026-04-29

### Added
- `prepare_variants_sharded.py` for memory-bounded preparation of sharded BIM/PVAR/VCF inputs.
- `sort_variants.py --drop-duplicates` to explicitly retain the first row for duplicate target identities after sorting.
- `apply_vmap_*` and `project_payload.py --retain-snp-id` to opt out of generated projected payload IDs.

### Changed
- Projected payloads now write corrected output variant IDs as `chrom:pos:a1:a2` by default.
- Reduced peak memory use in `intersect_variants.py` and `match_vmap_to_target.py`; the effect is largest when the first input to `intersect_variants.py` is small, or when the target table for `match_vmap_to_target.py` is small.

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
- `match_vmap_to_target.py` now correctly keeps searching duplicate source candidates when an early allele match has missing provenance.
- Fixed importer behavior for duplicate target identities (`chrom,pos,a1,a2`) by keeping first-seen rows and auditing dropped duplicates in QC output.

## [0.1.0] - 2026-03-23

### Added
- Initial repository layout
