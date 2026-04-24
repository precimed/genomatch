# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project is intended to follow Semantic Versioning.

## [Unreleased]

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
