# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project is intended to follow Semantic Versioning.

## [Unreleased]

### Changed
- Speed up `apply_vmap_to_bfile.py` at biobank-scale sample counts by applying sample-axis reconciliation directly on packed BED rows instead of per-sample decode/scatter/re-encode loops.

## [0.1.0] - 2026-03-23

### Added
- Initial repository layout
