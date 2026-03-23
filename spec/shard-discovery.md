# Filename shard discovery

This file is the single source of truth for filename-based `@` shard discovery in `match/`.

It applies to:

- raw importer input discovery for `import_bim.py`, `import_pvar.py`, and `import_vcf.py`
- source-prefix discovery for `apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py`

## Shared contract

When `@` is present in an accepted input path or source prefix, the tool treats it as a filename placeholder and searches for matching shards on disk.

Discovery rules:

- replace `@` by the supported exact shard tokens
- supported exact shard tokens are:
  - `1` through `22`
  - `X`, `Y`, `MT`
  - `chr1` through `chr22`
  - `chrX`, `chrY`, `chrM`
  - `23`, `24`, `25`, `26`
  - `XY`, `chrXY`
- `23`, `24`, `25`, and `26` are PLINK-family shard labels
- `XY` and `chrXY` are accepted literal aliases for PLINK split-X PAR shard filenames
- record the exact replacement token used for each discovered shard
- that exact replacement token is the shard label; it must not be normalized, canonicalized, or reinterpreted biologically
- reject discovery if zero shards are found

## Tool-family responsibilities

This file defines only shared filename discovery and exact shard-label capture.

The consuming specs remain responsible for tool-family-specific behavior:

- importers define provenance assignment, shard concatenation order, and `source_index` rules
- payload application defines component validation, exact lookup by stored `source_shard + source_index`, and output sharding rules
