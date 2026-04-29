# Filename shard selection

This file is the single source of truth for filename-based `@` shard selection in `match/`, including both discovered shards and explicit user-provided shard lists.

It applies to:

- raw importer input discovery for `import_bim.py`, `import_pvar.py`, and `import_vcf.py`
- wrapper-level raw shard selection for `prepare_variants.py` and `prepare_variants_sharded.py`
- source-prefix discovery for `apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py`

## Shared contract

When `@` is present in an accepted input path or source prefix, the tool treats `@` as a filename placeholder and resolves shard tokens in one of two modes.

### Mode A: explicit user-provided shards (`--shards`)

- interpret `--shards` as a comma-separated list of exact replacement tokens
- `--shards` requires the corresponding path/prefix argument to contain `@`
- reject `--shards` with empty tokens
- reject `--shards` with duplicate tokens
- for each token, construct the shard path by replacing `@` with that exact token
- fail clearly if any requested shard path does not exist
- do not run filename discovery first
- do not validate explicit tokens against the discovery token set below
- process shards in user-provided token order
- record the exact token as the shard label

### Mode B: discovered shards (no `--shards`)

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
- process discovered shards in deterministic lexical path order
- record the exact replacement token as the shard label
- reject discovery if zero shards are found

## Tool-family responsibilities

This file defines only shared filename shard selection and exact shard-label capture.

The consuming specs remain responsible for tool-family-specific behavior:

- importers define provenance assignment and `source_index` rules
- wrappers define any shard-grouping requirements beyond shared shard selection
- payload application defines component validation, exact lookup by stored `source_shard + source_index`, and output sharding rules
