# Workflow wrappers

This file is the authoritative wrapper-spec document for `prepare_variants.py`, `prepare_variants_sharded.py`, and `project_payload.py`.

These wrappers are orchestration only. They compose canonical tools into a common workflow and do not define alternate import, normalization, contig, metadata, liftover, matching, intersection, or payload-application semantics. Those semantics remain defined by the underlying canonical-tool specs, especially [importers.md](importers.md), [contigs-and-metadata.md](contigs-and-metadata.md), [variant-transforms.md](variant-transforms.md), [mapping.md](mapping.md), and [payload-application.md](payload-application.md).

## Common workflow

The common wrapper workflow ties together:

1. `prepare_variants.py` to prepare each raw input into a shared build and contig-naming convention
2. `intersect_variants.py` to derive one shared target `.vtable` from prepared inputs
3. `project_payload.py` to project each original payload into that shared target row set

In that workflow:

- `prepare_variants.py` emits one prepared `.vmap` per raw input, while retaining stage-by-stage wrapper outputs
- `prepare_variants_sharded.py` is the memory-bounded sharded-input variant of `prepare_variants.py`; it emits one final prepared `.vmap` while retaining per-shard stage outputs
- `intersect_variants.py` consumes prepared `.vmap` or `.vtable` inputs and emits one shared `.vtable`
- `project_payload.py` matches a prepared source `.vmap` to that shared target `.vtable` and rewrites the original payload into target-row order
- exact intersection, matching, and payload-application semantics remain defined by the canonical tools rather than by the wrappers

## `prepare_variants.py`

`prepare_variants.py` is the convenience wrapper for the preparation stage of the common workflow. It advances a retained stage-by-stage `.vmap` waterfall and then copies the last retained stage to `<output>.vmap`.

### Flow

```text
raw input
  |
  v
import_<format>.py
  -> <prefix>.imported.vmap
  |
  v
normalize_contigs.py        [only if current contig naming != --dst-contig-naming]
  -> <prefix>.normalized.vmap
  |
  v
guess_build.py          [in place; skipped under --resume when build is resolved]
  |
  v
restrict_build_compatible.py [--allow-strand-flips] [--norm-indels] [--sort when build already matches --dst-build]
  -> <prefix>.build_compatible.vmap
  |
  +--> liftover_build.py --target-build <dst-build>   [if build differs]
  |     -> <prefix>.lifted.vmap
  |
  +--> normalize_contigs.py --to plink_splitx         [only when requested final contig naming is plink_splitx and that final normalization was deferred until build resolution]
  |     -> <prefix>.splitx.vmap
  |
  +--> drop_strand_ambiguous.py                      [optional]
  |     -> <prefix>.strand_filtered.vmap
  |
  \--> restrict_contigs.py --chr2use/--contigs <value> [optional]
        -> <prefix>.contigs.vmap
  |
  \--> copy last retained stage
        -> <output>.vmap
```

### Contract

- `current_path` means the latest retained `.vmap` stage in the wrapper waterfall
- accept raw `--input` and `--input-format`
- for `--input-format sumstats`, `--input` may be omitted; when omitted, resolve the raw sumstats path from `path_sumStats` in `--sumstats-metadata` as `<directory of --sumstats-metadata>/<path_sumStats>`
- require `--sumstats-metadata` for `--input-format sumstats` and reject it for other importer formats
- accept optional `--id-vtable` only for `--input-format sumstats` and pass it through to `import_sumstats.py`
- accept optional `--max-allele-length` (positive integer), defaulting to `150`, and pass it through unchanged to the selected `import_*` tool
- when summary-stat metadata omits `CHR` / `POS`, require `--id-vtable`
- accept optional `--prefix`, defaulting it to `--output`
- require final `--output`
- treat `--output` as an output stem, not as a full filename
- the final prepared output artifact is exactly `<output>.vmap`
- omitted `--prefix` means the exact `--output` stem value
- default `--dst-build` to `GRCh38`
- default `--dst-contig-naming` to `ncbi`
- accept `--allow-strand-flips` / `--no-allow-strand-flips`, defaulting strand-flip allowance to enabled
- accept `--norm-indels` / `--no-norm-indels`, defaulting indel normalization to enabled
- dispatch to the appropriate `import_*` tool for the selected input format
- for `--input-format sumstats`, pass `--id-vtable` through unchanged when supplied
- pass `--max-allele-length` through unchanged when supplied
- skip `normalize_contigs.py` only when the current retained object already declares the requested `--dst-contig-naming`
- run `normalize_contigs.py` when the current retained object omits `contig_naming` or declares a different naming than requested
- exception: if `--dst-contig-naming=plink_splitx` and the current retained stage has `genome_build=unknown`, the wrapper must defer final `plink_splitx` normalization until build is known
- in that deferred `plink_splitx` case, if the current retained stage also omits `contig_naming`, the wrapper must first normalize to interim build-independent `plink` so that `guess_build.py` can run on declared contigs
- run `guess_build.py` in place on the current retained stage
- with `--resume`, skip `guess_build.py` when the current retained stage already has resolved `genome_build`
- run `restrict_build_compatible.py` before any optional liftover
- by default, pass both `--allow-strand-flips` and `--norm-indels` through to `restrict_build_compatible.py`
- if `--no-allow-strand-flips` is supplied to `prepare_variants.py`, omit `--allow-strand-flips` when invoking `restrict_build_compatible.py`
- if `--no-norm-indels` is supplied to `prepare_variants.py`, omit `--norm-indels` when invoking `restrict_build_compatible.py`
- if the current retained stage already matches `--dst-build`, pass `--sort` to `restrict_build_compatible.py`
- final `prepare_variants.py` output is therefore in declared coordinate order: either from `restrict_build_compatible.py --sort` in same-build cases or from `liftover_build.py` in build-conversion cases
- under `--norm-indels`, `prepare_variants.py` still runs exactly one retained build-compatible stage, even though `restrict_build_compatible.py` may internally make at most one `bcftools norm` call for branch 2 and at most one `bcftools norm` call for branch 3
- skip `liftover_build.py` when the current retained stage already matches `--dst-build`
- in the deferred `plink_splitx` case, after any required liftover and after build is known, run one final `normalize_contigs.py --to plink_splitx` stage before any optional `drop_strand_ambiguous.py` or `restrict_contigs.py`
- run `drop_strand_ambiguous.py` only when explicitly requested
- run `restrict_contigs.py` only when `--chr2use` / `--contigs` is supplied
- compute the full planned wrapper-managed stage graph before invoking canonical tools
- planned wrapper-managed outputs are the retained stage outputs under `<prefix>.*.vmap` plus the final `<output>.vmap`
- `<prefix>.build_compatible.vmap` remains the retained stage name regardless of whether `restrict_build_compatible.py` was invoked with `--allow-strand-flips`, `--norm-indels`, both, or neither
- `<prefix>.build_compatible.vmap` is the final same-build retained stage; `prepare_variants.py` does not retain a separate `.sorted.vmap` stage
- `<prefix>.splitx.vmap` is retained only for the deferred final `plink_splitx` normalization case
- after the requested preparation steps, copy the last retained `current_path` to `<output>.vmap`
- by default, fail if any planned output or retained intermediate for the invocation already exists
- `--resume` and `--force` are mutually exclusive
- if `--resume` is supplied, iterate through planned stages in order and skip a stage if and only if that stage output already exists
- if `--resume` is supplied, never overwrite any outputs
- if `--resume` is supplied and a later planned stage output exists while an earlier required stage output is missing, fail because the retained stage graph has a hole
- if final `<output>.vmap` already exists under `--resume`, nothing remains to do
- if `--resume` reruns a liftover step because the final lifted object is missing but retained `bcftools +liftover` output still exists, pass `--resume` through to `liftover_build.py` so it can reuse that retained bcftools output and rerun parse/QC/output logic
- if `--force` is supplied, delete all wrapper-managed outputs first and then rerun the full planned stage graph cleanly from scratch
- print the invoked subcommands
- retain intermediate outputs by default

### `prepare_variants.py`: output, preflight, `--resume`, and `--force`

- Treat `--output` as an output stem, not as a full filename. The final wrapper artifact is exactly `<output>.vmap`.
- If `--prefix` is omitted, default `--prefix` to `--output`.
- retained intermediates: `<prefix>.imported.vmap`, `.normalized.vmap`, `.build_compatible.vmap`, `.lifted.vmap`, `.splitx.vmap`, `.strand_filtered.vmap`, `.contigs.vmap`
- Before executing any canonical tool, the wrapper must determine the full planned stage graph and the full set of wrapper-managed outputs for the requested invocation.
- By default, the wrapper must fail if any planned wrapper-managed output already exists.
- `--resume` and `--force` are mutually exclusive and must not be accepted together.
- In `--resume` mode, the wrapper must iterate through the planned stage graph in order and skip a stage if and only if that stage output already exists.
- In `--resume` mode, the wrapper must never overwrite any output.
- In `--resume` mode, if a later planned stage output exists while an earlier required stage output is missing, this is an error. For example, if an `a -> b -> c` chain is planned and `a` and `c` exist but `b` is missing, the wrapper must fail rather than continue or overwrite.
- In `--resume` mode, if the final `<output>.vmap` already exists, there is nothing left to do.
- `--force` means: delete all wrapper-managed outputs for the planned invocation first, then run cleanly from scratch.
- for `prepare_variants.py`, `--force` deletes the full retained stage namespace under `<prefix>.*.vmap` plus the final `<output>.vmap`, even for optional stages that the current invocation would skip, including deferred final `plink_splitx` normalization

## `prepare_variants_sharded.py`

`prepare_variants_sharded.py` is a workflow wrapper for sharded raw genotype variant inputs. It runs the full `prepare_variants.py` stage graph independently for each discovered raw input shard, then concatenates and sorts the per-shard final `.vmap` outputs into one final `<output>.vmap`.

It is orchestration only. Stage semantics are exactly those of `prepare_variants.py` and the canonical tools it invokes.

### Contract

- accept the same user-facing preparation controls as `prepare_variants.py`, except as narrowed below
- require `--input` and require it to contain `@`
- support `--input-format` values `bim`, `pvar`, and `vcf`; reject `sumstats`
- require `--prefix` and require it to contain `@`
- require final `--output`, treat it as an output stem, and reject `@` in `--output`
- discover shards from the `--input` template using the shared raw-importer `@` discovery contract
- process discovered shards in deterministic lexical path order
- for each shard, invoke `prepare_variants.py` with a temporary one-shard `@` input template so the importer preserves the exact discovered shard token as `source_shard`
- for each shard token, replace `@` in wrapper `--prefix` with that exact token and use the resulting stem for both the per-shard retained prefix and per-shard final output
- pass all preparation controls through unchanged to each per-shard `prepare_variants.py` invocation, including build, contig naming, strand, indel, contig-filter, allele-length, `--resume`, and `--force` controls
- before final concatenation, validate that all per-shard final `.vmap` metadata agree on target `genome_build` and target `contig_naming`; mismatches fail clearly and no partial concatenated output is retained
- concatenate per-shard final `.vmap` rows in discovered shard order
- sort the concatenated rows with `sort_variants.py --drop-duplicates` into final `<output>.vmap`
- for the final `sort_variants.py` call, pass an explicit non-sharded sort scratch `--prefix` derived from wrapper `--prefix` by replacing each `@` with `all_targets` and then appending a sort-specific suffix such as `.sort_tmp`
- copy final metadata from the first shard's final `.vmap`, then set `derived_from` to the original `--input` `@` template
- leave per-shard QC sidecars in place under each per-shard retained prefix; do not concatenate them and do not write a top-level QC sidecar for the final `<output>.vmap`
- under `--resume`, if final `<output>.vmap` already exists, do nothing; otherwise skip any shard whose per-shard final `.vmap` already exists and redo final concatenation plus sort if any shard was processed
- under `--force`, first delete wrapper-managed per-shard variant objects and sidecars for all discovered shard tokens, plus final `<output>.vmap` and its metadata sidecar, then rerun cleanly

## `project_payload.py`

`project_payload.py` is the convenience wrapper for the post-prepare projection stage. `--target` always defines the target row set. If provenance must be supplied separately, the wrapper matches a source `.vmap` onto that target first, retains the matched mapping object, and then rewrites the original payload into the target row set.

### Flow 

- case 1: `--target` is `.vmap` and `--source-vmap` is omitted

  ```text
  target .vmap
      |
      \--> apply_vmap_*.py
            -> rewritten payload
  ```

- case 2: `--source-vmap` is supplied

  ```text
  prepared source .vmap         target .vtable/.vmap
            |                         |
            +-----------+-------------+
                        |
                        v
  match_vmap_to_target.py
    -> <prefix>.vmap
    |
    \--> apply_vmap_*.py
          -> rewritten payload
  ```

### Contract

- require raw `--input`, except for `--input-format sumstats` and `--input-format sumstats-clean` where `--input` may be omitted
- require `--input-format`
- require `--target`
- `--target` always defines the target row set
- `--target` must name a `.vtable` or `.vmap`
- accept optional `--source-vmap`
- if supplied, `--source-vmap` must name a prepared `.vmap`
- if `--target` names a `.vtable`, require `--source-vmap`
- if `--target` names a `.vmap`, `--source-vmap` is optional
- provenance comes from `--source-vmap` when supplied; otherwise, if `--target` is a `.vmap`, provenance comes from `--target`
- if `--source-vmap` is supplied, run exact `match_vmap_to_target.py` semantics first, with no extra normalization or relaxed allele handling
- if `--target` names a `.vmap` and `--source-vmap` is supplied, preserve exact `match_vmap_to_target.py` behavior: ignore target provenance and emit its warning
- accept optional `--prefix`, defaulting it to `--output`
- require final `--output`
- retain the matched `.vmap` from `match_vmap_to_target.py`
- the retained matched `.vmap` path is exactly `<prefix>.vmap`
- when `--target` is a `.vmap` and `--source-vmap` is omitted, no retained matched `<prefix>.vmap` is produced
- omitted `--prefix` means the exact `--output` value, except for PLINK output prefixes containing `@`, where the default retained prefix is derived by replacing each `@` with `all_targets`
- accept optional `--full-target`
- for `bfile` input, accept optional `--target-fam` and reject `--target-psam`
- for `pfile` input, accept optional `--target-psam` and reject `--target-fam`
- for `bfile` and `pfile` input, accept optional `--sample-id-mode {fid_iid,iid}` and pass it through unchanged to the underlying `apply_vmap_*` tool; default to `fid_iid`
- for `bfile` and `pfile` input, accept optional `--sample-axis union`
- reject `--sample-axis union` for `sumstats` and `sumstats-clean` input
- reject `--sample-axis union` together with an explicit `--target-fam` or `--target-psam`
- for `--input-format sumstats-clean`, accept optional `--fill-mode {column,row}` and `--use-af-inference` and pass them through unchanged to `apply_vmap_to_sumstats.py --clean`
- require `--sumstats-metadata` for `sumstats` and `sumstats-clean` input and reject it for `bfile` and `pfile`
- for `sumstats` and `sumstats-clean`, when `--input` is omitted, resolve input from `path_sumStats` in `--sumstats-metadata` as `<directory of --sumstats-metadata>/<path_sumStats>`
- support `sumstats`, `sumstats-clean`, `bfile`, and `pfile` input formats
- exception: if `--target` is a `.vmap` and `--source-vmap` is omitted, skip the match step and use `--target` directly as the mapping object for the apply step
- dispatch to `apply_vmap_to_sumstats.py` for `sumstats` input
- dispatch to the clean summary-stat apply path for `sumstats-clean` input
- for `sumstats-clean` input, pass `--fill-mode` and `--use-af-inference` through unchanged to the clean apply path
- dispatch to `apply_vmap_to_bfile.py` for `bfile` input and to `apply_vmap_to_pfile.py` for `pfile` input
- by default, pass `--only-mapped-target` to the underlying `apply_vmap_*` tool
- if `--full-target` is supplied, omit that wrapper-added `--only-mapped-target`
- if `--target-fam` or `--target-psam` is supplied, pass it through unchanged to the underlying `apply_vmap_*` tool
- if `--sample-id-mode` is supplied for `bfile` or `pfile` input, pass it through unchanged to the underlying `apply_vmap_*` tool
- for `sumstats` and `sumstats-clean` input, treat `--output` as the exact rewritten payload path and reject `@`
- for `bfile` and `pfile` input, treat `--output` as the PLINK output prefix and allow `@`
- `--prefix` must not contain `@`
- wrapper-managed outputs are the retained matched `<prefix>.vmap` plus the final rewritten payload outputs selected by `--output`
- exception: if `--target` is a `.vmap` and `--source-vmap` is omitted, wrapper-managed outputs do not include `<prefix>.vmap`
- if `--sample-axis union` synthesizes a target sample file, that synthesized file is a retained wrapper artifact written as `<prefix>.target_samples.fam` for `bfile` input or `<prefix>.target_samples.psam` for `pfile` input
- for `bfile` and `pfile` output prefixes containing `@`, determine wrapper-managed projected payload paths by replacing `@` with each target contig label present in `--target`; do not use glob-based inference
- `project_payload.py` does not support `--resume`
- `project_payload.py` supports `--force`
- by default, fail if any wrapper-managed output already exists
- `--force` means: delete all wrapper-managed outputs for the planned invocation first, including the retained matched `.vmap`, any retained synthesized target sample file, and final projected payload outputs, then run cleanly from scratch
- write the rewritten payload to `--output`
- print the invoked subcommands
- preserve the underlying `match_vmap_to_target.py`, summary-stat apply, `apply_vmap_to_bfile.py`, and `apply_vmap_to_pfile.py` semantics rather than defining a separate mapping or payload contract

### `project_payload.py`: wrapper-level `--sample-axis union`

- `--sample-axis union` is a wrapper-only convenience for `bfile` and `pfile` input
- it does not define alternate canonical payload semantics; instead, it synthesizes a target sample file and then invokes the canonical `apply_vmap_*` tool with that explicit target sample file
- when `--sample-axis union` is requested, participating source shards are the referenced retained mapped shards only
- a referenced retained mapped shard is a source shard that appears in the retained `.vmap` rows with `source_index != -1` after wrapper row-retention filtering: mapped-only by default, or full-target when `--full-target` is supplied
- shards discovered on disk but not referenced by retained mapped rows must not participate in the union
- if `--sample-axis union` is requested for non-`@` source input, warn that it is a no-op and behave as though it was omitted
- if `--sample-axis union` is requested for `@` source input but only one referenced retained mapped shard remains, warn that it is a no-op and behave as though it was omitted
- otherwise, `project_payload.py` must synthesize one explicit target sample file from the union of subject keys across the participating shards, using the selected `--sample-id-mode`
- the wrapper must retain that synthesized target sample file exactly as `<prefix>.target_samples.fam` for `bfile` input or `<prefix>.target_samples.psam` for `pfile` input
- the wrapper must then pass that retained synthesized target sample file to the canonical `apply_vmap_*` tool via `--target-fam` or `--target-psam`
- synthesized union subject order is deterministic first occurrence across participating shards in deterministic shard order
- for `bfile` input, the synthesized `.fam` must contain exactly `FID IID 0 0 sex -9`
- for `bfile` input, the synthesized `FID` and `IID` values are taken from the first occurrence of the selected subject key; father, mother, and phenotype fields are intentionally not reconciled and must be emitted as `0 0 -9`
- for `pfile` input, the synthesized `.psam` must contain `IID`
- for `pfile` input under `--sample-id-mode=fid_iid`, the synthesized `.psam` must additionally contain `FID`
- for `pfile` input under `--sample-id-mode=fid_iid`, all participating source `.psam` files must agree on whether the header contains `FID`; if they do not agree, fail clearly rather than synthesizing a mixed contract
- for `pfile` input under `--sample-id-mode=iid`, `FID` presence in participating source `.psam` files is ignored for subject matching
- for `pfile` input, the synthesized `.psam` may additionally include `SEX` when that value can be represented without conflict; arbitrary other source `.psam` columns are out of scope and must be omitted rather than heuristically merged
- in wrapper-level union synthesis, missing sex plus known sex for the same subject key is allowed; keep the known sex and warn
- in wrapper-level union synthesis, conflicting known male and known female codes for the same subject key must fail
- in wrapper-level union synthesis, duplicate subject keys within one source shard must fail
- overlap of subject keys across different participating shards is expected
- `project_payload.py` must treat retained synthesized target sample files as wrapper-managed outputs for default preflight existence checks and for `--force` cleanup

### `project_payload.py`: retained matched `.vmap`, output paths, and `--force`

- `project_payload.py` must retain the matched `.vmap` produced by `match_vmap_to_target.py`.
- exception: if `--target` is a `.vmap` and `--source-vmap` is omitted, `project_payload.py` skips `match_vmap_to_target.py` and does not write `<prefix>.vmap`.
- Accept optional `--prefix` for this retained intermediate. The retained matched mapping must be written as exactly `<prefix>.vmap`.
- `--prefix` must not contain `@` because the retained matched mapping is a canonical single-file `.vmap`.
- If `--prefix` is omitted, default `--prefix` to `--output`.
- If `--prefix` is omitted for `bfile` or `pfile` input and `--output` contains `@`, derive the default retained prefix by replacing each `@` in `--output` with `all_targets`.
- For `sumstats` and `sumstats-clean` input, treat `--output` as the exact rewritten payload path and reject `@`.
- For `bfile` and `pfile` input, treat `--output` as the PLINK output prefix and allow `@`.
- For `bfile` and `pfile` output prefixes containing `@`, determine the wrapper-managed payload outputs by replacing `@` with each target contig label present in `--target`.
- `project_payload.py` does not support `--resume`.
- `project_payload.py` supports `--force` to delete the retained matched `.vmap`, any retained synthesized target sample file, and final payload outputs first, then rerun cleanly from scratch.
- By default, `project_payload.py` must fail if any wrapper-managed output already exists, including the retained matched `.vmap`, any retained synthesized target sample file, and the final projected payload output.
- `--force` means: delete all wrapper-managed outputs for the planned invocation first, including the retained matched `.vmap`, any retained synthesized target sample file, and final projected payload outputs, then run cleanly from scratch.
- `project_payload.py` must print the invoked subcommands, just like `prepare_variants.py`.
