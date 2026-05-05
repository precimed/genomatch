# Contigs and metadata

## Contig normalization

Supported `contig_naming` values in v1:

- `ncbi`
- `ucsc`
- `plink`
- `plink_splitx`

In this toolkit, `ncbi` means NCBI-style contig naming such as `1`, `2`, `X`, `Y`, `MT`, rather than UCSC-style contig naming such as `chr1`, `chrX`, `chrM`. It does not imply provenance from NCBI as a data source.

Normalization rules retained from the old toolkit:

- optional `chr` prefix
- `M` / `chrM` -> `MT`
- `par`, `par1`, `par2`, `x_par1`, `x_par2`, `x_nonpar` -> `X`
- `XY` / `chrXY` -> `X`
- PLINK codes `23 -> X`, `24 -> Y`, `25 -> X`, `26 -> MT`

## Metadata and contigs

Metadata is authoritative once present. Importers may infer `contig_naming`, but they must not rewrite contig labels.

Importers infer the first supported naming convention compatible with all retained rows, using this fixed priority order:

- `ncbi`
- `ucsc`
- `plink`
- `plink_splitx`

If imported rows are mixed or invalid:

- importers preserve retained target rows as-is
- omit `contig_naming` from metadata
- emit a warning
- require `normalize_contigs.py` before downstream use

`normalize_contigs.py` is the only repair path. It may normalize target rows from mixed or inconsistent inputs, drops any row that cannot be mapped to the chosen target convention, preserves retained row order, emits the same artifact type as the input, and reports how many rows were normalized and how many were dropped as unresolved.

For `.vmap` input, if contig normalization collapses multiple retained rows to the same target `chrom:pos:a1:a2`, `normalize_contigs.py` must drop all colliding rows and emit `<output>.qc.tsv` with status `duplicate_target`.

`normalize_contigs.py` operates on the target side only. It must not rewrite `source_shard`.

`restrict_contigs.py` is the explicit target-filtering step after normalization. It accepts either a `.vtable` or a `.vmap`, filters target rows by explicit contig selection such as `--chr2use` / `--contigs`, preserves retained target order, and emits the same artifact type as the input. For `.vmap` input it preserves original source provenance. For `.vtable` input it remains provenance-free and emits `.vtable`.

All downstream tools fail if metadata omits `contig_naming`. Valid row labels under the declared `contig_naming` remain an object invariant, but downstream tools are not required to re-validate that invariant at every boundary; tools may fail if they encounter inconsistent stored contig labels while performing their own required work.

The intended cleanup workflow is:

1. import to `.vmap`
2. if import omitted `contig_naming` or declared a different naming than the desired downstream convention, repair target contig labels with `normalize_contigs.py`
3. set or confirm metadata on the current `.vmap` using `guess_build.py` or an explicit metadata edit
4. optionally filter retained target rows with `restrict_contigs.py`
5. continue with `restrict_build_compatible.py`, liftover, sorting, match, and apply using artifacts that already exclude unresolved contigs

Wrapper note: `prepare_variants.py` may need a special-case orchestration when the requested final naming is `plink_splitx` and build is still unknown. In that case it may first normalize to build-independent `plink`, then resolve build, and only then normalize to final `plink_splitx` after any required liftover. This is a wrapper sequencing rule, not a change to canonical `normalize_contigs.py` behavior.

Reference-aware tools such as `restrict_build_compatible.py` and `liftover_build.py` still convert valid declared inputs internally to UCSC naming for FASTA or chain use, but only after strict validation against the declared metadata. This support surface includes declared `ncbi`, `ucsc`, `plink`, and `plink_splitx` input naming.

## `guess_build.py` build-guess policy

`guess_build.py` accepts either a `.vtable` or a `.vmap` and emits the same type as the input. When the input is `.vmap`, build guessing applies to its target side only and updates target-side metadata only; it does not change source provenance.

In v1, `guess_build.py` evaluates only `GRCh37` and `GRCh38`.

For each build it computes an allele/reference compatibility rate against the configured UCSC internal FASTA for that build, after internal normalization to UCSC contig naming when needed.

By default, build guessing is computed on a random downsample of up to `10,000` target rows (`--sample-rows 10000`). Sampling is without replacement and deterministic for reproducibility. Use `--sample-rows 0` to disable downsampling and evaluate all rows.

The current fixed decision policy is:

- guess the best-scoring build only if its compatibility rate is at least `0.60`
- require the best-scoring build to exceed the second-best build by at least `0.10`
- report confidence `high` when the chosen build has compatibility rate at least `0.85`
- otherwise report confidence `medium` for a chosen build
- if those thresholds are not met, keep `genome_build = "unknown"` and report confidence `low`

These thresholds are intentionally fixed in code for v1 and are not currently user-configurable.

## `--chr2use` / `--contigs`

default: no filtering

support ranges like `1-22` or `1-22,X,Y,MT`

accept `chr` prefixes and PLINK codes

normalize duplicates while preserving first occurrence

reject invalid tokens with clear errors

Should work correctly with any target contig naming, e.g. `--chr2use chr1` should work if the target rows use NCBI naming.

`--contigs` is an alias for `--chr2use`.

In the intended v1 CLI surface, `--chr2use` / `--contigs` belongs only to:

- all `import_*` tools, where it filters retained target rows before writing a canonical `.vmap`
- `restrict_contigs.py`, where it filters target rows after contig cleanup

It is not part of `apply_vmap_*` tools and is not carried across the rest of the canonical `.vtable` / `.vmap` transformation surface.

## `plink_splitx`

Requirements:

- treat this as a contig-normalization operation only; do not liftover, reorder, or alter alleles or positions
- read genome build from the input metadata, not from a new CLI flag
- require build to be GRCh37/hg19, GRCh38/hg38, or T2T-CHM13v2.0/hs1; fail clearly if missing or unsupported
- use the uploaded schema/JSON file as the authoritative source of PAR coordinates and related definitions; do not hardcode PAR boundaries if they can be read from that file

For target `plink_splitx`:

- autosomes 1-22 unchanged
- X loci inside build-specific X PAR intervals become XY
- non-PAR X loci remain X
- Y remains Y
- MT remains MT

Use inclusive interval logic (`start <= pos <= end`).

Update output metadata to `contig_naming = "plink_splitx"` and preserve genome build.
