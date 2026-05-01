# Variant transforms

This file defines target-side row-transform semantics. Exact matching, intersection, and union semantics remain defined in [mapping.md](mapping.md).

Expected ploidy for coordinate-changing transforms is defined in [ploidy-model.md](ploidy-model.md).

## Reference access policy

- this policy exists to improve performance and support vectorized reference-backed operations across build-aware utilities
- any tool path that consults FASTA-backed reference bases is in scope, including `guess_build.py`, `restrict_build_compatible.py`, and `liftover_build.py`
- FASTA-consuming tools use bulk reference access by default when querying FASTA-backed reference bases
- reference access mode is controlled by `MATCH_REFERENCE_ACCESS_MODE`; accepted values are `BULK` and `LEGACY` case-insensitively, with default `BULK` when unset
- the legacy per-lookup FASTA fetch path is available only when `MATCH_REFERENCE_ACCESS_MODE=LEGACY` is set
- the legacy path is intended for tests and profiling, not as the default runtime behavior
- tools that use FASTA-backed reference access must preserve their existing validation and error semantics when switching between bulk and legacy access modes

**Implementation requirement:** bulk FASTA access should be cached per reference FASTA path and per contig, not per base position. Cache lifetime may extend for the duration of the process or tool invocation. GRCh37 and GRCh38 must use separate caches because they resolve to different FASTA paths. An implementation may eagerly load all primary contigs or lazily populate contig caches, but it must not fall back to repeated per-position `fetch()` calls in the default bulk mode.

**Implementation requirement:** in the default bulk mode, FASTA-consuming tools must not implement reference lookup as a per-variant loop around `fetch_reference_base()` or `fasta.fetch(contig, pos - 1, pos)`. They must batch by contig and answer lookups from contig-level cached reference data. The legacy per-variant fetch path is allowed only when `MATCH_REFERENCE_ACCESS_MODE=LEGACY` is set.

## Allele case normalization

- allele strings are canonical uppercase `A` / `C` / `G` / `T` sequences across target-side rows
- all `import_*` tools must normalize retained alleles to uppercase before writing `.vmap`
- when loading user-provided `.vtable` / `.vmap`, tools must normalize loaded allele strings to uppercase before downstream processing
- FASTA sequences loaded into reference caches must be normalized to uppercase at load time
- downstream transforms should operate on already-normalized uppercase alleles and should not rely on scattered fallback uppercasing in per-row hot paths

## Strand-flip scope

- `--allow-strand-flips` is not part of `match_vmap_to_target.py`
- The same-build contract is strong: `match_vmap_to_target.py` assumes source and target are already on the same build and strand convention
- Optional same-build strand-flip resolution belongs to `restrict_build_compatible.py`, which restricts same-build rows to those that are reference-compatible with the configured reference FASTA
- By default, reference-aware tools must not accept rows that are only reference-compatible after strand complementation
- If `--allow-strand-flips` is supplied to `restrict_build_compatible.py`, it may additionally consider strand-complemented alleles when evaluating same-build reference compatibility
- `liftover_build.py` expects source rows to already be source-build reference-compatible, for example after `restrict_build_compatible.py --allow-strand-flips`

## Variant transforms

- `drop_strand_ambiguous.py` drops target-side strand-ambiguous biallelic rows, where `reverse_complement(a1) == a2`
- `drop_strand_ambiguous.py` applies the strand-ambiguity rule to both SNP and non-SNP biallelic alleles
- `drop_strand_ambiguous.py` emits the same type as the input
- `.vmap` provenance preserves original `source_shard` and shard-local `source_index`
- `.vmap` metadata carries target-side metadata only
- `convert_vmap_to_target.py` materializes target rows and breaks provenance intentionally
- `sort_variants.py` accepts `.vtable` or `.vmap`, reads the input object in memory, sorts target rows into declared coordinate order, emits the same type as the input, and does not normalize contigs, validate against reference, liftover, change alleles, or change provenance
- `sort_variants.py` validates target rows against the declared `contig_naming` to compute declared coordinate order, but preserves stored `chrom` labels exactly and does not treat different contig spellings as equivalent target identities
- `sort_variants.py` must accept optional `--drop-duplicates`; when supplied, duplicate target identities are dropped after sorting by exact `(chrom, pos, a1, a2)`, ignoring `id` and any `.vmap` provenance fields
- `sort_variants.py --drop-duplicates` drops duplicate exact target identities after sorting; because exact duplicates have identical declared-coordinate sort keys, the retained duplicate representative is the first occurrence from the input
- `sort_variants.py --drop-duplicates` does not write a QC sidecar; it is intended for explicit whole-object duplicate cleanup rather than row-level import or transform auditing
- `liftover_build.py` resolves chain and FASTA assets from config + metadata only
- `liftover_build.py` accepts declared `ncbi`, `ucsc`, and `plink` input contig naming while using UCSC internal reference assets
- `liftover_build.py` rejects declared `plink_splitx` input clearly; users must normalize to a build-independent naming first, liftover, and then normalize back to `plink_splitx` afterward if desired
- `liftover_build.py` supports `--resume` to reuse retained `bcftools +liftover` output at `<output>.bcftools_output.vcf` and rerun parse/QC/output logic without rerunning bcftools
- `liftover_build.py` supports SNVs and source-build reference-anchored non-SNV rows where exactly one allele is a single-base `A` / `C` / `G` / `T` allele matching the source reference base (that single-base anchor may be either input allele) and the other allele is any non-empty `A` / `C` / `G` / `T` string; this is an input-compatibility condition, not a requirement that output canonical `a2` be single-base
- `liftover_build.py` may drop unsupported non-SNV rows with auditable QC; this includes rows that are not source-build reference-anchored in that sense and lifted outputs that are not representable as one canonical biallelic row with at least one single-base allele
- `liftover_build.py` may drop liftover results that land on unsupported non-primary target contigs with auditable QC
- `liftover_build.py` must compare expected source and lifted target `(male_ploidy, female_ploidy)` pairs using the shared ploidy model and must drop rows with auditable QC status `ploidy_class_changed` when that pair changes
- `liftover_build.py` may drop liftover results that collapse to duplicate target `chrom:pos:a1:a2` rows with auditable QC
- `liftover_build.py` writes QC rows for mapped and unmapped source variants using original `(source_shard, source_index)` provenance
- for `.vtable` input, where original provenance is absent, `liftover_build.py` writes QC rows using synthetic shard-local provenance `(".", row_index)`
- `liftover_build.py` emits declared coordinate order after liftover, using the shared declared-coordinate ordering rule
- `liftover_build.py` emits the same type as the input; it preserves original source provenance when the input is `.vmap`
- `restrict_build_compatible.py` resolves reference FASTA from config + metadata only
- `restrict_build_compatible.py` accepts declared `ncbi`, `ucsc`, `plink`, and `plink_splitx` input contig naming while using UCSC internal reference FASTA
- exception: `restrict_build_compatible.py --norm-indels` rejects declared `plink_splitx` input clearly; users must normalize to a build-independent naming first before running normalization-aware same-build restriction
- `restrict_build_compatible.py` supports both SNP and non-SNP rows when at least one allele is a single-base `A` / `C` / `G` / `T` allele matching the reference base (that single-base anchor may be either input allele) and the other allele is any non-empty `A` / `C` / `G` / `T` string; this is an input-compatibility condition, not a requirement that canonical `a2` be single-base
- `restrict_build_compatible.py --allow-strand-flips` applies that same rule after strand complementation as well
- `restrict_build_compatible.py` keeps only same-build reference-compatible rows
- `restrict_build_compatible.py` canonicalizes retained rows to reference-second target ordering: `a1=non-reference`, `a2=reference`
- `restrict_build_compatible.py` accepts optional `--sort`
- with `--sort`, `restrict_build_compatible.py` sorts retained output rows into declared coordinate order using the same contract as `sort_variants.py`
- `restrict_build_compatible.py` accepts optional `--drop-duplicates`; requires `--sort`
- with `--drop-duplicates`, after sorting, duplicate target identities are dropped by exact `(chrom, pos, a1, a2)`, retaining the first occurrence; for `.vmap` inputs, dropped duplicates are written to `<output>.qc.tsv` with status `duplicate_target`; for `.vtable` inputs, dropped duplicates are logged only (no QC sidecar, since `.vtable` carries no source provenance)
- if the reference allele is found in `a1`, `restrict_build_compatible.py` swaps the retained row and composes local `swap` into `.vmap allele_op`
- if `--allow-strand-flips` is supplied and the reference allele is found only after strand complementation, `restrict_build_compatible.py` emits local `flip` or `flip_swap` as needed while still ending with `a2=reference`
- `restrict_build_compatible.py` accepts optional `--norm-indels`
- without `--norm-indels`, `restrict_build_compatible.py` remains a same-build reference-aware filter and canonicalizer only
- with `--norm-indels`, `restrict_build_compatible.py` must partition retained rows into three branches by allele lengths after input validation
- branch 1: if both `a1` and `a2` are single-base alleles, do not call `bcftools norm`; use the existing same-build reference-aware logic only
- branch 2: if exactly one allele is longer than one base, first perform the existing same-build reference-aware compatibility logic, including optional strand-flip handling and explicit local `identity` / `swap` / `flip` / `flip_swap` tracking into `.vmap allele_op`, and then make one bulk `bcftools norm` call for that branch
- branch 3: if both `a1` and `a2` are longer than one base, this branch is available only when `--norm-indels` is supplied
- in branch 3, `--allow-strand-flips` does not apply
- in branch 3, the tool must not attempt to determine local `identity` or `swap` from sequence heuristics outside the branch-3 candidate-survival rule
- in branch 3, the tool must construct exactly two temporary VCF candidates per input row: one unswapped candidate with `REF=a2 ALT=a1`, and one swapped candidate with `REF=a1 ALT=a2`
- in branch 3, the tool must make one bulk `bcftools norm` call for that branch and map normalized output back to the originating input row and candidate orientation
- in branch 3, exactly one candidate must survive as one valid normalized biallelic row for that input row
- in branch 3, if only the unswapped candidate survives, emit local `identity`
- in branch 3, if only the swapped candidate survives, emit local `swap`
- in branch 3, if neither candidate survives, drop the input row with auditable QC
- in branch 3, if both candidates survive, drop the input row with auditable QC because local orientation is ambiguous
- `bcftools norm` is used only for normalization of already reference-compatible rows; it must not be used to discover strand flips, to repair REF/ALT orientation, or to decide `.vmap allele_op`
- with `--norm-indels`, `restrict_build_compatible.py` must still emit canonical target rows in reference-second order with `a1=non-reference` and `a2=reference`
- with `--norm-indels`, the tool may change `chrom`, `pos`, `a1`, and `a2` through reference-guided indel normalization, but those representation changes are not additional `allele_op` states
- with `--norm-indels`, any emitted local `identity` / `swap` / `flip` / `flip_swap` must come from one of the two explicit orientation-selection paths above: either the toolkit's pre-normalization compatibility logic or the branch-3 candidate-survival rule
- with `--norm-indels`, the tool must internally round-trip through temporary VCF representation with `REF=reference` and `ALT=non-reference`, then convert normalized rows back to canonical target ordering `a1=ALT`, `a2=REF`
- with `--norm-indels`, if normalization changes `chrom` or `pos`, the tool must compare expected source and normalized target `(male_ploidy, female_ploidy)` pairs using the shared ploidy model and must drop rows with auditable QC status `ploidy_class_changed` when that pair changes
- with `--norm-indels`, the `bcftools norm` call must use the configured reference FASTA for the declared build
- with `--norm-indels`, the `bcftools norm` call must be normalization-only: do not use multiallelic split or join modes, do not use atomization, and do not use `-c s`
- with `--norm-indels`, branch 2 must use `bcftools norm -c e`; a REF mismatch there is an implementation error because orientation was already resolved by the toolkit before normalization
- with `--norm-indels`, branch 3 must implement strict `bcftools norm -c x` filtering semantics; one candidate is expected to be reference-incompatible for many rows, so the implementation must exclude such candidates rather than abort the whole branch
- because upstream `bcftools` issue #2427 may segfault under `bcftools norm -c x`, the implementation may realize those semantics via a two-pass workaround: first `bcftools norm -c w`, then parse `REF_MISMATCH` warnings and rerun after excluding exactly those warned records
- with `--norm-indels`, retain the normalization VCFs next to `--output` as `<output>.bcftools_norm_input_e.vcf`, `<output>.bcftools_norm_output_e.vcf`, `<output>.bcftools_norm_input_x.vcf`, and `<output>.bcftools_norm_output_x.vcf`
- with `--norm-indels`, retain the `bcftools norm` invocation log next to `--output` as `<output>.bcftools_norm_output_e.log` for branch 2 and `<output>.bcftools_norm_output_x.log` for branch 3
- with `--norm-indels`, each invocation runs from scratch rather than supporting `--resume`; if any of those retained normalization VCFs already exist, delete them first and rerun cleanly
- with `--norm-indels`, each invocation also deletes any pre-existing retained normalization log files before rerunning
- input multiallelic rows are out of scope for `restrict_build_compatible.py`; the tool still operates on canonical biallelic `.vtable` / `.vmap` rows only
- with `--norm-indels`, if `bcftools norm` nevertheless introduces multiallelic output, those rows must not be emitted in canonical output
- with `--norm-indels`, rows whose normalized VCF output is not representable as one canonical biallelic `.vtable` row must be dropped with auditable QC rather than being split into multiple output rows
- with `--norm-indels`, the same internal contig-normalization model used by `liftover_build.py` applies here as well: convert declared input contigs to UCSC labels before calling `bcftools`, and convert normalized UCSC output contigs back to the declared naming afterward
- with `--norm-indels`, unsupported input contig labels still fail the whole invocation rather than being QC-filtered row by row
- with `--norm-indels`, "not representable as one canonical biallelic `.vtable` row" includes at least: multiallelic output, symbolic or missing alleles, non-`A/C/G/T` alleles, identical `REF` and `ALT`, multiple normalized output records for one input row, unsupported output contig labels, invalid position values, and any normalized row whose final canonical `a1` and `a2` are both longer than one base
- with `--norm-indels`, normalized rows whose final canonical `a1` and `a2` are both longer than one base must be dropped for now as unsupported complex indels; this keeps the same-build normalization contract aligned with the current `liftover_build.py` support surface
- with `--norm-indels`, duplicate-target detection still applies after normalization on final canonical `chrom:pos:a1:a2`
- if reference-aware canonicalization collapses multiple `.vmap` rows to the same target `chrom:pos:a1:a2`, `restrict_build_compatible.py` drops all colliding rows and writes `<output>.qc.tsv` with status `duplicate_target`
- with `--norm-indels`, `<output>.qc.tsv` must additionally audit rows dropped for normalization-specific reasons using explicit status codes
- the normalization-specific status vocabulary is:
- `norm_multiallelic`: normalized output is multiallelic
- `norm_not_atcg_alleles`: normalized output contains missing, symbolic, empty, or non-`A/C/G/T` alleles
- `norm_multiple_output_records`: one input row yields more than one normalized output record
- `norm_identical_ref_alt_alleles`: normalized `REF` and `ALT` are identical
- `norm_unsupported_complex_indel`: normalized output remains biallelic but both final canonical alleles are longer than one base
- `norm_ambiguous_orientation`: both branch-3 candidates survive
- `norm_ref_mismatch`: no valid normalized candidate survives because reference anchoring fails under normalization-only `bcftools norm`
- `norm_invalid_position`: normalized output position is invalid or non-positive
- `unsupported_target_contig`: normalized output contig cannot be converted back from internal UCSC naming into the declared contig naming
- `ploidy_class_changed`: normalized output changes the expected `(male_ploidy, female_ploidy)` pair relative to the input row
- `restrict_build_compatible.py` emits filtered output of the same type as the input; it preserves original source provenance when the input is `.vmap`
- `restrict_build_compatible.py` emits a JSON summary to stdout describing the input, declared metadata, reference normalization path, and retained / dropped row counts
- `restrict_build_compatible.py` does not perform liftover and is not a diagnostics-only validator
- `restrict_build_compatible.py` is the canonical same-build reference-normalization entrypoint even when `--norm-indels` is used; users are not required to run a separate pre-step before invoking `--norm-indels`

## Liftover allele operations

For biallelic SNVs, `liftover_build.py` should derive `allele_op` from `bcftools +liftover` output rather than inferring it only from the emitted alleles.

To make those states explicit and auditable, the toolkit should invoke `bcftools +liftover` with:

- `--flip-tag FLIP`
- `--swap-tag SWAP`

The canonical `allele_op` mapping must then be derived from those emitted annotations.

The intended mapping is:

- no `FLIP` flag and no `SWAP` annotation: `identity`
- `FLIP` flag and no `SWAP` annotation: `flip`
- no `FLIP` flag and `SWAP=1`: `swap`
- `FLIP` flag and `SWAP=1`: `flip_swap`

If liftover produces a case not representable by canonical v1 `allele_op` values, such as `SWAP=-1` for a new reference allele outside the original source allele pair, the toolkit should fail clearly or reject that variant with auditable QC rather than coercing it into a canonical op.
