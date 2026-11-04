# Variant transforms

This file defines target-side row-transform semantics. Exact intersection, restriction, and union semantics remain defined in [mapping.md](mapping.md).

Expected ploidy for coordinate-changing transforms is defined in [ploidy-model.md](ploidy-model.md).

Reference-aware restriction and liftover support the explicit biallelic sequence
model in [core-objects.md](core-objects.md#invariants) for both input and output,
including rows where both alleles are multibase.

## Tools accepting `.vmap` or `.vtable`

The following tools accept either `.vmap` or `.vtable` and apply the same
rules to their target columns:

- `guess_build.py`
- `normalize_contigs.py`
- `restrict_contigs.py`
- `drop_strand_ambiguous.py`
- `restrict_build_compatible.py`, including `--norm-indels`
- `liftover_build.py`
- `sort_variants.py`

Given the same ordered canonical target rows, equivalent target metadata,
options, and reference assets, these operations must produce the same target
rows and target metadata, with the same validation, filtering, allele transformations, ordering,
duplicate handling, warnings, and failure behavior for either input type.
Each operation emits the same object type as its input, using that type's
metadata layout.

On `.vmap`, preserve original `(source_shard, source_index)` and compose local
allele operations with the existing `allele_op`. On `.vtable`, apply the same
allele transformations without recording their history or creating provenance
columns. Missing or malformed required `.vmap` provenance remains an error;
equivalence of target-row behavior does not waive object-specific validation.

Duplicate policies belong to the operation, not the input object type; the
equivalence requirement also applies in documented duplicate-cleanup modes.

This rule does not broaden the source-input interfaces of `restrict_vmap.py`,
`assign_vmap_ids.py`, or payload-application tools, which require `.vmap`.
Set materialization and provenance removal retain their explicitly defined
output types.

### Duplicate-containing input

Duplicate-containing input has multiple rows with the same
[variant key](core-objects.md#variant-keys). Unless a tool explicitly supports
such input, its behavior for that input is unspecified. Tools may reject it. No global
duplicate-validation pass is required. Explicit cleanup, set-membership, and
ambiguity-check policies remain authoritative. This also applies to source
mappings of `.vmap`-only tools.

### Row-level QC

When an operation requires row-level QC and drops rows, it writes
`<output>.qc.tsv` for either input type, using the same exclusion statuses and
reporting the same exclusion counts. QC contains dropped rows only. If no rows
are dropped, do not emit a QC sidecar and remove any stale `<output>.qc.tsv`.
Operations explicitly defined without a QC sidecar, such as `sort_variants.py`,
retain that policy for both input types.

QC identifies a dropped row using the following QC row identifiers:

- For `.vmap`, use the original `(source_shard, source_index)`; the QC columns
  are `source_shard`, `source_index`, `source_id`, and `status`.
- For `.vtable`, use `row_idx` equal to the zero-based
  input row ordinal, assigned when reading the input and preserved through
  filtering, sorting, and external-tool round trips. The QC columns are
  `row_idx`, `source_id`, and `status`; omit `source_shard` and `source_index`.

For `.vtable`, the QC row identifier refers to this invocation's input object
only; it does not establish raw-payload provenance and is never added to the output
`.vtable`. Keep `row_idx` distinct from provenance fields throughout processing;
do not rename it to `source_index` or synthesize `.vmap` provenance for `.vtable`
rows. Other QC fields retain their tool-defined meanings.

## Reference access policy

- this policy exists to improve performance and support vectorized reference-backed operations across build-aware utilities
- any toolkit path that queries FASTA-backed reference sequence is in scope, including single-base queries and queries over an entire reference allele; this policy does not require a tool to recheck source reference compatibility when the caller is responsible for ensuring it
- FASTA-consuming tools use bulk reference access by default for both single-base and variable-length span queries
- reference access mode is controlled by `MATCH_REFERENCE_ACCESS_MODE`; accepted values are `BULK` and `LEGACY` case-insensitively, with default `BULK` when unset
- the legacy per-lookup FASTA fetch path is available only when `MATCH_REFERENCE_ACCESS_MODE=LEGACY` is set
- the legacy path is intended for tests and profiling, not as the default runtime behavior
- for identical queries and reference assets, bulk and legacy access must return equivalent reference sequences and preserve the consuming tool's output, QC, and failure semantics

Span queries must cover the requested allele length from the declared one-based
start position. Queries at the same contig and position with different lengths
must remain distinct. A missing contig or a span extending outside the contig
does not supply a complete reference sequence; partial or clipped spans must
not be accepted as matches. Reference characters other than `A`, `C`, `G`, or
`T` do not match canonical alleles. These rules apply in both access modes.

Reference I/O implementation constraints are defined in
[reference FASTA access](performance-contract.md#reference-fasta-access).

## Allele case normalization

- allele strings are canonical uppercase `A` / `C` / `G` / `T` sequences across target-side rows
- all `import_*` tools must normalize retained alleles to uppercase before writing `.vmap`
- when loading user-provided `.vtable` / `.vmap`, tools must normalize loaded allele strings to uppercase before downstream processing
- FASTA sequences loaded into reference caches must be normalized to uppercase at load time
- downstream transforms should operate on already-normalized uppercase alleles and should not rely on scattered fallback uppercasing in per-row hot paths

## Strand-flip scope

- Same-build set operations are exact `chrom:pos:a1:a2` operations over already-prepared target rows; they do not perform strand-flip or swap recovery.
- Optional same-build strand-flip resolution belongs to `restrict_build_compatible.py`, which restricts same-build rows to those that are reference-compatible with the configured reference FASTA
- Strand handling in `restrict_build_compatible.py` follows the permissions and priority in [same-build reference compatibility](#same-build-reference-compatibility).
- `liftover_build.py` relies on the [required source reference compatibility](#required-source-reference-compatibility)

## Variant transforms

- `drop_strand_ambiguous.py` drops target-side strand-ambiguous biallelic rows, where `reverse_complement(a1) == a2`
- `drop_strand_ambiguous.py` applies the strand-ambiguity rule to both SNP and non-SNP biallelic alleles
- `convert_vmap_to_target.py` materializes target rows and breaks provenance intentionally
- `sort_variants.py` reads the input object in memory, sorts target rows into declared coordinate order, and does not normalize contigs, validate against reference, liftover, or change alleles
- `sort_variants.py` validates target rows against the declared `contig_naming` to compute declared coordinate order, but preserves stored `chrom` labels exactly and does not treat different contig spellings as equivalent variant keys
- `sort_variants.py` must accept optional `--drop-duplicates`; when supplied, duplicate variant keys are dropped after sorting by exact `(chrom, pos, a1, a2)`, ignoring `id` and any `.vmap` provenance fields
- `sort_variants.py --drop-duplicates` drops duplicate variant keys after sorting; because exact duplicates have identical declared-coordinate sort keys, the retained duplicate representative is the first occurrence from the input
- `sort_variants.py --drop-duplicates` does not write a QC sidecar; it is intended for explicit whole-object duplicate cleanup rather than row-level import or transform auditing

## Same-build reference compatibility

Input filtering determines which rows are eligible. Reference-allele assignment
chooses the permitted allele order and strand for each eligible row.

### Input filtering

After allele-case normalization, drop rows where `a1 == a2` with status
`identical_alleles`, before input-length filtering and reference-allele assignment,
using the shared row-level QC contract.

The three input classes are defined by the supplied allele lengths, before
reference-allele assignment or normalization; multibase means longer than one base:

| Input class | Input alleles |
| --- | --- |
| 1 | Both single-base (SNVs) |
| 2 | Exactly one multibase |
| 3 | Both multibase |

For these controls, "indel" and "non-SNV input" refer to multibase input
(classes 2 and 3), including padded SNVs and MNVs. Eligibility depends on the
supplied representation, not the variant type after normalization.

`restrict_build_compatible.py` accepts `--include-indels`, disabled by default.
Without it, exclude classes 2 and 3 before reference-allele assignment or normalization.
With it, admit all three input classes, subject to the canonical validation
rules. This is an input filter; it imposes no allele-length restriction on
transformed output. Audit these input exclusions with status `non_snv_input`
under the shared row-level QC contract.

`--allow-indel-swaps`, `--allow-strand-flips`, and `--norm-indels` are independent
opt-in options, all disabled by default in this primitive. The multibase options
`--allow-indel-swaps` and `--norm-indels` apply to classes 2 and 3.
Enabling either without `--include-indels` must issue a warning explaining that the requested multibase
operations are no-ops because multibase input is excluded. Neither operation
option enables input inclusion.

### Reference-allele assignment

`restrict_build_compatible.py` compares each proposed reference allele with the
FASTA sequence of the same length starting at the row's declared position.
A matching first base alone is insufficient. A sequence containing characters
other than `A`, `C`, `G`, or `T`, extending outside the contig, or differing from
the proposed reference allele does not match.

For every admitted input class and imported format, choose the first permitted
operation whose proposed full REF matches FASTA, in this order:

| Local operation | Proposed REF | Permission |
| --- | --- | --- |
| `identity` | `a2` | Always |
| `swap` | `a1` | Class 1, or `--allow-indel-swaps` |
| `flip` | Reverse complement of `a2` | Class 1 and `--allow-strand-flips` |
| `flip_swap` | Reverse complement of `a1` | Class 1 and `--allow-strand-flips` |

Apply the chosen [allele operation](core-objects.md#allele-operations) to obtain
the REF/ALT allele assignment, with `REF=a2` and `ALT=a1`.
If no REF/ALT allele assignment allowed by the enabled operation controls matches
the reference FASTA, drop the row with status `no_reference_match`. This includes
cases where a matching assignment requires a disabled operation.
Report exclusion reasons and assignment counts by operation; audit dropped rows
under the shared row-level QC contract. No option retains incompatible rows.
Apply normalization and final output checks only to the selected REF/ALT
allele assignment; normalization must not select or repair allele order or
strand. Failure of the selected assignment does not authorize trying another
operation.

A retained row has FASTA-compatible `a2`. When both input alleles match,
preserving the supplied order is a compatibility policy, not proof that unknown
source REF/ALT roles have been recovered.

## Same-build restriction

- `restrict_build_compatible.py` resolves reference FASTA from config + metadata only
- `restrict_build_compatible.py` accepts declared `ncbi`, `ucsc`, `plink`, and `plink_splitx` input contig naming while using UCSC internal reference FASTA
- exception: `restrict_build_compatible.py --norm-indels` rejects declared `plink_splitx` input clearly; users must normalize to a build-independent naming first before running normalization-aware same-build restriction
- input eligibility and reference compatibility follow [same-build reference compatibility](#same-build-reference-compatibility)
- `restrict_build_compatible.py` canonicalizes retained rows to reference-second target ordering: `a1=non-reference`, `a2=reference`
- `restrict_build_compatible.py` preserves retained input order and accepts neither `--sort` nor `--drop-duplicates`; explicit sorting and keep-first cleanup belong to `sort_variants.py`
- if reference-aware canonicalization makes distinct input [variant keys](core-objects.md#variant-keys) equal on output, including normalized rows matching unchanged rows, `restrict_build_compatible.py` drops all affected rows with status `duplicate_target` under the row-level QC rules
- only rows passing output validation and the ploidy check are considered when checking for duplicate output variant keys; previously excluded rows must not cause otherwise valid rows to be dropped as `duplicate_target`
- input multiallelic rows are out of scope for `restrict_build_compatible.py`; the tool operates on canonical biallelic `.vtable` / `.vmap` rows only
- `restrict_build_compatible.py` emits a JSON summary to stdout describing the input, declared metadata, reference normalization path, and retained / dropped row counts
- `restrict_build_compatible.py` does not perform liftover and is not a diagnostics-only validator
- `restrict_build_compatible.py` is the canonical same-build reference-normalization entrypoint even when `--norm-indels` is used; users are not required to run a separate pre-step before invoking `--norm-indels`

### Normalization

- without `--norm-indels`, preserve coordinates and allele representation after applying the selected operation in every input class
- with `--norm-indels`, an unexpected REF mismatch after successful reference-allele assignment, including in normalized output, is an invariant failure that fails the invocation
- class 1 input does not go through `bcftools norm`
- with `--norm-indels`, the tool may change `chrom`, `pos`, `a1`, and `a2` through reference-guided indel normalization, but those representation changes are not additional `allele_op` states
- with `--norm-indels`, if normalization changes `chrom` or `pos`, the tool must compare expected source and normalized target `(male_ploidy, female_ploidy)` pairs using the shared ploidy model and must drop rows with auditable QC status `ploidy_class_changed` when that pair changes
- with `--norm-indels`, the `bcftools norm` call must use the configured reference FASTA for the declared build
- with `--norm-indels`, the `bcftools norm` call must be normalization-only: use `-c e`; do not use multiallelic split or join modes, do not use atomization, and do not use `-c x`, `-c w`, or `-c s`
- with `--norm-indels`, each invocation runs normalization from scratch rather than supporting `--resume`
- with `--norm-indels`, if no records require normalization, skip `bcftools norm`
- with `--norm-indels`, the same internal contig-normalization model used by `liftover_build.py` applies here as well: convert declared input contigs to UCSC labels before calling `bcftools`, and convert normalized UCSC output contigs back to the declared naming afterward
- with `--norm-indels`, unsupported input contig labels still fail the whole invocation rather than being QC-filtered row by row

Normalized output must be representable as one canonical biallelic target row
per retained input row. Drop rows with the following failures, using these
status codes under the shared row-level QC contract; do not split an input row
into multiple output rows.

| Status | Failure |
| --- | --- |
| `norm_ref_mismatch` | An input row submitted for normalization yields no output record; retain this legacy status without inferring that a REF mismatch occurred |
| `norm_multiallelic` | Normalized output is multiallelic |
| `norm_not_atcg_alleles` | Normalized output contains missing, symbolic, empty, or non-`A/C/G/T` alleles |
| `norm_multiple_output_records` | One input row yields more than one normalized output record |
| `norm_identical_ref_alt_alleles` | Normalized `REF` and `ALT` are identical |
| `norm_invalid_position` | Normalized output position is invalid or non-positive |
| `unsupported_target_contig` | Normalized output contig cannot be converted back from internal UCSC naming into the declared contig naming |

These row-level exclusions do not fail the invocation. Normalization output
that cannot be parsed as VCF records or associated with known input identifiers
fails the invocation rather than being silently discarded.

### Normalization diagnostics

All diagnostic VCFs and logs produced by normalization must use the namespace
`<output>.bcftools_norm*`, beside the restriction output. Names within that
namespace and partitioning across calls or batches are implementation-defined;
retain the files as produced, without requiring aggregation.

Before processing rows on every `restrict_build_compatible.py` run, remove all
pre-existing files matching `<output>.bcftools_norm*`, whether or not
`--norm-indels` is set. This includes files from batches no longer produced by
the rerun. If no normalization runs, the namespace must remain empty.

## Liftover

`bcftools +liftover` owns the build-conversion transformation, including indel
representation, reference/alternate allele changes, and chain-directed strand
changes. Genomatch accepts its results subject to the canonical-output and
payload-projection requirements below; it does not independently reselect allele
assignments, renormalize lifted alleles, or filter strand-ambiguous results.
Explicit strand filtering remains a separate operation.

- `liftover_build.py` resolves chain and FASTA assets from config + metadata only, using list-valued directed `liftover` edges keyed by `(source, target)`
- `liftover_build.py` accepts declared `ncbi`, `ucsc`, and `plink` input contig naming while using UCSC internal reference assets
- `liftover_build.py` rejects declared `plink_splitx` input clearly; users must normalize to a build-independent naming first, liftover, and then normalize back to `plink_splitx` afterward if desired
- `liftover_build.py` supports `--resume` to reuse retained `bcftools +liftover` output at `<output>.bcftools_output.vcf` and rerun parse/QC/output logic without rerunning bcftools
- `liftover_build.py` constructs source VCF records as `REF=a2, ALT=a1`. The restriction options do not control target-build `FLIP`/`SWAP` operations, which remain governed by [Liftover allele operations](#liftover-allele-operations).
- `liftover_build.py` audits non-retained input rows under the shared row-level QC contract
- `liftover_build.py` emits declared coordinate order after liftover, using the shared declared-coordinate ordering rule

### Required source reference compatibility

Callers must supply source rows whose entire `a2` sequence matches the
source-build FASTA starting at the declared position, as defined in
[same-build reference compatibility](#same-build-reference-compatibility).
`liftover_build.py` relies on callers to ensure this compatibility and does not
perform its own source-FASTA sequence check or reselect source allele order.
Callers who cannot guarantee it must first run `restrict_build_compatible.py`
with the required input-inclusion and operation options. Preparation wrappers
already run that stage before liftover.

The caller's responsibility applies only to source reference compatibility.
Input schema, required metadata, contig naming, and reference-asset validation still apply.
Inputs violating this requirement have no guaranteed detection or
repair behavior; errors or exclusions produced by the external tool follow
the normal liftover error and QC rules.

### Output acceptance

Each retained input row must yield one canonical biallelic target row, with valid
target coordinates, distinct non-empty `A` / `C` / `G` / `T` alleles in
reference-second order (`a1=ALT`, `a2=REF`), and `FLIP`/`SWAP` annotations that
map to a supported `allele_op` under the [liftover allele-operation rules](#liftover-allele-operations).
Both alleles may be multibase. Provenance and operation composition follow the
rules for `.vmap` and `.vtable` inputs. The final candidate must preserve the
expected `(male_ploidy, female_ploidy)` pair under the shared ploidy model.

Exclude results that cannot meet these conditions with an applicable QC status:

| Status | Reason |
| --- | --- |
| `unmapped` | No lifted result for the input row |
| `unsupported_representation` | The result cannot be represented as one canonical biallelic row with valid position and sequence alleles, or its `FLIP`/`SWAP` annotations do not map to a supported `allele_op` |
| `unsupported_target_contig` | The target contig is outside the supported primary contigs or cannot be expressed in the declared naming |
| `ploidy_class_changed` | The final candidate changes the expected ploidy pair |
| `duplicate_target` | Distinct input variant keys become equal after liftover; drop all affected rows |

Only candidates passing the other acceptance conditions participate in detecting
duplicate output variant keys; previously excluded rows must not cause otherwise
valid rows to be dropped as `duplicate_target`. No further internal processing
order is prescribed.

Row-level output exclusions do not fail the invocation. Invalid input schema
or required metadata, configuration/reference-asset errors, external-tool failures, and
output that cannot be parsed or associated unambiguously with input rows fail
clearly; do not fabricate provenance or allele operations to recover from such
invocation-wide failures.

## Liftover allele operations

`bcftools +liftover` reports strand reversal with `FLIP` and a REF/ALT exchange
with `SWAP=1`. Payload projection needs this allele correspondence to recode
genotypes and effects correctly; normalized allele strings alone are
insufficient to infer it.

`liftover_build.py` must request these annotations with
`--flip-tag FLIP --swap-tag SWAP` and decode them as a local allele operation:

- no `FLIP` flag and no `SWAP` annotation: `identity`
- `FLIP` flag and no `SWAP` annotation: `flip`
- no `FLIP` flag and `SWAP=1`: `swap`
- `FLIP` flag and `SWAP=1`: `flip_swap`

These are local operations, before composition with the input `.vmap`'s
existing `allele_op` under the [input-type rules](#tools-accepting-vmap-or-vtable).

Other annotation values, including `SWAP=-1` for a new reference allele outside
the source pair and `SWAP=2`, cannot be represented by the four supported
operations for payload projection. Drop such rows with QC status
`unsupported_representation` under the [liftover acceptance contract](#liftover);
do not coerce them into an operation or infer a replacement from allele strings.
