# Payload application

## Output variant IDs

This section applies to all `apply_vmap_*` tools.

- payload-application tools write corrected output variant IDs by default
- the default corrected output ID is generated from the `.vmap` target row as `chrom:pos:a1:a2`
- default output IDs are generated from the `.vmap` target rows that are emitted; `.vmap` `id` values are not authoritative for the default final payload IDs
- this rule affects only the payload output ID fields: summary-stat `SNP`, PLINK 1 `.bim` column 2, and PLINK 2 `.pvar` `ID`
- provenance lookup, allele operations, and row ordering continue to ignore `id` and operate on target coordinates, alleles, and stored source provenance
- all `apply_vmap_*` tools accept optional `--retain-snp-id`; when supplied, every output row uses the target-side `.vmap` `id` value instead of the generated corrected target ID

## Shared `apply_vmap_*` contract for genotype payloads

This section applies to `apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py`.

- both tools consume a genotype payload plus a `.vmap`
- both tools resolve source payload rows by exact `source_shard + source_index`; they do not normalize or reinterpret source provenance implicitly
- `source_shard` is exact stored provenance, not a biological chromosome field
- both tools use `.vmap` rows as the exact definition of output variant rows
- both tools preserve `.vmap` row order
- if the same source row is referenced by multiple `.vmap` rows, both tools must allow that reuse and emit each `.vmap` row independently in `.vmap` row order
- empty `.vmap` input must fail cleanly rather than emit an empty payload
- `@` templates are a payload-only convention; they do not apply to canonical `.vtable` or `.vmap` artifacts
- both tools support single-file source input and filename-based `@` source discovery
- both tools support target-side `@` output sharding
- when writing sharded output, replace `@` with the target-side contig label exactly as it appears in retained `.vmap` `chrom`
- do not emit empty output shards for contigs with no retained rows
- reject missing required source shards
- reject out-of-range shard-local `source_index`
- do not apply smart label resolution at lookup time; provenance lookup is exact
- both tools accept optional `--sample-id-mode {fid_iid,iid}`; the default is `fid_iid`
- `--sample-id-mode` defines subject-key matching only for explicit target-sample reconciliation, not for source-variant provenance lookup

## Sample-axis reconciliation for genotype payloads

- by default, if no explicit target sample file is supplied, both tools preserve the source payload sample axis
- with a single-file source payload, the output sample axis is therefore the source sample axis
- with an `@`-sharded source payload and no explicit target sample file, all referenced source shards in one invocation must still have identical source sample-file contents; this remains the default fail-fast behavior
- both tools accept an optional explicit target sample file which, when present, defines the full output sample axis exactly
- for BFILE payloads, the explicit target sample file is `--target-fam`
- for PFILE payloads, the explicit target sample file is `--target-psam`
- both tools accept optional `--sample-axis native`
- `--sample-axis native` cannot be combined with `--target-fam` or `--target-psam`
- in native sample-axis mode, every emitted output shard must have one unambiguous native sample axis
- in native sample-axis mode, all mapped rows emitted to one output shard must reference source shards whose sample metadata signatures are identical
- in native sample-axis mode, different output shards may use different sample axes
- in native sample-axis mode, an emitted output shard with no mapped source rows must fail because no native sample axis is defined
- in native sample-axis mode, each output sample file is copied from the native source sample file selected for that emitted output shard
- when an explicit target sample file is supplied, output subject order is exactly the target sample-file order
- when an explicit target sample file is supplied, the tool must copy that target sample file exactly to every emitted output payload
- when an explicit target sample file is supplied, subjects absent from a given referenced source shard must be emitted as missing for rows sourced from that shard
- if a target sample file includes subjects absent from every referenced source shard, those subjects must still be retained in output and represented as missing throughout
- duplicate subject keys within one source shard must fail
- duplicate subject keys in an explicit target sample file must fail
- overlap of subject keys across different source shards is expected and does not itself imply an error
- for BFILE payloads, subject keys are derived from `.fam` as `(FID, IID)` under `--sample-id-mode=fid_iid` and as `IID` under `--sample-id-mode=iid`
- for PFILE payloads, subject keys are derived from `.psam` as `IID` under `--sample-id-mode=iid`
- for PFILE payloads under `--sample-id-mode=fid_iid`, subject keys are `(FID, IID)`
- for PFILE payloads under `--sample-id-mode=fid_iid`, all referenced source `.psam` files and any explicit target `.psam` must agree on whether the header contains `FID`; if they do not agree, fail clearly before genotype processing
- for PFILE payloads under `--sample-id-mode=fid_iid`, `FID` presence is determined from the `.psam` header, not from whether `FID` values are empty or populated
- for PFILE payloads under `--sample-id-mode=fid_iid`, if `FID` is present then empty or missing `FID` values still participate in the subject key as stored values and do not trigger fallback to `IID`
- when an explicit target sample file is supplied, both tools must always summarize reconciliation-added missingness
- reconciliation-added missingness is defined on `.vmap` rows only
- reconciliation-added missingness is counted as output row / output subject cells that are missing only because the selected output sample axis contains a subject absent from that row's referenced source shard
- when an explicit target sample file is supplied, both tools must warn if any retained output subject has more than 50% reconciliation-added missingness across `.vmap` rows
- when an explicit target sample file is supplied, both tools must warn if any retained output variant has more than 50% reconciliation-added missingness across output subjects

## Abstract allele-operation level

- `identity`: retain the mapped source row in the same allele orientation relative to the target row
- `flip`: no genotype-index rewrite; only nucleotide interpretation differs
- `swap`: rewrite the genotype payload to swap the two alleles relative to target allele order
- `flip_swap`: same genotype-rewrite semantics as `swap`

Genotype rewriting is driven by allele ordering, not by nucleotide complementation itself.

### Implementation note: sample-axis plan

The text in this subsection is non-normative implementation guidance only.

One implementation pattern that fits this contract well is to construct one shared sample-axis plan before genotype-row rewriting begins.

That plan can hold:

- the final output subject order
- the exact output sample file path to copy or synthesize
- the resolved output sex vector used for BFILE haploid validation
- for each referenced source shard, a subject-key-to-local-index lookup
- for each referenced source shard, a local-index-to-output-index scatter map
- counters needed for reconciliation-added missingness summaries and threshold warnings

With that structure, both `apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py` can follow the same high-level row algorithm:

- load the mapped source row from its referenced source shard
- scatter source-sample values into the output row by the precomputed shard-specific scatter map
- apply any allele-operation rewrite required by `swap` or `flip_swap`
- for BFILE, run haploid validation against the resolved output sex vector when enabled

This subsection is implementation guidance only; any implementation that preserves the normative behavior above is acceptable.

## `apply_vmap_to_sumstats.py`

- `apply_vmap_to_sumstats.py` resolves source rows by exact `source_shard + source_index`
- `apply_vmap_to_sumstats.py` fails cleanly on out-of-range shard-local provenance
- for single-file imported summary-stat payloads, `source_shard` is `.`
- `apply_vmap_to_sumstats.py` supports single-file payloads only; it does not allow `@` for outputs
- for `apply_vmap_to_sumstats.py`, `--input` may be omitted; in that case, use `path_sumStats` from `--sumstats-metadata`
- when `apply_vmap_to_sumstats.py --input` is omitted, resolve `path_sumStats` as `<directory of --sumstats-metadata>/<path_sumStats>`
- `path_sumStats` is filename-only (no `/`) and therefore resolves only within the metadata directory
- `apply_vmap_to_sumstats.py` supports `--clean` to emit canonical cleaned summary statistics
- without `--clean`, `apply_vmap_to_sumstats.py` is a metadata-aware projection mode: it retains only metadata-declared payload/stat columns, excluding source-side variant columns and joined variant fields, while still interpreting metadata-declared fields needed to produce a target-side projection and allele-orientation-safe effect values
- projection mode does not derive new statistical fields from other statistical fields or from metadata constants; derivation and canonical cleaned-stat completion are exclusive to `--clean`
- `apply_vmap_to_sumstats.py` does not reuse source-side POS / SNP / effect-allele / other-allele values in the input payload as output variant values, even when defined by the metadata
- source-side variant columns and joined variant fields are input-only; they are not preserved, rewritten in place, or emitted in joined form
- metadata-declared input columns must be resolved robustly, including harmless surrounding whitespace in input headers; ambiguous matches must fail rather than choosing an arbitrary column
- a metadata-declared input column match is ambiguous when more than one input header matches the metadata value after trimming surrounding whitespace and applying case-insensitive comparison
- in projection mode, metadata-declared retained payload columns are emitted under the metadata-declared column name after trimming surrounding whitespace
- projection mode must fail cleanly before reading payload rows if two retained metadata-declared payload/stat columns produce the same output name after trimming surrounding whitespace
- projection mode must fail cleanly before reading payload rows if any retained metadata-declared payload/stat column would be emitted as `CHR`, `POS`, `SNP`, `EffectAllele`, or `OtherAllele`
- unrecognized input columns are not carried through projection or clean output
- `apply_vmap_to_sumstats.py` preserves `.vmap` row order throughout its pipeline
- `apply_vmap_to_sumstats.py` always emits explicit canonical output variant columns from the `.vmap` target rows, not from source payload values; these columns are named and derived as follows:
  - `CHR` from target `chrom`
  - `POS` from target `pos`
  - `SNP` from the shared output-ID rule above
  - `EffectAllele` from target `a1`
  - `OtherAllele` from target `a2`
- in `--clean` mode, a transformation of payload columns is performed according to `spec/sumstats-harmonization.md`,
  after the base payload-application step has assembled output rows from `.vmap`,
  but before effect-direction normalization is applied based on `allele_op`.
- in `--clean` mode, `--fill-mode {column,row}` and `--use-af-inference` are exposed to the user
- `Direction` is passed through unchanged by clean harmonization and unchanged by `swap` / `flip_swap`
- the clean harmonization logic is a missing-value and missing-column completion pipeline; it does not attempt to verify semantic consistency between overlapping fields that are already present
- existing non-missing values are retained unless an explicit transform or range rule says otherwise
- the following logic applies to handle `allele_op=swap` and `allele_op=flip_swap`:
  - negate signed effects (`BETA` and `Z`)
  - invert odds ratios (`OR`)
  - invert and swap lower and upper `OR` confidence intervals `ORL95` and `ORU95`
  - swapped alleles complement effect frequencies (`EAF`, `CaseEAF`, `ControlEAF`)
  - **Implementation requirement:** Must be implemented vectorized using boolean masks and column assignment, not row-by-row loops. Create a mask for rows where `allele_op in {"swap", "flip_swap"}`, then apply operations to masked rows using vectorized assignment (e.g., `df.loc[mask, column] = ...`). Row-by-row looping with individual cell assignment is not permitted for performance reasons.
- when swapped-allele numeric effect transforms cannot be applied because the payload value is non-numeric, non-finite, or non-invertible, `apply_vmap_to_sumstats.py` must emit a warning and set the field to missing
- `apply_vmap_to_sumstats.py` writes tab-delimited output in both projection mode and `--clean` mode
- every missing output value is emitted as an empty field in both projection mode and `--clean` mode
- consumers must parse output as delimiter-separated text with an explicit tab separator; whitespace-splitting readers are unsupported because they cannot preserve empty fields
- `apply_vmap_to_sumstats.py` writes compressed output when `--output` ends with `.gz`; otherwise it writes plain UTF-8 text
- order of output columns in both projection mode and `--clean` mode is: `CHR`, `POS`, `SNP`, `EffectAllele`, `OtherAllele`, followed by payload columns
- without `--clean`, payload columns include only metadata-declared retained payload/stat columns, excluding source-side variant columns and joined variant fields; their relative order is the same as in the input file
- with `--clean`, payload columns are restricted to the canonical cleaned-stat columns produced by `spec/sumstats-harmonization.md`
- `apply_vmap_to_sumstats.py` writes a YAML metadata sidecar at `<output>.meta.yaml` for both projection mode and `--clean` mode
- `src/genomatch/schemas/raw-sumstats-metadata.yaml` is the normative schema contract for the output sidecar as well as input summary-stat metadata
- the output sidecar records retained payload column mappings using the same top-level `col_*` keys as raw input metadata; for example, a retained projected beta column is recorded as `col_BETA: <emitted column name>`, not in a nested mapping structure
- `path_sumStats` in the output sidecar is the filename-only basename of `--output`, consistent with the filename-only contract for input `path_sumStats`
- the output sidecar records at least: output `path_sumStats`, `delimiter: "\t"`, `missing_value: ""`, target `genome_build` and `contig_naming` from the `.vmap`, canonical output variant column mappings (`col_CHR: CHR`, `col_POS: POS`, `col_SNP: SNP`, `col_EffectAllele: EffectAllele`, `col_OtherAllele: OtherAllele`), retained payload `col_*` mappings, and `clean`
- in `--clean` mode, the sidecar also records `fill_mode` and `use_af_inference`

**Implementation requirement:** In both projection mode and `--clean` mode, write output using vectorized DataFrame operations (e.g., pandas `to_csv()`), not row-by-row iteration. Output must be tab-delimited with missing values as empty fields.

## `apply_vmap_to_bfile.py`

`apply_vmap_to_bfile.py` is one payload-specific realization of the shared genotype-payload `apply_vmap_*` contract.

Expected ploidy, payload-validation rules, and `.ploidy` semantics are defined in [ploidy-model.md](ploidy-model.md).

- the payload type is PLINK 1 BED/BIM/FAM
- `apply_vmap_to_bfile.py` defines output `.bim` rows from `.vmap` rows, writes genetic-position / cM as `0`, and writes the SNP field according to the shared output-ID rule above
- `apply_vmap_to_bfile.py` accepts optional `--target-fam`
- if `--target-fam` is not supplied, `apply_vmap_to_bfile.py` must propagate the source payload `.fam` to every emitted output `.fam`
- if neither `--target-fam` nor `--sample-axis native` is supplied, all referenced source shards in one invocation must have identical `.fam` contents; implementations may enforce this as a single global precheck across all referenced shards
- if `--target-fam` is supplied, that file defines the output sample axis exactly and must be copied exactly to every emitted output `.fam`
- `identity` and `flip` leave BFILE genotype encoding unchanged
- `swap` and `flip_swap` are genotype-swapping operations
- when `@` is present in the source prefix, shard discovery is filename-based and each discovered shard prefix must have matching `.bim`, `.bed`, and `.fam` components
- if the output prefix does not contain `@`, emit one PLINK payload across all `.vmap` rows in `.vmap` row order
- if the output prefix contains `@`, emit one PLINK payload per target contig with rows in `.vmap` row order
- every emitted PLINK output must include `.bed`, `.bim`, and `.fam`
- `apply_vmap_to_bfile.py` emits `.ploidy` according to the shared ploidy-model contract
- `apply_vmap_to_bfile.py` emits `<output>.qc.tsv` when ploidy-validation incompatibilities are observed in output rows
- `<output>.qc.tsv` columns are `source_shard`, `source_index`, `id`, `status`, `n_affected`
- `status` is one of `haploid_het_incompatible` or `absent_nonmissing`
- `n_affected` is the number of output samples with that incompatibility for the row
- if no retained row has ploidy-validation incompatibilities, `<output>.qc.tsv` is not emitted
- when the output prefix contains `@`, this `.qc.tsv` behavior applies per emitted output shard
- if `--target-fam` is supplied, haploid validation uses the sex column from `--target-fam`
- if `--sample-axis native` is supplied, haploid validation uses the sex column from the native `.fam` selected for each emitted output shard
- if `--skip-ploidy-check` is supplied, genotype-content ploidy validation and `.qc.tsv` incompatibility reporting are skipped; `.ploidy` emission is unchanged
- `apply_vmap_to_bfile.py` follows the shared ploidy-model validation contract and does not redefine ploidy by rewriting offending genotype content

## Bounded-memory requirement for `apply_vmap_to_bfile.py`

BFILE payloads may be very large. Implementations must not require loading the full genotype matrix for all variants into memory at once.

Normative requirements:

- `apply_vmap_to_bfile.py` must read and write genotype data in bounded chunks
- this applies to both single-file and `@`-sharded source payloads
- for `@`-sharded source payloads, implementations may batch internally by source shard and source row ranges
- internal batching must not change output semantics
- exact `.vmap` row-order output semantics must still hold even when `.vmap` order interleaves rows from multiple source shards

## `apply_vmap_to_pfile.py`

`apply_vmap_to_pfile.py` is the PLINK 2 PFILE analogue of `apply_vmap_to_bfile.py`. It follows the same row-selection, provenance-resolution, allele-op, ordering, discovery, and sharding semantics as the shared genotype-payload `apply_vmap_*` contract, and differs only in payload-specific realization.

Expected ploidy and payload-validation rules are defined in [ploidy-model.md](ploidy-model.md).

- `apply_vmap_to_pfile.py` applies a `.vmap` to a PLINK 2 PFILE payload (`.pgen/.pvar/.psam`)
- every emitted PFILE output must include `.pgen`, `.pvar`, and `.psam`
- output `.pvar` rows are defined from `.vmap` rows, not copied from source `.pvar`, with `ID` written according to the shared output-ID rule above
- output allele columns in `.pvar` must match `.vmap` target-side `a1/a2`
- if a mapped row uses `swap` or `flip_swap`, the genotype payload must be rewritten to stay consistent with the target-side allele order encoded in output `.pvar`
- `apply_vmap_to_pfile.py` accepts optional `--target-psam`
- if `--target-psam` is not supplied, `.psam` is propagated from the source payload sample axis
- if neither `--target-psam` nor `--sample-axis native` is supplied, all referenced source shards in one invocation must have identical `.psam` contents; implementations may enforce `.psam` equality as one global precheck across all referenced shards
- if `--target-psam` is supplied, that file defines the output sample axis exactly and must be copied exactly to every emitted output `.psam`
- with `--target-psam`, `apply_vmap_to_pfile.py` must not attempt heuristic merging of arbitrary extra source `.psam` metadata columns; the explicit target `.psam` is the output metadata source of truth
- biallelicity of each referenced source row is determined from the source `.pvar` allele structure for that row
- if a referenced source row is non-biallelic by source `.pvar`, fail
- supported mapped content is: hardcalls, hardcall phase information when present, unphased dosages, and haploid hardcalls or haploid dosages under PLINK/PGEN allele-count conventions
- unsupported mapped content is: non-biallelic source rows, phased dosage preservation, and any mapped source content that cannot be represented by the chosen pgenlib read or write path
- fail on unsupported mapped inputs; do not silently degrade
- if `--sample-axis native` is supplied, ploidy validation uses the sex column from the native `.psam` selected for each emitted output shard
- if `--skip-ploidy-check` is supplied, genotype-content ploidy validation and incompatibility warnings are skipped
- `apply_vmap_to_pfile.py` follows the shared ploidy-model validation contract; hardcall phase flags do not affect ploidy compatibility

For mapped rows:

- `identity` and `flip`: copy hardcalls, hardcall phase flags, and dosages through unchanged; at payload-application time, `flip` is treated exactly like `identity`
- `swap` and `flip_swap`: swap the two biallelic hardcall allele codes for every nonmissing allele observation, preserve hardcall phase flags unchanged, rewrite every nonmissing unphased dosage value as `2 - dosage`, preserve missing hardcall alleles as missing, preserve missing dosage values as missing, and apply the same swap rule to haploid calls and haploid dosages under the shared ploidy-model dosage convention; at payload-application time, `flip_swap` is treated exactly like `swap`
- no imputation, hardcall derivation from dosage, dosage derivation from hardcalls, phase inference, or call repair is performed
- payload application does not complement nucleotides; it only follows the `.vmap` `allele_op` semantics established upstream

- `apply_vmap_to_pfile.py` must retain whatever supported genotype channels are present in the mapped input rows, rather than forcing one maximally rich output structure
- supported channels are hardcalls, hardcall phase information, and unphased dosages
- the tool must preserve available supported channels through to output
- if mapped rows are incoherent in which supported channels are present across variants, the tool should warn and still retain whatever supported channels are available for each mapped row; this inconsistency is not by itself a hard failure
- unsupported mapped content still causes failure
- mapped rows with unsupported source payload content must fail rather than degrade to missing; this includes non-biallelic mapped rows, mapped rows requiring unsupported preservation, and mapped rows that cannot be read or written under the chosen pgenlib path for supported channels

Discovery and output rules:

- when `@` is present in the source prefix, shard discovery is filename-based and each discovered shard prefix must have matching `.pgen`, `.pvar`, and `.psam` components
- if source-prefix `@` discovery finds zero shards, fail cleanly
- resolve each `.vmap` row to the unique discovered payload shard for its exact `source_shard`
- if the output prefix does not contain `@`, emit one PFILE payload across all `.vmap` rows in `.vmap` row order
- if the output prefix contains `@`, emit one PFILE payload per target contig with rows in `.vmap` row order
- every emitted output shard must include `.pgen`, `.pvar`, and `.psam`

## Bounded-memory requirement for `apply_vmap_to_pfile.py`

PFILE payloads may be very large. Implementations must not require loading the full genotype matrix for all variants into memory at once.

Normative requirements:

- `apply_vmap_to_pfile.py` must read and write genotype data in bounded chunks
- this applies to both single-file and `@`-sharded source payloads
- internal batching may group by source shard and source row ranges
- internal batching must not change output semantics
- exact `.vmap` row-order output semantics must still hold even when `.vmap` rows interleave multiple source shards
