# Payload application

## Shared `apply_vmap_*` contract for genotype payloads

This section applies to `apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py`.

- both tools consume a genotype payload plus a `.vmap`
- both tools resolve source payload rows by exact `source_shard + source_index`; they do not normalize or reinterpret source provenance implicitly
- `source_shard` is exact stored provenance, not a biological chromosome field
- both tools use retained target-side `.vmap` rows as the definition of output variant rows
- by default, both tools preserve full target-row order from the `.vmap`
- if the same mapped source row is referenced by multiple retained target rows, both tools must allow that reuse and emit each retained target row independently in retained target-row order
- if `--only-mapped-target` is supplied, both tools must drop rows with `source_index == -1` before writing output while preserving retained relative target-row order
- if `--only-mapped-target` leaves zero retained rows, both tools must fail cleanly rather than emit an empty payload
- `@` templates are a payload-only convention; they do not apply to canonical `.vtable` or `.vmap` artifacts
- both tools support single-file source input and filename-based `@` source discovery
- both tools support target-side `@` output sharding
- when writing sharded output, replace `@` with the target-side contig label exactly as it appears in retained `.vmap` `chrom`
- do not emit empty output shards for contigs with no retained rows
- reject missing required source shards
- reject out-of-range shard-local `source_index`
- do not apply smart label resolution at lookup time; provenance lookup is exact
- preserve the `.vmap` invariant that `source_index == -1` implies `source_shard == "."` and `allele_op == "missing"`
- both tools accept optional `--sample-id-mode {fid_iid,iid}`; the default is `fid_iid`
- `--sample-id-mode` defines subject-key matching only for explicit target-sample reconciliation, not for source-variant provenance lookup

## Sample-axis reconciliation for genotype payloads

- by default, if no explicit target sample file is supplied, both tools preserve the source payload sample axis
- with a single-file source payload, the output sample axis is therefore the source sample axis
- with an `@`-sharded source payload and no explicit target sample file, all referenced source shards in one invocation must still have identical source sample-file contents; this remains the default fail-fast behavior
- both tools accept an optional explicit target sample file which, when present, defines the full output sample axis exactly
- for BFILE payloads, the explicit target sample file is `--target-fam`
- for PFILE payloads, the explicit target sample file is `--target-psam`
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
- reconciliation-added missingness is defined on retained mapped rows only
- reconciliation-added missingness is counted as retained mapped output row / output subject cells that are missing only because the selected output sample axis contains a subject absent from that row's referenced source shard
- retained unmatched target rows do not contribute to reconciliation-added missingness summaries or warnings
- when an explicit target sample file is supplied, both tools must warn if any retained output subject has more than 50% reconciliation-added missingness across retained mapped rows
- when an explicit target sample file is supplied, both tools must warn if any retained output variant has more than 50% reconciliation-added missingness across output subjects

## Abstract allele-operation level

- `identity`: retain the mapped source row in the same allele orientation relative to the target row
- `flip`: no genotype-index rewrite; only nucleotide interpretation differs
- `swap`: rewrite the genotype payload to swap the two alleles relative to target allele order
- `flip_swap`: same genotype-rewrite semantics as `swap`
- `missing`: emit an unmatched or all-missing target row if unmatched rows are retained; otherwise drop it under `--only-mapped-target`

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

- initialize an all-missing output row on the full output sample axis
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
- `apply_vmap_to_sumstats.py` does not use source-side POS / SNP / effect-allele / other-allele columns in the input payload, even when defined by the metadata
- joined source variant fields (`CHR:POS_A1_A2` or `CHR:POS`) described in the positional layout documented by `schemas/raw-sumstats-metadata.yaml` are also ignored in `apply_vmap_to_sumstats.py`; output variant columns are populated from target `.vmap` rows
- source-side variant columns and joined variant fields are dropped from the output
- `apply_vmap_to_sumstats.py` preserves full target-row order by default throughout its pipeline; payload columns for the unmatched target rows are populated by missing values
- with `--only-mapped-target`, unmatched target rows are dropped before writing output
- with `--only-mapped-target` in `--clean` mode, rows whose `P` value remains missing after harmonization are also dropped before writing output; in this context, missing `P` is treated as an all-missing payload row
- `apply_vmap_to_sumstats.py` defines output variant columns from the `.vmap` target rows, not from source payload values; these columns are named and derived as follows:
  - `CHR` from target `chrom`
  - `POS` from target `pos`
  - `SNP` from target `id`
  - `EffectAllele` from target `a1`
  - `OtherAllele` from target `a2`
- in `--clean` mode, a transformation of payload columns is performed according to `spec/sumstats-harmonization.md`,
  after the base payload-application step has matched retained output rows to the target side of a `.vmap`,
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
- without `--clean`, `apply_vmap_to_sumstats.py` preserves the input file delimiter in output
- without `--clean`, `apply_vmap_to_sumstats.py` preserves legacy payload-field formatting, including legacy missing-value encodings such as `n/a` for synthetic numeric missing values
- with `--clean`, `apply_vmap_to_sumstats.py` writes tab-delimited output
- with `--clean`, every missing value is emitted as an empty field
- `apply_vmap_to_sumstats.py` writes compressed output when `--output` ends with `.gz`; otherwise it writes plain UTF-8 text
- order of output columns is: `CHR`, `POS`, `SNP`, `EffectAllele`, `OtherAllele`, followed by payload columns
- without `--clean`, the order of payload columns after dropping source-side variant columns and joined variant fields is the same as in the input file; with `--clean`, the order is as produced by 
  the output of `spec/sumstats-harmonization.md`.

**Implementation requirement:** In `--clean` mode, write output using vectorized DataFrame operations (e.g., pandas `to_csv()`), not row-by-row iteration. Output must be tab-delimited with missing values as empty fields.

## `apply_vmap_to_bfile.py`

`apply_vmap_to_bfile.py` is one payload-specific realization of the shared genotype-payload `apply_vmap_*` contract.

Expected ploidy, payload-validation rules, and `.ploidy` semantics are defined in [ploidy-model.md](ploidy-model.md).

- the payload type is PLINK 1 BED/BIM/FAM
- `apply_vmap_to_bfile.py` defines output `.bim` rows entirely from retained target-side `.vmap` rows and writes genetic-position / cM as `0`
- `apply_vmap_to_bfile.py` accepts optional `--target-fam`
- if `--target-fam` is not supplied, `apply_vmap_to_bfile.py` must propagate the source payload `.fam` to every emitted output `.fam`
- if `--target-fam` is not supplied, all referenced source shards in one invocation must have identical `.fam` contents; implementations may enforce this as a single global precheck across all referenced shards
- if `--target-fam` is supplied, that file defines the output sample axis exactly and must be copied exactly to every emitted output `.fam`
- `identity` and `flip` leave BFILE genotype encoding unchanged
- `swap` and `flip_swap` are genotype-swapping operations
- `missing` emits all-missing genotypes for retained unmatched rows
- by default, unmatched target rows are retained as missing-genotype rows
- if every retained `.vmap` row is unmatched (`source_index = -1`), fail cleanly rather than emit an all-missing PLINK payload
- when `@` is present in the source prefix, shard discovery is filename-based and each discovered shard prefix must have matching `.bim`, `.bed`, and `.fam` components
- if the output prefix does not contain `@`, emit one PLINK payload across all retained target rows in target-row order
- if the output prefix contains `@`, emit one PLINK payload per target contig with retained rows in target-row order after optional `--only-mapped-target` filtering
- every emitted PLINK output must include `.bed`, `.bim`, and `.fam`
- `apply_vmap_to_bfile.py` emits `.ploidy` according to the shared ploidy-model contract
- if `--target-fam` is supplied, haploid validation uses the sex column from `--target-fam`
- `apply_vmap_to_bfile.py` follows the shared ploidy-model validation contract and does not redefine ploidy by rewriting offending genotype content

## Bounded-memory requirement for `apply_vmap_to_bfile.py`

BFILE payloads may be very large. Implementations must not require loading the full genotype matrix for all variants into memory at once.

Normative requirements:

- `apply_vmap_to_bfile.py` must read and write genotype data in bounded chunks
- this applies to both single-file and `@`-sharded source payloads
- for `@`-sharded source payloads, implementations may batch internally by source shard and source row ranges
- internal batching must not change output semantics
- exact target-row output semantics must still hold even when target order interleaves rows from multiple source shards
- if `--only-mapped-target` is supplied, the same bounded-chunk requirement applies after filtering unmatched target rows

## `apply_vmap_to_pfile.py`

`apply_vmap_to_pfile.py` is the PLINK 2 PFILE analogue of `apply_vmap_to_bfile.py`. It follows the same row-selection, provenance-resolution, allele-op, ordering, filtering, discovery, and sharding semantics as the shared genotype-payload `apply_vmap_*` contract, and differs only in payload-specific realization.

Expected ploidy and payload-validation rules are defined in [ploidy-model.md](ploidy-model.md).

- `apply_vmap_to_pfile.py` applies a `.vmap` to a PLINK 2 PFILE payload (`.pgen/.pvar/.psam`)
- every emitted PFILE output must include `.pgen`, `.pvar`, and `.psam`
- output `.pvar` rows are defined from retained target-side `.vmap` rows, not copied from source `.pvar`
- output allele columns in `.pvar` must match retained target-side `a1/a2`
- if a mapped row uses `swap` or `flip_swap`, the genotype payload must be rewritten to stay consistent with the target-side allele order encoded in output `.pvar`
- `apply_vmap_to_pfile.py` accepts optional `--target-psam`
- if `--target-psam` is not supplied, `.psam` is propagated from the source payload sample axis
- if `--target-psam` is not supplied, all referenced source shards in one invocation must have identical `.psam` contents; implementations may enforce `.psam` equality as one global precheck across all referenced shards
- if `--target-psam` is supplied, that file defines the output sample axis exactly and must be copied exactly to every emitted output `.psam`
- with `--target-psam`, `apply_vmap_to_pfile.py` must not attempt heuristic merging of arbitrary extra source `.psam` metadata columns; the explicit target `.psam` is the output metadata source of truth
- biallelicity of each retained mapped source row is determined from the source `.pvar` allele structure for that row
- if a retained mapped source row is non-biallelic by source `.pvar`, fail
- supported retained mapped content is: hardcalls, hardcall phase information when present, unphased dosages, and haploid hardcalls or haploid dosages under PLINK/PGEN allele-count conventions
- unsupported retained mapped content is: non-biallelic retained source rows, phased dosage preservation, and any retained mapped source content that cannot be represented by the chosen pgenlib read or write path
- fail on unsupported retained mapped inputs; do not silently degrade
- `apply_vmap_to_pfile.py` follows the shared ploidy-model validation contract; hardcall phase flags do not affect ploidy compatibility

For retained mapped rows:

- `identity` and `flip`: copy hardcalls, hardcall phase flags, and dosages through unchanged; at payload-application time, `flip` is treated exactly like `identity`
- `swap` and `flip_swap`: swap the two biallelic hardcall allele codes for every nonmissing allele observation, preserve hardcall phase flags unchanged, rewrite every nonmissing unphased dosage value as `2 - dosage`, preserve missing hardcall alleles as missing, preserve missing dosage values as missing, and apply the same swap rule to haploid calls and haploid dosages under the shared ploidy-model dosage convention; at payload-application time, `flip_swap` is treated exactly like `swap`
- no imputation, hardcall derivation from dosage, dosage derivation from hardcalls, phase inference, or call repair is performed
- payload application does not complement nucleotides; it only follows the `.vmap` `allele_op` semantics established upstream

For unmatched retained rows:

- when an unmatched target row is retained in output, `apply_vmap_to_pfile.py` must emit an all-missing `.pgen` row rather than fail
- this applies to ordinary unmatched target rows retained by default and to explicit `allele_op = missing` rows
- if every retained `.vmap` row is unmatched (`source_index = -1`), fail cleanly rather than emit an all-missing PLINK 2 payload
- `apply_vmap_to_pfile.py` must retain whatever supported genotype channels are present in the retained mapped input rows, rather than forcing one maximally rich output structure
- supported channels are hardcalls, hardcall phase information, and unphased dosages
- the tool must preserve available supported channels through to output
- if retained mapped rows are incoherent in which supported channels are present across variants, the tool should warn, still retain whatever supported channels are available for each retained mapped row, and still emit unmatched or all-missing rows coherently with the output representation used by the implementation; this inconsistency is not by itself a hard failure
- unsupported retained mapped content still causes failure
- retained mapped rows with unsupported source payload content must fail rather than degrade to missing; this includes non-biallelic retained mapped rows, retained mapped rows requiring unsupported preservation, and retained mapped rows that cannot be read or written under the chosen pgenlib path for supported channels
- an inserted all-missing row must be coherent with the emitted output representation used by the implementation rather than requiring one canonical missing-row encoding
- for every channel the implementation emits for an inserted all-missing row, every sample must be represented as missing in that channel
- exact low-level representation of all-missing hardcall-phase state is implementation-specific

Discovery and output rules:

- when `@` is present in the source prefix, shard discovery is filename-based and each discovered shard prefix must have matching `.pgen`, `.pvar`, and `.psam` components
- if source-prefix `@` discovery finds zero shards, fail cleanly
- resolve each retained mapped row to the unique discovered payload shard for its exact `source_shard`
- if the output prefix does not contain `@`, emit one PFILE payload across all retained target rows in retained target-row order
- if the output prefix contains `@`, emit one PFILE payload per target contig with retained rows
- every emitted output shard must include `.pgen`, `.pvar`, and `.psam`

## Bounded-memory requirement for `apply_vmap_to_pfile.py`

PFILE payloads may be very large. Implementations must not require loading the full genotype matrix for all variants into memory at once.

Normative requirements:

- `apply_vmap_to_pfile.py` must read and write genotype data in bounded chunks
- this applies to both single-file and `@`-sharded source payloads
- internal batching may group by source shard and source row ranges
- internal batching must not change output semantics
- exact target-row output semantics must still hold even when retained target rows interleave multiple source shards
- if `--only-mapped-target` is supplied, the same bounded-chunk requirement applies after filtering unmatched target rows
