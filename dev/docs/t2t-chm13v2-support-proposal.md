# Proposal: T2T-CHM13v2.0 support

Status: draft proposal, not yet normative.

## Goal

Add explicit support for `T2T-CHM13v2.0` as a destination and same-build
reference genome in the existing `.vmap` / `.vtable` workflow.

This proposal is intentionally scoped to the named build
`T2T-CHM13v2.0`. It is not a general plugin model for arbitrary genome
assemblies. The implementation may hard-code this build as one additional
supported value wherever the current public surface enumerates genome builds.
T2T should be a first-class supported build across canonical tools and
workflows, not a special-case liftover-only target.

Initial supported conversion edges:

- `GRCh37` -> `T2T-CHM13v2.0`, using `ref/chain/hg19-chm13v2.chain.gz`
- `GRCh38` -> `T2T-CHM13v2.0`, using `ref/chain/grch38-chm13v2.chain.gz`
- `T2T-CHM13v2.0` -> `GRCh37`, using `ref/chain/chm13v2-hg19.chain.gz`
- `T2T-CHM13v2.0` -> `GRCh38`, using `ref/chain/chm13v2-grch38.chain.gz`

## Reference assets

The initial local assets are:

```text
ref/chain/hg19-chm13v2.chain.gz
ref/chain/grch38-chm13v2.chain.gz
ref/chain/chm13v2-hg19.chain.gz
ref/chain/chm13v2-grch38.chain.gz
ref/ucsc/T2T-CHM13v2.0/chm13v2.0.fa
```

Download sources:

- Reference FASTA: `https://42basepairs.com/browse/s3/human-pangenomics/T2T/CHM13/assemblies?file=chm13v2.0.fa`
- Chain files: `https://42basepairs.com/browse/s3/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo`
- Contig aliases: `https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/GCA_009914755.4.chromAlias.txt`

The current FASTA and chain assets use UCSC-style contig labels. The local
chain filenames above are the configured runtime filenames, not a strict
download-name contract. Upstream chain sources may provide plain `.chain` files
or differently named `.over.chain.gz` files; users may download, gzip, and name
them to match their `liftover` config entries.

For `T2T-CHM13v2.0` / UCSC `hs1`, the practical primary contigs are:

- UCSC labels: `chr1` through `chr22`, `chrX`, `chrY`, `chrM`
- NCBI-style aliases: `1` through `22`, `X`, `Y`, `MT`
- Accession-style aliases also exist and should remain metadata/reference
  information, not a new public `contig_naming` mode for this change.

`chrY` in `chm13v2.0.fa` is from HG002/NA24385 rather than the CHM13 cell
line. A `chm13v2.0_noY` FASTA is out of scope for this proposal because the
toolkit's primary-contig contract is currently hard-coded to include `Y`.

## Proposed config model

Replace the current split `references.ucsc`, `references.ncbi`, and top-level
magic-keyed `chain` blocks with two explicit concepts:

- `builds`: one required UCSC-style FASTA asset per supported build.
- `liftover`: a list of directed build-pair edges. Each edge is exactly the
  triple `source`, `target`, and `chain`.

The toolkit keeps the existing hard-coded primary-contig contract:

- UCSC reference labels are `chr1` through `chr22`, `chrX`, `chrY`, `chrM`.
- NCBI-style object labels are `1` through `22`, `X`, `Y`, `MT`.
- PLINK-style object labels remain derived from that same primary-contig set.

NCBI-style FASTA assets should be removed from the active config example. They
may exist locally for inspection or fixture generation, but reference-aware
tools should not resolve them.

Proposed `config.yaml` shape:

```yaml
builds:
  GRCh37:
    ucsc_fasta: ucsc/GRCh37/hg19.p13.plusMT.no_alt_analysis_set.fa
  GRCh38:
    ucsc_fasta: ucsc/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
  T2T-CHM13v2.0:
    ucsc_fasta: ucsc/T2T-CHM13v2.0/chm13v2.0.fa

liftover:
  - source: GRCh37
    target: GRCh38
    chain: chain/hg19ToHg38.over.chain.gz
  - source: GRCh38
    target: GRCh37
    chain: chain/hg38ToHg19.over.chain.gz
  - source: GRCh37
    target: T2T-CHM13v2.0
    chain: chain/hg19-chm13v2.chain.gz
  - source: GRCh38
    target: T2T-CHM13v2.0
    chain: chain/grch38-chm13v2.chain.gz
  - source: T2T-CHM13v2.0
    target: GRCh37
    chain: chain/chm13v2-hg19.chain.gz
  - source: T2T-CHM13v2.0
    target: GRCh38
    chain: chain/chm13v2-grch38.chain.gz
```

Notes:

- Config files may omit T2T entries when T2T workflows are not used. T2T FASTA
  and chain assets are required only for operations whose declared source or
  target build is `T2T-CHM13v2.0`.
- Chain files and reference FASTA files are always UCSC-style in the active
  runtime contract. Temporary VCFs used for `bcftools` are therefore also
  UCSC-style.
- The edge identity is the ordered pair `(source, target)`, not a derived key
  such as `hg19ToHg38`.
- Paths remain absolute or relative to the config file, matching the current
  `MATCH_CONFIG` behavior.
- The legacy config shape does not need a transition path. This is a breaking
  config change; users must update to `builds` and list-valued `liftover`.

## Reference resolution behavior

`reference_utils.py` should expose two resolution paths.

Same-build FASTA resolution:

- Input: `genome_build`.
- Resolve exactly one UCSC FASTA from `builds.<genome_build>.ucsc_fasta`.
- Fail clearly if no FASTA exists for the build, if the requested build is
  unsupported, or if the FASTA index does not contain the hard-coded UCSC
  primary-contig set required by the toolkit.
- In particular, `restrict_build_compatible.py` should fail clearly when the
  declared build is `T2T-CHM13v2.0` but the config has no T2T FASTA entry.

Liftover asset resolution:

- Input: `source_build`, `target_build`.
- Find exactly one `liftover` edge whose `source` and `target` match.
- Resolve the chain path from that edge.
- Resolve source and target UCSC FASTA paths from `builds`.
- Fail if the configured edge points to missing assets or to FASTA assets that
  do not satisfy the hard-coded UCSC primary-contig contract.
- In particular, `liftover_build.py` should fail clearly when either source or
  target build is `T2T-CHM13v2.0` but the config lacks the required T2T FASTA
  or directed T2T chain edge.

This keeps the current GRCh37/GRCh38 UCSC-chain workflow and makes the UCSC
dependency explicit instead of pretending NCBI-style FASTA assets are active
runtime alternatives.

## Tool behavior changes

### `restrict_build_compatible.py`

Current behavior always converts declared contigs to UCSC labels before FASTA
lookup. This should remain the required behavior.

Required behavior:

- Accept `T2T-CHM13v2.0` metadata for same-build restriction.
- Resolve only the configured UCSC FASTA for the declared build.
- Strictly convert declared object labels from `ncbi`, `ucsc`, `plink`, or
  `plink_splitx` to UCSC labels before reference lookup.
- Continue to emit rows in the input object's declared `contig_naming`.
- Fail early for contig labels outside the hard-coded primary-contig contract.
- Do not silently drop rows because a reference asset lacks an expected primary
  contig; missing reference surfaces are config errors.
- Support the `--norm-indels` path for T2T using the configured T2T UCSC FASTA
  with the same `bcftools norm` semantics as GRCh37/GRCh38.

### `liftover_build.py`

Current behavior derives a chain key from a hard-coded two-entry map and always
uses UCSC temporary VCF contigs. New behavior should derive the temporary VCF
surface from the selected liftover edge.

Required behavior:

- Accept `--target-build T2T-CHM13v2.0` and accept
  `T2T-CHM13v2.0` as a source build when the input metadata declares it.
- Resolve a directed edge from `(source_build, target_build)`.
- Support any ordered pair of distinct supported builds when a directed
  liftover edge exists in config: `GRCh37`, `GRCh38`, and
  `T2T-CHM13v2.0`.
- Reject `source_build == target_build` with the usual "liftover is not
  required" error.
- Continue to reject declared `plink_splitx` input. Users must normalize to
  build-independent `plink`, run liftover, and normalize back to `plink_splitx`
  afterward when needed.
- Convert input rows to UCSC labels for the temporary VCF.
- Convert lifted UCSC output back to the input object's declared
  `contig_naming`, preserving the current output-naming contract.
- Retain existing auditable QC statuses for unmapped variants, unsupported
  target contigs, duplicate targets, unsupported non-SNV rows, and
  `ploidy_class_changed`.

Initial implementation should support every configured directed edge among the
three supported builds: GRCh37<->GRCh38, GRCh37<->T2T, and GRCh38<->T2T.

### `prepare_variants.py` and `prepare_variants_sharded.py`

Required behavior:

- Document and validate `--dst-build T2T-CHM13v2.0` as a supported destination.
- Preserve the current default `--dst-build GRCh38`.
- Run `guess_build.py` only when source build metadata is `unknown`.
- Skip `guess_build.py` whenever source build metadata is already known, not
  only during `--resume`.
- T2T source data must carry explicit `genome_build=T2T-CHM13v2.0` metadata.
- Work with T2T as source or destination through the same canonical stage graph
  used for GRCh37/GRCh38.

### `guess_build.py`

`guess_build.py` should remain limited to `GRCh37` and `GRCh38`.

T2T build guessing is out of scope for this proposal. A future T2T classifier
would need a T2T-specific reference or calibration set, such as dbSNP lifted to
T2T, because common SNPs in euchromatic regions may be identical or nearly
identical between GRCh38 and T2T at the classifier's current scoring surface.
Without that calibration, T2T inputs risk being absorbed into the GRCh38 score.

### Contig normalization and ordering

The existing canonical primary-contig model already covers `1-22`, `X`, `Y`,
and `MT` and can represent T2T UCSC labels as `chr1-22`, `chrX`, `chrY`,
`chrM`. The implementation should not add a new public contig-naming mode for
T2T accessions.

Accession labels from the UCSC alias table should be handled, if needed, as
reference alias metadata or a future normalization enhancement. They are out of
scope for the initial public input contract.

### Ploidy model and `plink_splitx`

The current haploid-region schema is build-specific for `GRCh37` and
`GRCh38`. T2T support should rename
`src/genomatch/schemas/human_haploid_regions_grch37_grch38.json` to
`src/genomatch/schemas/human_haploid_regions.json`, update code references to
the new path, update the schema title, and add a `T2T-CHM13v2.0` build entry.

The T2T PAR intervals are:

- PAR1:
  - `chrX:1-2394410`
  - `chrY:1-2458320`
- PAR2:
  - `chrX:153925835-154259566`
  - `chrY:62122810-62460029`

Using the existing schema shape, T2T should define:

- `PAR1` and `PAR2` as diploid in both sexes.
- `X_nonPAR` as all `X` intervals outside PAR1 and PAR2, with male haploid and
  female diploid ploidy.
- `Y_MSY` as all `Y` intervals outside PAR1 and PAR2, with male haploid and
  female absent ploidy.
- `MT` as haploid in both sexes.
- A T2T source entry:
  ```json
  {
    "label": "T2T CHM13 v2.0 analysis set (PAR BED)",
    "url_hint": "marbl/CHM13 GitHub repository (analysis-set downloads, chm13v2.0_PAR.bed)"
  }
  ```

Initial behavior:

- `ncbi`, `ucsc`, and `plink` contig naming can be supported for T2T across
  reference-aware transforms after the haploid schema is extended.
- `normalize_contigs.py --to plink_splitx` should support
  `T2T-CHM13v2.0` metadata by applying the T2T X PAR intervals above.
- `liftover_build.py` should continue rejecting declared `plink_splitx` input.
- `liftover_build.py` and `restrict_build_compatible.py --norm-indels` must
  continue to preserve expected ploidy pairs. For T2T source or target rows on
  `X`, `Y`, or `MT`, this should use the T2T haploid-region definitions above.

## Documentation changes

Update these files when the implementation is made normative:

- `DOWNLOADS.md`: replace the current UCSC-internal reference model with the
  `builds` plus list-valued `liftover` model, add T2T download commands, and
  explain the full-FASTA `chrY` caveat.
- `config.example.yaml`: use the new config shape and include the T2T FASTA
  and four T2T liftover edges.
- `README.md`: add `T2T-CHM13v2.0` to the supported-build summary and wrapper
  option text, while preserving the primary-contig-only scope.
- `SPEC.md`: update the reference-model section after `DOWNLOADS.md` becomes
  normative for the new config.
- `spec/contigs-and-metadata.md`: replace the UCSC-internal note with
  an explicit UCSC-reference requirement for reference-aware tools.
- `spec/variant-transforms.md`: update `restrict_build_compatible.py` and
  `liftover_build.py` resolution rules.
- `spec/ploidy-model.md`: add the T2T PAR and haploid-region intervals above.
- Rename `src/genomatch/schemas/human_haploid_regions_grch37_grch38.json` to
  `src/genomatch/schemas/human_haploid_regions.json`, update code references,
  and add a `T2T-CHM13v2.0` build entry with PAR, X non-PAR, Y MSY, and MT
  regions.

Avoid repeating the full config example across all docs; keep the complete
shape in `DOWNLOADS.md` and link to it elsewhere.

## Implementation plan

1. Add tests for the new config parser and resolver behavior:
   - new `builds` / `liftover` shape
   - legacy shape rejected with a clear config-upgrade error
   - no magic chain-key derivation
   - duplicate or missing liftover edges fail clearly
   - update `write_match_config` in `tests/utils.py` to emit the new
     `builds` plus list-valued `liftover` schema for GRCh37/GRCh38-only test
     configs, with T2T entries omitted by default
2. Refactor `reference_utils.py`:
   - add UCSC FASTA resolution by build
   - add directed edge resolution from list-valued `liftover`
   - keep path resolution relative to `MATCH_CONFIG`
3. Update same-build reference restriction:
   - resolve the UCSC FASTA for the declared build
   - keep converting valid declared input contigs to UCSC labels
   - keep reporting the existing `internal_reference_naming=ucsc` in the JSON
     summary
4. Update liftover:
   - resolve edge assets from `(source_build, target_build)`
   - convert temporary VCF rows to UCSC labels
   - parse lifted contigs as UCSC labels
5. Add hard-coded build support where public CLIs or schemas enumerate builds:
   - `T2T-CHM13v2.0` as a valid source and destination build
   - update `SUPPORTED_GENOME_BUILDS` in `src/genomatch/vtable_utils.py`,
     currently `{"GRCh37", "GRCh38", "unknown"}`
   - no broad arbitrary-build acceptance
6. Keep `guess_build.py` limited to `GRCh37` and `GRCh38`.
7. Update `prepare_variants.py` and `prepare_variants_sharded.py` to skip
   `guess_build.py` whenever source build metadata is already known, not only
   during `--resume`.
8. Rename the haploid-region schema to
   `src/genomatch/schemas/human_haploid_regions.json`, update code references,
   and extend it with the T2T PAR, X non-PAR, Y MSY, and MT definitions.
9. Update docs and examples.
10. Run focused resolver, restrict, liftover, wrapper, ploidy, and split-X
   tests, then the full test suite in the `match-liftover` conda environment
   before committing code.

## Test plan

Unit tests:

- `reference_utils.py` resolves T2T FASTA and T2T liftover edges from the new
  config shape.
- `tests/utils.py::write_match_config` emits the new config shape and omits
  T2T entries by default.
- Unknown build, missing UCSC FASTA, missing chain, duplicate edge, and
  unsupported direction all fail with clear errors.
- Configs without T2T entries continue to work for GRCh37/GRCh38-only
  workflows.
- `restrict_build_compatible.py` fails clearly for T2T input when the config
  lacks a T2T FASTA entry.
- `liftover_build.py` fails clearly when source or target is T2T and the config
  lacks the T2T FASTA or the needed directed T2T chain edge.
- Legacy config shape is rejected with a clear config-upgrade error.
- `SUPPORTED_GENOME_BUILDS` includes `T2T-CHM13v2.0`.
- `guess_build.py` still evaluates only `GRCh37` and `GRCh38`.
- `prepare_variants.py` and `prepare_variants_sharded.py` skip `guess_build.py`
  whenever source build metadata is already known.
- Same-build restriction converts NCBI-style object labels to UCSC labels for
  reference lookup, without resolving an NCBI-style FASTA.
- Same-build T2T restriction uses the UCSC T2T FASTA for UCSC input.
- Same-build T2T restriction supports `--norm-indels`.
- Liftover resolves both upward edges to T2T and downward edges from T2T.
- Liftover rejects identical source and target builds with the existing
  same-build error.
- `expected_ploidy_pair()` returns the expected T2T PAR, X non-PAR, Y MSY, and
  MT ploidy pairs.
- `normalize_contigs.py --to plink_splitx` accepts `T2T-CHM13v2.0` metadata and
  applies the T2T X PAR intervals.
- `prepare_variants.py --dst-build T2T-CHM13v2.0` plans the expected
  restrict/liftover sequence from known GRCh37 or GRCh38 metadata.
- `prepare_variants.py --dst-build GRCh37` or `--dst-build GRCh38` plans the
  expected down-liftover sequence from explicit T2T metadata, without running
  `guess_build.py`.

Integration tests:

- Real `bcftools +liftover` from a small GRCh38 fixture to
  `T2T-CHM13v2.0`.
- Real `bcftools +liftover` from a small GRCh37 fixture to
  `T2T-CHM13v2.0`.
- Real `bcftools +liftover` from a small T2T fixture to `GRCh38`.
- Real `bcftools +liftover` from a small T2T fixture to `GRCh37`.
- T2T-dependent real-reference tests should use a `pytest.skip()` guard when
  the T2T FASTA or T2T chain assets are unavailable in the local reference
  tree.
- `tests/real_liftover_helpers.py` should extend `write_real_match_config` and
  `write_ucsc_subset_fasta` for T2T paths, adding new FASTA constants and new
  chain constants while following the existing GRCh37/GRCh38 pattern exactly.
- A `dbsnp_cleansumstat_reference_T2T_GRCh38.txt` fixture, or an equivalent
  fixture, is needed to drive entry selection for T2T integration cases. It
  should be generated the same way as the existing
  `dbsnp_cleansumstat_reference_GRCh38_GRCh37.txt` fixture.
- A row mapping to an unsupported target contig is audited rather than silently
  omitted.
- Duplicate lifted targets are audited with `duplicate_target`.

Manual validation:

- Confirm `samtools faidx ref/ucsc/T2T-CHM13v2.0/chm13v2.0.fa` exists before
  running reference-aware tools.
- Inspect chain headers to confirm the expected source and target naming
  surfaces are UCSC-style before enabling the edges in default examples.

## Open questions

- Should accession-style T2T labels be accepted by `normalize_contigs.py` in
  the initial release, or deferred until there is a general alias-table
  mechanism?
