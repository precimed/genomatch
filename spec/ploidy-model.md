# Ploidy model

This file defines the expected ploidy model used by coordinate-changing variant transforms and genotype-payload application.

## Scope

- The expected ploidy model is defined for ordinary human diploid samples with recorded sex.
- Aneuploidy, mosaicism, sample-specific structural ploidy changes, and other abnormal ploidy states are out of scope.
- The model is used only for target-side variant rows on primary reference contigs.

## Expected ploidy

- Expected ploidy is derived from target `chrom`, target `pos`, target `genome_build`, and sample sex.
- For payload application, sample sex is taken from the resolved output sample axis, whether that comes from propagated source `.fam` / `.psam` metadata or an explicit target sample file.
- The authoritative source of build-specific haploid-region definitions is [human_haploid_regions_grch37_grch38.json](../schemas/human_haploid_regions_grch37_grch38.json).
- The model must not depend on literal contig-naming strings such as `ncbi`, `ucsc`, `plink`, or `plink_splitx`; tools must normalize contig interpretation internally before applying ploidy logic.

For one retained row, the canonical row-level ploidy descriptor is the expected ploidy pair:

- `(male_ploidy, female_ploidy)`

Allowed ploidy values are:

- `2` for diploid
- `1` for haploid
- `0` for absent

The pair `(male_ploidy, female_ploidy)` is the normative object used for row-level comparison by coordinate-changing tools. This avoids ambiguity in sex-dependent regions such as X non-PAR and Y.

Terminology:

- `diploid in both sexes` means the expected ploidy pair is `(2, 2)`
- examples of rows diploid in both sexes include autosomes `1-22` and X PAR rows
- `non-diploid in at least one sex` means the expected ploidy pair is not `(2, 2)`
- examples of rows non-diploid in at least one sex include X non-PAR rows, Y rows, and MT rows

## Coordinate-changing transforms

Tools that can change a row's effective biological ploidy by changing `chrom` or `pos` must preserve the expected ploidy pair.

Normative rule:

- If a tool changes `chrom` or `pos` in a way that changes the expected `(male_ploidy, female_ploidy)` pair, that row must not be emitted.
- Such rows must be dropped with auditable QC status `ploidy_class_changed`.
- This comparison is applied to the final candidate row that would otherwise be emitted by the tool, after all coordinate-changing transform logic for that row is complete and immediately before final output/QC accounting.

This currently applies to:

- `liftover_build.py`
- `restrict_build_compatible.py --norm-indels` for normalization-aware branches that can change `chrom` or `pos`

This does not apply to tools that only relabel contigs or reorder rows without changing biological loci, such as:

- `normalize_contigs.py`
- `sort_variants.py`
- matching, intersection, union, and materialization tools

## Payload application

`apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py` share the same target-side ploidy model.

Normative rule:

- Payload application must validate emitted genotype content against expected target-side ploidy.
- This validation is secondary QC only. It must not redefine the upstream contract that coordinate-changing transforms preserve ploidy class.
- Validation warnings must be auditable and aggregated; payload-application tools must not emit one warning per offending genotype.
- Payload-application tools must not rewrite offending genotype content to missing as part of ploidy validation.

## Unknown sex handling

- Ploidy validation is performed subject by subject whenever that subject's sex is known.
- Unknown-sex subjects do not disable validation for other subjects.
- For sex-independent regions, validation still applies even when some subjects have unknown sex.
- For sex-dependent regions, subjects with unknown sex are skipped by ploidy validation and counted in aggregate warning/reporting as unvalidated if reporting is implemented.

## BFILE validation and `.ploidy`

For BFILE payload application:

- BED does not encode source ploidy explicitly; validation operates on observed emitted hardcall states.
- Observed hardcall states are `hom-ref (0)`, `het (1)`, `hom-alt (2)`, and missing.
- If target ploidy is diploid, observed `hom-ref (0)`, `het (1)`, `hom-alt (2)`, and missing are accepted.
- If target ploidy is haploid, observed `het (1)` is incompatible.
- If target ploidy is haploid, observed `hom-ref (0)`, `hom-alt (2)`, and missing are accepted.
- If target ploidy is absent, any nonmissing observed hardcall state is incompatible.
- If target ploidy is absent, only missing is accepted.
- Compatible and incompatible hardcall states are validated after sample-axis reconciliation and after any `allele_op`-driven swap.

`.ploidy` contract:

- `apply_vmap_to_bfile.py` emits `.ploidy` whenever retained target rows include any row that is non-diploid in at least one sex.
- `.ploidy` is a per-variant expected ploidy summary, not a per-sample validation report.
- Each `.ploidy` row stores exactly `(male_ploidy, female_ploidy)` for the corresponding emitted output row.
- `.ploidy` emission is independent of whether some samples have unknown sex.
- When the output prefix contains `@`, `.ploidy` emission is determined per emitted output shard: a shard emits `.ploidy` iff that shard includes at least one row that is non-diploid in at least one sex; shards containing only `(2, 2)` rows do not emit `.ploidy`.

## PFILE validation

For PFILE payload application:

- Validation uses the same expected target-side ploidy model as BFILE.
- Hardcall validation is based on emitted allele observations after sample-axis reconciliation and after any `allele_op`-driven swap.
- Hardcall phase flags do not affect ploidy compatibility.
- For validation, the emitted per-sample allele pair is classified as one of:
- fully missing: `(-9, -9)`
- haploid-like call: exactly one allele is `0` or `1`, and the other is `-9`
- diploid hom-ref: `(0, 0)`
- diploid het: `(0, 1)` or `(1, 0)`
- diploid hom-alt: `(1, 1)`
- malformed partial state: any other partially missing or otherwise malformed pair
- If target ploidy is diploid, fully missing, haploid-like, diploid hom-ref, diploid het, and diploid hom-alt are accepted; malformed partial states are incompatible.
- If target ploidy is haploid, fully missing, haploid-like, diploid hom-ref, and diploid hom-alt are accepted; diploid het and malformed partial states are incompatible.
- If target ploidy is absent, only fully missing is accepted; all other states are incompatible.
- If target ploidy is absent, any nonmissing unphased dosage is incompatible.

Dosage convention:

- For payload-application purposes, PFILE unphased dosages are treated as allele-count values on a `0..2` scale, including haploid chromosomes.
- Analysis-tool conventions that rescale haploid dosages to `0..1` are out of scope for payload-application semantics.
- Under this convention, `swap` does not require ploidy-aware dosage rewriting; the existing `2 - dosage` rule applies uniformly.

Unsupported content remains unchanged:

- Phased dosage preservation remains unsupported.

## Swap semantics

`swap` and `flip_swap` are payload-encoding transforms, not ploidy transforms.

- For BFILE hardcalls, swap semantics do not depend on ploidy.
- For PFILE hardcalls, swap semantics do not depend on ploidy.
- For PFILE hardcall phase flags, swap semantics do not depend on ploidy and phase flags are preserved unchanged.
- For PFILE unphased dosages, swap semantics follow the payload-application dosage convention above and do not branch on expected ploidy.
