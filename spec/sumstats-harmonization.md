# Column Specification

This document defines summary-statistic harmonization for the clean summary-stat application path.

It implements a module for `apply_vmap_to_sumstats.py` which applies after output rows are matched to the target side of a `.vmap`,
but before effect-direction normalization is applied based on `allele_op`.
Row retention and `--only-mapped-target` behavior are defined by [payload-application.md](payload-application.md), not by this module.
This module is a missing-value and missing-column completion pipeline. It does not attempt to verify semantic consistency between overlapping fields that are already present, for example `BETA` vs `OR`, `EAF` vs `OAF`, `P` vs `Z`, or row-level `N` vs metadata totals. Existing non-missing values are retained unless an explicit transform or range rule says otherwise.

This specification describes a pipeline, with each step updating the current stage in place for the purpose of the next step.
Each step accepts a dataframe (`sumstats`) and `metadata` describing its columns, and returns modified `sumstats`+`metadata`.

This module never drops or re-orders lines; each row of the output `sumstats` dataframe corresponds to the row of the input.
Rows containing only missing values are valid input and valid output for this module.
If the caller uses `--only-mapped-target`, any later dropping of rows with missing `P` is outside this module and belongs to the caller contract in [payload-application.md](payload-application.md).

Harmonization pipeline is a fixed sequence of transformation steps. CLI parsing, file I/O, and wrapper orchestration are described in [payload-application.md](payload-application.md).

Missing values may use any convenient internal representation.

User options for this operation:

- `fill-mode`: `column` or `row`; default: `column`
- `use-af-inference`: boolean

The implementation of this module must be contained in a separate file or module from the base non-clean apply path so the added harmonization logic does not pollute ordinary `apply_vmap_to_sumstats.py` behavior.

### Metadata

Metadata is described by `schemas/raw-sumstats-metadata.yaml`.
The cleaned output payload/stat column order is declared by `schemas/cleaned-sumstats.yaml`, excluding `CHR`, `POS`, `SNP`, `EffectAllele`, and `OtherAllele` which are outside this module.
`schemas/cleaned-sumstats.yaml` is guidance for canonical column naming and preferred ordering; it is not a normative validation contract for the output of this module.

This module does not use, drop, or otherwise process source-side `CHR`, `POS`, `SNP`, `EffectAllele`, `OtherAllele`, or joined source variant fields. Those are handled by base payload application.

`metadata` may define any subset of the following optional payload/stat columns:

- `col_BETA`
- `col_SE`
- `col_Z`
- `col_P`
- `col_OR`
- `col_ORL95`
- `col_ORU95`
- `col_N`
- `col_CaseN`
- `col_ControlN`
- `col_StudyN`
- `col_INFO`
- `col_Direction`
- `col_EAF`
- `col_CaseEAF`
- `col_ControlEAF`
- `col_OAF`
- `col_CaseOAF`
- `col_ControlOAF`

Only `col_*` keys actually present in `metadata` are required to resolve uniquely. Missing metadata entries are allowed; such columns may later be derived by the pipeline.
Rules below refer to column names without `col_` prefix, i.e. `BETA` refers to the column defined by `col_BETA` metadata when that metadata entry is present.

### metadata check and normalization

- normalize `metadata`'s model labels before stat selection:
  - `linear mixed-model` becomes `linear`
  - `logistic mixed-model` becomes `logistic`
- if the normalized model is `ordinal`, `cox`, or `other`, fail with `not implemented`
- if `stats_TotalN` is present and model is `logistic`, emit a warning and ignore `stats_TotalN`

### sumstats and metadata pre-processing step

- sumstat headers: remove every character not in `[A-Za-z0-9]`
- sumstat headers: lowercase
- sumstat headers: trim surrounding whitespace
- metadata string values used as column references: apply the same cleanup so `col_*` references match the cleaned headers
- arbitrary unrecognized columns may remain ambiguous after header normalization; they are later dropped
- raise error if any metadata-declared `col_*` reference cannot be identified uniquely in the `sumstats` dataframe after the full normalization above
- rename `sumstats` columns to canonical names as per `schemas/cleaned-sumstats.yaml`, and update `metadata` accordingly.
- if two metadata-declared input columns collapse to the same canonical cleaned column name after normalization and renaming, fail cleanly rather than choosing one
- drop any columns that are unrecognized (not described by the metadata)
- drop columns with all missing values, and remove them from metadata

### P-value normalization

Apply these rules to `col_P`:

- treat non-numeric `P` as missing
- if both `stats_neglog10P` and `stats_log10P` are true, raise an error
- if either `stats_neglog10P` or `stats_log10P` is true and `col_P` is missing, raise an error
- if `stats_neglog10P: true`, treat `col_P` as `-log10(P)` and convert it to plain `P`
- if `stats_log10P: true`, treat `col_P` as `log10(P)` and convert it to plain `P`
- if neither `stats_neglog10P` nor `stats_log10P` is true, keep `col_P` as already being plain `P`

### Restrict to numerically valid and within range values

- treat all non-numeric or non-finite values as missing, except for `Direction`
- treat `SE <= 0` (negative or zero) as missing
- treat `P < 0` or `P > 1` as missing
- treat `OR <= 0`, `ORL95 <= 0`, and `ORU95 <= 0` as missing
- treat `N < 0`, `CaseN < 0`, `ControlN < 0` as missing
- treat `INFO < 0` as missing
- treat `AF < 0` or `AF > 1` as missing for all allele frequency columns (`EAF`, `CaseEAF`, `ControlEAF`, and their `OAF` counterparts)
- `Direction` is passed through unchanged

### Effective AF definition from Other AF

- if `col_EAF` is missing, populate from `1-col_OAF`
- if `col_CaseEAF` is missing, populate from `1-col_CaseOAF`
- if `col_ControlEAF` is missing, populate from `1-col_ControlOAF`
- drop `OAF`, `CaseOAF`, `ControlOAF` if present, and remove from metadata

### `x_from_y` deriving rules

Apply these additional rules sequentially.

In `fill-mode=row`: each rule is valid only where the output field is missing for a given row and all required inputs are present for that row.
In `fill-mode=column`: each rule is valid only when the whole output column is absent from the dataframe and all input columns are present. A present column with some missing values is not partially filled in this mode.
An all-missing column that was dropped earlier counts as absent for later `fill-mode=column` derivation.

- `N`: if missing, use `stats_TotalN` ; apply only if model is `linear`
- `CaseN`: if missing, use `stats_CaseN` ; apply only if model is `logistic`
- `ControlN`: if missing, use `stats_ControlN` ; apply only if model is `logistic`
- `N`: if missing, use `Neff_from_Nca_and_Nco(CaseN, ControlN)` ; apply only if model is `logistic`.
- `EAF`: if missing, use `AF_from_CaseAF_ControlAF(CaseEAF, ControlEAF, CaseN, ControlN)` ; apply only if model is `logistic`.
- `BETA`: if missing, use `beta_from_oddsratio(OR)`; apply only if model is `logistic`.
- `SE`: if missing, use `se_from_ORu95_ORl95(ORU95, ORL95)`; apply only if model is `logistic`.

- `Z`: if missing, use `zscore_from_pval_beta(P, BETA)`; prefer this rule over BETA/SE rule, because in GWAS p-value is more robust
- `Z`: if missing, use `zscore_from_beta_se(BETA, SE)`
- `P`: if missing, use `pval_from_zscore(Z)`
- `SE`: if missing, use `se_from_zscore_beta(Z, BETA)`
- `BETA`: if missing, use `beta_from_zscore_se(Z, SE)`

If `use-af-inference` is enabled:
- `N`: if missing, use `N_from_zscore_beta_af(Z, BETA, EAF)`
- `BETA`: if missing, use `beta_from_zscore_N_af(Z, N, EAF)`
- `SE`: if missing, use `se_from_zscore_N_af(Z, N, EAF)`

After derivation, apply a second validation pass using the same non-numeric, non-finite, and range rules as above.

Historical note: the following rules from cleansumstats will not be implemented:
- `Z`: if missing, use `zscore_from_pval_beta_N(P, BETA, N)`
- `P`: if missing, use `pval_from_zscore_N(Z, N)`

### Output

Re-order payload/stat columns according to `schemas/cleaned-sumstats.yaml`, excluding `CHR`, `POS`, `SNP`, `EffectAllele`, and `OtherAllele` which are outside this module.
Return updated `sumstats` dataframe and updated `metadata`.

## Rule definitions

### Mathematical notation

- `abs(x)`: absolute value of `x`
- `sign(x)`: sign of `x`, taking values `-1`, `0`, or `1`
- `log(x)`: natural logarithm of `x`
- `pnorm(x)`: cumulative distribution function of the standard normal distribution
- `qnorm(p)`: quantile function of the standard normal distribution
- `pt(x, df)`: cumulative distribution function of the Student t distribution with `df` degrees of freedom
- `qt(p, df)`: quantile function of the Student t distribution with `df` degrees of freedom

### Formulas

If a formula encounters invalid math or produces a non-finite value, the output is missing.

- `pval_from_neglog10p(x) = 10^(-x)`
- `pval_from_log10p(x) = 10^(x)`
- `zscore_from_beta_se(BETA, SE) = BETA / SE`
- `zscore_from_pval_beta(P, BETA) = sign(BETA) * abs(qnorm(P / 2))`
- `zscore_from_pval_beta_N(P, BETA, N) = sign(BETA) * abs(qt(P / 2, N - 2))`
- `pval_from_zscore_N(Z, N) = 2 * pt(-abs(Z), N - 2)`
- `pval_from_zscore(Z) = 2 * pnorm(-abs(Z))`
- `beta_from_zscore_se(Z, SE) = Z * SE`
- `beta_from_zscore_N_af(Z, N, AF) = Z / sqrt(2 * AF * (1 - AF) * (N + Z^2))`
- `se_from_zscore_beta(Z, BETA) = BETA / Z`
- `se_from_zscore_N_af(Z, N, AF) = 1 / sqrt(2 * AF * (1 - AF) * (N + Z^2))`
- `N_from_zscore_beta_af(Z, BETA, AF) = Z^2 / (2 * AF * (1 - AF) * BETA^2) - Z^2`
- `beta_from_oddsratio(OR) = log(OR)`
- `se_from_ORu95_ORl95(ORu95, ORl95) = (log(ORu95) - log(ORl95)) / (2 * qnorm(0.975))`
- `Neff_from_Nca_and_Nco(Nca, Nco) = 4 / ((1 / Nca) + (1 / Nco))`
- `AF_from_CaseAF_ControlAF(CaseAF, ControlAF, Nca, Nco) = (CaseAF * Nca + ControlAF * Nco) / (Nca + Nco)`
