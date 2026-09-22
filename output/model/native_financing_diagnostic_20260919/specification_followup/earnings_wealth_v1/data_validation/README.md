# Direct four-year gross-household-earnings validation

Diagnostic only; no target or model specification is adopted.

Input: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/psid_income_md_extract.dta`; SHA-256 `147da916b58700ceb5dc1846921c7c6983e467d926bf946c7a52e85e9ce0c2a5`.
The input has 48,329 positive annual person-years for 6,769 persons after the existing extract filters.

## Definition

A block is the arithmetic sum of four observed consecutive annual `EARNINDRRC` values; the logged block mean is used for covariances because dividing every block by four changes only the log intercept. Only annual-era years 1984--1997 are eligible. Post-1997 biennial observations are never interpolated or treated as zeros.

For each offset `o = 0,1,2,3`, starts are `1984 + o + 4k`; therefore blocks do not overlap within an offset. Block weights are the geometric mean of the four annual `IW` values. Log block means are residualized by weighted least squares on block-start integer age and block-start calendar-year fixed effects. Covariances use geometric-mean block weights and only adjacent non-overlapping blocks at four-year lags 1 and 2.

The three moments `(gamma_0, gamma_1, gamma_2)` identify the three nonnegative parameters `(persistent variance, rho, iid variance)` absent sampling noise. The optimizer imposes `rho in [1e-8, 1-1e-8]` and both variances `>= 0`; unconstrained implied iid variance and boundary hits are reported rather than clipped silently.

## Point estimates and support

| offset | persons | blocks | lag-1 pairs | lag-2 pairs | rho(4yr) | persistent var | iid var | iid boundary | unconstrained iid var |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|
| 0 | 3932 | 7812 | 3836 | 1449 | 0.776145 | 0.481263 | 0.0272219 | False | 0.0272219 |
| 1 | 3982 | 7957 | 3929 | 1491 | 0.797071 | 0.47871 | 0.0423288 | False | 0.0423288 |
| 2 | 3911 | 7339 | 3389 | 1150 | 0.813485 | 0.465589 | 0.0473337 | False | 0.0473337 |
| 3 | 3294 | 5190 | 1896 | 0 | nan | nan | nan |  | nan |

## Frequency and concept comparison

The existing no-fixed-type annual candidate reports `rho_annual = 0.970303`, persistent variance `0.692827`, and transitory variance `0.344442`. The direct estimates above are four-year-frequency covariances from complete annual cells, so their `rho_4yr` and variances are not numerically comparable without an explicit aggregation map. The annual candidate also uses the same EARNINDRRC concept but a different residualization/moment schedule and an endpoint-plus-iid period approximation; this packet does not force annual AR(1) equivalence.

The existing fixed-effect annual packet's fitted annual `rho` is 0.886345, with fixed-effect variance 0.393053, persistent variance 0.331897, and transitory variance 0.309752. That fixed-effect decomposition is a separate annual observer and is not imposed in this direct four-year no-fixed-type fit.

All lag-0/1/2 empirical and fitted covariances, pair counts, and pair-person counts are in `block_covariances.csv`; the full optimizer receipt is in `fit.csv`. Existing annual fixed-effect and no-fixed candidate source values are preserved in `existing_annual_comparison.json` for concept-mismatch review.

`residual_variance_by_age_bin.csv` reports the weighted residual log variance by block-start age bins 25--34, 35--44, 45--54, and 55--60. The source extract contains no ages 18--24, so it cannot validate the model's entrant-age distribution; that is an external entry restriction rather than evidence of zero entrant risk.

The Stata variable labels identify `year` as survey year and `EARNINDRRC` as tax-year earnings. The diagnostic retains the survey-year labels and applies no unverified timing shift; a common shift would leave within-person covariance unchanged, but the survey/tax-year distinction limits level and age-profile interpretation.

`synthetic_observer_recovery.json` applies the same block construction, missingness shape, age/year FE residualization, weights, and covariance estimator to a known four-year AR(1)+iid process. Its finite-sample bias is reported rather than treated as a pass/fail calibration result. `block_construction_validation.json` confirms no overlapping blocks within any offset.

The deterministic bootstrap-weight test is recorded in `run_metadata.json`; it verifies that omitted persons receive zero frequency rather than the point-estimate default weight.

## Interpretation limits

The four-year process is estimated on complete annual cells only, so support is short and selected toward survivors with positive reported labor earnings. It is not an estimate of the post-1997 biennial population. `EARNINDRRC` is gross labor earnings, while the model budget is after-payroll-tax period household resources; the model's tax, pension, transfers, and entry mapping remain separate objects. This direct process must therefore be compared to the existing annual source as a measurement diagnostic, not substituted mechanically.

Person bootstrap draws requested: 199. Bootstrap files are present only when the run is invoked with a positive `--bootstrap-reps`.
