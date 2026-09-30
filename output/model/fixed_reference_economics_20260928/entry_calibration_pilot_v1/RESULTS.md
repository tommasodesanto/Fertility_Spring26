# Three entrant-wealth calibration pilots: terminal results

Torch array **18895422** completed all three arms with exit 0:0. All selected points passed two exact full-equilibrium repeats. These are experimental local improvements, not adopted calibrations.

The experiments replace the reference entrant-wealth mapping with three author-approved laws and credit limits: empirical five wealth/income ratios independently times current annual gross entrant earnings with renter floor −0.25; all-zero entry with floor zero; and negative five-bin ratios censored to zero with positive ratios scaled by 0.3632385158888715 to preserve mean, with floor zero. Censoring is applied to the five-bin approximation, not raw survey data. All arms use 120 assets × 9 income states and a common 2% annual interest rate. Earnings, taxes, pension, corrected homeowner rules, mortality repayment, targets and weights are held fixed. Price clears birth renewal, population clears absolute housing supply, and H0 and psi_child are fixed. Nine remaining coordinates are searched against ten scored targets; fourteen total rows are reported.

Baseline comparisons isolate these arm specifications at common parameter guesses. Baseline-to-selected changes measure the limited refitting within each arm. The old 160×15 / D=0.53 comparison has a different contract and is not used as a calibration-gain benchmark.

| Arm | Baseline loss | Selected loss | Reduction | Full GE | Lifecycle solves | Search elapsed |
|---|---:|---:|---:|---:|---:|---:|
| empirical_credit | 22.308 | 21.728 | 2.60% | 15.000 | 133.000 | 42.15 min |
| zero_wealth | 21.709 | 20.584 | 5.18% | 15.000 | 114.000 | 36.87 min |
| nonnegative_mean | 20.682 | 19.697 | 4.76% | 14.000 | 122.000 | 46.29 min |

Each baseline Jacobian has numerical rank nine. Its largest/smallest singular-value ratio is approximately empirical_credit: 63,390, zero_wealth: 42,112, nonnegative_mean: 43,106. Full numerical rank does not establish strong identification; no fresh selected-point Jacobian was computed. Each arm performed one derivative/Gauss–Newton round. The nonnegative arm preserved the final-repeat time reserve by omitting the second damped proposal. This is a time-bounded pilot, with no search convergence or fine-grid accuracy claim.

The main residual fit problems are early fertility, wealth/earnings and the first-birth housing response. Complete economic comparisons follow; the linked CSVs retain full precision and all weights, gaps and contributions.

| Moment | Target | Empirical selected | Zero selected | Nonnegative selected |
|---|---:|---:|---:|---:|
| initial_normalization | 2.100 | 2.100 | 2.100 | 2.100 |
| cps_childlessness | 0.198 | 0.201 | 0.201 | 0.202 |
| cps_exactly_one | 0.214 | 0.209 | 0.210 | 0.209 |
| nchs_mean_age | 25.976 | 25.933 | 26.021 | 25.963 |
| nchs_share30 | 0.249 | 0.222 | 0.228 | 0.225 |
| wealth_earnings | 6.927 | 6.262 | 6.231 | 6.308 |
| bequest_wealth | 0.007 | 0.007 | 0.007 | 0.007 |
| old_dispersion | 3.516 | 3.043 | 3.019 | 3.038 |
| mean_rooms | 5.729 | 5.836 | 5.846 | 5.825 |
| ownership_30_55 | 0.676 | 0.649 | 0.666 | 0.663 |
| first_birth_rooms | 1.465 | 1.632 | 1.623 | 1.634 |
| family_rooms | 0.385 | 0.354 | 0.358 | 0.361 |
| recent_parent_ownership | 0.128 | 0.118 | 0.118 | 0.118 |
| early_fertility | 0.810 | 0.535 | 0.527 | 0.533 |

## empirical_credit

Baseline price 0.799826477; selected price 0.795437316. Baseline population 1.013404037; selected population 1.006481161.

### Baseline

[Complete 14-target CSV](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/empirical_credit/000_baseline/phase_b_ge/selected_root/target_fit.csv) · [All 31 parameters, restrictions, bounds and bound flags](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/empirical_credit/000_baseline/phase_b_ge/selected_root/parameters.csv) · [Standard 17-plot folder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/empirical_credit/000_baseline/phase_b_ge/selected_root/standard_diagnostics)

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
|---|---:|---:|---:|---:|---:|---|
| initial_normalization | 2.100 | 2.100 | -1.042e-10 | — | — | normalization |
| cps_childlessness | 0.198 | 0.201 | 0.003 | 35532.304 | 0.329 | scored |
| cps_exactly_one | 0.214 | 0.209 | -0.005 | 26952.821 | 0.649 | scored |
| nchs_mean_age | 25.976 | 25.931 | -0.046 | 139.828 | 0.292 | scored |
| nchs_share30 | 0.249 | 0.222 | -0.028 | 0.000 | 0.000 | validation |
| wealth_earnings | 6.927 | 6.248 | -0.679 | 7.595 | 3.499 | scored |
| bequest_wealth | 0.007 | 0.007 | -1.998e-04 | 5165289.256 | 0.206 | scored |
| old_dispersion | 3.516 | 3.037 | -0.479 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.816 | 0.087 | 128.021 | 0.966 | scored |
| ownership_30_55 | 0.676 | 0.642 | -0.034 | 2339.362 | 2.743 | scored |
| first_birth_rooms | 1.465 | 1.632 | 0.167 | 137.565 | 3.850 | scored |
| family_rooms | 0.385 | 0.356 | -0.029 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.118 | -0.009 | 27055.823 | 2.248 | scored |
| early_fertility | 0.810 | 0.535 | -0.274 | 100.000 | 7.525 | scored |
### Selected

[Complete 14-target CSV](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/empirical_credit/011_gn_0.5/phase_b_ge/selected_root/target_fit.csv) · [All 31 parameters, restrictions, bounds and bound flags](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/empirical_credit/011_gn_0.5/phase_b_ge/selected_root/parameters.csv) · [Standard 17-plot folder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/empirical_credit/011_gn_0.5/phase_b_ge/selected_root/standard_diagnostics)

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
|---|---:|---:|---:|---:|---:|---|
| initial_normalization | 2.100 | 2.100 | -1.923e-11 | — | — | normalization |
| cps_childlessness | 0.198 | 0.201 | 0.003 | 35532.304 | 0.331 | scored |
| cps_exactly_one | 0.214 | 0.209 | -0.005 | 26952.821 | 0.658 | scored |
| nchs_mean_age | 25.976 | 25.933 | -0.044 | 139.828 | 0.266 | scored |
| nchs_share30 | 0.249 | 0.222 | -0.027 | 0.000 | 0.000 | validation |
| wealth_earnings | 6.927 | 6.262 | -0.665 | 7.595 | 3.359 | scored |
| bequest_wealth | 0.007 | 0.007 | -1.736e-04 | 5165289.256 | 0.156 | scored |
| old_dispersion | 3.516 | 3.043 | -0.473 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.836 | 0.107 | 128.021 | 1.455 | scored |
| ownership_30_55 | 0.676 | 0.649 | -0.027 | 2339.362 | 1.679 | scored |
| first_birth_rooms | 1.465 | 1.632 | 0.167 | 137.565 | 3.817 | scored |
| family_rooms | 0.385 | 0.354 | -0.031 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.118 | -0.010 | 27055.823 | 2.476 | scored |
| early_fertility | 0.810 | 0.535 | -0.274 | 100.000 | 7.532 | scored |

Both selected fertility dispersion parameters, `kappa_fert` and `kappa_fert_continuation`, are flagged near their lower bounds in this arm, as in all three arms. The flags measure distance relative to the wide [0.020, 50.000] search interval; exact estimates and bounds are in the linked parameter CSV.

## zero_wealth

Baseline price 0.786114877; selected price 0.789671405. Baseline population 0.994318684; selected population 1.000244437.

### Baseline

[Complete 14-target CSV](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/zero_wealth/000_baseline/phase_b_ge/selected_root/target_fit.csv) · [All 31 parameters, restrictions, bounds and bound flags](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/zero_wealth/000_baseline/phase_b_ge/selected_root/parameters.csv) · [Standard 17-plot folder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/zero_wealth/000_baseline/phase_b_ge/selected_root/standard_diagnostics)

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
|---|---:|---:|---:|---:|---:|---|
| initial_normalization | 2.100 | 2.100 | -5.526e-10 | — | — | normalization |
| cps_childlessness | 0.198 | 0.200 | 0.002 | 35532.304 | 0.127 | scored |
| cps_exactly_one | 0.214 | 0.211 | -0.003 | 26952.821 | 0.261 | scored |
| nchs_mean_age | 25.976 | 26.006 | 0.029 | 139.828 | 0.121 | scored |
| nchs_share30 | 0.249 | 0.227 | -0.022 | 0.000 | 0.000 | validation |
| wealth_earnings | 6.927 | 6.220 | -0.706 | 7.595 | 3.789 | scored |
| bequest_wealth | 0.007 | 0.007 | -1.565e-04 | 5165289.256 | 0.127 | scored |
| old_dispersion | 3.516 | 3.030 | -0.486 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.864 | 0.134 | 128.021 | 2.309 | scored |
| ownership_30_55 | 0.676 | 0.661 | -0.015 | 2339.362 | 0.557 | scored |
| first_birth_rooms | 1.465 | 1.626 | 0.161 | 137.565 | 3.567 | scored |
| family_rooms | 0.385 | 0.351 | -0.034 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.117 | -0.010 | 27055.823 | 2.969 | scored |
| early_fertility | 0.810 | 0.529 | -0.281 | 100.000 | 7.884 | scored |
### Selected

[Complete 14-target CSV](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/zero_wealth/011_gn_0.5/phase_b_ge/selected_root/target_fit.csv) · [All 31 parameters, restrictions, bounds and bound flags](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/zero_wealth/011_gn_0.5/phase_b_ge/selected_root/parameters.csv) · [Standard 17-plot folder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/zero_wealth/011_gn_0.5/phase_b_ge/selected_root/standard_diagnostics)

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
|---|---:|---:|---:|---:|---:|---|
| initial_normalization | 2.100 | 2.100 | -3.452e-07 | — | — | normalization |
| cps_childlessness | 0.198 | 0.201 | 0.003 | 35532.304 | 0.254 | scored |
| cps_exactly_one | 0.214 | 0.210 | -0.004 | 26952.821 | 0.354 | scored |
| nchs_mean_age | 25.976 | 26.021 | 0.044 | 139.828 | 0.275 | scored |
| nchs_share30 | 0.249 | 0.228 | -0.021 | 0.000 | 0.000 | validation |
| wealth_earnings | 6.927 | 6.231 | -0.695 | 7.595 | 3.671 | scored |
| bequest_wealth | 0.007 | 0.007 | -1.297e-04 | 5165289.256 | 0.087 | scored |
| old_dispersion | 3.516 | 3.019 | -0.496 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.846 | 0.116 | 128.021 | 1.727 | scored |
| ownership_30_55 | 0.676 | 0.666 | -0.010 | 2339.362 | 0.251 | scored |
| first_birth_rooms | 1.465 | 1.623 | 0.158 | 137.565 | 3.446 | scored |
| family_rooms | 0.385 | 0.358 | -0.027 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.118 | -0.010 | 27055.823 | 2.551 | scored |
| early_fertility | 0.810 | 0.527 | -0.282 | 100.000 | 7.968 | scored |

Both selected fertility dispersion parameters, `kappa_fert` and `kappa_fert_continuation`, are flagged near their lower bounds in this arm, as in all three arms. The flags measure distance relative to the wide [0.020, 50.000] search interval; exact estimates and bounds are in the linked parameter CSV.

## nonnegative_mean

Baseline price 0.797381848; selected price 0.799093457. Baseline population 1.008448108; selected population 1.011339610.

### Baseline

[Complete 14-target CSV](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/nonnegative_mean/000_baseline/phase_b_ge/selected_root/target_fit.csv) · [All 31 parameters, restrictions, bounds and bound flags](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/nonnegative_mean/000_baseline/phase_b_ge/selected_root/parameters.csv) · [Standard 17-plot folder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/nonnegative_mean/000_baseline/phase_b_ge/selected_root/standard_diagnostics)

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
|---|---:|---:|---:|---:|---:|---|
| initial_normalization | 2.100 | 2.100 | -1.045e-09 | — | — | normalization |
| cps_childlessness | 0.198 | 0.201 | 0.003 | 35532.304 | 0.273 | scored |
| cps_exactly_one | 0.214 | 0.209 | -0.004 | 26952.821 | 0.474 | scored |
| nchs_mean_age | 25.976 | 25.950 | -0.026 | 139.828 | 0.098 | scored |
| nchs_share30 | 0.249 | 0.225 | -0.025 | 0.000 | 0.000 | validation |
| wealth_earnings | 6.927 | 6.301 | -0.625 | 7.595 | 2.970 | scored |
| bequest_wealth | 0.007 | 0.007 | -2.360e-04 | 5165289.256 | 0.288 | scored |
| old_dispersion | 3.516 | 3.040 | -0.476 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.834 | 0.104 | 128.021 | 1.390 | scored |
| ownership_30_55 | 0.676 | 0.661 | -0.015 | 2339.362 | 0.549 | scored |
| first_birth_rooms | 1.465 | 1.635 | 0.170 | 137.565 | 3.984 | scored |
| family_rooms | 0.385 | 0.355 | -0.030 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.117 | -0.011 | 27055.823 | 3.059 | scored |
| early_fertility | 0.810 | 0.534 | -0.276 | 100.000 | 7.599 | scored |
### Selected

[Complete 14-target CSV](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/nonnegative_mean/011_gn_0.5/phase_b_ge/selected_root/target_fit.csv) · [All 31 parameters, restrictions, bounds and bound flags](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/nonnegative_mean/011_gn_0.5/phase_b_ge/selected_root/parameters.csv) · [Standard 17-plot folder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/nonnegative_mean/011_gn_0.5/phase_b_ge/selected_root/standard_diagnostics)

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
|---|---:|---:|---:|---:|---:|---|
| initial_normalization | 2.100 | 2.100 | -1.662e-09 | — | — | normalization |
| cps_childlessness | 0.198 | 0.202 | 0.004 | 35532.304 | 0.445 | scored |
| cps_exactly_one | 0.214 | 0.209 | -0.005 | 26952.821 | 0.593 | scored |
| nchs_mean_age | 25.976 | 25.963 | -0.013 | 139.828 | 0.023 | scored |
| nchs_share30 | 0.249 | 0.225 | -0.024 | 0.000 | 0.000 | validation |
| wealth_earnings | 6.927 | 6.308 | -0.619 | 7.595 | 2.909 | scored |
| bequest_wealth | 0.007 | 0.007 | -2.117e-04 | 5165289.256 | 0.232 | scored |
| old_dispersion | 3.516 | 3.038 | -0.478 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.825 | 0.095 | 128.021 | 1.165 | scored |
| ownership_30_55 | 0.676 | 0.663 | -0.013 | 2339.362 | 0.384 | scored |
| first_birth_rooms | 1.465 | 1.634 | 0.169 | 137.565 | 3.951 | scored |
| family_rooms | 0.385 | 0.361 | -0.024 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.118 | -0.009 | 27055.823 | 2.323 | scored |
| early_fertility | 0.810 | 0.533 | -0.277 | 100.000 | 7.672 | scored |

Both selected fertility dispersion parameters, `kappa_fert` and `kappa_fert_continuation`, are flagged near their lower bounds in this arm, as in all three arms. The flags measure distance relative to the wide [0.020, 50.000] search interval; exact estimates and bounds are in the linked parameter CSV.

## Verification and limitations

All 228 collected remote files have independently matching SHA-256 hashes; baseline and selected each contain the standard 17 PNGs. Each baseline and selected table has 14 targets and 31 parameters. The baseline repeat and both selected repeats reproduce target/parameter CSV and closure JSON bytes exactly. All baseline and selected plot hashes match their repeat receipts. All three input target fingerprints match the pinned plan. The frozen checkpoint/authenticator source receipts match across arms. Runtime source authentication was enforced by the pinned launcher; this collection does not claim a fresh source replay. No large model arrays were downloaded.

Native estate-accounting caveats remain: availability valuation is provisional, liabilities have an unresolved creditor rule, and the ledger does not certify counterparty or physical housing settlement. Numerical/ledger admissibility does not resolve those economic limitations. There has been no production adoption or transition validation.

[Collection verification](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/final_verification.json) · [Remote download hashes](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/remote_file_sha256.json)
