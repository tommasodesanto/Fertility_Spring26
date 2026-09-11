# Initial fit at the mapped starting point

**This is a complete diagnostic readout, not calibrated SMM.** The parenthood-only utility and balanced initial pensions have been evaluated at the mapped starting structural coordinates; those nine coordinates have not been re-estimated. One initial preference level is separately normalized to 2.1. No complete early weight profile is active, so every actual weight and loss contribution below is unavailable and no total loss is calculated.

The saved initial-loop receipt reports two repetitions and eight stationary solves. The fresh observation pass reads the same checkpoint (4afc7fc6f4db32a1bb220bb3b2b30c6e2c0b96822a152c294f6c73c228320852) and performs no new household or equilibrium solve. Initial market residual is 4.3850071e-06; the absolute scaled pension imbalance is 2.77249189e-10, below its 1e-6 gate. Recorded household budget-violating mass and occupied negative value steps are zero. These checks concern the supplied initial calculation: no perfect-foresight path or policy result is claimed.

## All 13 initial restrictions

Gap means model minus target, in the displayed units. Ownership and fertility shares are fractions. “Ref.” is an available inverse variance or inverse squared scale for reference only. C1/C2 and N1/N2 refer to the unadopted alternatives below. The normalization is separate from the twelve proposed scored restrictions.

| Restriction / provenance | Target | Model | Gap | Uncertainty | Ref. precision | Actual weight | Loss contribution |
|---|---:|---:|---:|---|---:|---:|---:|
| [Initial model completed fertility](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/OVERNIGHT_PLAN.md) | 2.1 | 2.10003286 | 3.2862891e-05 | No empirical SE | — | — | — |
| [Childless women, ages 40–44](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/fertility_target_contract.json) | 0.198278751 | 0.180585019 | -0.0176937325 | Pooled SE unchosen | C1 | — | — |
| [Exactly one child among mothers, ages 40–44](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/fertility_target_contract.json) | 0.213655325 | 0.233236067 | 0.0195807416 | Pooled SE unchosen | C2 | — | — |
| [Period mean first-birth age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/fertility_target_contract.json) | 25.9762639 | 26.5851408 | 0.608876934 | Process / discrepancy scale unchosen | N1 | — | — |
| [First births at age 30+](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/fertility_target_contract.json) | 0.249278013 | 0.265167566 | 0.0158895535 | Process / discrepancy scale unchosen | N2 | — | — |
| [Wealth / annual gross labor earnings](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 6.14586139 | 6.09687696 | -0.0489844317 | 0.362855153 SE | 7.59509847 | — | — |
| [Annual bequests / aggregate wealth](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 0.0088 | 0.012368035 | 0.00356803499 | 0.00044 synthetic | 5165289.26 synthetic | — | — |
| [Old wealth/income p90 / median, ages 76–84](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 3.51593509 | 3.72975572 | 0.213820638 | 0.306910785 SE | 10.6163615 | — | — |
| [Mean occupied rooms, capped at 9](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 5.56109738 | 6.24638234 | 0.685284961 | 0.0883812008 SE | 128.020702 | — | — |
| [Ownership, heads 30–55](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 0.648334034 | 0.574768359 | -0.0735656752 | 0.0206752729 SE | 2339.36237 | — | — |
| [First-birth room response, −1 to +3](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 0.720246262 | 0.438846663 | -0.281399599 | 0.0852600513 SE | 137.565275 | — | — |
| [Rooms: 3+ versus 1–2 resident children (model dependent proxy)](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 0.347066932 | 0.166215341 | -0.180851591 | 0.0597051546 SE | 280.528084 | — | — |
| [Recent-parent ownership gap](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 0.162895509 | — | — | 0.00607952468 SE | 27055.823 | — | — |

**The only unavailable restriction is recent-parent ownership.** ACS requires the oldest resident own child to be under four and compares against no resident own children. The count-only model cannot reconstruct those groups. The old any-dependent-versus-never-parent statistic is not substituted. The family-room value is explicitly a dependent-count proxy for the resident-child contrast; its numerical gap does not certify exact empirical comparability.

## CPS age-projection sensitivity

The main diagnostic assumes a parity transition is uniformly timed inside the four-year cell. The sensitivity holds post-birth parity fixed throughout each cell. Both use overlap of ages 40–44 with model cells starting at 38 and 42; both remain approximations to women in the CPS.

| Observation | Target | Uniform birth time | Gap | Constant post-cell | Gap |
|---|---:|---:|---:|---:|---:|
| Childless women, ages 40–44 | 0.198278751 | 0.180585019 | -0.0176937325 | 0.167372728 | -0.0309060233 |
| Exactly one child among mothers, ages 40–44 | 0.213655325 | 0.233236067 | 0.0195807416 | 0.209841702 | -0.0038136237 |

First-birth timing is identical across these stock projections: it is calculated from first-birth flows. The exactly-one gap changes sign across projections; approximation sensitivity must remain visible before scoring these moments.

## Two validation observations

| Observation / provenance | Target | Model | Gap | Empirical SE | Actual weight / contribution |
|---|---:|---:|---:|---:|---|
| [Ownership, heads 25–34](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 0.431158377 | 0.31360335 | -0.117555027 | 0.0210233346 | — / — |
| [Old wealth/income median, ages 76–84](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json) | 7.285793 | 6.40231025 | -0.88348275 | 0.522531159 | — / — |

## Interpretation of this starting point

Wealth relative to gross labor earnings is close: 6.096877 versus 6.145861. The old wealth-dispersion proxy is also comparatively close, although its denominator is modeled pension income rather than PSID total family income. These are descriptive matches, not identification results.

Housing fit is weaker. Mean capped rooms are 6.246382 versus 5.561097, while ownership at 30–55 is 57.4768% versus 64.8334%, a deficit of 7.3566 percentage points. Young ownership is 31.3603% versus 43.1158%, a deficit of 11.7555 points. The first-birth rooms response is 0.438847 versus 0.720246. The larger-family dependent-count proxy gives 0.166215 rooms versus the resident-child target 0.347067. These gaps show that average housing quantity alone does not describe the tenure and family-housing fit.

Annual bequests/wealth are 0.0123680 versus the external normalization 0.0088. First births occur about 0.609 years later in the flow-based mean and the share at 30+ is 1.5890 percentage points above the early target. Childlessness is below its target under both CPS projections; the exactly-one share requires the projection qualification above. No policy inference or claim of an optimum follows from this one mapped point.

## Complete structural coordinates and restrictions

The following nine values and all bounds/near-bound flags are copied from the supplied parameter receipt. They are inherited or mapped diagnostic coordinates, not new estimates. Only the initial preference level is normalized.

| Structural coordinate | Value | Lower | Upper | Transform | Near bound |
|---|---:|---:|---:|---|---|
| beta_annual | 0.995276579 | 0.94 | 0.9995 | discount | False |
| kappa_fert | 2.16817304 | 0.02 | 50 | log | False |
| kappa_fert_continuation | 1.77077059 | 0.02 | 50 | log | False |
| chi | 1.05372702 | 0.1 | 5 | log | False |
| H0 | 8.43484408 | 0.2 | 80 | log | False |
| theta0 | 0.57034989 | 0 | 8 | softzero | False |
| theta1 | 0.103723951 | 0.02 | 16 | log | True |
| first_birth_fixed_cost | 4.61973138 | 0 | 8 | softzero | False |
| h_P | 0.716850009 | 0.1 | 2.3 | log | False |

The bequest wealth-shift parameter `theta1` retains the receipt’s near-bound flag. `h_P` is the mapped parenthood housing requirement; the per-dependent-child slope is fixed to zero. These restrictions leave nine searched structural coordinates in a future calibration, plus the separate normalization.

| Normalization or other restriction | Value | Status from receipt |
|---|---:|---|
| hbar_child_rooms | 0 | zero restriction |
| psi_child | 0.282703567 | normalized to 2.1 |
| payroll_tax | 0.179 | externally fixed |
| pension_period | 2.04636139 | budget derived |
| housing_supply_elasticity | 0.63 | externally fixed |
| tenure_choice_kappa | 0.005 | externally fixed |
| alpha_cons | 0.733 | externally fixed |
| sigma | 2 | externally fixed |

All 17 rows, including original status and interpretation fields, are copied verbatim in [parameters.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_fit_readout/parameters.csv). No alternative bound criterion is introduced. The fixed equivalence-scale shape and other inherited lifecycle restrictions remain as documented in the approved plan; this table is the supplied candidate receipt, not a claim that every inherited primitive was newly estimated.

## Unadopted fertility precision alternatives

The next values are references, not selected weights. CPS uses official annual generalized-variance approximations and alternative pooled covariance assumptions. NCHS process precision differs from annual temporal dispersion and deterministic window sensitivity; none is automatically an SMM discrepancy scale.

| Row | Scale option | Scale | Inverse squared scale | Interpretation |
|---|---|---:|---:|---|
| C1 | pooled_GVF_zero_crossyear_covariance | 0.00375125495 | 71063.5398 | Pooled approximate sampling-scale candidate |
| C1 | pooled_GVF_correlation_one | 0.00530503574 | 35532.3042 | Pooled approximate sampling-scale candidate |
| C2 | pooled_GVF_zero_crossyear_covariance | 0.0043074188 | 53897.1501 | Pooled approximate sampling-scale candidate |
| C2 | pooled_GVF_correlation_one | 0.00609113027 | 26952.8208 | Pooled approximate sampling-scale candidate |
| N1 | conditional_multinomial_process | 0.00224989151 | 197549.915 | Process SE, not survey sampling SE |
| N1 | annual_temporal_dispersion | 0.0845673695 | 139.828068 | Annual temporal SD, not SE |
| N1 | leave_one_year_sensitivity | 0.0396911738 | 634.763746 | Deterministic window-sensitivity scale, not SE |
| N2 | conditional_multinomial_process | 0.000168243828 | 35328216.9 | Process SE, not survey sampling SE |
| N2 | annual_temporal_dispersion | 0.00849226186 | 13866.0654 | Annual temporal SD, not SE |
| N2 | leave_one_year_sensitivity | 0.0037096109 | 72668.0123 | Deterministic window-sensitivity scale, not SE |

CPS pooled candidates condition on observed base shares; cross-year covariance and random mother-denominator weights are not certified. ACS errors are metro-resampling errors, not official ACS design SEs. PSID marginal bootstrap errors do not supply a combined cross-block covariance. The bequest scale is synthetic. All actual weights and contributions remain null in every machine-readable file.

## Measurement and provenance limits

- Housing and wealth use the new explicit uniform within-age-cell overlap, not the old nearest-label age masks. Rooms are capped at nine before income aggregation; PSID birth-response rooms remain uncapped.
- National/model households are compared with 42 MET2013 housing cities; the new direct-city footprint differs from old admitted-PUMA geography. Ownership does not reproduce the empirical DUE structure filter.
- Resident child ages and adult-child residence are absent from the model. Family-room values use dependent counts by explicit request; recent-parent ownership remains unavailable.
- The old wealth sample overlaps [76,85), but its modeled pension-plus-lump-sum income is a proxy for PSID family income. The actual lump-sum transfer here is zero. The real $1,000 income cutoff, observed-child-history selection and within-wave wealth/income dates remain unverified.
- The preserved pooled PSID first-birth contrast has a non-flat-prepath and normalization caveat. The stationary model uses identical origin policies for the pair and allows destination continuation births; its risk-set weighting differs from the event study. Initial application assumes stability.
- CPS stock projections assume the model household reproductive member represents the maternal population. NCHS timing conditions on known first births; exact residence coverage is unresolved. Boundary collapse places 7.73136% of early first births below age 18 into the first model cell; it does not create pre-18 maternal states. Unknown birth order remains a separate source caution.
- The 2.1 completed-fertility normalization and top-bin representative are distinct from literal female period TFR and the household entry law. This readout neither changes nor newly certifies those identities.
- Existing contract documents include historical “not implemented” language. The saved probe now demonstrates the diagnostic implementations; their substantive approximation and weight-activation gates remain. No moment is dropped and no older 2023 target is used here.

## Files and checks

- [13-restriction fit CSV](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_fit_readout/target_fit.csv); [validation CSV](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_fit_readout/validation_fit.csv); [both CPS projections](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_fit_readout/cps_projection_sensitivity.csv).
- [all reference precision options](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_fit_readout/reference_precision_options.csv); [complete machine-readable summary](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_fit_readout/summary.json).
- [saved measurement receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_fit_readout/source_initial_measurement_summary.json); [saved initial-loop receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_fit_readout/source_initial_smoke_summary.json); [housing/wealth provenance](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json); [fertility provenance](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/fertility_target_contract.json).

Checked the common checkpoint hash, all 13+2 row assignments, exactly one missing model row, every gap, available housing numerator/denominator arithmetic, both CPS stock-ratio calculations, source snapshot hashes, and byte-identical parameter copying. No empirical builder, model or cluster job ran for this table-building task. The lead reports the final 23 compiled checks passed; their detailed log is not among these table inputs.
