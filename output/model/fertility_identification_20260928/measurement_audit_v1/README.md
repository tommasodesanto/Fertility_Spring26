# Income and fertility measurement audit

**2007 stationary reference — block0506, September 28 verified export.**
The common-primary loss remains **19.581**. This packet inspects existing
results and prepares a numerical experiment; it contains **zero model solves**
and adopts no new specification, target, parameter or reference.

The stationary economy approximates the 2007 distribution under the author's
deliberate replacement-fertility restriction. For each calibration proposal,
the outer loop adjusts the child-benefit parameter \(\psi\) until completed
fertility is 2.1 and checks demographic renewal. The subsequent transition is
toward 2023; re-estimation of that transition remains deferred. The inherited
export filename `lifecycle_2023.csv` does not make its stationary observations
2023 transition results.

## Reference identity and verification

The source is the immutable Torch staging tree
`/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project`,
bound inside the container to the historical local project path. It is not an
assertion that today's shared working checkout still has identical source.
The inspected export is `../resume_v1/selected_export/primary/`; its checkpoint
symlink resolves to `../resume_v1/repeat_0212_primary/case/initial_state.pkl.gz`.

The [identity receipt](reference_identity.json) records the full identities:

- Contract: `68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf`.
- Source manifest: `07d84336a3112b251afe505908113d9c00585b91c34bd50f0dee108435db496d`.
- Actual inspected checkpoint: `b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d`.
- Parameter table: `4d229fe18ee9c43a8a709413a10690a845d0cf7028efe57925d4f4b1d7c8d036`, identical to the original block0506 table.
- Primary objective: `04528ece6513e4da9436bb8b35b54cd7e9cecd3a6695d050a0a7aea5a71d5e60`.
- Target/weight fingerprint: `20e531855075b807886da8e936930dca9090b989b2ce381de39042cf206efdbd`.

Torch job **18735483** passed in 35 seconds: full source ancestry, checkpoint,
all 30 export artifact hashes, all 17 standard plot hashes, independent loss
arithmetic, and all saved pre/post fertility count distributions against the
checkpoint arrays. The 17 original plots were retained without re-rendering;
their earlier visual review remains the visual evidence. No checkpoint was
downloaded, copied into Git, or included in a new repository snapshot.

## Actual saved income process

The matrix stored in `P.Pi_z`, the saved solution's `income_transition`, and
the effective runtime matrix agree to \(5.6\times10^{-17}\). An independent
Rouwenhorst recursion using the approved external estimate reproduces the
saved 15-by-15 matrix to \(1.1\times10^{-16}\), the income grid to
\(5.4\times10^{-15}\), and the binomial stationary weights exactly.

| Income object, per four-year period | Approved | Saved matrix implies |
|---|---:|---:|
| Log persistence | 0.735 | 0.735 |
| Log innovation standard deviation | 0.484 | 0.484 |
| Stationary log variance | 0.508 | 0.508 |
| First-lag log covariance | 0.374 | 0.374 |
| Second-lag log covariance | 0.274 | 0.274 |
| Mean income multiplier | 1.000 | 1.000 |

The empirical second-lag covariance is 0.290: the approved single-process
approximation misses that held-out covariance by 5.354%. This is a known
restriction of the fitted process, distinct from the verified implementation.
The saved flags select a Markov process, disable fixed permanent income types,
and set retirement income dependence on the shock to zero. There is no separate
iid component in this 15-state process. Row sums and stationarity pass at
machine precision. The alternative reset-to-stationary matrix differs by up
to 0.598, so the constructor conjecture is resolved in favor of the approved
Rouwenhorst matrix.

[Full income receipt](income_audit.json), [actual matrix](saved_income_transition.csv),
[independent reconstruction](approved_rouwenhorst_transition.csv). The external
estimate is the existing `single_process_external_estimate.json` in
`output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/`.
This verifies the inherited configuration, not its empirical adequacy or grid
convergence. Entry wealth/income coupling and stationary income risk at entry
remain inherited approximations.

## Early fertility: motherhood versus additional children

Let \(s\) be the share who are mothers and \(k\) their mean capped number of
children. Then early fertility is \(e=s k\). All child counts here are children
ever born, capped at three, not children currently at home.

| Age-25 object | CPS | Model | Model minus CPS |
|---|---:|---:|---:|
| Mothers, percent | 45.725 | 45.011 | -0.714 pp |
| Children among mothers | 1.770 | 1.190 | -0.581 |
| Children per woman | 0.810 | 0.535 | -0.274 |

The symmetric exact decomposition is
\[
e_M-e_D=(s_M-s_D)\frac{k_M+k_D}{2}
       +(k_M-k_D)\frac{s_M+s_D}{2}.
\]
The motherhood component is -0.011 children; the conditional-count component
is -0.264, or **96.144%** of the total gap. These are arithmetic components,
not causal effects or a proof about attainable fit. Replacing the count target
with motherhood would remove information about the margin that currently
fails; it is not an innocuous relabeling or an adopted proposal.

The model age-25 shares with zero, one, two, and three-or-more children are
0.550, 0.365, 0.085, and 0.000. Entry at 18 and at most one birth in each
four-year cell rule out a third birth by this age. CPS has a weighted 2.916%
with **more than** three births, so its three-or-more share is at least that
large; the exact age-25 count distribution is not in the compact receipt.
Quantifying the total effect of omitted pre-18 or closely spaced births needs
birth-history measurement or a separately specified experiment.

The empirical source pools June 2004 and 2006 women at completed interview
age 25, valid `FREVER` 0–20 and positive fertility-supplement weights, with
1,774 observations. Pooling weights observations, not annual means. The raw
uncapped mean is 0.857; the capped-three mean is 0.810. Its reported bootstrap
SE is 0.028, based on person resampling within year; it is not a survey-design
variance estimate. Motherhood and conditional-count uncertainty is not yet
estimated. The source is `output/model/early_fertility_target_20260926/early_fertility_target.json`
and `code/data/cps_fertility/build_early_fertility_target.py`.

Completed age 25 means the interval \([25,26)\). The model measures its mean
using 0.125 of the pre-birth and 0.875 of the post-birth distribution in
\([22,26)\), assuming uniform birth timing. This was checked against the
serialized arrays, not inferred from the constructor. The two CPS samples
represent different birth cohorts, approximately mid-1978–mid-1979 and
mid-1980–mid-1981. They are not a period fertility rate or one cohort's panel.

[Full decomposition](early_fertility_decomposition.json) and
[checkpoint replay](fertility_checkpoint_replay.json).

## First-birth timing and the lifecycle

The author-requested [three-panel lifecycle figure](fertility_lifecycle_comparison.png)
shows children per woman, motherhood, and children among mothers over identical
five-year age windows. `plot_lifecycle.py` renders only the saved table on Torch;
`lifecycle_plot_qa.json` verifies all six plotted series against that table.
It is supplemental and leaves the standard 17 graphs intact. The author is
considering a wider-age fertility target; no replacement or averaging rule is
adopted. A candidate must retain information on additional births at young ages
and be checked for sensitivity to the first-/later-birth taste scales,
first-birth cost and child-benefit curvature before changing the ten-moment
system. Close late-age levels do not by themselves validate the early timing.

| First-birth age cell | NCHS share, percent | Model share, percent |
|---|---:|---:|
| 18–21, with empirical younger tail collapsed | 34.170 | 32.648 |
| 22–25 | 21.556 | 25.993 |
| 26–29 | 19.347 | 18.994 |
| 30–33 | 14.256 | 11.614 |
| 34–37 | 7.505 | 6.050 |
| 38–41 | 2.631 | 2.908 |
| 42–45, with empirical older tail collapsed | 0.536 | 1.793 |

The data are 6,611,269 pooled period first births in 2003–2006, not the
birth histories of the CPS age-25 women. The target maps ages 12–21 to midpoint
20 and ages 42–49 to midpoint 44. First births before 18 comprise 7.731% of the
data. The raw completed-integer-age mean is 25.161; the mapped mean used in
calibration is 25.976. The model mapped mean is 25.933. Thus close means coexist
with visible differences in the seven shares. Raw-builder residence filtering
and unknown birth-order coverage remain unresolved; no new claim of resident-only
coverage is made. [Exact shares](first_birth_age_cells.csv) and
[timing conventions](first_birth_timing.json).

The following newly aligned five-year windows use the existing pooled CPS
weighted count shares. Model windows integrate the same uniform pre/post
interpolation over the identical age intervals, retaining model age masses.
They are supplemental cross-sectional diagnostics, not a replacement for the
17 standard graphs or a cohort trajectory. This alignment differs from the
earlier four-year midpoint display in `morning_fits/`.

| Age window | CPS capped children | Model capped children | Gap |
|---|---:|---:|---:|
| 20–24 | 0.475 | 0.308 | -0.167 |
| 25–29 | 0.990 | 0.697 | -0.293 |
| 30–34 | 1.476 | 1.096 | -0.379 |
| 35–39 | 1.709 | 1.453 | -0.256 |
| 40–44 | 1.718 | 1.731 | 0.012 |

The discrepancy concentrates before age 40. The [complete lifecycle table](fertility_lifecycle_matched_windows.csv)
also reports motherhood, children among mothers, three-or-more shares, and
raw versus model-coded counts. The source age-profile table is
`output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/age_profile/age_profile_candidates.csv`.
There is no CPS observation for the full 42–45 model cell in this source;
do not silently treat ages 40–44 as that cell.

## Replacement stationarity and the Claude review

The imposed 2.1 exceeds CPS ages-40–44 uncapped fertility 1.878 by 0.222;
against CPS with the model's top-bin coding, 1.891, the difference is 0.209.
This quantifies part of the maintained stationary approximation. It is not
a discrepancy between the model and an empirical 2007 period fertility rate.
The exact saved three-or-more weight is 3.602359422009, authenticated in the
parallel frozen-reference manifest against the same checkpoint identity.

The model's **observed ages-40–44** shares are 0.201, 0.167, 0.331, 0.300,
versus CPS 0.198, 0.171, 0.344, 0.286. Its model-coded mean in that window is
1.912. Its terminal post-birth population instead has a three-or-more share
of 0.382 and the normalized model-coded mean 2.100. Births still occur in
the last fertile model cell. These populations cannot be substituted for one
another. Claude's approximately 0.43 three-or-more share follows from imposing
the replacement mean on the wrong age-window distribution; it is not the
saved result. [Full comparison](stationary_approximation_comparison.json).

Claude's 0.665 calculation is correct conditional on fixing the empirical
first-birth cell shares, identifying ages-40–44 motherhood with eventual
motherhood, and converting period shares into cohort lifetime probabilities.
The active target system does not impose all those conditions. Matching one
mean and adding up the seven shares leaves five free directions. For example,
changing shares at midpoints 20, 24, 28 by
\((\epsilon,-2\epsilon,\epsilon)\) preserves the mean and total share.
It changes the early-fertility ceiling. This disproves the claimed inference,
not the possibility that the complete active system is difficult to fit.

The verified arithmetic does not establish that “most of the gap is time
aggregation,” nor justify dropping or replacing the early count target. A
coherent period/cohort comparison still requires the correct female exposures,
age-specific total birth rates, and cohort birth histories. The observed
replacement/cohort mean difference alone cannot identify its effect on the
age-25 gap. None of these limitations reverses the author's 2007-to-2023
timeline or the deliberate replacement closure.

## Identification and a bounded numerical experiment

The [plain-English Jacobian readout](jacobian_readout.md) connects all ten
parameter directions to early fertility and first-birth age, distinguishes
verified profiles from linear predictions, and explains weak local directions.
The author's later two-birth suggestion is examined separately in the
[fixed-policy spacing replay](../two_births_v1/README.md), a diagnostic with
zero model solves, not an adopted specification or a new calibration.

The active system has ten scored moments and ten searched coordinates, plus
the separately normalized \(\psi\). The full/half-step Jacobians are formally
rank ten, but only seven half-step singular values exceed the step-difference
norm. Their condition numbers are about \(2.69\times10^6\) and
\(4.65\times10^5\). This is a numerical sensitivity warning, not a statistical
rank test or a proof of global identification. Working weights are not a
joint efficient SMM covariance matrix. Neither target count nor formal rank
certifies precise identification.

The saved Jacobian reproduces Claude's **unsolved** ridge-10 linear prediction
of loss 9.519. His displayed coordinates are rounded, not exact inputs. The
prepared [full-precision proposal](proposed_step_parameters.csv) uniformly
damps the log step to a maximum magnitude of 0.150. Full- and half-step
Jacobians predict primary losses 13.335 and 13.346, respectively. The early
moment is predicted to fall to 0.534: this proposal tests improvement in the
other moments, not a resolution of the early-fertility problem.
[All fourteen predicted rows](proposed_step_predictions.csv) retain roles,
weights and predicted contributions; they are not calibration results.

The [bounded experiment plan](bounded_experiment_plan.json) proposes two
evaluations of that identical parameter point. One uses the retained initial
normalization guess; the other uses the Jacobian-predicted guess for \(\psi\).
Both retain the same bracket step, price initialization, targets, weights,
boundaries, full 2.1 normalization, renewal, and every scientific gate. The
reference required six stationary solves, totaling 970.676 seconds. Warm-start
savings remain unmeasured.

Proposed changes are **diagnostic recalibration of the ten fitted parameters**
and the corresponding derived child-benefit level; none is adopted. Earnings,
entry wealth/income, timing, transfers and floors, preference functional forms,
targets and weights remain fixed. The only difference between the paired arms
is the initial normalization guess. Maximum scope: two objectives, 23 stationary
solves each, 1,800 seconds per objective, two Torch workers, 70 minutes total.
The plan defines additional cross-arm equivalence screens in advance and
requires a synthetic dispatch/failure/timeout smoke, immutable opt-in source,
complete fit/bounds tables and the unchanged 17-plot packet. Stop after the pair
or an integrity/unknown failure; no automatic extension, search revival or
promotion. Final exact-repeat verification would be separate work.

**No experiment was launched.** The normalization-start adapter and its launch
smoke have not been implemented. Current task authority is used for this
read-only audit and concrete experiment design. Presentation work and
transition re-estimation remain deferred.

## Reproduction and retained failures

`audit_saved.py` runs with `sbatch run.sh` inside the pinned Torch stage.
`compare_and_plan.py` runs in the same container with one core, 4 GiB and a
five-minute cap; it reads only compact tables and writes this packet.
`inputs/` on Torch contains exact copies of the approved income receipt, early
fertility receipt, pooled age-profile table, NCHS counts, Jacobian/SVD tables,
and the same-checkpoint frozen-reference manifest. Their original paths are
given above and their staged hashes are in `small_table_input_hashes.json`.
The scripts require a Slurm job and refuse Mac execution.

The first inspection, 18735365, stopped after checkpoint loading because the
small external-income JSON was absent from the snapshot. Inspection 18735393
completed the income comparison but its final fertility assertion treated a
named-share dictionary as an array. Table job 18735474 had the same schema
mistake. Failed scripts and logs are retained under `failed_attempt_1/` and
`failed_attempt_2/`. Explicit key extraction fixed the reporting scripts;
18735483 and 18735484 then passed. Final small-table job18735908 passed in four
seconds using the exact saved top-bin weight rather than its three-decimal
display. A subsequent wording-only clarification distinguishes timing from
count coding; no numeric result changed. An independent read-only review
checked the decomposition, age-window integration and Claude corrections.
These were audit/staging
repairs, with no model or parameter changes and no equilibrium solves.

## Complete reference fit and parameter restrictions

The following tables reproduce the authenticated primary export. Display is
limited to three decimals; CSVs preserve full precision. Normalization,
scored targets and zero-weight validation rows are separated. The full 31-row
parameter table includes all ten searched coordinates and their bounds,
normalized/derived quantities and external restrictions. The two fertility
scale flags use the retained raw-range near-bound rule; they do not mean the
optimizer is sitting on the lower endpoint.

**Separate normalization**

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Replacement completed fertility | 2.100 | 2.100 | -1.663e-06 | — | — |

**Ten scored targets**

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Childlessness, ages40–44 | 0.198 | 0.201 | 0.003 | 35532.304 | 0.275 |
| Exactly one child among mothers, ages40–44 | 0.214 | 0.209 | -0.004 | 26952.821 | 0.483 |
| Mapped mean first-birth age | 25.976 | 25.933 | -0.043 | 139.828 | 0.264 |
| Net worth / mean annual working earnings | 6.927 | 6.326 | -0.600 | 7.595 | 2.737 |
| Annual positive estates / wealth | 0.007 | 0.007 | -2.336e-04 | 5.165e+06 | 0.282 |
| Mean rooms | 5.729 | 5.848 | 0.119 | 128.021 | 1.798 |
| Ownership, ages30–55 | 0.676 | 0.655 | -0.021 | 2339.362 | 1.076 |
| First-birth rooms response | 1.465 | 1.622 | 0.157 | 137.565 | 3.395 |
| Recent-parent ownership gap | 0.128 | 0.120 | -0.008 | 27055.823 | 1.758 |
| Children ever born at25, capped3 | 0.810 | 0.535 | -0.274 | 100.000 | 7.513 |

**Three zero-weight validation rows**

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| First births at30+, share | 0.249 | 0.224 | -0.026 | 0.000 | 0.000 |
| Older wealth/income p90 / median | 3.516 | 3.069 | -0.447 | 0.000 | 0.000 |
| Rooms gap,3+ versus1–2 children | 0.385 | 0.353 | -0.032 | 0.000 | 0.000 |

**All parameters and restrictions**

| Parameter | Estimate | Lower | Upper | Near bound | Role / restriction |
|---|---:|---:|---:|---|---|
| H0 | 6.294 | 0.200 | 80.000 | False | free in evening DUE calibration |
| beta_annual | 0.963 | 0.940 | 0.990 | False | free in evening DUE calibration |
| chi | 1.094 | 0.100 | 5.000 | False | free in evening DUE calibration |
| first_birth_fixed_cost | 0.621 | 0.000 | 8.000 | False | free in evening DUE calibration |
| kappa_fert | 0.176 | 0.020 | 50.000 | True | free in evening DUE calibration |
| kappa_fert_continuation | 0.332 | 0.020 | 50.000 | True | free in evening DUE calibration |
| theta0 | 0.125 | 0.000 | 8.000 | False | free in evening DUE calibration |
| delta_alpha_jump | 0.135 | 0.000 | 0.250 | False | free in evening DUE calibration |
| child_benefit_curvature | 0.061 | 0.000 | 0.800 | False | free in evening DUE calibration |
| tenure_choice_kappa | 0.012 | 0.001 | 0.100 | False | free in evening DUE calibration |
| psi_child | 0.136 | — | — | — | normalized to completed fertility 2.1 |
| child_benefit_CRRA_coefficient | 0.127 | — | — | — | derived from normalized one-child benefit |
| theta1 | 0.008 | — | — | — | fixed external restriction |
| sigma | 2.000 | — | — | — | fixed |
| alpha_cons | 0.733 | — | — | — | fixed CEX childless expenditure share |
| delta_alpha | 0.000 | — | — | — | fixed zero later-child loading |
| h_P | 0.000 | — | — | — | no housing floor |
| utility_reference_rent | 0.110 | — | — | — | fixed substantive utility normalization |
| q_annual | 0.020 | — | — | — | author-retained 2% annual real rate |
| financed_share | 0.800 | — | — | — | inherited credit contract |
| housing_supply_elasticity | 0.630 | — | — | — | fixed provisional external mapping |
| payroll_tax | 0.080 | — | — | — | derived from adopted pension ratio |
| pension_period | 0.918 | — | — | — | balanced PAYGO |
| annual_depreciation | 0.014 | — | — | — | adopted |
| period_depreciation | 0.055 | — | — | — | compounded |
| annual_property_tax | 0.011 | — | — | — | adopted |
| period_property_tax | 0.042 | — | — | — | linear period convention |
| selling_cost | 0.060 | — | — | — | retained |
| rental_cap | 6.000 | — | — | — | retained provisional |
| wealth_grid_nodes | 160.000 | — | — | — | retained exact grid |
| income_states | 15.000 | — | — | — | retained B15 |

Full precision: [target fits](../resume_v1/selected_export/primary/target_fit.csv), [parameter bounds/restrictions](../resume_v1/selected_export/primary/parameters.csv).
