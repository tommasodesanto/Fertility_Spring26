# Native fixed-price financing diagnostic

This is a full-lifecycle partial-equilibrium diagnostic from the frozen September-14 native checkpoint. It holds checkpoint/schema inspection and one receipt plus compressed native arrays for every completed arm. Every arm applies the native sequential-fertility and current-choice period evaluator to the saved beginning-of-period mass. It does not calculate a stationary endpoint, a transition, a market-clearing price, or a fiscal-clearing transfer.

Run two independent `baseline` cases first. Both require all saved policy arrays to reproduce at `atol=1e-10, rtol=0`; only then run `mortgage_only`, `unsecured_only`, and `both` as independently timed processes.

## Submission receipt (September 19, 2026)

Torch job **18034069** (`native_finance`) was submitted from
`code/cluster/submit_e5f_native_financing_diagnostic.sh` to remote root
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a`.
The run uses the fixed native checkpoint and exact native parameters, ψ,
pension, rebate, baseline price, and growth rate; it does not perform GE
recalibration or a transition. The five full solves are two baseline controls
followed by `mortgage_only` (φ=1, λ=0), `unsecured_only` (φ=0.8, λ=5), and
`both` (φ=1, λ=5). The first error halts the batch; there are no automatic
retries. Job **18034069** completed in **3m24s** with exit 0. Both baseline
controls reproduced saved policies, current mass, and births at
`atol=1e-10, rtol=0`; all five cases passed zero budget-excess mass and
occupied wealth-value-drop gates.

Current-calibration four-year birth flows at common initial population and
prices were baseline **0.1155903878306699**, mortgage-only
**0.11638809385408035** (**+0.690114496872%**), unsecured-only
**0.12377387987346156** (**+7.079734047419%**), and both
**0.12451625195091737** (**+7.721977828574%**). These are birth flows, not
completed fertility, and are partial-equilibrium diagnostics; they do not
measure GE offsets or revised-calibration effects. Mortgage-only relaxes both
the deposit and collateral limit. The supplementary measured report is
[report/report.md](report/report.md).

The credit cap λ=5 means five times age-specific after-tax four-year income
as constructed by `build_debt_caps`, not five annual incomes; the cap is zero
after retirement. The 82/86 taper is retained only as an old experimental
specification. The PE report uses the same baseline pre-choice mass and native
sequential mapper. The diagnostic sequence is intended to establish the
current-calibration mechanism before any author-approved structural change or
recalibration; no full modern BGM earnings specification is included.

Retained native full-fit references are
[complete target fit](../paper_baseline_sep14/replay_20260917/native_output/selected_target_fit.csv)
and [complete parameter table](../paper_baseline_sep14/replay_20260917/native_output/selected_parameters.csv).
The earnings variant and recalibration/GE closure remain outstanding; no
automatic launch follows.

## Income pilot receipts

Stationary job **18040896** completed in 8m28s with exit 0; cohort job
**18041070** completed in 2m16s with exit 0. The stationary pilot held all nine structural parameters fixed and solved
the new earnings equilibrium with normalized
\(\psi=0.2429719621740803\), yielding 2.1000283897. Its six GE solves and 17
plots passed. The objective is **651.9418197098323**, versus retained
**179.2984242480252** under the identical objective. The cohort exercise is a
different fixed-preference exercise: explicit birth flows are
**1.87241076 → 1.66643427**, first-birth probability **0.82046004 →
0.76890885**, and mean grid first-birth age **24.095 → 25.702**. Mortgage
effects are small; total births and first births are reported separately.
Do not read cohort birth flows as completed fertility or call them TFR.

The stationary target-fit table below reports all 13 rows, including the
unweighted normalization row. Displayed numbers are rounded; the linked CSVs
retain full precision. `—` denotes a null weight or loss contribution.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.10003 | 2.83897e-05 | — | — |
| Childless women, ages 40–44 | 0.198279 | 0.233382 | 0.0351036 | 35532.3 | 43.7851 |
| Exactly one child among mothers, ages 40–44 | 0.213655 | 0.198902 | -0.0147534 | 26952.8 | 5.86666 |
| Period mean first-birth age | 25.9763 | 27.4911 | 1.51479 | 139.828 | 320.847 |
| First births at age 30+ | 0.249278 | 0.307885 | 0.0586071 | 13866.1 | 47.6271 |
| Wealth / annual gross labor earnings | 6.14586 | 6.69313 | 0.547266 | 7.5951 | 2.27473 |
| Annual bequests / aggregate wealth | 0.0088 | 0.00700983 | -0.00179017 | 5.16529e+06 | 16.5533 |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.55287 | 1.03693 | 10.6164 | 11.415 |
| Mean occupied rooms, capped at 9 | 5.5611 | 6.17848 | 0.61738 | 128.021 | 48.7962 |
| Ownership, heads 30–55 | 0.648334 | 0.490673 | -0.157661 | 2339.36 | 58.1499 |
| First-birth room response, −1 to +3 | 0.720246 | 1.29942 | 0.579172 | 137.565 | 46.145 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.350358 | 0.00329151 | 280.528 | 0.00303925 |
| Recent-parent ownership gap | 0.162896 | 0.119701 | -0.0431943 | 27055.8 | 50.4793 |

The full 17-parameter diagnostic table reports estimate, bound, near-bound
flag, and status. Beta is **0.99** at the active restricted upper bound;
the raw generic table's upper bound is **0.9995** and should not be confused
with the active pilot bound.

| Parameter | Estimate | Active bounds | Near active bound | Status |
|---|---:|---:|:---:|---|
| beta_annual | 0.99 | 0.94–0.99 | True | Retained value; fixed in this pilot |
| kappa_fert | 0.337734 | 0.02–50.0 | True | Retained value; fixed in this pilot |
| kappa_fert_continuation | 0.397856 | 0.02–50.0 | True | Retained value; fixed in this pilot |
| chi | 1.04965 | 0.1–5.0 | False | Retained value; fixed in this pilot |
| H0 | 8.1121 | 0.2–80.0 | False | Retained value; fixed in this pilot |
| theta0 | 0.081051 | 0.0–8.0 | False | Retained value; fixed in this pilot |
| theta1 | 0.0852036 | 0.02–16.0 | True | Retained value; fixed in this pilot |
| first_birth_fixed_cost | 0.265765 | 0.0–8.0 | False | Retained value; fixed in this pilot |
| h_P | 2.3 | 0.1–2.3 | True | Retained value; fixed in this pilot |
| hbar_child_rooms | 0 | fixed | False | zero restriction |
| psi_child | 0.242972 | normalized | False | normalized to 2.1 |
| payroll_tax | 0.179 | fixed | False | externally fixed |
| pension_period | 2.04636 | derived | False | budget derived |
| housing_supply_elasticity | 0.63 | fixed | False | externally fixed |
| tenure_choice_kappa | 0.005 | fixed | False | externally fixed |
| alpha_cons | 0.733 | fixed | False | externally fixed |
| sigma | 2 | fixed | False | externally fixed |

Source pins: [target-fit CSV](income_stationary_18040896/target_fit_comparison.csv),
[parameter CSV](income_stationary_18040896/parameters_comparison.csv),
[stationary receipt](income_stationary_18040896/collection_receipt.json), and
[cohort receipt](income_cohort_18041070/job_receipt.json). Bounded refit
preparation is in progress elsewhere; no full recalibration has been launched.

## Rental-access interaction and earnings candidate

Grouped native job **18037585** completed in 28 seconds and rental-access job
**18037587** completed in 3m13s, both with exit 0. At common prices,
pre-population, and fertility preferences, expanding the native rental room
cap from 6 to 10 rooms reduced the mortgage-only birth-flow increment by
**84.1755512677%** and the first-birth-flow increment by **88.4941545740%**.
This supports a rental space-access role; it does not identify mediation or
establish GE, recalibration, or a sign proof. All three rental cases passed
zero budget/mass/value violations, the baseline control reproduced exactly,
and each case wrote 17 standard PNGs. Group reconciliation was below `1e-12`.
The grouped mortgage comparison attributes 98.89% of the mortgage first-birth increment
to initial renters, who are 96.53% of eligible mass; renter Q3 first-birth
rate rose 0.412 percentage points and mean rooms rose 0.342. See
`grouped/metadata.json`, `grouped/groups.csv`, and
`rental_access/comparisons.csv`.

The persistent-plus-iid-transitory candidate uses the existing nested
13-moment, 3-parameter fit: objective `39.5085474236852` versus
`14.3178085110745` for the full fixed-effect fit. Flag: materially worse
earnings fit; no full-model recalibration has run. See
`earnings_candidate/README.md` and `earnings_candidate/candidate.json`.

## Reporting verification

Report-only job **18036306** completed in 11 seconds with exit 0. All four
comparison rows are available. First births use the exact loss of childless
mass at the fertility stage; housing uses realized current-tenure mass.
The four-panel figure is supplementary; the full standard policy-function
gallery and a separate fertility-probability normalization audit remain
outstanding. These results do not isolate housing mediation or completed
fertility.

The native solver stores tenure probabilities in float32. Their largest
occupied-state sum discrepancy is 6.69e-8, within the explicit reporting
tolerance 2e-11 plus one storage epsilon (1.1923e-7 total). Location
probabilities sum exactly to one on inspected states. This reporting check
never changes probabilities or the model's existing acceptance checks.
Earlier reporting jobs failed on precision assumptions and serialization;
none reran the household model. See `review_receipt.json`.

To regenerate the supplementary packet on Torch (Anaconda 2025.06), run the
following on a compute node, using the retained experiment and native source:

```sh
experiment=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a
native=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches
python "$experiment/code/model/tools/build_e5f_native_financing_report.py" --input "$experiment/run" --output "$experiment/report" --checkpoint "$native/baseline_replay_20260917/replay/case/evaluation/raw/repetition_02/initial_state.pkl.gz" --source-root "$native/final_night_20260913/corrected_initial_source_v2/code/model"
```

## Earnings cohort launch and stop

Job **18037924**, submitted with a one-hour allocation, failed after 20 seconds
at `pre-mass changed at age 0` in the baseline cohort. No candidate earnings
comparison or recalibration was produced. The feasibility projection changed
the standardized entrant distribution; correction and a cluster smoke test
are required before resubmission. See `income_cohort_18037924/status.json`.
Launcher: `code/cluster/submit_e5f_native_income_followup.sh`.

Job **18040822** replaces the independent redraw with the native conditional
entry rule for both earnings processes, following gate-only job **18040625**.
The original native entrant has zero gate-induced L1 change; the independent
redraw changes by 0.0023731346586869503. All feasibility checks remain active.
The entry asset distribution may change with income. See
`income_cohort_18040822/submission.json`; this is a submitted diagnostic,
not a completed recalibration.

Job **18040822** subsequently stopped after both old-income cohorts completed:
the first candidate solve reached the standard plot routine with an age-zero-only
distribution, causing division by zero for empty older ages. No completed
candidate comparison is claimed. The corrected plotting packet must use the
normalized full lifetime cohort cross-section.

Independent stationary pilot **18040896** passed candidate-specific frozen
wrapper preflight and is running. It keeps nine structural parameters fixed,
solves the new income specification from fresh policies/distributions, and
normalizes fertility to 2.1 under the unchanged target contract. It permits
at most eight stationary solves (1800 seconds native, 2100 seconds wrapper).
See `income_stationary_18040896/submission.json` and
`earnings_candidate/calibration_plan.remote.json`. It is not a structural search.

Latest cohort submission: **18041070**. The graph packet now uses a normalized
full lifetime cohort distribution after all ages complete, preserving the
original standard plot set. Stationary pilot **18040896** has entered its first
native equilibrium solve. Both jobs have one-hour allocations; no completed
new-income comparison or recalibration is claimed at submission.
