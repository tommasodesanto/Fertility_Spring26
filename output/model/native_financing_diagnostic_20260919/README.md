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
