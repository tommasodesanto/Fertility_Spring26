# Implementation Status: Intergenerational Housing Fertility

Updated: 2026-06-07

## Rule

Do not leave simplifications implicit. If code uses a shortcut, record it here
with one of these labels:

- `INTENDED`: part of the current target model.
- `SIMPLIFICATION`: deliberately simpler than a richer target.
- `NOT IMPLEMENTED`: target object not coded yet.
- `DIAGNOSTIC ONLY`: produced for inspection, not a calibration target.

## Current Pass

The current code is a runnable first-pass scaffold, not a calibrated
quantitative model.

- `INTENDED`: one aggregate housing-services market. The old multi-market
  choice dimension from the workhorse code is kept with `I=1`, so it is
  mechanically present but economically degenerate.
- `INTENDED`: Coven-style user-cost closure in stationarity:
  \[
  q=(r+\delta+\tau^p)P.
  \]
- `INTENDED`: aggregate housing services clear against an upward-sloping supply
  curve:
  \[
  H^D(q)=H^S(P)=H_0\left(\frac{q}{\bar q}\right)^\eta.
  \]
- `INTENDED`: renter housing is continuous up to `hR_max`; owner housing is a
  discrete rung choice on `H_own`.
- `INTENDED`: tenure/rung choices include a small Type-I-EV taste-shock
  smoothing term, `tenure_choice_kappa=0.01`. This is the model's current
  discrete-choice smoothing device and is also needed for reliable market
  clearing with discrete owner rungs.
- `INTENDED`: owner adjustment includes transaction/sale wedges and the
  workhorse down-payment constraint. In this code, `phi` is the financed share,
  so the down-payment threshold is \((1-\phi)Ph\).
- `INTENDED`: new owner choices face a payment-to-income screen. The current
  implementation blocks new purchases/adjustments when
  \[
  (q_{\text{int}}\phi+\tau^p)Ph > \psi^{PTI} y_a.
  \]
  Incumbent owners who keep the same owner rung do not requalify each period.
- `INTENDED`: lifecycle income profile \(y_a\) enters the budget and the
  payment-to-income screen.
- `INTENDED`: persistent Markov income/productivity state \(z_t\). The default
  grid is `z_grid=[0.70, 1.00, 1.30]` with stationary entry weights
  `z_weights=[0.30, 0.40, 0.30]`. The transition matrix is
  `Pi_z = rho_z I + (1-rho_z) z_weights`, with default
  `rho_z=0.85`. The Bellman equation takes expectations over
  \(z_{t+1}\sim\Pi_z(z_{t+1}\mid z_t)\), and the discrete-time forward
  distribution propagates mass across \(z_{t+1}\) with the same matrix.
  Working-age income is \(y_{a,z}=z y_a\), and the same type-scaled income
  enters the budget and payment-to-income screen.
- `INTENDED`: 4-year decision periods. Defaults are `period_years=4`,
  `age_start=22`, `J=16`, and `J_R=11`, so agents enter at age 22, first retire
  at age 66, and die after the last age-82 decision period.
- `INTENDED`: one dependent-child state with stochastic age-out duration
  `stage_durations=[A_m/period_years]=[4.5]`. With `A_m=18`, this gives an
  expected child-at-home duration of 18 years, so child maturation into the
  adult entrant pool is mapped to the 4-year decision period rather than
  happening after one period.
- `INTENDED`: period-unit flows. Discounting, interest, depreciation, property
  tax, income, child goods costs, and minimum consumption are scaled to the
  4-year decision period in the default parameter file.
- `INTENDED`: fertility accounting follows the old workhorse convention.
  `mean_completed_fertility` is the household/unitary-agent mean parity; the
  reported and targeted TFR moment is `tfr = 2 * mean_completed_fertility`.
  This factor-of-two conversion is not a demographic transition change; it is
  the project convention for mapping unitary household parity into the
  paper-facing fertility scale.

## Simplifications

- `SIMPLIFICATION`: the Markov income grid and transition probabilities are
  smoke-test defaults, not a calibrated discretization of an estimated earnings
  process.
- `SIMPLIFICATION`: `tenure_choice_kappa=0.01` is not calibrated. With
  deterministic tenure/rung choices, the one-market price solver can bracket
  housing excess but still stop at a small residual because aggregate demand
  jumps at tenure thresholds.
- `SIMPLIFICATION`: retirement income is the common pension flow. It is not
  currently indexed to the household's realized income history or current
  \(z_t\).
- `SIMPLIFICATION`: the Markov-income path currently uses full Bellman solves
  in each price iteration. Howard policy-evaluation acceleration is still
  available only in the non-Markov inherited path.
- `INTENDED`: the active Markov-income forward path computes the old
  event-study housing response statistics by propagating birth cohorts under
  the same tenure, savings, child-aging, and income-transition rules used by
  the main forward distribution. The price iteration skips these event
  statistics and recomputes them once at the accepted price.
- `SIMPLIFICATION`: fertility remains the workhorse one-shot completed-family
  choice for childless fertile households. It is not a sequential parity hazard.
- `SIMPLIFICATION`: the single dependent-child state is not age-resolved.
  Maturation is geometric with mean 18 years, not deterministic child ages
  \(0,\ldots,17\). This keeps the state space small while preserving the
  adult-maturation timing implied by 4-year periods.
- `SIMPLIFICATION`: the model uses a collateral-constrained user-cost shortcut:
  \(qh\) enters the flow budget and \((1-\phi)Ph\le b\) enters the new-owner
  feasibility screen. The down payment is not subtracted as a separate asset
  purchase and housing equity is not a separate continuous state.
- `SIMPLIFICATION`: no mortgage coupon, amortization, maturity, refinancing, or
  present-value mortgage lock-in state.
- `SIMPLIFICATION`: the payment-to-income screen uses a simple interest-plus-tax
  payment. It is not yet Coven's full amortized mortgage payment formula.
- `SIMPLIFICATION`: the housing supply shifter `H0` and rent anchor `r_bar` are
  smoke-test normalizations, not calibrated targets.
- `SIMPLIFICATION`: old retention currently comes only through the inherited
  workhorse owner state, transaction/sale wedge, and bequest utility. A clean
  policy-created old-retention wedge is not yet implemented.

## Not Implemented

- `NOT IMPLEMENTED`: formal calibration, counterfactual tables, and production
  parameter search. The current random-search tool is diagnostic only.
- `NOT IMPLEMENTED`: estate-tax counterfactuals, inheritance kernels, bequest
  principal adding-up, and estate-revenue rebates.
- `NOT IMPLEMENTED`: mortgage-rate lock-in from coupon gaps.
- `NOT IMPLEMENTED`: transition dynamics.
- `NOT IMPLEMENTED`: landlord balance sheets by house size. The current model
  clears aggregate housing services; rents are pinned by user cost.

## Verification Targets

Every nontrivial change should pass:

```bash
python -m compileall -q code/model/intergen_housing_fertility
cd code/model && .venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
```

For model inspection, also run:

```bash
cd code/model && .venv/bin/python -m intergen_housing_fertility.cli diagnostics --fixed-prices --quiet
```

## Diagnostic Packet

Current diagnostic packets include aggregate lifecycle moments, market-clearing
plots, relative housing-market residual plots, outcome bars by income state,
age profiles by income state, childless-renter wealth distributions at selected
ages, and childless-renter policy functions at two fertile-window ages. The
policy plots show consumption, post-tenure housing services, expected
fertility, and owner-entry policy over liquid wealth with separate lines for
each \(z_t\).
Policy plots mask value-function states below \(-10^9\), since logit
probabilities at all-infeasible grid points are not economically meaningful.

## Diagnostic Calibration

- `DIAGNOSTIC ONLY`: `python -m intergen_housing_fertility.cli calibrate-small`
  runs a checkpointed random search. Its default target set is
  `old_nonlocation`, which uses the old workhorse targets that remain defined
  in the one-market code: TFR, childlessness, mean age at first birth,
  ownership, family ownership gap, child-linked housing increments, young
  liquid wealth relative to income, old-age ownership, and old
  parent-childless ownership gap. It explicitly excludes the old targets that
  require multiple markets: gradients, center/periphery shares, migration, and
  inversion targets.
- `DIAGNOSTIC ONLY`: the legacy `core` target set is retained only for quick
  smoke calibration of ownership, young ownership, old ownership, completed
  fertility, and childlessness. It is not the requested old-target objective.
- `DIAGNOSTIC ONLY`: 2026-06-07 fixed the old-target `tfr` extractor to use
  the workhorse convention `tfr = 2 * mean_completed_fertility`. Earlier
  2026-06-05 and 2026-06-06 one-market diagnostic losses evaluated `tfr` on the
  raw household-parity scale and are therefore not comparable as objective
  values. Their market-clearing and event-statistic smoke information remains
  useful.
- `DIAGNOSTIC ONLY`: 2026-06-07 local one-case objective smoke
  `output/model/intergen_housing_fertility_age_mapping_smoke/calibrate_one_case/`
  verifies the corrected age and fertility accounting: `period_years=4`,
  `stage_durations=[4.5]`, expected child-at-home duration `18` years,
  household mean parity `0.699`, reported `tfr=1.399`, and market residual
  \(8.09\times 10^{-5}\).
- `DIAGNOSTIC ONLY`: 2026-06-07 Torch smoke job `10574606`, run tag
  `intergen_old_nonlocation_age_tfr_smoke_20260607`, completed a 2-case array
  with the corrected age and `tfr` accounting. It wrote `metadata.json`,
  `cases.jsonl`, `best.json`, and `summary.json` under
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_intergen_housing_fertility_intergen_old_nonlocation_age_tfr_smoke_20260607/`.
- `DIAGNOSTIC ONLY`: 2026-06-07 Torch overnight screening job `10574658`, run
  tag `intergen_old_nonlocation_age_tfr_overnight_20260607`, was launched with
  `64` tasks, `48` cases per task, `J=16`, `Nb=70`, `n_house=6`,
  `max_iter_eq=60`, and target set `old_nonlocation`. This is an overnight
  diagnostic random search, not a formal production calibration.
- `DIAGNOSTIC ONLY`: 2026-06-05 Torch run
  `intergen_old_nonlocation_20260605` used the new one-market code and the
  old non-location target subset. It completed `48` tasks and `2,304` valid
  cases. The scalar best had loss `119.007`, but it is not economically useful:
  TFR `0.117`, childlessness `0.891`, ownership `0.937`, and old ownership
  `1.000`. Across all cases, both child-linked housing increment diagnostics
  were exactly zero because the active Markov-income statistics path did not
  compute those event-study moments yet.
- `DIAGNOSTIC ONLY`: 2026-06-06 local smoke
  `output/model/intergen_housing_fertility_old_nonlocation_smoke_after_eventstats_fast/`
  verifies that the Markov-income event-study statistics are no longer
  mechanical zeros. Baseline smoke moments include
  `housing_increment_0to1=-0.201` and `housing_increment_1to2=0.044`.
- `DIAGNOSTIC ONLY`: 2026-06-06 Torch run
  `intergen_old_nonlocation_eventstats_20260606` used the new one-market code
  and the old non-location target subset after patching the Markov-income
  event-study statistics. It completed `48` tasks and `2,304` valid cases.
  Scalar best loss was `145.253`, but the scalar-best point remains an
  economically poor corner: TFR `0.117`, childlessness `0.891`, ownership
  `0.937`, old ownership `1.000`, `housing_increment_0to1=2.610`, and
  `housing_increment_1to2=0.308`. The best more interpretable cases still have
  late births, weak family ownership gaps, and poor old-age parent-childless
  gaps. The next step is to inspect policy functions and tighten/search around
  economically admissible regions rather than treating the scalar best as a
  candidate calibration.
- `DIAGNOSTIC ONLY`: 2026-06-04 local run
  `output/model/intergen_housing_fertility_small_calibration/` used
  `J=12`, `Nb=40`, four owner rungs, and 24 checkpointed cases. Best case:
  market residual \(3.37\times 10^{-5}\), ownership `0.686`, young ownership
  `0.255`, old ownership `0.995`, completed fertility `0.847`, childless rate
  `0.298`. The associated plot packet is
  `output/model/intergen_housing_fertility_small_calibration_best_diagnostics/`.
- `DIAGNOSTIC ONLY`: a targeted 24-case probe around higher entry wealth and
  looser rental caps did not improve the diagnostic loss. Current evidence is
  that the owner/renter block creates a sharp tradeoff: lowering the rental cap
  can generate young ownership, but it tends to push old ownership too high and
  depress fertility.
