# Implementation Status: Intergenerational Housing Fertility

Updated: 2026-06-24

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
- `INTENDED`: `chi` is the reduced-form owner utility/service premium from the
  May 29 quantitative slides. It scales residual owner housing services after
  the family-space requirement is subtracted, \(\chi[H-\bar h(n,cs)]\), rather
  than making owner physical rooms count as \(\chi H\) before the floor. Thus
  ownership is more attractive in utility, but the family-space floor is still
  measured in physical rooms.
- `INTENDED`: the default model does not impose a payment-to-income screen:
  `use_pti_constraint=False`. If that optional screen is manually enabled, new
  owner purchases/adjustments must satisfy the simple payment test
  \[
  q_{\text{int}}D+\tau^p Ph \le \psi^{PTI} y_a,
  \]
  where \(D=\max\{Ph-\text{cash},0\}\) is actual transaction debt, not the
  maximum allowed LTV debt \(\phi Ph\). Equivalently, the code raises the
  required cash threshold only when the PTI-implied minimum cash exceeds the
  collateral down payment \((1-\phi)Ph\). Incumbent owners who keep the same
  owner rung do not requalify each period.
- `INTENDED`: lifecycle income profile \(y_a\) enters the budget and, only if
  manually enabled, the optional payment-to-income screen.
- `INTENDED`: persistent Markov income/productivity state \(z_t\). The default
  grid is `z_grid=[0.70, 1.00, 1.30]` with stationary entry weights
  `z_weights=[0.30, 0.40, 0.30]`. The transition matrix is
  `Pi_z = rho_z I + (1-rho_z) z_weights`, with default
  `rho_z=0.85`. The Bellman equation takes expectations over
  \(z_{t+1}\sim\Pi_z(z_{t+1}\mid z_t)\), and the discrete-time forward
  distribution propagates mass across \(z_{t+1}\) with the same matrix.
  Working-age income is \(y_{a,z}=z y_a\), and the same type-scaled income
  enters the budget and the optional payment-to-income screen when enabled.
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
  the main forward distribution. During price iteration, the code solves the
  same Bellman problem and propagates the same full distribution, but computes
  only the market-clearing statistics needed to update price. Full reporting
  moments, including event-study housing responses, are recomputed once at the
  accepted price.
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
- `SIMPLIFICATION`: if enabled, the optional payment-to-income screen uses a
  simple interest-plus-tax payment. It is not Coven's full amortized mortgage
  payment formula.
- `SIMPLIFICATION`: the housing supply shifter `H0` and rent anchor `r_bar` are
  smoke-test normalizations, not calibrated targets.
- `SIMPLIFICATION`: old retention currently comes only through the inherited
  workhorse owner state, transaction/sale wedge, and bequest utility. A clean
  policy-created old-retention wedge is not yet implemented.
- `DIAGNOSTIC ONLY`: `estate_tax_rate` and `estate_tax_exemption` implement a
  narrow terminal bequest-tax wedge by reducing positive terminal resources
  before evaluating bequest utility. This is useful for proof-of-concept
  comparisons, but it is not a full estate-tax reform: there is no government
  budget, no rebate, no inheritance-transfer kernel, and no bequest-principal
  adding-up.

## Not Implemented

- `NOT IMPLEMENTED`: formal calibration, counterfactual tables, and production
  parameter search. The current random-search tool is diagnostic only.
- `NOT IMPLEMENTED`: full estate-tax counterfactuals with inheritance kernels,
  bequest-principal adding-up, and estate-revenue rebates. The implemented
  terminal bequest-tax wedge is diagnostic only.
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

- `DIAGNOSTIC ONLY`: 2026-06-08 calibration-structure ledger:
  `docs/model/intergen_housing_fertility_calibration_structure_20260608.md`.
  The note classifies the current parameter inventory, target list, missing
  multi-market moments, replacement non-location moments, and the smallest
  disciplined production calibration vector. Read it before launching another
  calibration search.
- `DIAGNOSTIC ONLY`: 2026-06-08 added
  `python -m intergen_housing_fertility.cli informed-smoke`, a deterministic
  parameter-led smoke panel. It varies internal economic blocks
  (`beta`, `alpha_cons`, `b_entry_fixed`, Stone-Geary terms, child costs,
  fertility taste/scale, housing needs, `chi`, and bequest terms), while
  keeping finance, supply, and menu objects fixed. It supports `--labels`,
  `--case-limit`, and `--fixed-price` so short local panels are recoverable.
- `DIAGNOSTIC ONLY`: 2026-06-08 informed-smoke records:
  `output/model/intergen_housing_fertility_informed_smoke_20260608/` contains
  three GE records before the interactive run was stopped for runtime;
  `output/model/intergen_housing_fertility_informed_smoke_fixedprice_20260608/`
  contains the six-case fixed-price sensitivity panel. The GE baseline has
  `tfr=1.430`, childlessness `0.381`, first-birth age `34.70`,
  prime-age ownership `0.312`, young liquid wealth/income `0.192`,
  old-age ownership `0.792`, old parent-childless gap `0.262`,
  \(H_{01}=0.604\), and \(H_{12}=0.227\). Patient `beta` raises wealth and
  ownership but did not clear tightly in the short GE panel; impatient `beta`
  clears and lowers wealth/ownership. In the fixed-price panel, higher
  `alpha`-housing-share/Stone-Geary pressure raises fertility and the
  first-child housing response but does not clear at the fixed price; higher
  `chi` hits prime-age ownership locally but drives old ownership near one and
  collapses the old parent-childless gap. Treat these as directional smoke
  evidence, not candidate calibrations.
- `DIAGNOSTIC ONLY`: 2026-06-08 added
  `python -m intergen_housing_fertility.cli local-panel`, a bounded multicore
  diagnostic panel. It varies only the implemented internal economic
  parameters in the current parameter ledger, keeps finance/supply/menu inputs
  fixed, ranks cases excluding `mean_age_first_birth`, checkpoints
  `cases.jsonl`, writes panel plots, and re-solves the best cases for the
  diagnostic packet. The first run
  `output/model/intergen_housing_fertility_local_panel_nb60_nz5_20260608/`
  used `Nb=60`, five Markov income states, six worker processes, and a
  30-minute submission budget. It completed `142` cases. The no-age-rank best
  (`case=33`) improves TFR and childlessness but is not a candidate
  calibration: `tfr=1.839`, childlessness `0.207`, prime-age ownership
  `0.747`, family ownership gap `0.377`, young liquid wealth/income `0.310`,
  old ownership `0.993`, old parent-childless gap `0.010`, \(H_{01}=0.641\),
  \(H_{12}=0.312\), and first-birth age `33.73`. The best full-loss cases
  are pulled toward very high fertility (`tfr` above `2.6`) and very low
  childlessness, confirming that the timing target is changing the search
  region rather than simply adding a harmless extra moment. Policy plots for
  the top records show sharp ownership transitions and nonmonotone
  childless-renter owner-entry policies; audit the tenure/finance block before
  treating these points as production calibration candidates.
- `DIAGNOSTIC ONLY`: 2026-06-08 added the target set
  `candidate_no_timing_v0` for a documented 13-moment trial. It excludes
  `mean_age_first_birth`, keeps parity composition diagnostic under the current
  `tfr = 2 * mean_completed_fertility` convention, and adds candidate targets
  for midlife liquid wealth/income (`1.20`), aggregate housing user-cost share
  (`0.24`), prime-age childless renter median rooms (`4.0`), and prime-age
  childless owner median rooms (`6.0`). These added moments are reactivated
  candidate targets, not a finalized empirical target system.
- `DIAGNOSTIC ONLY`: 2026-06-08 prepared
  `code/cluster/submit_intergen_housing_fertility_twohour_panel.sh` for a
  two-hour Torch panel test of `candidate_no_timing_v0`. The launcher uses one
  worker per array task, a 115-minute internal panel budget inside a 2:10 SLURM
  walltime, `J=16`, `Nb=60`, five income states, six owner rungs, and no
  diagnostic packets during the timed run. Task 1 includes deterministic anchor
  cases; other tasks use `--random-only` to avoid repeating those anchors.
  Collect outputs with
  `python tools/collect_intergen_panel_results.py --results-dir <RESULTS_DIR>`.
- `DIAGNOSTIC ONLY`: 2026-06-09 added
  `code/cluster/submit_intergen_housing_fertility_global_de.sh` and
  `python -m intergen_housing_fertility.cli global-de-panel`. This is a
  different global proposal algorithm, not a different model: each array task
  starts from broad Latin-hypercube bounds over the same 13 internal economic
  parameters and then applies differential-evolution proposals. It keeps
  `target_set=candidate_no_timing_v0`, `J=16`, `Nb=60`, five income states,
  six owner rungs, and the same checkpointed `cases.jsonl` output format.
- `DIAGNOSTIC ONLY`: 2026-06-09 final-best pathology audit:
  `docs/model/intergen_housing_fertility_pathology_audit_20260609.md`.
  It re-solves the final global-DE best under `Nb=60/120` and two owner-ladder
  choices. Main finding: the final best is not a production calibration;
  owner-entry policies remain highly nonmonotone, lifecycle ownership is
  wrong, and the temporary `linspace(2,10,6)` owner ladder materially
  contaminates the owner-room fit.
- `DIAGNOSTIC ONLY`: 2026-06-09 added
  `code/model/tools/run_intergen_policy_poc.py`, a fixed-theta policy
  comparison runner from the final global-DE toy best. The default cases are
  baseline, parent-targeted credit relief with `parent_dp_waiver_phi=0.95`,
  a property-tax increase from 1 percent to 2 percent annually with prices
  re-cleared, and a 30 percent terminal bequest-tax wedge. Output is in
  `output/model/intergen_policy_poc_20260609/`; the summary note is
  `docs/model/intergen_policy_poc_20260609.md`.
- `DIAGNOSTIC ONLY`: 2026-06-09 added
  `code/model/tools/audit_intergen_parent_credit_margin.py`, a treated-margin
  audit for the weak fertility response to parent-targeted credit relief. It
  evaluates feasibility and policy responses on baseline fertile childless
  renter states. Main finding: the first-child family-rung PTI screen is slack
  and LTV relief expands down-payment feasibility, but owner-entry and
  fertility probabilities barely move at the final toy theta. Output is in
  `output/model/intergen_parent_credit_margin_audit_20260609/`; the summary
  note is `docs/model/intergen_parent_credit_margin_audit_20260609.md`.
- `DIAGNOSTIC ONLY`: 2026-06-24 added
  `code/model/tools/build_intergen_mechanics_packet.py`, a non-production
  mechanics packet builder for saved intergen theta records. It preserves the
  active model logic, re-solves the one-market Markov-income model, and writes
  standard diagnostics, full target/model/gap tables, prime-age childless
  room-bin and owner-rung shares, lifecycle profiles, childless-renter
  owner-entry threshold tables, and optional fixed-theta policy proof-of-
  concept cases. Policy cases from this tool are mechanism diagnostics, not
  quantitative policy estimates.
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
- `DIAGNOSTIC ONLY`: 2026-06-08 status for job `10574658`: all `64` tasks
  finished and all `3,072` cases produced valid records. `2,868` cases cleared
  the market to residual \(\le 5\times 10^{-3}\), and `2,650` cleared to
  residual \(\le 10^{-4}\). The scalar-best converged case has loss `286.451`,
  residual \(4.57\times 10^{-5}\), `tfr=3.548`, childlessness `0.003`,
  ownership `0.937`, old ownership `1.000`, family ownership gap essentially
  zero, and first-birth age `30.48`. This is not an economically usable
  calibration point; the random search is still pulled toward a high-fertility,
  high-ownership corner. The best basic-filter candidate with
  `1.4<=tfr<=2.0`, `0.05<=childlessness<=0.30`, `0.40<=ownership<=0.75`, and
  `0.55<=old ownership<=0.90` has loss `480.752`, `tfr=1.962`,
  childlessness `0.269`, ownership `0.735`, old ownership `0.867`, family gap
  `0.529`, and first-birth age `32.18`; it is more interpretable but still far
  from the age-at-first-birth and wealth targets.
- `DIAGNOSTIC ONLY`: 2026-06-08 audit found and fixed a target-mapping bug in
  `calibration.py`: the old workhorse maps target `own_rate` to prime-age
  ownership `own_rate_3055`, but the one-market diagnostic extractor had been
  using aggregate ownership. The extractor now stores target `own_rate` from
  `own_rate_3055` and saves aggregate ownership separately as
  `aggregate_own_rate`. Re-ranking the completed Torch records with this fix
  does not rescue the scalar best; the dominant failure remains the
  first-birth-age target under the one-shot fertility architecture. Full audit:
  `docs/model/intergen_housing_fertility_audit_20260608.md`.
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
