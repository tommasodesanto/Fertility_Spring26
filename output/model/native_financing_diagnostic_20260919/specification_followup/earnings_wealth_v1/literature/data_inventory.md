# Earnings and household-income data inventory

Status: bounded provenance review, 2026-09-21. This memo is read-only evidence
for the earnings specification choice. It does not estimate a new process, run a
model, acquire data, or adopt a calibration.

## Operative model object

The quantitative model has four-year periods (`PERIOD_YEARS = 4.0`), starts at
age 18, retires at age 66, and has 17 model ages. The lifecycle and period
constants are in
[`calibration.py:21-27`](../../../../../../code/model/intergen_eqscale_seq_optimized/calibration.py:21)
and the active timing map is constructed in
[`calibration.py:1230-1245`](../../../../../../code/model/intergen_eqscale_seq_optimized/calibration.py:1230).

The model's resource income is a period flow. At working ages,
`P.income[i,j] = 4(1-tau_pay) w_i a_j` when period-flow scaling is on; the
income-state multiplier is applied in `income_at_state`. Retirement replaces it
with a balanced stationary pension. See
[`parameters.py:818-831`](../../../../../../code/model/intergen_eqscale_seq_optimized/parameters.py:818)
and [`parameters.py:834-859`](../../../../../../code/model/intergen_eqscale_seq_optimized/parameters.py:834).
The default payroll rate is `tau_pay = 0.179` and the pension is PAYGO-derived,
not imported from the PSID earnings variable
([`parameters.py:288-301`](../../../../../../code/model/intergen_eqscale_seq_optimized/parameters.py:288)).

Working-age model income is therefore after-payroll-tax and period-scaled. The
wealth targets use an annual gross denominator by dividing period income by four
and grossing up by `1/(1-tau_pay)`; the helper includes the model's lump-sum
property-tax transfer unless that transfer is explicitly removed. See
[`solver.py:327-335`](../../../../../../code/model/intergen_eqscale_seq_optimized/solver.py:327)
and the wealth-moment implementation's definition at
[`solver.py:6826-6833`](../../../../../../code/model/intergen_eqscale_seq_optimized/solver.py:6826).
This is a material contract: the empirical denominator must be gross labor
earnings if this model helper is retained, and a fiscal rebate must not be
silently described as labor earnings.

The age profile is an external piecewise profile
`[0.565, 0.838, 1.0, 0.985, 0.935]` over breaks `[18,25,35,45,55]`, normalized
to mean one over working ages. It is not estimated from the current earnings
packet. See [`parameters.py:862-881`](../../../../../../code/model/intergen_eqscale_seq_optimized/parameters.py:862).
There is no equivalence-scale adjustment in the model earnings state. The
current candidate has 5 persistent states crossed with 3 iid transitory nodes,
no permanent type, and a mean-one process; the constructor and its state
contract are in [`build_persistent_transitory_income_candidate.py:30-86`](../../../../../../code/model/tools/build_persistent_transitory_income_candidate.py:30).

## Existing US household-income builders

### July 27 PSID gross-labor packet: current candidate source

The active candidate uses the PSID shelf's `EARNINDRRC`, described in the source
builder as RP/spouse combined real 2022-dollar gross labor earnings. The narrow
extract keeps `ID`, calendar `year`, `RELTOHEAD_`, `AGEREP`, `DEATHYEAR`,
`EARNINDRRC`, and `IW`; it restricts to reference persons (`RELTOHEAD_ == 10`),
1984--2019, alive in the observation year, ages 25--60, positive earnings, and
positive survey weight. It checks one row per person-year. See
[`extract_psid_income_md.do:11-20`](../../../../../../code/data/psid_followup_mar2026/extract_psid_income_md.do:11)
and [`build_psid_income_between_within.R:95-119`](../../../../../../code/data/psid_followup_mar2026/build_psid_income_between_within.R:95).

The full fixed-effect/AR(1)/transitory packet confirms the concept, sample and
estimator: log earnings are residualized on integer age and year fixed effects
with `IW` weights; the process is

\[
u_{it}=a_i+p_{it}+e_{it},\qquad p_{it}=\rho p_{i,t-1}+\eta_{it}.
\]

It uses exact calendar-year autocovariances at lags
`0,1,2,4,6,8,10,12,16,20,24,28,32`, geometric-mean pair weights, and 199
person-cluster bootstrap draws with 10% covariance-correlation shrinkage. These
conventions are documented at
[`psid_income_fixed_effect_md_20260727/README.md:3-24`](../../../../../../code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/README.md:3)
and implemented at
[`build_psid_income_between_within.R:278-353`](../../../../../../code/data/psid_followup_mar2026/build_psid_income_between_within.R:278).

The unrestricted fixed-effect fit has objective `14.3178` (13 moments, 4 free
parameters). The no-fixed-effect nested fit has objective `39.5085` (13 moments,
3 free parameters), a gap of `25.1907`; the fit summary has the exact counts
`89,936` person-years and `12,508` persons
([`md_fit_summary.csv:1-3`](../../../../../../code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/md_fit_summary.csv:1)).
The fixed-effect fit estimates variance `0.3931`, persistent variance `0.3319`,
transitory variance `0.3098`, and annual `rho = 0.8863`, with fixed-effect
variance 38.0% of fitted residual variance
([`md_parameter_estimates.csv:1-11`](../../../../../../code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/md_parameter_estimates.csv:1)).

The current no-permanent-type candidate reuses the no-fixed nested fit directly:
annual `rho = 0.970303`, persistent variance `0.692827`, transitory variance
`0.344442`, and 15 states (5 x 3). Its metadata explicitly says that removing
the permanent component materially worsens fit and that the artifact is
diagnostic only
([`candidate.json:1-45`](../../../../../../output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json:1)).
The candidate's empirical selection is therefore **not** a proven coding bug;
it is a substantive restriction whose fit cost is measured in the same packet.

Important measurement limits:

- The positive-earnings restriction drops zero, missing, or nonpositive
  observations; there is no imputation or explicit unemployment/transfer state
  in this builder. This makes the process conditional on observed positive
  RP/spouse labor earnings.
- `EARNINDRRC` is gross labor earnings, not disposable household income. The
  builder has no tax, transfer, equivalence, or payroll adjustment. The model
  applies its own flat payroll wedge after importing the risk process
  ([`earnings_candidate/README.md:3-16`](../../../../../../output/model/native_financing_diagnostic_20260919/earnings_candidate/README.md:3)).
- The sample is a reference-person household concept, but there is no
  equivalence-scale transformation. Spouse aggregation changes the object from
  an individual wage to exogenous household labor earnings; it is not the same
  object as a one-earner wage process in Sommer, Floden--Linde, or a
  productivity process before tax.
- The candidate's calendar-year lag fit is not an exact four-year *sum* of
  annual household earnings. Its period transitory variance is
  `log(1+(exp(V_e,annual)-1)/4)`, matching the first two level moments of four
  independent mean-one iid shocks, while the persistent state is the four-year
  endpoint. The source README labels this a coarse approximation
  ([`earnings_candidate/README.md:18-35`](../../../../../../output/model/native_financing_diagnostic_20260919/earnings_candidate/README.md:18)).

### July 16 PSID `INCFAMR` process: separate, older family-income builder

`estimate_intergen_income_entry_targets.R` is a different pipeline. It uses
`INCFAMR`, described as annual family income, for reference persons aged 25--60,
1984--2019, with `INCFAMR > 1000`, alive-year and positive-`IW` filters. It
residualizes `log(INCFAMR)` on age dummies and year dummies, computes exact
calendar-year pair covariances at lags `1,2,4,6,8`, and fits persistent-plus-
transitory and AR(1)-only models by pair-count-weighted nonlinear least squares.
See [`estimate_intergen_income_entry_targets.R:154-273`](../../../../../../code/data/psid_followup_mar2026/estimate_intergen_income_entry_targets.R:154).

This process mixes annual and biennial observation regimes: lag-one pairs exist
only when the later observation is in the annual era (`second year <= 1997`),
while longer lags pool annual-era and biennial-era pairs and report the split
([`intergen_income_entry_targets_20260716/README.md:11-26`](../../../../../../code/data/psid_followup_mar2026/output/intergen_income_entry_targets_20260716/README.md:11)).
The current PSID gotcha is consequential: after 1997 PSID observations are
biennial, so a four-year average cannot be reconstructed by silently inserting
calendar-year values. This older builder does not do that insertion.

The saved `INCFAMR` estimates are very different from the current candidate:
`rho_annual = 0.97494`, `sigma_eta = 0.17688`, and `sigma_e = 0.41280`, based on
99,215 person-years and 13,028 persons with 499 bootstrap draws
([`block1_income_process_estimates.csv:1-16`](../../../../../../code/data/psid_followup_mar2026/output/intergen_income_entry_targets_20260716/block1_income_process_estimates.csv:1)).
The difference is not an automatic correction: `INCFAMR` is a broader family-
income object than `EARNINDRRC`, and the samples, moment lags, bootstrap design,
and source definitions differ. It must not be spliced into the current
RP/spouse-gross candidate without an explicit measurement decision.

### Existing PSID entry-wealth builder

The same July 16 pipeline constructs entrant wealth ratios from `NETWORTH2R /
INCFAMR`, childless reference persons, renters (`HOMEOWN == 2`), and positive
income. It reports weighted quintile-bin means and restricts the primary sample
to ages 18--24 because wealth supplements limit support
([`intergen_income_entry_targets_20260716/README.md:28-39`](../../../../../../code/data/psid_followup_mar2026/output/intergen_income_entry_targets_20260716/README.md:28)).
The active model turns those annual-income ratios into entrant liquid wealth by
multiplying them by model annual gross income at the entry state
([`solver.py:472-484`](../../../../../../code/model/intergen_eqscale_seq_optimized/solver.py:472)).
The wealth ratio's denominator is therefore currently `INCFAMR`, whereas the
earnings risk candidate is based on `EARNINDRRC`. This cross-file concept
inconsistency is a specification issue to close before a refit, not evidence
that either builder is individually broken.

## What is and is not available beyond PSID

The targeted active-code inventory found:

- `code/data/cps_fertility/` contains a CPS fertility-target builder and outputs,
  not a longitudinal CPS earnings-process estimator. CPS age profiles could
  validate cross-sectional levels or dispersion, but CPS repeated cross-sections
  cannot identify household serial covariance without an explicit rotation-panel
  design and compatible linkage.
- Active ACS builders under `code/data/Stata/` and
  `code/data/mms_center_periphery/` provide cross-sectional income, housing and
  geography objects. They can validate age profiles and household composition,
  but cannot identify a persistent/transitory serial process.
- No active SIPP earnings-process builder or local SIPP input was found in the
  targeted `code/data` inventory. No data-access workaround or acquisition was
  attempted.
- The prioritized published references are parameter/architecture checks rather
  than drop-in estimates. Sommer (2016, JME) supplies an annual age-profile plus
  persistent and iid wage-risk architecture (`rho=.95`, persistent innovation
  SD `.21`, iid SD `.17`); De Nardi (2004, ReStud) estimates after aggregating
  PSID income into complete five-year cells; Bick (2016, JEEA) constructs
  three-year period income before estimating; Boar--Gorea--Midrigan use a
  disposable-income concept with taxes and transfers. The project ledger records
  these distinctions and verifies that no prioritized paper supplied an exact
  four-year US stochastic-income implementation
  ([`ACTIVE_DECISION_LEDGER.md:246-273`](../../../../../../docs/model/ACTIVE_DECISION_LEDGER.md:246)).

Published values can validate broad persistence/risk ranges or motivate a
sensitivity arm. They cannot identify this model's household gross-versus-net
concept, spouse aggregation, missingness, age profile, or entry joint
distribution without remeasurement.

## Why the earlier PSID process was judged weak

There are three separate findings, and they should not be conflated.

1. **Measured fit weakness.** The no-fixed-effect nested fit is objectively much
   worse than the fixed-effect fit (`39.51` versus `14.32`), and the live
   no-permanent-type candidate documents that gap. Long lags provide information
   that a highly persistent AR(1) is standing in for permanent heterogeneity
   ([`psid_income_fixed_effect_md_20260727/README.md:9-24`](../../../../../../code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/README.md:9)).
   This is evidence against *that restriction on that sample*, not proof of a
   bad optimizer or a bad PSID file.
2. **Concept mismatch.** `EARNINDRRC` is positive RP/spouse gross labor earnings,
   while the model budget is after-payroll-tax household resources and some
   wealth denominators use `INCFAMR` family income. The July 22 reconciliation
   explicitly flags the open risk concept (pre-tax productivity versus after-tax
   household income), source choice, and state-count adequacy
   ([`eqscale_calibration_reconciliation_20260722.md:647-669`](../../../../../../docs/model/eqscale_calibration_reconciliation_20260722.md:647)).
3. **Historical feasibility interaction.** The July 17 memo found that realistic
   risk interacted with a positive Stone--Geary consumption/housing requirement,
   creating cash shortfalls unless a transfer, support restriction, or resource
   margin was supplied. That diagnosis concerns the old floor architecture; it
   does not show that the exact current persistent-plus-iid candidate fails
   ([`intergen_income_risk_feasibility_decision_memo_20260717.md:7-20`](../../../../../../docs/model/intergen_income_risk_feasibility_decision_memo_20260717.md:7)).

The code's annual-to-four-year AR(1) formulas are already mathematically correct:
`rho_4 = rho_a^4` and
`sigma^2_{eta,4}=sigma^2_{eta,a}\sum_{k=0}^3 rho_a^{2k}`. The iid conversion
matches the first two level moments of an average of four independent shocks.
The unresolved issue is that an endpoint persistent state multiplied by an
averaged iid proxy is not the exact block average along the persistent path. The
ledger records this as an approximation-validation problem, not an established
missing-conversion bug
([`ACTIVE_DECISION_LEDGER.md:220-240`](../../../../../../docs/model/ACTIVE_DECISION_LEDGER.md:220)).

## Feasible own-estimation and validation route

The lowest-cost defensible route is a two-stage data/model measurement check:

1. Freeze the income concept first: either the current RP/spouse gross labor
   object (`EARNINDRRC`) or a deliberately broader disposable/family-income
   object. State whether the model's flat payroll tax is the only tax treatment,
   whether transfers enter income, and whether the wealth-ratio denominator is
   changed to the same object. Do not combine `EARNINDRRC` risk with `INCFAMR`
   wealth denominators by default.
2. Using the existing narrow extract, form true four-year **sums** (or means with
   an explicit factor of four) only where four annual observations exist. Keep
   annual-era cells separate from later biennial endpoint pairs; never fill an
   unobserved year with zero or interpolation. At minimum report complete-cell
   counts, positive/missing shares, and the annual/biennial composition.
3. Estimate the same age/year residualization and person-cluster uncertainty on
   the four-year object. Compare the full covariance vector and tails, not only
   `rho` or a variance. For the model, generate synthetic observations, apply
   the identical empirical measurement operation, and score the same rows. This
   directly tests whether the candidate process survives its own data observer.
4. Keep a small validation matrix: (a) July 27 `EARNINDRRC` gross labor,
   (b) July 16 `INCFAMR` family income, and (c) a published annual reference
   such as Sommer, with each mapped consistently through the budget and tax
   contract. Report process choice, tails and entrant composition separately.
5. Hold entrant wealth fixed while testing the income process. The current
   income-grid diagnostic shows that changing income-grid resolution changes the
   entrant wealth distribution materially (L1 distances `0.355--0.396`), so a
   household-fit change is otherwise not attributable to earnings risk alone
   ([`income_grid_cohort_v2/comparison.md:3-13`](../../../../../../output/model/native_financing_diagnostic_20260919/specification_followup/income_grid_cohort_v2/comparison.md:3)).

Approximate local cost, without model solves: constructing four-year cells and
the point estimates from the existing narrow extract should be minutes; a
199-draw person bootstrap is likely on the order of tens of minutes, using the
existing checkpointed R design as the cost benchmark. The current full 499-draw
July 16 pipeline already exists and reports its estimates, so a first pass can
avoid repeating it. The direct four-year route is limited by the short annual
PSID era and by positive-earnings selection; later biennial data can validate
four-year endpoints but cannot identify four annual sums. If the annual support
is too small, retain the published-process arm as a clearly labeled sensitivity
rather than silently treating a proxy as an exact four-year estimate.

No calibration, regression, new data acquisition, or model run was performed in
this inventory.
