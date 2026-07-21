# Intergen Numerical And Economics Issues Ledger

Date: 2026-06-25; urgent update 2026-07-21

This is a working ledger of issues to think about in the one-market
intergenerational housing-fertility model. It is not a task checklist and not a
calibration target sheet. The goal is to organize open questions so the model
can be made economically well behaved before serious calibration or policy
work.

## Reserved Economics Question: Normality or Inferiority of Children

**Status: open; expand later.** Determine whether children are normal or
inferior in the current model, using a precise income/wealth comparative
static for fertility rather than informal language. Separate the response of
completed fertility from birth timing, and distinguish a genuine preference
income effect from changes operating through housing access, tenure, prices,
borrowing constraints, and the bequest motive. Identify the parameter and
state-space conditions under which the sign can change, then test the result
in both partial and general equilibrium.

## Absolute Urgency: Funded Property-Tax Baseline

**Priority: P0, for July 21--22.** The current calibrated baseline charges the
1% property tax but does not return the revenue. A fixed-M5-theta test in the
parity-verified optimized solver closes the stationary government budget with
an equal lump-sum rebate. The unrebated M5 TFR is `2.0346`; rebating the
baseline 1% revenue raises it to `2.3485` (`+15.43%`). Relative to that funded
baseline, a rebated 2% tax raises TFR by `8.04%` and lowers the house price by
`17.41%`; using part of the revenue for the targeted `0.4` purchase grant gives
`+7.98%` and `-17.17%`. This is strong fixed-parameter evidence that returning
the previously burned tax revenue increases fertility in the current model.

It is not yet a paper-ready policy result. Fiscal closure moves the fixed-theta
loss from about `9` to `219.7`, so the existing calibration cannot simply be
relabelled as funded. Before relying on or updating the circulated policy
numbers: (1) move the tested fiscal closure into the production model; (2)
recalibrate under the rebated 1% baseline; (3) rerun the rebated 2% tax and the
funded grant package; and (4) repeat the population-adjusted estimand used in
the paper. Exact test artifacts are in
`output/model/intergen_funded_property_tax_test_20260721/`; the driver is
`code/model/tools/run_intergen_funded_property_tax_test.py`.

## Current State

The active model is `code/model/intergen_housing_fertility/`, the June 2026
one-market/no-location strand.

Recent code corrections:

- PTI is off by default. The active borrowing restriction is collateral/down
  payment only.
- Optional PTI was fixed so that, if enabled later, it uses actual transaction
  debt rather than max-LTV debt.
- `chi` was corrected to be an owner utility/service premium on residual owner
  housing services after physical family-space needs. It is not a physical room
  multiplier.
- First-look plots and density plots were corrected to distinguish policy
  branches, realized choices, liquid wealth, total wealth, tenure, and child
  status more cleanly.

The latest 4-hour Torch run completed cleanly with 12,398 valid records under
`candidate_replacement_young_old_roomgap_v1`, but it did not solve the
economics.

## Main Frontier Fact

The model can hit pieces of the target system, but not jointly.

- Room separation and owner room levels come with almost no young owners.
- Young ownership plus lower old ownership collapses the owner-renter room gap.
- High young ownership usually keeps old ownership too high.
- The latest run found zero records with young ownership at least 0.25, old
  ownership at most 0.85, and owner-renter room gap at least 1.5.

Interpretation: this is probably not just insufficient search. It is evidence
of a mechanism or model-object problem.

## Things To Think About

### 1. Math-Code Alignment

- Is every object in code exactly the object in the draft/slides?
- Are owner housing services, physical rooms, and utility premium separated
  correctly?
- Does the down-payment/collateral constraint implement the intended rule: a
  house costing `pH` can be financed up to `phi * pH`, with cash need
  `(1 - phi) * pH`?
- If PTI exists as an optional feature, is it cleanly off by default and
  economically interpretable when on?
- Are tenure taste shocks present where we think they are, and are they
  entering the intended choice margin?

### 2. Lifecycle Ownership Mechanism

- Why does the model struggle to get young ownership without also generating
  very high old ownership?
- Is old ownership mainly inherited mechanically from prime-age ownership
  rather than separately disciplined by bequest/retention motives?
- Do we need an explicit downsizing, maintenance, moving, health, liquidity, or
  old-age housing-adjustment mechanism?
- Is young access controlled too much by the same objects that control old
  retention and room separation?

### 3. Owner-Renter Room Separation

- Why does the owner-renter room gap appear only when young ownership is nearly
  zero or ownership composition becomes distorted?
- Is `chi` doing too many jobs: ownership utility, owner-renter space
  separation, young access, and old retention?
- Do we need separate primitives for tenure utility and owner size/services?
- Are owner rungs and renter continuous housing comparable empirical objects?
- Are renters too large because the rental cap/menu/support is too loose, or
  because housing preferences/floors push them there?

### 4. Housing Ladder And Grid Use

- Are owner rungs being used in a sensible way, or are choices bunching on too
  few rungs?
- Is the upper rung unused because it is genuinely unattractive or because of
  grid/support/price artifacts?
- Is the wealth grid mostly empty relative to the ergodic distribution?
- Should policy plots be read against liquid wealth, total liquidated wealth,
  or both?
- Are apparent nonmonotonicities economic, or numerical artifacts from
  interpolation, grid coarseness, or branch-switching?

### 5. Fertility-Housing Interaction

- Why do some points with better ownership profiles have poor fertility or
  excessive childlessness?
- Are child utility, child goods costs, and family-space floors jointly
  overloading the same variation?
- Does the one-shot completed fertility architecture map cleanly enough to the
  housing-response moments?
- Are housing increments for first and second child being interpreted as
  housing demand moments, not sequential birth hazards?

### 6. Ergodic Mass And Policy Interpretation

- Policy plots should always be read together with where agents actually are in
  the stationary distribution.
- Aggregate density should come first, then useful splits by tenure and parent
  status.
- Conditional policy branches should not be confused with realized active
  housing after tenure choice.
- Representative policy plots should include both low/high income states and a
  small number of ages, rather than every slice at once.

### 7. Target-Object And Measurement Questions

These are not the main focus today, but they matter because some apparent model
failures may be target-object failures.

- Is the old-age ownership target the right object for bequest/retention, or
  should old wealth, downsizing, or retention be used instead?
- Are renter and owner room targets measured consistently with model room
  services?
- Is the owner-renter room gap better than owner median rooms for this
  discrete-rung model?
- Is aggregate housing user-cost share comparable to model housing cost
  objects?
- Which moments identify `beta`, `theta0`, `theta_n`, `chi`, and `h_bar`
  objects cleanly?

### 8. Numerical Reliability And Run Workflow

- We should not rerun the model just to regenerate plots when a saved solution
  object exists.
- Every serious run should produce a standard readout: status, stderr, record
  count, best scalar point, frontier screens, top loss contributions, and
  representative records.
- Plots should be reproducible from saved result records and solution caches.
- Run speed should be logged with enough detail to separate laptop throttling,
  Slurm variation, compilation, and genuine solver difficulty.
- A frontier ledger is more useful than only best-loss tables.

## Representative Point Types To Keep

For each major run, keep a small set of representative records:

- Best scalar-loss point.
- Best high-room-gap point.
- Best young-access point.
- Best low-old-ownership point.
- Best soft-joint point, even if the scalar loss is bad.
- Any pathological but informative corner point, such as nearly universal
  ownership or nearly zero ownership.

For each representative point, generate the same diagnostic packet so
comparisons are visual and mechanical rather than ad hoc.

## Current Interpretation

The current evidence points toward a missing or conflated mechanism in the
tenure lifecycle block. The model needs to separately control young owner
access, old owner retention/exit, and owner-renter space separation. Right now,
movement that helps one of these margins tends to damage another.

Before another broad calibration search, use the existing frontier to decide
whether to adjust mechanisms, target objects, or diagnostics.

## Audit Update: Concrete Issue Ledger

The June 25 objective-level audit sharpened the ledger. The same-renter
consumption and housing drops are real objects in the saved policy functions,
but the dense \(b'\) probe does not support the interpretation that the
golden-section optimizer is missing the global optimum. The current bottleneck
is therefore less "rewrite the Bellman optimizer" and more "separate the
economic mechanisms and clean up the numerical support."

| Priority | Issue | Current Read | Why It Matters |
|---|---|---|---|
| `P0` | Same-renter \(c,h\) dips | Dense \(b'\) audit finds the stored renter policy is essentially globally optimal. \(b'(b)\) is monotone, but the marginal saving rate can exceed one, so current \(c\) and \(h\) fall as wealth rises. | Need to distinguish genuine threshold-saving toward lumpy tenure/fertility choices from grid-amplified kinks before reading policy plots structurally. |
| `P0` | Young ownership near zero | `own_rate_2534` is about `0.0005` versus target `0.341`. This appears to be a value-ranking problem, not just affordability: young households can afford small owner rungs but prefer renting. | The model lacks a robust young-owner advantage under the current user-cost rent closure. |
| `P0` | Old ownership too high | Old-age ownership is about `0.916` versus target `0.764`. Owner status is too sticky through sale costs, owner service premium, user-cost wedge, gross-estate bequest, and no old-age exit channel. | The lifecycle tenure profile is not economically well behaved. |
| `P1` | Renter cap versus room targets | `hR_max=6` helps create owner-renter room separation, but makes renter share with rooms \(\ge 6\) nearly mechanically unreachable. | The current support can fit the owner-renter room gap and fail the large-renter target by construction. |
| `P1` | Bequest child-level penalty | With \(\sigma \simeq 2\), child-weighted CRRA bequest utility lowers terminal utility at fixed estate for households with more children. | This is a cardinal-normalization problem that can distort fertility and saving incentives. |
| `P1` | Grid and feasibility hygiene | Large unreachable grid regions and \(-10^{10}\) penalty cliffs are latent interpolation hazards. Owner residual services are floored at `1e-10` rather than treated as explicit infeasibility. | Could contaminate value interpolation even if current headline policies are not optimizer failures. |
| `P2` | Fertility interpretation | Fertility is one-shot completed family size, not sequential parity progression. | Do not interpret \(1\to2\) objects as literal second-birth hazards without changing the state space. |
| `P2` | Documentation mismatch | `CALIBRATION_STATUS.md` and `code/model/intergen_housing_fertility/README.md` are current. Parent `code/model/README.md` still describes the older center-periphery path. | Future work can start from stale orientation unless this is cleaned up. |

## Audit Update: Diagnostics And Candidate Fixes

This is the current ordered to-do list. The first column is the diagnostic to
run; the second column is the class of solutions it can discipline. These are
not calibration searches.

| Priority | Diagnostic | Potential Solution |
|---|---|---|
| `1` | Fixed-price grid refinement at \(p=0.8144\), with `Nb` in `{60,120,240}`. Compare same-renter \(b'\), \(c\), and \(h\) slices. | If dips shrink, refit the wealth grid around reachable support. If stable, treat the dips as real threshold-saving economics. |
| `2` | Young owner value-gap table for renters at ages 26--34 across wealth, owner rungs, and `chi`. | Add a genuine young-owner advantage: rental markup, tax wedge, expected capital gains, starter-rung redesign, or early parental down-payment transfers. |
| `3` | Old owner retention decomposition: sale cost, `chi`, user-cost wedge, bequest convention, and gross-versus-net estate. | Add an old-age exit, downsizing, maintenance, or liquidity mechanism; test the net-estate bequest convention. |
| `4` | Renter cap sweep: `hR_max` in `{6,8,10,inf}` at fixed theta. | Decide whether the large-rental target is off-model, or replace the hard cap with economic rental supply or price wedges. |
| `5` | Bequest ablations: `theta_n=0`, zero-estate-normalized bequest, and net estate \(b+(1-\psi)pH\). | Normalize bequest utility so children affect marginal bequest motives without adding a fixed terminal utility penalty. |
| `6` | Feasibility-mask audit for owner residual services and infeasible continuation values. | Replace the `1e-10` owner service floor and \(-10^{10}\) cliffs with explicit masks; rerun a moment-invariance check. |
| `7` | Tenure-smoothing sweep: `kappa_t` in `{0,0.005,0.01,0.02,0.05}` with dense excess-demand curves. | Document whether market clearing relies on smoothing and decide whether smoothing is structural or numerical. |
| `8` | Budget and mass-conservation invariant audit on occupied branches. | Fix any residual above numerical tolerance before calibration. |
| `9` | Documentation cleanup for the active strand. | Update parent `code/model/README.md` or add a prominent pointer to `intergen_housing_fertility/README.md`. |

Recommended first move: run the fixed-price grid-refinement diagnostic. It
does not change primitives, does not launch a calibration search, and directly
resolves whether the same-renter policy dips are economically stable or
grid-amplified.

## Parallel Diagnostic Sprint Results

The parallel diagnostic sprint is summarized in
`docs/model/intergen_parallel_debug_synthesis_20260625.md`. The underlying
diagnostic outputs are under
`output/model/intergen_parallel_audits_20260625/`.

Main updates:

- Young ownership is primarily a value-ranking failure, not a static
  down-payment-feasibility failure. Owner rungs are usually feasible for young
  childless renter-origin mass, but renting dominates in branch value.
- Old ownership is near-absorbing same-rung ownership. Among old owners,
  roughly `0.976` stay on the same owner rung next period, while only `0.002`
  sell to rent and `0.004` owner-downsize.
- Same-renter \(c,h\) dips are not just a coarse-wealth-grid artifact. Fixed
  price refinement through `Nb=240` reduces the flagged age-30 dip but does not
  eliminate dips in high-mass slices. \(b'(b)\) remains monotone.
- The renter \(h_R\ge 6\) moment has a likely definition problem in the Markov
  income extractor. Full 7D cached distribution mass gives weak
  \(h_R\ge6\) around `0.124`, close to the `0.138` target, while the saved
  implemented moment reports `0.013` after income-state policy collapse.
- The owner ladder has no clean starter-owner product: 2 rooms is utility-dead,
  4 rooms is a weak/floor rung for family states, and 6 rooms is effectively
  the first robust owner product. But starter-rung probes alone do not fix
  young ownership.
- Bequest utility has a confirmed child-level penalty under the current
  \(\sigma=2\) CRRA form, and terminal bequests use gross \(b+pH\) rather than
  net liquidation value \(b+(1-\psi)pH\).
- Accounting and KFE invariants are mostly clean, but the wealth grid is
  overwide and purchase/down-payment cliffs are reachable near occupied mass.
- Tenure smoothing affects choices as a real logit mixture, but local market
  clearing at the saved price does not mechanically depend on smoothing:
  hard-argmax \( \kappa_t=0 \) nearly clears locally.
- Fertility accounting is internally consistent as one-shot completed family
  size with `tfr = 2 * mean_completed_fertility`; the main risk is terminology
  that invites a sequential hazard interpretation.

Updated first moves:

1. Audit and, if needed, fix the renter large-room moment definition before
   changing primitives or re-running calibration.
2. Inspect the renter Bellman/intratemporal kink regions behind persistent
   same-renter \(c,h\) dips.
3. Decompose young owner value gaps and old owner retention as separate
   mechanism failures.
4. Only then test structural changes: young-entry transfer or rent wedge,
   old-age exit/downsize shock, large-rental price wedge, bequest
   normalization, and starter-owner menu redesign.
