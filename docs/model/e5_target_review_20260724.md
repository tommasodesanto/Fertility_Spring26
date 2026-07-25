# E5 target review ledger (2026-07-24)

Row-by-row author review of the E5 target system. Each row is registered here
when signed off in discussion. Format: target, source, parameter connection.
This file becomes the E5 contract table when the review completes.

## Standing decisions (author-confirmed today)

- Fertility cohort anchor: women 40--44 in the June 2024 CPS fertility
  supplement, i.e. birth cohorts 1979--84 -- the youngest essentially-complete
  cohort. Mild completeness caveat (a few percent of births still ahead at
  40--44) disclosed in the note.
- Environment (prices, incomes, entrant wealth, housing menus) stays measured
  on 2019--2023 cross-sections; the seam between the fertility cohort and the
  modern environment is one disclosed assumption (preferences stable across
  adjacent cohorts).
- Top-bin convention: the CPS public file topcodes children at five, so the
  top-bin mean among 3+ families is the CAPPED mean (~3.602, builder to
  confirm). Adopted as measured, truncation disclosed; no uncapping
  adjustment. Replaces the provisional 3.4. The model's completed-fertility
  moment must use the same capped weight so model and data means share one
  convention.
- Stocks builder: `code/data/cps_fertility/` owned by the M agent
  (download+checksum manifest, sample definition, bootstrap SEs and
  covariance, window sensitivities). The E strand consumes the CSV.
- Timing rows: built by the E strand from NCHS cohort fertility tables
  (birth-certificate based) for the same 1979--84 cohorts.
- Period synthetic-cohort package (TFR ~1.62, synthetic childlessness, mafb
  ~27.3, share 30+) = declared robustness bracket at freeze time, not in E5.

## Registered rows

### Saving and bequests block [SIGNED OFF 2026-07-24]

- Aggregate wealth / gross labor earnings = 6.8731 -- primary row for beta.
  PROJECT-MEASURED AND COMPLETE: PSID 2005-2019, aggregate NETWORTHR ages
  18-85 over RP/spouse gross labor earnings (EARNINDRRC) ages 18-65, IW
  weights; 49,550 family-years, 10,432 persons; person-clustered bootstrap
  SE 0.3988 (interval [6.19, 7.69]). Documented in
  docs/model/intergen_wealth_target_beta_audit_20260723.md. The borrowed
  6.90 chain (De Nardi-Yang via Hendricks) stays retired; the numerical
  closeness is coincidental. NOTE for the E-package port: the matched model
  statistic uses the GROSS earnings denominator (no payroll wedge) -- port
  the current matched definition, not the July-23 after-tax version.
- Annual bequest flow / wealth = 0.0088 -- primary row for theta_0.
  Borrowed (Gale-Scholz via De Nardi-Yang), labeled as borrowed; model side
  is the at-death post-saving flow under repaired timing (B1 port).
- Late-life wealth dispersion p90/p50 = 3.4481 -- primary row for theta_1.
  Project-measured: PSID living reference persons 76-84, waves 1984-2019,
  bootstrap SE 0.1325 (499 person-clustered draws). Vintage footnote: the
  dispersion row pools 1984-2019 while the level row uses 2005-2019; the
  ratio is scale-free, so pooling is defensible, disclosed once.
- Registered caution from E4: with the flow unanchored, theta_0/theta_1
  collapsed and old-age ownership hit 0.945; these rows plus deterministic
  tenure are the guard. If old-age ownership still lands near 0.95 in E5,
  that is the health/long-term-care limitation conversation, not a weight
  tweak.
- A one-page unified wealth-target source ledger (M side) is endorsed; this
  entry cites the fragments meanwhile.

### Tenure and supply block [SIGNED OFF 2026-07-24]

- Ownership rate 30-55 (0.5755, ACS) -- primary row for chi_O, the owner
  utility premium. One-to-one.
- Aggregate occupied rooms (5.780, ACS) -- primary row for H-bar, the supply
  scale. One-to-one.
- kappa_T: DROPPED. Tenure choice is deterministic (the author's standing
  preference, reaffirmed; E4 drove kappa_T to ~0 at no cost; honest income
  risk generates tenure churn without taste noise). The Brier row dies with
  it. The four-year tenure switch rate becomes an untargeted validation
  check of the income-risk architecture; if it misses, that is evidence
  about the income process, not a dial to turn.
- Presentation rule (author): parameters fixed at zero (kappa_T, theta_n)
  do NOT appear in the assigned-parameters tables -- footnote treatment
  only ("tenure choice is deterministic"; "bequest utility does not vary
  with the number of children"). The strategy note's tables get one
  coherent revision pass when this review completes.
- Free parameters now: 10 (beta, delta_jump, delta_a, psi, kappa_E,
  kappa_C, chi_O, H-bar, theta_0, theta_1).

### Row 8 -- Consumption share alpha_0: FIXED EXTERNALLY [SIGNED OFF 2026-07-24]

- Decision: alpha_0 = 0.733 (CEX childless-cash-renter expenditure slope,
  2019-2023) becomes an external, leaving the free list: 11 free parameters.
  Wired as a contract default recorded in run metadata and verified by the
  collector -- NOT a CLI flag.
- Why this was controversial (May 2026 history, from the 05-28 transcript):
  (i) a May 27 override-flag implementation produced a record whose metadata
  alpha (0.93) did not match the effective alpha of its saved moments,
  invalidating a slide baseline; (ii) under Stone-Geary, alpha was a
  marginal share with no clean external counterpart (observed shares include
  the floors). Both are resolved: the E branch has no floors, so the
  childless-renter expenditure share equals 1 - alpha_0 exactly; and
  externals are now contract defaults with metadata checks.
- Cross-check kept for the note: E4's freely estimated alpha_0 = 0.754,
  three percent from the CEX value it never saw.

### Rows 6-7 -- Expenditure-tilts block [SIGNED OFF 2026-07-24]

- Row 6: first-child rooms increment, 0.664. PSID event-study response,
  horizon-0 diff-in-diff on the active Markov path. Identifies the SUM
  delta_jump + delta_a (a new parent gets the jump plus one per-child step).
- Row 7: rooms gap 3+ vs 1-2 children ages 30-55, 0.368, literal bins (L4).
  ACS build. The jump cancels within parents, so identifies delta_a alone;
  delta_jump is the residual of row 6. Triangular block, solved in order.
- delta_jump = one-time housing-share shift at parenthood; delta_a =
  per-additional-child shift. Registered tension: at E4 both rows overshoot
  at tiny tilts (scale e(n) already lifts family housing demand) while the
  family gap wants larger tilts; where E5 settles this tug-of-war tests
  whether the tilt device serves the space and tenure stories at once.

### Row 5 -- Family ownership gap [SIGNED OFF 2026-07-24]

- Target: 0.168, TARGETED in the hard loss. Parent minus childless ownership
  rate, ages 30--55. Model: same difference in the stationary population.
- Source: project-measured matched ACS build; SE status confirmed at the
  weights pass.
- Connection: no dedicated parameter -- the system's overidentifying
  restriction. For the draft: most closely connected to the composition of
  the childless margin (who stays childless, governed by the fertility
  noise scales) and the expenditure tilts jointly. Decision trail: July 23
  put it in the loss; a validation-only alternative was considered and
  rejected today because untargeted runs deliver 0.04--0.07 against 0.168.
  The overall hard-row count (aspiration: 13 rows for 12 parameters) is
  audited together at the end of this review.

### Fertility block summary [rows 1-4 complete]

Three parameters (psi, kappa_E, kappa_C) on four hard rows. Primary
connections for the draft: psi -- completed fertility; kappa_E -- first-birth
timing (mean and 30+ share); kappa_C -- childlessness. Overidentified by
roughly one moment (the timing pair is correlated, disclosed); the residual
pattern is a test of the sequential structure, and the E5 Jacobian reports
actual loadings. Birth-order progression rates stay untargeted diagnostics.

### Row 4 -- Share of first births at 30+ [SIGNED OFF 2026-07-24; value pending build]

- Target: pending the NCHS build; among mothers of the 1979--84 cohorts, the
  share with first birth at age 30 or later. Denominator = mothers (pure
  timing; no overlap with childlessness). Model: first births at ages 30+
  over all first births (computed; denominator convention verified at
  wiring).
- Source: same NCHS cohort-table builder as row 3; covariance included.
- Connection: kappa_E through the tail of the timing distribution -- stops
  the model from matching the mean with a wrong-shaped distribution.

### Row 3 -- Mean age at first birth [SIGNED OFF 2026-07-24; value pending build]

- Target: pending the NCHS build; mean age at first birth for the 1979--84
  birth cohorts (expected mid-26s). Model: mass-weighted mean age over
  realized first births in the stationary population (already computed).
- Source: NCHS cohort fertility tables (birth-certificate based;
  age-specific first-birth rates by mother's birth cohort; complete
  coverage by construction). Built by the E strand, committed script +
  manifest.
- Connection: mostly connected with kappa_E, the entry noise scale. Entry
  noise spreads the age at which women start trying and the fecundity
  schedule converts late starts into late or missing births, so the
  location of the first-birth age distribution is kappa_E's footprint.
  This row's absence is what let kappa_E reach its bound in E4.

### Row 2 -- Completed childlessness [SIGNED OFF 2026-07-24]

- Target: 0.188 (0.188180 reproduced). Share of women ages 40--44 with zero
  children ever born. Model: parity-0 share over post-fertile ages; the
  chosen/clock decomposition stays an untargeted diagnostic.
- Source: same CPS June 2024 builder as row 1; SE and covariance with row 1
  from the same bootstrap.
- Connection: mostly connected with kappa_C, the continuation noise scale,
  through the option value of entry: quiet continuation makes entering
  parenthood a committed multi-child track, which moves the share who never
  enter (in the v3 scan, kappa_C from 15.9 to 0.3 doubled childlessness at
  fixed entry noise). Draft states the primary connection; the fertility
  block remains jointly identified and the E5 Jacobian will report the
  actual loadings.

### Row 1 -- Completed fertility [SIGNED OFF 2026-07-24]

- Target: 1.918 (1.918425 reproduced from microdata). Mean children ever born
  per woman at ages 40--44. Model: mean parity over post-fertile ages, 3+ bin
  at the measured capped mean from the same builder.
- Source: CPS June 2024 supplement, women 40--44 (cohorts 1979--84); builder
  + SEs in progress (M agent).
- Connection: primary parameter psi (per-child utility flow). psi raises the
  value of every birth attempt, entry and continuation alike, so it shifts
  the level of completed childbearing; mean completed parity is the level
  statistic. Block-identified with (kappa_E, kappa_C) against the four
  fertility rows; completed fertility is psi's primary assignment.
  Cross-loads disclosed: fecundity (external), expenditure tilts, house
  prices.

## OPEN ITEM: top-group family-size convention (registered 2026-07-24)

The top family-size group carries a mean of 3.602 children. That number is
the CPS June 2024 mean among women with three or more children ever born,
computed on a public file that topcodes children ever born at five, so it is
a CAPPED mean. The model's completed-fertility moment applies the same cap,
so model and data means share one convention -- but both are then slightly
below the true mean of the group.

Author decision (2026-07-24): the note states the number and its source and
does not narrate the cap, which reads as confusing detail. VERIFY BEFORE THE
PAPER: (i) reconfirm the builder computes the capped mean as intended;
(ii) quantify how far the cap moves the mean (upper bound from an uncapped
external source, e.g. NSFG or the CPS internal file); (iii) decide whether to
correct or to keep the capped convention and disclose it in a footnote.

## OPEN ITEM: income-process reference (registered 2026-07-24)

Current external: persistence 0.9136 and innovation variance 0.0426 from
Floden and Linde (2001, RED 4(2), Table IV, GMM on PSID 1988-92 heads,
measurement error separately identified and excluded), with the innovation
scaled by (1 - tau), tau = 0.181, from Heathcote, Storesletten and Violante
(2017). Provenance is verified and recorded in
`code/model/intergen_eqscale_seq_optimized/externals.py`.

Author question (2026-07-24): Floden-Linde is a 2001 estimate on a 1988-92
panel and may not be the best modern reference. PREFERRED OPTION, raised by
the author: estimate the process ourselves on the PSID. The data are already
in the project for the wealth targets, and an own estimate would match the
model's exact sample and income definition (household earnings over our age
range) rather than borrowing an estimate on individual wages from a different
period. Fallback if we do not estimate: compare against Storesletten, Telmer
and Yaron (2004), Guvenen (2009), Kaplan and Violante (2010), Blundell,
Pistaferri and Preston (2008), and the processes used in the closest
structural housing papers (Sommer and Sullivan 2018; Sommer 2016). Also
decide whether the persistent / transitory split should be modeled rather
than folded into one persistent component. Any change here moves every wealth
moment, so it is a recalibration, not an edit.

INCONSISTENCY TO RESOLVE (flagged by the author 2026-07-24): the model uses
two different tax objects, and their values are nearly identical by
coincidence, which obscures the difference.
- `tau_pay = 0.179` is a flat AVERAGE RATE (OECD, via DUE). It sets the LEVEL
  of disposable income; `P.income` is stored net of it.
- `tau = 0.181` is the HSV PROGRESSIVITY parameter, a curvature. It is used
  only to scale the innovation s.d. by (1 - tau), i.e. to set how much RISK
  survives the tax system.
A consistent specification would apply the HSV function
`y_net = lambda * y^(1 - tau)` to both, calibrating lambda to reproduce the
average rate, instead of a flat rate for the level and HSV curvature for the
risk. Decide before the paper; this is a model change, hence a
recalibration.

## OPEN ITEM: fecundity constants are not fitted (registered 2026-07-24)

The fecundity schedule that gates every fertility decision in E5,
`1 - pi_a = 0.02 * exp(0.134 * (a - 18))` with a hard zero from age 45, is
hand-set. The code labels it "ILLUSTRATIVE fecundity schedule; not
calibrated" (`run_seq_smoke.py`). The FUNCTIONAL FORM is Sommer (2016, JME),
who fits it to Trussell-Wilson (1985) shares permanently infertile by age but
reports the fit graphically, not numerically -- so the constants could not be
copied and were chosen here. Leridon (2004/2005) is NOT our anchor; it
belongs to the de la Croix-Pommeret (2021) branch, which uses a bounded
logistic instead. Known discrepancy with Sommer: our schedule still gives
~0.35 conception probability at 44 and then jumps to zero, while his climbs
smoothly to full infertility at 45.

Consequence to keep in view: the two new timing targets (mean age at first
birth 25.31, share of first births at 30+ 0.270) are matched through this
schedule, so kappa_E's estimate is conditional on it. A refit that moves the
schedule may move the timing fit without any change in preferences.

Next step (not tonight): least-squares fit of (w1, w2) to the Trussell-Wilson
schedule, cross-checked against Menken-Trussell-Larsen (1986), documented in
`docs/model/fecundity_schedule_note_20260720.tex`, then a judgment on whether
the calibration must be re-run. The strategy note cites Sommer for the form
and does not display the constants as estimated.

## LAUNCH RECORD (2026-07-24 evening)

E5 launched on the signed system: smoke `14736164` (47/47 exact-loop
evaluations, arm E5, 10 free / 12 targets, timing rows live), production
array `14736281_[1-8]` (225-minute chains, strict twice-repeated winners),
dependent collector `14736282`. Seed: certified E4 winner with the
continuation scale carried over. Weight rule: (gap/SE)^2, measured SEs
where they exist, declared SEs elsewhere, documented in `e5_profile.py`.
Only the certified collector output is reportable.
