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

EVIDENCE FROM SOMMER (2016), the closest family paper (verified 2026-07-24
by reading the FEDS working-paper text, Section 4.1): Sommer does NOT use
Floden-Linde. He sets rho = 0.95 and the persistent innovation s.d. to 0.21,
with a separate transitory component of s.d. 0.17, citing Card (1991),
Hubbard, Skinner and Zeldes (1995), Storesletten, Telmer and Yaron (1998),
and Meghir and Pistaferri (2004); seven-state Tauchen for the persistent
component, two-state i.i.d. for the transitory. Note the consilience: our
GROSS innovation s.d. from Floden-Linde is 0.206, essentially Sommer's 0.21.
Our persistence 0.9136 is below his 0.95. We fold the transitory component
away entirely. (Same reading also confirms the fecundity attribution: Sommer
fits an exponential in age to Trussell-Wilson (1985) permanent-infertility
point estimates, with cumulative infertility set to 1 at 45 against 97
percent in the data.)

TWO TAX OBJECTS -- design rationale, not an accident. The July-22
reconciliation memo already states the intent: "the flat model payroll tax
shifts levels only and leaves log-income risk untouched; the HSV wedge stands
in for the progressivity the model does not have."
- `tau_pay = 0.179` is a flat AVERAGE RATE (OECD, via DUE); it sets the LEVEL
  of disposable income and `P.income` is stored net of it.
- `tau = 0.181` is the HSV PROGRESSIVITY parameter, a curvature; it scales
  the innovation s.d. by (1 - tau), i.e. it sets how much RISK survives.
The near-equality of the two numbers is coincidence and reads as confusing,
so the strategy note states the scaling and does not narrate the two roles.
STILL WORTH DECIDING before the paper: whether to replace the pair with a
single HSV function `y_net = lambda * y^(1 - tau)`, lambda calibrated to
reproduce the average rate. That is a model change, hence a recalibration.

## FECUNDITY FIT RUN 2026-07-24 (supersedes the "not fitted" framing below)

Corrected object definition. pi_a is a FOUR-YEAR CONCEPTION PROBABILITY, not
a permanent-infertility share. Because the model period is four years, the
July-20 memo's reading is the right one: Leridon (2004) four-year figures
(0.91 / 0.84 / 0.64 at ages 30 / 35 / 40) are DIRECT empirical counterparts
of pi_a. Sommer (2016) supplies the exponential FORM. The strategy note's
footnote now states this and reports the model values (0.90 / 0.80 / 0.62).

Least-squares fit to the three Leridon points (lead-run, numpy grid with a
closed-form step in w1):
  current (0.02,    0.134  )  SSR = 0.001788
  fitted  (0.01331, 0.14960)  SSR = 0.000189   (9.5x better)

DECISION-RELEVANT FINDING: the refit does NOT fix the E5 fertility misses.
The fitted schedule is slightly MORE generous at every age below 42 and
essentially IDENTICAL at 44 (0.3491 fitted vs 0.3482 current). Late-age
conception stays near 0.35 at 44 because Leridon's data really are that
generous (0.64 within four years at age 40). So the model's forgiving late
biology is empirically correct, not a calibration error, and the earlier
hypothesis -- that a proper fecundity fit would simultaneously repair the
late-first-birth, childlessness, and completed-fertility rows -- is WRONG
and is retracted here.

Consequence: the E5 fertility misses must be explained by the choice side,
not biology. See the E5 verdict section.

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

## E5B CERTIFIED (2026-07-26 morning) -- THE OPERATIVE WINNER

Jobs `14756965_[1-8]` + collector `14756966`; 6/8 chains eligible (two
failed the exact tight-repeat rule); winner loss 385.14 vs E5's 485.50 and
the restricted probe's 444.31. Report:
`output/model/eqscale_seq_e5b_recalibration_20260725/report/`. The E5
winner of 2026-07-24 is superseded for all reporting.

Winner (annual beta 0.9911): delta_jump 0.0239, delta_a 0.0122, psi 0.5850,
kappa_E 3.552, kappa_C 1.254, chi 1.0380, H0 8.1624, theta0 0.1684, theta1
0.0500, all interior (registered psi bound +-3 kept).

Fit: EIGHT of twelve rows within ~1.3 SE -- tfr 1.904 (-0.5 SE), mafb 25.64
(+1.3), own_rate 0.564 (-0.4), own_family_gap 0.166 (-0.2), rooms increment
(-0.2), 3+ rooms gap (+0.2), bequest flow (+0.6), aggregate rooms (+1.0).
The four failures are exactly the structural trio of the probe verdicts:
childlessness 0.080 (-14.1 SE; WORSE than E5, as predicted -- lower noise
removes the only childlessness device), share30+ 0.332 (+7.8, improved),
wealth 5.455 (-3.6, improved), p90/p50 2.068 (-10.4, unchanged). The
current architecture now fits everything it can fit; the remaining misses
are the mechanism questions in the 2026-07-25 memo, not calibration slack.

Policy packet at this winner launched as job `14796061` (--arm e5,
verification gate at 1e-6 against the certified table before any output).

## OVERNIGHT PROBE VERDICTS (2026-07-25 morning; both runs complete, 10/10 and 64/64)

Artifacts: `output/model/eqscale_seq_e5_probe_20260725/report/frontier.csv`
(collector defect noted: its `beta_annual` column holds the PERIOD beta;
annual = value^(1/4)) and
`output/model/eqscale_seq_e5_wealth_slice_20260725/task_*/wealth_slice.csv`.

1. THE E5 WINNER IS DOMINATED. Probe cell kappa_E = 4 reaches tight loss
   444.3 vs the certified E5 winner's 485.5 despite one FEWER free
   parameter. The E5 chains were stuck in the high-noise basin (echo of the
   single-basin finding of 2026-07-11). At kappa_E = 4: annual beta 0.9914,
   psi 0.62, kappa_C 0.63, chi 1.04, ownership 0.589 (target 0.575), wealth
   5.58, tfr 2.015, childlessness 0.096, mafb 26.1, share30+ 0.361.

2. TIMING IS BIMODAL IN kappa_E. At kappa_E <= 1 the model collapses to
   enter-at-18 (mafb 18.2, childlessness 0.000, tfr 2.72): the DETERMINISTIC
   economics favors immediate entry, so all observed postponement at larger
   kappa_E is noise, not prices. Between kappa_E = 1 and 2 the solution jumps
   straight from degenerate-early to over-dispersed (mafb 26+, share30+
   0.36+). NOWHERE on the frontier is the data pair (25.31, 0.270)
   approached: whenever the mean is near 26, the 30+ share is ~0.36. Under
   logit entry + current fecundity, the timing PAIR looks structurally
   unmatchable, and the mechanism warning stands: at current parameters the
   down-payment channel does not delay births; noise does.

3. CHILDLESSNESS HAS A CEILING OF ~0.118 on the entire frontier (target
   0.188). With entry noise the only childlessness device, the model cannot
   come close; this is now the binding fertility failure and points to a
   mechanism (e.g. permanent heterogeneity in psi or a career margin), not a
   dial.

4. THE PSI CORNER DISSOLVES: with the bound widened to 6, no cell pins
   (max |psi| = 3.21 at kappa_E = 36.3; 0.33-0.92 at low kappa_E). The E5
   cap at 3.0 was binding but is not the deep problem.

5. WEALTH 6.87 IS REACHABLE, OWNERSHIP VIA beta/theta0 IS NOT. Slice: wealth
   crosses 6.87 only at annual beta >= ~0.99 (6.93 at beta_a 0.9981/theta0
   0.5; 7.39 at theta0 1; 8.03 at theta0 2), consistent with the M-side
   beta-r tension. But ownership never exceeds 0.447 anywhere on the 64-node
   grid (chi and fertility fixed at the E5 winner), while the kappa_E = 4
   probe cell (chi re-optimized, beta_a 0.991) puts ownership at 0.589: the
   ownership target needs the discount factor AND chi jointly, not beta
   alone.

6. p90/p50 IS FLAT AT ~2.0-2.2 EVERYWHERE in both runs (target 3.45), even
   at theta0 = 8 where the flow overshoots (0.022). theta1 is not generating
   dispersion at any point visited: the late-life tail needs a mechanism
   (echoes the old "theta1 inert" finding), not reweighting.

DECISIONS FOR THE AUTHOR (none taken): (a) an E5b re-run of the full
10-parameter system seeded at the kappa_E = 4 probe winner, same contract,
to replace the dominated winner; (b) the parked weights discussion now has
content: the timing pair and the childlessness level appear structurally
unmatchable under the current architecture, so the choice is between
retargeting (e.g. mean-age only) and a mechanism change (entry
heterogeneity); (c) a policy packet at the kappa_E = 4 point to measure
whether prices move births at all once noise is moderate.

## OVERNIGHT PROBES LAUNCHED (2026-07-24 ~23:45, author-approved "proceed with all")

Motivated by the E5 verdict (psi at its 3.0 cap with fertility overshooting;
kappa_E only 49->36; wealth 4.35 vs 6.87 with ownership 0.41). Target system
and weights UNCHANGED in both runs; these are diagnostics, not calibrations.

1. kappa_E-profile probe, job `14741878_[1-10]` (cpu_short, 3:55): fixes
   kappa_fert per cell at {0.5, 1, 2, 4, 8, 12, 18, 24, 30, 36.339}, widens
   psi_child to [-6, 6] (E5_PSI_BOUND=6), re-optimizes the other nine
   parameters from the certified E5 winner seed, 200-minute chains, strict
   twice-repeated tight winners. Arm `E5_PROBE_KE`, 9 free / 12 targets.
   Smoke: 13/13 strict exact-loop evaluations, metadata verified.
   Collector: `collect_e5_probe.py` on
   `output/model/eqscale_seq_e5_probe_20260725/production` (run at collection).
   Questions: is there an interior kappa_E that reaches the timing rows, at
   what cost elsewhere, and does freeing psi resolve the fertility tension?

2. beta-theta0 wealth slice, job `14741879_[1-8]` (cpu_short, 2:00): 64
   tight solves (Nb=120, 40/2.5e-5) at the E5 winner, beta_annual in
   {0.94..0.9981} x theta0 in {0.05..8}, all else fixed; per-node CSV
   checkpoints, resume-safe. Smoke: 2 nodes + resume-skip verified.
   Output: `output/model/eqscale_seq_e5_wealth_slice_20260725/task_*/`.
   Question: is wealth/earnings 6.87 reachable at all in this structure, and
   what happens to ownership and the bequest rows where it is closest?

Infrastructure commit `01f5cb9`; default paths bit-identical (suite 123
passed with the gates off).

## LAUNCH RECORD (2026-07-24 evening)

E5 launched on the signed system: smoke `14736164` (47/47 exact-loop
evaluations, arm E5, 10 free / 12 targets, timing rows live), production
array `14736281_[1-8]` (225-minute chains, strict twice-repeated winners),
dependent collector `14736282`. Seed: certified E4 winner with the
continuation scale carried over. Weight rule: (gap/SE)^2, measured SEs
where they exist, declared SEs elsewhere, documented in `e5_profile.py`.
Only the certified collector output is reportable.

## E6A FECUNDITY-TAIL GATE (2026-07-27)

The handoff's cheapest-first step is complete at the certified E5b parameter
vector. The existing two-parameter fit to the Leridon four-year conception
probabilities does not repair timing: it moves the 38+ first-birth share only
from `0.0963` to `0.0942` and raises the signed loss from `385.14` to
`402.27`.

A default-off terminal-decay gate was therefore added to the shared fecundity
function. It preserves the fitted conception schedule through age 40 and
decays success from 40 to the hard close at 45. The preferred external reaches
`0.03` immediately before age 45, matching the terminal evidence used in the
handoff; 5 and 10 percent terminal-success schedules are retained as
sensitivities. The gate is used identically by the household problem and
population flow, adds no free calibration parameter, and leaves the default
path bitwise unchanged.

Strict fixed-winner results (all market residuals below `2.0e-5`):

| Schedule | Signed loss | Mean first-birth age | Share 30+ | Share 38+ | Share 42+ |
|---|---:|---:|---:|---:|---:|
| E5b current | 385.14 | 25.637 | 0.332 | 0.0963 | 0.0348 |
| Leridon two-parameter | 402.27 | 25.597 | 0.330 | 0.0942 | 0.0336 |
| Terminal success 10% | 355.72 | 25.397 | 0.320 | 0.0851 | 0.0254 |
| Terminal success 5% | 338.19 | 25.272 | 0.314 | 0.0792 | 0.0199 |
| Terminal success 3% | **331.40** | 25.196 | 0.311 | 0.0756 | 0.0165 |

Verdict: the terminal tail clears the E6a refit gate, but biology alone does
not generate the missing 25--30 hump. At fixed parameters the 3 percent
schedule improves timing and childlessness but moves the family ownership gap
from `0.166` to `0.137`; the full refit must determine whether that cost is
recoverable. Exact outputs:
`output/model/eqscale_seq_e6a_fecundity_tail_frontier_20260727/`.

## E6A REFIT LAUNCH (2026-07-27)

The exact-loop Torch smoke `14843503_[1-2]` completed successfully. Each
chain produced 13/13 strict cases and identified arm `E6A`, the unchanged
signed twelve-row target system, ten free parameters, and the three standing
external restrictions. Both the household problem and population flow used
the registered terminal-tail schedule.

Production array `14843676_[1-8]` is now running with a 225-minute and
1,000-case cap per chain. Chains 1--2 start at the certified E5b winner,
chains 3--5 use `start_mix = 0.10`, and chains 6--8 use
`start_mix = 0.25`. Dependent collector `14843677` runs only after all eight
array tasks succeed. Based on the smoke throughput, the planned maximum is
about 21.6 CPU-hours and 3.75 wall-clock hours under full concurrency. Only
the collector's twice-repeated tight winner will be treated as reportable.
