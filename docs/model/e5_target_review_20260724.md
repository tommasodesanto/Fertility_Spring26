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
