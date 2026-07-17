# M5 recalibration contract (drafted 2026-07-16, pre-launch)

Purpose: recalibrate the one-market intergenerational model conditional on
externally disciplined inputs, replacing the three provisional constructions
identified by the M4 audit (`docs/model/intergen_m4_calibration_audit_20260716.md`):
the reverse-engineered income process, the circular entrant-wealth condition,
and the grid-discrete nonhousing-median target. M4
(`output/model/intergen_standard_bequest_recalibration_20260716/final_report/`)
is the warm start and comparison baseline, not the paper object.

Status legend: [PENDING-DATA] = filled from
`code/data/psid_followup_mar2026/output/intergen_income_entry_targets_20260716/`
before wiring; [DECISION] = Tommaso's call before launch.

## External inputs (new)

1. **Income process** — persistent component of a persistent+transitory
   decomposition estimated from PSID residual log family income (reference
   persons 25–60, 1984–2019, age and year effects removed), fit by weighted
   minimum distance on the autocovariance vector, person-bootstrap SEs.
   Annual (rho, sigma_eta) = [PENDING-DATA Block 1]. Fed through the existing
   verified conversion (`local_panel.py:1049-1094`: rho^4, time-aggregated
   innovation, 5-state Rouwenhorst). The single-AR(1) restricted fit is
   reported alongside for the paper. [DECISION D1 default: persistent-only
   feed.] The matched constants (0.9602, 0.0645) are retired everywhere; the
   `local_panel` module defaults and driver constants collapse to one source.
2. **Entrant wealth at age 18** — quintile-bin ratio distribution for
   childless renters aged 18–24 (fallback 18–26 if < 500 family-years),
   same construction as the retired 25–35 object, honest provenance string.
   Nodes/weights = [PENDING-DATA Block 2]. The 25–35 mean (0.17922556)
   remains a hard target — now a genuine lifecycle test, circularity removed.

## Target system (v2)

14 moments / 13 free parameters (fork: 15/14 if D4 frees tenure kappa):

- The 12 established moments, unchanged values and weights.
- `old_total_estate_wealth_to_annual_income_median_7684` = 6.50131577436537,
  weight 18.585767349158665 (unchanged; identifies theta0 — nested-zero
  evidence: dropping theta0 costs 303 loss points through this row alone).
- **NEW** `old_nonhousing_ge_1x_income_share_6575` = [PENDING-DATA Block 3],
  weight 1/SE^2 — replaces the grid-discrete nonhousing median (Jacobian row
  exactly zero; node spacing ~3 SEs). The median stays as a diagnostic.
- [DECISION D4 fork] If tenure kappa is freed: add
  `old_age_own_rate_6575` = [PENDING-DATA Block 4] with its fresh person-
  bootstrap 1/SE^2 weight as the 15th moment. Rationale: kappa's identifying
  variation is the ownership age path; the moment was demoted at M1 only
  conditional on kappa=0.

Untargeted diagnostics reported every readout: estate p90/p50, family-size
estate gap, housing/nonhousing decomposition 76–84, ownership path 62–78,
young liquid wealth by age, nonhousing median 65–75.

## Parameters

- Free (13): the 11 clean-frontier coordinates (M4 bounds/transforms
  unchanged) + theta0 [0, 8] softzero + theta1 [0.02, 16] log.
- External: theta_n = 0 (Hurd 1989; Kopczuk–Lupton 2007; DFJM 2025; own PSID
  gap 0.101, SE 0.563); tenure_choice_kappa = 0 unless D4 frees it;
  `normalize_bequest_utility = True`; SSA survival; no owner-LTV taper;
  standing conventions (phi=0.80, psi=0.06, sigma=2, q, delta, eta_supply).
- The complete free/fixed table is shown to Tommaso before launch. No
  restriction travels silently inside a phrase.

## Search design

- 8 chains, walltime 3:55 (cpu_short cap), <= 1,000 evals each — ~3x M4's
  effective budget; search evaluator (10, 1e-4); two tight repeats
  (40, 2.5e-5), bit-identical required.
- Warm starts: M4 winner; M1 winner; exact theta0=0 nested seed (collector
  fails if the free search ends above it); dispersed theta1 starts
  {0.1, 0.25, 0.42, 1.0, 4.0, 12.0} — the sub-0.42 region M4 never explored.
- Dependent jobs behind the strict collector:
  (a) two-sided theta1 profile, predeclared ladder spanning [0.1, 2.0], other
      coordinates re-optimized per cell — the interior-optimum question is
      testable this time;
  (b) tenure-kappa fixed-theta sweep {0, 0.01, 0.02, 0.05} at the M5 winner
      (skipped if D4 already freed kappa);
  (c) fresh Jacobian at the winner.
- Winner verified at Nb=240 before any promotion (standing convention that M4
  skipped).

## Ex ante gates (all encoded in acceptance_criteria.csv, including #7)

1. Strict, bit-identical two tight repeats.
2. Estate median within 1 bootstrap SE (|gap| <= 0.232).
3. Composition share within 2 bootstrap SEs of [PENDING-DATA].
4. Established-12 summed loss <= 4.5.
5. Young liquid wealth (25–35) within +/-0.15 of 0.179 — the watch moment:
   realistic income risk pushes it up, honest entry wealth pushes it down;
   the July-10 rejection of realistic risk was contaminated by the old entry
   condition.
6. Free-theta0 winner weakly below the theta0=0 nested seed.
7. Identification: Jacobian ranks, condition number, and weakest-direction
   composition reported AND recorded as a criterion row; a flat or
   bound-pinned theta1 is reported as set-identification (profile band from
   dependent job (a)), never as a point estimate.

Failure disposition: if gates 2/3/5 cannot hold jointly under the estimated
income process, that is evidence about a missing mechanism (medical/LTC
expense risk is the predeclared leading candidate, externally calibrated) —
not a license to shrink income risk again.

## Documentation duties

Run README written before launch; CALIBRATION_STATUS updated at launch and at
collection; all code committed before launch; this contract updated in place
when [PENDING-DATA] and [DECISION] slots resolve, with the resolution dated.
