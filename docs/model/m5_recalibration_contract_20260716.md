# M5 recalibration contract (drafted 2026-07-16, pre-launch)

Purpose: recalibrate the one-market intergenerational model conditional on
externally disciplined inputs, replacing the three provisional constructions
identified by the M4 audit (`docs/model/intergen_m4_calibration_audit_20260716.md`):
the reverse-engineered income process, the circular entrant-wealth condition,
and the grid-discrete nonhousing-median target. M4
(`output/model/intergen_standard_bequest_recalibration_20260716/final_report/`)
is the warm start and comparison baseline, not the paper object.

DECISIONS RESOLVED 2026-07-16 night (Tommaso): D1 = literature AR(1),
Sommer–Sullivan `rho=0.90`, `sigma=0.20` annual ("unfortunately in the past we
found out it was very hard to keep sigma at 0.2 but let's try"); D4 = tenure
kappa FREED with the ACS old-age ownership moment ("for ownership, acs is
fine"). Contract is therefore **15 moments / 14 free parameters**.

## External inputs (new)

1. **Income process** — literature AR(1): annual `rho = 0.90`,
   `sigma = 0.20` (Sommer–Sullivan anchor), fed through the verified
   conversion (`local_panel.py:1049-1094`; rho^4 = 0.6561, stationary log
   variance 0.2105, 5-state Rouwenhorst z ≈ [0.36, 2.26]). Own PSID estimates
   (Block 1: persistent+transitory rho=0.9749, sigma_eta=0.1769, sigma_e=0.4128;
   AR(1)-only rho=0.9436, sigma=0.2846) are appendix evidence that the
   literature value is conservative — family income without education/fixed-
   effect controls is even more dispersed. The matched constants (0.9602,
   0.0645) are retired everywhere. Known risk, accepted with eyes open: the
   July-10 sigma=0.20 test blew up the young liquid-wealth moment (1.90 vs
   0.179), but under the circular entry condition and fixed parameters; gate 5
   adjudicates it honestly this time.
2. **Entrant wealth at age 18** — 18–24 childless-renter quintile ratio
   distribution (Block 2: 1,835 family-years; bin means −2.2225/−0.0526/
   0.1042/0.3516/3.1035, mean 0.2589, median 0.0985), replacing the circular
   25–35 injection. The 25–35 mean (0.17922556) remains a hard target — now a
   genuine lifecycle test.

## Target system (v3, resolved)

15 moments / 14 free parameters:

- The 12 established moments, unchanged values and weights.
- `old_total_estate_wealth_to_annual_income_median_7684` = 6.50131577436537,
  weight 18.585767349158665 (unchanged; identifies theta0 — nested-zero
  evidence: dropping theta0 costs 303 loss points through this row alone).
- **NEW** `old_nonhousing_ge_1x_income_share_6575` = 0.608333139649131,
  weight 9435.18732291246 (1/SE^2, SE 0.0102949617302331, Block 3) — replaces
  the grid-discrete nonhousing median (Jacobian row exactly zero; node
  spacing ~3 SEs). The median stays as a diagnostic.
- **RESTORED** `old_age_own_rate` = 0.76426097 (ACS-sourced, legacy weight
  160.0 — kept ACS for consistency with the other ownership targets; the
  PSID reference-person alternative 0.834, SE 0.0077, is documented in Block
  4). Identifies the freed `tenure_choice_kappa`: its variation is the
  ownership age path, and the moment was demoted at M1 only conditional on
  kappa = 0.

Untargeted diagnostics reported every readout: estate p90/p50, family-size
estate gap, housing/nonhousing decomposition 76–84, ownership path 62–78,
young liquid wealth by age, nonhousing median 65–75.

## Parameters

- Free (14): the 11 clean-frontier coordinates (M4 bounds/transforms
  unchanged) + theta0 [0, 8] softzero + theta1 [0.02, 16] log +
  tenure_choice_kappa [0, 0.12] softzero (the production-profile search range;
  starts include 0 since every prior estimate sat there).
- External: theta_n = 0 (Hurd 1989; Kopczuk–Lupton 2007; DFJM 2025; own PSID
  gap 0.101, SE 0.563);
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
- Nested seed: exact (theta0 = 0, tenure_choice_kappa = 0) at the strict M1
  coordinates — the no-bequest, deterministic-tenure submodel; the collector
  fails if the free winner ends above it.
- Dependent jobs behind the strict collector:
  (a) two-sided theta1 profile, predeclared ladder spanning [0.1, 2.0], other
      coordinates re-optimized per cell — the interior-optimum question is
      testable this time;
  (b) fresh Jacobian at the winner (kappa sweep obsolete — kappa is free).
- Winner verified at Nb=240 before any promotion (standing convention that M4
  skipped).

## Ex ante gates (all encoded in acceptance_criteria.csv, including #8)

1. Strict, bit-identical two tight repeats.
2. Estate median within 1 bootstrap SE (|gap| <= 0.232).
3. Composition share within 2 bootstrap SEs (|gap| <= 0.0206).
3b. Old-age ownership |gap| <= 0.03 vs the ACS 0.764.
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
