# RECEIPT — selling-cost dose-response + lock-in check (experimental; NOT a calibration; nothing adopted)

## Outcome
- Part 1: all 8 dose-response solves done (P.psi in {0.00, 0.02, 0.04, 0.06}
  x financed share {0.8, 1.0}), all need1
  (hbar_child_rooms 1.0, hbar_first_child_jump 1.49169824815624), fixed price
  0.6744838540900874 AND fixed H0 7.288573389887633, no GE root, no
  recalibration. Labels psi06/psi04/psi02/psi00 x phi08/phi10.
- `results.csv` written (18 rows x 8 cases + 4 within-psi financing effects
  phi10-phi08). Rows: first/second/third birth flows; completed fertility
  (TFR); childlessness 40-44; ownership_30_55; own_rate_25_34;
  first_birth_rooms; renter-parent cap share m=1,2,3; matched entrant
  (age 18-21) first-birth prob and ownership; plus diagnostics
  (gates_status, renewal/housing/PAYGO residuals, lifecycle solve seconds).
  No economic interpretation per TASK.
- Identity gate PASSED bit-for-bit: psi06 arms reproduce B1 first_birth_flow
  (diff 0.0 both phis); psi00 arms reproduce B2 first_birth_flow (diff 0.0
  both phis). Threshold asked: exact.
- Part 2: lock-in tables written (`part2_tables/`) from the saved psi06/psi00
  arrays (both phis), zero extra solves. Childless (n=0, m/cs=0), age cells
  0-4, by beginning tenure state: at-risk mass + share, attempt prob, pi_j,
  birth prob. Post-birth tenure transitions (keep/upsize/downsize/to-rent):
  NOT available from saved arrays (verdict recorded in part2_summary.md).
- No repo source file edited. Torch never attempted per TASK.

### Case parameter objects (exact; all else = reference point bit-identical via deepcopy)
- Reference point (10 free coords): beta_annual 0.9678020936224091,
  chi 1.0912063362521283, child_benefit_curvature 0.09799189459865237,
  first_birth_fixed_cost 0.39721835190447047, h_P 2.49169824815624,
  kappa_fert 0.12761433940425662, kappa_fert_continuation 0.34800333637865355,
  psi_child 0.1742874677877052, tenure_choice_kappa 0.015482897400866402,
  theta0 0.10903797455048245.
- Fixed in all cases: price 0.6744838540900874, H0 7.288573389887633,
  purchase_saving_fraction 0.25 (quarter rule), unsecured credit d_bar 0.0,
  entry `nonnegative_mean`, grid 120x9, H_own (2,4,6,8,10) n_house 5,
  population scale 1.0. Only P.psi varies (None = reference 0.06; else
  0.04/0.02/0.0 with `selling_cost` parameter row re-flagged) and the
  financed-share vector (0.8/1.0). `params_used.json` pins this
  machine-readably; each report dir holds its own 31-row `parameters.csv`
  (need1 hbar reporting accommodation as in the previous probes).

### Solve times and flags
- Lifecycle solve seconds: 7.2-7.6 per case (all 5-rung).
- `fp.gates` (stationary=True): PASSED in all 8 cases; 17 standard plots each.
- Measured fixed-price diagnostics (not gates; no root by design): renewal
  residuals -0.086 to -0.061; absolute housing residuals +0.025 to +0.218;
  PAYGO residuals <= 3.2e-14 (all cases). Full per-case values in
  `results.csv`; closures/gates in each report dir.

### Part 2 definitions (exact; see part2_tables/part2_summary.md)
- At-risk mass g_pre reconstructed exactly as the observer does: entrant
  cohort at j=0 from saved entry_by_loc; advance of saved post-fertility
  g_beginning_distribution otherwise (saved loc/tenure/saving policy, rebuilt
  location/tenure maps). Saved g_beginning_distribution is POST-fertility
  (fertility is applied to it in place in the KFE loop); weighting attempt
  probs by it would understate at-risk means by birth selection, so g_pre is
  used. Fecundity pi_j (j=0..6): 0.98, 0.965817, 0.941576, 0.900144,
  0.82933, 0.708298, 0.501436; 0 after. Engine flags recorded in the summary
  (sequential births; stochastic aging; survival; due-stayer; entry-censor
  with 0.0 censored mass; no readiness gate; no joint nesting).

## Verification
- Engine/source identity: `run_dose_probe.py` is `run_barrier_probe.py` with
  only docstring/budget/CASES/psi-status-line changed; all
  `setup_imports`/`build_base`/solve/observe logic identical, so all
  `source_pins.json`, `engine_pins.json`, packet `manifest.json` SHAs and the
  frozen-source overlay were asserted at startup before every solve (same
  executed quarter engine). Smoke case psi04_phi08 passed standalone first;
  the full run then skipped it via case_results.json and completed the rest.
- Identity: psi06 vs B1 and psi00 vs B2 first_birth_flow diffs all 0.0.
- Part 2 exactness (part2_tables/part2_verification.csv, all 4 cases):
  re-applying the KFE sequential fertility logic to reconstructed g_pre
  reproduces saved post-fertility arrays to L1 ~1.4e-16; implied aggregate
  first births equal observer first_birth_flow to 0.0 (one case 6.9e-18);
  reconstructed entry mass equals saved entry rate; censored mass 0.0.
  Cross-check: j=0 renter birth prob 0.22375257... equals the observer
  matched-entrant first-birth prob bit-for-bit.
- Post-birth tenure verdict verified by stage key inventory (87 keys; only
  tenure_choice/tenure_probs, both pre-tenure policies by beginning state;
  identical across the 4 cases).

## Artifacts (all inside this folder; no repo source file edited)
- `run_dose_probe.py` (driver), `build_dose_results.py` (table builder),
  `part2_lockin.py` (Part 2 reconstruction + tables), `results.csv`,
  `params_used.json`, `reference_point.json`, `floor_semantics.json`.
- `part2_tables/` (mass CSV, first-birth CSV, verification CSV, summary MD).
- `case_results.json`, `latest_completed.json`, `completed.json`,
  `probe_full.log`, `setup/`, `part2_setup_tmp/` (auth + utility verification).
- Per case `<label>/phase_b_ge/<label>/`: `closure.json`, `target_fit.csv`,
  `parameters.csv`, `observers.json`, `extra_moments.json`, `gates.json`,
  `lifecycle_2023.csv`, `young_ownership_age_measurement.csv`,
  `market_quantity_units.json`, `policy_array_summary.json`,
  `standard_diagnostics/` (17 PNGs), `stage/` (solution arrays).
- Job IDs: none (local only; Torch never attempted per TASK). Local run:
  background pid 82213 (plus the earlier psi04_phi08 smoke process),
  single thread via `code/model/.venv/bin/python` with
  `NUMBA/OMP/OPENBLAS/MKL_NUM_THREADS=1`.

## Unresolved
- A first draft of Part 2 weighted attempt probs by saved
  g_beginning_distribution (post-fertility) mass; that understates at-risk
  means (e.g. j=0 birth prob 0.1409 vs the engine-true 0.2237). Caught by the
  entrant-cell cross-check before delivery; the shipped tables use the
  exactly-reconstructed at-risk mass (L1 ~1e-16). No stale table remains.
- Choice-specific tenure values and post-birth-branch tenure policies are not
  saved (engine returns VH/tcj/probs only); keep/upsize/downsize/to-rent
  shares after a first birth cannot be computed -- stated in tables, not
  resolved.
- Variant-arm renewal/housing residuals are measured fixed-price diagnostics
  by design (no root, no renormalization); reported, not resolved. `fp.gates`
  passed for all arms.

## Reported cost
- Wall time ~15 minutes total, within the 75-minute TASK budget: inspection +
  driver authoring; 8 lifecycle solves (~1 min total); observers + 8 x 17
  plots (remainder); Part 2 reconstruction (minutes). Single-core local run
  as above. Torch: 0 s (never attempted per TASK).
