# RECEIPT — per-child housing need probe (experimental mechanism test; NOT a calibration)

## Outcome
- Step 0 gate PASSED: fixed-price reproduction of the chain54 quarter point is bit-identical.
  - `first_birth_flow`: probe 0.049745516141897185 vs saved 0.049745516141897185 (diff 0.0).
  - `ownership_30_55`: probe 0.6457219408551765 vs saved 0.6457219408551765 (diff 0.0).
- All 4 solves done (2 need arms x 2 financed shares), fixed price AND fixed H0, no GE root, no recalibration.
- `results.csv` written (22 rows x 4 cases + financing effects within need + DiD). No economic interpretation per TASK.

### Child-floor semantics (engine/shared.py:263-290)
- `nk` = number of children currently at home (child-state axis `cs`), because
  `child_state_mode == "independent_count"` (verified on the built P object).
  Under the shared-clock mode it would instead be parity (ever born); that mode is off.
- Reference values: `hbar_first_child_jump = 2.49169824815624`, `hbar_child_rooms = 0.0`
  (== `h_P`; floor at one child at home = 2.49169824815624 rooms; childless floor 0).
- need1: `hbar_child_rooms = 1.0`, `hbar_first_child_jump = 1.49169824815624`
  (one-child floor unchanged at 2.49169824815624; two children +1 room; three +2).

### Case parameter objects (exact; all else = reference point bit-identical via deepcopy)
- Reference point (10 free coords): beta_annual 0.9678020936224091, chi 1.0912063362521283,
  child_benefit_curvature 0.09799189459865237, first_birth_fixed_cost 0.39721835190447047,
  h_P 2.49169824815624, kappa_fert 0.12761433940425662,
  kappa_fert_continuation 0.34800333637865355, psi_child 0.1742874677877052,
  tenure_choice_kappa 0.015482897400866402, theta0 0.10903797455048245.
- Fixed in all cases: price 0.6744838540900874, H0 7.288573389887633,
  purchase_saving_fraction 0.25 (quarter rule), unsecured credit d_bar 0.0,
  entry `nonnegative_mean`, grid 120x9, population scale 1.0.
- Per case: need0_phi08 (jump 2.49169824815624, rooms 0.0, phi 0.8);
  need0_phi10 (same floor, phi 1.0); need1_phi08 (jump 1.49169824815624, rooms 1.0, phi 0.8);
  need1_phi10 (same need1 floor, phi 1.0). `params_used.json` pins this machine-readably;
  each report dir holds its own 31-row `parameters.csv` (h_P row reports the first-child jump).

### Solve times and flags
- Lifecycle solve seconds: need0_phi08 7.60, need0_phi10 7.51, need1_phi08 7.17, need1_phi10 7.34.
- `fp.gates` (stationary=True): PASSED in all 4 cases; 17 standard plots each; no exceptions.
- Measured fixed-price diagnostics (not gates): renewal residuals
  (-6.95e-10, -0.00491, -0.0816, -0.0864); absolute housing residuals at fixed H0
  (0.0, 0.0626, 0.0371, 0.0948); PAYGO residuals ~1e-14 (all cases);
  adult-entry relative gaps recorded per case in `case_results.json`.
- Reporting accommodation (need1 only): the frozen `actual_parameters` adapter requires
  `hbar_child_rooms == 0`, so need1 `parameters.csv` rows were produced via an identical copy
  with only that field zeroed; all other estimates exact, h_P row = need1 first-child jump.
  Flagged as `hbar_reporting_accommodation: true` in `case_results.json`/closure.

### Row definitions used in results.csv
- Birth flows: sums over `parity_birth_flows_by_age` columns (uniform_birth_time observer);
  g-array cross-check (`g_post - g_pre` by parity) agrees to <1e-13 (separate rows).
- Completed fertility: `chain.extract_moments` TFR == target_fit `initial_normalization` model.
- Childlessness 40-44 = target_fit `cps_childlessness` model (observer `childless_rate_40_44`).
- `own_rate_25_34` = `production_whole_nodes` row of `young_ownership_age_measurement.csv`
  (whole-node `age_to_index` aggregation; two uniform-cell alternatives stored alongside).
- Renter cap shares: renters (tenure 0) with `hR_pol >= hR_max - 1e-9`, by children at home
  m = 1, 2, 3 (m = 0 memo row); shares are of renter parents with that m.
- Entrant cell: age cell j = 0 (ages 18-21); entrant mass identical in all cases by
  construction and verified all parity-0 (`entrant_nonzero_parity_mass` = 0.0 in all 4);
  first-birth prob = j0 parity 0->1+ flow / entrant mass; ownership = j0 owner / j0 mass.

## Verification
- Engine identity: executed engine `purchase_rules_overnight_v1/engines/quarter` SHA256 matches
  the TASK-named `code/model/experiments/quarter_saving_solvency/source` engine
  (shared.py `1fbe36f7…`, solver.py `19e6f924…`); all `engine_pins.json`, `source_pins.json`,
  and packet `manifest.json` SHAs asserted at startup.
- Frozen-source overlay identical to local chain54 runs (`local_runtime/bootstrap.py`
  replica: `frozen_sources/e5f_exact_policy_cache.py` digest `d51bbd13…` asserted).
- Reference arm equals the saved GE point bit-for-bit (see gate diffs above; renewal
  residual reproduces `-6.952464159937222e-10`; housing residual 0.0 at fixed H0).
- Selected-point identity: chain54 `search/completed.json` selected
  (loss 51.55603604909936, price 0.6744838540900874, H0 7.288573389887633).

## Artifacts (all inside this folder; nothing outside it edited)
- `run_probe.py` (driver), `build_results.py` (table builder), `results.csv`, `params_used.json`.
- `reference_point.json`, `floor_semantics.json`, `case_results.json`, `latest_completed.json`,
  `completed.json`, `setup/` (auth + utility verification).
- Per case `need*_phi*/phase_b_ge/need*_phi*_selected/`: `closure.json`, `target_fit.csv`,
  `parameters.csv`, `observers.json`, `extra_moments.json`, `gates.json`,
  `lifecycle_2023.csv`, `young_ownership_age_measurement.csv`, `market_quantity_units.json`,
  `policy_array_summary.json`, `standard_diagnostics/` (17 PNGs), `stage/` (solution arrays).

## Unresolved
- None blocking. Variant-arm renewal/housing residuals are measured fixed-price diagnostics
  by design (no root, no renormalization); they are reported, not resolved.
- `fp.gates` passed for all arms, so no convergence investigation was needed.

## Reported cost
- Wall time well within the 75-minute TASK budget: 4 lifecycle solves (~30 s total),
  observers + 4 x 17 plots (minutes); single-core local run with
  `NUMBA_NUM_THREADS=OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`,
  `code/model/.venv/bin/python`.
