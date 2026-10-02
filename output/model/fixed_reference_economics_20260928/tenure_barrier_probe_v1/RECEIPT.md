# RECEIPT — tenure-barrier probe (experimental mechanism test; NOT a calibration; nothing adopted)

## Outcome
- Part A tables written (`partA_tables/`: per-(m, age-cell) CSVs for need1_phi10 and
  need1_phi08, pooled CSV, summary MD). Zero model solves.
- Part B: all 12 counterfactual solves done (6 arms x phi 0.8/1.0), fixed price
  0.6744838540900874 AND fixed H0 7.288573389887633, no GE root, no recalibration.
- `results.csv` written (23 rows x 12 cases + 6 within-arm financing effects + 5
  arm-minus-B1 DiDs). No economic interpretation per TASK.
- B1 identity gate PASSED: B1 first_birth_flow reproduces previous probe need1 arms
  bit-for-bit — need1_phi08 diff 0.0, need1_phi10 diff 0.0 (threshold asked ~1e-10).
- B3/B4/B6 proceeded: the engine supports an arbitrary-length owner grid (see
  Verification). Nothing was skipped.

### Case parameter objects (exact; all else = reference point bit-identical via deepcopy)
- Reference point (10 free coords): beta_annual 0.9678020936224091, chi 1.0912063362521283,
  child_benefit_curvature 0.09799189459865237, first_birth_fixed_cost 0.39721835190447047,
  h_P 2.49169824815624, kappa_fert 0.12761433940425662,
  kappa_fert_continuation 0.34800333637865355, psi_child 0.1742874677877052,
  tenure_choice_kappa 0.015482897400866402, theta0 0.10903797455048245.
- Fixed in all cases: price 0.6744838540900874, H0 7.288573389887633,
  purchase_saving_fraction 0.25 (quarter rule), unsecured credit d_bar 0.0,
  entry `nonnegative_mean`, grid 120x9, population scale 1.0.
- Per arm: B1 need1 (jump 1.49169824815624, rooms 1.0); B2 = B1 + P.psi 0.06 -> 0.0
  (selling cost; expected `selling_cost` row updated, status flagged);
  B3 = B1 + H_own (2,...,10), n_house 5 -> 9; B4 = B3 + P.psi 0.0;
  B5 = B1 + tenure_choice_kappa -> 0.001 (plan.json search-range lower bound 0.001
  verified); B6 need0 (jump 2.49169824815624, rooms 0.0) + finer grid.
  `params_used.json` pins this machine-readably; each report dir holds its own
  31-row `parameters.csv` (need1 arms use the same hbar reporting accommodation as
  the previous probe).

### Part A headline numbers (pooled; full cells in partA_tables/)
- Screen applied is the operative income-augmented one: with
  `native_purchase_income=True`, renter buying rung tn is feasible iff
  b >= (1-phi)*p*H - y/Rg AND b - p*H >= max(-phi*p*H - y/Rg, b_grid[0]).
- Feasibility is ~100% on every rung in both arms (phi 0.8 and 1.0): mean beginning
  wealth ~1.5-1.6 with mean current income ~3.6-4.0 covers every down payment.
- Tenure choice still puts 36-60% of this population on rent (t0); rung H=2 gets
  ~0 probability in all cells; remainder concentrates on rungs 8/10 (m>=2) and 6/8/10.
- Mean income state z-index 4.0-4.7 of 0-8; zero mass at either wealth-grid bound
  (grid -12 to 3000). Choice-specific tenure values are not saved, so no
  rent-vs-best-rung value gap is reported (see Unresolved).

### Solve times and flags
- Lifecycle solve seconds: 7.0-7.8 (5-rung arms), 13.3-13.8 (9-rung arms).
- `fp.gates` (stationary=True): PASSED in all 12 cases; 17 standard plots each.
- Measured fixed-price diagnostics (not gates; no root by design): renewal residuals
  -0.092 to -0.004 (B6 need0-fine ~-0.004/-0.010); absolute housing residuals
  -0.149 to +0.638 (B5 low-kappa largest); PAYGO residuals <= 2.2e-12 (all cases).
  Full per-case values in `results.csv`; closures/gates in each report dir.

### Row definitions used in results.csv
- Birth flows: `parity_birth_flows_by_age` columns (uniform_birth_time observer).
- Completed fertility: `chain.extract_moments` TFR == target_fit `initial_normalization` model.
- Childlessness 40-44 = target_fit `cps_childlessness` model.
- `own_rate_25_34` = `production_whole_nodes` row of `young_ownership_age_measurement.csv`.
- Renter cap shares: realized renters (tenure 0) with `hR_pol >= hR_max - 1e-9`, by
  children at home m = 1, 2, 3; shares are of renter parents with that m.
- Ownership of parents by m: owner share of `g_current` mass at child-state m
  (m = 0 memo row included).
- Entrant cell: age cell j = 0 (ages 18-21); first-birth prob and ownership as in
  the previous probe.

## Verification
- Engine/source identity: `run_barrier_probe.py` reuses the previous probe's
  `setup_imports`/`build_base` unchanged, so all `source_pins.json`,
  `engine_pins.json`, packet `manifest.json` SHAs and the frozen-source overlay
  were asserted at startup before every solve (same executed quarter engine).
- Reference-arm identity: B1 first_birth_flow diff 0.0 vs saved need1 arms (both phis).
- Owner-grid support (bounded static search + empirical smoke): core loops derive the
  tenure count from shapes (`nt = 1 + P.n_house`): shared.py:330/333-334/504/529/533/550/588,
  household.py:154/179-180, kernels.py tenure kernels take `nt` from `Vd.shape`,
  joint_nested.py:223; no hard-coded rung count found in engine or native
  pilot/runner. B3 smoke then passed solve + 14-row target fit + gates + plots, so
  B3/B4/B6 ran as specified.
- Part A screen replication is exact: with the income-augmented screen, positive
  choice probability on infeasible (b, rung) cells is 0.0 in all 96 occupied
  (m, age) cells. (A beginning-wealth-only screen does NOT replicate: it flags up
  to 0.0012 mass with large buy probabilities at b ~ 0-0.3 in high income states,
  which the engine feasibly serves out of current income.)
- Part A weights: saved `g_beginning_distribution x tenure_probs[...,0]` summed over
  origins reproduces realized renter mass `g[:,0]` to ~3e-10 in both arms, so
  `g_beginning_distribution` is the saved pre-tenure (pre-choice) distribution used
  for weights. `g_post_fertility` is not among the 87 saved stage arrays.
- Tenure-choice-specific values (rent vs best rung) are not saved: the kernels
  compute per-rung values in a local array and return only VH/tcj/probs
  (kernels.py:412-434); the stage npz key inventory (87 keys) contains no per-rung
  value array. Value gap skipped as instructed.

## Artifacts (all inside this folder; no repo source file edited)
- `run_barrier_probe.py` (driver, adapted from previous `run_probe.py`),
  `build_barrier_results.py` (table builder), `partA.py` (Part A decomposition),
  `results.csv`, `params_used.json`, `reference_point.json`, `floor_semantics.json`.
- `partA_tables/` (2 x 51-row cell CSVs, pooled CSV, summary MD).
- `case_results.json`, `latest_completed.json`, `completed.json`, `probe_full.log`,
  `partA_setup_tmp/` (auth + utility verification).
- Per case `<label>/phase_b_ge/<label>/`: `closure.json`, `target_fit.csv`,
  `parameters.csv`, `observers.json`, `extra_moments.json`, `gates.json`,
  `lifecycle_2023.csv`, `young_ownership_age_measurement.csv`,
  `market_quantity_units.json`, `policy_array_summary.json`,
  `standard_diagnostics/` (17 PNGs), `stage/` (solution arrays).
- Small check scripts: `check_arrays.py`, `check_screen.py`, `check_flags.py`,
  `check_infeas.py`.
- Job IDs: no new Torch job was submitted (see Unresolved). Pre-existing production
  array 19040483 (24 purchase-rule calibration tasks) was verified RUNNING at task
  start and was left untouched. Local Part B run: background pid 60572 (plus B3
  smoke pid before it), single thread via `code/model/.venv/bin/python`.

## Unresolved
- Torch side not used: `torch.sh status` passed at ~12:2x EDT, then fresh SSH
  connections failed (`Permission denied`, exit 255) for both direct `ssh torch`
  and `torch.sh status` at ~12:40 EDT. Per TASK fallback, all 12 solves ran locally
  (~7 s each at 5 rungs, ~14 s at 9 rungs); the NEW dated folder under the Torch
  scratch tree was NOT created (no transport). Nothing was submitted to the cluster.
- `g_post_fertility` is not in the saved stage arrays; Part A weights use saved
  `g_beginning_distribution` (pre-choice mass), with the ~3e-10 reproduction check
  above as the equivalence evidence.
- Choice-specific tenure values are not saved; the rent-vs-best-rung value gap and
  its distribution are skipped (instructed fallback).
- Variant-arm renewal/housing residuals are measured fixed-price diagnostics by
  design (no root, no renormalization); reported, not resolved. `fp.gates` passed
  for all arms, so no convergence investigation was needed.

## Reported cost
- Wall time ~35 minutes total (12:12-12:46 EDT), within the 2.5-hour TASK budget:
  Part A (imports + build_base + tables) minutes; 12 lifecycle solves (~2 min total);
  observers + 12 x 17 plots (remainder). Single-core local run with
  `NUMBA_NUM_THREADS=OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`,
  `code/model/.venv/bin/python`. Torch jobs <= 1 hour: none submitted (0 s).
