# Intergen Quantitative Code SOTA Index

Date: 2026-06-22

Status: codebase-orientation and speedup-planning note. This does not report
calibration progress and does not change model logic.

## Source-Of-Truth Read

Required startup files read for this index:

- `AGENTS.md` and `CLAUDE.md`
- `memory/AGENT_MEMORY.md`
- latest daily note, `memory/daily/2026-06-22.md`
- `CALIBRATION_STATUS.md`
- `code/model/intergen_housing_fertility/README.md`
- `code/model/intergen_housing_fertility/IMPLEMENTATION_STATUS.md`
- `docs/model/intergen_three_day_calibration_log_20260619.md`
- `git status -sb`

There is a live orientation discrepancy. `CALIBRATION_STATUS.md`, the intergen
package README, and the June 19-22 log identify the June 2026 active
quantitative strand as the one-market/no-location intergenerational model under
`code/model/intergen_housing_fertility/`. Some older memory text and
`code/model/README.md` still describe the center-periphery `dt_cp_model` path as
"active." For this strand, treat `dt_cp_model` as inherited/reference code, not
the current calibration object.

The worktree was already dirty before this note, mostly in older
center-periphery code, LaTeX outputs, and untracked empirical/prompt files. This
index does not touch those files.

## Active Model Package

Active package:

- `code/model/intergen_housing_fertility/`

Core files:

- `README.md`: package entry note and standard commands.
- `IMPLEMENTATION_STATUS.md`: live implementation ledger and shortcut list.
- `parameters.py`: one-market defaults, 4-year periods, Markov income process,
  owner ladder, renter cap, financing/PTI constraints, and supply closure.
- `solver.py`: GE, fixed-price, Markov-income Bellman, forward distribution,
  statistics, and legacy non-Markov reference path.
- `calibration.py`: diagnostic target sets, weights, moment extraction, and
  small diagnostic calibration driver.
- `local_panel.py`: bounded local-panel and global-DE diagnostic search stack,
  including the 13-parameter internal vector and 5-state income override.
- `diagnostics.py`: diagnostic packet writer.
- `cli.py`: public command line surface.
- `kernels.py`: Numba kernels for savings, policy evaluation, tenure choice,
  interpolation, and forward distribution.

Usual local runtime:

```bash
cd code/model
.venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
```

The active quantitative object is one market, no location, persistent Markov
income, one dependent-child stage, one-shot completed fertility, owner rungs,
continuous renter housing, down-payment constraints, PTI screening, and a
stationary aggregate housing-services market. It is not a production
calibration.

## Active Solver, Calibration, And Diagnostic Files

Solver and model inspection:

- `code/model/intergen_housing_fertility/solver.py`
- `code/model/intergen_housing_fertility/parameters.py`
- `code/model/intergen_housing_fertility/diagnostics.py`
- `code/model/tools/audit_intergen_final_best_pathologies.py`
- `code/model/tools/audit_intergen_parent_credit_margin.py`
- `code/model/tools/run_intergen_policy_poc.py`

Diagnostic calibration/search:

- `code/model/intergen_housing_fertility/calibration.py`
- `code/model/intergen_housing_fertility/local_panel.py`
- `code/model/tools/collect_intergen_panel_results.py`
- `code/model/tools/audit_intergen_sensitivity_jacobian.py`
- `code/model/tools/audit_intergen_mechanism_grid.py`

Active cluster launchers:

- `code/cluster/submit_intergen_housing_fertility_twohour_panel.sh`
- `code/cluster/submit_intergen_housing_fertility_global_de.sh`
- `code/cluster/submit_intergen_sensitivity_jacobian.sh`
- `code/cluster/submit_intergen_mechanism_grid.sh`
- `code/cluster/submit_intergen_mechanism_repair_grid.sh`
- `code/cluster/submit_intergen_housing_fertility_calibration.sh`
- `code/cluster/collect_intergen_shutdown_snapshot.sh`
- `code/cluster/pull_intergen_shutdown_snapshots_local.sh`

The mechanism-grid launchers are diagnostic tools, not new SMM target systems.
They preserve frontier parameter points and run deterministic perturbation
blocks.

## Active Target Sets

The active target inventory is in
`code/model/intergen_housing_fertility/calibration.py`. None of these should be
described as a final production SMM objective.

Important current sets:

- `candidate_no_timing_v0`: old non-location moments excluding first-birth
  timing, plus midlife wealth, housing user-cost share, and childless
  renter/owner room medians.
- `candidate_no_timing_core_feasibility_v1`,
  `candidate_no_timing_cost_test_v1`,
  `candidate_no_timing_oldage_test_v1`, and
  `candidate_no_timing_roomcost_test_v1`: June 17 frontier diagnostic subsets.
- `candidate_replacement_v1`: ACS/PSID replacement target trial with room
  means/shares, old nonhousing wealth, and childless renter/owner room objects.
- `candidate_replacement_nh_median_v1`: replaces old nonhousing mean wealth
  with old nonhousing median wealth.
- `candidate_replacement_due_lifecycle_v1` and
  `candidate_replacement_due_lifecycle_soft_v1`: add lifecycle ownership-slope
  pressure.
- `candidate_replacement_due_lifecycle_owngap_v1`: uses old
  parent-childless ownership gap instead of old parent-childless nonhousing
  wealth gap.
- `candidate_replacement_old_retention_v1`: direct old ownership and old
  ownership-gap pressure.
- `candidate_replacement_young_old_own_v1`: adds young ownership.
- `candidate_replacement_young_old_roomgap_v1`: adds owner-renter room-gap
  pressure.
- `candidate_replacement_total_median_v1` and
  `candidate_replacement_total_due_lifecycle_v1`: total-wealth robustness
  variants.

Identification rule: the search vector has 13 internal parameters,
\(\{\beta,\alpha,b_0,\bar c_0,\bar c_n,\bar h_0,\bar h_{\mathrm{jump}},
\bar h_n,\psi_{\mathrm{child}},\kappa_n,\chi,\theta_0,\theta_n\}\). Any target
revision must keep at least 13 informative moments or fix parameters
externally. Dropping a difficult target is not valid unless a replacement
moment identifies the same parameter block.

## Active Output And Result Locations

Local output/result locations:

- `output/model/intergen_globalde_final_best_diagnostics/`: June 9 final
  global-DE toy best diagnostic packet. It is a reference diagnostic point, not
  a calibrated benchmark.
- `output/model/cluster_pulls/intergen_overnight_frontier_20260617/`: June 17
  frontier target-set ensemble pulls.
- `output/model/intergen_sensitivity_jacobian_20260618/`: June 18 sensitivity
  and Jacobian audit outputs.
- `output/model/cluster_pulls/intergen_replacement_cluster_wave_20260618/`:
  June 18 replacement-target cluster pull.
- `output/model/intergen_repair_best_diagnostics_20260622/`: diagnostic packet
  for the best joint-distance focused repair-grid case.
- `output/model/intergen_repair_best_loss_with_roomgap_diagnostics_20260622/`:
  diagnostic packet for a focused repair case with better scalar loss and room
  gap.
- `output/model/intergen_speed_audit_20260622/`: local timing audit from this
  note.

Torch scratch locations named by the June 19-22 log:

- `/scratch/td2248/projects/Fertility_Spring26_20260617_fast/code/cluster`
- `/scratch/td2248/projects/Fertility_Spring26_20260617_fast/code/cluster/audit_points_20260621_frontier/`
- `/scratch/td2248/projects/Fertility_Spring26_20260617_fast/output/model/intergen_frontier_jacobian_20260621_*`
- `/scratch/td2248/projects/Fertility_Spring26_20260617_fast/output/model/intergen_mechanism_grid_20260622/`
- `/scratch/td2248/projects/Fertility_Spring26_20260617_fast/output/model/intergen_mechanism_repair_grid_20260622/aggregate_repair_readout.{txt,json}`

## Stale, Reference, And Archival Paths

Reference but not current intergen calibration object:

- `code/model/dt_cp_model/`: inherited center-periphery model. It is useful for
  solver patterns, especially the existing Howard policy-evaluation cycle, but
  not the June 2026 one-market intergen calibration object.
- `code/model/README.md`: currently overstates the center-periphery path as the
  active Python implementation. It should eventually be clarified.

Historical/archival:

- `calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/`
- `calibration_archive/legacy_matlab_2026-05-07/`
- older intergen runs from 2026-06-05 and 2026-06-06 before the corrected TFR
  extraction and Markov event-statistic path.

Do not compare losses across target systems, room-unit normalizations,
geographies, or objective definitions without checking comparability.

## Current Best Diagnostic Fits

None is final. The useful current points are mechanism diagnostics:

- June 9 global-DE toy best:
  `output/model/intergen_globalde_final_best_diagnostics/source_record.json`
  has stored rank loss `11.503191936648555` under `candidate_no_timing_v0`.
  It is not economically acceptable: prime-age ownership is `0.222`, young
  ownership is `0.004`, old ownership is `0.783`, TFR is `1.660`, and the owner
  room fit is contaminated by the earlier owner-ladder choice.
- June 21 Wave 1-3 campaign scalar best: `old_retention_w3` loss `11.578`,
  TFR `1.736`, childlessness `0.287`, young ownership `0.006`, old ownership
  `0.908`, and room gap `1.389`. It gets a good scalar score by emptying the
  young-owner pipeline.
- June 22 mechanism grid best scalar source-target loss: scalar old-retention
  point with `phi=0.90`, loss `9.386`, TFR `1.743`, childlessness `0.285`,
  young ownership `0.006`, old ownership `0.862`, and room gap `1.363`. It is
  not an economic fix because young ownership remains essentially zero.
- June 22 focused repair best joint-distance case:
  `hR4.2_chi0.82_b3.5_phi0.90_fert`, rerun in
  `output/model/intergen_repair_best_diagnostics_20260622/`, has TFR `1.675`,
  childlessness `0.194`, young ownership `0.232`, old ownership `0.844`,
  renter rooms `4.192`, owner rooms `6.247`, room gap `2.055`, old NH median
  wealth `3.203`, and old parent-childless NH wealth gap `-0.582`. This is a
  real allocation improvement, but it misses young ownership, overstates old
  ownership relative to target, and has the wrong-signed old-wealth gap.
- June 22 focused repair lower-loss/room-gap case:
  `hR4.0_chi0.90_b1.5_phi0.90_fert`, rerun in
  `output/model/intergen_repair_best_loss_with_roomgap_diagnostics_20260622/`,
  has rerun rank loss `22.543`, TFR `1.962`, childlessness `0.177`, young
  ownership `0.248`, old ownership `0.872`, and room gap `2.214`. It still
  misses old ownership and has a negative old parent-childless NH wealth gap.

Conclusion: the scalar-best and mechanism-screen bests are diagnostic evidence,
not production calibration results.

## Known Economic Failures

- Young ownership, old exit, owner-renter room separation, fertility, and old
  wealth do not currently fit jointly.
- Lowering old ownership and generating room separation often empties the
  young-owner pipeline.
- Preserving young ownership and old exit tends to compress the owner-renter
  room gap because renters become too large.
- Positive old parent-childless nonhousing wealth gaps exist, but the cases
  with positive gaps tend to have too much childlessness and still-high old
  ownership.
- Old nonhousing median wealth loads mostly on `beta`; it is not a clean
  bequest-block identifier for `theta0`/`theta_n`.
- The current one-shot completed-fertility architecture cannot be interpreted
  as a sequential parity hazard.
- Target count and local rank are not enough. Several target systems are full
  rank locally, but the available levers create a sharp economic tradeoff.

## Known Numerical And Speed Bottlenecks

- The Markov-income GE path dispatches through
  `solve_markov_income_equilibrium` when `uses_markov_income(P)` is true
  (`solver.py` around lines 670-679).
- `solve_markov_income_equilibrium` calls `solve_markov_income_at_prices` in
  every damped price iteration, in scalar refinement, and once more for final
  full reporting (`solver.py` around lines 750-792).
- `solve_markov_income_at_prices` always calls
  `solve_bellman_full_markov_income` (`solver.py` around lines 964-979).
- Markov fixed-price timings explicitly report `bellman_mode="full_only"`,
  `n_full=1`, and `n_eval=0` (`solver.py` around lines 991-999).
- The non-Markov inherited path already has a Howard cycle using `stored_bp`
  and `solve_bellman_eval` (`solver.py` around lines 1287-1330 and
  1699-1706).
- Markov income adds a \(z\) dimension to policy/value arrays and requires
  continuation over \(\Pi_z\): `solve_bellman_full_markov_income` allocates
  arrays of shape `(Nb, nt, I, J, Nz, npar, ncs)` and computes
  \(E[V_{j+1}(z')|z]\) (`solver.py` around lines 1725-1756 and 1810-1819).
- Scalar refinement is expensive because every bracketing/root evaluation is
  currently a full Bellman solve.
- Distribution/statistics time is also material, especially for the final
  full-statistics pass.

## Speed Audit: 2026-06-22

Command used:

```bash
cd code/model
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python - <<'PY'  # inline wrapper around existing solver
```

The inline wrapper did not alter source code. It called the existing
`run_model_cp_dt` with:

- `base_overrides(J=16, Nb=60, n_house=5, max_iter_eq=3)`
- `income_process_overrides(5)`
- default baseline parameters, not a search candidate
- scalar refinement left on, matching the current calibration stack default

Output file:

- `output/model/intergen_speed_audit_20260622/markov_ge_j16_nb60_nz5_h5_maxiter3_scalar_refine_baseline.json`

Results:

| Metric | Value |
|---|---:|
| Grid | `J=16`, `Nb=60`, `Nz=5`, `n_house=5`, `H_own=[2,4,6,8,10]` |
| `max_iter_eq` | `3` |
| Wall time | `40.495` sec |
| Aggregate Bellman full time | `21.189` sec |
| Aggregate distribution/statistics time | `19.175` sec |
| Fixed-price elapsed time | `40.442` sec |
| Fixed-price calls | `9` |
| Full Bellman solves | `9` |
| Eval Bellman solves | `0` |
| Distribution passes | `9` |
| Fast-stat calls | `8` |
| Full-stat calls | `1` |
| Scalar refinement | Used, Brent, bracket found, `3` iterations |
| Final price | `0.8184156354` |
| Final market residual | `4.28e-06` |
| Final timing block | `bellman_mode="full_only"`, `n_full=1`, `n_eval=0` |

The final timing block alone understates total cost because it covers only the
accepted final fixed-price solve. Aggregating all fixed-price calls shows the
actual short-run cost of the GE loop plus scalar refinement.

## Howard Speedup Audit

The speedup is the right first in-repo path because the non-Markov model
already implements exactly the relevant pattern: solve the full Bellman
occasionally, store the savings policy `bp_pol`, and use a cheaper Bellman
policy-evaluation branch between full optimizations. This avoids repeatedly
running golden-section savings optimization when the GE price update is small.
It also preserves the existing exact dynamic-programming structure and can be
guarded with a full-only parity flag.

It likely dominates a neural-net solver at this stage because the current
bottleneck is repeated exact Bellman optimization inside a known state space,
not an unknown high-dimensional function class. A neural approximation would
introduce approximation error, training instability, extra constraint handling
for borrowing/PTI/tenure feasibility, and a new validation burden before the
calibration target system is economically settled. Howard evaluation is
deterministic, already proven in the inherited path, and directly testable
against full solves.

Important nuance: current diagnostic searches often use `max_iter_eq=3`. The
inherited non-Markov Howard loop performs the first three GE iterations as full
solves, so simply copying that loop and leaving scalar refinement full-only
would not speed those exact short runs. The Markov implementation should either
wire policy evaluation into scalar refinement or expose a conservative
Markov-specific warmup flag while always recomputing the accepted final solution
with a full Bellman solve.

## Markov Howard Implementation Plan

Add or extend functions:

- Add `solve_bellman_eval_markov_income(stored_bp, r_hat, p_hat, P, b_grid, SD)`.
- Preferably refactor the Markov Bellman into a shared
  `solve_bellman_markov_income_core(..., stored_bp=None, eval_mode=False)` so
  the full and eval branches stay parallel.
- Extend `solve_markov_income_at_prices(..., stored_bp=None, eval_mode=False)`
  to dispatch to full or eval mode and to return consistent timings.
- Extend `solve_markov_income_equilibrium` to keep a Markov `stored_bp`, count
  `n_full` and `n_eval`, and expose `force_full_bellman`, `howard_freq`, and
  a Markov warmup flag such as `markov_howard_warmup`.
- Extend `refine_one_market_markov_income` to accept a stored policy for cheap
  price evaluations, but force a full Bellman solve at the accepted/refined
  price before reporting final moments.

Policy objects to store:

- Main stored object: `bp_pol` with Markov shape
  `(Nb, nt, I, J, Nz, npar, ncs)`.
- Optional precomputed interpolation indices/weights for stored `bp_pol`, with
  the same shape, analogous to the non-Markov `stored_idx`/`stored_wt` branch.
- Do not freeze tenure, location, fertility, or child-aging probabilities.
  Recompute those objects under current prices after evaluating the stored
  savings policy.
- Preserve \(z\)-state indexing. For each current \(z\), continuation must use
  \[
  V^e_{j+1}(z)=\sum_{z'}\Pi_z(z,z')V_{j+1}(z')
  \]
  before child-aging and tenure/fertility choice. A non-Markov eval core cannot
  be reused without this continuation step.

Wiring:

- In each Markov GE iteration, do a full solve when `force_full_bellman` is
  true, no stored policy exists, the warmup condition holds, the Howard
  frequency triggers, or the market iteration stalls.
- Otherwise call Markov eval mode with `stored_bp`.
- After every full solve, update `stored_bp = bp_pol.copy()`.
- For scalar refinement, use eval mode for intermediate price probes only if a
  stored policy exists and the price is within a guarded relative distance from
  the full-solve price; otherwise use full mode. Always run one full Markov
  solve at the accepted refined price before final packing.
- Keep `force_full_bellman=True` as the parity/disable flag.

Tests:

1. Compile:

   ```bash
   python -m compileall -q code/model/intergen_housing_fertility
   ```

2. Basic smoke:

   ```bash
   cd code/model
   .venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
   ```

3. Fixed-price Markov parity:

   - Run a full Markov fixed-price solve at a representative price.
   - Re-run Markov eval mode at the same price using the full solve's `bp_pol`.
   - Require no invalid probabilities, no nonfinite policies, and key moments
     within `1e-5` absolute or `1e-4` relative. Market demand and distribution
     mass should be within `1e-6` to `1e-5`.

4. Full-only GE parity:

   - Run a small Markov GE case with `force_full_bellman=True`.
   - Run the new code with Howard effectively disabled (`howard_freq=1` or
     equivalent).
   - Require bit-level or near-bit-level equality up to JSON/float noise.

5. Howard GE parity:

   - Run the current representative short stack:
     `J=16`, `Nb=60`, `Nz=5`, `n_house=5`, `max_iter_eq=3`.
   - Compare full-only versus Howard-enabled results after the final accepted
     full solve.
   - Acceptable tolerances for this diagnostic scale: price relative error
     `<=1e-3`, market residual change `<=5e-5`, hard-target moment absolute
     differences `<=2e-3`, and objective loss difference `<=1e-2`.

6. Timing comparison:

   - Re-run the speed audit above.
   - Require `n_eval>0` in the aggregate timing record for configurations where
     Howard is intended to operate.
   - Require lower full Bellman count and lower wall time at the same grid,
     while the final accepted solution is recomputed in full mode.

Do not start a long calibration after this implementation. The next run after
approval should be a parity/timing packet, not a new production search.

## Staged Cleanup Plan, Not For This Pass

No folders should be moved as part of this SOTA task. A later cleanup can be
done in stages:

1. Clarify `code/model/README.md` so it distinguishes the older
   center-periphery package from the June 2026 one-market intergen package.
2. Refresh `memory/AGENT_MEMORY.md` to remove or qualify stale
   center-periphery "current working model" language.
3. Add an output index for the intergen diagnostic result folders, separating
   current June 2026 diagnostics from pre-fix and historical runs.
4. Archive only after review: stale drafts, obsolete generated diagnostics, and
   old exploratory outputs should be moved with a manifest, not ad hoc.
