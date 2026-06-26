# Handoff: intergen one-market model — audit fixes + calibration request

Date: 2026-06-25. For: the calibration agent. This is self-contained.

## TL;DR / what you're being asked to do

The one-market intergenerational housing-fertility model under
`code/model/intergen_housing_fertility/` was audited and **three measurement/code
bugs were found and fixed** (all committed + pushed to GitHub `main`). The model's
policies and accounting are now verified sound, and the targeted moments are
trustworthy. **Please run a 12-hour, 24-node calibration on Torch using two
different optimizer families.** Details in the last section. Pull `main` first and
re-sync the Torch scratch copy — the fixes change the loss, so prior incumbents
are stale.

## Model identity / canonical files

- Active model: `code/model/intergen_housing_fertility/` (one market, no location,
  Markov income, 4-year periods, J=16). One-shot completed-fertility (NOT
  sequential births — do not reinterpret `1to2` as a birth hazard).
- venv: `code/model/.venv/bin/python`
- Live config (set by `code/model/run_intergen_model.py` overrides, not
  `parameters.py` defaults): `Nb=60`, `income_states=5`,
  `H_own=[2,4,6,8,10]`, `hR_max=6.0`, `J=16`, `max_iter_eq=10`,
  `use_pti_constraint=False`, `tenure_choice_kappa=0.01`, `interp_method=linear`.
- Active target set: `candidate_replacement_young_old_roomgap_v1` (in
  `calibration.py`).
- Warm-start theta (13 internal params): 
  `output/model/intergen_room_distribution_current_best_20260623/summary.json`.
- Canonical status: `CALIBRATION_STATUS.md` (updated today). Implementation
  ledger: `code/model/intergen_housing_fertility/IMPLEMENTATION_STATUS.md`.

## Findings from the audit (so you know what's solid)

1. **The solver is NOT broken.** The renter savings policy looked jagged/
   non-monotone; a dense-grid global-optimum sweep proved the stored
   golden-section policy IS the global optimum (max value gap < 5e-3 over 992
   occupied nodes). `b'(b)` is monotone; the c/h "dips" are genuine
   target-saving toward future fertility/ownership thresholds and **persist under
   monotone-cubic interpolation and a 2x finer grid**, so they are real economics,
   not numerics. They also largely average out at the population (moment) level.
   **Do not change the optimizer; do not chase the policy shapes.**

2. **Markov room-moment bug — FIXED (commit `7219f64`).**
   `compute_markov_statistics` collapsed the income (z) dimension BEFORE applying
   nonlinear (>=threshold / cap / median) operators to the renter policy — a
   Jensen error (`1{E_z[h]>=t}` instead of `E_z[1{h>=t}]`). This made
   `prime30_55_childless_renter_share_rooms_ge6` read **0.013** when the true value
   is **0.124** (≈ target 0.138), and corrupted `prime_childless_renter_median_rooms`
   and the renter cap shares. Means and owner-rung moments were unaffected. Fixed
   via `markov_renter_room_moments()` (recompute on the full income-resolved
   distribution). **Implication: any prior cluster run / incumbent ranking that
   weighted renter room-share or median moments was optimizing against corrupted
   moments — recompute before trusting.**

3. **`housing_increment` event-study bug — FIXED (commit `8ade95f`).**
   The model differenced housing over a raw 12-year (3-period) window with NO
   control, so `housing_increment_0to1` read **1.467** — mostly lifecycle drift
   (ownership ramps 0->63% over ages 38-46), not the birth response. The PSID
   target (0.664) is a CONTROLLED Sun-Abraham room response ~3 years post-birth.
   Redefined (Markov path) as a difference-in-differences: birth cohort minus a
   no-birth control cohort over a configurable horizon `housing_event_horizon`
   (default **0** = the birth period, ≈ 3 years in a 4-year model). Now reads
   **0.477** vs target 0.664 (a modest under-shoot vs the old wild overshoot).
   `housing_increment_1to2` similarly moved to 0.246. Loss contribution of the
   0to1 moment drops from ~9.0 to ~0.5. (Horizon is tunable: 0->0.48, 1 (4yr)->0.40,
   2 (8yr)->-0.05.)

4. **Childless room moments are NOT a bug (verified).** The ACS "childless"
   targets condition on no children currently in household (`ni==0`,
   "childless-in-household"), which MATCHES the model's `current_child_bin_dt`
   (empty-nesters counted as childless in both). So `owner_share_rooms_ge6`
   (0.964 vs 0.596) and the room means are **genuine economic misses** (owners
   bunch in big units on a coarse ladder), not measurement. Keep as targets.

5. **Accounting is clean** (budgets, mass conservation, parity/TFR all to machine
   precision). No other Jensen siblings. `own_rate_2534≈0` is correctly measured
   (genuine economic miss). One latent, dormant code defect: the `H=2` owner rung
   has net service <=0 even childless (floored to 1e-10, ~1e-7 mass now) — fragile
   if `chi` falls or `h_bar_0` rises; not urgent, flag only.

## New infrastructure added (commit `8f97ed6`) — optional for you

- **Configurable wealth grid** (`make_grid`: `b_min/b_max/b_core_lo/b_core_hi/
  b_mid_hi/b_frac_*`; defaults reproduce the old grid). The current grid wastes
  ~24/60 nodes (mass lives in `b∈[-2.3, 9.1]`; deep-negative and far-positive are
  empty). A resized grid (`b_min≈-8, b_max≈25`, dense core) improves moment
  smoothness for derivative-based search but does NOT change the (genuine) policy
  dips. Diagnostic: `code/model/tools/grid_interp_lab.py --diagnose`.
- **`interp_method` switch** (`linear` default / `monotone_cubic` PCHIP).
  `monotone_cubic` routes off the compiled kernel to the slower Python path.
  **Keep `linear` for the calibration** (fast); cubic is for experiments only.

## Corrected moments snapshot (current point, after all fixes)

Confirmed by re-solve: aggregate loss **29.85** (down from 38.05; almost the
entire drop is the `housing_increment` fix, 9.01 -> 0.49). The dominant REMAINING
misses are all **economic**, and are what calibration should target:

| moment | target | model | note |
|---|---:|---:|---|
| `own_rate_2534` (young own) | 0.341 | ~0.000 | biggest miss; structural (see below) |
| `own_rate` | 0.575 | 0.371 | |
| `old_age_own_rate` | 0.764 | 0.916 | near-absorbing owner state |
| `owner_share_rooms_ge6` | 0.596 | 0.964 | owners bunch in big units |
| `young_liquid_wealth_to_income` | 0.179 | 0.610 | |
| `housing_increment_0to1` | 0.664 | **0.477** | FIXED (was 1.467) |
| `housing_increment_1to2` | 0.488 | **0.246** | FIXED |
| `prime30_55_childless_renter_share_rooms_ge6` | 0.138 | **0.124** | FIXED (was 0.013) |
| `tfr` | 1.700 | 1.710 | hit |
| `childless_rate` | 0.150 | 0.284 | |

Economic context (do NOT fix via mechanism changes — out of scope this round):
young ownership ≈0 and old ownership ≈0.92 are the dominant misses and are
structural (rent is pinned at user cost, so owning has no flow advantage for the
young beyond the small `chi` premium; old ownership is near-absorbing). Calibration
alone likely cannot fully close these without a new mechanism (parental transfers /
old-age exit shock), which is a separate future decision. **For this run, just find
the best SMM fit given the current structure.**

## Calibration request

Please run a **12-hour, multi-node SMM calibration on Torch**:

- **24 nodes**, **two different optimizer families** (e.g. global differential
  evolution AND the seeded random/global + local-adaptive search) so we get
  algorithmic variation — split the 24 nodes between them as you see fit.
- Target set `candidate_replacement_young_old_roomgap_v1`; grid `J=16, Nb=60,
  income_states=5, H_own=[2,4,6,8,10], hR_max=6.0, max_iter_eq=10`,
  `interp_method=linear`.
- Warm-start / seed from
  `output/model/intergen_room_distribution_current_best_20260623/summary.json`
  (13 internal params), but also let the global searches explore broadly.
- Discipline (project long-run-search rules): smoke-test the exact loop first;
  write progress + checkpoints at least every case / 5 min; keep best-so-far and
  latest-case artifacts; record the output dir, target set, starting point, and
  exact launch command in `CALIBRATION_STATUS.md`.
- **Before trusting any saved incumbent, recompute its loss with the fixed stats**
  — prior results that weighted renter room-share/median or `housing_increment`
  are stale.

Verification before launch:
```bash
python -m compileall -q code/model/intergen_housing_fertility
cd code/model && .venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
```
Confirm `git log --oneline -5` shows commits `7219f64`, `8f97ed6`, `8ade95f` on
`main`, and re-sync the Torch scratch copy to that commit.
