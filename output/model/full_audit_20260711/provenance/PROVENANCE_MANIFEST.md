# Provenance Manifest — July 10 Overnight Combined Recalibration

Audit date: 2026-07-11. Specialist: provenance. Companion machine-readable file: `hashes.json` (same directory).

## (a) Git state

- HEAD: `fe092ca99030e54fb256206725e66d9694dc539d`, 2026-07-10 10:26:47 -0400, Tommaso De Santo, "Add repaired equilibrium figures to preliminary note".
- Branch: `main`.
- Working tree: DIRTY. The dirty tree is the audited object; the overnight ran from a torch snapshot of it (`/scratch/td2248/projects/Fertility_Spring26_20260710_combined_spec`). **Torch-snapshot diff is BLOCKED: no SSH to the cluster during this audit.** Local-tree ≡ snapshot is asserted by the lead's shared context, not verified byte-for-byte.

### Modified tracked files (model/calibration-relevant)

`code/model/intergen_housing_fertility/`: calibration.py, cli.py, local_panel.py, parameters.py, solver.py.
`code/model/dt_cp_model/`: direct_calibration.py, objective.py (dt_cp strand, not load-bearing for the intergen candidate).
`code/model/tools/`: audit_intergen_final_best_pathologies.py, audit_intergen_mechanism_grid.py, audit_intergen_parent_credit_margin.py, audit_intergen_solver_accuracy.py, audit_intergen_wealth_units.py, build_intergen_mechanics_packet.py, calibrate_direct_geometry.py, make_direct_fit_slide_note.py, make_full_dispersion_diagnostic_note.py, plot_ownership_spatial_lifecycle_diagnostic.py, run_intergen_policy_poc.py.
`code/model/README.md` (J 16→17 in examples).
`code/cluster/`: submit_intergen_housing_fertility_{calibration,global_de,local_polish,twohour_panel}.sh, submit_intergen_mechanism_grid.sh, submit_intergen_mechanism_repair_grid.sh, submit_intergen_sensitivity_jacobian.sh, submit_python_direct_geometry_overnight.sh, torch.sh — all mechanical: J default 16→17 plus new pass-through env vars `INTERGEN_INCOME_PROCESS/_ANNUAL_RHO/_INNOVATION_SD` (default `current`, 0.90, 0.20).

### Untracked files (model/calibration-relevant)

- `code/model/intergen_housing_fertility/tests/test_bequest_normalization.py`, `tests/test_combined_specification.py`
- `code/model/tools/run_intergen_combined_recalibration.py` (the repo driver for the combined spec)
- `code/model/tools/compare_intergen_bequest_normalization.py`, `compare_intergen_combined_specification.py`, `compare_intergen_income_processes.py`
- `code/cluster/submit_intergen_combined_recalibration.sh`, `submit_intergen_combined_spec_test.sh`
- ALL of `tmp/overnight_combined_20260710/` (the entire `tmp/` directory is untracked; the actual overnight drivers live only here and on the unreachable snapshot)
- `output/model/combined_recalibration/overnight_20260710_report/` (git-ignored via `output/` in .gitignore): current_bound_best{.json,_fit.csv,_parameters.csv}, housing_relaxed_best{...}, summary.json — the result artifacts themselves are outside version control.
- Other untracked tools (policy counterfactuals, phi sweeps, redfin/house-size plots, dk_lm tests, `code/model/intergen_surrogate_calibration/`) are listed in `hashes.json`-adjacent git status but not hashed; they are not on the candidate's evaluation path.

Clean and tracked (important negatives): `production_profile.py` (last commit 2b988b3 "Restore exact July calibration runtime"), `solver.py`-adjacent `kernels.py`, `diagnostics.py`, `utils.py`, `__init__.py`, `tests/test_code_repair_20260709.py` — all unmodified vs HEAD.

## (b) SHA256 hashes

See `hashes.json` (37 files: full intergen package incl. tests, the repo driver, the three compare_intergen_* tools, both combined submit scripts, all 7 files in `tmp/overnight_combined_20260710/`, and the 7 overnight report artifacts under `output/model/combined_recalibration/overnight_20260710_report/`). Each entry records sha256, byte size, and git status (tracked_clean / modified_vs_HEAD / untracked / ignored).

## (c) What changed vs HEAD in the modified tracked files

Load-bearing for the 14.780 candidate (all in `code/model/intergen_housing_fertility/`):

- **local_panel.py** (+147/-?): (1) NEW Rouwenhorst branch in `income_process_overrides(income_states, process, annual_innovation_sd, annual_rho)` — builds period rho `rho_annual**PERIOD_YEARS`, period innovation sd via `innovation_sd * sqrt(sum(rho^{2k}))`, binomial weights, mean-one z grid, Rouwenhorst Pi_z. (2) **Override merge-order flip in `run_local_panel_case`** (line ~847-852): HEAD had `**income,` before `**(extra_overrides or {})`; dirty tree has `**(extra_overrides or {}), **income,` — income now wins over the production profile's stale 5-point z_grid/Pi_z. Economically substantive: at HEAD the profile grid would silently override Rouwenhorst. (3) `run_local_polish` gains `normalize_bequest_utility`, `additional_search_bounds`, `additional_targets`, `additional_weights`, `fixed_model_overrides` kwargs — the entire H0-search + rooms-target + fixed-spec plumbing the overnight used. (4) Audit fingerprint `income_process_fingerprint` and `Pi_z` logged in headers. All substantive.
- **solver.py** (+11/-2, `bequest_utility_vec` lines 5608-5620): adds `normalize_bequest_utility` switch subtracting the b=0 utility level: log case `utility - log(theta1)`, CRRA case `utility - theta1**(1-sigma)/(1-sigma)`. Economically substantive (changes bequest motive level; overnight ran with it True).
- **calibration.py** (+4): `extract_moments` adds `"aggregate_mean_occupied_rooms_18_85": aggregate_housing_demand / max(total_mass, 1e-12)` — the new 15th moment. Substantive.
- **parameters.py** (+3/-1): default `P.normalize_bequest_utility = False`; `apply_overrides` reshapes overridden `H0` to 1-D (`np.asarray(...).reshape(-1)`) — needed because H0 is now a searched scalar. Substantive-supporting.
- **cli.py** (+20): new CLI flags `--income-process/--income-annual-rho/--income-innovation-sd/--normalize-bequest-utility` threaded to the panel/DE/polish entry points. Mechanical plumbing.

Not load-bearing for the candidate:

- **dt_cp_model/objective.py** (+56): adds `_attach_housing_expenditure_share_moments` (user-cost housing share) — dt_cp strand. **dt_cp_model/direct_calibration.py** (+16): optional alpha_cons calibration with bounds [0.80,0.90].
- **tools/**: audit_intergen_final_best_pathologies.py, audit_intergen_mechanism_grid.py, audit_intergen_solver_accuracy.py — J 16→17 one-liners. audit_intergen_wealth_units.py — reports new entry-wealth objects. audit_intergen_parent_credit_margin.py (+462) — rewritten to the loss6 record, optional matplotlib, fixed-price PE mode, J=17. build_intergen_mechanics_packet.py — adds `--combined-corrected-spec` flag that hard-codes the same combined spec constants (rooms target 5.779970481941968 w=6, Rouwenhorst 0.9601845894041878/0.06453733259357768, q=(1.02)^4-1). run_intergen_policy_poc.py — J=17 default, target-set arg, new LTV95 policy cases. calibrate_direct_geometry.py / make_direct_fit_slide_note.py / make_full_dispersion_diagnostic_note.py / plot_ownership_spatial_lifecycle_diagnostic.py — dt_cp housing-share target, record-vs-fresh moment checks, plot styling. None enter the intergen candidate's evaluation path, but build_intergen_mechanics_packet.py and run_intergen_policy_poc.py will be used to interpret it downstream.

## (d) Repo driver vs overnight drivers

Repo driver `code/model/tools/run_intergen_combined_recalibration.py` (untracked) vs `tmp/overnight_combined_20260710/run_overnight_combined.py` and `run_wave2_lateral.py` (untracked). Shared: identical constants (MATCHED_ANNUAL_RHO 0.9601845894041878, MATCHED_ANNUAL_INNOVATION_SD 0.06453733259357768, ROOMS_TARGET 5.779970481941968, weight 6.0, H0 in [1,10], q=(1.02)^4-1, delta=1-(1-0.011)^4, eta_supply=[1.75], normalize_bequest_utility=True), same target set PRODUCTION_TARGET_SET, J=17, Nb (PRODUCTION_SEARCH_NB=120 vs literal 120), n_house=5, max_iter_eq=PRODUCTION_MAX_ITER_EQ=10, income_states=5.

Differences in controls:

| control | repo driver | run_overnight_combined.py (wave 1) | run_wave2_lateral.py (wave 2) |
|---|---|---|---|
| seed source | `--seed-theta` JSON file | `select_seed()`: best strict `best.json` under `--seed-root` | same `select_seed()` from wave-1 tree |
| evals cap | `min(args.max_evals, 160)`; default env COMBINED_MAX_EVALS=60 | `--max-evals` default 1600, launched at 850 (run_array.sh), NO hard cap | 60 LHS global draws + `--local-evals` 720 polish |
| minutes cap | `min(args.minutes, 60.0)`; default env 220 (so 60 effective) | `--minutes` default 690, launched at 330, NO cap | 330 total; global stage deadline `min(100, minutes*0.28)` min; polish gets remainder minus 2 |
| bounds | production GLOBAL_DE_BOUNDS unchanged + [("H0",1,10)] | `bounds_for_arm()` per arm; **mutates module global `lp.GLOBAL_DE_BOUNDS = list(bounds)`**; relaxed arms: preference_relaxed widens c_bar_0 (0.04,1.80), psi_child (0,0.60), kappa_fert (0.20,12), theta_n (0,2.50); housing_relaxed widens h_bar_0 (0.25,6), h_bar_jump (0.05,3.50), h_bar_n (0.02,3), chi (0.40,1.60); metadata marks relaxed arms `diagnostic_only: true` | same `bounds_for_arm` via arm map {global_current→baseline_pattern, global_preference→preference_relaxed, global_housing→housing_relaxed}; **no diagnostic_only marker in wave-2 metadata** |
| jitter | none | `--jitter` per task: unit-space N(0, jitter) clipped to [0,1]; slots (0, 0.01, 0.025, 0.05, 0.09, 0.15) | none (LHS draws instead; `units[0]` replaced by the wave-1 seed) |
| method | pattern, initial_step arg (env, default 0.04), min_step 0.003, shrink 0.5 | pattern or nelder-mead by arm; initial_step 0.04 (NM) / 0.06 (pattern), min_step 0.0015, shrink 0.60 | pattern only, initial_step 0.08, min_step 0.0015, shrink 0.60 |
| seeds | env COMBINED_SEED, default 20260710; submit script sets 20260710+1000*TASK_ID | `--task-seed` = 2026071100+TASK_ID (run_array.sh) | 2026071200+TASK_ID (run_wave2.sh); polish seed = task_seed+100000; smoke 2026071190+TASK_ID |
| Nb | PRODUCTION_SEARCH_NB (=120) via env COMBINED_NB | literal 120 | literal 120 |
| profile | PRODUCTION_PROFILE_NAME unless env COMBINED_PROFILE!=1 | PRODUCTION_PROFILE_NAME always | PRODUCTION_PROFILE_NAME for polish; global stage calls `run_local_panel_case` directly with `extra = production_profile_overrides()` + fixed spec |
| failure gate | requires best strict_converged | `strict_record()` (also checks timings.strict_converged and residual<=tol) | same strict gate on global stage and final |

Launch controls from the shell scripts (all `cpu_short`, account torch_pr_570_general, 1 CPU, threads pinned to 1, PYTHONPATH=$SNAPSHOT/code/model):

- `run_array.sh`: wave 1, 05:45:00 walltime, 4G; arms = (baseline_pattern, baseline_nelder, preference_relaxed, housing_relaxed) × 6 jitter slots ⇒ 24 array tasks; MAX_EVALS=850, MINUTES=330 (env-overridable); smoke mode (OVERNIGHT_SMOKE=1): 2 evals, 12 min, jitter 0.
- `run_wave2.sh`: wave 2, 05:45:00, 4G; arms = (global_current, global_preference, global_housing) × 8 ⇒ 24 tasks; --global-draws 60 --local-evals 720 --minutes 330; PYTHONPATH additionally prepends $ROOT so wave 2 can import run_overnight_combined.
- `run_wave2_smoke.sh`: 00:20:00; 3 tasks, seeds 2026071190+ID, 2 draws / 2 evals / 12 min.
- `reduce.sh`: 00:20:00, runs reduce_results.py; reducer takes cross-arm min over strict records only, writes cross_arm_summary.json; no model call.
- Repo submit script `code/cluster/submit_intergen_combined_recalibration.sh`: 04:00:00, 8G, COMBINED_MAX_EVALS=50, COMBINED_MINUTES=235 (but driver caps at 60 min), step ladder (0.020…0.090) by task — this is the earlier short seeding run, NOT the overnight.

Notable design flags for the search-design specialist: (i) relaxed-bounds arms are declared diagnostic-only in wave 1 but wave 2 reuses them without that marker, and `housing_relaxed_best` is reported next to `current_bound_best` in the report artifacts; (ii) `lp.GLOBAL_DE_BOUNDS` is mutated at module level by both tmp drivers; (iii) wave 2's `--minutes 330` on a 05:45 walltime leaves ~15 min slack; (iv) the repo driver could never have produced an 850-eval overnight run (hard 160/60 caps) — the tmp drivers are materially different programs, not copies.

## (e) Verdict: reproducible from COMMITTED source alone?

**No.** From `git HEAD fe092ca` alone the 14.780 candidate cannot be reproduced, for four independent reasons:

1. **The evaluation code is uncommitted.** Five working-tree modifications to tracked files are load-bearing: `local_panel.py` (Rouwenhorst branch, `run_local_polish` extended signature, override merge order), `solver.py` (`normalize_bequest_utility`), `calibration.py` (`aggregate_mean_occupied_rooms_18_85` moment), `parameters.py` (`normalize_bequest_utility` default, H0 reshape), `cli.py` (plumbing). At HEAD, `run_local_polish` does not even accept `additional_search_bounds` / `additional_targets` / `fixed_model_overrides` / `normalize_bequest_utility` — the drivers would crash with TypeError.
2. **Even the parts that would run would compute a different economy at HEAD.** At HEAD `run_local_panel_case` merges `**income` BEFORE `**(extra_overrides)`, so the production profile's stale 5-point z_grid/Pi_z (`production_profile.py:136-141`, which is tracked and clean) would override any Rouwenhorst income; the dirty tree flips this order so Rouwenhorst wins. The lead verified the candidate is only reproduced under the dirty-tree order.
3. **The drivers are untracked.** `code/model/tools/run_intergen_combined_recalibration.py`, both combined submit scripts, and the entire `tmp/overnight_combined_20260710/` directory (`tmp/` is wholly untracked) exist only in the working tree and on the unreachable torch snapshot.
4. **The result artifacts are git-ignored.** `output/` is ignored, so `output/model/combined_recalibration/overnight_20260710_report/*` (current_bound_best.json etc.) is outside version control; only the SHA256 hashes in `hashes.json` pin them.

Reproducibility from the CURRENT DIRTY TREE is separately established (lead reproduced loss 14.780020699972585 exactly via `output/model/full_audit_20260711/scripts/reproduce_candidate.py`). The correct remediation is to commit the five modified package files, the two new tests, the repo driver + submit scripts, and archive copies of the three tmp drivers + four launch scripts, then tag that commit as the provenance of the July 10 overnight.

Additional observation for other specialists: the merge-order flip was applied only to `run_local_panel_case` (local_panel.py:849-853). A second solve site, `write_best_case_diagnostics` (local_panel.py:1303-1307), merges only `base_overrides + income + theta` — no profile extra_overrides and no fixed spec (q/delta/eta/normalize_bequest_utility). It is called from `run_local_panel` (line 264), not from `run_local_polish`, so it did not touch the overnight loss; but any diagnostics packet generated through the panel path would be solved under a different specification than the ranked record.

Blocked check: byte-identity of the torch snapshot `/scratch/td2248/projects/Fertility_Spring26_20260710_combined_spec` vs this working tree (SSH unavailable). Until verified, the possibility of snapshot-side divergence remains open, though the exact-loss local reproduction makes material divergence on the evaluation path unlikely.
