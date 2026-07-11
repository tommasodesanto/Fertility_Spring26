# Calibration-search integrity worksheet — July 10 overnight combined spec

Auditor: search-integrity specialist, 2026-07-11. Local evidence only; SSH to torch blocked.
Audited objects: `tmp/overnight_combined_20260710/` (drivers), `code/model/intergen_housing_fertility/local_panel.py`,
`output/model/combined_recalibration/overnight_20260710_report/`.
Check script: `output/model/full_audit_20260711/scripts/check_search_integrity.py`
(output: `search_audit/check_search_integrity_output.json`). No model solves were run.

## 1. Bounds patching

Both drivers monkey-patch the module-level bounds at runtime:

- `run_overnight_combined.py:116` `lp.GLOBAL_DE_BOUNDS = list(bounds)`
- `run_wave2_lateral.py:52` `lp.GLOBAL_DE_BOUNDS = list(bounds)`

`bounds_for_arm("baseline_pattern")` returns EXACTLY `PRODUCTION_SEARCH_BOUNDS` — verified by
value comparison in the check script (`baseline_equals_production_bounds: true`). The construction
(`run_overnight_combined.py:61-78`) only substitutes entries listed in `replacements`, which is
empty for the baseline arm.

Relaxed boxes (all other entries identical to production; H0 in [1,10] appended everywhere):

| arm | c_bar_0 | psi_child | kappa_fert | theta_n | h_bar_0 | h_bar_jump | h_bar_n | chi |
|---|---|---|---|---|---|---|---|---|
| production | [0.08, 1.28] | [0, 0.35] | [1, 12] | [0, 1.5] | [1, 6] | [0.05, 2.5] | [0.02, 2] | [0.4, 1.15] |
| preference_relaxed | [0.04, 1.80] | [0, 0.60] | [0.20, 12] | [0, 2.50] | prod | prod | prod | prod |
| housing_relaxed | prod | prod | prod | prod | [0.25, 6] | [0.05, 3.50] | [0.02, 3.00] | [0.40, 1.60] |

Metadata fields under patching (`local_panel.py` `run_local_polish`):

- `meta["bounds"]` (lines 641-644) = `search_bounds` = patched `GLOBAL_DE_BOUNDS` + H0 → the ACTUAL
  search box. Correct.
- `meta["source_controlled_bounds"]` (lines 645-648): `for name, lo, hi in GLOBAL_DE_BOUNDS` — under
  patching this equals the patched box, so the label "source_controlled" is FALSE for the relaxed
  arms (the source-controlled box is `PRODUCTION_SEARCH_BOUNDS`). MISLABELED.
- `meta["production_profile_spec"]["bounds"]` (line 634 → `production_profile.py:116-119`) = unpatched
  `PRODUCTION_SEARCH_BOUNDS`.

So relaxed-arm `metadata.json` on scratch carries two contradictory boxes plus one mislabel.
Severity: MINOR for reproducibility — the actual box is correctly and unambiguously recorded in
`meta["bounds"]`; the contradiction only bites a reader who trusts `source_controlled_bounds`.

Report tables: `housing_relaxed_best_parameters.csv` lists the PATCHED (relaxed) bounds
(chi upper 1.6, h_bar_0 lower 0.25, h_bar_jump upper 3.5, h_bar_n upper 3.0) and
`current_bound_best_parameters.csv` lists production bounds + H0. i.e. the report lists the actual
per-arm search bounds. Correct.

Note the beta row: estimate 0.8411423842686492 with bounds [0.94, 0.995] looks contradictory but is
the annual/period transform: 0.9576732996982973^4 = 0.8411423842686492 exactly (verified);
`bound_scale_estimate` is annual. `theta_from_global_unit` (`local_panel.py:939-940`)
`theta["beta"] = value**PERIOD_YEARS`.

## 2. Wave-1 design (job 13313018)

`run_array.sh`: no `#SBATCH --array` directive — the range was supplied at submit time (array=1-24
appears in July-10 transcripts; consistent with `ARM_INDEX=(TASK_ID-1)/6`, 4 arms x 6 slots).
Arms `(baseline_pattern baseline_nelder preference_relaxed housing_relaxed)` (line 34); jitters
`(0 0.01 0.025 0.05 0.09 0.15)` per slot (line 31); `MAX_EVALS` default 850, `MINUTES` default 330
(lines 29-30); task seed `2026071100 + TASK_ID` (line 41) — distinct per task.

`jitter_theta` (`run_overnight_combined.py:81-94`) uses `np.random.default_rng(task_seed)`;
verified deterministic, distinct across seeds, and identity at jitter=0. Consequence: the four
jitter=0 tasks (one per arm) start at the SAME theta (modulo box clipping); distinct wave-1 starts
= 1 shared seed + 20 jittered = 21 points, all within jitter<=0.15 of one seed.

Seed provenance: `select_seed(args.seed_root)` (lines 45-58) recursively globs
`**/task_*/best.json` and takes the strict min of `rank_loss`. The root is `$OVERNIGHT_SEED_ROOT`,
set at submit time and NOT recorded in any local file — BLOCKED (scratch `overnight_launch.json`
per task records `source_seed_path`/`source_seed_loss`, unreachable). July-10/11 transcripts show a
prior strict tree `stage3_onehour_from_17_5912` (strict bests ~17.14+) and an incumbent
16.100969 quoted while wave 1 ran; the wave-1 root plausibly pointed at that stage-3 tree.

`baseline_nelder` reuses the baseline box (`run_overnight_combined.py:114`), method nelder-mead;
all other arms use pattern (line 118).

## 3. Wave-2 design (job 13314033)

`run_wave2.sh`: 3 arms x 8 tasks (ARM_INDEX=(TASK_ID-1)/8, line 21), task seed 2026071200+TASK_ID,
`--global-draws 60 --local-evals 720 --minutes 330` (lines 29-30).

Global stage (`run_wave2_lateral.py:78-117`): LHS over the full arm cube (14 dims incl. H0),
`units[0]` replaced by the wave-1 seed's unit image (line 83). Deadline
`min(100, 330*0.28)*60 = 92.4 min` (line 86), checked only for `idx > 0` (line 89), so the seed
always evaluates. Budget arithmetic: 60 draws need ~37 min at 37 s/solve, ~74 min at 74 s/solve —
fits inside 92.4 min even on 2x-slower nodes; truncation occurs only if the average case exceeds
~92 s (plausible for non-convergent random LHS points that burn `max_iter_eq=10` + scalar refine).
Whether any task truncated is recorded in `lateral_summary.json:global_completed` on scratch —
BLOCKED. Then pattern polish (initial_step 0.08) from the global-stage best (lines 121-152).

Critical silent-projection mechanic: `global_unit_from_theta` (`local_panel.py:946-966`) ends with
`return np.clip(unit, 0.0, 1.0)` — a seed theta OUTSIDE the arm box is silently clipped onto the
box boundary, never rejected (verified: `out_of_box_seed_silently_clipped: true`; `None` is
returned only for missing keys). `select_seed` picks the min across ALL wave-1 arms, including the
two arms wave-1 launch metadata itself marks `diagnostic_only: true`
(`run_overnight_combined.py:123`). So `global_current`'s `units[0]` can be (and, per Section 5,
almost surely is) a relaxed-arm point projected onto the production box.

## 4. Strictness

- `record_selection_loss` (`local_panel.py:1151-1155`): `return math.inf` unless
  `record.get("strict_converged")`; `is_better_record` (1132-1138) requires a finite selection loss.
  `strict_converged` is set in `run_local_panel_case` (863-867) as solver-strict AND
  `market_residual <= tol_eq`. So no non-strict record can become `best` in-run, regardless of its
  rank_loss — it is IGNORED, not selected and not a crash.
- `select_seed` filters on `strict_record` (`run_overnight_combined.py:53`), which additionally
  checks `timings["strict_converged"]`, finite residual, residual<=tol. Same in the reducer
  (`reduce_results.py:13-23`, `eligible` at 45). Reducer glob is one level (`*/task_*/best.json`),
  matching task-level best.json written by both waves.
- Crash-not-mask: RuntimeError fires only when NO strict record exists at a checkpoint
  (`run_wave2_lateral.py:116-117` global stage, `155-156` final; `run_overnight_combined.py:168-169`
  wave 1). Correction to the lead's framing: a strictly-better non-strict record does NOT crash the
  task; it is silently excluded. The crash path is "zero strict records".
- Latent MINOR bug: `run_wave2_lateral.py:154` `lp.is_better_record(polish_best, best)` with
  `polish_best=None` (polish produced no strict record) raises AttributeError inside
  `record_selection_loss(None)` (`local_panel.py:1152` `record.get`) before the intended
  RuntimeError at line 156. Fail direction is still a crash, so no result contamination; it can
  only kill a task whose polish seed (already strict at the global stage) fails to re-converge.

## 5. Shared-coordinate analysis / effective breadth

From `overnight_20260710_report/summary.json`: the two reported candidates share EXACTLY (16 digits)
beta=0.8411423842686492, c_bar_0=1.28, c_bar_n=0.4559787014242588, tenure_choice_kappa=0.0,
theta0=0.1318301350511569; and H0 differs by exactly 1.44 = 2 x (0.08 step x range 9). Verified
decompositions: alpha_cons diff = 1 x 0.08 x 0.55; kappa_fert diff = 1 x (0.08x0.6^4) x 11 — i.e.
pattern-polish coordinate steps at wave-2's initial_step 0.08 with shrink 0.6.

The two candidates come from different arms and different tasks (global_current/task_1 polish line
362; global_housing/task_17 polish line 207) with different task seeds, so their LHS draws are
independent — 16-digit coordinate equality is impossible for two independent LHS points. Therefore
in BOTH winning tasks the global-stage best was `units[0]`, the common wave-1 seed, and both
reported candidates are pattern-polish descendants of that single seed. The 14-D LHS stage
(<= 24 tasks x 59 = 1,416 dispersed draws) never produced a strict point better than the incumbent
lineage in any winning arm.

Box-boundary fingerprint of the seed: current_bound_best sits AT c_bar_0=1.28 (prod upper), chi=1.15
(prod upper), h_bar_0=1.0 (prod lower), while housing_relaxed_best has chi=1.2006 (>1.15) and
h_bar_0=0.921 (<1.0) interior to the relaxed box. h_bar_0<1 is representable only in the
housing_relaxed box, so the shared wave-1 seed almost surely came from the wave-1 housing_relaxed
arm — an arm whose own launch metadata says `diagnostic_only: true` — and `current_bound_best` is
that seed silently clipped onto the production box (Section 3 mechanic) and polished within it. Its
three at-bound parameters are exactly the clipped coordinates. SUSPECTED (decisive check on scratch:
`lateral_summary.json:source_seed_path` per wave-2 task).

Effective search breadth: <= ~39k evaluations capacity (wave 1: 24x850=20,400; wave 2:
24x(60+720)=18,720); `summary.json` reports 23,830 strict records. Genuinely dispersed points:
21 wave-1 starts (all within jitter 0.15 of one seed) + <=1,416 wave-2 LHS draws ~= 4% of all
evaluations; everything else is local descent along ONE lineage. 60 LHS draws per task in 14
dimensions is negligible coverage (2^14 = 16,384 orthants). The "global" stage functioned as a
sanity check that random points are worse, not as a global search; the report's implicit claim that
these are global optima over the stated boxes is NOT supported. The candidates remain valid,
strictly-converged local optima with correctly recorded losses (lead reproduced 14.780020699972585
exactly).

## 6. Override merge order (independent re-verification)

`run_local_panel_case` (`local_panel.py:848-853`):
```
overrides = {
    **base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq),
    **(extra_overrides or {}),
    **income,
    **theta,
}
```
- Wave-2 GLOBAL stage: `extra = production_profile_overrides()` then `extra.update({q, delta,
  eta_supply, normalize_bequest_utility})` (`run_wave2_lateral.py:68-76`); `income` is passed as the
  separate `income` argument (line 99 of the call at 92-103). Inside the merge, `**income` follows
  `**extra_overrides`, so the Rouwenhorst z_grid/Pi_z/persistence overwrite the profile's stale
  5-point grid, and `**theta` (incl. H0) comes last. Confirms the lead's finding for this
  differently-built path.
- Polish path (`run_local_polish`, lines 565-579): `profile_extra_overrides =
  production_profile_overrides()` then `["normalize_bequest_utility"]=...` then
  `.update(fixed_model_overrides)`; passed as `extra_overrides` with `income` separate (670-681).
  Same ordering, same conclusion. Fixed spec keys (q, delta, eta_supply) do not collide with income
  keys; theta keys do not collide with either.

## 7. Smoke runs

- Wave-1 smoke path exists: `run_array.sh:21-25` `OVERNIGHT_SMOKE=1` → 4 tasks (one per arm),
  MAX_EVALS=2, MINUTES=12, JITTER=0 — exercises the exact `run_local_polish` loop, strict gating,
  and metadata writing, per arm.
- Wave-2 smoke script exists: `run_wave2_smoke.sh` — 3 tasks (one per arm), `--global-draws 2
  --local-evals 2 --minutes 12`, seeds 2026071190+ID, 20-min walltime. It would validate the LHS +
  seed-injection + strict gate + polish + final-best plumbing, but NOT the 92.4-min global deadline
  or truncation behavior.
- Whether either smoke actually RAN: cluster-side (slurm logs / smoke outdirs on scratch) — BLOCKED.
  No smoke job id or completion line for the July-10 combined wave was found in local transcripts
  (the "24-worker smoke completed 0:0" line in the 2026-07-10 transcript is the July-9 repair wave,
  different snapshot).

## Blocked items (SSH unavailable)

1. `$OVERNIGHT_SEED_ROOT` values for both waves (what seeded wave 1; whether wave-2 root covered
   wave-1 results only or also stage-3).
2. Per-task `overnight_launch.json` / `lateral_summary.json` (`source_seed_path`,
   `source_seed_loss`, `global_completed` truncation check).
3. `cases.jsonl` files (distinct-vector counts, LHS-vs-seed loss distributions).
4. `cross_arm_summary.json` from the reducer (13314034).
5. Smoke-run execution evidence (slurm logs).
