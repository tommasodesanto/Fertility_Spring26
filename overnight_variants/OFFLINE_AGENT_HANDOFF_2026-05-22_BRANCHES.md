# Offline Agent Handoff: Overnight Branches and V5 Global Search

This document is meant to be self-contained for an agent that will not have the
chat transcript. It summarizes what was done, what was decided, which files
matter, and how to resume from the live global search.

Repository:

```text
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
```

Active isolated work folders:

```text
overnight_variants/2026-05-22_income_mortgage_risk/
overnight_variants/2026-05-22_developer_missing_middle/
```

Do not edit live model code unless the user explicitly asks:

```text
code/model/dt_cp_model/
```

## Required Startup for Next Agent

Before making claims or changes, read these files:

```text
AGENTS.md
memory/AGENT_MEMORY.md
CALIBRATION_STATUS.md
latex/model_writeup.tex
overnight_variants/REPORT_COMPARISON.md
overnight_variants/AGENT_HANDOFF_2026-05-22_BRANCHES.md
overnight_variants/OFFLINE_AGENT_HANDOFF_2026-05-22_BRANCHES.md
```

Then read branch-specific reports:

```text
overnight_variants/2026-05-22_income_mortgage_risk/README.md
overnight_variants/2026-05-22_income_mortgage_risk/REPORT_V5_HANK_Z_OUTSIDE_CLOSURE.md
overnight_variants/2026-05-22_income_mortgage_risk/REPORT_V5_HANK_Z_OUTSIDE_CLOSURE_NZ5.md
overnight_variants/2026-05-22_income_mortgage_risk/REPORT_V5_SPEED_AUDIT.md
overnight_variants/2026-05-22_income_mortgage_risk/REPORT_V5_NZ5_DIRECTIONAL_CALIBRATION.md
overnight_variants/2026-05-22_income_mortgage_risk/REPORT_V5_NZ5_GLOBAL_SEARCH_PLAN.md
overnight_variants/2026-05-22_developer_missing_middle/README.md
overnight_variants/2026-05-22_developer_missing_middle/REPORT_V3_GE.md
```

## Git and Dirty-Tree Warning

There are pre-existing dirty files outside `overnight_variants/`. They were
not created by the overnight branch work and should not be touched, staged, or
reverted unless the user explicitly asks.

At the time this document was written, known unrelated dirty paths included:

```text
code/model/docs/architecture.html
code/model/docs/tour.html
code/model/dt_cp_model/solver.py
code/model/tools/tour_template.html
latex/april_20_project_presentation.pdf
latex/april_20_project_presentation.tex
latex/main_note.pdf
latex/main_note.tex
latex/model_writeup.pdf
latex/model_writeup.tex
latex/two_branch_extension_memo_2026_05_22.pdf
latex/two_branch_extension_memo_2026_05_22.tex
```

The live global-search output directory is intentionally untracked while the
run is active:

```text
overnight_variants/2026-05-22_income_mortgage_risk/global_search_v5_nz5/
```

Do not commit partial global-search outputs unless the user asks. When the run
finishes, ask whether to commit the final search outputs or only a summary.

## High-Level Decision

Recommendation at handoff:

```text
Branch 1 HANK-z income risk should be the priority, but only after calibration
evidence improves and an Nz=7 validation run confirms the economics.
```

Reasons:

- Branch 1 is the canonical one-state income-risk extension: one structural
  \(z\) state, standard Rouwenhorst process, and the paper-facing
  outside-option closure.
- Branch 1 now clears strict GE and responds to calibration levers.
- Branch 2, the developer missing-middle branch, is mechanically useful but
  does not solve the room-target or ownership-gradient problem.
- Do not add the structural mortgage-account state \(\mu\) yet. The current
  evidence says recalibrating the HANK-\(z\) economy is the next bottleneck.

## Branch 1: Income Risk / HANK-z / Outside Closure

Folder:

```text
overnight_variants/2026-05-22_income_mortgage_risk/
```

Status:

```text
Mechanically working; economically yellow; active Nz=5 global search running.
```

### What Branch 1 Is

The branch is an isolated copied-solver prototype. It is not live code.

Implemented in the copied solver:

- A structural idiosyncratic earnings state \(z\).
- Household state is \((b,d,i,a,n,s,z)\).
- Value functions, policies, fertility probabilities, location probabilities,
  tenure choices, and forward distributions all carry \(z\).
- The \(z\) process is a Rouwenhorst Markov chain.
- Serious validation grid: \(N_z=7\), \(\rho_z=0.95\), unconditional
  \(\sigma_z=0.35\).
- Exploratory/calibration grid: \(N_z=5\), same \(\rho_z\) and \(\sigma_z\).

Not implemented:

- There is no structural mortgage-account state \(\mu\) in the Bellman
  recursion.
- The older \(\mu\) objects in the branch are diagnostic only.

### Closure Audit and Decision

The user asked whether the V5 outside-option closure should use an outer
normalization pass or impose the benchmark normalization directly in the GE
loop.

The relevant paper equation is:

\[
S E_0(p) = q^E(p) [M + S B_0(p)].
\]

Benchmark normalization is \(S=1\). Therefore, in the benchmark loop, the
correct implementation is:

\[
M = E_0(p)/q^E(p) - B_0(p).
\]

Decision:

- The benchmark should impose \(S=1\) inside each GE iteration.
- At the current candidate equilibrium/composition, compute \(E_0(p)\) and
  \(B_0(p)\).
- Calibrate the outside value \(\bar W^E\) to target \(q^E=0.9\).
- Set \(M\) residually from the expression above.
- Hold final \(M\) and \(\bar W^E\) fixed for counterfactuals.
- The earlier outer preliminary solve and re-solve was removed.

This was implemented only in the isolated copied solver under:

```text
overnight_variants/2026-05-22_income_mortgage_risk/dt_cp_model/solver.py
overnight_variants/2026-05-22_income_mortgage_risk/run_income_mortgage_risk_v5_hank_z_outside_closure.py
```

Important solver concepts added:

- `outside_option_benchmark_normalized`
- `BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES`
- `uses_benchmark_normalized_outside_closure`
- `benchmark_normalized_outside_population_scale`

### Branch 1 Important Files

Scripts:

```text
run_income_mortgage_risk_v5_hank_z_outside_closure.py
run_v5_speed_audit.py
plot_income_mortgage_risk_v5_hank_z_outside_closure_nz5.py
run_v5_nz5_directional_calibration.py
run_v5_nz5_global_search.py
```

Reports:

```text
REPORT_V5_HANK_Z_OUTSIDE_CLOSURE.md
REPORT_V5_HANK_Z_OUTSIDE_CLOSURE_NZ5.md
REPORT_V5_SPEED_AUDIT.md
REPORT_V5_NZ5_DIRECTIONAL_CALIBRATION.md
REPORT_V5_NZ5_GLOBAL_SEARCH_PLAN.md
```

Results and diagnostics:

```text
results_income_mortgage_risk_v5_hank_z_outside_closure.csv
diagnostics_income_mortgage_risk_v5_hank_z_outside_closure.csv
diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_closure.csv
diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_trace.csv
results_income_mortgage_risk_v5_hank_z_outside_closure_nz5.csv
diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_nz5.csv
diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_nz5_closure.csv
diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_nz5_trace.csv
```

Figures:

```text
figures_v5_hank_z_outside_closure_nz5/HANK_Z_OUTSIDE_CLOSURE_NZ5_FIGURE_PACKET.pdf
```

### V5 Nz=7 Validation Result

Command:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v5_hank_z_outside_closure.py --quiet --nb 30 --nz 7 --rho-z 0.95 --sigma-z 0.35 --kappa-entry 1000000 --baseline-max-iter-eq 35 --max-iter-eq 60 --normalization-passes 2 --scale-target-tol 1e-5 --tol-eq 5e-4
```

Result:

```text
accepted: True
strict converged: True
convergence reason: strict_tol
iterations: 17
runtime: 620.134208930991s
final GE error: 0.00022218249381993278
prices: [0.5195066598297923, 0.6004039929449273]
S: 1.0000000000000009
qE: 0.9000000000000001
outside probability: 0.1
M: 0.005780247614275978
outside value: -3327810.8839415493
SMM loss: 311.4039737325341
```

Moment table:

| Moment | Target | Nz=7 V5 |
|---|---:|---:|
| `tfr` | 1.700 | 1.562 |
| `childless_rate` | 0.150 | 0.351 |
| `mean_age_first_birth` | 26.000 | 35.280 |
| `tfr_gradient` | 0.133 | -0.194 |
| `own_rate` | 0.627 | 0.609 |
| `own_gradient` | 0.170 | -0.129 |
| `own_family_gap` | 0.110 | 0.353 |
| `prime_childless_renter_median_rooms` | 4.000 | 6.002 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 |
| `housing_increment_0to1` | 0.664 | 0.588 |
| `housing_increment_1to2` | 0.566 | 1.489 |
| `young_liquid_wealth_to_income` | 0.600 | 1.866 |
| `center_share_nonparents` | 0.494 | 0.233 |
| `center_share_newparents` | 0.416 | 0.371 |
| `migration_rate` | 0.032 | 0.033 |
| `old_age_own_rate` | 0.863 | 0.753 |
| `old_age_parent_childless_gap` | 0.070 | 0.263 |
| `inv_pop_share_C` | 0.450 | 0.385 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.156 |

Read:

- Mechanical closure is good.
- Economics are not good yet.
- \(N_z=7\) is too slow for exploratory calibration.

### Runtime Audit

Important distinction:

```text
Warm means computational warm-up in the same Python process, mainly Numba
kernels already compiled for the same array shapes. It does not mean economic
warm start. No policies, prices, values, or distributions were reused.
```

Cold accepted `Nb=30,Nz=7`:

```text
620.13s
```

Same-process warm accepted `Nb=30,Nz=7`:

```text
200.82s
```

Standalone `Nb=30,Nz=5` accepted plotting run:

```text
83.10s including solve plus plot packet
```

Parallel 8-worker global search observed first-wave solve times:

```text
about 230-550s per candidate under contention
```

Expected 8-hour global-search size:

```text
roughly 800-1000 completed full Nz=5 evaluations
```

### V5 Nz=5 Plotting Result

Command:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python plot_income_mortgage_risk_v5_hank_z_outside_closure_nz5.py --quiet --nb 30 --nz 5 --rho-z 0.95 --sigma-z 0.35 --kappa-entry 1000000 --max-iter-eq 60 --tol-eq 5e-4
```

Result:

```text
accepted: True
strict converged: True
convergence reason: strict_tol
iterations: 12
runtime: 83.10237065298134s
final GE error: 0.0003464318147695515
prices: [0.511963296326644, 0.6124693819761028]
S: 0.9999999999999998
qE: 0.8999999999999999
outside probability: 0.09999999999999998
M: 0.005806434956255813
outside value: -1504137.9225856075
loss: 311.29977174646774
```

Moment table:

| Moment | Target | Nz=5 V5 |
|---|---:|---:|
| `tfr` | 1.700 | 1.558 |
| `childless_rate` | 0.150 | 0.351 |
| `mean_age_first_birth` | 26.000 | 34.793 |
| `tfr_gradient` | 0.133 | -0.239 |
| `own_rate` | 0.627 | 0.575 |
| `own_gradient` | 0.170 | -0.082 |
| `own_family_gap` | 0.110 | 0.383 |
| `prime_childless_renter_median_rooms` | 4.000 | 5.767 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 |
| `housing_increment_0to1` | 0.664 | 0.591 |
| `housing_increment_1to2` | 0.566 | 1.492 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 |
| `center_share_nonparents` | 0.494 | 0.316 |
| `center_share_newparents` | 0.416 | 0.393 |
| `migration_rate` | 0.032 | 0.034 |
| `old_age_own_rate` | 0.863 | 0.709 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 |
| `inv_pop_share_C` | 0.450 | 0.409 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.196 |

Read:

- Same qualitative failures as \(N_z=7\).
- Good for plots and exploratory calibration.
- Do not make substantive economic claims from \(N_z=5\) alone.

### Directional Calibration Audit

Script:

```text
run_v5_nz5_directional_calibration.py
```

Reports and outputs:

```text
REPORT_V5_NZ5_DIRECTIONAL_CALIBRATION.md
v5_nz5_directional_calibration.csv
v5_nz5_directional_calibration_moments.csv
v5_nz5_directional_calibration_round2.csv
v5_nz5_directional_calibration_round2_moments.csv
```

Best probe:

```text
case: alpha_high_fertility_mid_finance
alpha_cons: 0.80
kappa_fert: 4.0
phi: 0.90
accepted: True
strict converged: True
iterations: 16
runtime: 107.72640849198797s
loss: 166.7224914373928
final GE error: 0.00041013185652335614
S: 1
qE: 0.9
outside probability: 0.1
M: 0.004853356989326867
outside value: -1504146.3877423254
```

Best-probe moment table:

| Moment | Target | Baseline Nz=5 | Best Probe |
|---|---:|---:|---:|
| `tfr` | 1.700 | 1.558 | 1.673 |
| `childless_rate` | 0.150 | 0.351 | 0.310 |
| `mean_age_first_birth` | 26.000 | 34.793 | 34.445 |
| `tfr_gradient` | 0.133 | -0.239 | -0.251 |
| `own_rate` | 0.627 | 0.575 | 0.644 |
| `own_gradient` | 0.170 | -0.082 | 0.046 |
| `own_family_gap` | 0.110 | 0.383 | 0.184 |
| `prime_childless_renter_median_rooms` | 4.000 | 5.767 | 5.092 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 | 5.400 |
| `housing_increment_0to1` | 0.664 | 0.591 | 0.693 |
| `housing_increment_1to2` | 0.566 | 1.492 | 1.375 |
| `young_liquid_wealth_to_income` | 0.600 | 1.215 | 1.322 |
| `center_share_nonparents` | 0.494 | 0.316 | 0.342 |
| `center_share_newparents` | 0.416 | 0.393 | 0.418 |
| `migration_rate` | 0.032 | 0.034 | 0.034 |
| `old_age_own_rate` | 0.863 | 0.709 | 0.727 |
| `old_age_parent_childless_gap` | 0.070 | 0.291 | 0.199 |
| `inv_pop_share_C` | 0.450 | 0.409 | 0.431 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.196 | 1.182 |

Interpretation:

- The model is movable.
- `alpha_cons` is important; fixed `alpha_cons=0.70` appears too restrictive
  once HANK-\(z\) risk is active.
- `kappa_fert` and `phi` are useful near-term search parameters.
- Lowering `beta` helps wealth and timing but damages ownership; do not lower
  it aggressively without an offsetting tenure channel.
- Geography shifters are powerful but can destabilize rent ratio and gradients;
  treat geography carefully.
- Remaining hard problems are first-birth age, childlessness, fertility
  gradient, and young liquid wealth.

## Live V5 Nz=5 Global Search

Status at this document's creation:

```text
run tag: v5_nz5_global_20260522_215756
completed: 36
active: 8
submitted: 44
best objective: 166.7224914373928
best eval id: 2
finished: false
```

This snapshot will be stale as the run continues. Always read `status.json`.

Run directory:

```text
overnight_variants/2026-05-22_income_mortgage_risk/global_search_v5_nz5/v5_nz5_global_20260522_215756/
```

Screen sessions:

```text
v5_nz5_global_20260522_215756
v5_nz5_global_20260522_215756_caffeinate
```

Monitoring:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
screen -ls
cat global_search_v5_nz5/v5_nz5_global_20260522_215756/status.json
cat global_search_v5_nz5/v5_nz5_global_20260522_215756/heartbeat.txt
tail -f global_search_v5_nz5/v5_nz5_global_20260522_215756/driver.log
sed -n '1,180p' global_search_v5_nz5/v5_nz5_global_20260522_215756/best_summary.md
```

Attach and detach:

```bash
screen -r v5_nz5_global_20260522_215756
```

Detach with `Ctrl-a d`.

Files written by the global search:

```text
config.json
bounds.json
evaluations.jsonl
evaluations.csv
latest.json
best.json
best_summary.md
status.json
heartbeat.txt
driver.log
driver.pid
driver.screen
caffeinate.screen
```

Search setup:

- `workers=8`
- `budget_sec=28800`
- `max_evals=1000000`
- `seed=20260523`
- `nb=30`
- `nz=5`
- `rho_z=0.95`
- `sigma_z=0.35`
- `kappa_entry=1000000`
- `max_iter_eq=35`
- `tol_eq=5e-4`
- `global_prob=0.75`
- Objective is full 19-moment weighted SMM loss plus penalties for GE failure
  or invalid outside-option scale accounting.

Search parameters:

```text
beta
b_entry_fixed
psi_child
h_bar_jump
h_bar_n
c_bar_n
kappa_fert
chi
kappa_loc
mu_move
theta0
theta_n
h_bar_0
E_C
r_bar_C
alpha_cons
phi
hR_max
h_own_max
```

Bounds are in:

```text
global_search_v5_nz5/v5_nz5_global_20260522_215756/bounds.json
```

How to summarize after the global search finishes:

1. Read `status.json`.
2. Confirm `"finished": true`.
3. Read `best_summary.md` and `best.json`.
4. Count completed evaluations:

```bash
wc -l global_search_v5_nz5/v5_nz5_global_20260522_215756/evaluations.jsonl
```

5. Inspect the top candidates:

```bash
python - <<'PY'
import csv, math
path = "global_search_v5_nz5/v5_nz5_global_20260522_215756/evaluations.csv"
rows = list(csv.DictReader(open(path)))
rows = [r for r in rows if r.get("objective_loss") not in ("", "inf", "Infinity")]
rows.sort(key=lambda r: float(r["objective_loss"]))
for r in rows[:20]:
    print(r["objective_loss"], r["moment_loss"], r["eval_id"],
          r["tfr"], r["childless_rate"], r["mean_age_first_birth"],
          r["tfr_gradient"], r["own_rate"], r["own_gradient"],
          r["young_liquid_wealth_to_income"], r["inv_pop_share_C"],
          r["inv_rent_ratio_C_over_P"])
PY
```

6. If a candidate beats `166.72`, rerun it cleanly at `Nz=5`. The global
   search script does not currently include a one-command replay helper, so
   the agent should construct a replay script from `best.json` by applying
   the `parameters` vector to the same V5 setup. Use `run_v5_nz5_global_search.py`
   functions if convenient.
7. If the clean `Nz=5` rerun is acceptable, validate with `Nz=7`.
8. Only then discuss moving Branch 1 into live code.

Important:

- Do not judge the search from a candidate with large penalty.
- A candidate with `accepted=False`, `strict_converged=False`, invalid scale,
  or negative `outside_entry_flow` should not be treated as viable even if
  moment loss looks low.
- The global search best objective includes penalties; sort by
  `objective_loss`, not just `moment_loss`.

## Branch 2: Developer Missing-Middle Supply

Folder:

```text
overnight_variants/2026-05-22_developer_missing_middle/
```

Status:

```text
Mechanically useful; economically yellow; do not merge.
```

Main script:

```text
run_developer_missing_middle_v3_ge.py
```

Main report:

```text
REPORT_V3_GE.md
```

What Branch 2 does:

- Full type-price developer GE prototype.
- Updates all six owner prices \(p_{iq}\) and rents \(r_{iq}\) for
  \(i \in \{C,P\}\) and \(q \in \{S,M,L\}\).
- Updates entry shares.
- Prices owner rungs by \(p_{iq(k)}\).
- Prices renters by realized room interval \(r_{iq(h^R)}h^R\).
- Fixed costs \(F_{iq}\) are reported through entry-threshold diagnostics.
  They are not active as discrete type shutdown conditions.

Why V3 mattered:

- Earlier renter block priced by \(\bar h_n\), mismatching the demand
  diagnostic by realized \(h^R\).
- V3 fixes this by solving renter choices piecewise over `S/M/L` intervals.
- This makes type-price GE mechanically more coherent.

Key branch 2 results from comparison:

| Moment | Target | Developer V3 GE |
|---|---:|---:|
| `tfr` | 1.700 | 1.899 |
| `childless_rate` | 0.150 | 0.148 |
| `mean_age_first_birth` | 26.000 | 33.698 |
| `tfr_gradient` | 0.133 | 0.110 |
| `own_rate` | 0.627 | 0.578 |
| `own_gradient` | 0.170 | 0.018 |
| `own_family_gap` | 0.110 | 0.070 |
| `prime_childless_renter_median_rooms` | 4.000 | 6.500 |
| `prime_childless_owner_median_rooms` | 6.000 | 8.200 |
| `housing_increment_0to1` | 0.664 | 0.418 |
| `housing_increment_1to2` | 0.566 | 0.252 |
| `young_liquid_wealth_to_income` | 0.600 | 0.815 |
| `center_share_nonparents` | 0.494 | 0.412 |
| `center_share_newparents` | 0.416 | 0.394 |
| `migration_rate` | 0.032 | 0.035 |
| `old_age_own_rate` | 0.863 | 0.767 |
| `old_age_parent_childless_gap` | 0.070 | -0.044 |
| `inv_pop_share_C` | 0.450 | 0.447 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.178 |

Read:

- Branch 2 clears mechanically.
- It does not solve the room target problem.
- Prime childless renters still choose too much housing.
- Prime childless owners are far above the room target.
- Ownership gradient is almost flat.
- Middle owner mass is too weak.

Decision:

- Do not merge Branch 2 into live code yet.
- If resumed, run a narrow diagnostic on room-bin/type mapping and owner rung
  ladder before broadening the solver.

Rerun:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_v3_ge.py --quiet --nb 30 --iterations 12 --price-damp 0.25 --entry-damp 0.25
```

## Cross-Branch Comparison

Current live benchmark closure in `CALIBRATION_STATUS.md` is
`renewal_valve_calibrated`. The outside-option closure is a paper-closure
prototype in the isolated V5 branch, not the current live benchmark closure.

High-level comparison:

- Branch 1 HANK-\(z\) adds structural income risk and the paper-facing
  outside-option benchmark normalization.
- Branch 2 missing-middle adds type-price developer supply and type-specific
  rents/prices.
- Branch 1 is the recommended next live direction, conditional on calibration.
- Branch 2 should remain a diagnostic/prototype until room-bin and owner-rung
  mechanics are better disciplined.

Do not compare losses across incompatible closures/grids without stating the
target system and grid. The current global search is `Nb=30,Nz=5`; final claims
require validation at `Nz=7`.

## Recent Commits Relevant to This Work

```text
74286b6 Normalize V5 outside closure inside benchmark loop
730e739 Add V5 HANK-z speed audit
1dd6178 Add V5 Nz5 equilibrium figure packet
5b10eef Add V5 Nz5 directional calibration audit
40eb8ff Add V5 Nz5 parallel global search
7866d3e Add overnight branch agent handoff
```

This offline handoff itself should be committed after creation.

## Suggested Next-Agent Work Plan

If the global search is still running:

1. Check `status.json`.
2. Check `screen -ls`.
3. Do not kill the run unless the user asks or it is clearly unhealthy.
4. If no heartbeat/status update appears for more than 30 minutes, investigate
   with `screen -r`, `driver.log`, and `ps`.

If the global search has finished:

1. Summarize best candidate and top 10-20 candidates.
2. Compare best to directional best (`166.72`).
3. Reject candidates with penalties or invalid scale.
4. Rerun best cleanly at `Nz=5`.
5. If clean rerun holds, run `Nz=7` validation.
6. Update `REPORT_COMPARISON.md`, branch README, and this handoff or a new
   post-search report.
7. Ask before committing large raw search outputs; a compact summary may be
   preferable.

If the user asks whether to merge:

```text
Do not merge yet. Branch 1 must first produce an acceptable calibrated candidate
and pass Nz=7 validation. Branch 2 is not merge-ready.
```

If the user asks what is most promising:

```text
Branch 1 with alpha_cons, kappa_fert, and phi in the search is promising.
The model responds. The hard remaining targets are first-birth timing,
childlessness, fertility gradient, and young liquid wealth.
```
