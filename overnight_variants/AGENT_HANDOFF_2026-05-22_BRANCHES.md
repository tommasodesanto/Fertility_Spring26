# Agent Handoff: 2026-05-22 Overnight Branches

This is the short handoff for the two isolated overnight branches. Work stayed
inside `overnight_variants/`; do not edit or revert live code under
`code/model/dt_cp_model/`.

## Read First

- Root instructions: `AGENTS.md`
- Live calibration status: `CALIBRATION_STATUS.md`
- Variant comparison: `overnight_variants/REPORT_COMPARISON.md`
- Branch 1 folder:
  `overnight_variants/2026-05-22_income_mortgage_risk/`
- Branch 2 folder:
  `overnight_variants/2026-05-22_developer_missing_middle/`

There are pre-existing dirty files outside `overnight_variants/`, including
live model docs/code and LaTeX files. They were not touched by this work and
should not be staged or reverted unless the user explicitly asks.

## Branch 1: HANK-z Income Risk / Outside Closure

Folder:
`overnight_variants/2026-05-22_income_mortgage_risk/`

Status: **mechanically working, economically yellow, currently calibrating**

What was implemented:

- A true structural idiosyncratic earnings state \(z\) in the copied solver,
  with state \((b,d,i,a,n,s,z)\).
- Rouwenhorst income process. Serious validation grid so far:
  \(N_z=7\), \(\rho_z=0.95\), unconditional \(\sigma_z=0.35\).
- No structural mortgage-account state \(\mu\). Existing mortgage/default
  objects in this branch are diagnostic only.
- Paper-facing outside-option closure in V5.

Main closure decision:

- The benchmark should impose \(S=1\) directly inside the GE loop.
- At each candidate benchmark composition, compute \(E_0(p)\) and \(B_0(p)\),
  calibrate \(\bar W^E\) to target \(q^E\), and set
  \[
  M = E_0(p)/q^E(p)-B_0(p).
  \]
- The earlier outer normalization/re-solve pass was removed as unnecessarily
  complicated for the benchmark. The final \(M\) and \(\bar W^E\) are the
  counterfactual objects.

Key V5 `Nz=7` validation run:

- Command report: `REPORT_V5_HANK_Z_OUTSIDE_CLOSURE.md`
- Accepted strict GE in 17 iterations.
- Runtime: `620.13s` cold.
- Final GE error: `0.00022218249381993278`
- \(S=1.0000000000000009\)
- \(q^E=0.9000000000000001\)
- outside probability: `0.1`
- \(M=0.005780247614275978\)
- outside value: `-3327810.8839415493`
- Loss: `311.4039737325341`

Runtime decision:

- `Nz=7` is validation-only for now.
- `Nz=5` is the exploratory/calibration grid.
- Same-process warm-up in the speed audit meant Numba/shape warm-up only, not
  an economic warm start.

Key V5 `Nz=5` plotting run:

- Report: `REPORT_V5_HANK_Z_OUTSIDE_CLOSURE_NZ5.md`
- Figure packet:
  `figures_v5_hank_z_outside_closure_nz5/HANK_Z_OUTSIDE_CLOSURE_NZ5_FIGURE_PACKET.pdf`
- Accepted strict GE in 12 iterations.
- Runtime including plot packet: `83.10s`
- Final GE error: `0.0003464318147695515`
- \(S=1\), \(q^E=0.9\), outside probability `0.1`
- \(M=0.005806434956255813\)
- outside value: `-1504137.9225856075`
- Loss: `311.29977174646774`

Directional calibration audit:

- Report: `REPORT_V5_NZ5_DIRECTIONAL_CALIBRATION.md`
- Script: `run_v5_nz5_directional_calibration.py`
- Best probe so far:
  `alpha_high_fertility_mid_finance`, with
  `alpha_cons=0.80`, `kappa_fert=4.0`, and `phi=0.90`.
- Accepted strict GE in 16 iterations.
- Runtime: `107.73s`
- Loss: `311.30 -> 166.72`
- Improved 16 of 19 target gaps.
- Good movements: TFR near target, ownership level near target, ownership
  gradient becomes positive, room moments improve, geography improves modestly.
- Remaining failures: first-birth age far too late, childlessness too high,
  fertility gradient wrong-signed, young liquid wealth too high/worse.

Current global search:

- Script: `run_v5_nz5_global_search.py`
- Plan: `REPORT_V5_NZ5_GLOBAL_SEARCH_PLAN.md`
- Objective: full 19-moment weighted SMM loss plus GE/scale penalties.
- Search dimension: 19 parameters, including original structural calibration
  block, `E_C`, `r_bar_C`, `alpha_cons`, `phi`, `hR_max`, and `h_own_max`.
- Workers: 8 parallel worker processes, one numerical thread each.
- Run tag:
  `v5_nz5_global_20260522_215756`
- Run directory:
  `overnight_variants/2026-05-22_income_mortgage_risk/global_search_v5_nz5/v5_nz5_global_20260522_215756/`
- Screen sessions:
  `v5_nz5_global_20260522_215756` and
  `v5_nz5_global_20260522_215756_caffeinate`

Monitoring commands:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
cat global_search_v5_nz5/v5_nz5_global_20260522_215756/status.json
tail -f global_search_v5_nz5/v5_nz5_global_20260522_215756/driver.log
sed -n '1,140p' global_search_v5_nz5/v5_nz5_global_20260522_215756/best_summary.md
screen -ls
```

Snapshot at handoff:

- completed: `16`
- active: `8`
- submitted: `24`
- best objective: `166.7224914373928`
- best eval id: `2`
- best is still the seeded `alpha_cons=0.80`, `kappa_fert=4.0`, `phi=0.90`
  case; random/global candidates had not yet beaten it.

Expected next action:

1. Let the 8-hour global search run.
2. After completion, inspect `best_summary.md`, `best.json`, and
   `evaluations.csv`.
3. If a candidate materially beats `166.72`, rerun it cleanly at `Nz=5`, then
   validate at `Nz=7`.
4. Do not move this branch into live code until there is an acceptable
   calibrated candidate and an `Nz=7` validation run.

## Branch 2: Developer Missing-Middle Supply

Folder:
`overnight_variants/2026-05-22_developer_missing_middle/`

Status: **mechanically useful, economically yellow, do not merge**

What was implemented:

- Full type-price developer GE prototype in
  `run_developer_missing_middle_v3_ge.py`.
- Updates all six \(p_{iq}\) and \(r_{iq}\) objects for
  \(i\in\{C,P\}\) and \(q\in\{S,M,L\}\), plus entry shares.
- Owner rungs are priced by \(p_{iq(k)}\).
- Renter choices now use realized room interval pricing
  \(r_{iq(h^R)}h^R\), fixing the earlier mismatch that priced renters by
  \(\bar h_n\).
- Fixed costs \(F_{iq}\) are diagnostic through entry thresholds only; they
  are not active as discrete shutdown conditions.

Key run:

- Report: `REPORT_V3_GE.md`
- Script:
  `run_developer_missing_middle_v3_ge.py`
- Accepted full type-price GE on the coarse `Nb=30` run.
- Runtime category: `expensive`.
- Moment comparison is in `overnight_variants/REPORT_COMPARISON.md`.

Branch 2 read:

- The corrected full-GE prototype no longer fails mechanically.
- Middle demand is non-degenerate.
- It still does not solve the room target problem:
  prime childless renter median rooms are `6.500` versus target `4.000`,
  prime childless owner median rooms are `8.200` versus target `6.000`, and
  middle owner mass remains essentially too weak.
- The ownership gradient is nearly flat (`0.018` versus target `0.170`).

Branch 2 decision:

- Do not merge into live code yet.
- Next step, if resumed, should be a narrow diagnostic on room-bin/type mapping
  and the owner rung ladder, not another broad solver merge.

Rerun command:

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_v3_ge.py --quiet --nb 30 --iterations 12 --price-damp 0.25 --entry-damp 0.25
```

## Cross-Branch Decision

Recommendation remains: **Branch 1 first, but only as a controlled live branch
after calibration evidence improves.**

Reasons:

- Branch 1 is the cleaner canonical income-risk extension: one structural
  \(z\) state, standard Rouwenhorst process, and a paper-consistent
  outside-option benchmark normalization.
- Branch 1 now clears strict GE and responds to calibration parameters.
- Branch 2 is useful but has not fixed the room/ownership-gradient problem.
- Do not add structural \(\mu\) yet. The current evidence says the next problem
  is recalibrating HANK-\(z\), not expanding the state space.

## Recent Commits

- `74286b6` Normalize V5 outside closure inside benchmark loop
- `730e739` Add V5 HANK-z speed audit
- `1dd6178` Add V5 Nz5 equilibrium figure packet
- `5b10eef` Add V5 Nz5 directional calibration audit
- `40eb8ff` Add V5 Nz5 parallel global search

The live global-search output directory is untracked while the run is active.
Do not commit partial global-search outputs until the run finishes or the user
asks for an intermediate checkpoint commit.
