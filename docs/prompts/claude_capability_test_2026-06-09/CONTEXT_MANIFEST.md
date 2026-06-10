# Context Manifest

Claude should treat this repository as a research project, not as a generic
software package. The goal is to assess the theory draft, the active Python
model, the calibration system, and the project organization.

## Mandatory Startup

Read these first, in this order:

1. `memory/AGENT_MEMORY.md`
2. Latest `memory/daily/YYYY-MM-DD.md`, which is currently expected to be
   `memory/daily/2026-06-09.md`
3. `CALIBRATION_STATUS.md`
4. `SESSION_DIARY.md` only if historical detail is needed

Also read `CLAUDE.md` and `AGENTS.md`. If they conflict with this harness,
follow the repository-level instructions and report the conflict.

Before any substantive work, run:

```bash
git status -sb
```

Do not revert, delete, or overwrite user changes. The worktree may be dirty.

## Theory Artifacts

Primary current theory artifact:

- `latex/intergenerational_housing_fertility_v4.tex`
- `latex/intergenerational_housing_fertility_v4.pdf`

Candidate figures:

- `latex/figures/fig6_ce_planner_wedge.tex`
- `latex/figures/fig6_ce_planner_wedge.pdf`
- `latex/figures/fig7_entry_fertility_decomposition.tex`
- `latex/figures/fig7_entry_fertility_decomposition.pdf`

Historical comparison only:

- `latex/archive/intergenerational_housing_fertility_v3_20260609/`

Specific theory objects that must not disappear without explicit diagnosis:

- warm-glow bequests
- outside option and endogenous entry
- competitive equilibrium
- constrained planner problem
- one or two precise efficiency propositions
- distinction between fertility intensive margin and entry/extensive margin

## Active Codebase

Start from:

- `code/model/README.md`
- `code/model/dt_cp_model/parameters.py`
- `code/model/dt_cp_model/solver.py`
- `code/model/dt_cp_model/objective.py`
- `code/model/dt_cp_model/direct_calibration.py`
- `code/model/dt_cp_model/evaluate.py`
- `code/model/dt_cp_model/kernels.py`
- `code/model/dt_cp_model/theta.py`
- `code/model/tools/calibrate_direct_geometry.py`
- `code/model/tools/collect_direct_geometry_results.py`
- `code/cluster/submit_python_direct_geometry_overnight.sh`

Use `rg --files` and `rg` for discovery. Do not orient by broad scans of
generated output or archived folders. Avoid `calibration_archive/` except when
you need historical reference parity.

## Calibration And Output Context

Treat `CALIBRATION_STATUS.md` as the live calibration source of truth unless
local evidence shows it is stale. If it is stale, report the exact discrepancy.

High-value calibration checks:

- target values in `code/model/dt_cp_model/parameters.py`
- objective weights and active moment list in `code/model/dt_cp_model/objective.py`
- active direct-geometry parameterization in
  `code/model/dt_cp_model/direct_calibration.py`
- latest benchmark or best-record outputs named by `CALIBRATION_STATUS.md`
- output records under `output/model/`, only after identifying which run family
  is live and comparable

Do not compare losses across target systems, room scales, geographies, closure
rules, or objective definitions unless you have verified comparability.
