# Stationary-state sandbox

**Package-version notice (2026-09-16):** this sandbox currently imports
`intergen_eqscale_seq_optimized` from `main`. The retained September 13 2007
`corrected_initial` state was solved on a snapshot of branch
`codex/balanced-social-security` at commit `70abd4a8` plus a corrected
`solver.py` (sha256 `2992412586b81cef...`); `main`'s package differs from
that snapshot in 8 files (`kernels.py`, `solver.py`, `parameters.py`,
`utils.py`, plus three modules `main` lacks). A `psi_mode: joint` residual
trace confirms the economy is not at the retained equilibrium even at the
retained `(psi, transfer)`, consistent with a different solver rather than a
different root method. The snapshot is being fetched into
`output/model/e5f_final_night_20260913/corrected_initial_source_fetched/`.
A `package_root` spec/CLI option lets `run_ss.py` import that snapshot
instead of `main`'s package -- see below.

Change one assumption, re-solve the 2007 general-equilibrium stationary
state, and read moments-versus-targets and policy plots from one command.

## Commands

```
cd code/model
PYTHONPATH=. python sandbox/run_ss.py --spec baseline
PYTHONPATH=. python sandbox/run_ss.py --spec kappa_h_zero --fast
PYTHONPATH=. python sandbox/compare.py baseline kappa_h_zero
make sandbox SPEC=kappa_h_zero
make sandbox-check
make sandbox-test
```

`run_ss.py` writes exactly four files to `output/model/sandbox/<spec>/`
(or `--out DIR`): `summary.md`, `moments.csv`, `parameters.csv`, `graphs.pdf`.
Nothing else -- no receipts, hashes, JSON sidecars, or per-date folders.

## What it reproduces

The retained September 13 2007 stationary calibration:
`output/model/e5f_final_night_20260913/corrected_initial/`
(`parameters.csv`, `candidate_result.json`, `target_fit.csv`, 13 target
rows). `--spec baseline` loads that calibration's 9 free parameters and its
derived psi_child intercept, re-solves the same GE stationary state via the
production package's own solver (`run_e1_chain.py` / `solver.py`) and root
helpers (`audit_closed_reproductive_closure.py`,
`run_e5f_transition_calibration.py:solve_old_steady_state`), and should
reproduce it closely (see Deviations below for what does and does not match
to 1e-8).

`parameters.csv` and `candidate_result.json` report the 9 free coordinates
under **display names** (`beta_annual`, `h_P`) that differ from the literal
override keys the solver accepts (`beta`, `hbar_first_child_jump`).
`run_ss.py::display_theta_to_overrides` undoes that renaming; see its
docstring for the exact formulas and the two source citations that
established them.

## Spec format

`sandbox/specs/<name>.yaml` is a flat mapping (a tiny hand-rolled parser in
`spec_io.py` reads it; PyYAML is not installed in `code/model/.venv`, and
every spec here is a flat scalar mapping, so a real YAML parser was not worth
adding):

```yaml
# comment
overrides:
  tenure_choice_kappa: 0.0
fix_psi: false          # optional; true keeps psi_child fixed instead of
                          # re-deriving it to hit completed fertility 2.1
psi_child: 0.15          # optional; only read when fix_psi: true
```

Any key under `overrides:` is passed straight through as a `run_model_cp_dt`
override (parameter name -> value), except the three mechanism-switch keys
below, which are validated by `mechanisms.check_switches_supported` first.

Five specs are provided: `baseline`, `kappa_h_zero`
(`tenure_choice_kappa: 0.0`), `s1_concave_benefit`
(`child_benefit_form: log`), `s2_ssk_weighting` (`scale_weighting: multiply`),
and `s3_child_penalty` (`child_earnings_penalty: 0.10` -- **not runnable**,
see below).

## Mechanism switches

Implemented as runtime overrides in `sandbox/mechanisms.py`, monkeypatching
`solver.precompute_shared` only for the duration of one solve (restored
immediately after, even on exception). No file under
`code/model/intergen_eqscale_seq_optimized/` or `code/model/tools/` is ever
edited.

- **child_benefit_form** in `{linear, log, power}` with
  `child_benefit_curvature`: replaces `psi*m` (solver.py:2254) by
  `psi*log(1+m)` or `psi*((1+m)**(1-eps)-1)/(1-eps)`.
- **scale_weighting** in `{deflate, multiply}`: changes the equivalence-scale
  exponent at solver.py:2278 from `sigma-1` to `sigma`, for the
  `eqscale_form in {power, sqrt}` branches (the only branches that apply a
  sigma-dependent exponent in production). Under `eqscale_form == linear`
  (the `gamma_e` branch), `multiply` is a documented no-op: production
  applies no sigma exponent there at all, at any sigma, so there is no
  well-defined "multiply" analogue without inventing a new functional form.
- **child_earnings_penalty**: **not implemented**. `income_at_state`
  (solver.py:226) has no child-state argument, and `P.income`
  (parameters.py:725, `set_income_given_w_and_pension`) is shaped `(I, J)`
  with no child-count dimension; the Bellman-loop call sites that do have
  the child state in scope call `income_at_state` without it. Applying
  `(1 - tau_c(m))` from sandbox code alone would need either a wider
  `income_at_state` signature threaded through ~10 call sites, or an
  `(I, J, n_child_states)` income array -- both are package edits this
  sandbox may not make. Running `--spec s3_child_penalty` raises
  `NotImplementedError` with this explanation rather than silently no-op'ing
  or producing wrong numbers.

Each implemented switch has a bitwise-off unit test in
`sandbox/tests/test_mechanisms.py`: at the default value
(`child_benefit_form=linear`, `scale_weighting=deflate`),
`sandboxed_precompute_shared` returns the untouched package function's
output with no rebuild at all, so it is bitwise-identical to unpatched
production by construction, checked against a real (tiny-grid) parameter
object built the same way `run_ss.py` builds one.

## What the four output files contain

- **summary.md**: spec name and overrides, wall time per Bellman/root
  evaluation and total, number of root evaluations, warm-start note, the
  scalar loss, the 13-row target table (target/model/gap/weight/loss
  contribution, with a `source` column -- see Deviations), and the parameter
  table (estimate/bounds/near-bound flag).
- **moments.csv**: one row per target-table moment (`moment`, `model`,
  `source`).
- **parameters.csv**: one row per free/fixed parameter plus a
  `_solved_price` bookkeeping row (the equilibrium price vector, used only
  to warm-start a later spec's price root -- not a model parameter).
- **graphs.pdf**: the package's own standard 17-graph diagnostic packet
  (`intergen_eqscale_seq_optimized.diagnostics.write_diagnostics`, unchanged
  and unredesigned) plus a supplemental renter-housing-policy-by-age page,
  as one multi-page PDF.

## Warm start

The price root is warm-started from a prior baseline run's saved price
(`p_init_override`, an existing hook in `run_model_cp_dt`, solver.py:1098)
when `--out` is not `baseline` and a baseline run exists; `summary.md`
reports how many evaluations that saved. **Bellman/value-function warm start
is not implemented**: `run_model_cp_dt` has no V-init override hook
analogous to `p_init_override` (the value function is always cold-started
inside the package), so this half of the requested warm start cannot be done
from sandbox code without a package edit; `summary.md` states this
explicitly on every run.

## Regression gate and measured timing

`make sandbox-check` runs `--spec baseline` at the full production grid
(Nb=120, J=17) and asserts the sandbox's model moments match the retained
`target_fit.csv` to 1e-8, and that the pension/rebate residuals meet the
solver's own 2e-4 gate. See "Deviations from the literal spec" below for
which of the 13 rows this check actually covers and why.

Measured locally (see the assistant's final report for the actual run's
numbers): expect roughly 5-8 minutes per full-grid Bellman/root evaluation
(the retained calibration's own two verification repetitions took 752s/2 =
376s each), and 1 evaluation for `--spec baseline` (the retained psi is
already within the 1e-6 fertility tolerance, so the psi-root search
short-circuits at the first evaluation). `--fast` uses Nb=40 instead of 120.

## Deviations from the literal spec (read before trusting a number)

1. **The 13-row target table is only exactly reproducible for a subset of
   rows.** The retained `target_fit.csv`'s 12 scored rows were computed by a
   specialized dated period/cohort-timing measurement pipeline
   (`transition_cross_section_moments` / `cohort_timing_moments` in
   `run_e5f_transition_calibration.py`), not by the plain stationary
   `chain.extract_moments(sol, P)` this sandbox reuses. `sandbox/target_table.py`
   maps each retained label to the closest `extract_moments` key and tags it
   `exact_stationary_key` (5 rows, including the completed-fertility
   normalization: matches or should match very closely),
   `approximate_stationary_analogue` (7 rows: e.g. mean occupied rooms
   uncapped rather than capped at 9, since `extract_moments` has no
   capped-rooms measure), or `unavailable_in_extract_moments` (1 row,
   "exactly one child among mothers 40-44": no analogue exists at all).
   `make sandbox-check` gates only the `exact_stationary_key` rows to 1e-8;
   the rest are reported for inspection, not gated, and `summary.md` states
   this caveat on every run.
2. **Bellman/value-function warm start is not implemented** (see Warm start
   above) -- no such hook exists in the package.
3. **`child_earnings_penalty` is not implemented** (see Mechanism switches
   above) -- no clean sandbox-only hook exists; `--spec s3_child_penalty`
   fails loudly rather than silently.
4. **`scale_weighting: multiply` only changes anything under
   `eqscale_form in {power, sqrt}`.** The retained calibration uses
   `eqscale_form=power` (via the E5 profile), so `s2_ssk_weighting.yaml`
   does have an effect; a spec that also sets `eqscale_form: linear` would
   see no effect from `scale_weighting: multiply`, by construction (see
   point 2 under Mechanism switches).
