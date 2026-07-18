# Codex worker task: means-tested transfer floor (probe-only, default-off) + probe driver

## Goal

Two deliverables in `code/model/`:

(A) A minimal, **default-off** means-tested transfer floor in the
intergenerational solver, markov-income Bellman path only.

(B) A local probe driver `code/model/tools/run_intergen_transfer_floor_probe.py`
that runs 8 fixed-theta GE solves and appends per-cell results to disk.

Nothing may change any result when the floor parameters are zero (the default).
Bitwise-identical behavior with the floor off is an acceptance criterion.

## Scope — exclusive file ownership

- `code/model/intergen_housing_fertility/parameters.py` (2-line addition)
- `code/model/intergen_housing_fertility/solver.py` (precompute_shared, markov
  Bellman call sites + fallback, core-path guard, census)
- `code/model/intergen_housing_fertility/kernels.py` (the two `full_*_block_kernel`s only)
- `code/model/tests/test_transfer_floor.py` (new)
- `code/model/tools/run_intergen_transfer_floor_probe.py` (new)

## Do not touch

- `tenure_choice_kernel` / `tenure_logit_kernel` (transfer must NOT fund down payments)
- the eval kernels `eval_renter_block_kernel` / `eval_owner_block_kernel` (core-path only; guarded off instead)
- KFE mass-moving code, `_censor_entry_dead_mass`, bequest utilities (`Vbq`), `income_at_state`
- `calibration.py`, `local_panel.py`, `production_profile.py`, `cli.py`
- anything under `output/`
- the debt/rollover floor logic inside the kernels (lines using `b_grid[b]` / `current_b` / `current_unsecured`) — the transfer must not expand borrowing capacity

## Economic specification (normative)

Per-period means-tested transfer, DNFJ-style (De Nardi–French–Jones 2010 eq. 10
shape), benefit-capped:

    Rv      = R_gross * b + y                (period cash-on-hand, existing)
    G(n,s)  = G0 + Gn * nk   if children present (same kp logic as c_bar)
            = G0             otherwise
    T       = min( max(0, G(n,s) - Rv), G(n,s) )     # cap: deep-debt states get at most G — no debt bailout
    Rv_eff  = Rv + T

Units: **period (4-year) flows**, same convention as `c_bar_0` (`P.c_bar_0 =
0.10 * P.period_years` style). G0, Gn default to 0.0 → T ≡ 0 → all arithmetic
reduces to today's exactly.

The floor must not: expand debt capacity, count toward down payments, enter
estates, or change measured gross income (`income_at_state` untouched).

## A. Solver changes, site by site

### A1. parameters.py — in `setup_parameters()`, immediately after the
`P.c_bar_0 ... P.h_bar_n` block (~lines 49–53):

```python
# Means-tested floor guarantee (period units, like c_bar_0); 0 = no floor.
P.transfer_floor_G0 = 0.0
P.transfer_floor_Gn = 0.0
```

### A2. solver.py — `precompute_shared()` (lines 2010–2067)

Read `g0 = float(getattr(P, "transfer_floor_G0", 0.0))`, `gn = float(getattr(P,
"transfer_floor_Gn", 0.0))`. Build `g_bar = np.zeros((P.n_parity,
P.n_child_states))` inside the existing nn/cs loop, mirroring `c_bar`:

- if `kp`: `g_bar[nn, cs] = g0 + gn * nk`
- else: `g_bar[nn, cs] = g0`

Add to the returned SimpleNamespace: `g_bar=g_bar`,
`gb_flat=g_bar.reshape(1, nc, order="F")`.

### A3. kernels.py — `full_renter_block_kernel` (line 753)

- New argument `gb_v` **immediately after `psi_v`** (line 761), comment `# (nc,)`.
- In the column loop, hoist after `psic = psi_v[c]` (line 794): `gc = gb_v[c]`.
- Immediately after `Rvb = Rv1d[b]` (line 802) insert:

```python
if gc > 0.0:
    Tb = gc - Rvb
    if Tb > 0.0:
        if Tb > gc:
            Tb = gc
        Rvb = Rvb + Tb
```

Everything downstream (hi at 810, eval_renter_scalar calls, surplus at 862, ct
at 869) reads the topped-up `Rvb` automatically. Do NOT touch `current_b`
/ `rollover_floor` (804–806).

### A4. kernels.py — `full_owner_block_kernel` (line 886)

Same: new arg `gb_v` after `psi_v` (line 894); hoist `gc = gb_v[c]` after
`psic = psi_v[c]` (line 921); same 6-line block right after `Rvb = Rv1d[b]`
(line 929). Do NOT touch `current_unsecured` (932) — it must stay on
`b_grid[b] - bf`.

### A5. solver.py — markov Bellman `solve_bellman_full_markov_income` (2080–2437)

- Where `cb_v` / `hb_v` / `psi_v_flat` are prepared as contiguous 1-d `(nc,)`
  arrays for the kernels, build `gb_v` identically from `SD.gb_flat`.
- Pass `gb_v` in the `full_renter_block_kernel` call (2216–2221) and the
  `full_owner_block_kernel` call (2273–2278).
- Non-numba fallback branch: right after `d_nc = SD.cb_flat + ri * SD.hb_flat`
  (2224) compute once:

```python
Rv_eff_nc = Rv + np.clip(SD.gb_flat - Rv, 0.0, SD.gb_flat)   # (Nb, nc)
```

  (`np.clip(x, 0.0, 0.0) == 0.0`, so with the floor off this equals `Rv`
  broadcast — values identical.) Then in the fallback only, replace:
  - `hi` (2235): `Rv[:, 0]` → `Rv_eff_nc[:, c]`
  - `golden_renter` Rv argument (2237): `Rv[:, 0]` → `Rv_eff_nc[:, c]`
  - `surplus_nc` (2243): `Rv` → `Rv_eff_nc`
  - `ct_cap` (2248): `Rv` → `Rv_eff_nc`
  - owner fallback: `hi` (2290) and `golden_owner` Rv argument (2292):
    `Rv[:, 0]` → `Rv_eff_nc[:, c]`; `co_nc` (2297): `Rv` → `Rv_eff_nc`

  Note the owner fallback sits in a different loop than the renter fallback;
  `Rv_eff_nc` can be computed once per (i, z) block right after `Rv` (2209) if
  that is cleaner — it must be identical for both tenure branches.

### A6. solver.py — `solve_bellman_core` (line 2438), first statements of the body:

```python
if float(getattr(P, "transfer_floor_G0", 0.0)) != 0.0 or float(getattr(P, "transfer_floor_Gn", 0.0)) != 0.0:
    raise NotImplementedError("transfer floor: markov-income Bellman path only")
```

The core-path calls of the two full kernels (2602, 2674) must pass a zeros
array for `gb_v` (build `gb_zero = np.zeros(nc)` once near the top). Core
behavior stays byte-identical.

### A7. solver.py — `_dead_mass_census_at_age` (3382–3450), after `resources = ...` (3416):

```python
gG = float(SD.g_bar[nn, cs]) if hasattr(SD, "g_bar") else 0.0
tr = min(max(gG - resources, 0.0), gG)
resources = resources + tr
```

and add `"transfer": tr` to the census dict. (Descriptive only — deadness
detection is V-based at 3404; this keeps the census truthful under the floor.)

## B. Probe driver `code/model/tools/run_intergen_transfer_floor_probe.py`

Mirror the M5 tight-solve contract from `tools/run_intergen_bequest_exit_chain.py`
EXACTLY (this is how gate0 reproduced the Torch winner locally):

- Build an `argparse.Namespace` with: `arm="M5"`, `J=PRODUCTION_J`,
  `Nb=PRODUCTION_SEARCH_NB`, `max_iter_eq=40`, `tol_eq=2.5e-5`,
  `ltv_terminal=0.4`, `theta1=0.25`, `seed_theta0=0.30`, `seed_theta_n=0.75`,
  `seed_kappa=0.0`, `fixed_theta0=None` (add any other fields
  `arm_contract`/`survival_schedule`/`common_overrides` read).
- `active, fixed, mechanism = arm_contract(args)`;
  `base = common_overrides(args, mechanism)`;
  `tight = {**base, "max_iter_eq": 40, "tol_eq": 2.5e-5}`.
- `theta = load_theta(ROOT / "output/model/intergen_income_disciplined_recalibration_20260716/report/results.json", seed_arm="M5")`.
- Import from `tools.run_intergen_bequest_exit_chain`: `load_theta`,
  `target_system`, `arm_contract`, `common_overrides`, `survival_schedule` (as
  needed). Import from the package: `run_model_cp_dt`, `InfeasibleThetaError`,
  `income_at_state`, `income_transition_values`, `extract_moments`,
  `diagnostic_loss`, `get_target_set`, `income_process_overrides`.
- Per cell: `cell = {**tight}`; if the cell changes sigma:
  `cell.update(income_process_overrides(5, "rouwenhorst", sigma_annual, 0.9601845894041878))`.
  Then `cell["transfer_floor_G0"] = G0_period`, `cell["transfer_floor_Gn"] = Gn_period`,
  and `theta_cell = {**theta}` with `c_bar_0` replaced where the cell says so.
  Solve `sol, P, p_eq = run_model_cp_dt({**cell, **theta_cell}, verbose=False)`.
- CELLS (names use ANNUAL units; pass PERIOD units = annual × 4). Annual→period:
  c_bar_0 0.10→0.40, 0.14→0.56; G0 0.13→0.52, 0.16→0.64; Gn 0.10→0.40.

| # | name | sigma_annual | c_bar_0 (period) | G0 (period) | Gn (period) | expectation |
|---|------|--------------|------------------|-------------|-------------|-------------|
| 0 | base_repro | matched (leave `base` income overrides untouched — do NOT re-apply income_process_overrides) | theta's own (≈1.2596) | 0.0 | 0.0 | must reproduce gate0 |
| 1 | floor_only | matched (untouched) | theta's own | 0.52 | 0.40 | ≈ baseline; isolates floor |
| 2 | lowc_nofloor | 0.20 | 0.40 | 0.0 | 0.0 | does low c_bar_0 alone restore feasibility? |
| 3 | s12_lowc | 0.12 | 0.40 | 0.52 | 0.40 | main grid |
| 4 | s20_lowc | 0.20 | 0.40 | 0.52 | 0.40 | main grid |
| 5 | s12_midc | 0.12 | 0.56 | 0.64 | 0.40 | main grid |
| 6 | s20_midc | 0.20 | 0.56 | 0.64 | 0.40 | main grid |
| 7 | gate_check_smallG | 0.20 | theta's own (≈1.2596) | 0.52 | 0.40 | MUST raise InfeasibleThetaError (floor below bundle must not mask infeasibility) |

- Wrap each solve in `try/except InfeasibleThetaError as exc` → record
  `{"status": "infeasible", "stage": exc.stage, "dead_mass": exc.dead_mass,
  "census_head": exc.census[:5]}` and continue to the next cell.
- Per-cell record (JSON): cell params (annual and period), status, elapsed
  seconds, `market_residual = float(sol.best_max_abs_rel_excess)`, the FULL
  15-row target-fit table (moment, target, model, gap, weight,
  loss_contribution — replicate the exit-chain `target_fit()` logic with
  `target_system()`), `diagnostic_loss`, the complete `extract_moments(sol, P)`
  dict, `entry_censored_share = float(getattr(sol, "entry_censored_share", -1.0))`,
  and the receipt statistics below.
- Receipt statistics (post-processing inside the probe script, NO solver
  changes): rebuild `g_bar` from P with the same formula; get `z_grid` from
  `income_transition_values(P)`; with `g = np.asarray(sol.g)` of shape
  `(Nb, nt, I, J, Nz, npar, ncs)` and `b_grid = np.asarray(sol.b_grid)`,
  compute vectorized (broadcast/einsum, no per-node Python loops):
  `Rv[b,ten,i,j,z,nn,cs] = P.R_gross * b + income_at_state(P, i, j, z)`
  (income_at_state depends on (i, j, z) only — precompute a (I, J, Nz) array),
  `T = clip(g_bar[nn,cs] - Rv, 0, g_bar[nn,cs])`. Report:
  `receipt_mass_share` (mass with T>0 / total mass), `outlays_total`
  (sum g*T), `outlays_over_income` (sum g*T / sum g*y), receipt splits: renter
  vs owner, b<0 vs b>=0, age bands 18–33 / 34–65 / 66+.
- After cell 0: load
  `output/model/intergen_m5_draft_refresh_20260717/gate0/gate0_result.json`,
  compare loss and market_residual at full precision, print MATCH/MISMATCH
  with both values; continue either way (the lead adjudicates).
- Output dir `output/model/transfer_floor_probe_20260718/`: append each cell to
  `results.jsonl` (flush per cell), write final `results.json`, print one
  progress line per cell (name, status, loss, elapsed).
- `--cells` CLI arg (comma list of cell indices, default all), `--outdir`
  override. `if __name__ == "__main__": main()`.
- Script must run with `PYTHONPATH=code/model` from the repo root (absolute
  imports `from intergen_housing_fertility...` / `from tools...` like the other
  tools scripts).

## C. New test `code/model/tests/test_transfer_floor.py` (fast, no GE solve)

1. `precompute_shared` g_bar: tiny P via `setup_parameters()` + overrides;
   check g_bar for a kp column (= G0 + Gn*nk) and a non-kp column (= G0), and
   that gb_flat matches the F-order reshape.
2. Kernel equivalence off: call `full_renter_block_kernel` on a tiny grid
   (Nb≈5, nc=1) twice — once with the new `gb_v = np.zeros(1)`, once comparing
   against expected hand-built values — and assert the zeros call reproduces
   the SAME outputs as a pre-change reference computed by the same formula
   with no floor (use `np.array_equal` on Vo, bp, co, ho for gb_v=0 vs a
   duplicate call; and hand-check one interior cell).
3. Floor revives a dead node: construct a cell with `Rv < dc` (dead: Vo ≤ -1e9)
   under `gb_v=0`, then with `gb_v = [G]`, `G > dc + margin`, assert
   `Vo > -1e9` and `co ≥ c_bar` at that node.
4. Benefit cap: a deep-debt node (`Rv < 0`) with `G > 0`: assert the implied
   top-up never exceeds G — e.g. a node with `Rv + G < dc` must STAY dead.
5. Owner kernel: repeat (2) and (3) minimally for `full_owner_block_kernel`.

## Verification (run these; paste output verbatim)

- `python3 -m py_compile` on the three edited files + two new files.
- `PYTHONPATH=code/model python3 -m pytest code/model/tests/test_transfer_floor.py -q`
- `PYTHONPATH=code/model python3 -m pytest code/model/tests/test_entrant_feasibility.py -q` (must stay green)
- Import check only for the driver:
  `PYTHONPATH=code/model python3 -c "import tools.run_intergen_transfer_floor_probe"`.
- Do NOT run any GE solve or the probe cells — the lead runs those.

## Stop and report if

- Any line anchor above does not match the current code (drift) — stop, report
  the mismatch, change nothing else.
- Any previously green test goes red.
- The change would need to touch a file in the Do-not-touch list.
- Ambiguity about which Rv site feeds which kernel.

## Required output

Files changed with line ranges, the verbatim test/verification output, and any
deviation from this spec (there should be none without a stop-and-report).
