# Addendum v2: debt-blind means test (fixes the leveraged-owner eligibility artifact)

Round-1 probe finding: means-testing on cash-on-hand `Rv = R_gross*b + y`
makes every mortgage-leveraged owner (b<0 secured) eligible — 79% receipt,
outlays 18.9% of income, floor becomes a mortgage subsidy. Fix: the means
test must use **debt-blind resources**

    x_test = y + R_gross * max(b, 0)        # income + positive liquid assets; debt never qualifies you
    T      = min( max(0, G(n,s) - x_test), G(n,s) )
    Rv_eff = Rv + T                          # budget still pays debt service out of Rv as before

`x_test >= 0` always, so the cap can bind only via G itself; keep it as a belt.
Everything else from tmp/transfer_floor_probe_spec_20260718.md stands
(default-off, markov-only, no borrowing-capacity/down-payment/estate/income-
measurement contamination).

## Changes

### kernels.py — both `full_renter_block_kernel` and `full_owner_block_kernel`

- New argument `Rvt1d` immediately AFTER `Rv1d` (comment `# (Nb,) debt-blind means-test resources`).
- Replace the round-1 insert (after `Rvb = Rv1d[b]`) with:

```python
if gc > 0.0:
    Tb = gc - Rvt1d[b]
    if Tb > 0.0:
        if Tb > gc:
            Tb = gc
        Rvb = Rvb + Tb
```

### solver.py — markov path `solve_bellman_full_markov_income`

- After `Rv = Rg * b + yj` (2209) add:

```python
Rv_test = Rg * np.maximum(b, 0.0) + yj
```

  and next to `Rv1d_full = np.ascontiguousarray(Rv[:, 0])` add
  `Rvt1d_full = np.ascontiguousarray(Rv_test[:, 0])`.
- Pass `Rvt1d_full` right after `Rv1d_full` in both full-kernel calls.
- Fallback: replace the `Rv_eff_nc` line with

```python
Rv_eff_nc = Rv + np.clip(SD.gb_flat - Rv_test, 0.0, SD.gb_flat)
```

  (downstream fallback uses of `Rv_eff_nc` unchanged).

### solver.py — core path

Guard already blocks the floor there; pass `Rv1d_full` again for the new
`Rvt1d` slot in the two core-path kernel calls (inert because `gb_zero`).

### solver.py — census `_dead_mass_census_at_age`

Replace the round-1 transfer lines with:

```python
gG = float(SD.g_bar[nn, cs]) if hasattr(SD, "g_bar") else 0.0
x_test = float(P.R_gross) * max(b_now, 0.0) + y_now
tr = min(max(gG - x_test, 0.0), gG)
resources = resources + tr
```

("transfer": tr stays in the census dict.)

### tools/run_intergen_transfer_floor_probe.py

- `receipt_statistics`: compute `test_resources = R_gross * np.maximum(b_grid, 0)[:, None, ...] + income[...]`
  and `transfer = np.clip(guarantee - test_resources, 0.0, guarantee)`
  (replace the old `resources`-based transfer; keep the rest).
- CELLS: change cells 5 and 6 `G0_period` from 0.64 to **0.72** (0.18 annual —
  top of the measured childless range incl. housing assistance; round-1
  showed the actual bundle cost c_bar_0 + r*h_bar_0 at the solved rent is
  ~0.165 annual, so 0.16 was knife-edge and died).
- Update cell expectations: cell 1 -> "receipt confined to bottom-income
  states and poor retirees; moments near baseline"; cells 5/6 -> "main grid".

### tests

- `test_transfer_floor.py`: update direct kernel calls for the new `Rvt1d`
  argument. For b>=0 nodes pass `Rvt1d = resources` (same vector) — there
  max(b,0)=b so test resources equal cash-on-hand when the test constructs
  resources as R*b+y with b>=0. Replace the round-1 deep-debt cap test with
  the sharper regression this bug demands:
  - **mortgage-blind test**: a node with resources < dc (deep debt) but
    `Rvt1d` (income-side) ABOVE G gets NO transfer: outputs identical to the
    gb=0 call (np.array_equal) — leverage is not subsidized;
  - **income-poor test**: a node with `Rvt1d < G` and G > dc is revived
    (V > -1e9) even at negative b, with the top-up never exceeding G
    (assert a node with Rv + G < dc stays dead).
- `test_entrant_feasibility.py`: the two direct kernel calls get the new
  `Rvt1d` argument right after the resources vector — pass the same
  `resources` array (inert: gb_v is zeros there).

## Verification (same discipline as round 1)

- py_compile all touched files.
- `PYTHONPATH=code/model code/model/.venv/bin/python3 -m pytest code/model/intergen_housing_fertility/tests/test_transfer_floor.py code/model/intergen_housing_fertility/tests/test_entrant_feasibility.py -q`
  (use this exact interpreter; it has numba+pytest).
- Import check the driver. No GE solves.

## Stop and report if

Anchors drifted (they will have moved by the round-1 diff — locate by the
round-1 inserted lines), any green test goes red, or any Do-not-touch file
would need changes.
