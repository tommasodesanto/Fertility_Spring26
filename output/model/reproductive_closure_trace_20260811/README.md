# Fixed-price trace of the reproduction and clearing schedules (smoke grade)

Date: 2026-08-11 (overnight session). Companion to
`latex/reproductive_closure_theory_note.tex`, which defines every object used
here. Diagnostic only — not a calibration, not a policy run.

## What this is

The reproductive closure needs two schedules: mature entrant households per
incumbent household, `B(r)`, against the replacement requirement `E`; and the
population supported by housing clearing, `S(r) = H^S(r)/D(r)`. Both are
computable from fixed-price solves of the stationary model (one Bellman pass
plus one distribution pass per price; no equilibrium loop). This packet traces
them on 37 log-spaced asset prices from 0.02 to 2.40 (benchmark 0.7761).

## Provenance

- Source parameterization: `output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json`,
  `winners.E1.theta` (child-room-floor arm; exact strict repeat verified by the loader).
- Policy held fixed across all prices at the funded baseline: 1% annual
  property tax, lump-sum transfer `0.2091064453125` (its balanced value at the
  benchmark price only).
- Solver: `intergen_eqscale_seq_optimized` package via
  `run_e5_repaired_policy_with_entry.load_repaired_profile(smoke=True)` and
  `solve_mode="pe"` with `p_fixed`; smoke grid `Nb=40`, full lifecycle `J=17`
  (required for the demographic accounting).
- Reproduce: `code/model/.venv/bin/python trace_reproduction_schedule.py trace_points.csv <prices...>`
  then `build_trace_figure.py`. Runtime ~1.4s per price point locally.

## Headline readings (smoke grade)

- `E = 0.061733456` at every point and every price — the replacement
  requirement is a demographic constant (theory note, Lemma 1).
- `B(r)/E` is monotonically decreasing in the price everywhere on the grid
  (theory note, Assumption in Prop. 1 verified) and nearly flat: benchmark
  `0.8573` (production-grade value `0.85738`), rising only to a ceiling of
  `0.8764` as the price falls to 2.6% of benchmark. Local elasticity at the
  benchmark ~`0.026`.
- Lifetime births per household: `1.813` at the benchmark, ceiling `1.853`
  under nearly free housing; replacement in household units is `2.115`
  (conversion 0.5, maturation leak `0.0544`). The schedule never reaches
  replacement: at this parameterization there is no closed-economy reproductive
  steady state (theory note, Prop. 5 dichotomy — the economy is sustained by
  the entrant inflow).
- Through Prop. 2 of the theory note, the trace implies an effective housing
  share of the marginal child's cost of ~`2.2%`, against a reproduction gap of
  `14.3%`: the stabilization condition `s_H >= gap` fails by a factor of ~6.5.
  Per the note's Section 5, this flatness is exactly the object that stationary
  fertility targets do not discipline — treat it as the current
  parameterization's implied slope, not as an established fact about fertility
  and housing costs.
- `D(r)` falls from 9.73 to 2.91 rooms per household across the grid, so the
  clearing scale `S(r)` is strictly increasing (0.001 to 16.7 relative to
  benchmark supply units).

## Caveats

- Smoke grade: at the benchmark price the implied clearing scale is `1.0098`
  rather than 1.0000, a ~1% level error from the coarse wealth grid. Shapes
  and demographic ratios are reliable (benchmark `B/E` matches production to
  4 digits); levels of `D` and `S` are not production numbers.
- The transfer is not re-balanced away from the benchmark price, so the fiscal
  budget is unbalanced off-benchmark by construction.
- One arm only; the certified maturation-repair arm sits lower
  (`B/E = 0.788` at its benchmark) and would trace lower everywhere.
