# Developer Missing-Middle Prototype

This folder is an isolated overnight branch test. It contains a copied
`dt_cp_model/` package and does not modify the live implementation under
`code/model/dt_cp_model`.

The smoke driver is:

```bash
cd overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_smoke.py --quiet
```

The second-pass driver is:

```bash
cd overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_v2.py --quiet --iterations 5 --nb 50
```

The full-equilibrium type-price driver is:

```bash
cd overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_v3_ge.py --quiet --nb 30 --iterations 12 --price-damp 0.25 --entry-damp 0.25
```

## Implemented Objects

- Two unit types:
  - `S`: ordinary units
  - `M`: family-capable middle units, with the boundary at roughly 5+ rooms
- Owner rung mapping \(q(k)\): owner rungs with \(H_k\ge 5\) are treated as
  middle units in the two-type smoke test.
- Rental mapping \(q(h^R)\): renter choices are priced by a two-type
  approximation based on the state housing need in the Bellman step and by
  realized \(h^R\) in diagnostics.
- Type-specific owner budgets, owner maintenance, down-payment constraints, and
  sale proceeds in the copied solver.
- Developer inverse supply:
  \[
  r_{iq}=\tau_{iq}+\nu_{iq}S_{iq}^{1/\eta_{iq}}.
  \]

## Clearing Protocol

The first smoke test pins ordinary prices to the current benchmark scalar
prices and clears only middle-unit prices in a small outer loop. Fixed costs
\(F_{iq}\) are included in reported entry-threshold diagnostics but are not used
to shut down a type during this first clearing loop.

This is a partial test by design. `REPORT.md` states whether the partial
type-clearing loop stabilized and whether the room-screen diagnostics moved in
the right direction.

## Second-Pass V2

V1 failed partly because every owner rung above 5 rooms was priced as middle
housing. V2 adds a third large-unit type and bounds the middle type to the
5--6.5 room interval. It also first measures scalar-price baseline demand by
type, then calibrates \(\nu_{iq}\) so the type supply curve passes through that
baseline demand at the initial price. This avoids treating an arbitrary initial
stock split as the missing-middle wedge.

V2 is documented in `REPORT_V2.md`, `results_developer_missing_middle_v2.csv`,
and `diagnostics_developer_missing_middle_v2.csv`.

## Third-Pass V3 Full Equilibrium

`run_developer_missing_middle_v3_ge.py` supersedes the pinned-price tests for
branch-decision evidence. It updates all six location/type prices
\(p_{iq}\), \(q\in\{S,M,L\}\), and the two entry shares in an outer fixed-point
loop. Ordinary prices are not pinned.

The V3 correction also fixes the main household-side mismatch from the earlier
prototype: renter budgets are solved piecewise over realized room intervals, so
the budget uses \(r_{iq(h^R)}h^R\) rather than assigning a rent type from
\(\bar h_n\). Owner budgets, maintenance, down-payment constraints, and sale
proceeds continue to use \(p_{iq(k)}\).

The coarse `Nb=30` run accepted after 4 price iterations, with final relative
market excess 0.0228 and elapsed time 123.30 seconds. The verdict is yellow:
the type-specific GE loop clears and middle demand is non-degenerate, but the
room screen does not improve and the ownership gradient remains far below the
target. See `REPORT_V3_GE.md`,
`results_developer_missing_middle_v3_ge.csv`, and
`diagnostics_developer_missing_middle_v3_ge.csv`.
