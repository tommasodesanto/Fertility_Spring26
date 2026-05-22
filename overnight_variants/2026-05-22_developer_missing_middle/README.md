# Developer Missing-Middle Prototype

This folder is an isolated overnight branch test. It contains a copied
`dt_cp_model/` package and does not modify the live implementation under
`code/model/dt_cp_model`.

The smoke driver is:

```bash
cd overnight_variants/2026-05-22_developer_missing_middle
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_developer_missing_middle_smoke.py --quiet
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
