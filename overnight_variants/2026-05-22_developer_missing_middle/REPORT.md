# Developer Missing-Middle Smoke Report

Verdict: **red**

## What Ran

The copied solver was modified so owner budgets use \(p_{iq(k)}\), owner
maintenance uses \(p_{iq(k)}\), sale/down-payment accounting uses the same
type price, and renter choices use a two-type rent approximation based on the
state's housing need. A small outer loop pinned ordinary prices to the current
benchmark and updated only middle-unit prices from the developer inverse supply
schedule.

This is therefore a **partial** type-specific market-clearing test, not a full
2-location by 2-type general-equilibrium fixed point.

## Developer Block

- active unit types: `S` ordinary and `M` family-capable middle
- middle boundary: `5+ rooms`
- \(F_{iq}\): documented and included in entry-threshold diagnostics, but not
  used to shut types off in the first clearing loop
- operational wedges: \(\tau_{iM}\), \(\nu_{iM}\), and \(\eta_{iM}\)
- ordinary prices pinned: `True`

## Clearing Status

- middle-price iterations: `3`
- last relative middle-price step: `0.0170414`
- final max type excess relative to demand: `1.09283e+06`
- elapsed seconds: `10.79`
- runtime category: `moderate-to-expensive` for partial PE type loop; full GE
  type clearing remains `expensive`
- SMM loss on partial smoke moments: `61.043`

## Key Moment Movement

| Moment | Target | Benchmark | Variant |
|---|---:|---:|---:|
| childless renter median rooms | 4.000 | 6.365 | 6.552 |
| childless owner median rooms | 6.000 | 6.800 | 6.800 |
| H01 | 0.664 | 0.441 | 0.740 |
| H12 | 0.566 | 0.192 | 0.314 |
| ownership | 0.627 | 0.643 | 0.000 |

## Interpretation

The branch is red in this smoke test because the partial type loop did not
deliver a usable market allocation: ordinary prices were pinned, type excess
remained large, and ownership collapsed. The diagnostic CSV reports type
prices, stocks, demands, middle wedges, room bins, renter cap mass,
lower/middle owner mass, and whether the H2/starter-family margin is alive.
