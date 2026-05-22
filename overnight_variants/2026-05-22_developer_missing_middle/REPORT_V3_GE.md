# Developer Missing-Middle V3 Full Equilibrium Report

Verdict: **yellow**

## What This Is

This supersedes the V2 partial test. V3 updates every location/type price
\(p_{iq}\) and the entry shares in a copied fixed-point loop. Ordinary prices
are no longer pinned. The household solve at each iteration uses the current
type-price vector, charges renters by realized room interval \(q(h^R)\), then
developer inverse supply updates all \(S/M/L\) rents.

## Equilibrium Status

- accepted: `True`
- iterations: `4`
- final relative market excess: `0.0227967`
- final relative price step: `0.00807004`
- final entry-share step: `0.00280192`
- elapsed seconds: `133.51`
- runtime category: `expensive`
- SMM loss: `47.961`
- \(F_{iq}\): reported in entry thresholds but not used for discrete shutdown

## Type Markets

| Market | p | r | supply | demand |
|---|---:|---:|---:|---:|
| S loc 0 | 0.495 | 0.035 | 0.000 | 0.000 |
| M loc 0 | 0.514 | 0.036 | 0.547 | 0.541 |
| L loc 0 | 0.498 | 0.035 | 3.537 | 3.531 |
| S loc 1 | 0.584 | 0.041 | 0.234 | 0.231 |
| M loc 1 | 0.650 | 0.046 | 0.505 | 0.517 |
| L loc 1 | 0.587 | 0.041 | 2.233 | 2.237 |

## Moment Movement

| Moment | Target | Benchmark | V3 GE |
|---|---:|---:|---:|
| childless renter median rooms | 4.000 | 6.365 | 6.500 |
| childless owner median rooms | 6.000 | 6.800 | 8.200 |
| H01 | 0.664 | 0.441 | 0.418 |
| H12 | 0.566 | 0.192 | 0.252 |
| ownership | 0.627 | 0.643 | 0.578 |

## Read

This is now a genuine type-price fixed point rather than a pinned ordinary-price
test, and renter budgets now use realized room intervals. The verdict is yellow
only if the status above accepts; it is not green because the room medians,
ownership gradient, and \(H01\) remain away from the target discipline.
