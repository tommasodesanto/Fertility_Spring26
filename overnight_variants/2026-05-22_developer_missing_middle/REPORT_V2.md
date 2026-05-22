# Developer Missing-Middle V2 Report

Verdict: **yellow**

## What Changed Relative To V1

V1 treated every owner rung above 5 rooms as middle housing. That made the
middle price act like a tax on almost all ownership and collapsed the owner
sector. V2 uses three unit types:

- `S`: non-middle ordinary units, below 5 rooms;
- `M`: missing-middle units in the 5--6.5 room interval;
- `L`: large units above 6.5 rooms.

The V2 loop first measures scalar-price baseline demand by type. It then
chooses \(\nu_{iq}\) so each type's supply equals that baseline demand at the
initial type price. This avoids interpreting arbitrary initial stocks as a
developer wedge.

## Status

- active type price iterations: `1`
- last active relative price step: `6.31161e-05`
- final active type excess relative to demand: `0.000206536`
- ordinary price pinned: `True`
- \(F_{iq}\): reported in entry thresholds but not used for discrete shutdown
- elapsed seconds: `8.25`
- SMM loss: `48.5582`

## Moment Movement

| Moment | Target | Benchmark | V2 |
|---|---:|---:|---:|
| childless renter median rooms | 4.000 | 6.365 | 6.232 |
| childless owner median rooms | 6.000 | 6.800 | 8.200 |
| H01 | 0.664 | 0.441 | 0.382 |
| H12 | 0.566 | 0.192 | 0.259 |
| ownership | 0.627 | 0.643 | 0.606 |

## Read

This is a better engineering test than V1 because the mapping is no longer
taxing the whole owner ladder as middle housing. It is still not green: the
type-clearing mechanics work, but owner median rooms move the wrong way and
\(H01\) falls relative to the current benchmark. Ordinary prices are also
pinned, the household state is unchanged, and the type mapping is piecewise by
room bins rather than estimated from ACS unit-type prices.
