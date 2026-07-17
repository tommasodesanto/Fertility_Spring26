# Intergen Room-Gap 14-Moment Rescore

Date: 2026-06-23

Status: diagnostic rescore of existing cases. This did not rerun the model and
does not change solver or calibration code.

## Source

Input packet:

`output/model/cluster_pulls/intergen_full_global_roomgap_16h_readout_20260623/all_cases_compact.csv.gz`

Output packet:

`output/model/cluster_pulls/intergen_full_global_roomgap_14moment_rescore_20260623/`

Files:

- `rescore_summary.json`
- `top200_14moment_cases.csv`
- `screen_stats_14moment.csv`
- `best14_contributions_with_diagnostics.csv`
- `README.md`

The rescore uses all `263,282` retained finite case records from the full
global room-gap search.

## Proposed 14 Hard Moments

The aim is an \(n+1\)-style diagnostic objective: 14 hard moments for the 13
searched parameters.

| Moment | Target | Weight | Best 14-moment model |
|---|---:|---:|---:|
| TFR | 1.700 | 20.0 | 1.764 |
| childless rate | 0.150 | 20.0 | 0.276 |
| overall ownership | 0.575 | 100.0 | 0.490 |
| parent-child ownership gap | 0.168 | 45.0 | 0.291 |
| housing increment 0 to 1 child | 0.664 | 14.0 | 1.110 |
| housing increment 1 to 2 children | 0.488 | 8.0 | 0.331 |
| young liquid wealth-income | 0.179 | 12.0 | 0.593 |
| old nonhousing wealth-income median | 2.230 | 0.8 | 2.059 |
| old parent-childless nonhousing wealth gap | 1.007 | 2.0 | 0.291 |
| childless renter mean rooms | 3.805 | 6.0 | 4.560 |
| childless owner-renter mean-room gap | 2.419 | 12.0 | 1.493 |
| childless owner share rooms >= 6 | 0.596 | 25.0 | 0.762 |
| old ownership | 0.764 | 160.0 | 0.877 |
| young ownership, ages 25-34 | 0.341 | 80.0 | 0.175 |

Demoted diagnostics:

| Moment | Target | Best 14-moment model |
|---|---:|---:|
| childless owner mean rooms | 6.224 | 6.053 |
| childless renter share rooms >= 6 | 0.138 | 0.100 |
| old parent-childless ownership gap | 0.083 | 0.146 |

This preserves identification coverage for fertility, tenure/access, housing
space, lifecycle saving, old retention, and old wealth composition while
avoiding hard-targeting every correlated room and old-family object.

## Rescore Result

The best record under the 14-moment objective is the same record as under the
17-moment objective:

- source: DE wave 4, task 3, case 492
- 14-moment loss: `26.529963`
- original 17-moment loss at same point: `26.902664`

The three demoted moments contributed only `0.373` loss at the best original
point, so dropping them changes the scalar level but not the selected basin.

Largest hard-moment loss contributions at the best 14-moment point:

| Moment | Contribution | Direction |
|---|---:|---|
| childless owner-renter mean-room gap | 10.293 | too small |
| childless renter mean rooms | 3.417 | renters too large |
| housing increment 0 to 1 child | 2.783 | too large |
| young ownership | 2.200 | too low |
| young liquid wealth-income | 2.054 | too high |
| old ownership | 2.041 | too high |
| old parent-childless nonhousing wealth gap | 1.027 | too small |

## Screens Under The 14-Moment Objective

The joint feasibility result is unchanged because it is about the records
themselves, not the objective weights.

| Screen | Count | Best 14-loss | Young own | Old own | Room gap | Childless |
|---|---:|---:|---:|---:|---:|---:|
| young >= 0.20, old <= 0.90, gap >= 1.50 | 0 | | | | | |
| young >= 0.25, old <= 0.88, gap >= 1.50 | 0 | | | | | |
| young >= 0.20, old <= 0.90 | 16,157 | 46.367 | 0.206 | 0.897 | 1.319 | 0.169 |
| young >= 0.20, gap >= 1.50 | 14 | 42.912 | 0.412 | 1.000 | 1.511 | 0.252 |
| old <= 0.90, gap >= 1.50 | 27,217 | 29.012 | 0.000 | 0.872 | 2.017 | 0.283 |
| gap >= 2.00 | 8,445 | 29.012 | 0.000 | 0.872 | 2.017 | 0.283 |
| young >= 0.25 | 82,493 | 36.674 | 0.265 | 0.979 | 1.242 | 0.246 |
| old <= 0.85 | 120,628 | 31.574 | 0.000 | 0.843 | 1.865 | 0.318 |

## Interpretation

Reducing the objective from 17 to this 14-moment system is conceptually cleaner,
but it does not reveal a hidden good point in the existing 263k-case evidence.
The demoted moments were not the binding source of the current scalar best or
of the joint frontier failure.

The hard remaining failures are still:

- young ownership is too low,
- old ownership is too high,
- the owner-renter room gap is too small at the best scalar point,
- room-gap/old-exit cases nearly eliminate young ownership,
- young-ownership cases keep old ownership too high and shrink the room gap.

This does not rule out an optimizer run directly using the 14-moment objective:
the DE chains evolved under the 17-moment loss. But because the local/random
component was broad and the all-case screen is empty, another blind global run
under this specific 14-moment loss is unlikely to change the qualitative
frontier by itself.

Recommended next step: if this 14-moment system is the intended diagnostic
objective, add it explicitly as a target set and use it for a small focused
verification or active-learning run, not another unconstrained overnight search.

