# Intergen Full Global Room-Gap Search Readout

Date: 2026-06-23

Status: diagnostic readout, not a production calibration.

## Source

This note summarizes the 30-task split global search launched from the active
June 2026 one-market intergenerational package,
`code/model/intergen_housing_fertility/`. The old center-periphery
`dt_cp_model` path is not the calibration object here.

Cluster source root:

`/scratch/td2248/projects/Fertility_Spring26_20260617_fast/code/cluster`

Local compact output packet:

`output/model/cluster_pulls/intergen_full_global_roomgap_16h_readout_20260623/`

Important output files:

- `frontier_summary.json`: machine-readable all-case summary.
- `all_cases_compact.csv.gz`: compact table for all retained case records.
- `top200_cases.csv`: top scalar-loss records.
- `screen_grid.csv`: threshold-grid screen counts and best records.
- `frontier_bins.csv`: representative frontier bins.
- `best_loss_contributions.csv`: loss decomposition for the best scalar record.

## Run Design

Target set: `candidate_replacement_young_old_roomgap_v1`.

Grid and solver stack:

- `J=16`
- `Nb=60`
- `income_states=5`
- `n_house=5`
- `max_iter_eq=3`

Search design:

- 15 differential-evolution array tasks per wave.
- 15 local/random-panel array tasks per wave.
- 8 chained two-hour waves were submitted.
- Waves 1-7 completed for both algorithms.
- Wave 8 was cancelled by user request near the wall-clock end, after writing
  partial records.

Retained records:

- Total finite case records: `263,282`.
- DE records: `116,253`.
- Local/random records: `147,029`.
- All retained records have finite rank loss and `status == "ok"`.

This was a broad diagnostic search. It does not establish calibration progress
unless a retained case is re-solved and validated under the final intended
target/solver settings.

## Key Targets

The current frontier tension is mainly over:

| Moment | Target |
|---|---:|
| TFR | 1.700 |
| childless rate | 0.150 |
| young ownership, ages 25-34 | 0.341 |
| old ownership, ages 65-75 | 0.764 |
| childless owner-renter mean-room gap | 2.419 |
| old parent-childless nonhousing wealth gap | 1.007 |
| old nonhousing wealth-income median | 2.230 |

## Best Scalar Record

The best scalar-loss record is from DE wave 4, task 3, case 492.

| Moment | Model | Target |
|---|---:|---:|
| rank loss | 26.903 | |
| TFR | 1.764 | 1.700 |
| childless rate | 0.276 | 0.150 |
| young ownership | 0.175 | 0.341 |
| old ownership | 0.877 | 0.764 |
| childless owner-renter room gap | 1.493 | 2.419 |
| old parent-childless nonhousing wealth gap | 0.291 | 1.007 |
| old nonhousing wealth-income median | 2.059 | 2.230 |

Largest weighted loss contributions at this point:

| Moment | Contribution | Direction |
|---|---:|---|
| childless owner-renter room gap | 10.293 | too small |
| childless renter mean rooms | 3.417 | renters too large |
| first-birth housing increment | 2.783 | too large |
| young ownership | 2.200 | too low |
| young liquid wealth-income | 2.054 | too high |
| old ownership | 2.041 | too high |
| old parent-childless nonhousing wealth gap | 1.027 | too small |

This point is not economically satisfactory. It improves the scalar objective
relative to many random cases, but it still fails the central joint allocation
test.

## All-Case Screens

These counts use all `263,282` retained case records, not only per-task bests.

| Screen | Count | Best loss | Young own | Old own | Room gap | Childless | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|
| young own >= 0.20, old own <= 0.90, room gap >= 1.50 | 0 | | | | | | |
| young own >= 0.25, old own <= 0.88, room gap >= 1.50 | 0 | | | | | | |
| young own >= 0.20, old own <= 0.90, room gap >= 2.00 | 0 | | | | | | |
| young own >= 0.20, old own <= 0.90 | 16,157 | 47.618 | 0.206 | 0.897 | 1.319 | 0.169 | -0.827 |
| young own >= 0.20, room gap >= 1.50 | 14 | 47.207 | 0.247 | 0.994 | 1.537 | 0.207 | 0.489 |
| old own <= 0.90, room gap >= 1.50 | 27,217 | 29.204 | 0.000 | 0.872 | 2.017 | 0.283 | 0.079 |
| room gap >= 2.00 | 8,445 | 29.204 | 0.000 | 0.872 | 2.017 | 0.283 | 0.079 |
| young own >= 0.25 | 82,493 | 42.622 | 0.276 | 0.979 | 1.189 | 0.291 | 0.579 |
| old own <= 0.85 | 120,628 | 32.378 | 0.000 | 0.843 | 1.865 | 0.318 | 0.316 |

The sharp diagnostic result is stronger than the earlier task-best readout:
there are no retained cases with young ownership above 0.20, old ownership
below 0.90, and a room gap above 1.50. There are also no retained cases with
young ownership above 0.10, old ownership below 0.90, and a room gap above 1.50;
the only case with old ownership below 0.90 and room gap above 1.50 has young
ownership about 0.081 and very high scalar loss.

## Frontier Interpretation

The search confirms the existing frontier rather than finding a hidden basin.

The model can separately produce pieces of the desired allocation:

- room gap above 2.0 with old ownership below 0.90,
- young ownership above 0.25,
- old ownership below 0.85,
- positive old nonhousing wealth gaps.

But these pieces do not survive jointly:

- Large room-gap cases almost always empty young ownership.
- High young-ownership cases keep old ownership too high and shrink the
  owner-renter room gap.
- Cases with young ownership and room separation have old ownership near one.
- Old nonhousing wealth gaps improve in some regions but do not repair the
  young-ownership/old-exit/room-gap conflict.

The search therefore supports the current diagnosis: under the current
one-market Markov-income model, grid, and target system, the main failure is
not an obvious missed global-search basin. The unresolved issue is either a
model mechanism failure, a target-object mismatch, or a rough-solve artifact
that needs targeted verification rather than another blind search.

## Recommended Next Steps

1. Re-solve a small set of frontier records at higher verification quality:
   the best scalar point, the best `young >= 0.20 and old <= 0.90` point, the
   best `young >= 0.20 and room gap >= 1.50` point, the best `old <= 0.90 and
   room gap >= 1.50` point, and the closest joint-distance point.
2. Refit the surrogate/ML diagnostic on this current room-gap target system,
   using the compact all-case table, to map smooth local frontiers cheaply.
3. Audit the economic mechanism around the tenure/access and old-retention
   blocks before changing targets: down-payment/PTI thresholds, owner rung
   usage, renter cap mass, old downsizing, and bequest/nonhousing wealth
   composition.
4. Do not drop or demote a hard moment without an identification-preserving
   replacement for the same parameter block.

