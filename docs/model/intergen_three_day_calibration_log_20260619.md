# Intergen Three-Day Calibration Log

Date started: 2026-06-19

## Operating Rules

- Diagnostic exploration only; no solver/model-code changes without explicit
  approval.
- Do not delete cluster or local outputs.
- Preserve SMM identification: if a hard target is removed or weakened, add a
  replacement moment for the same parameter block or state the system is
  underidentified.
- Record large misses as mechanism evidence, not just search failures.

## Overnight DUE-Lifecycle Suite

Scratch root:
`/scratch/td2248/projects/Fertility_Spring26_20260617_fast/code/cluster`

All jobs below completed with exit code `0:0`.

| Run | Slurm | Target set | Cases | Best loss | Key read |
|---|---:|---|---:|---:|---|
| Pulse | 11108589 | `candidate_replacement_due_lifecycle_v1` | 2,400 | 26.216 | Better than local partial run, but high ownership and slope remain |
| Strict lifecycle | 11108949 | `candidate_replacement_due_lifecycle_v1` | 11,200 | 20.422 | Best strict run fits old median wealth but keeps old ownership too high |
| Soft lifecycle | 11108950 | `candidate_replacement_due_lifecycle_soft_v1` | 11,200 | 17.569 | Best scalar fit, but it ignores the lifecycle slope |
| Old ownership gap | 11108951 | `candidate_replacement_due_lifecycle_owngap_v1` | 11,200 | 23.214 | Parent-childless ownership gap can be matched, old ownership still high |
| Total wealth | 11108952 | `candidate_replacement_total_due_lifecycle_v1` | 10,976 | 24.310 | Total-wealth robustness does not fix lifecycle ownership |

### Current Best By Own Objective

The soft-lifecycle run has the lowest scalar loss, `17.569`, but this is not a
convincing economic fit:

| Moment | Model | Target/benchmark |
|---|---:|---:|
| TFR | 1.723 | 1.700 |
| Childlessness | 0.253 | 0.150 |
| Ownership 30--55 | 0.514 | 0.575 |
| Young ownership 25--34 | 0.068 | 0.341 |
| Old ownership 65--75 | 0.979 | 0.764 |
| Old-minus-young ownership | 0.910 | 0.423 |
| Old NH wealth median/income | 1.830 | 2.230 |
| Parent-childless old NH wealth gap | 0.715 | 1.007 |
| Childless renter rooms | 4.434 | 3.805 |
| Childless owner rooms | 5.482 | 6.224 |

The strict lifecycle run is economically more informative despite higher loss:

| Moment | Model | Target/benchmark |
|---|---:|---:|
| TFR | 2.048 | 1.700 |
| Childlessness | 0.190 | 0.150 |
| Ownership 30--55 | 0.801 | 0.575 |
| Young ownership 25--34 | 0.369 | 0.341 |
| Old ownership 65--75 | 0.975 | 0.764 |
| Old-minus-young ownership | 0.606 | 0.423 |
| Old NH wealth median/income | 2.288 | 2.230 |
| Parent-childless old NH wealth gap | 0.465 | 1.007 |
| Childless renter rooms | 4.166 | 3.805 |
| Childless owner rooms | 5.331 | 6.224 |

Interpretation: once young ownership is close to target, old ownership remains
far too high. This points to an old-retention/downsizing mechanism problem, not
only a broad lifecycle-slope weighting problem.

## June 19 Wave 1

Launched after reading the overnight results.

| Slurm | Run tag | Target set | Purpose |
|---:|---|---|---|
| 11136888 | `intergen_old_retention_wave1_20260619` | `candidate_replacement_old_retention_v1` | Directly target old ownership and old parent-childless gap |
| 11136889 | `intergen_young_old_own_wave1_20260619` | `candidate_replacement_young_old_own_v1` | Separately target young and old ownership |
| 11136890 | `intergen_young_old_roomgap_wave1_20260619` | `candidate_replacement_young_old_roomgap_v1` | Add owner-renter room-gap pressure |

All three arrays started and produced finite initial losses. Each uses
`16` tasks with `%8` concurrency, `J=16`, `Nb=60`, five income states,
`n_house=5`, `max_iter_eq=3`, `700` max evaluations per task, and a `115`
minute task budget.

## Failure Modes To Track

1. Old ownership remains too high even when \(\theta_0\) is low or old
   ownership is directly targeted.
2. Young ownership and old ownership move in opposite directions through
   tenure/access parameters, making a lifecycle slope alone too coarse.
3. Owner-renter room separation remains too small: owners are too small and
   renters too large relative to ACS room targets.
4. Parent-childless old wealth gaps remain below target, so \(\theta_n\) is
   still weakly identified by current old wealth moments.
