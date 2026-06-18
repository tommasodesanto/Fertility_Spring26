# Intergen Replacement Cluster Wave

Date: 2026-06-18

This wave tests whether the worse fit under `candidate_replacement_v1` is a
search problem, a target-object problem, or a model-mechanism problem. It keeps
the SMM count discipline: each tested target set has 13 target moments for the
13 varied internal parameters.

## Source Code

Committed local source state: `d92174c`.

Torch scratch copy:
`/scratch/td2248/projects/Fertility_Spring26_20260617_fast`

The scratch copy is not a Git checkout. The relevant source files and the
Jacobian source records were synced with `rsync -R` before submission.

## Target Sets

All runs use the one-market/no-location intergenerational model with
`J=16`, `Nb=60`, `income_states=5`, `n_house=5`, and `max_iter_eq=3`.

| Target set | Purpose | Slurm job |
|---|---|---:|
| `candidate_replacement_v1` | Current candidate: PSID old nonhousing mean wealth level plus parent-childless nonhousing gap | `11095203` |
| `candidate_replacement_nh_median_v1` | Same 13-count system, but old nonhousing level is the PSID median matched to the model median statistic | `11095204` |
| `candidate_replacement_total_median_v1` | Same 13-count system, but old wealth uses total-wealth median and total parent-childless wealth gap | `11095205` |

Each global-DE job is an 8-task Slurm array on `cpu_short`, with
`INTERGEN_GLOBAL_EVALS_PER_TASK=300`, `INTERGEN_GLOBAL_POP_SIZE=24`, and
`INTERGEN_MINUTES=115`.

## Jacobian Audit

Slurm job `11095206` runs `code/model/tools/audit_intergen_sensitivity_jacobian.py`
under `candidate_replacement_v1`.

It audits two explicit source records from the local replacement panel:

| Point | Reason |
|---|---|
| `replacement_best_draw_0055` | Best scalar fit under `candidate_replacement_v1`: high ownership, low old nonhousing wealth |
| `replacement_high_oldwealth_draw_0064` | Highest old nonhousing wealth basin in the local panel: old wealth close to target, ownership collapses |

The goal is to compare the local Jacobian in the high-ownership basin against
the high-old-wealth basin. The key question is whether the old wealth, ownership,
and room-separation moments load on distinct parameter directions or collapse
onto the same few directions.

## Evaluation Plan

When jobs finish:

1. Collect each global-DE array with `code/model/tools/collect_intergen_panel_results.py`.
2. Compare best and frontier records across all three target sets, not just
   scalar loss.
3. Check whether the median old-wealth target produces a basin with moderate
   ownership, old wealth near target, and credible owner-renter room separation.
4. Inspect Jacobian rank, condition number, top moment sensitivities, and column
   collinearities for `11095206`.
5. If no target set can jointly fit old wealth, ownership, and room separation,
   report the mechanism before proposing any target replacement.

This wave is diagnostic. It does not define a production SMM target system.

## Results

Collected locally under:
`output/model/cluster_pulls/intergen_replacement_cluster_wave_20260618/`

All Slurm jobs completed with exit code `0:0`. The three global-DE arrays each
produced 2,400 finite records. The Jacobian job completed with empty stderr.

| Target set | Best rank loss | Main read |
|---|---:|---|
| `candidate_replacement_v1` | 26.553 | Better than the local panel best `34.009`, but old nonhousing mean wealth remains far below target |
| `candidate_replacement_nh_median_v1` | 9.093 | Large scalar improvement; median old-wealth target is much more compatible with the model |
| `candidate_replacement_total_median_v1` | 13.430 | Total-wealth median can be matched, but ownership and first-birth housing response deteriorate |

Best-case moment summaries:

| Moment | Mean old NH target | NH median target | Total median target |
|---|---:|---:|---:|
| Rank loss | 26.553 | 9.093 | 13.430 |
| TFR | 1.651 | 1.882 | 1.599 |
| Childless rate | 0.283 | 0.248 | 0.345 |
| Ownership | 0.537 | 0.509 | 0.393 |
| Family ownership gap | 0.219 | 0.362 | 0.249 |
| First-birth room increment | 1.069 | 0.956 | 1.334 |
| Second-birth room increment | 0.231 | 0.439 | 0.643 |
| Young liquid wealth/income | 0.223 | 0.213 | 0.267 |
| Old nonhousing wealth/income mean | 2.590 | 1.704 | 2.590 |
| Old nonhousing wealth/income median | 2.516 | 1.601 | 2.516 |
| Old total wealth/income median | 5.205 | 4.121 | 5.364 |
| Childless renter mean rooms | 4.780 | 4.052 | 3.948 |
| Childless owner mean rooms | 5.361 | 5.556 | 5.781 |
| Owner-renter room gap | 0.581 | 1.505 | 1.833 |

Frontier slice:

- Mean old nonhousing target: only 19 of 2,400 records have ownership above
  0.50 and old nonhousing mean wealth within 25 percent of the target, but their
  best scalar fit has a negative owner-renter room gap. No record also has an
  owner-renter room gap above 1.5.
- Nonhousing median target: 414 records have ownership above 0.50 and old
  nonhousing median wealth within 25 percent of target. The best scalar fit has
  a room gap of 1.50 but old nonhousing median wealth is still low at 1.60
  versus target 2.23. No record jointly satisfies ownership above 0.50, old
  median wealth within 25 percent, and room gap above 1.5.
- Total median target: 751 records have ownership above 0.50 and old total
  median wealth within 25 percent of target. The best scalar fit has ownership
  0.623, total median wealth 4.51 versus target 5.26, and room gap 1.48. This
  is close on the housing distribution but still misses fertility timing and
  first-birth room growth.

## Jacobian Results

The Jacobian job `11095206` audited two basins under `candidate_replacement_v1`.

| Point | Rank | Condition number | Interpretation |
|---|---:|---:|---|
| `replacement_best_draw_0055` | 13 | 9.41e3 | Formally full rank but ill-conditioned |
| `replacement_high_oldwealth_draw_0064` | 12 | 4.46e8 | Locally rank deficient; high-old-wealth basin has extremely weak identification |

The strongest local load for old nonhousing wealth is `beta`, not the bequest
parameters. At `replacement_best_draw_0055`, the scaled derivative of old
nonhousing wealth with respect to `beta` is 2.19; for `theta0` it is only
0.069. The parent-childless old wealth gap is moved by `beta`, `c_bar_0`, and
only weakly by `theta_n`. This means the current old-wealth targets are not
cleanly identifying the intended bequest block.

The strongest collinearities are also economically important:

- `alpha_cons` and `chi` have cosine 0.947 at the best basin;
- `c_bar_0` and `c_bar_n` have cosine 0.936;
- `alpha_cons` and `h_bar_0` have cosine -0.929;
- in the high-old-wealth basin, `alpha_cons` and `h_bar_0` have cosine -0.976.

Interpretation: the replacement moments improve the target system, especially
when old wealth is treated as a median object, but they do not yet give a clean
identified SMM calibration. The room-distribution moments, fertility
increments, and old-wealth block are still loading on overlapping parameter
directions.
