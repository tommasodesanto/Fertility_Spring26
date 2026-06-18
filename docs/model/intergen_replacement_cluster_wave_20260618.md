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
