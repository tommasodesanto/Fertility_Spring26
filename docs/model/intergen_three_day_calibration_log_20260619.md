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

## Three-Day Campaign Plan

The campaign is a bounded diagnostic search, not a production calibration. Each
active SMM specification must preserve at least one informative target for each
of the 13 varied internal parameters:

\[
\{\beta,\alpha,b_0,\bar c_0,\bar c_n,\bar h_0,\bar h_{\mathrm{jump}},
\bar h_n,\psi_{\mathrm{child}},\kappa_n,\chi,\theta_0,\theta_n\}.
\]

The monitor should use the current wave first, then decide whether the next
wave should be a longer search on an existing identified target set or a new
13-moment diagnostic variant. A moment can be moved only if the log records
which parameter block it was identifying and what replacement moment carries
that identification. Do not treat frontier slices with fewer than 13 hard
moments as calibrated SMM specifications.

The leading mechanism question is whether the one-market model can jointly
deliver:

1. young ownership near the ACS/PUMS target,
2. old ownership below the current high-retention basin,
3. owner-renter room separation near the ACS room gap,
4. fertility and childlessness near their completed-fertility targets, and
5. old-age nonhousing wealth and parent-childless old-wealth gaps near PSID.

If a wave cannot improve that joint frontier, record why in economic terms:
search failure, target-measurement mismatch, weak identification, or a missing
model mechanism.

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

### 2026-06-19 09:38 EDT Partial Readout

The three arrays were running normally at the first hourly check: tasks `1--8`
active for each array and tasks `9--16` pending behind the `%8` throttle. This
is an early partial readout, mostly initial population records, so it should not
be interpreted as the final frontier.

| Run | Finite cases | Best partial loss | Key warning |
|---|---:|---:|---|
| `old_retention` | 268 | 27.277 | Old ownership still `0.965` against target `0.764` |
| `young_old_own` | 242 | 30.401 | Young ownership `0.039` and old ownership `0.959` both miss badly |
| `young_old_roomgap` | 256 | 40.182 | Room gap improves to `1.405`, but target is `2.419` and old ownership remains `0.955` |

The early mechanism signal is the same as the overnight readout: ownership
among the old remains very high even when directly targeted, and increasing the
owner-renter room gap does not by itself fix the lifecycle ownership slope.
Keep monitoring to distinguish a search-stage artifact from a genuine old
retention/downsizing mechanism failure.

### 2026-06-19 09:40 EDT Frontier Slice

The arrays were still running normally, but the growing JSONL files made a more
useful partial frontier possible. The best scalar points were:

| Run | Finite cases | Best partial loss | TFR | Own 25--34 | Old own | Room gap |
|---|---:|---:|---:|---:|---:|---:|
| `old_retention` | 680 | 26.297 | 1.722 | 0.046 | 0.961 | 0.333 |
| `young_old_own` | 605 | 30.401 | 1.872 | 0.039 | 0.959 | 0.854 |
| `young_old_roomgap` | 630 | 40.182 | 1.959 | 0.213 | 0.955 | 1.405 |

The frontier slices sharpen the mechanism. Low-old-ownership points do exist:
`old_retention` finds old ownership `0.772`, `young_old_own` finds `0.755`,
and `young_old_roomgap` finds `0.742`. But those points generally have almost
zero young ownership, poor fertility, or weak room separation. Requiring
young ownership above `0.25` and old ownership below `0.85` leaves only a
handful of high-loss candidates, and the best such candidates either collapse
the room gap or fertility. No partial run yet has a candidate with owner mean
rooms at least `6` and renter mean rooms at most `4`.

Current interpretation: the model can mechanically reduce old ownership, but
not yet while also preserving young access, fertility, and owner-renter space
separation. That is stronger evidence for a joint mechanism/target tension than
for a simple weighting failure, though the wave is still too young for a final
claim.

### 2026-06-19 09:45 EDT Partial Update

The first batch of tasks was still active after about 11 minutes; tasks `9--16`
had not yet started. The `old_retention` target set found a much better scalar
partial point, but the mechanism is not yet acceptable:

| Run | Finite cases | Best partial loss | TFR | Own 25--34 | Old own | Room gap | Old NH median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 811 | 18.389 | 1.665 | 0.002 | 0.889 | 1.323 | 2.974 |
| `young_old_own` | 716 | 30.401 | 1.872 | 0.039 | 0.959 | 0.854 | 0.458 |
| `young_old_roomgap` | 753 | 40.182 | 1.959 | 0.213 | 0.955 | 1.405 | 0.915 |

This improves the scalar objective for `old_retention`, and old ownership falls
from the previous best around `0.96` to `0.889`. But young ownership is almost
zero, while the target is `0.341`; the model is again lowering old ownership by
nearly emptying the young owner pipeline. Frontier restrictions confirm the
tradeoff: requiring old ownership below `0.85` and young ownership above `0.25`
leaves only three `old_retention` candidates so far, with the best at high loss
`47.253`, TFR `0.926`, childlessness `0.541`, and room gap `0.070`.

Current mechanism note: separate old-retention pressure is useful, but the
model still has trouble combining (i) young access, (ii) lower old retention,
and (iii) owner-renter room separation. The next full readout should inspect
whether later DE generations can repair this, or whether the next wave needs a
diagnostic that conditions on young access while varying old-retention
mechanisms.

### 2026-06-19 10:01 EDT Partial Frontier

The first batch of tasks was still active after about 28 minutes, with
tasks `9--16` pending behind the `%8` array throttle. All arrays were healthy
in `squeue`/`sacct`.

| Run | Finite cases | Best partial loss | TFR | Own 25--34 | Old own | Room gap | Old NH median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 1,817 | 17.600 | 1.995 | 0.011 | 0.960 | 1.073 | 2.059 |
| `young_old_own` | 1,623 | 26.373 | 1.846 | 0.024 | 0.969 | 1.084 | 1.601 |
| `young_old_roomgap` | 1,767 | 40.182 | 1.959 | 0.213 | 0.955 | 1.405 | 0.915 |

The economically useful frontier did not improve in the scalar-best points.
The best `old_retention` point now fits some fertility and old-wealth objects
better, but it still has almost no young ownership. The room-gap run found one
candidate with owner mean rooms above `6` and renter mean rooms below `4`
(`6.087` versus `3.822`), but that candidate has young ownership `0.000`,
old ownership `0.001`, TFR `1.372`, and high loss `146.721`. This is evidence
that room separation can be forced only by leaving the tenure distribution
economically empty, not yet by a plausible owner-renter allocation.

### 2026-06-19 10:16 EDT Partial Frontier

The first batch of tasks was still running normally after about 42 minutes.
Scalar-best points remained economically weak, but the constrained frontier
moved in an informative direction. In the `young_old_roomgap` run, the best
candidate satisfying old ownership below `0.85` and young ownership above
`0.25` improved to loss `100.158`, with young ownership `0.435` and old
ownership `0.778`. That point is not a good fit: TFR is `1.511`,
childlessness is `0.327`, and the owner-renter room gap is only `0.447`.

Interpretation: the model can now find points with plausible young access and
lower old retention, but those points erase the owner-renter space separation.
This is the converse of the high-room-gap frontier, which achieves room
separation only in near-zero-ownership corners. The live failure is therefore
the joint interaction between tenure timing and space separation, not just a
failure to find either margin separately.

### 2026-06-19 10:20 EDT Partial Frontier

The first batch of tasks was still running after about 46 minutes; tasks
`9--16` remained pending behind the array throttle. No Slurm failure was visible
in `squeue` or `sacct`.

| Run | Finite cases | Best partial loss | TFR | Own 25--34 | Old own | Room gap | Old NH median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 2,984 | 17.600 | 1.995 | 0.011 | 0.960 | 1.073 | 2.059 |
| `young_old_own` | 2,662 | 26.373 | 1.846 | 0.024 | 0.969 | 1.084 | 1.601 |
| `young_old_roomgap` | 2,939 | 38.021 | 1.480 | 0.009 | 0.953 | 1.683 | 1.830 |

The best constrained frontier points are still economically far from target.
The best `old_retention` point with old ownership below `0.85` has old
ownership `0.818` and room gap `1.531`, but young ownership is zero and
childlessness is `0.348`. The best point with both old ownership below `0.85`
and young ownership above `0.25` remains high loss: in `young_old_roomgap`,
young ownership is `0.435` and old ownership is `0.778`, but the room gap is
only `0.447`, TFR is `1.511`, and childlessness is `0.327`.

Current decision: let Wave 1 complete before launching another array. The next
wave should target this joint frontier directly rather than simply increase
scalar search budget on a generic objective.

### 2026-06-19 10:23 EDT Health Check

The same first-batch tasks were still running normally after about 49 minutes.
All visible stderr files were zero bytes. The case counts increased to
`3,196`, `2,856`, and `3,150` finite records across `old_retention`,
`young_old_own`, and `young_old_roomgap`, respectively.

There was no meaningful improvement in the joint frontier. No current candidate
satisfies the relaxed screen:

\[
\text{old ownership}\le 0.85,\qquad
\text{young ownership}\ge 0.25,\qquad
\text{owner-renter room gap}\ge 1.5.
\]

This reinforces the current mechanism read: the model can move each component
in isolation, but the searched parameter region has not yet produced young
access, old exit, and owner-renter space separation together.

## Failure Modes To Track

1. Old ownership remains too high even when \(\theta_0\) is low or old
   ownership is directly targeted.
2. Young ownership and old ownership move in opposite directions through
   tenure/access parameters, making a lifecycle slope alone too coarse.
3. Owner-renter room separation remains too small: owners are too small and
   renters too large relative to ACS room targets.
4. Parent-childless old wealth gaps remain below target, so \(\theta_n\) is
   still weakly identified by current old wealth moments.
