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

### 2026-06-19 10:55 EDT User Departure Baseline

The three arrays were still healthy when the user left for the three-day
campaign. First-batch tasks remained active, one `young_old_roomgap` task had
completed, and `young_old_roomgap` task `9` had started. All visible stderr
files were still zero bytes. The hourly heartbeat automation is active with
`72` hourly checks.

Current scalar-best points:

| Run | Cases | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 5,243 | 15.961 | 1.853 | 0.261 | 0.000 | 0.853 | 1.470 | 4.449 | 5.919 | 2.517 |
| `young_old_own` | 4,785 | 22.653 | 1.835 | 0.269 | 0.136 | 0.900 | 0.538 | 4.957 | 5.495 | 2.517 |
| `young_old_roomgap` | 5,228 | 34.920 | 1.853 | 0.246 | 0.144 | 0.956 | 1.507 | 4.108 | 5.615 | 1.373 |

Current constrained-frontier points:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Old ownership `<=0.85` | `old_retention` | 18.949 | 1.900 | 0.196 | 0.008 | 0.753 | 0.977 | 4.742 | 5.719 |
| Old `<=0.85`, young `>=0.25` | `old_retention` | 47.253 | 0.926 | 0.541 | 0.316 | 0.847 | 0.070 | 5.264 | 5.334 |
| Room gap `>=2` | `young_old_roomgap` | 51.042 | 1.603 | 0.329 | 0.000 | 0.651 | 2.123 | 4.702 | 6.825 |
| Owner rooms `>=6`, renter rooms `<=4` | `young_old_roomgap` | 146.721 | 1.372 | 0.349 | 0.000 | 0.001 | 2.264 | 3.822 | 6.087 |
| Joint softer: old `<=0.88`, young `>=0.20`, gap `>=1.0` | `young_old_roomgap` | 71.021 | 2.850 | 0.041 | 0.249 | 0.676 | 1.100 | 4.474 | 5.574 |

The economic interpretation is unchanged but sharper. The model can reduce old
ownership, can generate owner-renter room separation, and can produce some
young ownership. It has not yet produced those three margins together with
credible fertility and old-wealth moments. The immediate next action is to let
Wave 1 finish and then launch a bounded Wave 2 that targets the joint frontier,
not a blind scalar-loss continuation.

### 2026-06-19 10:58 EDT Pulse

Wave 1 remained healthy. Several first-batch tasks had completed and replacement
tasks had started; visible stderr files were still zero bytes. Because the
arrays were still running, no Wave 2 was launched.

The `old_retention` scalar best improved to loss `14.778`, with TFR `1.786`,
childlessness `0.254`, aggregate ownership `0.510`, young ownership `0.023`,
old ownership `0.936`, and owner-renter room gap `1.105`. This is a better
scalar point but not a better economic point: young ownership remains nearly
empty, and old ownership remains too high.

The constrained frontier was essentially unchanged:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Old ownership `<=0.85` | `old_retention` | 18.949 | 1.900 | 0.196 | 0.008 | 0.753 | 0.977 | 4.742 | 5.719 |
| Old `<=0.85`, young `>=0.25` | `old_retention` | 47.253 | 0.926 | 0.541 | 0.316 | 0.847 | 0.070 | 5.264 | 5.334 |
| Room gap `>=2` | `young_old_roomgap` | 51.042 | 1.603 | 0.329 | 0.000 | 0.651 | 2.123 | 4.702 | 6.825 |
| Joint softer: old `<=0.88`, young `>=0.20`, gap `>=1.0` | `young_old_roomgap` | 71.021 | 2.850 | 0.041 | 0.249 | 0.676 | 1.100 | 4.474 | 5.574 |

Adding a fertility screen makes the tension sharper: there is still no
candidate satisfying TFR in `[1.55,1.90]`, childlessness below `0.28`, old
ownership below `0.88`, young ownership above `0.20`, and room gap above `1.0`.
The current evidence says the search can find each margin separately, but not
the joint lifecycle-tenure-space-fertility object.

### 2026-06-19 11:00 EDT Monitor-Only Pulse

Wave 1 was still running cleanly, with replacement tasks beginning in the
second batch and all visible stderr files at zero bytes. No new wave was
launched.

The `old_retention` and `young_old_own` scalar bests were unchanged from the
previous pulse. The useful new movement was in `young_old_roomgap`: among
candidates with old ownership below `0.85`, the best point improved to loss
`41.961`, old ownership `0.850`, and room gap `1.819`. The miss remains large:
young ownership is `0.000`, TFR is only `1.353`, childlessness is `0.356`, and
the old parent-childless nonhousing wealth gap is negative. This is not yet a
joint solution; it is another example of improving old exit and room separation
by emptying young ownership and damaging fertility.

This reinforces the current mechanism read: the model can move each component
in isolation, but the searched parameter region has not yet produced young
access, old exit, and owner-renter space separation together.

### 2026-06-19 11:06 EDT Monitor Pulse

Wave 1 remained active and healthy. Several long first-batch tasks were still
running near the `90` minute mark, and some second-batch tasks remained pending
behind the array throttle. Visible stderr files were all zero bytes, so no
recovery action was needed and no Wave 2 was launched.

The only substantive movement was in the `young_old_roomgap` low-old-ownership
frontier. Among candidates with old ownership below `0.85`, the best point
improved to loss `38.022`, old ownership `0.829`, and room gap `1.599`.
This is informative but still not a joint fit: young ownership is only `0.012`,
TFR is `1.477`, childlessness is `0.311`, and the old nonhousing wealth median
is `2.059`. The result reinforces the same mechanism diagnosis: the model can
combine lower old ownership and more room separation, but it does so by nearly
emptying the young-owner pipeline and weakening fertility.

### 2026-06-19 10:25 EDT Frontier Check

The same first-batch tasks were still running after about 51 minutes; second
batch tasks remained pending behind the `%8` array throttle and stderr files
were still zero bytes. Finite records reached `3,311`, `2,968`, and `3,272`.

Best scalar points were unchanged. A softer joint screen was also empty:

\[
\text{old ownership}\le 0.88,\qquad
\text{young ownership}\ge 0.20,\qquad
\text{owner-renter room gap}\ge 1.0.
\]

Decision: do not launch an overlapping wave. Wave 1 remains healthy and has not
yet reached the second batch. The next useful action is to let the current
arrays roll forward, then reassess whether the completed frontier calls for a
longer search on an existing identified target set or a new 13-moment
diagnostic variant.

### 2026-06-19 10:27 EDT Partial Improvement

The first-batch tasks were still healthy after about 53 minutes. Finite records
reached `3,431`, `3,073`, and `3,391`. The `young_old_own` target set found a
new scalar-best point:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Old NH median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `young_old_own` | 22.653 | 1.835 | 0.269 | 0.136 | 0.900 | 0.538 | 2.516 |

This is a genuine improvement in the tenure-timing direction relative to the
previous `young_old_own` scalar best, but it does not solve the joint problem.
Old ownership remains too high, young ownership is still below the ACS target
of about `0.341`, and owner-renter room separation is far below the ACS gap of
about `2.419`. The relaxed and softer joint screens remain empty.

Decision unchanged: do not launch overlapping jobs. Let Wave 1 finish.

### 2026-06-19 10:28 EDT Scalar/Mechanism Divergence

The first-batch tasks were still running after about 55 minutes, with the second
batch pending. Stderr files remained zero bytes. The `old_retention` target set
found a slightly better scalar point:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Old NH median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 17.267 | 1.859 | 0.237 | 0.025 | 0.979 | 1.103 | 1.830 |

This lowers the scalar loss but worsens the old-retention mechanism: old
ownership rises further above target, while young ownership remains almost
empty. The best low-old-ownership frontier point in `old_retention` remains
young-ownership-zero, and the best young-and-old ownership point still erases
the room gap. The softer joint screen remains empty.

Decision unchanged: no overlapping launch while Wave 1 is healthy and still in
its first batch.

### 2026-06-19 10:30 EDT Room-Gap Frontier Improvement

The first-batch tasks were still running after about 57 minutes. The
`young_old_roomgap` target set found a materially better high-room-gap frontier
point:

| Screen | Loss | Own 25--34 | Old own | Room gap | TFR | Childless | Renter rooms | Owner rooms |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `roomgap >= 2` | 51.041 | 0.000 | 0.651 | 2.123 | 1.603 | 0.329 | 4.702 | 6.825 |

This point is useful because it separates two problems. The model can generate
large owner-renter room separation and low old ownership at the same time, but
only by emptying the young-owner pipeline. The best point with both young
ownership above `0.25` and old ownership below `0.85` still has a room gap of
only `0.447`. The joint and softer joint screens remain empty.

Decision unchanged: no overlapping launch. Wait for the second batch before
designing the next wave.

### 2026-06-19 10:33 EDT Low-Old-Ownership Frontier Improvement

The first-batch tasks were still running after about 60 minutes; the second
batch remained pending behind the array throttle and stderr files were still
zero bytes. The `old_retention` target set improved its best low-old-ownership
frontier point:

| Screen | Loss | Own 25--34 | Old own | Room gap | TFR | Childless | Renter rooms | Owner rooms |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `old <= .85` | 20.937 | 0.000 | 0.843 | 0.945 | 1.931 | 0.238 | 4.953 | 5.898 |

This is useful but not enough. It shows the search can lower old ownership
without destroying fertility, but it does so by eliminating young ownership and
with a room gap still far below the ACS target. The best point with both young
ownership above `0.25` and old ownership below `0.85` still has a near-zero
room gap. The live failure remains the three-way interaction among young access,
old exit, and owner-renter space separation.

Decision unchanged: no overlapping launch while the first batch is still
healthy.

### 2026-06-19 10:38 EDT Health And Best-So-Far

The first-batch tasks were still running after about 65 minutes. Tasks
`9--16` remained pending behind the `%8` array throttle, and all visible stderr
files were still zero bytes. Finite cases reached `4,237`, `3,773`, and
`4,198` for `old_retention`, `young_old_own`, and `young_old_roomgap`.

The current scalar best by each objective is:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Old NH median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 17.267 | 1.859 | 0.237 | 0.025 | 0.979 | 1.103 | 1.830 |
| `young_old_own` | 22.653 | 1.835 | 0.269 | 0.136 | 0.900 | 0.538 | 2.516 |
| `young_old_roomgap` | 34.919 | 1.853 | 0.246 | 0.144 | 0.956 | 1.507 | 1.373 |

The `young_old_roomgap` scalar best improved relative to the prior partial
readout, but the economic frontier is still not solved. A single softer joint
candidate now exists in `young_old_roomgap`, with young ownership `0.249`, old
ownership `0.676`, and room gap `1.100`, but it has TFR `2.850`,
childlessness `0.041`, old nonhousing median wealth `0.229`, and high loss
`71.021`. The strict joint screen remains empty.

Decision unchanged: let Wave 1 finish. The next wave should be designed around
the joint frontier, not around the scalar best alone.

### 2026-06-19 10:43 EDT Old-Retention Scalar And Frontier Move

The first-batch tasks were still running after about 69 minutes, with tasks
`9--16` pending behind the `%8` throttle and all visible stderr files still at
zero bytes. Finite records reached `4,516`, `4,071`, and `4,496` for
`old_retention`, `young_old_own`, and `young_old_roomgap`.

The `old_retention` target set found a better scalar point:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Old NH median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 15.961 | 1.853 | 0.261 | 0.000 | 0.853 | 1.470 | 2.516 |

It also improved the best low-old-ownership frontier:

| Screen | Loss | Own 25--34 | Old own | Room gap | TFR | Childless | Renter rooms | Owner rooms |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `old <= .85` | 20.085 | 0.000 | 0.813 | 1.066 | 1.513 | 0.308 | 4.693 | 5.759 |

The movement is informative but still not a solved fit. The scalar-best point
now nearly reaches the old-ownership target, but it does so with effectively no
young ownership. The low-old-ownership frontier similarly lowers old ownership
while missing fertility and leaving the young-owner pipeline empty. The best
point with both young ownership above `0.25` and old ownership below `0.85`
still has a near-zero room gap, and the strict joint screen remains empty.

Decision unchanged: continue Wave 1 rather than launching an overlapping
search. The next wave should condition more directly on preserving young access
while forcing old exit and room separation, or it should diagnose why that
frontier is infeasible in the current mechanism.

### 2026-06-19 11:12 EDT Second-Batch Monitor

Wave 1 remained healthy. The first eight `old_retention` tasks completed with
exit code `0:0`; its second batch was running. The `young_old_own` and
`young_old_roomgap` arrays still had late first-batch tasks running and a few
second-batch tasks pending behind the `%8` throttle. All visible stderr files
were still zero bytes.

The current scalar best by each objective was:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Old NH median |
|---|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 14.778 | 1.786 | 0.254 | 0.023 | 0.936 | 1.105 | 2.516 |
| `young_old_own` | 22.653 | 1.835 | 0.269 | 0.136 | 0.900 | 0.538 | 2.516 |
| `young_old_roomgap` | 31.956 | 1.598 | 0.316 | 0.004 | 0.952 | 1.760 | 1.830 |

The `young_old_roomgap` scalar fit improved from the earlier `34.919` point,
mostly by moving the renter room level closer to target and raising the room
gap. It is still not an economic solution: young ownership is essentially zero,
old ownership is still far above target, and childlessness is too high.

The constrained frontier still shows the same mechanism tension. Points with
low old ownership exist; points with large owner-renter room gaps exist; and
points with positive young ownership and low old ownership exist. But no run
yet has a candidate satisfying the strict joint screen

\[
\text{old ownership}\le 0.85,\qquad
\text{young ownership}\ge 0.25,\qquad
\text{owner-renter room gap}\ge 1.5,
\]

and no candidate satisfies the softer fertility-ok joint screen. The only
softer joint candidate remains high-fertility and low-childlessness, with
second-child housing growth negative. Do not launch an overlapping wave while
Wave 1 is healthy; the next wave should start only after the current arrays
clear.

### 2026-06-19 11:17 EDT High-Room-Gap Frontier Move

Wave 1 was still running cleanly, with no visible stderr failures. The scalar
best points did not change, and the strict joint screen remained empty. The
only material movement was in the `young_old_own` high-room-gap frontier:

| Screen | Loss | Own 25--34 | Old own | Room gap | TFR | Childless | Renter rooms | Owner rooms |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `roomgap >= 2` | 72.170 | 0.000 | 0.645 | 2.001 | 0.999 | 0.528 | 5.213 | 7.213 |

This improves the previous high-room-gap frontier under the `young_old_own`
objective, but it reinforces the same mechanism diagnosis. The model can
generate owner-renter room separation and low old ownership only by emptying
the young-owner pipeline and missing fertility badly. No current candidate
combines young ownership, old exit, room separation, and fertility in the
targeted range.

### 2026-06-19 11:25 EDT Three-Day Chain Setup

Wave 1 was still healthy when the next chain was queued. Remaining Wave 1 tasks
were running with zero visible stderr failures and no nonzero Slurm exit codes.
The partial finite-record counts were `7,274`, `6,869`, and `7,400` for
`old_retention`, `young_old_own`, and `young_old_roomgap`, respectively.

Current scalar-best points:

| Run | Cases | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 7,274 | 14.778 | 1.786 | 0.254 | 0.023 | 0.936 | 1.105 | 4.515 | 5.620 | 2.516 |
| `young_old_own` | 6,869 | 22.653 | 1.835 | 0.269 | 0.136 | 0.900 | 0.538 | 4.957 | 5.495 | 2.516 |
| `young_old_roomgap` | 7,400 | 31.956 | 1.598 | 0.316 | 0.004 | 0.952 | 1.760 | 3.949 | 5.709 | 1.830 |

Current frontier screens:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Old ownership `<=0.85` | `old_retention` | 18.949 | 1.900 | 0.196 | 0.008 | 0.753 | 0.977 | 4.742 | 5.719 |
| Old `<=0.85`, young `>=0.25` | `old_retention` | 47.253 | 0.926 | 0.541 | 0.316 | 0.847 | 0.070 | 5.264 | 5.334 |
| Room gap `>=2` | `young_old_roomgap` | 51.042 | 1.603 | 0.329 | 0.000 | 0.651 | 2.123 | 4.702 | 6.825 |
| Owner rooms `>=6`, renter rooms `<=4` | `young_old_roomgap` | 146.721 | 1.372 | 0.349 | 0.000 | 0.001 | 2.264 | 3.822 | 6.087 |
| Joint softer: old `<=0.88`, young `>=0.20`, gap `>=1.0` | `young_old_roomgap` | 71.021 | 2.850 | 0.041 | 0.249 | 0.676 | 1.100 | 4.474 | 5.574 |

Decision: do not change the repository model or target code. Queue two
dependent follow-up waves that preserve the same 13-moment identified target
systems but vary differential-evolution search pressure.

| Wave | Slurm jobs | Dependency | Target sets | Search pressure |
|---|---|---|---|---|
| Wave 2 broad | `11146362`, `11146363`, `11146364` | after Wave 1 jobs `11136888`, `11136889`, `11136890` | old retention, young/old ownership, young/old/room-gap | 24 tasks each, `%8`, `pop_size=40`, mutation `1.05`, crossover `0.80`, 1,100 requested evals/task, 115 minute budget |
| Wave 3 tight | `11146365`, `11146366`, `11146367` | after Wave 2 jobs | same three target sets | 24 tasks each, `%8`, `pop_size=24`, mutation `0.55`, crossover `0.95`, 1,100 requested evals/task, 115 minute budget |

The purpose is to separate search-depth and search-geometry failures from an
economic impossibility result. If these waves still cannot produce even a
relaxed joint candidate with young ownership, old exit, and room separation,
the next move should be either a scratch-only experimental 13-moment target
variant or a written mechanism diagnosis, not another blind scalar-loss run.

### 2026-06-19 11:29 EDT Decomposition Screens

Wave 1 was still running normally, with all visible stderr files at zero bytes.
Wave 2 and Wave 3 remained pending behind the intended Slurm dependencies.
Finite-record counts rose to `7,544`, `7,148`, and `7,691`.

The scalar bests and strict joint screens were essentially unchanged. Two
additional decomposition screens are more informative:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `young_old_roomgap` | 43.148 | 1.758 | 0.275 | 0.000 | 0.730 | 1.553 | 4.812 | 6.366 | 2.288 |
| Young ownership + old exit + fertility, no room-gap restriction | `old_retention` | 43.466 | 1.605 | 0.284 | 0.233 | 0.813 | 0.587 | 5.516 | 6.103 | 0.458 |

This separates the failure into two corners. The model can combine fertility,
low old ownership, and some owner-renter room separation, but only by nearly
emptying young ownership. It can also combine young ownership, low old
ownership, and fertility, but the owner-renter room gap collapses and old
nonhousing wealth is too low. The missing object is still the joint
young-access/old-exit/space-separation allocation, not a single isolated
moment.

### 2026-06-19 11:33 EDT Decomposition Frontier Tightens

Wave 1 was still in its second batch. Wave 2 and Wave 3 remained pending on
the intended dependencies, and no nonzero stderr files were visible. The scalar
best points did not move, but one decomposition frontier improved sharply.

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention` | 25.588 | 1.815 | 0.279 | 0.000 | 0.779 | 1.757 | 4.739 | 6.496 | 1.601 |
| Young ownership + old exit, no room-gap restriction | `young_old_own` | 41.633 | 1.983 | 0.245 | 0.462 | 0.771 | -0.221 | 5.883 | 5.662 | 2.516 |

The first row is a real improvement relative to the previous decomposition
screen: the model can now jointly hit plausible fertility, low old ownership,
and meaningful room separation. But it does so with zero young ownership. The
second row shows the opposite corner: young ownership and low old ownership can
coexist, but fertility is high and the owner-renter room gap turns negative.
This is stronger evidence that the missing margin is not merely old exit or
space separation in isolation; it is preserving the young-owner pipeline while
maintaining owner-renter space separation.

## Failure Modes To Track

1. Old ownership remains too high even when \(\theta_0\) is low or old
   ownership is directly targeted.
2. Young ownership and old ownership move in opposite directions through
   tenure/access parameters, making a lifecycle slope alone too coarse.
3. Owner-renter room separation remains too small: owners are too small and
   renters too large relative to ACS room targets.
4. Parent-childless old wealth gaps remain below target, so \(\theta_n\) is
   still weakly identified by current old wealth moments.

### 2026-06-19 11:40 EDT Three-Day Campaign Commitment

The user authorized a three-day unattended diagnostic campaign with hourly
check-ins. The rules are unchanged: do not alter main model/solver code, do
not delete outputs, and do not weaken the target system without recording the
replacement moment or external restriction that preserves identification.

Wave 1 remained healthy. Tasks `1--8` for all three arrays had completed with
exit code `0:0`; tasks `9--16` were running. Wave 2 broad and Wave 3 tight
were still pending on their intended dependencies. No nonzero stderr files were
visible.

Best scalar points:

| Run | Cases | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 8,265 | 14.778 | 1.786 | 0.254 | 0.023 | 0.936 | 1.105 | 4.515 | 5.620 | 2.517 |
| `young_old_own` | 7,923 | 22.653 | 1.835 | 0.269 | 0.136 | 0.900 | 0.538 | 4.957 | 5.495 | 2.517 |
| `young_old_roomgap` | 8,487 | 31.956 | 1.598 | 0.316 | 0.004 | 0.952 | 1.760 | 3.949 | 5.709 | 1.830 |

Best diagnostic decomposition points:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention` | 25.588 | 1.815 | 0.279 | 0.000 | 0.779 | 1.757 | 4.739 | 6.496 | 1.601 |
| Young ownership + old exit + fertility, no room-gap restriction | `young_old_own` | 41.633 | 1.983 | 0.245 | 0.462 | 0.771 | -0.221 | 5.883 | 5.662 | 2.517 |

Current interpretation: the scalar best is not the economic best. The best
economic evidence is the two-corner decomposition above. The model can combine
fertility, low old ownership, and room separation only by emptying young
ownership. It can combine young ownership and low old ownership only by
collapsing the owner-renter room gap. The main campaign question is whether
Wave 2/3 can find the missing joint allocation or whether this is a structural
mechanism/target mismatch.

### 2026-06-19 11:48 EDT No-Young Frontier Improves

Wave 1 remained healthy in its second batch. Wave 2 broad and Wave 3 tight were
still pending on the intended dependencies. No nonzero stderr files or Slurm
failures were visible.

Finite records reached `26,212` across the three Wave 1 target sets. Scalar
bests remained unchanged:

| Run | Cases | Loss | TFR | Childless | Own 25--34 | Old own | Room gap |
|---|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 8,754 | 14.778 | 1.786 | 0.254 | 0.023 | 0.936 | 1.105 |
| `young_old_own` | 8,435 | 22.653 | 1.835 | 0.269 | 0.136 | 0.900 | 0.538 |
| `young_old_roomgap` | 9,023 | 31.956 | 1.598 | 0.316 | 0.004 | 0.952 | 1.760 |

The useful frontier movement was in the decomposition screen that ignores young
ownership:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention` | 18.472 | 1.798 | 0.286 | 0.000 | 0.831 | 1.904 | 4.682 | 6.586 | 1.830 |

This is a genuine improvement over the previous no-young frontier loss of
`25.588`: the model now finds plausible fertility, lower old ownership, and a
larger owner-renter room gap together. But young ownership remains effectively
zero, so the central joint failure is unchanged. The missing allocation is
still young access plus old exit plus owner-renter space separation.

### 2026-06-19 11:53 EDT Room-Gap Scalar Improves, Joint Failure Persists

Wave 1 was still healthy in its second batch, with Wave 2 and Wave 3 pending on
the intended dependencies. No stderr or Slurm failures were visible.

Finite records reached `27,160` across Wave 1. The `young_old_roomgap` scalar
best improved from loss `31.956` to `28.854`:

| Run | Cases | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `young_old_roomgap` | 9,348 | 28.854 | 1.686 | 0.280 | 0.037 | 0.963 | 1.711 | 4.175 | 5.886 | 1.830 |

This is a scalar improvement, not a solution to the main mechanism failure.
The improved point has plausible fertility and a moderate room gap, but old
ownership is still too high and young ownership is still nearly empty. The
fixed screens are unchanged: no `joint_core` candidate, no fertility-ok soft
joint candidate, and the best young-ownership/old-exit candidate still has a
negative owner-renter room gap.

### 2026-06-19 11:59 EDT Young/Old/Fertility Frontier Improves

Wave 1 remained healthy. The first eight tasks in each array had completed
with exit code `0:0`; the remaining Wave 1 tasks were running. Wave 2 broad
and Wave 3 tight remained pending on the intended dependencies. No nonzero
stderr files were visible.

Finite records reached `28,499` across Wave 1. Scalar bests were essentially
unchanged:

| Run | Cases | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention` | 9,489 | 14.778 | 1.786 | 0.254 | 0.023 | 0.936 | 1.105 | 4.515 | 5.620 | 2.517 |
| `young_old_own` | 9,217 | 22.653 | 1.835 | 0.269 | 0.136 | 0.900 | 0.538 | 4.957 | 5.495 | 2.517 |
| `young_old_roomgap` | 9,793 | 28.854 | 1.686 | 0.280 | 0.037 | 0.963 | 1.711 | 4.175 | 5.886 | 1.830 |

The useful movement was in the decomposition screen that keeps young ownership,
old exit, and fertility in range while ignoring the room-gap restriction:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Young ownership + old exit + fertility, no room-gap restriction | `young_old_own` | 35.068 | 1.695 | 0.260 | 0.228 | 0.727 | 0.417 | 5.657 | 6.074 | 2.745 |

This improves the previous best in that screen, but it does not solve the
joint frontier. The candidate preserves fertility and lowers old ownership
while keeping some young ownership, but renters are too large and the
owner-renter room gap remains far below the ACS gap. The missing object remains
the joint allocation with young access, old exit, and owner-renter space
separation.

### 2026-06-19 15:21 EDT Wave 1 Complete, Wave 2 Broad Running

Wave 1 completed cleanly: all `48` array tasks across the old-retention,
young/old-ownership, and young/old/room-gap target sets exited with code `0:0`.
Wave 2 broad was running its second batch of tasks, with tasks `17--24` pending
behind the array throttle. Wave 3 tight remained pending on the intended Wave 2
dependency. No nonzero stderr files were visible.

The complete Wave 1 record has `33,600` finite cases. Wave 2 already had
`23,046` partial finite cases. The best scalar point moved from Wave 1's
`14.778` to a partial Wave 2 loss of `14.197`, but the economic miss did not
change:

| Run | Cases | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w2` | 7,707 | 14.197 | 1.643 | 0.357 | 0.029 | 0.922 | 1.097 | 4.650 | 5.746 | 1.830 | 0.833 |

This is a scalar improvement, not a better calibration. It lowers the loss by
raising the parent-childless old nonhousing wealth gap and keeping second-child
housing growth close to target, but it leaves young ownership almost empty,
old ownership too high, childlessness too high, and the owner-renter room gap
far below the ACS target.

The constrained frontier remained split:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Joint core: old `<=0.85`, young `>=0.25`, gap `>=1.5` | none | -- | -- | -- | -- | -- | -- | -- | -- |
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention_w1` | 18.472 | 1.798 | 0.286 | 0.000 | 0.831 | 1.904 | 4.682 | 6.586 |
| Young ownership + old exit + fertility, no room-gap restriction | `young_old_own_w1` | 35.068 | 1.695 | 0.260 | 0.228 | 0.727 | 0.417 | 5.657 | 6.074 |
| Wave 2 young + old, no fertility/gap restriction | `old_retention_w2` | 38.596 | 1.485 | 0.307 | 0.285 | 0.836 | 0.426 | 5.394 | 5.820 |

Current interpretation: the broader Wave 2 search is not merely stuck at the
Wave 1 scalar point, but the frontier geometry is unchanged. Points that
preserve old exit and room separation empty the young-owner pipeline; points
with young ownership and old exit produce large renters and little or no
owner-renter room separation. Let Wave 2 finish and then compare it with Wave 3
tight before launching a new experimental target system.

### 2026-06-19 16:22 EDT Wave 2 Soft-Joint Frontier Moves

Wave 2 broad remained healthy. Tasks `1--8` had completed with exit code
`0:0`; tasks `9--16` were running; tasks `17--24` were pending behind the
array throttle. Wave 3 tight remained pending on the intended Wave 2
dependency. No nonzero stderr files were visible.

The aggregate record reached `67,688` finite cases: `33,600` from completed
Wave 1 and `34,088` partial cases from Wave 2. The scalar best remained
Wave 2's `old_retention` point at loss `14.197`, still with almost no young
ownership:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w2` | 14.197 | 1.643 | 0.357 | 0.029 | 0.922 | 1.097 | 4.650 | 5.746 | 1.830 | 0.833 |

Wave 2 did, however, improve the relaxed joint frontier. The best candidate
with old ownership below `0.88`, young ownership above `0.20`, and an
owner-renter room gap above `1.0` moved from Wave 1's loss `71.021` to
Wave 2's loss `51.180`:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `young_old_own_w2` | 51.180 | 2.559 | 0.100 | 0.201 | 0.852 | 1.019 | 5.194 | 6.213 | 3.432 | -0.449 |

This is useful frontier information, not a satisfactory fit. The candidate is
just inside the soft ownership/room screen, but it generates too many births,
too little childlessness, renters that are far too large, and a wrong-signed
old parent-childless nonhousing-wealth gap. The stricter joint screen remains
empty:

\[
\text{old ownership}\le 0.85,\qquad
\text{young ownership}\ge 0.25,\qquad
\text{owner-renter room gap}\ge 1.5.
\]

The earlier decomposition also remains: fertility plus old exit plus room
separation still comes with essentially zero young ownership, while young
ownership plus old exit and fertility still collapses the room gap.

### 2026-06-19 19:25 EDT Wave 2 Complete, Wave 3 Started

Wave 2 broad completed cleanly. All `72` Wave 2 array tasks exited with code
`0:0`, and no nonzero stderr files were visible. Wave 3 tight started its
first batch of tasks across the same three identified target sets.

The complete Wave 2 record contains `62,270` finite cases. It did not solve the
joint allocation problem. The best scalar point remained the old-retention
Wave 2 candidate:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w2` | 14.197 | 1.643 | 0.357 | 0.029 | 0.922 | 1.097 | 4.650 | 5.746 | 1.830 | 0.833 |

The strict joint screen remained empty even after the complete broad search:

\[
\text{old ownership}\le 0.85,\qquad
\text{young ownership}\ge 0.25,\qquad
\text{owner-renter room gap}\ge 1.5.
\]

The complete Wave 2 decomposition confirms the same two corners:

| Screen | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH gap |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention_w2` | 27.898 | 1.684 | 0.300 | 0.000 | 0.777 | 1.652 | 4.955 | 6.607 | 0.401 |
| Young ownership + old exit + fertility, no room-gap restriction | `old_retention_w2` | 36.703 | 1.779 | 0.207 | 0.213 | 0.829 | 0.526 | 5.563 | 6.089 | -0.207 |
| Young ownership + old exit, no fertility/gap restriction | `old_retention_w2` | 38.596 | 1.485 | 0.307 | 0.285 | 0.836 | 0.426 | 5.394 | 5.820 | -0.134 |

The first Wave 3 partial read had only `6,078` finite cases, so it is too early
to interpret tightly. Its best scalar point was `old_retention_w3` with loss
`19.319`, TFR `1.905`, young ownership `0.010`, old ownership `0.877`, and
room gap `1.191`. It did not yet improve the strict or soft joint screens.

Current interpretation: Wave 2's broad search materially increased coverage,
but it reinforced the same mechanism failure rather than finding a hidden
joint basin. The model finds old exit and room separation only by emptying
young ownership, and it finds young ownership plus old exit only with large
renters and a small or negative owner-renter room gap. Let Wave 3 tight finish
before deciding whether the next step is a new identified 13-moment
experimental target system or a written mechanism diagnosis.

### 2026-06-19 20:25 EDT Wave 3 Scalar Best Improves, Frontier Unchanged

Wave 3 tight was running normally in its first batch. Tasks `1--8` were active
for each of the three target sets, with tasks `9--24` pending behind the array
throttle. No nonzero stderr files were visible.

The partial Wave 3 record reached `16,828` finite cases. It produced a new
global scalar best, but not an economically better joint allocation:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w3` | 12.795 | 1.744 | 0.282 | 0.019 | 0.919 | 1.391 | 4.498 | 5.889 | 2.288 | 0.381 |

This point improves the scalar loss because it balances several targeted
margins better than the prior scalar best, including old median nonhousing
wealth and the room gap. But it leaves young ownership essentially empty, old
ownership too high, childlessness too high, and the owner-renter room gap well
below the ACS target.

The strict joint screen is still empty:

\[
\text{old ownership}\le 0.85,\qquad
\text{young ownership}\ge 0.25,\qquad
\text{owner-renter room gap}\ge 1.5.
\]

The partial Wave 3 decompositions continue to show the same tradeoff:

| Screen | Best Wave 3 run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH gap |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention_w3` | 22.189 | 1.873 | 0.255 | 0.001 | 0.846 | 1.798 | 4.493 | 6.290 | 0.362 |
| Young ownership + old exit + fertility, no room-gap restriction | `old_retention_w3` | 42.748 | 1.780 | 0.241 | 0.546 | 0.827 | -0.052 | 5.404 | 5.352 | 0.121 |
| Young ownership + old exit, no fertility/gap restriction | `old_retention_w3` | 28.682 | 1.512 | 0.359 | 0.543 | 0.848 | 0.174 | 5.294 | 5.468 | 0.435 |

Interpretation: the tight search can lower the scalar objective, but so far it
does so by staying in the same corner: poor young access, high old retention,
and insufficient owner-renter space separation. This is useful evidence
against a pure search-depth explanation. Continue Wave 3 before launching a
new target variant.

### 2026-06-19 21:27 EDT Wave 3 Continues; Scalar Best Moves Again

Wave 3 remains healthy. Tasks `1--8` have completed for each of the three tight
target sets, tasks `9--16` are running, and tasks `17--24` are queued behind the
array throttle. Completed tasks report exit code `0:0`, and no nonzero stderr
files are visible.

The partial Wave 3 record now contains `28,041` finite cases. The scalar best
moved again:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w3` | 12.197 | 1.846 | 0.268 | 0.030 | 0.898 | 1.349 | 4.451 | 5.800 | 2.288 | 0.362 |

This improves the scalar score but still fails the economic joint allocation
screen. Young ownership remains essentially absent, old-age ownership remains
too high, and the owner-renter room gap is below the target.

The strict joint screen remains empty:

\[
\text{old ownership}\le 0.85,\qquad
\text{young ownership}\ge 0.25,\qquad
\text{owner-renter room gap}\ge 1.5.
\]

The current Wave 3 frontier still splits into the same two corners:

| Screen | Best Wave 3 run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH gap |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention_w3` | 18.086 | 1.828 | 0.259 | 0.000 | 0.849 | 1.582 | 4.760 | 6.342 | 0.205 |
| Young ownership + old exit + fertility, no room-gap restriction | `old_retention_w3` | 38.014 | 1.705 | 0.264 | 0.346 | 0.824 | 0.442 | 5.035 | 5.477 | -0.362 |
| Young ownership + old exit, no fertility/gap restriction | `old_retention_w3` | 28.682 | 1.512 | 0.359 | 0.543 | 0.848 | 0.174 | 5.294 | 5.468 | 0.435 |

Interpretation: Wave 3 is usefully lowering the objective, so it should be
left running. But the new best point does not change the mechanism diagnosis:
the current target geometry still has not found a candidate with old exit,
young ownership, and owner-renter room separation at the same time.

### 2026-06-19 22:25 EDT Wave 3 Partial Frontier Softens, Not Solved

Wave 3 is still running normally. Tasks `9--16` are active at roughly
`1:37` elapsed, tasks `17--24` are still queued behind the array throttle, and
no nonzero stderr files are visible.

The partial record has grown to `38,813` finite Wave 3 cases and `134,683`
finite cases across Waves 1--3. The scalar best improved again:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w3` | 11.717 | 1.908 | 0.227 | 0.004 | 0.912 | 1.363 | 4.482 | 5.845 | 2.517 | 0.327 |

The new scalar best is numerically close to the earlier reference diagnostic
loss range, but it is still an ownership-corner point: young ownership is nearly
zero and old-age ownership remains high.

One new soft joint candidate appears under the looser screen
\[
1.55\le \mathrm{TFR}\le 1.9,\quad
\mathrm{childless}\le 0.28,\quad
\mathrm{old\ own}\le 0.88,\quad
\mathrm{own}_{25--34}\ge 0.20,\quad
\mathrm{room\ gap}\ge 1.0.
\]

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w3` | 24.434 | 1.750 | 0.247 | 0.214 | 0.869 | 1.036 | 4.482 | 5.518 | 2.974 | -0.008 |

This is progress relative to the earlier soft frontier, but it is not a
solution. The strict joint screen remains empty:

\[
\text{old ownership}\le 0.85,\qquad
\text{young ownership}\ge 0.25,\qquad
\text{owner-renter room gap}\ge 1.5.
\]

The main tradeoff remains:

| Screen | Best Wave 3 run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH gap |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention_w3` | 15.014 | 1.697 | 0.280 | 0.000 | 0.842 | 1.656 | 4.447 | 6.102 | 0.107 |
| Young ownership + old exit + fertility, no room-gap restriction | `young_old_own_w3` | 26.971 | 1.630 | 0.288 | 0.207 | 0.838 | 0.796 | 4.793 | 5.589 | 0.031 |
| Young ownership + old exit, no fertility/gap restriction | `old_retention_w3` | 25.450 | 1.947 | 0.213 | 0.326 | 0.849 | 0.396 | 5.253 | 5.650 | 0.082 |

Interpretation: Wave 3 is not merely repeating the first batch; it has found a
softer compromise. But the compromise still gets there by compressing the
owner-renter room gap and producing an almost zero parent-childless old wealth
gap. Continue the current arrays; do not launch a new diagnostic wave while
tasks are still active.

### 2026-06-19 23:25 EDT Wave 3 Final Batch Running; Strict Screen Still Empty

Wave 3 is in its final batch. Tasks `1--16` have completed for each array with
exit code `0:0`; tasks `17--24` are running at roughly `42` minutes elapsed.
No nonzero stderr files are visible.

The partial record has reached `49,559` finite Wave 3 cases and `145,429`
finite cases across Waves 1--3. The scalar best is now close to the earlier
reference diagnostic range:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w3` | 11.578 | 1.736 | 0.287 | 0.006 | 0.908 | 1.389 | 4.430 | 5.819 | 2.059 | 0.168 |

This is a better scalar score, but it remains an ownership-corner point: young
ownership is essentially zero, old-age ownership is too high, and the
owner-renter room gap is still below the empirical gap.

The strict joint screen remains empty:

\[
\text{old ownership}\le 0.85,\qquad
\text{young ownership}\ge 0.25,\qquad
\text{owner-renter room gap}\ge 1.5.
\]

The soft joint candidate from the previous pulse still exists, but no better
soft candidate has replaced it. It has TFR `1.750`, childlessness `0.247`, young
ownership `0.214`, old ownership `0.869`, room gap `1.036`, and old
parent-childless nonhousing wealth gap `-0.008`.

The partial decomposition remains:

| Screen | Best Wave 3 run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH gap |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fertility + old exit + room gap, no young-ownership restriction | `old_retention_w3` | 15.014 | 1.697 | 0.280 | 0.000 | 0.842 | 1.656 | 4.447 | 6.102 | 0.107 |
| Young ownership + old exit + fertility, no room-gap restriction | `young_old_own_w3` | 26.971 | 1.630 | 0.288 | 0.207 | 0.838 | 0.796 | 4.793 | 5.589 | 0.031 |
| Young ownership + old exit, no fertility/gap restriction | `old_retention_w3` | 25.450 | 1.947 | 0.213 | 0.326 | 0.849 | 0.396 | 5.253 | 5.650 | 0.082 |

Interpretation: the search is still productive in scalar terms, but it is not
yet weakening the core mechanism failure. The current model/target geometry can
approach the scalar benchmark by giving up young ownership, or it can get young
ownership plus old exit by compressing the owner-renter room gap. Let the final
Wave 3 batch complete before designing a next bounded diagnostic wave.

### 2026-06-21 Final Wave 1--3 Readout After Auth Recovery

Torch access was restored on June 21. All nine Wave 1--3 arrays completed with
Slurm exit code `0:0`, and all visible stderr files for jobs `11136888`,
`11136889`, `11136890`, `11146362`, `11146363`, `11146364`, `11146365`,
`11146366`, and `11146367` were zero bytes. The completed campaign produced
`158,405` finite candidate records.

Finite cases by run:

| Run | Cases |
|---|---:|
| `old_retention_w1` | 11,200 |
| `young_old_own_w1` | 11,200 |
| `young_old_roomgap_w1` | 11,200 |
| `old_retention_w2` | 20,662 |
| `young_old_own_w2` | 20,348 |
| `young_old_roomgap_w2` | 21,260 |
| `old_retention_w3` | 20,979 |
| `young_old_own_w3` | 20,624 |
| `young_old_roomgap_w3` | 20,932 |

The final scalar best is the same point seen in the partial final-batch pulse:

| Run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `old_retention_w3` | 11.578 | 1.736 | 0.287 | 0.006 | 0.908 | 1.389 | 4.430 | 5.819 | 2.059 | 0.168 |

This is not an economic solution. It gets a good scalar score by nearly
emptying the young-owner pipeline. The completed strict joint screen remains
empty:

\[
\text{own}_{25--34}\ge 0.25,\qquad
\text{old ownership}\le 0.85,\qquad
\text{owner-renter room gap}\ge 1.5.
\]

The best completed frontier screens are:

| Screen | Count | Best run | Loss | TFR | Childless | Own 25--34 | Old own | Room gap | Renter rooms | Owner rooms | Old NH med | Old NH gap |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Old ownership `<=0.85` | 42,564 | `old_retention_w3` | 14.135 | 1.688 | 0.285 | 0.055 | 0.845 | 1.309 | 4.622 | 5.931 | 2.288 | -0.072 |
| Young `>=0.25`, old `<=0.85` | 509 | `old_retention_w3` | 25.450 | 1.947 | 0.213 | 0.326 | 0.849 | 0.396 | 5.253 | 5.650 | 2.745 | 0.082 |
| Soft joint: young `>=0.20`, old `<=0.88`, gap `>=1.0` | 8 | `old_retention_w3` | 24.434 | 1.750 | 0.247 | 0.214 | 0.869 | 1.036 | 4.482 | 5.518 | 2.974 | -0.008 |
| Fertility-ok soft joint | 3 | `old_retention_w3` | 24.434 | 1.750 | 0.247 | 0.214 | 0.869 | 1.036 | 4.482 | 5.518 | 2.974 | -0.008 |
| Fertility + old exit + room gap, no young restriction | 1,794 | `old_retention_w3` | 14.609 | 1.880 | 0.241 | 0.001 | 0.874 | 1.579 | 4.485 | 6.064 | 2.517 | 0.198 |
| Young ownership + old exit + fertility, no room-gap restriction | 297 | `young_old_own_w3` | 16.522 | 1.888 | 0.245 | 0.223 | 0.867 | 0.807 | 4.847 | 5.654 | 2.288 | 0.299 |
| Young ownership + room gap + fertility, no old restriction | 645 | `young_old_own_w3` | 24.339 | 1.693 | 0.271 | 0.209 | 0.946 | 1.045 | 4.289 | 5.334 | 2.059 | 0.305 |
| Room gap `>=2` | 4,576 | `old_retention_w3` | 13.673 | 1.808 | 0.300 | 0.000 | 0.879 | 2.015 | 4.252 | 6.267 | 2.517 | 0.350 |

Final mechanism read: the completed Wave 1--3 campaign does not solve the
joint target system. The model can generate old exit plus room separation, but
then young ownership is essentially zero. It can generate young ownership plus
old exit and tolerable fertility, but then the owner-renter room gap collapses.
It can generate young ownership plus a positive room gap and tolerable
fertility, but then old ownership remains too high. The soft joint candidates
are real progress, but they still have a room gap near one room rather than the
ACS target gap and essentially no parent-childless old nonhousing wealth gap.

Do not interpret this as permission to drop targets. Any formal next target
revision must preserve identification: the current evidence says which
moments/blocks are in tension, not that they can be removed.

### 2026-06-21 Frontier Jacobian Audit Launched

After the completed Wave 1--3 frontier showed persistent three-way tension, a
bounded local identification audit was launched on Torch. This is not a new
calibration search: it holds the 13-parameter vector fixed by target-set point
and computes finite-difference derivatives of the active target moments around
selected frontier candidates.

The selected points are stored under
`/scratch/td2248/projects/Fertility_Spring26_20260617_fast/code/cluster/audit_points_20260621_frontier/`.
They cover the scalar old-retention best, the soft joint compromise, the old
exit plus room-gap/no-young-ownership corner, the young-ownership plus old-exit
/no-room-gap corner, the young-ownership plus room-gap/no-old-exit corner, and
the room-gap target set's scalar best.

Six one-point Jacobian jobs were submitted, each with `J=16`, `Nb=60`,
`income_states=5`, `n_house=5`, `max_iter_eq=3`, relative step `0.01`, and
27 solves per point (`1+2*13`). The total audit is 162 solves, parallelized
across independent `cpu_short` jobs:

| Job | Point | Target set |
|---:|---|---|
| `11384924` | `scalar_best_old_retention_w3` | `candidate_replacement_old_retention_v1` |
| `11384925` | `soft_joint_old_retention_w3` | `candidate_replacement_old_retention_v1` |
| `11384926` | `old_gap_no_young_old_retention_w3` | `candidate_replacement_old_retention_v1` |
| `11384927` | `young_old_no_gap_young_old_own_w3` | `candidate_replacement_young_old_own_v1` |
| `11384928` | `young_gap_no_old_young_old_own_w3` | `candidate_replacement_young_old_own_v1` |
| `11384929` | `roomgap_scalar_young_old_roomgap_w3` | `candidate_replacement_young_old_roomgap_v1` |

Initial queue check: all six jobs were running, and all six Slurm stderr files
were zero bytes. Outputs will be written under
`/scratch/td2248/projects/Fertility_Spring26_20260617_fast/output/model/intergen_frontier_jacobian_20260621_*`.
