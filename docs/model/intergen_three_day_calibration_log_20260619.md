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

## Failure Modes To Track

1. Old ownership remains too high even when \(\theta_0\) is low or old
   ownership is directly targeted.
2. Young ownership and old ownership move in opposite directions through
   tenure/access parameters, making a lifecycle slope alone too coarse.
3. Owner-renter room separation remains too small: owners are too small and
   renters too large relative to ACS room targets.
4. Parent-childless old wealth gaps remain below target, so \(\theta_n\) is
   still weakly identified by current old wealth moments.
