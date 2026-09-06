# Evening housing-response follow-up

Author renewed authorization at 22:07 UTC September 5 and will return no earlier
than 21:00 Eastern (01:00 UTC September 6). The question is whether the first-birth
housing response can improve with the same twelve targets, eleven parameters,
original bounds and numerical gates. Production and the completed candidate
policy comparison are retained. A new bounded design follows inspection of the
recent panel and historical decisions. Exact smoke array `17010311` completed
two identical full-history evaluations. The 23-case coordinate array `17010431`, twelve joint cases `17014282`, and two final exact repeats `17015163` completed with all gates passing. The unchanged-bound winner is evening `round2/task_010`, loss 21.7980329171, first-birth rooms 0.4766454140 versus target 0.7202462624, and mean rooms 6.3325296653 versus 5.7799704819. Complete fit and parameter tables are in that case folder. Production remains the separate September 4 task_010.

## Predeclared experiment and budget

Use the verified repaired scientific snapshot
`/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/`, bundle
`33167d84113e2bd38d9ee48dcd9ab0403790348610d998d4032fb8c1797ad3e3`.
The existing observation adapter now accepts only the original and this verified
bundle, and pins the exact chosen bundle in every launch and validation. Five
pure guard tests pass. No scientific model file changes.

1. Two exact complete-history repeats of the reviewed morning candidate using
   its original generator. All twelve fit rows, parameters and 253 numeric
   historical entries must reproduce, with checkpoints and seventeen graphs.
2. After smoke collection passes, 23 coordinate cases at radius 0.0025 in the
   existing transformed bounded coordinates, around the reproduced candidate.
   This is a new center and half the morning coordinate radius. All eleven
   parameters and the twelve-row objective are retained; no Jacobian inversion.
3. After all 23 pass, at most twelve predeclared joint cases from the best point.
   Nine combinations set the first-child intercept to its existing maximum 0.5,
   either preserving the three-child floor, preserving the two-child floor, or
   retaining the per-child slope, each at H0 multipliers 1, 0.97, and 0.85.
   Three further cases combine the three-child-preserving endpoint with the
   individually improving non-housing coordinate moves at half, one, and two
   step lengths. Exclude H0 and both housing floors from this latter direction;
   clip only to the unchanged search bounds, and omit duplicate points.
   For m>0, floor(m)=jump+m*slope; preserving the m-child floor therefore uses
   slope_new=slope_old-(0.5-jump_old)/m. Reject any floor pattern outside bounds.
4. Two exact repeats of the cross-stage lowest-loss point, using its original
   generator. If all new points worsen fit, retain the starting candidate.

Maximum 39 complete historical evaluations, one CPU/8GB each, six concurrent,
30 minutes hard limit per case. Recent full historical evaluations took about
7--11 minutes: approximately 4.6--7.2 CPU-hours and 50--80 minutes of numerical
work at six-way concurrency, plus queues, preparation and reporting. Each case
contains old-equilibrium normalization and five historical market dates. The
reported dated-model solve count and old-normalization stationary solve count
are kept separately; no complete count of inner Bellman iterations is claimed.

New calibration launches stop at 00:10 UTC September 6. Reserve through 01:00 UTC
for final collection/readout. Do not launch a dependent stage after a source,
target, accounting, market, population, probability or feasibility gate failure.
No automatic numerical retry. One-minute heartbeats and completed-case artifacts;
investigate thirty minutes without progress. Keep latest-completed and best-so-far
summaries during each stage by partial collection. Stop after these two rounds.
This is a same-target candidate search, not a claim of a globally optimal fit or
unreachable empirical target. Production stays retained and no policy runs are
included in this predeclared calibration budget.

The existing 0.5 bound is held fixed pending the historical evidence review.
No uncertainty interval or robust policy estimate follows from this small panel.

## Historical reconciliation checked by the lead

The bounded worker report is preserved as `history_worker_report.md`. Its E5b
analogy needs a precise qualification: the August 9 transcript at lines 580--655
uses first-child and additional-child expenditure-share tilts (`delta_alpha_jump`,
`delta_alpha`), not today's housing-floor intercept and slope. Its older target
and fit are not comparable to this experiment.

Git commit `f6ab6704a389789a05a0281a14371b648ffbb0eb` on August 18 introduces
`hbar_first_child_jump` and its estimated [0, 0.5] domain in the transition
driver. The contemporaneous status, first 140 lines at that commit, describes
a fixed diagnostic grid 0, .1, .2, .25, .35, .5 with a best point at .2, and a
joint estimate .199446. This explains why the bound was not then binding; it
does not establish an empirical or theoretical justification for 0.5. The
August 18 transcript manifest contains no sessions. No direct author choice of
that endpoint was recovered in the bounded search. Earlier assertions of
misspecification are prior assessments, not proofs of an unreachable target.

The active solver constructs floor(m)=jump+m*slope for m>0. The lead checked
this in `precompute_shared`. The dated target compares two branches from the
same successful-first-birth risk set, permits continuation births only in the
treated branch and holds the control childless. `analyze.py` reconciles each
saved response to treated minus control housing. Its cross-candidate level
decomposition is neither causal nor fixed-composition. The supplemental fit
scatterplot leaves the seventeen standard diagnostic graphs unchanged.

```sh
python3 output/model/e5f_evening_housing_refinement_20260905a/analyze.py
```

To regenerate the unchanged diagnostic packet without solving, run inside the
frozen Torch project (substitute the chosen stage/case and a new output folder):

```sh
python3 output/model/e5f_evening_housing_refinement_20260905a/analyze.py --graphs smoke/task_001 --graphs-out output/model/e5f_evening_housing_refinement_20260905a/graph_replay_smoke
```

The demonstrated replay reproduced all seventeen original PNG hashes exactly,
with zero model solves. It verifies the saved case, checkpoint, source bundle
and graph-helper hash before reading the packet. The example destination already
exists; choose another new directory to repeat it.

## Contingent diagnostic after the exact bounded repeats

Only if the complete, twice-reproduced bounded winner has jump=0.5, a separate
conditional profile can distinguish an inherited search ceiling from a lack
of housing response in the model. This is outside the 39-case search budget
and is not a recalibration or a changed production bound. The explicit fixed-
jump option has existed in the scientific driver since August 18.

The prepared `fixed_profile.py` requires the bounded final repeats first. It
would run two exact fixed-option smokes at the same physical winner, three
points at jump .6, .75 and min(.9, jump_old+3*(slope_old-.1)), preserving the
three-child floor, then two exact repetitions of the best conditional point.
Apart from this compensating per-child slope, the other nine estimated
coordinates remain at the bounded winner and are checked in collection; the old fertility
intercept is re-normalized by the unchanged driver. All twelve targets and
weights, scientific source, markets, population and numerical gates stay fixed.
No identifying row is dropped. The diagnostic jump is explicitly external;
the remaining ten-coordinate domain is not re-estimated in this profile.

Maximum seven additional histories, at most three concurrent, one CPU/8GB and
30 minutes hard limit each. The complete coordinate panel took about 8--12
minutes per case. At 23:18 UTC, before inspecting any joint result or making a conditional launch,
the contingency budget was revised prospectively to 0.9--1.4 CPU-hours and
24--36 minutes for three stages, plus queues. Its launch cutoff is 00:20 UTC;
do not start this contingent follow-up after 23:50 UTC. The original 39-case
search keeps its separate 00:10 cutoff. This preserves the seven-history cap
and reserves the remainder of the author's pre-21:00 window for review.
Stop dependent stages on any
failure; no retry. Each case has the same heartbeat, checkpoint and seventeen
graphs, with separate profile latest-completed and best-so-far summaries.

`fixed_jump_adapter.py` is an output-owned copy of the observation adapter, with
only explicit fixed-option routing, ten-coordinate classification and profile
guards changed. The lead inspected the complete diff. Eight pure plan guards
and the shell syntax check pass; the two-case fixed-option smoke completed as array `17015491`. Its plan pins the unchanged scientific bundle and target fingerprint. Both original bounded final repeats reproduced all twelve fit rows, parameters and 253 numeric history entries exactly.
The existing search adapter and scientific code are untouched by this preparation.

## Targeted value-screen check

Coordinate case 8 (lower owner-service premium) has one occupied adjacent-wealth
value drop at age 30, renter, income state 12, parity/dependent children 3/3,
wealth -0.81395 to -0.67442. The drop is 0.0044611; exposed pre-choice mass share
is 2.28775e-7. This point has worse fit and is not the provisional best.

`inspect_saving_flag.py` reuses the already tested paired local/exhaustive-saving
loop, with two fixed-price Bellman evaluations at this single saved 2023 state.
This is additional numerical diagnosis, not another historical calibration case.
First require exact local policy/distribution reproduction from the saved
pre-choice distribution; only then evaluate exhaustive saving, require occupied
value dominance, and report changes in births, rooms and ownership. Inspect the
full source diff against the prior six-solve checker. Stop on failure, no retry.
Budget one CPU/8GB, expected 2--4 minutes, hard 10 minutes, no following market
run. Both evaluations write the seventeen graphs and a checkpoint, with latest-
completed and diagnostic best-so-far summaries. This cannot establish an error
bound for the historical SMM fit or the matched first-birth housing response.

## Local sensitivity and quantile reconstruction

The complete 23-case panel supports direct evaluation but gives a fragile
finite-difference Jacobian: singular values range from 726.413 to 0.544012,
condition number 1,335.29. Algebraic rank is eleven; ten singular values exceed
the predeclared relative 1e-3 screen. The weakest right direction is dominated
by the bequest shifter theta1. This is a local conditioning diagnostic under the
chosen weights, transforms and step, not global identification or a proof of
underidentification. No row or parameter is dropped and no Jacobian is inverted.

Upward and downward derivatives disagree principally because of the living-old
wealth-to-income p90/p50 statistic. The floor columns have one-sided cosine above .999.
The model uses the inverse weighted empirical CDF of its discrete state support
(`utils.py`, `weighted_quantile`), with pre-transaction asset mass and gross
housing value. A nine-checkpoint read-only reconstruction via `analyze.py
--quantile-check` reproduces every stored ratio using the unchanged definition.
For example, a small discount-factor reduction shifts p90 from 27.89536 to
27.58747 while p50 stays near 8.443: the ratio moves from 3.30359 to 3.26755.
The saved CDF brackets document the percentile's change of support value.
This identifies a source of nonsmooth objective responses; it is not a coding
failure or a reason to silently smooth, drop or reweight the empirical row.

The paired value check completed as job `17013496` in 34.76 seconds. Local
arrays and distributions reproduce exactly; exhaustive saving removes the sole
occupied value drop and weakly dominates the saved occupied values (minimum
difference -4.44e-16). Fixed-state differences are +0.000233535% births,
+0.0000691684% rooms and -0.0000215806 ownership percentage points. These are
conditional quantity differences, not historical SMM or event-response error bounds.

## Recovering the exact scientific source

All 51 fingerprinted scientific files are byte-identical to Git commit `c63821a6e027e90c85db1f8c137e617d4f7b5492`; `source_recovery.json` maps every bundle key to its repository path and SHA-256. The local input `output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json` matches the frozen remote input hash. This permits recovery independently of the temporary checkout. Observation adapters and plans are saved with this experiment and pin their own hashes.

After all completed histories finish, `analyze.py --verify-artifacts` validates every completed stage and original receipt on Torch and writes the full artifact manifest. After collection, `analyze.py --verify-collected` verifies the local copies and explicitly inventories checkpoints retained only on Torch. These commands do not solve the model.

Both fixed-option smokes (`17015491`) exactly reproduced all twelve fit rows, parameters and 253 numeric historical entries. Conditional array `17015738` evaluated jumps 0.6, 0.75, 0.9 with slopes 0.2038538348, 0.1538538348, 0.1038538348. Every point preserves the three-child floor 1.2115615043 and the other nine estimated coordinates; the original old-intercept normalization is rerun. These are completed conditional diagnostics, not revised production bounds.

## Completed conditional result and economical stop

All three extended points completed and passed all existing gates. Jumps 0.6, 0.75, 0.9 give losses 22.0891683, 30.6794904, 49.4962663 and first-birth responses 0.4995381, 0.5396778, 0.5800638. The three-plus-versus-one/two-child room gap deteriorates from 0.3727079 at the anchor to 0.3477174, 0.3104515, 0.2745712 against target 0.3676996. The parent ownership gap also widens too much. Full twelve-row fits and parameter ledgers are in `profile_grid/report/`; `analyze.py --profile-summary` reproduces the branch/floor ledger and supplemental three-panel plot. No occupied value drops occur at any of these points.

The previously twice-reproduced fixed-option anchor remains best. A final two-case repeat plan was prepared but deliberately not submitted: it would repeat exactly the same already verified point. `profile_stop.json` records this reduction from the seven-history cap to five completed histories; it does not change a numerical gate or an accepted result. The final manifest requires both exact fixed-option smoke receipts and verifies that none of the extended points improves fit. There is no new calibrated point outside the original bounds and no production promotion. The diagnostic does not exclude an improvement between 0.5 and 0.6 or under joint re-estimation.

## Final separate follow-up: preserve the observed room gap

At 00:26 UTC, after the original profile was closed, the lead declared one new, bounded hypothesis. Preserving the three-child preference floor compressed the *observed* 3+-versus-1/2-child room gap. The stable central floor columns instead imply physical derivatives -0.0190214 with respect to the first-child term and +0.6904211 with respect to the per-child slope at the morning center. This suggests evaluating a different direction; it is not a full Jacobian inversion or a prediction accepted as a solve.

`empirical_room_gap.py` declares four fixed first-child terms .525, .55, .575, .6, with slopes .2306218583, .2313106198, .2319993812, .2326881427 chosen to close the observed room-gap residual to first order. All other nine estimated coordinates stay at the bounded winner and are checked before collection; old fertility is re-normalized by the same scientific driver. All twelve actual target rows and weights remain active. The exact same fixed-option loop has already passed two complete-history smokes, including checkpoints, all seventeen original PNGs and 253 history entries. No scientific, adapter or numerical-gate code changes.

New maximum: four independent histories, at most four concurrent, followed by two exact repeats only for a newly improving point. At observed 9--10 minute solves this is about 0.9--1.0 CPU-hours and 20--30 minutes plus queue/reporting; one CPU/8GB, hard 30 minutes each. Preparation must begin before 00:35 UTC; all new launches stop at 00:55 UTC. No automatic retry or further stage. Stop on any source, target, market, population, accounting or feasibility failure. One-minute heartbeat, checkpoint and 17 graphs per case; separate latest/best summaries live in `empirical_room_gap/`. This new, six-history cap does not retroactively change the completed 39-case search or five-case floor-preserving diagnostic. Production remains retained.

The four follow-up solves completed 0:0. Initial collection stopped on a reporting KeyError: the scientific driver stores the externally fixed first-child term in `summary.best_candidate.theta`, not in the ten-row free-parameter table. `collect_empirical_room_gap.py` corrects only that lookup, retains all nine held-coordinate and both housing-coordinate checks, and calls the unchanged full collector. Original preparation, centers, numerical plans and all results remain unchanged. Its SHA is recorded in each collection receipt and pinned in any repeat plan. Use this corrected collector for the follow-up; do not rerun the original preparation script's collector action. No model solve was retried.

## Final verified result and collection

All numerical work is complete: 39 bounded histories, five original conditional
histories, six observed-room-gap histories, plus two fixed-price Bellman checks.
All jobs completed 0:0; allocated time totals 8.048611 CPU-hours. No further
stage is planned. The final diagnostic at `empirical_room_gap/grid/task_004`
has loss 20.6952742796, first-birth rooms 0.5157782240 versus 0.7202462624,
and family-size room gap 0.3675410210 versus 0.3676995588. Its fixed jump is
0.6, above the old 0.5 ceiling, and its slope is 0.2326881427. Mean rooms
remain high at 6.3491145755 versus 5.7799704819, and the parent ownership gap
is 0.1757242866 versus 0.16766167. All nine other coordinates remain at the
bounded winner. This is a conditional diagnostic, not a production promotion.

Repeat array `17019294` reproduces all twelve target-fit rows, every physical
parameter, 253 numeric historical entries and all seventeen standard PNGs
exactly in both cases. `empirical_room_gap/repeats/plan.json` pins its original
generator, scientific source, target fingerprint and corrected collector.
The complete fits, weights, losses and parameter restrictions are in each case
and the stage report ledgers. The selected point lies at the endpoint of the
small tested interval; no global frontier or unreachable-target claim follows.

`artifact_manifest.json` verifies 1,050 original artifact hashes across fifty
histories, with 884 standard graphs including the two fixed-price packets.
`collection_receipt.json` verifies 2,742 local files against that manifest;
52 large checkpoints remain on Torch. Saved-state graph replay reproduces all
seventeen PNGs exactly for the smoke, largest original extension, and final
selected diagnostic, without solving the model. The lead visually reviewed all
seventeen graphs for the bounded winner, the largest original extension, and
the selected observed-gap diagnostic. Old-age ownership saturation and uneven
tenure thresholds remain validation limitations.

The review is `docs/model/e5f_evening_housing_refinement_review.md`, rendered as
`output/pdf/e5f_evening_housing_refinement_review.pdf`. It includes full fit and
parameter tables and clearly labels earlier policy results at the morning
candidate. Rebuild with:

```sh
python3 code/model/tools/build_e5f_independent_audit_pdf.py --source docs/model/e5f_evening_housing_refinement_review.md --output output/pdf/e5f_evening_housing_refinement_review.pdf --heading 'FERTILITY / EVENING CALIBRATION REVIEW' --no-source-index
```

The economic next step is joint housing and tenure refinement under the full
twelve-row objective, with an explicit decision about extending the first-child
bound. The empirical housing target is retained. Production remains the separate
September 4 `task_010` and no new policy results are attributed to tonight's points.
