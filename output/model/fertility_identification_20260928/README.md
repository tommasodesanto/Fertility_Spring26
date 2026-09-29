# Current-point fertility identification and reoptimization

Experiment index: see [EXPERIMENT_REGISTER.md](EXPERIMENT_REGISTER.md) for the frozen reference, completed experiments, and authorized experiments awaiting budgets and source certification.

Author-authorized September28 follow-up. Main18716710 submitted after six accepted smokes; see launch.json and approval_v1.json.
Frozen reference is overnight selected block0506, primary loss19.581310760.
No model, earnings, entry distribution, target value, grid, bound, closure or
gate changes. Three untargeted checks remain visible. Child benefit is always
renormalized to completed fertility2.100, retaining the positive-benefit gate.

Plan: six exact full-loop reference smokes (two per free-weight lane) with lead acceptance; forty normalized
central differences (all10coordinates, full and half steps); six population
searches of at most80objectives each; two fresh repeats per selected lane and17
standard plots. Six-hour controller cap, search stops at five hours, repeats
stop at5h50, exports within six hours. Maximum538 objectives,24single-thread
Torch workers,192GiB;1800-second objective cap and unchanged23SS normalization
cap. Median prior objective769s,90th percentile1114s:538/24*769=4.8worker-wall
hours idealized, with six-hour hard cap rather than guaranteed completion.
No Mac model work. Checkpoint every case/at most5minutes;30min stale diagnosis.

Lanes: original primary weights;10x and100x earlyfertility weights with all10
coordinates free; fixed continuation-scale profiles at half,double,quadruple
reference with other9free, each using10x earlyweight. At least two starting
points (reference plus successful earlyfertility-oriented overnight point).
Compare all under common original primary weights, raw14moment gaps, and
other-moment primary loss excluding earlyfertility. Fixedparameter profiles
are not exact-target constraints. No target swapping or demotion is adopted.

Current Jacobian must report step-size sensitivity, derivative of normalized
child benefit, scaled10-by10 scored-moment SVD and remaining full14rows.
No global identification, unreachable-target or optimum claims from this finite
exercise. Candidate acceptance requires original scientific gates and exact
repeats, with all failures and unrun cases retained.

## September 28, 11:22 EDT — launch dependency repair
Importcheck18716108 COMPLETED/PASS, exact pinned controller-helper import and
contract verification, zero solves. Replacement smoke18716129 submitted using
immutable run_v2.sh SHA02eb80711e5b18104d8f73ff8d7d4a80b73dd266005a9efd90043f05d30b50b5.
Outputs use smoke_v1 (previous failure never created it); original failure log
identification_smoke_v1.log remains intact; new log identification_smoke_launcher_v2.log.
Main remains unapproved; use run_v2.sh for subsequent approved main submission.


Smoke18715827 failed before clock creation or any model solve: dynamic import
of the pinned recovery scheduler could not find its sibling
`e5f_utility_comparison_design`. No search approved. Original launcher and failure
log retained. Versioned `run_v2.sh` adds the existing portable helper directory
to PYTHONPATH; no source, model, objective or contract change. Torch importcheck
18716108 checks pinned contract verification and the exact failing helper import
with zero solves. Await its PASS before a replacement smoke submission.
No numerical clock exists yet; the six-hour cap still starts at smoke execution.


## September 28, 11:44 EDT — verification passed; main submitted

All six exact-loop smokes pass. Independent Torch review18716599 confirms source
pins, all14 physical target rows/all31 parameters against anchor, three exact
within-lane table pairs, correct1x/10x/100x early weights, and all17 PNGs identical
across six cases. Lead viewed all17 plots; inherited high-wealth ownership, age30
housing and retirement-profile caveats remain. Search anchor has no PNGs;
comparison is cross-smoke, not an asserted anchor-image comparison.
Approval_v1.json written; main18716710 submitted ONCE using run_v2.sh,24CPU192GiB.
Original clock retained: search16:23:14 EDT, repeats17:13:14, end17:23:14.
Current Jacobian then six bounded search lanes; no model or target changes.

## September 28 13:36 EDT — authenticated continuation submitted

Original job18716710 is terminal FAILED with75 records:59 success (40 Jacobian,
19 search),14 censored timeouts and2 known nonpositive-benefit rejections
misclassified by the controller. All original evidence remains immutable.
Authentication job18721565 completed exit0 and verified all75 records under
original pins, with zero solves. Lead checked manifest and inherited clock.
Resume job **18721946** submitted ONCE via resume_v1.sh; pending Priority at13:36.
Do not duplicate. No attempted case is rerun; original538 total cap and16:23 search /
17:23 hard end remain unchanged. Model, targets, bounds and gates unchanged.
Continue10minute monitoring until running, then30minutes. No final certification yet.


September28 13:43 EDT: resume18721946 RUNNING on cs677,24 actual evaluators.
Compute-node heartbeat verifies75 imported records,24active,37pending and unchanged
clock, no stop reason. Login-node directory metadata briefly stale; use compute-node
receipts if needed. Monitoring restored30minutes. No new completed cases yet.

September28 14:15 EDT monitoring: resume18721946 healthy,24 actual workers;
100 records =75 imported+25 new (11success,14timeouts); total70success,28timeouts,
2 preserved known rejections, zero new fatal errors. Heartbeat1sec/checkpoint37sec;
MaxRSS126971456KiB within192GiB allocation. Best points unchanged, no intervention.
Original deadlines and gates unchanged. Last plot review13:27; refresh next check.

September28 15:16 EDT: resume18721946 running24workers,136 records (88success,45timeout,
3inadmissible), no new fatal; generation1 active, heartbeat fresh, MaxRSS121.1GiB/192GiB.
Primary best unchanged. Half-scale profile now0133, not final certified. Saved renderer
18724418 PASS zero solves: four lanes all17 hashes unchanged, half0073 all17 inspected;
market residual7.11e-8, highwealth ownership/age30housing/retirement caveats remain.
Latest half0133 postdates rendered snapshot. No model/gate/deadline changes.

September28 15:46 EDT:161 records,99success/59timeouts/3inadmissible;24 actual
workers, fresh heartbeat and408second checkpoint age. MaxRSS125.4GiB/192GiB.
No new fatal errors or best-point changes. No intervention or deadline extension.

September28 16:16 EDT:186 records (105success,78timeout,3inadmissible),24workers,
fresh heartbeat; checkpoint410seconds, memory125.4GiB/192GiB. Best unchanged.
Verified original controller freezes selection at16:23 and requests two repeats
per available lane before17:13, then exports before17:23. No extension.
Saved-only renderer18731455 queued v4 to refresh half0133; zero objectives.

## Local Git protection and review caveat — September 28

Author authorized prevention of local memory overload. Added *.pkl.gz ignore
rule and verified an example checkpoint path is ignored. Repository-local
gc.auto=0 and maintenance.auto=false disable automatic maintenance; pack.threads=2,
pack.windowMemory=256m and pack.deltaCacheSize=128m bound specific packing caches
and concurrency, not total process RAM. No existing checkpoint data or Git objects
deleted; the oversized retained snapshot still needs separate recovery/cleanup.

Claude review's 0.665 early-fertility ceiling conditions on the empirical first-birth
cell shares, not merely their targeted mean. Matching one mean does not fix those
shares. Therefore this is a conditional timing-distribution bound, not proof that
the current target system is infeasible. The observer does confirm age25 interpolation
with post weight0.875; coarse birth spacing merits measurement checks. No target
swap or economic change adopted. Predicted Gauss-Newton gains remain untested.

September28 16:46 EDT: search stopped at cutoff. Six repeats successful (primary,
early10,early100 pairs), four profile repeats active; final equality/export pending.
218 records:113success,102timeouts,3inadmissible. v4 saved-only renderer passes;
all17 half0133 plots inspected, other four lanes17hashes unchanged. Half0133 housing
residual1.65e-6; existing wealth/age30housing/retirement caveats retained.

## September28 17:17 EDT — experiment complete and authenticated

Torch18721946 COMPLETED exit0. All10 repeats pass, with exact14target/31parameter
comparisons and17PNG pair identity. Independent export hash check passes for all
five available lanes; each17PNG packet matches visually reviewed saved_review_v4.
Local compact complete tables: resume_v1/selected_export/<lane>/{target_fit,parameters}.csv.
No binary model checkpoints downloaded. Main primary remains19.581310760; no baseline
improvement. Early100: early0.558505, primary88.614507; doubled continuation:early0.687888,
primary807.454086, first-birth age23.693. Half continuation:early0.505300, primary152.816537.
Quadruple lane unavailable. Total228 attempts includes6smokes+40Jac+172search+10repeats.
Final records117success,102timeouts (includes cutoff censoring),3inadmissible.
Do not equate these timeouts with economic nonexistence. No new model/target/gate
changes; no optimum or unreachable-target claim. Search and all numerical work
finished within original window. Next author-authorized questions: verify actual
saved income process and decompose lifecycle/age25 first versus subsequent births.
Claude's0.665 ceiling is conditional on cell shares, not implied by a targeted mean.
