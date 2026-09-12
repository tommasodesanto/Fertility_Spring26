# Beta-restricted calibration searches

**Maintained author instruction: beta is estimated with a 0.99 ceiling, not fixed.**
The author clarified this after the fixed-beta recovery was launched. The original
lower bound remains 0.94. All nine structural coordinates remain free; the
unchanged scorer's original beta upper bound of 0.9995 is tightened externally
by the pinned search plan and enforced before each solve and after each score.
The target/weight fingerprint is unchanged; the new parameter restriction is
explicit in the plan and reporting metadata. This is not a claim of identification.

**Capped search 17425504 is RUNNING on 18 cores (September 11, 20:32 EDT).**
`run_capped_beta.py`, `test_run_capped_beta.py`, `plan_capped_beta_099.json`,
`submit_capped_beta.sh` and `capped_seed_manifest.json` define this separate batch.
It starts from the saved boundary candidate r0_joint_01, whose five source/score
receipts are pinned. Two fresh exact seed repetitions must match all numerical
cells before the search. Up to three adaptive rounds each use at most 18 feasible
derivative probes and 12 joint constrained proposals, retaining all nine
directions; beta uses an inward derivative at its upper bound. The actual first
stage has 16 probes because beta and h_P are at their respective upper bounds.

Budget: 18 CPUs / 96 GiB; at most 90 search trials plus two seed repetitions and
two final repetitions (93 case calls / 94 repetitions / 752 maximum stationary
solves). Expected 65–100 minutes excluding queue, with a three-hour hard limit,
7800-second search budget and 3000-second verification reserve. Every stage checks
its worst-case remaining time. The numerical source, targets, weights, utility,
pension accounting and all strict acceptance gates remain unchanged. Known
inadmissible trials receive no score; unknown errors, recurrent mass failures,
a majority of housing failures or a missing derivative direction stop for review.

Nine controller tests passed locally and on Torch: complete three-round search
with an interior beta optimum, all nine directions, cap rejection, exact replay,
rejected-trial handling, deadline reserves, and actual bounds in output tables.
The source/seed preflight passed, including the actual feasible beta probe
0.9897990204537143. Actual seed repetitions remain pending at launch. Results
save every target/parameter row, 30-second heartbeats, per-case latest/best
records, and the unchanged 17 diagnostic graphs. Remote batch:
`/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/capped_beta_099_20260911/`.

**Boundary comparison: fixed beta 0.99 recovery job 17424705 remains running.**
Its original launch status follows; it is no longer the maintained estimation
restriction. Its completed trials remain useful comparisons for the capped search.
The earlier request to push harder on beta 0.99 was initially implemented as a
fixed-beta profile; the subsequent clarification above supersedes that interpretation.
`plan_recovery_beta_099.json` and `submit_recovery.sh` define this separate batch;
the original remote batch and its sources remain intact. The recovery pins 87
evidence files in `recovery_manifest.json`, validates all 15 successful first-round
derivatives, and excludes the failed theta1-positive trial. It must freshly
reproduce every numerical cell of the original fixed-beta seed before using
those derivatives. The missing derivative side uses the verified negative side;
all eight parameter directions remain present.

The continuation evaluates 12 joint proposals, then one further round of 16
derivatives and 12 joint proposals around the best valid candidate. Eight workers,
at most 40 new search trials, one seed replay, and two final repetitions:
42 new case calls / 43 repetitions / at most 344 stationary equilibrium solves.
Expected duration is 45–75 minutes excluding queue; the hard cap remains three
hours with 7800 seconds for search and 3000 seconds reserved for final verification.
Thirteen controller tests passed locally and on Torch, including complete recovery
and exact-repetition handling. Remote import/proposal preflight passed. Actual
seed replay is still pending; no improved result is yet claimed.

The population tolerance remains 1e-8. Only a narrowly recognized, preflighted
sequential-age advancement mass failure with relative gap in (1e-8, 2e-8] can be
recorded as an **ineligible** trial; it receives no score and cannot enter either
the Jacobian or selection. The original failure counts toward a total allowance
of one, so any recurrence stops the search after the current batch for review.
Unknown failures also stop. All failure artifacts remain available. Inspection
suggests unnormalized float32 tenure probabilities in calendar advancement as a
possible cause, not a proven diagnosis. Neither solver code nor gates were changed.
Results, heartbeat, latest and best summaries are under the remote batch's
`results/`; the selected result includes the full fit/parameter tables and the
unchanged 17 diagnostic graphs. This is an initial-state calibration profile,
not historical preference estimation or a certified new policy result.

**Original run:** array 17403262 finished. Task 0 (annual beta 0.98) completed
with loss 222.524449 and two exact repetitions; task 1 (0.99) stopped after its
first derivative stage at the mass gate, best unrepeated loss 168.201074.
Original run plans below retain their original pins and are historical receipts;
the updated local controller is launched with the new recovery plan only.

Author approved these profiles on September11. Each profile reoptimizes
all eight other structural coordinates against the same12scoredmoments and
separate completed-fertility normalization2.1. No target, weight, original source,
observer, numerical tolerance or economic equation is changed. The unrestricted
selected calibration is initialization and a reference only: it cannot be
selected as a profile result. These profiles do not impose a permanent global
beta ceiling or establish identification.

Each of two independent Slurm tasks uses8CPUs/64GiB, with a3h cap,7800-second
search budget and3000seconds reserved for final verification. At most two rounds
of16coordinate derivatives and12joint proposals; up to3fixed-beta seed attempts
(the two H0±2%alternatives only if the central seed fails the exact known housing
gate), then two exact final repetitions. Maximum60case evaluations/61repetitions
and488stationarysolves perprofile. Prior case times300–490seconds imply roughly
70–120minutes excluding queue, but conservative per-stage guards may stop earlier.
Each two-wave search stage requires4200seconds remaining under the wrapper cap.
No time/case budget is an accuracy or global-optimality guarantee.

The controller pins the original wrapper/objective/source and the exact selected
seed receipts. It verifies every candidate's fixed beta and all other parameter
values. The unchanged raw scorer continues to describe its original nine-free
contract; separate profile tables/metadata explicitly mark eightfree parameters
and beta as an additional fixed restriction. Every target remains in the saved
full fit table, and selected states retain the stable17diagnostic graphs.
Known housing-equilibrium rejections do not stop unrelated valid candidates;
unexpected source, observer, accounting or code errors stop for review.

`run_profile.py` and `test_run_profile.py` implement/test the actual adaptive loop.
`plan_beta_098.json`, `plan_beta_099.json` and `submit.sh` define the independent
cluster tasks. The existing exact scored-model-loop smoke is pinned. A complete
stub-based adaptive-loop smoke tests both rounds, both beta profiles and final
repetition handling; actual fixed-beta seed solves must pass the complete
source/normalization/equilibrium/accounting/observer gates before any search.

Literature context: De Nardi, French and Jones (2010), *Why Do the Elderly Save?
The Role of Medical Expenses*, JPE118(1), sectionVIII.C,p69, reports beta0.99 in
the endogenous-medical-spending extension and0.97 in its benchmark. This provides
an example of0.99, not a universal upper bound or an identification argument for
this housing/fertility model. Primary fulltext:
https://users.nber.org/~denardim/research/De_Nardi_French_Jones_JPE_2010.pdf
