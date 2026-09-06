# Full joint overnight calibration

Tommaso authorized autonomous overnight calibration at about 22:15 Eastern on
September 5, renewed after the app update. He is too tired to review tonight.
This task must deliver the best verified fit obtained within its stated budget
and a short readout, rather than another housing-only diagnostic.

The earlier evening search varied eleven parameters within the old bounds and
ended at loss 21.798033. Its final loss-20.695274 point was a conditional
housing diagnostic with nine coordinates fixed. This overnight search instead
allows every one of the eleven estimated parameters to move. It retains all
twelve empirical targets and weights, including first-birth rooms 0.720246.
The young-ownership row remains validation-only. Production is not promoted
automatically, and no policy result will be relabeled as belonging to a new fit.

## Frozen specification and one enlarged search interval

Use `/scratch/td2248/projects/Fertility_Spring26_full_joint_overnight_20260906a/`.
The existing first-child housing-intercept upper search bound is expanded from
0.5 to 2.0 explicitly. This is a search restriction, not a new equation or
external estimate. The historical review found no empirical justification for
the old ceiling, and the recent diagnostic improves beyond it. The wider
interval permits searching that direction; it does not assert that its endpoint
is economically validated. Every other parameter bound is unchanged.

`source_contract.json` and `search_bound.patch` pin the entire change relative
to scientific commit c63821a6e027e90c85db1f8c137e617d4f7b5492. Of 51 fingerprinted
files, only `run_e5f_transition_calibration.py` changes: it adds an explicit
upper-bound argument with the old 0.5 default and validates it. Equations,
preferences, objective, target measurement, market solver, population, fiscal
and numerical gates remain fixed. The active working scientific code is not
edited. The new scientific bundle is
`ce38de90a85de7102f4d462bd1f2618fbad17f649c5c57d840d5152e5917dff6`.
The source-input SHA and target fingerprint are inherited unchanged and pinned
in `contract.json`, every stage plan, runtime and collection.

Dated housing-supply elasticity is externally fixed at 0.63 and tenure taste
dispersion at 0.005. Old initialization retains elasticity 1.75, fertility
normalization 2.1, outside-entry share .169, the household renewal convention,
and the historical household/age bridge. This run does not settle their policy
interpretation. Grids remain 120 wealth nodes, 17 ages and 15 income states.
Twelve target rows for eleven estimated parameters preserve the
maintained target count; the earlier weak bequest direction and nonsmooth old
wealth-to-income quantile remain identification/numerical limitations.

## Predeclared search and budget

1. Two identical full-history evaluations through the actual concurrent launch,
   observation, checkpoint, graph, collection and failure-handling loop. Each
   bridges the selected fixed-housing result to the new estimated domain:
   exact target/weight identities, fit/history differences at most 1e-10 and
   physical-parameter differences at most 1e-12 allow only transform rounding.
   The two new executions must then agree exactly, including all 17 PNGs.
2. A 32-member initial population around the best point, with all eleven
   coordinates variable. The first-child term spans initial values .5 to1.2;
   the other transformed coordinates initially vary by +/- .01. Actual search
   bounds are the complete declared intervals, not these initial radii.
3. At most eight deferred differential-evolution generations of32 trials.
   A mixture of current-to-best/1 (80%) and rand/1 (20%), mutation scale
   uniform[.5,.9], and crossover .9 proposes simultaneous parameter changes.
   Every proposed objective is evaluated with the original complete history.
   Bounds are repaired toward the parent, never changed. Selection accepts
   only strictly lower actual twelve-row loss. Stop early after three tiny
   improvements only if the population's largest unit-coordinate range is<.002.
4. Up to two local refinement rounds: all22 signed coordinate probes, followed
   by up to12 joint directions combining all/top-three/top-six improving
   coordinates at four step lengths. Rank coordinates by actual objective gain.
   Initial radius .0025 halves in the second round. Duplicate joint proposals
   are omitted. This remains the full objective with all eleven parameters.
5. Two exact repetitions of the best point using its ORIGINAL center generator,
   not a transformed/re-encoded fitted summary. Require all twelve fit rows,
   physical parameters, numeric historical entries and 17 PNGs to match exactly.

Maximum histories: 2+32+8*32+2*(22+12)+2 =360. At the observed9--12 minutes per
history, this is54--72 CPU-hours of case work. Up to12 cases run concurrently,
1CPU per child, with8GB per concurrent slot and a hard30-minute case timeout.
The full Slurm allocation requests12CPUs/96GB and at most8.5 hours (102 allocated
CPU-hours maximum); efficient use depends on stage synchronization and case
lengths. Recent timings imply roughly5--7 hours including launch/collection,
plus queueing. The smoke allocation uses2CPUs/16GB and35 minutes maximum.
Search uses at most7.5 hours from launch, reserves verification time, and stops
all work by08:00 Eastern September6 even if scheduling starts late. A stage is
skipped when it cannot fit its remaining time/case allowance. No automatic
numerical retry, no repeated expansion, and no run beyond these caps.

Every completed case has the unchanged17-graph packet, dated checkpoint,
complete12-row fit, all parameter estimates/bounds and a source/gate receipt.
Children heartbeat every60seconds; the coordinator writes state every60seconds
and after each completed batch of futures. Keep `latest_completed_case.json`,
`best_so_far.json`, and the stage/full ledgers. Failures or timeouts stop all
owned children and dependent stages; existing verified results remain on disk.
Thirty minutes without heartbeat requires investigation. Occupied value-screen
flags remain visible and require diagnosis before a winner is recommended;
they are not silently converted into economic conclusions.

## Implementation review and verification

A single bounded worker supplied an initial controller and pure tests. The lead
found missing final repetitions, inadequate failure cancellation/budget gates,
and a polishing comparison against an already-updated incumbent. The lead
replaced these sections before any model launch. The final nine pure checks
cover determinism, all-coordinate mutation, bounds, stale-loss exclusion,
complete evaluated selection, budget reserves, correct polishing directions,
original-generator repetition, failure cancellation and the complete simulated search-to-repetition control flow. The unchanged default
search interval, explicit2.0 interval, invalid argument rejection and full
scientific fingerprint are also checked. The exact cluster smoke must pass
before the long allocation is submitted.

`contract.json` pins the final controller, adapter, planner, scientific bundle,
input, targets, initial generator, seed and budgets. `prepare.py` constructs
these artifacts without solving. Controller commands use
`code/cluster/submit_e5f_joint_overnight_search.sh`; the existing Torch wrapper
was used to probe the connection/queue, and its documented raw-submit route is
used because this new full-history family is not one of its legacy families.


## Launch and continuation

Exact-loop smoke `17022832` completed0:0 in9m22s using2CPUs. Both runs exactly
match all12 fit rows, physical parameters,253 history entries and17PNG hashes.
The bridge from the fixed-housing source differs by at most4.7827e-12 in fit
entries and2.1622e-11 in history entries, within predeclared rounding tolerances.
All11 parameter-domain rows and their expanded/retained bounds are verified.

Full job `17023057` uses cs/cpu48,12CPUs/96GB. The first full submission was
rejected before starting because cpu_short allows only6hours; live QOS metadata
confirmed the correct longer CPU route. No model, controller, contract or gate
changed. `routing_preflight.json` retains both launcher hashes and the rejected
resource request. The two-stage smoke and full source contracts are unchanged.

The finite local collector is `code/model/tools/collect_e5f_joint_overnight_search.py`.
It never submits or retries a model job. It periodically pulls compact evidence,
checks available original hashes, collects the winner's standard graph packet,
and prepares `MORNING_REVIEW.md` plus
`output/pdf/e5f_full_joint_overnight_review.pdf` when the controller reaches a
terminal state. The rendered PDF and new selected graphs require visual review;
its automatic receipt explicitly records that pending check. Computation runs
independently on Torch; local collection resumes when the laptop is awake.
Read `search/best_so_far.json`, `search/search_state.json` and
`search/MORNING_SUMMARY.md` for the authoritative live/completed result. A saved
best point is not labeled verified unless the final repetition receipt passes.


### Final routing update before any full-search evaluation

The cs/cpu48 job17023057 remained pending without a start estimate and was
canceled in that state. No full-search evaluation was restarted. To use the
night productively, job17023172 requests a six-hour cpu_short allocation with
job QoS cpu48, matching the completed smoke's actual QoS. An explicit cpu_short
job-QoS request was rejected before submission; this distinction is recorded in
routing_preflight.json. No credentials, policy or resource limit were bypassed.

The immutable `short_queue_contract.json` SHA is
`4a1baefe628b95b9204f6522a22b86c2a0f667d206a0bf11e0c17c0bd0a46341`.
It changes only the effective finish time and a routing note relative to the
exact-smoked contract: finish08:47:19UTC (04:47 Eastern), with the existing
one-hour verification reserve. Controller, adapter, scientific source, seed,
search proposals, all eleven parameters, targets, weights and numerical gates
remain identical. Maximum360 histories is unchanged; time may truncate the
search. The six-hour allocation caps reserved CPU time at72CPU-hours. The local
collector now follows17023172. The earlier contracts and queue evidence are
retained; the final submission receipt identifies the operative contract.

## Market-gate stop and separately recorded recovery

Monitoring found job17023172 FAILED after14m56s and14 completed search trials
(16 histories including the initial two smokes). Candidate17's housing-market
residual6.118e-4 exceeded2e-4 after both existing bracket attempts. No completed
trial improved the seed. The failure is retained in search/failure.json and
initial_population/task_017/adapter_failure.json; all completed fits, parameters
and dated checkpoints remain. It is not evidence of an unreachable target.

See [the separate recovery design](recovery_01/README.md). Smoke17023921 tests
smaller full-objective coordinate and joint steps with all eleven parameters
free, source and gates unchanged, a210-history maximum and the original04:47
Eastern stop. Full recovery awaits its four-case smoke receipt. The app heartbeat
monitor is `monitor-full-joint-fertility-calibration`, every15minutes, in the
current calibration task; the older numerical-task automation stays paused.
