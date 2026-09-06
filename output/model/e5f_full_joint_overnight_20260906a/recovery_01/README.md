# Smaller full-objective recovery

The broad overnight attempt 17023172 stopped at candidate 17: housing-market
residual 6.118e-4 exceeded the unchanged 2e-4 gate, including the existing
second-bracket attempt. Fourteen completed trials did not improve the seed.
Other children were canceled, and the old source, output and failed point are
preserved. This is a numerical failure at a tested parameter vector, not an
identification result or a proof that any empirical target is unreachable.

This continuation changes the search method: probe a smaller neighborhood of
the verified seed, update using actual full-objective improvements, and repeat
across all eleven parameters. It does not rerun candidate 17. Every one of the
51 scientific source files, the case adapter, parameter domains, targets and
weights is identical to attempt A. No tolerance, model or measurement changes.
Ten pure controller checks pass, including every coordinate in six rounds and
the final exact original-generator repetitions.

## Predeclared budget and method

1. Smoke the actual launch/observation/collection loop with two identical seed
   histories and two joint proposals. Seed histories must exactly reproduce
   each other, including 17 PNGs, and satisfy the existing fixed-domain bridge.
   Joint smoke proposals perturb all eleven transformed coordinates by an
   alternating signed 0.0003125; all source and numerical gates must pass.
2. Start from the best of these four fully evaluated points. At each of up to
   six rounds, probe all 22 signed coordinate directions at radius 0.00125.
   Use actual objective gains against that round's original anchor to combine
   all/top-three/top-six improving directions at scales 0.5, 1, 2 and 4.
   Evaluate every combination as a complete history with the same twelve rows.
   The best point is retained only for a strictly smaller actual objective.
3. If a round gains less than 0.1%, halve the radius down to 0.0003125. Stop
   when the minimum radius gives an absolute gain below 1e-8. The radius affects
   proposals; the full eleven-parameter bounds remain unchanged.
4. Reserve two original-generator exact repetitions of the best point. Match
   all twelve fit rows, physical parameters, 253 numeric history entries and
   seventeen PNGs. Any failed numerical/source/measurement gate stops all
   owned children and dependent stages. No automatic model retry.

Maximum 4+6*(22+12)+2 = 210 histories. The 14 completed original search trials
had a median runtime around 7.6 minutes; recovery budget estimates use 8--12
minutes per history and 2--3 batches per round. Allow approximately 3--4 hours
including the two sequential smoke batches and verification, subject to queueing.
Smoke requests 2 CPUs/16 GB for 35 minutes; full recovery requests 12 CPUs/96 GB
for at most five hours. The unchanged absolute finish is 2026-09-06 08:47:19 UTC
(04:47 Eastern), and one hour remains reserved for verification. A stage that
cannot fit its declared remaining case/time allowance is skipped. No further
recovery is preauthorized by this design if this one also fails.

## Reproduce and monitor

`../prepare.py --local-recovery` creates a new immutable snapshot from attempt
A, copies the unchanged launch script, and writes this contract. It refuses to
overwrite an existing recovery. `contract.json` SHA:
580723c437640677038966d5f5bcfb1638037446fc319bfb443e160f7832bd6e.
Runtime controller SHA:
6086ecb6549fe8883f4b306c1f976d7154fa5b10480a5cd5190fcf985f4d3f88.
Scientific bundle remains ce38de90a85de7102f4d462bd1f2618fbad17f649c5c57d840d5152e5917dff6.

Remote snapshot:
`/scratch/td2248/projects/Fertility_Spring26_full_joint_overnight_20260906b`.
The relative output path is this folder's repository path. Smoke job 17023921
must pass before submitting the full recovery. `smoke/smoke_verification.json`
needs status pass, exact_standard_pngs 17 and joint_probe_cases 2. The full
controller independently rechecks all four completed cases and source contracts
before any new proposal. Queue submission reuses
`code/cluster/submit_e5f_joint_overnight_search.sh`, with cpu_short/cpu48,
12 CPUs, 96 GB, five-hour Slurm cap, this contract/hash, and E5F_JOINT_MODE=search.
Save search_submission.json with the exact command and returned job ID.

After submitting, start the finite local collector:

```sh
python3 code/model/tools/collect_e5f_joint_overnight_search.py --recovery --watch --job-id JOB_ID
```

It collects only compact evidence and the selected standard graph packet.
It never submits or retries a model job. The distinct output is
`output/pdf/e5f_full_joint_recovery_review.pdf`; visual review remains required.
`search/search_state.json`, `best_so_far.json`, `latest_completed_case.json`,
`all_cases.csv`, final verification and complete fit/parameter tables remain the
live evidence. Heartbeats occur every 60 seconds; investigate after 30 minutes
without progress. The app monitor checks every 15 minutes, stays quiet on
unchanged state, and must pause after completion/review or the final cutoff.

Production remains September 4 task_010. No policy runs, production promotion
or author-owned manuscript edits are part of this continuation.
