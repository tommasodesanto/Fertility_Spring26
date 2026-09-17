# Successive permanent-surprise refit

## Final recovery assessment — September 16, 2026

**No fitted shock has been accepted. The user queue is empty.** Recovery
17858740 ended after 09:37:22 with a preserved candidate requiring continuation.
The first candidate passes the 104-date finite market/fiscal root but misses
the fertility target. The second remains unconverged after 11 mappings; later
historical stages and policy were therefore not dispatched.

| First-shock candidate | Preference level | Period TFR | Target | Numerical/fit status |
|---|---:|---:|---:|---|
| Recovered first candidate | 0.1289153 | 1.953246 | 1.974875 | Finite root verified; fertility gap −0.021629 exceeds 0.005 tolerance. |
| Recovered second candidate | 0.1339153 | 1.990199 | 1.974875 | TFR provisional: finite root fails, so no accepted fit gap/bracket. |

All historical targets remain visible below; structural parameters were not
re-estimated. Their unchanged full [target-fit table](../../e5f_final_night_20260913/corrected_initial/target_fit.csv)
and [parameter table](../../e5f_final_night_20260913/corrected_initial/parameters.csv)
belong to the initial calibration, not to a newly fitted historical path.

| Surprise date | Target observation years | TFR target | Final stage status |
|---|---|---:|---|
| 2007 | 2008–2011 | 1.974875 | Two candidates attempted; neither fitted and accepted. |
| 2011 | 2012–2015 | 1.861000 | Not started. |
| 2015 | 2016–2019 | 1.755375 | Not started. |
| 2019 | 2020–2023 | 1.645750 | Not started. |

The second candidate's maximum housing/PAYGO/rebate relative errors are
0.020276%/0.003178%/0.020425%; its scaled joint score is 0.0408493 versus
0.0002. It reproduces exactly, but reproduction does not establish equilibrium.
The controller stopped before a round that could not fit its remaining
candidate budget. Checkpoints survived. Both stationary endpoints passed;
the first converged forecast still has a terminal household-mass gap of 1.742%
and distribution L1 distance of 9.446%. Horizon robustness is unverified.

The separate announced-path solver diagnostic 17865027 also failed after
12 mappings/06:01:49, even with its explicitly relaxed fiscal gate and four
terminal dates omitted from the trimmed acceptance score. Those diagnostic
criteria are not production criteria. Four tenure-dispersion probes completed
but reused the saved candidate and retained the baseline terminal endpoint;
none solved a new equilibrium. Raising dispersion has not been demonstrated
to fix convergence. The earlier untrimmed diagnostic 17860152 was cancelled.

Evidence: [source receipts and scheduler records](recovery_final_readout_20260916.json),
with remote paths, SHA-256 hashes, contracts and native-row residual maxima.
No new jobs, code changes, target/gate changes or slide updates were made by
this readout. The sections below preserve launch and earlier inspection history.

Smoke job 17732511 passed native validation and automatically submitted
first-stage job 17732824 (initially PENDING, Priority). Each stage dispatches
its successor only after its acceptance checks pass. The passed smoke is
not a completed fit. The smoke and dispatch receipts are saved here.

| Stage | Surprise | Observed birth years | Forecast decisions |
|---|---:|---|---:|
| 0 | 2007 | 2008–2011 | 104 |
| 1 | 2011 | 2012–2015 | 103 |
| 2 | 2015 | 2016–2019 | 102 |
| 3 | 2019 | 2020–2023 | 101 |
| Policy | 2023 | Unexpected permanent 2% tax | 100 |

All forecasts use the 2423 boundary. Households expect only the latest
preference to remain permanent; later shocks are surprises. Four fitted
levels, fixed calibrated 2007 structure, original birth-vintage law, no
immigration/rescaling, fixed asset-price housing supply, equal tax rebates,
and balanced PAYGO. The last accepted forecast supplies its own inherited
2023 households and baseline tail. No extra 1% policy solve is needed.

Each preference candidate receives its own verified stationary endpoint.
The policy endpoint is recomputed using the fitted final preference. The
measured derivative profiles seed a scaled-step Broyden root, with the
unchanged market/fiscal and household gates. Derivatives outside the measured
lag range are zero in this approximate initial matrix; the native equilibrium
evaluation remains the criterion for acceptance.

Remote batch:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/long_successive_refit_20260914a`.

Status: `output/smoke/{passed,failure,dispatch}.json`; then
`output/stage_N_YEAR/{latest_completed,best_so_far,accepted,failure,dispatch}.json`.
Per-candidate round folders contain native rows, fertility, residual receipts,
terminal distances, plots and standard diagnostic galleries. Each accepted
stage writes `accepted_next_state.pkl.gz` plus a verified hash and round-trip
state/queue comparison. Policy writes the matching baseline rows and fertility.

Stop caps are 30 minutes for smoke, one hour per endpoint, ten hours per
candidate, 24 hours per fitted stage, 12 hours for policy and seven calendar
days overall. Each candidate has at most four eight-mapping root rounds;
each stage at most six trial values. These caps do not predict completion.
Finite convergence is not terminal convergence or horizon adequacy; those
remain explicitly separate, and all results remain provisional until checked.

Sources: `code/cluster/run_e5f_long_successive_refit.py` and
`code/cluster/prepare_e5f_long_successive_refit.py`. Detailed specification:
`docs/model/e5f_long_horizon_successive_surprise_plan.md`.

The September 14 requested first-candidate transition graph and comparison of
the four-announced-shock tax policy with the frozen slide result are in
`readout/README.md`. These are separate diagnostic artifacts; the deck is unchanged.

Final job inspection after login renewal is in `final_job_readout_20260915.json`.
The 24-hour refit ended without an accepted shock. A later first-candidate
iterate reached the residual threshold but timed out before final verification;
it is saved for recovery, not accepted. The separate four-shock policy ended
at its evaluation budget and remains unconverged. `CALIBRATION_STATUS.md`
contains the readable final assessment. No jobs were relaunched by this inspection.

## Authorized recovery, September 15

`recovery_manifest_20260915.json`, `recovery_submission_20260915.json` and
`recovery_source_manifest_20260915.json` describe the separate recovery batch
`long_successive_refit_recovery_20260915a`. Smoke job 17858190 was observed
RUNNING. The old experiment and its output are untouched.

The smoke first validates both cached roots and freshly re-audits both terminal
states, then runs the existing native loop smoke. Success automatically submits
stage 0. Its first two candidates reuse prior roots and terminal states; subsequent
candidates use the retained scalar bracket/secant search. Only certified native
roots contribute signed fertility gaps, and only a fitted root may dispatch the
next historical stage. Policy still requires the fitted inherited 2023 state.

Expected first recovery: two full 104-date mappings at observed 2,665–3,413
seconds each, approximately 1.5–2.2 hours with reporting. Budgeting uses a
conservative 3,900 seconds per mapping plus 900 seconds for artifacts. The root
reserves its final evaluation within the chosen mapping count. If the next round
cannot fit initial plus final evaluations, the candidate record is saved; the
stage reports that continuation is required. No tolerance or scientific
specification changes are part of this operation.

All original caps remain: up to four rounds per candidate, six candidates per
stage, 10 hours per candidate, 24 hours per stage and seven days for the chain.
Every completed native mapping writes rows, fertility, graphs and best/latest
receipts; the driver writes a one-minute heartbeat. A missing heartbeat for
30 minutes requires investigation. Queue waiting time is outside solve estimates.

Ten local routing and recovery tests, remote manifest/coordinate validation and
source compilation passed before submission. Native smoke and long recovery
remain subject to their saved gate receipts; scheduler RUNNING is not acceptance.
