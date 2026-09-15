# Successive permanent-surprise refit

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
