# Authorized final-night work

Author requested immediate launch of the plan in
`docs/model/e5f_two_closure_overnight_plan.md`. This folder records actual
submission and evidence, which must not be confused with the full planned matrix.

Remote batch:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/`.

| Job | Work | Status at submission readout |
|---|---|---|
| 17592167 | Saved-policy demographic operator | Completed, two identities pass |
| 17592542 | Superseded initial smoke source | Cancelled; not a calibration result |
| 17592728 /17593512 /17593865 | Slow nested initial method and dependent searches | Cancelled/superseded |
| 17594990 | Joint price/preference/rebate initial smoke | Completed, full numerical/scoring/accounting gates passed |
| 17596347 | Nine-coordinate rebated initial search, six workers | Running |
| 17595967 | A0/A+ ×6/24/100 fitted histories and rebated-tax policies | All six array tasks running |
| 17597260 | Experimental initial fertility-age profile, two workers | Running; source/target preflight passed |
| 17597052 | Refits after an exactly reproduced initial improvement | Depends on search success; checks improvement and remaining time |
| 17597285 | Cluster-side five-minute receipt collector | Submitted, independent of laptop |

The array and search coexist with the age pilot within fourteen numerical workers, below the shared
eighteen-worker cap. Each worker has one numerical thread. The initial search
has a three-hour limit with35minutes reserved for exact verification; history
tasks have twelve-hour Slurm limits. There are at most36 original-objective
search proposals and24 preference trials plus one alternative start per history.
Numerical candidate failures are isolated; proven contract corruption or three
matching systemic failures stop the search branch. The accepted joint smoke took
about six minutes; the first6-date forecast mapping took165seconds. These are
observed operation times, not forecasts of full convergence.

Full initial smoke tables: `joint_initial_smoke/selected_target_fit.csv` and
`joint_initial_smoke/selected_parameters.csv`; the complete score is unchanged.
This is a validated seed, not a newly optimized result.

Forecasts retain the original twelve scored targets, all weights, nine free
coordinates and annual beta cap0.99, plus separate initial fertility2.1.
Property-tax revenue is returned equally and PAYGO pensions balance separately.
Four successive surprises are fitted; current preferences are expected to
persist and remain fixed after2023. Historical A head-age conditioning remains
an imposed bridge. The migration experiment applies to the post2023 future
segment, including forecasts formed before2023.

The finite-boundary implementation evaluates remaining-lifetime household
values at constant boundary conditions. Boundary prices, pensions and transfers
are jointly rooted with dated values, using actual carried households. No
stationary population or hidden entry is imported. Fiscal feasibility after the
boundary is not established; horizon comparisons remain required.

B0/B+ full histories have not been launched: the formation-unit restriction and
B+ person-migration allocation remain unresolved. The demographic preflight
reports the conversion that would keep one saved household distribution's
mass constant only as a diagnostic; it is not an empirical formation estimate.

Verification before submission: 19 history/boundary tests passed on Torch;
six initial-wrapper tests passed on Torch;77household-kernel source files were
identical across the pinned initial and historical snapshots. The full initial native smoke passed. Native history
root and fit gates still govern admission. Full2023observations are saved as native
snapshots; complete authoritative table reconstruction remains a separate task.

Source entry points are `code/model/tools/run_e5f_final_rebated_history.py`,
`run_e5f_rebated_initial_overnight.py`, `run_e5f_rebated_initial_search.py`, and
`code/cluster/prepare_e5f_final_night_manifest.py`.

The follow-up controller automatically submits six separate history/policy refits
if the original-objective search improves its seed and reproduces twice. It
preserves the first array, caps combined numerical workers at14, counts queue
delay against the fixed September13 18:00UTC deadline, and skips refits if fewer
than three hours remain. The experimental age objective is not automatically
promoted into those original-objective histories. Cluster collection writes
`monitor_summary.json` every five minutes and discovers the follow-up array.
The app heartbeat checks meaningful changes every30minutes when the app is
available; the cluster jobs, follow-up submission and collector do not need it.

Pilot17597051/17597204 failed test setup before numerical work because local
fixture paths were unavailable on Torch.17597260 uses the three portable tests,
with the complete source/target and native observer checks retained by the
controller; the local fixture tests already passed. Its source is frozen under
`age_source_v2`. No result is promoted from the failed test launches.

The displayed initial parameter table overlays the enforced annual beta upper
bound0.99 on the frozen scorer metadata, which still prints0.9995. The original
scorer parameter rows are preserved in selected_parameters_raw.csv. This is a
reporting correction only; the candidate/search bound was already enforced.
