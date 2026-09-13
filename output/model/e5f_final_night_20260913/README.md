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

## Failure recovery and exact policy reuse

Search17596347 stopped after three matching candidate numerical-budget failures.
The explicit traceback was the20-evaluation initial root limit; valid candidates
remain in the ledger. Recovery17603133 re-evaluates the verified best twice in
`initial_search_recovered`, preserving the stopped search. Handoff17605279 will
use that independently verified result; the earlier handoff17597052 is obsolete.

Native cache probe17598785 reproduced every comparison exactly (including all
dated values, policies, final household distribution, economic rows and audit
residuals):183.9926seconds uncached versus43.5013seconds cached, with1 actual
dated solve and11 hits in this constant-input six-date mapping. The cache is
bounded and its key covers the complete solver arguments, including shared
state and continuation values; hits return fresh arrays. Unsupported argument
states bypass the cache. This is a measured mapping speedup, not a claim that
every long solve is4.23times faster. Probe17598460 failed its harness time-reserve
setup before a model mapping;17598785 used a zero-policy-reserve probe manifest.

The pending refit source allows24 root evaluations instead of8, since the short
roots were still converging at the earlier limit. The market/fiscal/reproduction
gates are unchanged. It shifts numerical guesses one date when carrying a
realized state. Six refits have a maximum576 primary mappings per track before
one bounded alternative start, constrained by the shared deadline and policy
reserve. Cache bounds are6GiB for6/24dates and24GiB for100dates; refit jobs request
64GiB memory. The worker ceiling remains18 and the planned overlap is14 numerical
workers. Native cache proof and source hashes govern admission.

Collector17607147 supersedes17597285 and reads the recovery/cache receipts and
discovers the follow-up array automatically. Weekly allowance at this check was
81% remaining; subsequent healthy monitoring should stay compact.
# Latest verified continuation

Recovery17603133 passed two exact numerical repetitions. Selected original loss:
179.2984252281 versus182.6491468669 for the rebated seed. The search itself stopped
after repeated numerical evaluation limits. Complete13-row fit and17-row parameter
tables, accounting gaps and caveats: `initial_search_recovered/README.md`.

Handoff17605279 completed and submitted array17608564: six cached A0/A+ refits
at6/24/100dates, root cap24, unchanged tolerances, fixed18:00UTC deadline.
Original17595967 remains independent evidence. Cluster collector17607147
automatically includes the refit array. A complete accepted historical path and
policy results remain outstanding. Age pilot17597260 has exact matching numerical
results but failed its broad checkpoint-hash comparison; recover the report into
a separate folder without another solve or changing the main target system.

Resource-only queue replacement: pending17608564 was cancelled before execution and split into17613033 (6dates,16GiB),17613034 (24dates,32GiB),17613035 (100dates,64GiB), two cases per array. Commands, output folders, source pins, tolerances and deadline unchanged. See `resource_resubmission.json`. Age report recovery passed without solves; full tables and limitations in `age_pilot_recovered/README.md`.

After the new6/24date paths passed their first complete mapping audits, superseded uncached24/100tasks17595967_{1,2,4,5} were cancelled to free memory for100date replacements. Old6date tasks0/3 remain as comparison; all outputs preserved. `collect_e5f_final_history_readout.py` is a no-solve extractor prepared for accepted2019/2023 native snapshots. Compilation and schema checks passed; native extraction remains pending an accepted final history window.
