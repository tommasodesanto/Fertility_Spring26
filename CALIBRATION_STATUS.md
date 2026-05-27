# Calibration Status

Updated: `2026-05-26 18:10 EDT`

This is the single live calibration and model-status note for the current
discrete-time center-periphery fertility model. Historical MATLAB status notes,
handoffs, and planning memos have been archived; use them only as background.

## Active Codebase

The active model implementation is the Python port:

- model setup and targets:
  `code/model/dt_cp_model/parameters.py`
- equilibrium solver:
  `code/model/dt_cp_model/solver.py`
- SMM objective:
  `code/model/dt_cp_model/objective.py`
- direct-geometry calibration logic:
  `code/model/dt_cp_model/direct_calibration.py`
- local/cluster calibration entry point:
  `code/model/tools/calibrate_direct_geometry.py`
- result collector:
  `code/model/tools/collect_direct_geometry_results.py`
- active Torch launcher:
  `code/cluster/submit_python_direct_geometry_overnight.sh`

Archived MATLAB reference code:

- `calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/`
- `calibration_archive/legacy_matlab_2026-05-07/code_matlab/`

## Current Closure

Benchmark closure:

- `population_closure = outside_option_benchmark_normalized`
- The model is stationary but not fixed-population in counterfactuals.
- At the benchmark, the stationary scale identity is imposed mechanically with
  \(S=1\):
  \[
  S E_0(p)=q^E(p)\left[M+S B_0(p)\right],
  \qquad
  M=\frac{E_0(p)}{q^E(p)}-B_0(p).
  \]
- The empirical outside-margin normalization is an outside-origin entrant
  share, not \(q^E\) itself. Baseline uses
  \(s^{E,\mathrm{out}}=0.169\), based on ACS 2012--2023 cumulative
  across-CBSA arrival rates from ages 18--22 in the MMS metro sample.
- The code's required city-entry probability is then
  \(q^{E,*}=(1-s^{E,\mathrm{out}})/(B_0/E_0)\). Using the previous
  near-calibrated \(B_0/E_0\simeq 0.93\) gives \(q^{E,*}\simeq 0.89\).
  Because the benchmark equilibrium is invariant to \(q^{E,*}\), the exact
  \(q^{E,*}\) should be recomputed after the new best fit to hit
  \(s^{E,\mathrm{out}}\) exactly.
- `outside_entry_flow` and `outside_value` are benchmark accounting objects,
  not calibrated SMM parameters under the default direct-geometry setup.
- `implied_total_population` is not an SMM target under this closure.
- Counterfactuals should hold benchmark \(M\), \(\bar W^E\), and \(\kappa_E\)
  fixed, then let \(S=q^E(p)M/[E_0(p)-q^E(p)B_0(p)]\) move.

Finite-scale guardrail:

- Candidates with \(E_0(p)-\rho B_0(p)\le 0\), negative outside entry flow, or
  nonfinite scale diagnostics should be rejected before objective comparison.

## Current Target System

Ownership targets now use ACS/MMS household heads rather than person records,
so adult children living in parent-owned homes are not counted as owners. The
audit script is `code/data/mms_center_periphery/audit_ownership_targets.R`; its
diagnostic packet is written to
`code/data/mms_center_periphery/output_ownership_audit/`.

Current ownership targets in `parameters.py`:

| Moment | Target | Source |
|---|---:|---|
| `own_rate` | `0.575` | ACS/MMS heads, ages 30--55, DUE housing restrictions |
| `own_gradient` | `0.186` | ACS/MMS heads, ages 30--55, periphery minus center |
| `own_family_gap` | `0.168` | ACS/MMS heads, ages 30--55, new parents minus childless |
| `old_age_own_rate` | `0.764` | ACS/MMS heads, ages 65--75, DUE housing restrictions |
| `old_age_parent_childless_gap` | `0.070` | PSID completed-children target; ACS co-resident `NCHILD` is not the right object |

The lifecycle ownership slope, ages 65--75 minus ages 25--34, and the
prime-age childless renter/owner median room levels are diagnostic only. They
are reported when available but are not hard SMM targets.

Tenure segmentation now allows small owner units while preserving the renter
cap:

- benchmark owner ladder:
  `H_own = [2.0, 4.0, 6.0, 8.0, 9.5, 11.0]`
- fast owner ladder:
  `H_own = [2.0, 4.0, 6.0, 11.0]`
- renter cap remains `hR_max = 8.0`

## Current Overnight Reduced-Target Runs

Launched: `2026-05-27 00:13 EDT`

Purpose: recalibrate after dropping the lifecycle ownership slope and
prime-age childless renter/owner median room levels from the hard SMM target
system. Local `main` and Torch scratch active calibration files were checksum
matched at commit `27c87f3`.

Verified target counts before launch:

- base targets: `15`
- direct targets: `17`, including `inv_pop_share_C` and
  `inv_rent_ratio_C_over_P`
- direct parameters: `15`

The current Python worker has one optimizer family, a seeded random/global plus
local adaptive proposal search. To get algorithmic variation without changing
model code, the overnight launch uses two proposal regimes per room
specification:

| Room spec | Regime | Run tag | Slurm waves |
|---|---|---|---|
| `hR_max=8.0` baseline | default proposal mix | `py_direct_reduced_targets_hR8_default_overnight_20260527` | `9670117`, then `9670118` |
| `hR_max=8.0` baseline | global-heavy proposal mix | `py_direct_reduced_targets_hR8_globalheavy_overnight_20260527` | `9670119`, then `9670120` |
| `hR_max=6.0` diagnostic | default proposal mix | `py_direct_reduced_targets_hR6_default_overnight_20260527` | `9670121`, then `9670122` |
| `hR_max=6.0` diagnostic | global-heavy proposal mix | `py_direct_reduced_targets_hR6_globalheavy_overnight_20260527` | `9670123`, then `9670124` |

Update after launch: the first 8 tasks from each cold-start first wave were
left running. The pending cold-start tasks `9--16` and second waves were
canceled after confirming that the saved incumbents were much better under the
reduced target system than the generic seed-bank starts. Incumbent losses
computed from saved moments under the reduced target system:

- `hR_max=8.0` incumbent: `16.097`
- `hR_max=6.0` incumbent: `35.230`

Additional incumbent-seeded arrays:

| Room spec | Regime | Run tag | Slurm job |
|---|---|---|---|
| `hR_max=8.0` baseline | default proposal mix | `py_direct_reduced_targets_hR8_incumbent_default_overnight_20260527` | `9670313` |
| `hR_max=6.0` diagnostic | default proposal mix | `py_direct_reduced_targets_hR6_incumbent_default_overnight_20260527` | `9670315` |
| `hR_max=8.0` baseline | tight local proposal mix | `py_direct_reduced_targets_hR8_incumbent_tight_overnight_20260527` | `9670412` |
| `hR_max=6.0` diagnostic | tight local proposal mix | `py_direct_reduced_targets_hR6_incumbent_tight_overnight_20260527` | `9670413` |

The first global-heavy incumbent arrays, `9670314` and `9670316`, were canceled
after the default/global-heavy starts accumulated roughly 60 evaluations per
run without improving the incumbents. The replacement tight local arrays use
`DT_DIRECT_GLOBAL_PROB=0.02`, `DT_DIRECT_INITIAL_SCALE=0.05`,
`DT_DIRECT_MIN_SCALE=0.003`, `DT_DIRECT_SHRINK=0.60`, and
`DT_DIRECT_STALL_WINDOW=6`.

Two-hour checkpoint, `2026-05-27 02:35 EDT`:

- best `hR_max=8.0`: run
  `py_direct_reduced_targets_hR8_incumbent_tight_overnight_20260527`,
  worker `7`, evaluation `295`, loss `10.697`, \(TFR=1.734\), ownership
  `0.337`, center share `0.467`, rent ratio `1.087`
- best `hR_max=6.0`: run
  `py_direct_reduced_targets_hR6_incumbent_default_overnight_20260527`,
  worker `3`, evaluation `274`, loss `15.474`, \(TFR=1.930\), ownership
  `0.622`, center share `0.432`, rent ratio `1.180`
- Status: both incumbent-local specifications improved relative to the saved
  incumbents, so the default and tight incumbent arrays were left running.

Launch settings:

- setup: `benchmark`
- bounds: `global`
- closure: `outside_option_benchmark_normalized`
- geography weight: `100`
- array size: `1-16%8` per run, four first-wave arrays active together
- chained waves: second wave submitted with `afterany` dependency on the first
  wave and the same run tag, so workers resume through existing `best.json` and
  `evaluations.jsonl`
- per-worker budget: `13,500` seconds per wave; Slurm wall time `03:55:00`
- default proposal: `DT_DIRECT_GLOBAL_PROB=0.12`,
  `DT_DIRECT_INITIAL_SCALE=0.18`
- global-heavy proposal: `DT_DIRECT_GLOBAL_PROB=0.30`,
  `DT_DIRECT_INITIAL_SCALE=0.30`

Pre-launch checks:

- local and Torch target-count checks both returned `15` base targets,
  `17` direct targets, and `15` direct parameters
- active model and cluster launcher checksums matched between local `main` and
  Torch scratch
- one-evaluation Slurm smoke jobs completed for both room specs:
  `py_direct_reduced_targets_hR8_smoke_20260527` and
  `py_direct_reduced_targets_hR6_smoke_20260527`
- the obsolete old-target array `9662249` was canceled before launch

## Latest Cluster Search

Active household-head ownership / small-owner-ladder outside-option search:

- Slurm job: `9662249`
- Run tag: `py_direct_outside_headown_smallown_global_4h_20260526`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_headown_smallown_global_4h_20260526`
- Purpose: first broad recalibration after correcting ownership targets to
  household heads/reference persons and allowing small owner units while
  preserving the renter cap.
- setup: `benchmark`
- bounds: `global`
- workers: `40` on Torch `cpu_short`, submitted as `1-40%32`
- internal worker budget: `13,500` seconds; Slurm wall time `03:55:00`
- population closure: `outside_option_benchmark_normalized`
- owner ladder: `H_own=[2.0,4.0,6.0,8.0,9.5,11.0]`
- renter cap: `hR_max=8.0`
- Pre-launch checks:
  - local one-evaluation direct-worker smoke completed
  - `check_population_closure.py` passed
  - Torch two-worker Slurm smoke job `9662202` completed and wrote
    `config.json`, `evaluations.jsonl`, `best.json`, and `status.json`
- Initial production check: `32` worker directories wrote configs/status/best
  records within the first minute; tasks `33--40` were pending on the array
  concurrency limit.
- On 2026-05-26 evening, while this run was still active, the array throttle
  was lowered to `24` and the eight weakest active workers by current best loss
  were canceled to free CPU slots for the `hR_max=6` diagnostic below. Their
  partial JSON outputs remain in the results directory.
- Snapshot used for the current fit note, pulled at `2026-05-26 22:45 EDT`:
  worker `29`, evaluation `440`, loss `19.662`, \(TFR=1.772\),
  ownership `0.313`, childless-renter median rooms `5.974`, center share
  `0.449`, and rent ratio `1.064`. The run was still active when this
  snapshot was pulled.
- Current fit note:
  `latex/current_fit_slide_diagnostics_20260526.pdf`; generated source is
  `code/model/tools/make_direct_fit_slide_note.py`.

Parallel renter-cap diagnostic:

- Slurm smoke job: `9664816`
- Slurm diagnostic job: `9664883`
- Run tag: `py_direct_outside_headown_smallown_hr6_30m_20260526`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_headown_smallown_hr6_30m_20260526`
- Purpose: test whether lowering only the renter cap from `hR_max=8.0` to
  `hR_max=6.0` helps the childless-renter room moment and prime-age ownership
  without treating it as a benchmark change.
- setup: `benchmark`
- bounds: `global`
- requested workers: `16` on Torch `cpu_short`, submitted as `1-16%16`
- active split at launch: `24` running workers remain on the `hR_max=8`
  production run and `8` workers started on the `hR_max=6` diagnostic; tasks
  `9--16` were pending on the user CPU limit.
- internal worker budget: `1,800` seconds; Slurm wall time `00:40:00`
- population closure: `outside_option_benchmark_normalized`
- owner ladder unchanged:
  `H_own=[2.0,4.0,6.0,8.0,9.5,11.0]`
- renter cap override: `hR_max=6.0`
- Smoke job `9664816` completed and verified written configs with
  `hR_max=6.0`.
- Final collection succeeded with `16` worker directories. Best diagnostic
  point: worker `12`, evaluation `59`, loss `55.106`, \(TFR=1.931\),
  ownership `0.629`, childless-renter median rooms `6.000`, old-age ownership
  `0.985`, center share `0.413`, and rent ratio `1.033`.
- Read: `hR_max=6` mechanically binds the childless-renter room moment but
  looks too distortionary relative to the live `hR_max=8` search: ownership is
  pushed high, late-life ownership is nearly universal, and the lifecycle
  ownership slope is too steep.

Corrected DUE-common-support diagnostic pulse:

- Slurm job: `9546318`
- Run tag: `py_direct_outside_commonH_owner0_tfrlt2_20m_20260525_003535`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_commonH_owner0_tfrlt2_20m_20260525_003535`
- Purpose: diagnostic run with the owner ladder lowered to zero and rental
  access extended to the owner maximum:
  `H_own=[0.0, 2.2, 4.4, 6.6, 8.8, 11.0]`,
  `hR_max=max(H_own)=11.0`.
- Existing property tax `tau_H=0.01`, transaction cost `psi=0.06`, and
  financed share `phi=0.80` remained active.
- setup: `benchmark`
- bounds: `global`
- workers: `32` on Torch `cpu_short`
- internal worker budget: `1,200` seconds; Slurm wall time `00:30:00`
- population closure: `outside_option_benchmark_normalized`
- launcher city-entry probability: `q^{E,*}=0.89`
- hard TFR cap: candidates with \(TFR\ge 2.0\) receive loss `1e6`.
- Smoke job `9546293` completed first with
  `DT_DIRECT_H_OWN_MIN=0` and `DT_DIRECT_HR_MAX=owner_max`; configs verified
  the lowered owner ladder and `hR_max=11.0`.
- Final collection succeeded with `32` worker directories, `29--49`
  evaluations per worker, and empty `squeue` for job `9546318`.
- Partial best: worker `26`, evaluation `24`, loss `82.2035`,
  \(TFR=1.893\), childless rate `0.074`, mean age first birth `29.68`,
  ownership `0.194`, childless-renter median rooms `6.872`,
  childless-owner median rooms `6.6`, center share `0.437`, rent ratio `1.538`.
- Read: the corrected common-support/owner-floor-zero diagnostic still does
  not repair the small-renter-unit miss. Lowering the owner floor reduces the
  owner-room median but leaves childless renters around `6.9` rooms and
  ownership far below target.

Earlier common-support diagnostic with owner grid left fixed:

- Slurm job: `9543538`
- Run tag: `py_direct_outside_commonH_tfrlt2_global_6h_20260525_000050`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_commonH_tfrlt2_global_6h_20260525_000050`
- Purpose: diagnostic run with perfect overlap between the renter support and
  owner support, leaving the existing owner grid fixed:
  `H_own=[4.0, 5.4, 6.8, 8.2, 9.6, 11.0]`,
  `hR_max=max(H_own)=11.0`.
- This does not add a new economic mechanism. The existing property tax
  `tau_H=0.01`, transaction cost `psi=0.06`, and financed share
  `phi=0.80` remain active.
- setup: `benchmark`
- bounds: `global`
- workers: `32` on Torch `cpu_short`
- wall time: `6:00:00`
- internal worker budget: `20,700` seconds
- population closure: `outside_option_benchmark_normalized`
- launcher city-entry probability: `q^{E,*}=0.89`
- hard TFR cap: candidates with \(TFR\ge 2.0\) receive loss `1e6`.
- Cluster smoke job `9543502` completed first with
  `DT_DIRECT_HR_MAX=owner_max`; configs verified `hR_max=11.0`.
- Stopped manually with `scancel` after about 20 minutes total runtime, by
  design. This was a diagnostic pulse, not a production calibration.
- Final partial collection succeeded with `32` worker directories, `15--36`
  evaluations per worker, and empty `squeue` for job `9543538`.
- Partial best: worker `7`, evaluation `13`, loss `83.5964`, \(TFR=1.902\),
  ownership `0.441`, childless-renter median rooms `6.700`,
  childless-owner median rooms `8.2`, center share `0.223`, rent ratio `1.143`.
- Read: common rental/owner support does not mechanically repair the
  small-rental-unit miss. Even with `hR_max=11.0`, the best diagnostic point
  keeps childless renters around `6.7` rooms and has much weaker ownership and
  geography than the previous `hR_max=8` outside-option search.

Active outside-option recalibration:

- Slurm job: `9516170`
- Run tag: `py_direct_outside_sout169_tfrlt2_global_6h_20260524_151347`
- Results directory:
  `/scratch/td2248/projects/Fertility_Spring26/code/cluster/results_python_direct_geometry_py_direct_outside_sout169_tfrlt2_global_6h_20260524_151347`
- setup: `benchmark`
- bounds: `global`
- workers: `32` on Torch `cpu_short`
- wall time: `6:00:00`
- internal worker budget: `20,700` seconds
- population closure: `outside_option_benchmark_normalized`
- empirical outside-origin entrant-share normalization:
  \(s^{E,\mathrm{out}}=0.169\)
- launcher city-entry probability: \(q^{E,*}=0.89\)
- hard TFR cap: candidates with \(TFR\ge 2.0\) receive loss `1e6` and cannot
  be selected as best.
- first collection succeeded with `32` worker directories; initial best was
  worker `1`, evaluation `1`, loss `353.731`, \(TFR=0.891\),
  ownership `0.121`, \(N=1.000\), \(pop_C=0.393\), rent ratio `1.276`.
- The prior uncapped run
  `py_direct_outside_sout169_global_6h_20260524_145525` was canceled after the
  best-so-far had \(TFR>2\); use the capped run above for live results.

The latest completed cluster search below was run under the prior
`renewal_valve_calibrated` closure. It is useful as a pre-outside-option
benchmark, but it is not a recalibration of the current live closure.

Results directory:

- `code/cluster/results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506`

Torch arrays:

- `8113159`, `8113160`, `8113161`
- `120` total tasks on `cpu_short`
- all tasks completed with exit code `0:0`

Search settings:

- setup: `benchmark`
- bounds: `global`
- parameters: `15`
- `geo_weight = 100`
- `renewal_retention = 1.0`
- `scale_target = 1.0` imposed mechanically
- per-stage time budget: `13,500` seconds

Best collected candidate:

- worker `37`, evaluation `1425`
- loss `8.1275`
- solve time `44.2` seconds
- GE accepted, but not strict: final equilibrium error `0.003006`,
  convergence reason `soft_tol_10x`
- prices: \(p_P=0.49766\), \(p_C=0.58848\)

Closure diagnostics:

| Diagnostic | Value |
|---|---:|
| `entry_per_unit_scale` | `0.0166667` |
| `mature_cityborn_per_unit_scale` | `0.0155334` |
| `outside_entry_flow` | `0.0011333` |
| `population_scale_denominator` | `0.0011333` |
| `implied_total_population` | `1.000` |
| `scale_factor` | `1.000` |
| `stationary_entry_relative_residual` | `0` |

## Target Fit

Targets come from `build_calibration_setup()` in the active Python
`parameters.py`. The model column below is from the latest completed
pre-ownership-correction renewal-valve run, so it is a historical benchmark
until a new outside-option calibration is collected. The last two rows are
direct-geometry targets disciplined by the geometry-weight block.

| Moment | Target | Model |
|---|---:|---:|
| `tfr` | `1.700` | `1.898` |
| `childless_rate` | `0.150` | `0.145` |
| `mean_age_first_birth` | `26.000` | `33.535` |
| `tfr_gradient` | `0.133` | `0.119` |
| `own_rate` | `0.575` | `0.643` |
| `own_gradient` | `0.186` | `0.139` |
| `own_family_gap` | `0.168` | `0.114` |
| `housing_increment_0to1` | `0.664` | `0.441` |
| `housing_increment_1to2` | `0.566` | `0.192` |
| `young_liquid_wealth_to_income` | `0.600` | `0.527` |
| `center_share_nonparents` | `0.494` | `0.405` |
| `center_share_newparents` | `0.416` | `0.382` |
| `migration_rate` | `0.032` | `0.035` |
| `old_age_own_rate` | `0.764` | `0.947` |
| `old_age_parent_childless_gap` | `0.070` | `0.062` |
| `inv_pop_share_C` | `0.450` | `0.441` |
| `inv_rent_ratio_C_over_P` | `1.140` | `1.182` |

Current read:

- Under the prior renewal-valve search, benchmark scale was imposed
  mechanically, the outside flow was positive, and the stationary entry
  residual was zero.
- Remaining misses are concentrated in first-birth timing, childless renter
  baseline rooms, the first- and second-birth housing responses, and old-age
  ownership.
- Do not compare this loss mechanically to pre-Python, pre-renewal-valve, or
  MATLAB-inversion losses.

## Reporting Rules

When reporting calibration results, always show targets next to model moments.
Do not report only losses or only model moments when target values are
available.

Do not treat `parity_progression_1to2` as a calibration target unless the
fertility architecture has been changed to support a sequential second-birth
hazard.

Use `CALIBRATION_STATUS.md` first. Historical notes are archived here:

- `docs/archive/calibration_docs_2026-05-07/`
- `docs/archive/root_notes_2026-05-07/`
- `calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/status_notes/`
