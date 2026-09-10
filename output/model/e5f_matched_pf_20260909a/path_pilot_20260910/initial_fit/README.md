# Initial-economy diagnostic: ready-to-run readout, no new SMM

**The smallest informative run loads the saved normalized old stationary economy and evaluates existing observers. It needs zero stationary GE solves and zero Bellman solves.** It can immediately show which inherited parameters already approximate the early observations, while leaving missing or mismatched observation rules visible. It cannot honestly produce a new 13-row calibration objective yet.

The complete proposed restriction-to-observer map is `initial_restriction_observer_map.csv`: all 13 rows, every data value, available uncertainty, source, exact observer function, approximation and weighting status. The initial level is the author's chosen **2.1**, not the superseded 2.0605 proposal. `inherited_parameters.csv` contains all eleven inherited search coordinates and four normalized/derived/fixed rows, estimates and bounds. No parameter was re-estimated.

## What can be read now

1. **Existing observers:** annual wealth/gross earnings and bequest flow; old wealth dispersion; ownership; first-birth period timing; stationary first-birth housing response. These are computationally available but retain the listed age, income, geographic or time-stability approximations.
2. **Proxy only:** the existing completed-parity 2.1 normalization, age-42 childlessness, uncapped rooms and old family-group comparisons. These are not reported as exact matches to female exposure, CPS ages40–44, capped rooms or the new resident-family groups. The diagnostic writes proxy values separately and leaves their fit gaps blank.
3. **No existing named observer:** exactly-one-child among mothers ages40–44. The underlying parity masses are saved, so an explicitly tested age/pool observer is straightforward; the readout does not quietly substitute pooled terminal parity.

`read_initial_baseline.py` loads the checkpoint only after validating the complete historical source/input contract. It evaluates the supplied **pre-announcement stationary policy on `old.stationary_g_pre`**, uses the stationary supply primitive to report its residual, and verifies zero solver calls. It does not use the separately age-reweighted `old.initial_state.g_pre` as if that were the stationary distribution. It runs `period_fertility_diagnostics`, `first_birth_housing_response` and `transition_cross_section_moments`, then writes all thirteen rows with explicit comparison status. It copies the complete inherited parameter table only after checking every value against the loaded state. It never constructs weights or a scalar objective.

## Exact cache and source

The immutable contract is `output/model/e5f_matched_pf_20260909a/horizon100_root_contract.json`; its SHA is pinned in `manifest.json`.

Normalized checkpoint:

```
/scratch/td2248/projects/Fertility_Spring26_matched_pf_normalized_20260909b/output/initial_01/sequential/normalized_old.pkl.gz
```

SHA-256: `dfeb34f077e7a16804119eabd43cb8ee0611b0c3bb0180ec8af8ad7dc51b94e4`.

The pickle stores the keys `old`, `demographics`, `contract_sha256`. `old` is `e5f_matched_pf_initial_state.NormalizedOldState`, containing `parameters`, `b_grid`, `solution`, `policy`, `shared`, `stationary_g_pre`, `initial_state`, `historical_conditioning`, `supply_rule`, `years`, `psi_path`, and `diagnostics`. This is verified from the producer and loader source, **not by loading large arrays during this task**. The local receipt directory contains only the small contract/summary/heartbeat JSONs, not the checkpoint.

`run_e5f_matched_pf_baseline.load_normalized` verifies checkpoint, summary, parent-contract hashes and selected-input/arm identity before returning these objects. Use source F at `/scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a`, or an exactly matching immutable checkout.

## Execution

Stage this small diagnostic folder and the existing `horizon100_root_contract.json` on Torch. Then run, with a fresh output directory:

```sh
python /path/to/initial_fit/read_initial_baseline.py \
  --source-root /scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a \
  --contract /path/to/horizon100_root_contract.json \
  --output /path/to/fresh_initial_readout
```

The `--validate-only` option checks the 13-row ledger, parameter receipt, complete source fingerprint and input contract without importing model modules or loading the checkpoint. The real readout needs only one checkpoint load, a forward stationary-policy evaluation and existing statistics/branch observers. No empirical data reread, equilibrium search or anticipated historical path is needed. Runtime and peak memory of this observer-only readout have not been measured; use one core and the existing16GB checkpoint allocation with a bounded15-minute timeout for the first execution. A failure should be diagnosed, not converted into missing rows silently.

Expected outputs: `all_13_initial_rows.csv`, `parameters.csv`, `supplemental_existing_observers.json`, and `validation.json`. No charts are generated by this diagnostic; no user-requested graph redesign is implied.

## Genuine prerequisites to a calibration loop

- **Units of 2.1:** the source normalizes `extract_moments(...)[tfr]`, which is topcode-adjusted completed parity over post-fertile model ages. The period helper divides birth flows by adult-household age mass, not female person exposure. The author selected the initial scalar; the observational and renewal meanings still need to be documented consistently. Retain existing2.1 normalization and entry gates during the initial diagnostic.
- **CPS age/cohort:** the current dated stock helper rounds42 to the model cell labelled42, covering42–45 if interpreted as four-year spans. CPS is women40–44. A overlap-weighted40–44 observer would use parts of cells38 and42, with an explicit within-cell assumption. The stationary initial distribution cannot recover the older CPS cohorts' actual histories. Exactly-one conditional motherhood is independent information; mean fertility conditional on motherhood is algebraically redundant with unconditional mean and childlessness.
- **Room censoring:** use `min(h,9)` before aggregation for the new ACS room rows, on income-resolved realized states. Capping the average or income-collapsed renter policy is wrong. The reviewed PSID birth response remains uncapped.
- **Family groups:** independent-count child states record the number of dependents, not their exact ages. Existing `own_family_gap` compares any dependent-child parent with never-parent, whereas the new ACS row uses recent parents (eldest child<4) versus no resident own child. A flow/first-birth branch can help form a recent-parent proxy, but exact equivalence and age treatment are not supplied by the current stock observer. The family room gap compares model dependents with empirical resident own children under an under18 sample signal; disclose or repair that mapping.
- **Age boundaries:** `age_to_index` rounds. Current ownership30–55 includes model starts30 through54, spanning30–57; young25–34 includes26 through34, spanning26–37. Old-tail code directly tests starts76–84, selecting78 and82, spanning78–85. These are existing approximations, not new defects caused by the early targets. The diagnostic exposes them rather than labels every row exact.
- **Population weighting:** the normalized old object includes both a stationary demographic distribution and an ACS2007 age-reweighted inherited distribution. Choosing which supplies the initial calibration moments matters for levels and must be explicit; reweighting cannot be silently treated as unchanged stationary equilibrium.
- **Weights and identification:** CPS/timing uncertainty or justified scales, a new source/target/weight fingerprint, and an informative local parameter Jacobian are still required. Available rows and parameter counts do not establish identification. Existing later objective rows are not removed or reweighted by this diagnostic.
- **Candidate plumbing:** `initialize_normalized_old_state` checks the selected parameter coordinates and allows only `psi_child` inside its normalization adapter. The existing matched observer hard-pins the old12-row transition contract. A new initial candidate loop must explicitly accept a ten-coordinate proposal, rebuild and renormalize its initial economy, observe the new contract, and preserve numerical/accounting gates; it is not an existing CLI switch. The old supply elasticity1.75 versus dated.63 role remains explicit.

## Validation completed and limits

Python syntax and ledger checks pass. All503 source fingerprints match the immutable worktree HEAD. Current working-tree validation correctly rejects one concurrent change in `run_e5f_matched_pf_baseline.py`, made by the preference-path pilot agent; it does not alter this readout's source pins. `local_validation.json` records that expected rejection and the matching immutable HEAD hash. The script can run against original remote F. Alternatively, a newly reviewed pilot contract must explicitly pin its changed source; the guard must not be bypassed.

No checkpoint, large array, model observer or model solve was executed here. No jobs submitted. This evidence packet is an executable diagnostic preparation, not completed initial fit or a new calibration.
