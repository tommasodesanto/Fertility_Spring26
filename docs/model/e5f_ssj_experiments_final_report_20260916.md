# Sequence-space and transition-root experiments: final report

September 16, 2026. All work below is experimental. No production solver,
kernel, gate, parameter, population law, or residual definition was changed.
Every numerical run lives in an isolated Torch batch under
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/`
with pinned input hashes, and its local mirror under
`output/model/e5f_sequence_space_prototype_20260913/`. The running narrative
with all intermediate evidence is `docs/model/e5f_sequence_space_prototype.md`.

## 1. Questions and answers

**Why does the transition root stall?** Because the retained Broyden solver
starts from a diagonal Jacobian and the true one is strongly non-diagonal.
Measured at the 2007 stationary economy, dated housing demand responds to the
next period's price (\(+1.03\)) almost as strongly as to the current one
(\(-1.91\)); the rebate row inherits the same structure (\(+205\) and
\(-57\) in the \(\times200\) fiscal units). With four-year periods the user
cost is a small difference between \(q_t\) and discounted \(q_{t+1}\), so
demand reacts to the slope of the price path. Seven rank-one updates cannot
learn that. Verified by a predictive check along the reference path: the
measured matrix predicts realized residual changes with relative error
0.13--0.45 against 0.43--1.38 for the diagonal start.

**Does a measured initial Jacobian fix it?** Yes for the smooth problem. On
the ten-period shocked root, identical start, controls, endpoint and
eight-mapping budget, the best score fell from 1.346 (diagonal start, job
17700926) to 0.02246 (job 17714834). Continued for one more eight-mapping
run from the exactly reproduced best, the unmodified retained solver
certified the root at \(2.69\times10^{-5}\) (job 17730694); a copy with a
direction-preserving step certified at \(9.59\times10^{-5}\) (job 17730695).
The two certified roots agree to \(4\times10^{-7}\) relative. Cost: seven
derivative mappings once (722 s) plus 13--14 root mappings, versus eight
unconverged mappings before.

**Why does the 104-date announced root not converge?** Not the Jacobian. At
its best iterate (0.977 after eight mappings, job 17711519) the rebate
residual sits at isolated dates (43--44, 68--71) plus a smooth tail, and the
native rows show ownership jumping 1.2--1.7 percentage points between
adjacent dates there; a log-price step of 0.001 moves the jump to the next
date. The frozen tenure-smoothing scale is 0.005, the June 28 lower search
bound (issue M19), numerically an argmax: whole grid nodes flip tenure
together, the property-tax base moves by about 0.1 percent, and the rebate
row jumps by 0.2--1.5 against a gate of \(2\times10^{-4}\) (\(10^{-6}\)
relative). No smooth root exists across a straddled threshold at that
tolerance; over 104 dates a few dates always straddle one. Rescue arms
(jobs 17732268, 17732269) and the loosened-gate arms (17860152, 17865027)
all confirmed this: the extrapolated Jacobian drives the smooth residual to
zero in one step (housing from \(2.6\times10^{-3}\) to \(10^{-16}\) away from
the threshold dates), the threshold pair then dominates, and the safeguard
eventually discards the Jacobian. Under a trimmed-score acceptance rule the
smooth part reached a trimmed score of 1.56 (housing tail within 1.6 times its
gate) with raw residual 0.55 before the threshold coordinates again tripped
the safeguard at mapping 6.

**Does more tenure smoothing remove the flips?** Yes. One fixed-path mapping
per scale (jobs 17866783, 17866789, 17866790, 17866791; the 0.005 control
reproduced the checkpoint residual to all digits):

| scale | jump at date 25 | date 44 | date 71 | typical date-to-date change |
|---|---|---|---|---|
| 0.005 (frozen) | 1.25 pp | 1.33 pp | 1.67 pp | 0.16 pp |
| 0.02 | 0.70 | 0.81 | 0.98 | 0.18 |
| 0.05 | 0.45 | 0.22 | 0.50 | 0.11 |
| 0.1 | 0.33 | 0.09 | 0.26 | 0.06 |

At 0.05 the profile is smooth apart from a 2.5-point re-sorting at date 1,
which is an artifact of the 2007 initial state having been solved at 0.005.
The ownership level on the fixed path is non-monotone in the scale and must
be read from a re-solved stationary equilibrium, not from this probe.

## 2. Complete inventory of changes

### 2.1 New files (all committed on `main`; none imported by production code)

| file | purpose |
|---|---|
| `code/model/tools/e5f_ssj_toeplitz_jacobian.py` | pure-numpy assembly of a block-Toeplitz Jacobian from measured lag profiles; horizon extension |
| `code/cluster/run_e5f_ssj_toeplitz_jacobian_root.py` | two-stage driver: seven-mapping derivative stage, then the retained ten-period root with `initial_jacobian`; options `jacobian_source`, `step_rule`, `warm_start_receipt` |
| `code/model/tools/e5f_ssj_scaled_step_root.py` | copy of `tmp/e5f_matched_pf/.../e5f_matched_pf_path_root.solve_price_path` with three diagnostic options (below) |
| `code/cluster/run_e5f_ssj_announced_rescue.py` | warm restart of the 104-date announced root through the announced batch's own `run_path`; options `jacobian_mode`, `step_rule`, `warm_start_checkpoint`, `mapping_budget`, `housing_gate`, `fiscal_gate_scaled`, `trim_count`, `skip_mapping_plots` |
| `code/cluster/run_e5f_ssj_tenure_smoothing_probe.py` | one fixed-path mapping with `tenure_choice_kappa` overridden on a deep copy of the initial-state parameters |
| `code/model/tools/collect_e5f_ssj_toeplitz_jacobian.py` | packet summariser (`summary.md`, reference-path prediction check) |
| `code/model/tools/test_e5f_ssj_toeplitz_jacobian.py`, `test_run_e5f_ssj_toeplitz_jacobian_root.py`, `test_e5f_ssj_scaled_step_root.py`, `test_run_e5f_ssj_announced_rescue.py`, `test_run_e5f_ssj_tenure_smoothing_probe.py` | 24 pure tests, no native model |
| `code/model/tools/e5f_sequence_space_prototype.py`, `test_e5f_sequence_space_prototype.py`, `code/cluster/check_e5f_sequence_space_native.py` | the earlier worker's interface prototype and two-date native smoke (committed as found) |

### 2.2 Deviations from the retained root, by experiment

Every deviation is also written into each batch's `submission.json` and the
`experiment_contract.json` the job writes at startup.

| experiment | jobs | deviations from the retained root |
|---|---|---|
| Ten-period measured start | 17714834 | `initial_jacobian` = measured block-Toeplitz matrix (seven native mappings at the 2007 stationary economy, log step \(10^{-5}\), date 5); afternoon batch deadline (expired) replaced by the handoff's two-hour cap; otherwise identical |
| Ten-period scaled step | 17717654 | same Jacobian reused; solver copy with uniform step scaling to max coordinate 0.2 instead of componentwise clipping |
| Ten-period continuations | 17730694, 17730695 | warm start from each arm's exactly reproduced best with its learned Broyden matrix (frozen operator refuses more than eight mappings, so continuation replaces a larger budget) |
| 104-date rescue | 17732268 (cancelled as duplicate), 17732269 | warm start from the announced root's certified best; Jacobian = its final Broyden matrix (arm A, retained solver) or extrapolated ten-date Toeplitz (arm B, scaled step); own 16200 s deadline |
| 104-date loosened gate | 17860152 (cancelled after mapping 2) | warm start from the uncertified rescue checkpoint (0.447); extrapolated Toeplitz Jacobian; scaled step; budget 12 (frozen validation still receives 8); fiscal rows gated at 0.5 scaled instead of \(2\times10^{-4}\), housing rows unchanged, via a per-coordinate tolerance vector in the solver copy; per-mapping plots skipped |
| 104-date loosened gate, trimmed acceptance | 17865027 | as above plus: best-point selection and the 1.5 worsening safeguard use the fifth-largest normalized coordinate (`trim_count=4`); certification still requires every coordinate inside its gate |
| Tenure-smoothing probe | 17866783, 17866789, 17866790, 17866791 | `tenure_choice_kappa` in {0.005, 0.02, 0.05, 0.1} on the fixed warm-start path; no root; terminal continuation left at the 0.005 endpoint |

### 2.3 What is unchanged in every experiment

Native household kernels, population law, four-vintage birth queue at
births/2.1, no immigration, equal property-tax rebate, PAYGO at 0.179, all
mass, policy-reproduction and feasibility gates, the verified endpoint, the
preference path, the residual definitions and the \(\times200\) fiscal
scaling, and the final exact-replay certification.

### 2.4 Documentation and status

`docs/model/e5f_sequence_space_prototype.md` (full narrative),
`CALIBRATION_STATUS.md` (September 13 entry), `memory/AGENT_MEMORY.md` (SSJ
gotcha), `output/model/e5f_sequence_space_prototype_20260913/README.md`
(packet index). Commits df1e72f1 through 22f38204 plus this report.

## 3. Decisions that remain the author's

1. **Tenure-smoothing scale.** The frozen 0.005 is a bound, not an estimate.
   Around 0.05 the threshold flips that set the long-root floor disappear.
   Adopting it requires re-solving the stationary equilibrium and checking
   the calibration fit at that scale; nothing here measures the level effect.
2. **Fiscal gate for long roots.** If the scale stays at 0.005, the rebate
   rows cannot be balanced date by date at \(10^{-6}\) relative along 104
   dates; a gate near one threshold jump (0.5--2 scaled) is what the model can
   meet.
3. **Solver start.** A measured initial Jacobian, built once per horizon and
   reused across policy roots, is what certified the ten-period root. For 104
   dates the extrapolated ten-date matrix solves the smooth part in one step;
   a directly measured 104-date matrix costs about seven long mappings.
4. **Step rule and safeguard.** Uniform step scaling and a trimmed acceptance
   rule are both small changes to the solver copy; whether either belongs in
   production is a separate verification task.

## 4. Job 17865027 close-out

(filled in when the job ends)
