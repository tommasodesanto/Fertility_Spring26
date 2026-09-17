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

Twelve mappings, 21638 s, not certified (`evaluation_budget`), final replay
exact. Full and trimmed normalized scores and raw max-abs residual by mapping:

| mapping | 1 | 2 | 3 | 4 | 5 | 6\* | 7 | 8\* | 9\* | 10\* | 11 | final |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| full score | 13.2 | 37.3 | 20.3 | 17.6 | 13.8 | 42.9 | 12.0 | 32.9 | 28.3 | 27.3 | 11.5 | 11.5 |
| trimmed score | 9.41 | 5.26 | 2.40 | 1.56 | 2.11 | 2.38 | 1.07 | 3.65 | 3.33 | 3.16 | 1.03 | 1.03 |
| raw max-abs | 0.447 | 1.518 | 0.829 | 0.694 | 0.552 | 1.694 | 0.475 | 1.291 | 1.113 | 1.074 | 0.456 | 0.456 |

\* safeguard fired (restored best, damping halved, Jacobian reset); final
damping 0.0625.

At the accepted best (mapping 11) every PAYGO and rebate row is inside the
0.5 gate, including the threshold dates (rebate \(-0.456\) at date 70, 0.91
of tolerance; \(-0.242\) at date 44). The coordinates that fail are the
**housing** rows at the threshold dates: relative imbalance
\(2.3\times10^{-3}\) at date 70, \(1.7\times10^{-3}\) at 71 and
\(1.2\times10^{-3}\) at 44 (11, 8 and 6 times the unchanged
\(2\times10^{-4}\) housing gate), plus 1.0--1.1 times the gate at dates 45
and 98--99. Every other housing residual is numerically zero. The tenure flip
therefore moves housing demand by about 0.2 percent at the dates where it
occurs, so the discreteness floor is in the housing rows as well as the
fiscal rows, at the size of one flip. Loosening only the fiscal gate is not
sufficient at the frozen smoothing scale; at 0.05 the flips shrink three- to
four-fold (Section 1), which would bring these housing residuals to roughly
\(6\times10^{-4}\), and at 0.1 to about the gate.

## 5. Prepared, not launched: the smoothed-scale transition test

Author-approved on September 16 for launch when Torch is available. One
command from the repository root:

```bash
code/cluster/submit_e5f_ssj_smoothed_transition.sh 0.05
```

It stages `code/cluster/run_e5f_ssj_smoothed_transition.py` (with the copied
solver and rescue helpers) in a new batch
`announced_original_queue_20260913c_ssj_smoothed_k0.05_<tag>`, writes the
pinned manifest and sbatch through
`code/cluster/prepare_e5f_ssj_smoothed_transition.py`, submits, and mirrors
the contracts locally. The job then:

1. re-solves the stationary equilibrium at the frozen scale (control) and at
   0.05 with the frozen `solve_terminal`, and writes `fit_table.md`
   (stationary moments side by side; the full SMM target table still needs
   the calibration collector);
2. re-solves the terminal at 0.05 and the final announced preference, warm
   started from the frozen endpoint coordinates;
3. runs the 104-date announced root with the inherited 2007 state replaced by
   the 0.05 stationary state (removes the date-1 artifact), the 0.05 terminal,
   the saved 104-date checkpoint as start, the extrapolated ten-date Toeplitz
   `initial_jacobian`, the direction-preserving step, a 12-mapping budget, no
   trimming, and the **retained** gates (max-abs \(2\times10^{-4}\)), so a
   certification would mean the model meets its own tolerance.

Budget: 31000 s numerical deadline, 540-minute Slurm cap, 3600 s per
stationary solve. Every deviation is written into the job's
`experiment_contract.json`. The announced preference path stays as fitted at
0.005, which is a stated approximation. Pass a different scale as the first
argument to test 0.1.

## 6. September 16 evening: smoothing-scale experiments (author-approved, copies only)

**0.05 transition run (jobs 17904418, 17904647).** The probe-scale stationary
solve from the 0.005 coordinates exhausted 16 evaluations (price 0.703 to
0.576, rebate 0.187 to 0.140, residual still 0.029); relaunched from that
iterate with the 24-evaluation budget it verified in 6 evaluations, and the
0.05 terminal verified in 17. Fixed-theta stationary comparison, 0.05 versus
0.005: price \(-17.9\) percent, housing stock \(-11.7\) percent, ownership
0.566 to 0.600, entry flow \(-3.3\) percent. The 104-date root at 0.05 was
running at the time of writing.

**Track A recalibration panel at 0.05 (job 17905502; first attempt 17905057
failed on the live panel size of 21).** Center: the September 4 baseline best
(loss 30.48 at 0.005). Anchor at the same theta with scale 0.05: loss 379.17,
of which 272.7 is the family ownership gap collapsing from 0.161 to 0.029
(target 0.168). Best single-coordinate move: 316.86 (child room floor down),
gap 0.031. The only move that restores the gap (fertility dispersion
2.17 to 0.53, gap 0.1665) destroys the fertility timing moments (loss 1436).
Mean rooms and wealth moments improve at 0.05; the room-response moments
worsen. Reading: at 0.05 the tenure taste noise is comparable to the utility
stakes that sort parents into ownership, so the mechanism the model is built
around is smoothed away. A micro recalibration cannot recover it; 0.05 is not
a usable scale. Anchors at 0.02 and 0.01 (jobs 17907117, 17907118) were
launched to locate where the gap survives. Packet:
`output/model/e5f_sequence_space_prototype_20260913/track_a_recal_kappa050/`.
Track B (the frozen packet's objective) was sized but not launched: it needs
a re-derived objective contract because the scale is pinned inside the
template's fixed-parameter source.

**Anchors at 0.01 and 0.02 (jobs 17907118, 17907117), same theta and objective:**

| scale | loss | family ownership gap (target 0.168) |
|---|---|---|
| 0.005 (frozen) | 31.85 | 0.161 |
| 0.01 | 101.81 | 0.121 |
| 0.02 | 299.31 | 0.062 |
| 0.05 | 379.17 | 0.029 |

The gap falls continuously and steeply with the tenure taste scale; there is
no value that removes the threshold flips while keeping the mechanism. The
smoothing route to a certifiable long transition is therefore closed on
economic grounds. What remains: (i) gates commensurate with one flip at the
frozen scale, or (ii) a numerical treatment of the flip that does not add
taste noise, namely interpolating the tenure indifference point inside a
wealth-grid cell so the mass moves continuously with the price without
changing any household's preferences.

**Own-rent stakes from the saved stationary solution (job 17907789,
`stakes_diagnostic_20260916a/`).** With the six-product logit at
\(\kappa_H=0.005\), \(\Delta V=\kappa_H[\log\sum_{\text{own}}p_j-\log p_{\text{rent}}]\)
is recoverable from the saved float32 tenure probabilities wherever they are
not saturated (\(|\Delta V|<0.44\)); 24 percent of the mass is saturated and
firmly settled. State-weighted by the pre-decision distribution:

| group | mass | own share | \(p_{25}\) | median \(\Delta V\) | \(p_{75}\) | share with \(|\Delta V|<0.02\) | \(<0.05\) |
|---|---|---|---|---|---|---|---|
| all | 1.00 | 0.555 | \(-0.002\) | 0.028 | 0.042 | 0.16 | 0.54 |
| childless (\(n=0\)) | 0.31 | 0.136 | \(-0.054\) | \(-0.024\) | 0.000 | 0.16 | 0.39 |
| parents (\(n\ge1\)) | 0.69 | 0.740 | 0.019 | 0.032 | 0.048 | 0.16 | 0.61 |
| current renters | 0.49 | 0.094 | \(-0.057\) | \(-0.024\) | \(-0.002\) | 0.22 | 0.35 |
| current owners | 0.51 | 0.993 | 0.027 | 0.035 | 0.053 | 0.10 | 0.73 |

The typical stake is a few hundredths of a utility unit: six times the frozen
noise scale, so tenure is effectively deterministic, but below a scale of
0.05, at which more than half of households would be mixed toward 50/50
whatever their family status. This is the mechanism behind the loss table
above. It also says where the flips come from: the 5 percent of households
within 0.005 of indifference are the grid nodes that switch when the price
path moves, and their mass is concentrated among current renters with
children (16 percent of that group within 0.005).

## 7. Overnight September 16--17: a receipted 104-date path under declared gates

Job 17907712 (batch `announced_original_queue_20260913c_ssj_declared_gate_20260916a`,
local packet `announced_declared_gate_20260916a/`). Warm start from the
exactly replayed best of job 17865027; measured-Jacobian start; scaled step;
no trimming; declared gates housing \(2.5\times10^{-3}\), fiscal 0.5 scaled.
The initial mapping reproduced the warm start exactly (raw 0.4557, normalized
0.919) and the final replay certified it: `converged=True`, reproduction gap
0, two mappings, 3857 s including per-mapping plots. This is a finite-horizon
market/fiscal equilibrium of the announced four-shock transition at the frozen
tenure scale **up to one tenure flip**, documented as such:

- housing residual above the retained \(2\times10^{-4}\) only at dates 44
  (\(-1.2\times10^{-3}\)), 70 (\(-2.3\times10^{-3}\)), 71
  (\(+1.7\times10^{-3}\)) and marginally 98--99 (\(-2.2\times10^{-4}\));
- rebate residual above 0.1 only at dates 44 (\(-0.24\)), 70 (\(-0.46\)),
  71 (\(+0.34\)); PAYGO at most \(1.3\times10^{-2}\);
- every other date meets the retained gates;
- terminal-distance diagnostic unchanged in kind (population 10 percent
  above the endpoint at date 104), i.e. the terminal approach is unverified
  exactly as in the announced contract.

Accepted price path (every 13th date): 0.695, 0.601, 0.471, 0.388, 0.334,
0.300, 0.280, 0.264; rebate path 0.183 to 0.096. Full rows, fertility path,
figures and the standard diagnostics are in the batch. The receipt states
its gates explicitly; it is not a certification under the retained contract
and is not production eligible.

**0.05 transition root close-out (job 17904647).** Twelve mappings, 20719 s,
not certified (`evaluation_budget`), final replay exact. Started from the
0.005 checkpoint prices (about 18 percent too high for the 0.05 economy), so
the initial residual was 71, concentrated at the last dates; the
measured-Jacobian start (extrapolated from the 0.005 economy) brought it to
0.370 by mapping 10 (71, 70, 59, 41, 19, 4.3, 2.9, 1.17, 0.60, 0.37), then
mapping 11 tripped the safeguard. At the best point the residual above the
retained gates sits at the unconverged tail (dates 98--103, housing up to
\(1.9\times10^{-3}\), rebate 0.37) and at one remaining flip at date 34
(housing \(-1.1\times10^{-3}\), rebate \(-0.20\)); the 0.005-economy flips at
dates 44 and 70--71 are gone, as the fixed-path probe predicted. So at 0.05
the flips are reduced, not eliminated, and the budget rather than a floor
ended the run. Given the Track A result the scale is not usable, so no
continuation was launched. Terminal distance: population 21 percent above the
0.05 endpoint at date 104 (terminal approach unverified, as in every arm).

## 8. Status at the close of September 17 morning

No job is running. Nothing in production changed. The ledger carries active
point 11 (WAIT) with the evidence trail; the design note for the interpolated
tenure boundary is `docs/model/e5f_tenure_boundary_interpolation_design_20260916.md`.
