# Sequence-space prototype for the original household queue

September 13, 2026. This is a small, isolated adapter to the official
[`sequence-jacobian`](https://github.com/shade-econ/sequence-jacobian) nonlinear
path-update convention. It is not a package migration and does not alter any
scientific module or running transition job.

## Exact contract

For a horizon \(T\), the native root coordinates are

\[
 U = (\log q_0,\ldots,\log q_{T-1},\ b_0,\ldots,b_{T-1},\
 r_0,\ldots,r_{T-1})\in\mathbb R^{3T},
\]

where \(q_t\) is the house price and \(b_t,r_t\) are pension and equal
property-tax rebate levels. The exact native mapping returns

\[
 H(U)=(h_t,\ 200 f^{PAYGO}_t,\ 200 f^{rebate}_t)_{t=0}^{T-1}.
\]

Here \(h_t\) is the existing housing relative imbalance and both fiscal
entries are the existing relative imbalance times 200. The prototype does not
redefine any residual, scalar normalization, or convergence gate.

Each evaluation carries the native level-valued state
\((g^{pre}_t,e_t^{1:4},\tilde e_t^{1:4})\): the household distribution,
adjusted-birth entry queue, and raw-birth entry queue. The bridge rejects an
evaluation missing any of these objects. It neither normalizes mass to one nor
overwrites the carried state with the terminal stationary state. Births still
enter after the original four-vintage delay at births/2.1, with no immigration.

## What uses Sequence-Jacobian

`factor_ssj` builds the package's `JacobianDict` and `FactoredJacobianDict`
from an explicitly supplied \(3T\times3T\) native residual Jacobian. The
nonlinear update then is exactly the official package's signed step
\(\Delta U=-H_U^{-1}H(U)\), while the residual is freshly evaluated by the
native household/population operator. This follows the package's
`Block.solve_impulse_nonlinear` design, which factors a Jacobian once and
reuses it while re-evaluating nonlinear residuals. The current local Python
environment lacks the package. An attempted isolated install to
`/private/tmp/e5f_ssj_20260913` could not resolve PyPI because this machine has
no package-network DNS access; no global environment was changed. The adapter
therefore has pure local interface tests only.

A package factorization does **not** make a dense finite-difference Jacobian a
fast-news Jacobian. A complete fast implementation would require derivatives
of: (i) terminal and dated household policies including discrete feasibility
thresholds; (ii) the distribution forward operator; (iii) both queue update
rules and aggregate population scale; (iv) births/2.1 and housing demand; and
(v) PAYGO/rebate accounts and the housing-supply curve. Those objects were not
constructed here.

## Evidence

Run locally from `code/model/tools`:

```bash
/Users/tommasodesanto/miniconda3/bin/python -m unittest -v test_e5f_sequence_space_prototype.py
```

These tests are interface-only: a toy affine native evaluator confirms a
zero-residual path and a central directional derivative; another test proves
that a result missing either queue is rejected. They are not native model
checks and do not establish an equilibrium or speedup.

## Recoverable ten-period design

The ten-period native root has \(3\times10=30\) coordinates. A dense central
finite-difference Jacobian would need 60 perturbed mappings, plus a baseline,
which exceeds the scope and is not a credible acceleration experiment. The
current observed native ten-period mapping is about three minutes; serial dense
construction is therefore roughly 183 minutes before any nonlinear replay.
Parallelism cannot be assumed safe because each mapping uses large caches.

The next admissible native experiment is instead one separately named,
horizon-two directional smoke with at most six full evaluations: baseline,
two central directional evaluations, and up to three native replays. It uses
one CPU, 24 GiB, all numerical-library threads one, a 10-minute wall limit,
and the frozen `afternoon_original_queue_20260913a/spec.json` hash contract.
It preserves the native finite residual tolerance \(2\times10^{-4}\), fresh
replay tolerance \(2\times10^{-10}\), and all household/feasibility/mass gates.

Exact command (prepare/submit through the existing Torch batch mechanism; do
not execute on the login node):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 /share/apps/anaconda3/2025.06/bin/python code/model/tools/e5f_sequence_space_prototype.py --native-smoke --spec /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/afternoon_original_queue_20260913a/spec.json --horizon 2 --max-evaluations 6 --wall-minutes 10 --cpus 1 --mem-gib 24 --account torch_pr_570_general
```

This command records the bounded contract but deliberately does not submit a
job. A ten-period benchmark is worth running only if that smoke reproduces the
native zero-shock and directional checks and a subsequent explicit derivative
construction is shown to be cheaper than its mapping cost. No speed claim is
made by this prototype.

## September 13 continuation: measured block-Toeplitz initial Jacobian

Handoff `docs/prompts/HANDOFF_claude_ssj.md`. New isolated files:
`code/model/tools/e5f_ssj_toeplitz_jacobian.py` (pure numpy assembly),
`code/cluster/run_e5f_ssj_toeplitz_jacobian_root.py` (two-stage driver),
`code/model/tools/e5f_ssj_scaled_step_root.py` (diagnostic solver copy, one
step-rule change), `code/model/tools/collect_e5f_ssj_toeplitz_jacobian.py`
(summary), and their tests. Packets:
`output/model/e5f_sequence_space_prototype_20260913/toeplitz_jacobian_10/` and
`.../toeplitz_jacobian_10_scaled/`. No production file, kernel, gate, parameter,
population law, or residual definition changed. The afternoon batch's absolute
deadline (00:07 UTC) had effectively expired; the diagnostic contract records
that it applied the handoff's two-hour cap instead, with every pinned
scientific input unchanged.

### Why the retained root stalls

In the ten-period shocked root (job 17700926) the rebate block dominated the
score at every one of the eight mappings (37.0, 30.7, 21.9, 20.2, 12.8, 4.5,
1.35) while the housing block reached 0.008 and PAYGO 0.09. The retained
Broyden start is diagonal (\(-1.63\) for housing, \(-200\) for both fiscal
rows) and has no cross-block or cross-date entries.

### Stage 1: seven native mappings at the 2007 stationary economy (job 17714834)

Ten dates, perturbation at date 5, central step \(10^{-5}\) in the log of
each coordinate (the retained solver differentiates every coordinate in logs).
Wall times: 41 s for the baseline and 110--115 s per perturbed mapping, 722 s
in total. The baseline scaled residual is \(7.6768\times10^{-5}\) (gate
\(2\times10^{-4}\)); stationary drift is at most \(6.8\times10^{-7}\) (limit
\(10^{-5}\)); every mapping passed the unchanged mass, policy-reproduction and
feasibility gates. Measured lag profiles \(\partial R_t/\partial\log u_s\)
with lag \(t-s\):

| block | lag \(-1\) | lag 0 | lag \(+1\) | lag \(+2\) | lag \(+3\) | lag \(+4\) |
|---|---|---|---|---|---|---|
| housing \(\leftarrow\) log price | \(+1.03\) | \(-1.91\) | \(-0.15\) | \(-0.10\) | \(-0.06\) | \(-0.04\) |
| rebate \(\leftarrow\) log price | \(+205\) | \(-56.8\) | \(-31.3\) | \(-20.4\) | \(-15.0\) | \(-12.2\) |
| PAYGO \(\leftarrow\) log price | 0 | \(+0.12\) | \(+0.39\) | \(+0.88\) | \(+2.18\) | \(+4.17\) |
| PAYGO \(\leftarrow\) log pension | 0 | \(-200.0\) | 0 | 0 | 0 | 0 |
| rebate \(\leftarrow\) log rebate | \(+1.0\) | \(-197.9\) | \(+1.0\) | \(+0.8\) | \(+0.6\) | \(+0.5\) |
| rebate \(\leftarrow\) log pension | \(+4.0\) | \(+7.8\) | \(+2.8\) | \(+0.8\) | \(+0.4\) | \(+0.3\) |

Economics. With a four-year period the housing user cost is a small difference
between \(q_t\) and discounted \(q_{t+1}\), so dated demand responds to the
*slope* of the price path: a one percent higher next-period price raises
current demand by about one percent, while a one percent higher current price
lowers it by 1.9 percent. Because property-tax revenue is levied on current
owner housing, the rebate residual inherits the same anticipation structure
(\(+205\) versus \(-57\), in the \(\times200\) fiscal units). The PAYGO row
responds to *past* prices with a profile that grows with the lag (fertility
and the four-vintage entry queue), which is the one block where zero-filling
lags beyond \(+4\) is a visible approximation. The own-date two-date smoke
numbers (\(-2.00\), \(-73.5\)) differ from the ten-date ones because the
two-date horizon truncates the anticipation channel.

Predictive check. Along the reference Broyden path the measured matrix
predicts the realized residual changes with relative error 0.13--0.45 across
iterations, against 0.43--1.38 for the diagonal default; the reference's final
Broyden matrix never learned the \(+1.03\) superdiagonal or the \(+205\) entry
after seven rank-one updates (`reference_path_prediction_check.json`).

### Stage 2: identical shocked ten-period root with the measured start (job 17714834)

Same start guess, controls (slope 1.63, damping 1, componentwise log-step
clip 0.2, worsening factor 1.5), verified endpoint and eight-mapping budget as
job 17700926; only the Broyden `initial_jacobian` differs. Score (max scaled
residual) by mapping:

| mapping | 1 | 2 | 3 | 4 | 5 | 6 | 7 | final replay |
|---|---|---|---|---|---|---|---|---|
| reference (diagonal start) | 37.04 | 30.68 | 21.90 | 20.23 | 12.80 | 4.525 | 1.346 | 1.346 |
| measured Toeplitz start | 37.04 | 47.21 | 47.92 | 14.43 | 1.753 | 0.2644 | 0.02246 | 0.02246 |

Neither run meets the \(2\times10^{-4}\) gate within eight mappings; the
measured start ends 60 times lower, with housing \(1.8\times10^{-4}\) (inside
the gate), PAYGO \(2.1\times10^{-3}\), rebate \(2.2\times10^{-2}\), and an
exact final replay (reproduction gap 0). Accepted prices agree with the
reference best to three digits (0.683 to 0.391), so both roots sit in the same
basin. Root wall time 1222 s versus 1409 s (about 172 s versus 195 s per
mapping, a node difference, not an algorithmic saving); the derivative stage
added 722 s once. Terminal distances are unchanged from the reference (carried
mass 163 percent above the endpoint), which is the known ten-period truncation
issue, not a root property.

The first two iterates worsened (47.2, 47.9) because the Newton direction
asked for log-price and log-rebate increases far above 0.2 at dates 2--9 and
the componentwise clip cut the price coordinates while the rebate coordinates
moved as planned, producing rebate residuals near \(-45\). Once the clipped
coordinates stopped binding the Broyden updates recovered.

### Second arm: direction-preserving step with the same Jacobian (job 17717654)

`e5f_ssj_scaled_step_root.solve_price_path_scaled` is a byte-for-byte copy of
the retained solver except that the damped Newton step is scaled by one common
factor so its largest coordinate equals 0.2, instead of clipping each
coordinate. It reused the pinned Jacobian from job 17714834 (no new
derivative mappings). Scores by mapping: 37.04, 31.13, 18.62, 12.69, 1.567,
0.3328, 0.04671, final replay 0.04671 (exact reproduction), monotone
throughout; 959 s root wall time at about 135 s per mapping on a faster node.
The clipped-step arm ended lower (0.0225) because its two wasted iterates
happened to leave Broyden with a better-updated matrix; the scaled arm never
lost ground. Neither arm meets the \(2\times10^{-4}\) gate in eight mappings;
both were contracting the rebate block by roughly a factor of eight per
mapping at the end, so about two more mappings would be needed.

### Assessment

- The retained root's slow progress is a Jacobian-initialization problem, not
  a mapping-cost problem. Seven native mappings (12 minutes at ten dates) buy a
  start that reaches in eight mappings what the diagonal start would need well
  beyond its budget to reach.
- The measured structure is a near-differencing operator in prices
  (\(-1.91\) own, \(+1.03\) next date), inherited by the rebate row
  (\(-57\), \(+205\)). Any solver treating dates as independent will
  misjudge steps; componentwise step clipping is actively harmful with a
  coupled Jacobian.
- This is not a fake-news Jacobian and no speedup is claimed for the
  104-date production root. For that horizon the same seven-mapping
  measurement costs seven long mappings (about 3.3 hours at 28 minutes each),
  and lags beyond the measured window would be zero-filled; the growing PAYGO
  lag profile says the window must be wider than \(\pm5\). A comb
  perturbation (several dates spaced wider than the lag decay) in one central
  pair per block would measure a wider window at the same cost, but the
  stationary-Toeplitz assumption is untested along a 400-year transition.
- Recommended next decision (author): (i) allow the retained solver an
  `initial_jacobian` built once per horizon from a measured comb and reused
  across policy roots at that horizon; (ii) decide whether to adopt uniform
  step scaling in the production solver, which is a numerics change and needs
  its own verification; (iii) raise the mapping budget from eight to ten when
  a measured start is supplied. None of these was applied to production.

### September 14, 01:23 UTC: certified ten-period root by warm-started continuation

The frozen root operator refuses budgets above eight mappings, so each arm was
continued in a second eight-mapping run warm-started from its exactly
reproduced best iterate with its learned Broyden matrix (driver option
`warm_start_receipt`; jobs 17730694 clipped, 17730695 scaled). The scaled arm
converged: scores 0.04671, 0.009937, 0.0009657, \(9.587\times10^{-5}\), final
replay exact (reproduction gap 0), five mappings in 561 s. Housing
\(5.5\times10^{-7}\), PAYGO \(7.7\times10^{-6}\), rebate \(9.6\times10^{-5}\).
Total cost of the certified root: seven derivative mappings once, then
\(8+5=13\) root mappings, against eight unconverged mappings for the diagonal
start (job 17700926, best 1.346). Accepted prices 0.6825 to 0.3906, pensions
2.046 to 1.882, rebates 0.178 to 0.077. The terminal-distance diagnostic is
unchanged (carried mass 163 percent above the endpoint), so the root is
finite-horizon market/fiscal converged with the terminal approach unverified,
exactly as the reference contract labels it. Packet:
`output/model/e5f_sequence_space_prototype_20260913/toeplitz_jacobian_10_scaled_cont/`.
The clipped arm (retained solver, unmodified) also converged in its
continuation: 0.02246, 0.005601, 0.0007447, 0.0002399, \(2.690\times10^{-5}\),
final replay exact, six mappings in 946 s (job 17730694), so the certified
ten-period root costs \(8+6=14\) retained-solver mappings after the one-time
seven-mapping measurement. The two certified roots agree to
\(4\times10^{-7}\) relative in every coordinate.

### September 14, 03:45 UTC: the 104-date announced root and discrete tenure thresholds

The announced four-shock job 17711519 ended at its eight-mapping budget with
best score 0.977 (housing \(5.4\times10^{-3}\), PAYGO 0.12, rebate 0.977),
13819 s total; its second mapping tripped the worsening safeguard, which
halved damping to 0.5 and reset the Jacobian to diagonal. Two experiment-only
rescue arms (jobs 17732268 and 17732269, driver
`code/cluster/run_e5f_ssj_announced_rescue.py`) restarted from its exactly
reproduced best iterate: a pure continuation with the learned Broyden matrix
and the unmodified solver, and the extrapolated ten-date Toeplitz matrix with
the scaled step. Both first steps worsened (1.37 accepted; 1.50 tripped the
safeguard, which discards the supplied matrix).

The reason is not the Jacobian. At the warm start the rebate residual is
concentrated at isolated dates (43--44: \(-0.23,+0.36\); 68--70:
\(-0.16,+0.27,-0.98\)) plus a smooth tail rise to 0.71 at dates 100--102.
The native rows show that the isolated spikes coincide with discrete jumps in
the owner rate of 1.2--1.7 percentage points between adjacent dates (dates 25,
44, 71) on an otherwise flat ownership profile whose typical date-to-date
change is 0.16 points. A log-price step of 0.001 moved the jump from date 71
to 72 and from 44 to 45 and the residual spike moved with it. These are the
deterministic tenure-choice thresholds on the housing grid: a mass of
households flips tenure when the price path crosses a threshold, the
property-tax base moves by about 0.1 percent, and the rebate row, which is 200
times the relative imbalance, jumps by 0.2 or more. The fiscal gate
(\(2\times10^{-4}\) scaled, i.e. \(10^{-6}\) relative) is three orders of
magnitude tighter than that jump, so no smooth root exists across a straddled
threshold; the ten- and six-date roots pass because their few dates happen not
to straddle one. Along a 104-date path a handful of dates always will.

Implications (author decisions, nothing applied): (i) the announced 104-date
"unconverged" label is a discreteness floor, not a search failure; the housing
block at \(5\times10^{-3}\) and the PAYGO block are already near what the
mapping can deliver; (ii) a certifiable long root needs either a fiscal gate
commensurate with the threshold jump (of order \(10^{-3}\) relative, i.e. a
scaled gate near 0.2--0.5), or a smooth tenure margin (the previously rejected
Frechet smoothing), or a rebate rule that is not evaluated date by date at
\(10^{-6}\); (iii) the measured-Jacobian start remains the right tool for the
smooth part of the problem and cannot fix a discontinuity.

### September 14 close-out of the rescue arms

Arm A (job 17732268) was cancelled by the experiment owner after its
safeguard reset made it a duplicate of arm B (identical best score
0.6681649545882107). Arm B (job 17732269) reached mapping 6 at score 0.447
(housing \(2.6\times10^{-3}\), PAYGO \(3.4\times10^{-2}\), rebate 0.447 at
date 101) under quarter damping, then hit the diagnostic's own 16200 s
numerical deadline before its final replay; no receipt was written and the
result is a checkpoint, not a certified root. The residual is still the
threshold pair at dates 43--44 (\(-0.13,+0.43\)) and 69--70
(\(+0.07,-0.42\)) plus the smooth tail (0.32 to 0.45 at dates 99--101). The
practical floor of the current 104-date contract under a deterministic tenure
margin is therefore of order 0.4--1.0 in scaled units. No further arm was
launched.

### September 15: author-approved experiments in copies (no production change)

The author confirmed the frozen `tenure_choice_kappa = 0.005` is the June 28
lower search bound, not an estimate (issue M19), and approved testing larger
values. Two experiment families were run, all in copied code and isolated
batches; every deviation from the retained root is listed in each batch's
`submission.json` under the local packet.

**Loosened-gate 104-date root (jobs 17860152, 17865027).** Copied solver
gains a per-coordinate tolerance vector (housing rows \(2\times10^{-4}\),
fiscal rows 0.5 scaled, i.e. \(2.5\times10^{-3}\) relative) and a 12-mapping
budget; warm start from the rescue checkpoint (0.447); extrapolated Toeplitz
initial Jacobian; direction-preserving step. Job 17860152's first step was
decisive evidence: the extrapolated Jacobian drove the housing residual from
\(2.6\times10^{-3}\) at most dates to essentially zero everywhere except
dates 43--44 and 70--71 (the tenure-threshold pair, \(+1.52/-1.15\) in the
rebate row) and the tail, but the max-abs worsening safeguard then discarded
the Jacobian and reset to a half-damped diagonal. The job was cancelled and
relaunched as 17865027 with one further copied-solver change: best-point
selection and the safeguard use a trimmed score that ignores the four worst
normalized coordinates, while certification still requires every coordinate
inside its gate. Result pending at the time of writing.

**Tenure-smoothing probe (jobs 17866783, 17866789, 17866790, 17866791).**
One native mapping of the fixed warm-start path per scale, everything else
frozen; the 0.005 control reproduced the checkpoint residual to all digits.
Ownership jumps between adjacent dates (percentage points):

| scale | date 25 | date 44 | date 71 | median date-to-date change | mean ownership on the fixed path |
|---|---|---|---|---|---|
| 0.005 (frozen) | 1.25 | 1.33 | 1.67 | 0.160 | 0.732 |
| 0.02 | 0.70 | 0.81 | 0.98 | 0.182 | 0.683 |
| 0.05 | 0.45 | 0.22 | 0.50 | 0.107 | 0.678 |
| 0.1 | 0.33 | 0.09 | 0.26 | 0.058 | 0.710 |

At 0.05 the mid-path flips fall by a factor of three to four and the profile
becomes smooth except at dates 1--2 (a 2.5-point re-sorting of the 2007
initial distribution, which was solved at 0.005; an initial-condition
artifact that disappears once the initial state is re-solved at the same
scale). The cost is the level: on the fixed price path ownership is about
five points lower at every date, so a re-solved root would sit at different
prices and the residuals on the fixed path are not meaningful. Mapping times
were 1741--2471 s. At 0.1 the mid-path flips are 0.1--0.3 points and the
typical change 0.06 points, i.e. the profile is smooth, while the date-1
initial re-sorting grows to 5.5 points (same artifact). Mean ownership on the
fixed path is non-monotone in the scale (0.732, 0.683, 0.678, 0.710), so the
level effect must be read off a re-solved root, not this probe.
