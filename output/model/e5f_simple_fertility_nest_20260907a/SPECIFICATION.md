# Experimental simultaneous choice with fertility nests

Status: implementation and bounded retained-parameter verification, September 7, 2026.
Production remains sequential. The author rejected the alternative construction
that adds outcome-dependent subnests to reproduce the sequential formulas.

## Economic object

A complete plan specifies whether to attempt conception and which housing
product to use after each possible biological outcome. All plan tastes are
observed together; conception remains independent uncertainty. There are six
housing products (renting and five owner sizes). With conception probability
\(0<p<1\), the menu contains six wait plans and 36 attempt plans.

Let \(q_{0h}\) and \(q_{1h}\) denote optimized continuation-inclusive values
without and with a successful conception, after selecting housing product
\(h\). The original budgets, saving choices and product sizes determine these
values. A wait plan has value \(q_{0h}\); an attempt plan has value
\((1-p)q_{0h_0}+p q_{1h_1}-pC\), where \(C\) is the successful first-birth
cost (zero for later births).

A single within-nest housing scale \(\kappa\) and outer fertility scale
\(\sigma_F\geq\kappa>0\) define a two-nest GEV distribution. Its inclusive
values are
\[
 I_W=\kappa\log\sum_h e^{q_{0h}/\kappa},\qquad
 I_A=\kappa\log\sum_{h_0,h_1}
 e^{[(1-p)q_{0h_0}+p q_{1h_1}-pC]/\kappa}.
\]
The expected value is \(\sigma_F\log(e^{I_W/\sigma_F}+e^{I_A/\sigma_F})\).
Errors use the mean-zero marginal normalization; no Euler constant is added.
Conditional probabilities are the corresponding softmax expressions. The
Cartesian attempt sum factorizes for efficient evaluation; this introduces no
additional nests and requires no numerical integration.

Conditional on attempting, failure housing probabilities are proportional to
\(e^{(1-p)q_{0h}/\kappa}\), and success housing probabilities to
\(e^{p q_{1h}/\kappa}\). This changes economic behavior: a less likely outcome
has more dispersed housing choices. There are correlated plan shocks, not two
independent additive scalar shocks. Housing can differ after success/failure.

At zero or certain conception, the unused housing coordinate is removed. Near
these endpoints, the extra contingent alternatives can still affect menu
entropy; exact endpoint collapse is discontinuous in that sense. This is an
explicit limitation of the simple complete-plan menu, not silently repaired
with additional correlations. Attempts retain the original age/parity eligibility.

## Measurement and verification

The original 12 moments, weights, 11 coordinates and bounds are retained.
Housing scale remains 0.005, supply elasticity 0.63, and the first-child room
jump upper bound 0.5. The old-state fertility normalization to 2.1 is recomputed
as a maintained derived normalization. No parameter search is authorized here.

The first-birth housing comparison retains production measurement: the same
selected origin-state birth mass is transported under successful conception
and under a childless control using the ordinary wait housing rule. Failed-
attempt housing is a different object and is not substituted for that control.

`code/model/intergen_eqscale_seq_optimized/fertility_nested.py` owns the choice
operator and population factorization. The existing joint lifecycle plumbing
passes both success and failure product probabilities. The experiment uses the
previously audited exhaustive saving solver and strict interpolation support.
A sequential control uses the same exhaustive saving solver; it retains its
existing logit interpolation and probability storage, so this control is not
an assertion that every floating-point operation is identical.

Run `python -m unittest discover -s code/model/tools -p
"test_e5f_simple_fertility*.py" -v` for independent plan enumeration,
probability/value derivatives, endpoint behavior, full-product Bellman/cohort
integration and calibration-contract checks (23 tests).

The bounded lifecycle driver is `code/model/tools/run_e5f_simple_fertility_probe.py`.
It checks ten-array exact default-off reproduction, then old sequential,
exhaustive sequential and new nested fixed-price solutions. These are not
market-clearing fits. Its receipts precede the two full historical objective
cases in `output/model/e5f_simple_fertility_nest_20260907a/`.

Each historical case is capped at two hours, with one-minute heartbeats,
unchanged numerical gates, saved dated checkpoints and complete target-fit and
parameter tables. Stop after the retained-parameter comparison. No policies,
search, automatic monitor or unsolicited figures are included.
