# Prepared two-node conditional-branch audit

**Prepared for lead review; not submitted or numerically certified.** The driver
performs no full Bellman, stationary equilibrium, historical path, calibration,
or population-distribution solve. It reads only the exact `new_balanced`
fixed-preference checkpoint identified by the four-cell receipt. Runtime is
bounded externally at five minutes, one CPU and 16 GB; no more than 30 global
scalar branch optimizations are allowed, with 17 expected.

Files: `run_branch_audit.py`, `input_contract.json`, and
`prelaunch_validation.json`. The contract pins the driver, immutable source
snapshot, original parent contract and checkpoint. Current input-contract SHA256:
`1c63d2549fbb0113b787bcbdaa0d4e06da919f02330c9de8ad09173782899394`.
The driver verifies all 625 parent source pins and the checkpoint before load.
It requires a new output directory and never overwrites the original artifacts.

## The question and exact scope

At age 30, combined income index 8, renter tenure, location 0, parity 0 and
child-state 0, owner-entry probability falls from 0.279660255 to 0.266427934 as
wealth moves from grid index 52, 0.8604651163, to index 53, 1. These two nodes
do not cross a simple down-payment threshold. The six-room owner branch accounts
for almost the entire decrease; rental housing is six rooms at both nodes.

The scientific question is whether the saved conditional branch actions are
suboptimal under their saved continuation values, whether interpolation changes
relative branch values, or whether the decline is a valid relative-value
comparison under the discrete housing menu. The audit distinguishes these
objects; a completed diagnostic does not automatically certify the model.

The stationary checkpoint retains pension 2.0463613896121218, preference intercept
0.2900515293650047, sigma 2, the parenthood-only housing requirement
0.7168500088318734 and zero per-child slope. Shared family primitives are checked
to have zero consumption and housing floors, zero current child reward and
equivalence scale one at the inspected no-child state. All structural parameters,
prices, debt restrictions, continuation values and the shock scale remain saved
values. The owner service premium `chi` is retained.

## Critical numerical distinction: owner entry uses interpolation

The current net financial wealth of a household that purchases size \(H\) is
\(b^o=b-PH\). The entry restriction is
\(b\geq(1-\phi)PH\), equivalently \(b^o\geq-\phi PH\), with \(\phi=0.8\).
The production tenure operator interpolates the **conditional owner value
function solved on its wealth grid** at \(b^o\). It does not re-optimize the
owner's saving at this off-grid value before comparing tenures.

For six rooms, the two entry wealth states map to -2.9395064959 and
-2.7999716121. Their interpolation brackets are owner wealth nodes 24–26:
-3.0465116279, -2.9069767442 and -2.7674418605. The two- and four-room products
are also affordable, so reconstructing the complete softmax uses their three
nearby bracket nodes each. Eight and ten rooms fail entry affordability and
receive the unchanged infeasible value.

Thus the driver audits two on-grid renter states and nine unique on-grid owner
states, then six exact off-grid owner objectives as supplemental interpolation
diagnostics: **17 global scalar maximizations**. Domain grid knots and the
rental-cap kink are internal candidates within each optimization. The off-grid
optima are not substituted into the production softmax or presented as a
repaired model.

## Objective and continuation for line-by-line review

For the chosen age and income state, the saved next-period inclusive value
function is integrated over the saved income transition row, combined with the
death/bequest value, and then passed through the exact child-aging value map:

\[
\bar V_j=\mathcal A_j\left[s_j\sum_{z'}\Pi(z,z')V_{j+1}(z')
+(1-s_j)V^{bq}\right].
\]

The bequest argument is next liquid wealth plus the full value of the retained
owner house (zero housing wealth for renters). The existing pure bequest helper
preserves completed-child dependence, estate taxes, the wealth shift and any
normalization. The saved beta multiplies this combined continuation once, in
the scalar objective. The driver records the survival probability actually used.
No fertility probability is multiplied into the conditional tenure objective:
the saved next-age inclusive value already includes its future options, while
the inspected current tenure choice conditions on the post-fertility childless
state.

Resources at conditional branch wealth \(x\) are
\(R x+y_j(z)+\min\{\max[g-(R\max(x,0)+y_j(z)),0],g\}\), retaining the saved
debt-blind means test. The independent scalar evaluator uses

\[
u(c,h)=\frac{e\,[c^\alpha h^{1-\alpha}]^{1-\sigma}}{1-\sigma}
+\psi_m+\beta\bar V_j(b'),
\]

where \(c,h\) are surpluses after the saved floors. For renters, the
intratemporal optimum splits residual spending by \(\alpha\), until housing
reaches its saved maximum; after that, remaining spending goes to consumption.
For owners, effective housing is
\(\chi[H-\texttt{owner\_h\_bar\_scale}\,\bar h]\), and consumption subtracts
the saved maintenance/property-tax/size cost. The initial house purchase is
already accounted for in \(x=b-PH\); it is not subtracted again.

The lower saving bound uses the exact saved age-specific debt-taper and debt-cap
helpers. Owner secured debt is separated before applying the unsecured rollover
rule. The upper bound is the production resource bound with its 1e-6 margin.
Continuation interpolation is clipped linear interpolation at the saved wealth
grid endpoints.

Unlike the production analytical first-order candidate formula, this driver
uses **bisection of the scalar derivative on each concave segment**. It considers
all segment endpoints, every continuation-grid knot inside the feasible interval,
the renter housing-cap kink, and the unique interior derivative root when one
exists. Each segment combines concave flow utility and linear continuation.
The global maximum across these candidates remains valid when continuation
slopes are not globally monotone. This is an independent maximizer for the same
fixed-continuation objective, not an independent verification of the entire
continuation function.

## Reported evidence and interpretation

`audit.json` will retain the branch domain, saved saving, independent optimum,
saved-action value, best-minus-saved gain, reconstructed consumption/housing,
and feasible cross-evaluations of other saved actions in the same tenure. It
separately records exact off-grid versus interpolated owner values.

At each original renter node, the reconstructed conditional branch values
produce a full tenure softmax. The pairwise check is

\[
\kappa\log\frac{p_{6\mathrm{room}}}{p_{\mathrm{renter}}}
=V^{\mathrm{entry}}_{6\mathrm{room}}-V_{\mathrm{renter}}.
\]

Probabilities are saved in float32, whereas values are float64; assess residuals
at their recorded precision. The driver reports mismatches and objective gains
without silently declaring a pass or changing any acceptance gate. Material
best-minus-stored gains identify a stored-action problem under fixed continuation.
Exact on-grid optima and log-odds reproduction with a different direct off-grid
value instead isolate an interpolation issue. If both checks match and the
relative-value decrease remains, that supports the inspected two-node economic
comparison but not a universal ownership gradient or global model certificate.

## Existing source grounding and checks completed

Immutable c6dd3508 source inspected:

- `intergen_eqscale_seq_optimized/solver.py`: full Markov Bellman, lines
  2447–2735; borrowing helpers, lines 123–162; child-aging value map, lines
  7056–7094; completed fertility and bequest helper, lines 7253–7310.
- `intergen_eqscale_seq_optimized/kernels.py`: scalar renter/owner objectives,
  lines 55–130; tenure-logit operator, lines 678–764; exhaustive saving and
  conditional borrowing/allocation rules, lines 768–1150.
- `tools/run_e5f_independent_numerical_audit.py`: previously reviewed analytical
  segment oracle, lines 282–335, and fixed-continuation reconstruction, lines
  356–409. This gave the bounded audit template; the new maximizer does not
  call its production-equivalent analytical optimizer.

Local Python syntax checks and four synthetic segment-optimization checks
passed, covering rental caps, owners and both increasing and nonmonotone
continuation slopes. Independent maxima weakly dominate a 10,001-point dense
grid in each synthetic case. The first synthetic fixture omitted the `pc`
constant; setting its intended zero value fixed the fixture without changing
the driver. These checks do not load or certify the model checkpoint. Runtime
source/parameter/reproduction checks and the actual two-node result remain
pending lead review and submission.

