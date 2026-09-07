# An explicit finite neighborhood in the genuine mixed economy

September 7, 2026. Bounded read-only theory pass. Main note, model conventions,
utilities, value functions, quantitative work, and author decisions are unchanged.

## Result

The original mixed example admits a **fully quantified infinite-horizon,
two-stage transition**, with substantial renting and positive child goods costs.
The certified changes are extremely small. This closes the missing explicit
radius for one mixed economy; it does not establish economically large reforms
or broad primitive conditions.

Use exactly the Section 6 mixed example at taste scale \(\sigma_\xi=4\):
\(\chi=3/20\), \(\tau^p=467/9250\), \(\vartheta_0=141/400\), and
\(\phi_0=4/5\), with all its other original primitives and its fixed taste
location. Initially \(Y=O=P=1\), and ownership is \(11/21\).
Then, for
\[
 0<\vartheta_0-\vartheta_1\le10^{-11},\qquad
 0<\phi_1-\phi_0\le10^{-8},                              \tag{1}
\]
the unexpected preference decline has a converging baseline with lower initial
fertility and lower terminal population. At **every** baseline date
\(t_p\ge1\), the unexpected credit reform has a converging continuation
from the same actual inherited state, with higher policy-date fertility and
higher terminal population than continuing the baseline.

All original individual inequalities hold throughout both infinite paths.
Owner probabilities stay in \((0.52378068,0.52383838)\), so roughly 47.6%
rent throughout. Positive costs and material renting are fixed features of this
economy, not parameters being sent to zero. Neither taste parameter changes
following either shock.

Uniqueness concerns the certified neighborhoods and the stated strict branches.
The result implies neither monotone adjustment nor all-date fertility ordering.
No policy welfare or planner implementation claim is made.

## Conditions and exact construction

The finite certificate retains the exact six-variable system
\(Z=(P,u,Y,O,M,R)\), its two forward unknowns \(v=(u',Y')\), and original
residuals \(F(Z,v;\vartheta,\phi)=0\), update \(G\), and actual-old boundary.
The maintained branches are restrictive owner down payments, positive saving,
a slack owner physical cap, both rental caps binding, and slack old-owner
retention and estate constraints. Their **uniform verification**, rather than
stationary feasibility alone, appears below.

The checker freezes a rational invertible matrix \(Q\), supplied in full in
the script and report. Its original numerical suggestion is not a premise:
all subsequent matrix identities, inverses, and inequalities use rational
arithmetic. Let \(J\) be the existing exact derivative at the stationary
anchor, let \(\mathsf U\) be the lower-right two-by-two block of
\(Q^{-1}JQ\), and introduce proof coordinates
\[
 \zeta=Q^{-1}(Z-Z^*)=(\zeta^s,\zeta^u)\in\mathbb R^4\times\mathbb R^2.
\]
These coordinates do not replace the model's variables or value functions.
In particular \(u\) in \(Z\) remains the original housing service cost.

The original forward solve is validated before being used. With
\(C_0=F_v(Z^*,v^*)^{-1}\), the map
\(v\mapsto v-C_0F(Z,v;\vartheta,\phi)\) is a uniform contraction on an
explicit rectangle about \(v^*\). Derivative bounds and the exact zero anchor
residual prove that it maps that rectangle into itself for every state and
parameter in the specified boxes. This supplies an actual smooth function
\(g(Z;\vartheta,\phi)\), not an assumed forward adjustment rule.

The two regimes use the following bounds. Parameter radii here allow both
\(\vartheta\) and \(\phi\) to vary; the experiment (1) is a subset.

| Bound | Baseline | Policy continuation |
|---|---:|---:|
| Proof-coordinate radius \(\|\zeta\|_\infty\) | \(10^{-8}\) | \(10^{-5}\) |
| Parameter radius about original anchor | \(10^{-11}\) | \(10^{-8}\) |
| Auxiliary forward radius for each component of \(v\) | \(10^{-6}\) | \(10^{-3}\) |
| Inner derivative bound \(\|I-C_0F_v\|_\infty\) | <0.000006082 | <0.006089 |
| Inner forcing / forward radius | <0.091470 | <0.091793 |
| Infinite interior contraction bound | <0.565005 | <0.915096 |
| Interior forcing / coordinate radius | <0.051348 | <0.051684 |

The forcing plus contraction bounds are strictly below one. The relatively
large auxiliary forward rectangles are proof bounds, not actual path ranges.
Actual adjacent equilibrium states remain in their coordinate balls.

## Proof

**Infinite boundary-value contraction.** Write
\(\widetilde g(\zeta;\lambda)=Q^{-1}[g(Z^*+Q\zeta;\lambda)-Z^*]\), with
\(\lambda=(\vartheta,\phi)\). For fixed inherited claims, the initial
boundary is affine in \(Z_0\). Solving its four restrictions for
\(\zeta^s_0\) gives \(C(\mathcal I)\zeta^u_0+d(\mathcal I)\).
Define an operator on **infinite** bounded sequences by
\[
 (\mathcal T\zeta)^s_0=C(\mathcal I)\zeta^u_0+d(\mathcal I),\qquad
 (\mathcal T\zeta)^s_{t+1}=\widetilde g_s(\zeta_t;\lambda),
\]
\[
 (\mathcal T\zeta)^u_t=\zeta^u_t+
 \mathsf U^{-1}[\zeta^u_{t+1}-\widetilde g_u(\zeta_t;\lambda)].          \tag{2}
\]
Only the two-by-two block \(\mathsf U\) is inverted. The four other
coordinates propagate forward; no assumption that the full six-dimensional
forward map contracts is made.

The checker bounds
\[
 \|D\widetilde g_s\|_\infty,\qquad
 \|[0,I_2]-\mathsf U^{-1}D\widetilde g_u\|_\infty
                  +\|\mathsf U^{-1}\|_\infty,\qquad
 \|C(\mathcal I)\|_\infty.
\]
Together with the forcing bounds, these establish a self-map and a strict
contraction of the closed sequence ball. A unique infinite fixed point exists.
Equation (2) at a fixed point is exactly \(Z_{t+1}=g(Z_t;\lambda)\), with
the original boundary.

**The stationary endpoint and complete tail.** On constant sequences, use
only the two interior updates in (2). They too form a contraction and preserve
the finite-dimensional coordinate ball. Their fixed point is a stationary
state \(\zeta^*(\lambda)\); both updates together imply
\(g(Z^*(\lambda);\lambda)=Z^*(\lambda)\).

For the infinite path, put
\(e_t=\sup_{j\ge t}\|\zeta_j-\zeta^*(\lambda)\|_\infty\).
The forward rows are bounded by \(k e_{t-1}\); the backward rows by
\(k e_t\), with the appropriate certified \(k<1\). Thus
\(e_t\le\max\{k e_{t-1},k e_t\}\), implying
\[
 e_t\le k e_{t-1}\le2r k^t.
\]
This proves exponential convergence of the entire infinite path. Prices,
choices, tenure, and conditional old states inherit this rate through the
uniformly smooth original formulas. Positive stationary demography gives
\(\bar n^*=1/\nu\) and \(Y^*=O^*\).

**Actual claims and surprise timing.** For inherited owner mass \(B\),
financial claims \(A\), purchased titles \(H\), and fixed cohort masses
\(Y_0,O_0\), the exact boundary used is
\[
 Y=Y_0,\quad O=O_0,\quad R=a(O_0-B),\qquad
 M=L\left[A+P\left(H+\frac{q\tau^p\bar H B}{Y_0+O_0}\right)\right].  \tag{3}
\]
It is linear in current \(Z\), conditional on those inherited data. The old
financial claims include the original mortgage repayment. Current price and
rebates revalue the original titles through (3).

For any later baseline date, the checker bounds the preceding baseline young
cohort, its owner probability, and both conditional original saving choices.
It forms exactly
\(B=Y_{t_p-1}\pi_{t_p-1}^O\),
\(A=Y_{t_p-1}\pi_{t_p-1}^O a_{t_p}^O\), and
\(H=Y_{t_p-1}\pi_{t_p-1}^O h_{t_p-1}^O\).
Saving and mortgage repayment use the **baseline** forecasts and original
\(\phi_0\). The current cohort is bounded from the baseline next state.
No old asset or forecast is recalculated under \(\phi_1\).

The resulting common inherited-state set has
\(\|C(\mathcal I)\|_\infty<0.053193\) and
\(\|d(\mathcal I)\|_\infty/r_P<0.008389\). Its four-by-four boundary
inverse is validated by a Neumann bound with error below
\(1.138\times10^{-8}\). Hence the same policy ball works at every
\(t_p\ge1\). The baseline continuation lies in that larger ball, so at
zero policy change uniqueness identifies the policy comparison's starting
point with the actual continuing baseline.

**Finite signs.** Differentiate (2), holding inherited claims fixed. The
result is an infinite linear contraction. The coefficient bounds are evaluated
at actual adjacent-state balls, tighter than the auxiliary forward rectangles.
Forcing divided by one minus the contraction constant bounds every derivative
coordinate. The code starts every derivative coordinate in that interval and
makes 64 outward interval substitutions over 32 retained dates. The omitted
far tail **always retains the original unrestricted norm bound**; no stationary
terminal value is imposed. A separate 100-step interval iteration encloses the
constant-state derivative.

The resulting uniform rational bounds imply
\[
 \frac{\partial Y_1}{\partial\vartheta}>1.694,\qquad
 \frac{\partial N^*}{\partial\vartheta}>4.797
 \quad\text{on the baseline family},
\]
\[
 \frac{\partial Y_{t_p+1}}{\partial\phi}>0.143,\qquad
 \frac{\partial N^*}{\partial\phi}>2.22
 \quad\text{on every policy continuation}.                \tag{4}
\]
These inequalities are exact assertions in the checker, not signs read from
rounded output. Integrating (4) over the finite intervals (1) establishes the
claimed comparisons. At either surprise, its inherited \(Y\) is fixed, so
\(\partial\bar n=\partial Y_{\mathrm{next}}/(\nu Y)\).
Here \(N^*=Y^*+O^*\) counts adult households.

## Verification and original household margins

The runnable artifact is
[/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/verify_simplified_olg_mixed_finite.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/verify_simplified_olg_mixed_finite.py).
Run it with `PYTHONDONTWRITEBYTECODE=1`. It imports the unchanged repository
helper for its **exact** anchor matrices; its own interval arithmetic,
automatic differentiation, finite bounds, and infinite-sequence certificate
are explicit. The report includes frozen rational \(Q\), exact rational bound
values and derivative endpoints, as well as readable decimals.

All arithmetic assertions use rational endpoints rounded outward at 110 bits.
Square-root bounds use integer square roots; logarithms use the convergent
atanh series with a proved remainder; exponentials use a Taylor remainder
bound. The original quadratic fertility solution, renter future rent,
property-tax rebates, and both value functions enter the interval Jacobians.
The ownership formula uses exact changes in \(W^R-W^O\) from the anchor,
so the irrational taste location is fixed without floating-point rounding.
The point check requires the interval derivative to contain the existing
exact \(F_v,J,Q_\phi\), and share. Separate fraction calculations establish the exact stationary
budgets, first-order conditions, household inequalities, housing clearing,
rebates and replacement; interval containment alone is not used to infer
that a residual equals zero.

Across the full policy forward rectangle, lower margins include owner saving
0.9736609, restrictive purchase valuation gap 0.1147181, physical-cap slack
0.9999915, old retention slack 0.5421830, and estate slack 0.1812352. Renter
saving exceeds 0.8084324 and its old-cap valuation gap exceeds 0.7579979.
At every possible later surprise, reoptimized initial-owner consumption exceeds
0.7999949, retention slack exceeds 0.5430582, and estate slack exceeds 0.1830382;
old renters also remain feasible and strictly capped. All log arguments and
current and future prices remain positive. These are box inequalities, not
sampled household checks. These inequalities verify the original first-order
and complementary-slackness conditions; concavity of each conditional
household problem establishes optimality, while the original logistic rule
handles tenure choice.

Sources:

- [Original mixed equations and local theorem](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/mixed_transition_proof.md:25).
- [Exact mixed anchor and matrices](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/verify_simplified_olg_mixed_transition.py:96).
- [Original household budgets and future rent](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:111).
- [Inherited claims and surprise forecasts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:685).

## Limits and next decision

This is an explicit positive radius at a material mixed equilibrium, rather
than unspecified continuity from the all-owner limit. Its shock bounds are
numerically tiny and should not be advertised as economically substantial.
The interval estimates deliberately discard correlations and use a generous
auxiliary forward rectangle; those losses, not a demonstrated economic failure,
limit this certificate. Widening the validated region or replacing it with
shorter primitive restrictions remains open.

The main note can retain its readable local claim. This certificate belongs
in supporting proof material until wider finite ranges justify a useful
illustration. No policy-relevant magnitude, monotone path, comparison across
binding branches, welfare result, or author adoption is inferred.
