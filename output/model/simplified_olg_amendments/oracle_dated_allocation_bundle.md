Please treat the attached follow-up as the active task, superseding the previous stationary-cohort welfare question. The individual housing comparison is now primary; derive the aggregate result separately. Use the exact household-model attachment already provided in this chat.

### File: docs/prompts/oracle_simplified_olg_dated_allocation.md
```md
# Follow-up: one dated housing allocation, not stationary-cohort welfare

Please refocus the analysis. Use the household model in my previous attachment
in this chat. Your previous answer compares stationary cohort welfare and
obtains a general gain from permanent old-to-young transfers. That is not the
question I want answered. I want a simple analytical statement that a planner
would assign MORE HOUSING TO THE CURRENT YOUNG than they receive in equilibrium.
The author now prefers an INDIVIDUAL result as the main statement, with an
aggregate statement alongside it when justified.

## Individual first, aggregate separately

The primary target compares the SAME young household in the market and full
dated planner allocations:
\[
h_i^{y,SP}>h_i^{y,eq}\qquad\text{for }i\in\mathcal I.
\]
Identify the economically relevant set \(\mathcal I\), preferably through
primitive wealth/income restrictions that also establish binding young
financing, rather than by assuming the desired marginal-utility inequality.
Try for a broad result but do not require every young household to receive
strictly more housing: a household already at its retained physical cap
cannot do so, and redistribution across heterogeneous young households can
give some less. State the set's positive mass from the type distribution.

Then separately seek
\[
H_Y^{SP}>H_Y^{eq}.
\]
More total young housing does not imply more housing for every young household.
Conversely, more housing for a subset does not establish the aggregate sign.
Either derive the additional aggregation conditions or identify that part as
unproved. If all young weakly gain housing and a positive mass strictly gains,
aggregation follows; otherwise account for the households whose housing falls.

Do not confuse the primary comparison with young-versus-old housing INSIDE
the planner allocation. Both claims above concern the full planner optimum,
not just a feasible local transfer or its welfare derivative. A local
individual result can be reported as a fallback, but must be labeled as such.

## The precise comparison

Take one date of a positive stationary competitive equilibrium of the existing
model as the reference cross-section. Stationarity supplies the current state;
the welfare comparison itself is entirely one-date. Fertility, population,
current tenure and future real allocations/estate payments are fixed. The
planner chooses BOTH current consumption and current housing for all young
and old households. It can relax individual financial restrictions and make
balanced transfers. Total current goods and housing are conserved. All actual
tenure-specific housing-size caps remain in the planner problem and may bind.

Give each currently living household weight one on its remaining utility.
Because future allocations and estates are preserved, the changing part is
exactly
\[
\int_Y[\log(c_i-\chi n_i)+\alpha\log(h_i-\kappa n_i)]\,di
+\int_O[\log c_j^2+\gamma\log h_j^2]\,dj.
\]
Fertility and ownership-taste terms are also constant. Do not replace this
objective by stationary cohort utility or a discounted sum of future cohorts.
This is a redistributive utilitarian comparison, not a Pareto-improvement claim.
The full choice set HERE means all current consumption/housing allocations
conditional on the specified fixed objects, not unrestricted dynamic choices.

Equal cohort masses and tenure shares are consequences of this stationary
reference, not extra assumptions. Tenure is chosen young and retained when old;
the stationary endowment law and choice rule give both ages the same joint
type/tenure distribution. This does NOT imply equal housing quantities.
Outside stationarity the shares need not coincide; that extension is deferred.

## Financial feasibility already checked; verify it without changing the task

For owner housing changes set
\[
\Delta a'_i=-P_{t+1}\Delta h_i^y,\qquad
\Delta a_j^e=-qP_{t+1}\Delta h_j^o.
\]
These preserve young continuation resources and old estate payments. Renters'
corresponding financial positions remain fixed. With
\(u_t=(1+q\tau^p)P_t-qP_{t+1}\), current transfers are
\(\Delta c_k+u_t\Delta h_k\). They balance because total current consumption
and housing are unchanged; rental-intermediary financing adjusts consistently.
Individual balance sheets may change while future net real opportunities do
not. Old nonnegative financial saving is relaxed by the direct planner, as is
young borrowing. Do not turn it into a physical housing restriction.

We do not need to choose outside estate recipients versus internal inheritance
for this variation if exactly the same estate payments are preserved. Check
this assertion. If a genuine additional accounting restriction is necessary,
identify precisely the affected obligation; do not reopen an unrestricted
stationary estate-optimization problem.

## Verified static progress

Let Q be the common normalized type distribution, N each age-group mass, and
\(\bar h_i\) the type's tenure-specific cap. For a common housing resource
multiplier \(\lambda>0\), the full static optimum has
\[
h_i^{y*}=\min\{\bar h_i,\kappa n_i+\alpha/\lambda\},\qquad
h_i^{o*}=\min\{\bar h_i,\gamma/\lambda\}.
\]
Consumption is optimized jointly; its separability with fixed fertility means
that doing so does not alter these housing choices. If \(\alpha\ge\gamma\),
\(n_i>0\), and the stock does not exhaust everyone's cap,
\[
\bar H<2N\int\bar h_i\,dQ,
\]
then \(H_Y^{SP}>\bar H/2\). The proof is typewise: young housing is at least
old housing, strictly wherever the old household is uncapped. If every old
were capped, every young would be capped too, contradicting the stock bound.
Thus \(H_Y^{eq}\le H_O^{eq}\) is a sufficient market-allocation condition for
\(H_Y^{SP}>H_Y^{eq}\). This allows individual caps to bind.

This is a useful intermediate lemma, NOT the final equilibrium theorem.
The missing task is to derive the individual comparison and, separately, the
aggregate ordering from economically intelligible primitives.
Do not merely ASSUME that young marginal utility is higher or that equilibrium
lies below the planner target. Also do not infer the full aggregate housing
direction from a positive local welfare derivative: with heterogeneity and
caps that implication can fail.

## A conservative primitive route to improve

The existing model uses current cash \(w=y^y+b+T\), old income \(v=y^o+T\),
\(d_p=1-q+q\tau^p\), \(d_L=1-\phi+q\tau^p\),
\(a_p=1+q\tau^p\), \(K=1+\gamma+\omega_B\), and
\(J=\gamma+\omega_B\). Both old-estate regimes, including physical caps, give
\[
h_o^d=\min\{h_d^{\max},g_d z/P\},\qquad
g_R=\frac{\gamma}{Kd_p},\qquad
g_O=\frac1K\min\{\gamma/d_p,J/a_p\}.
\]
The original young cash/mortgage constraints imply
\[
\frac vw\ge A_R=K/\gamma \quad\Rightarrow\quad h_o^R\ge h_y^R,
\qquad
\frac vw\ge A_O=\frac{g_O^{-1}-1+\phi/q}{d_L}
\quad\Rightarrow\quad h_o^O\ge h_y^O.
\]
For owners use \(Ph_y<w/d_L\) and
\(z\ge v+(1-\phi/q)Ph_y\); for renters use \(Ph_y<w/d_p\) and \(z\ge v\).
Clipping old housing at the common cap preserves the ordering.

Let \(A=\max\{A_R,A_O\}\). At zero tax the all-type inequality
\(y_i^o/(y_i^y+b_i)\ge A\) contains only primitives. This is an explicitly
identified subcase, not authorization to drop property taxes from the model.
For \(0\le\tau^p<2\), the existing bound is
\[
\bar T=\frac{\tau^p\int(y^y+b+qy^o)\,dF}{(1-q)(2-\tau^p)}.
\]
Then \(y_i^o\ge A(y_i^y+b_i)+(A-1)\bar T\) is a fully primitive sufficient
replacement. These bounds can require excessively high old income. They are
starting points to sharpen, not economically mild assumptions already accepted
by the author. They also do not independently establish that young finance
binds or isolate its contribution from ordinary redistribution.

## What I need from this run

1. Verify the dated planner and feasibility, preserving the question above.
2. Make the INDIVIDUAL housing comparison primary. Determine its scope and give
   explicit, interpretable parameter/endowment-distribution inequalities under
   which the specified young households receive larger homes in the full
   planner allocation. Then derive the aggregate result separately where
   possible. State whether the conclusions hold in every reference competitive
   equilibrium satisfying the restrictions. Conditional existence must be
   distinguished from an equilibrium-existence theorem.
3. Try to improve substantially on the conservative income bounds: exploit
   household optimality, permit aggregate rather than all-type restrictions,
   and retain both old-estate regimes and meaningful binding caps. The
   condition \(\alpha\ge\gamma\) is a candidate preference restriction, not
   a model implication. Keep genuine income/wealth heterogeneity.
4. Explain what the financial and tenure restrictions contribute. An equal-
   weight planner can redistribute even without a market imperfection; do not
   label all such redistribution a borrowing distortion. A stronger claim
   specifically attributable to a binding young constraint requires proof.

Deliver ONE main proposition, with an individual part and a separately
qualified aggregate part if proved, its proof, and a short interpretation of the
assumptions. Show the exact condition before discussing how mild it might be.
State clearly what the model implies versus what the theorem requires. If a
stronger conjecture fails, give the precise obstruction and the best valid
static result. Do not substitute a numerical point plus an open neighborhood
for analytical restrictions. Do not change utility, mortgages, inheritance or
tenure timing without explicitly separating the proposed change. Do not solve
fertility policy, transitions, money creation, or stationary-cohort welfare
in this run. Aim for a compact problem-set-style argument, with at most a short
technical appendix, rather than another menu of unrelated theorems.
```
