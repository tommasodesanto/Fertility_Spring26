# Simplified OLG re-audit and analytical map to the quantitative model

**Status:** review memo; no theory or paper source was changed
**Date:** September 3, 2026

## Bottom line

The simplified OLG package is now internally coherent. I do not find a remaining
algebraic error in the household problems, the capped-renter correction, the
steady-state system, the finite transition, the fertility comparative static,
or the reallocation ledger. The numerical example solves the model on the
branches stated in the text. The two important qualifications are substantive,
not cosmetic:

1. The reallocation proposition is a Pareto-improvement result for households
   alive at the reallocation date **conditional on bridge finance and
   household-specific transfers**. Without those instruments it is a
   marginal-value decomposition, not an efficiency theorem.
2. The transition proposition establishes a locally determined finite system
   for a fixed terminal date. The numerical horizon check is strong, but it is
   not an infinite-horizon existence or uniqueness theorem.

The full quantitative model also admits useful analytical results. The most
valuable ones are not closed-form policy functions. They are:

- an exact stationary person-law formula linking births per head, demographic
  renewal, fixed net migration, and population scale;
- an exact decomposition of how housing changes the fertility-attempt hazard;
- an implicit-function characterization of long-run price, transfers, and
  population;
- exact cohort and housing-demand accounting along a transition; and
- a money-metric reallocation test that respects the model's discrete housing
  products.

These objects can be stated compactly. A useful analytical appendix for the
quantitative model should be about three pages, not another long theory model.

## 1. What was checked

The review covers:

- `latex/simplified_olg_paper_theory_section.tex`;
- `latex/simplified_olg_paper_theory_appendix.tex`;
- `code/model/tools/build_simplified_olg_mixed_tenure_theory.py` and its tests;
- the active household solver in
  `code/model/intergen_eqscale_seq_optimized/`;
- the current person cohort law in
  `code/model/demographic_transition/person_cohort_law.py`;
- the terminal person-household fixed point and perfect-foresight drivers; and
- the live model and transition status in `CALIBRATION_STATUS.md`.

The old lifecycle analysis in
`latex/dynamic_intergenerational_housing_fertility_model.tex` contains several
correct ideas, but its numerical discussion and scalar renewal closure are not
the source of truth for the current quantitative exercise. The derivations
below are written against the current person law.

## 2. Re-audit of the simplified OLG model

### 2.1 Result-by-result verdict

| Object | Verdict | Essential qualification |
|---|---|---|
| Rental user cost | Correct | Rent and the property tax are end-of-period flows and enter household budgets at date-\(t\) present value. |
| Young renter budget | Correct | The old-renter cap requires the separate saving load now used in the appendix and code. |
| Young owner budget and mortgage timing | Correct | The down payment is a liquidity requirement, not an additional resource cost. |
| Old renter allocations | Correct | Both the interior and capped branches satisfy the old-age budget and Euler equation. |
| Old owner allocation | Correct | The displayed solution requires positive retention cost and slack retention and estate-composition constraints. |
| Lifetime reduction | Correct on the stated branches | The owner formula assumes the marginal purchased unit is sold while alive next period. |
| Binding-cap fertility equation | Correct | The equation is strictly decreasing over its feasible domain and has one interior root. |
| Fertility derivative | Correct | Its sign is conditional; relaxing a cap need not raise fertility if the goods-resource effect dominates. |
| Logistic tenure share | Correct | The taste draw smooths aggregate tenure without changing conditional policies. |
| Steady-state accounting | Correct | Average fertility is replacement at every positive closed steady state. |
| Steady-state existence result | Correct as a sufficient result | The contraction and boundary-sign hypotheses are assumptions that the numerical example checks. |
| Finite transition result | Correct | It is branchwise, local, and for fixed terminal date J. |
| Reallocation proposition | Correct under its implementation assumption | It is not an unconditional welfare theorem and does not compare different populations. |
| Numerical construction | Correct | All asserted branches remain active and the equilibrium residual is below 7e-10. |

### 2.2 Timing and user costs

A rental intermediary buys a unit for \(P_t\), collects \(r_t\), pays
\(\tau^pP_t\), and ends with a unit worth \(P_{t+1}\). Zero profit is

\[
q(r_t+P_{t+1})=P_t+q\tau^pP_t,
\]

or

\[
r_t=(q^{-1}+\tau^p)P_t-P_{t+1}.
\]

This is exactly the relation in the note. Since \(r_t\) is an end-of-period
flow, its current cost is \(q r_t\). At a constant price it becomes
\(r=(q^{-1}-1+\tau^p)P\).

For a young owner, the current budget contains the down payment and the present
value of the property tax. Mortgage principal reappears as next-period debt.
When the marginal unit is sold in old age, substituting the old-age resource
constraint into the current budget cancels the mortgage inflow and repayment.
The lifetime resource cost of one unit is therefore

\[
u_t^O=q(r_t+\ell_{t+1}),
\]

while for a renter it is \(u_t^R=q r_t\). The anticipated gains tax belongs in
the owner's user cost because the interior owner branch sells the marginal unit
next period. If the owner instead retained every unit until death, this branch
and its user cost would have to be replaced.

### 2.3 The capped old renter

This was the serious defect identified in the hostile audit and is now fixed.
When the old renter is unconstrained, old-age consumption, housing, and estate
spending have total log weight \(A=1+\gamma+\omega_B\). When the old rental cap
binds, housing is predetermined, so only consumption and the estate remain
interior. The correct saving load and effective lifetime resources are

\[
k^R=\beta(1+\omega_B),
\qquad
\widetilde w^R=w-q^2r_{t+1}h_R^{\max}.
\]

The young lifetime budget is consequently

\[
(1+k^R)x+\chi n+u_t^Rh=\widetilde w^R,
\]

and saving is

\[
a'=\beta(1+\omega_B)x-qT_{t+1}
   +q^2r_{t+1}h_R^{\max}.
\]

The current code and appendix use these expressions. Direct substitution into
the old-age budget and the saving Euler equation leaves residuals at numerical
precision.

### 2.4 Housing access and fertility

On a branch with binding gross-housing limit \(\widehat h\), define

\[
x=\frac{\widetilde w-u\widehat h-\chi n}{1+k},
\qquad s=\widehat h-\kappa n.
\]

The fertility condition is

\[
\frac{\vartheta}{n}=\frac{\chi}{x}+\frac{\alpha\kappa}{s}.
\]

Its substituted left side tends from \(+\infty\) to \(-\infty\) over the
feasible interval and has derivative

\[
-\frac{\vartheta}{n^2}
-\frac{\alpha\kappa^2}{s^2}
-\frac{\chi^2}{(1+k)x^2}<0.
\]

The interior root is therefore unique. Implicit differentiation gives

\[
\frac{dn}{d\widehat h}
=\frac{\alpha\kappa/s^2-\chi u/[(1+k)x^2]}
{\vartheta/n^2+\chi^2/[(1+k)x^2]+\alpha\kappa^2/s^2}.
\]

This formula is useful precisely because it does not force a positive sign.
More space reduces the space cost of children, but paying for that space also
reduces goods available to the family. The condition in the note correctly
states when the first force dominates.

### 2.5 Positive steady states

The accounting is simple. Since today's old are yesterday's young,
\(O^*=Y^*\). A constant young cohort requires

\[
\nu\bar n(P^*,T^*;\vartheta)=1.
\]

Constant prices imply \(\ell^*=0\). Housing clearing and the rebated property
tax then give

\[
Y^*=\frac{\bar H}{\bar h(P^*,T^*;\vartheta)},
\qquad
T^*=\frac{q\tau^pP^*}{2}\bar h(P^*,T^*;\vartheta).
\]

Thus a permanent change in \(\vartheta\) does not change steady-state average
fertility: housing prices adjust until fertility returns to replacement. It can
change tenure, conditional fertility, housing per household, and the scale of
the population. That replacement property is a feature of the closed model,
not a generic property of the quantitative model with fixed net migration.

The existence proposition is valid. It uses a contraction only to solve the
fiscal equation for \(T\) at each candidate \(P\), then applies the intermediate
value theorem to the replacement residual. The Jacobian in the appendix is the
correct Jacobian of those two equations.

### 2.6 Transition result

For a fixed horizon \(J\), the model has \(2(J+1)\) equilibrium conditions in
the price and transfer paths. Conditional on a strict appreciation pattern and
strict household active sets, this is a smooth finite system. A nonsingular
Jacobian gives a locally unique path on that branch for a small permanent
shock. The directional price signs then verify that the branch agrees with the
statutory positive-part gains tax.

This is a sound result. It deliberately does not prove that an infinite path
exists, that only one appreciation branch can be self-consistent, or that the
finite paths converge as \(J\to\infty\). The numerical horizon comparison is
evidence about the reported construction, not a theorem about all primitives.

One permanent shock can nevertheless generate gains for several cohorts. The
reason is the exact cohort recursion

\[
Y_t=Y_0\prod_{s=0}^{t-1}\nu\bar n_s.
\]

The first fertility response changes the next cohort; that cohort later becomes
old and produces another cohort. With a fixed housing stock, this state
propagation can move prices for more than two dates without introducing a
second exogenous shock.

### 2.7 The welfare statement

The safest verbal description is: **a Pareto-improving housing reallocation for
households alive at date \(t\), conditional on bridge finance and two-date
transfers**.

For a constrained young buyer and an interior old incumbent,

\[
MV_{it}^{Y}=q r_t+q\ell_{t+1}+\zeta_{it}^{O,F},
\qquad
MV_{jt}^{O}=q r_t-\ell_t.
\]

The private marginal-value difference is

\[
MV_{it}^{Y}-MV_{jt}^{O}
=\zeta_{it}^{O,F}+\ell_t+q\ell_{t+1}.
\]

The title-transfer ledger in the appendix is consistent. The incumbent is
compensated by the after-tax sale proceeds. The buyer pays for the unit and its
property tax, receives current services, and sells the marginal unit next
period. The government receives the two gains taxes that are absent when the
incumbent retains the unit until a stepped-up-basis estate sale. Market-rate
lenders have zero net value. After subtracting the real title-and-occupancy
cost \(K_{ijt}\), total surplus is

\[
\zeta_{it}^{O,F}+\ell_t+q\ell_{t+1}-K_{ijt}.
\]

If this is positive, transfers can divide the surplus while keeping both
affected households weakly above baseline. The result does **not** say that the
competitive allocation is inefficient relative to the original borrowing
technology. Bridge finance expands the feasible set. It also does not compare
policies that change the number or identity of future households.

### 2.8 Numerical verification

The focused OLG test suite passes all four tests. It checks the rental timing,
the capped-renter reduction and Euler equation, the owner reduction and Euler
equation, and the endpoint market, fiscal, and replacement conditions.

Additional direct checks give:

- maximum scaled transition residual: \(6.50\times10^{-10}\);
- analytical versus centered-difference fertility derivative error:
  \(3.30\times10^{-12}\);
- maximum difference between horizons 28 and 32 over the first twenty dates
  and the reported variables: \(1.89\times10^{-15}\); and
- economically material gains-tax dates: 0--3 and 7--8, as stated in the
  appendix.

The source currently reports the conservative bound \(3.11\times10^{-15}\)
for the 28-versus-32 comparison, while the saved verification file records a
24-versus-28 comparison over nine dates. The claim remains true, but the build
script should eventually store the exact comparison quoted in the prose.

## 3. Style and structure

The note is substantially clearer than the versions reviewed earlier. It now
defines the renter, owner, and old-age problems before using their solutions;
explains the rental user cost; treats goods and housing needs symmetrically;
uses ordinary aggregate notation; and puts one sentence of economics around
the important equations. It no longer reads like a record of internal model
development.

Representative papers reinforce this direction. Menzio typically states the
environment and equilibrium in ordinary sentences, gives the economic argument
around each result, uses narrow titles that describe what a proposition does,
and sends calculations that do not advance the argument to the appendix.
Fern\'andez and Fogli use the same discipline in an empirical setting: define
the object, display the specification, and immediately say what its coefficient
means. The useful lesson here is modest: do not announce a proof strategy for a
lecture-note calculation, and do not leave an equation to explain itself. See
[Menzio, *A Theory of Partially Directed Search*](https://web-facstaff.sas.upenn.edu/~gmenzio/linkies/PDS.pdf)
and [Fern\'andez and Fogli, *Fertility: The Role of Culture and Family Experience*](https://archive.nyu.edu/bitstream/2451/26107/2/5-14.pdf).

The remaining issues are choices rather than defects:

1. The main section occupies six pages in the current package, not five. It can
   lose roughly one page without losing a result by shortening the transition
   definition and moving the steady-state local-determination sentence wholly
   to the appendix.
2. “Conditionally locally efficient” is correct but heavier than the result.
   In exposition, “Pareto-improving reallocation for current households” is
   more direct; the formal definition can remain in the appendix.
3. The appendix is long because it proves branch formulas, existence, finite
   transition regularity, and a numerical example. That is acceptable for a
   standalone note. In a paper appendix, the contraction proof and some of the
   branch catalog can be compressed to short derivations.
4. Both figures are defensible in the standalone note. The transition figure
   shows the state propagation that prose alone does not establish. The wedge
   figure is simpler, but it makes the decomposition immediately visible. If
   only one figure fits in the paper, retain the transition figure and report
   the wedge decomposition in one line or a small table.
5. A single institutional sentence could make explicit that households do not
   rent out retained owner housing. The choice set already imposes this, so it
   is not a mathematical omission.

No stylistic revision should be made before the author finishes the current
read. There is no correctness emergency left that requires silently changing
the note.

## 4. What carries to the full quantitative model

### 4.1 Map of the two models

| Simplified-model object | Quantitative-model counterpart | What can be established |
|---|---|---|
| Continuous renter housing cap | Continuous renter choice with `hR_max` and child-space floor | Exact KKT shadow value on a fixed branch |
| Continuous owner housing | Discrete owner-size menu | Adjacent-product value differences or compensating variation, not a housing FOC |
| Logistic owner type | Positive logit over rent and owner products | Exact probability derivatives when alternatives remain feasible |
| Continuous fertility choice | Sequential wait/attempt logit with conception risk | Exact derivative of the attempt hazard with respect to the value gap |
| Two-cohort replacement condition | Annual person law plus headship and fixed net migration | Exact renewal ratio and population-scale formula |
| Two-period distribution law | Full wealth-age-income-tenure-parity distribution | Exact forward accounting; no closed-form invariant distribution |
| Gains tax and stepped-up basis | No corresponding object in the active code | Cannot claim a quantitative lock-in tax wedge |
| Real reallocation cost | Six-percent seller transaction wedge and other moving costs | Can enter a quantitative money-metric reallocation test |
| Finite transition IFT | Price-transfer path residual system | Branchwise local derivative if the path Jacobian is nonsingular |

The analytical bridge therefore has to preserve the mechanism while changing
the mathematical object. In particular, it would be wrong to paste the
simplified owner's marginal housing FOC into a model in which owner housing is
a discrete product.

### 4.2 Exact stationary person law

This is the cleanest new result for the quantitative model.

Let \(\mathbf p\) collect persons by sex and single year of age. Let \(A\) be
the annual survival-and-aging operator, \(\mathbf f\) the surviving-newborn
vector generated by one birth, \(\mathbf d\) the headship-rate vector, and
\(\mathbf m\) fixed annual net migration. Once the terminal household fixed
point has been solved, let \(b\) be annual births per household head. The
person law is

\[
\mathbf p'=A\mathbf p+\mathbf f\,b\,\mathbf d'\mathbf p+\mathbf m.
\]

Define the lifetime person response to one annual birth and the response to
fixed migration by

\[
\mathbf v=(I-A)^{-1}\mathbf f,
\qquad
\mathbf w=(I-A)^{-1}\mathbf m.
\]

Two scalar objects summarize the stationary demographic feedback:

\[
L=\mathbf d'\mathbf v,
\qquad
M_H=\mathbf d'\mathbf w.
\]

Here \(L\) is lifetime household heads generated by one annual birth under the
fixed survival and headship schedules, while \(M_H\) is the stock of heads
supported by fixed net migration absent endogenous births. The renewal ratio is

\[
R=bL.
\]

At a stationary point, write annual births as \(B=bS\), where
\(S=\mathbf d'\mathbf p\) is household-head mass. Then

\[
\mathbf p=\mathbf w+B\mathbf v,
\qquad
B=b(M_H+BL).
\]

Therefore, when \(R<1\) and \(M_H>0\),

\[
\boxed{
B^*=\frac{bM_H}{1-R},
\qquad
S^*=\frac{M_H}{1-R},
\qquad
\mathbf p^*=\mathbf w+\frac{bM_H}{1-R}\mathbf v.
}
\]

These are the formulas implemented in the current person-law solver. Net
migration can be signed by cell; existence additionally requires the implied
stationary person vector to be nonnegative.

The saved August 26 terminal roots provide a direct arithmetic check, although
they are not a recalibration of the September 3 quantitative profile. The
one-percent-tax root has \(b=0.016993\), \(R=0.556841\), and
\(S=0.666119\); the two-percent-tax root has \(b=0.017064\),
\(R=0.559181\), and \(S=0.669656\). Both imply
\(L=32.768867\) and \(M_H=0.295197\), and substituting those two objects into
\(S=M_H/(1-R)\) reproduces the saved head masses.

In a closed economy, \(\mathbf m=0\). A nonzero stationary population then
requires

\[
R(P,T,g)=1.
\]

The housing market determines its scale. If \(R<1\), the only closed stationary
population is zero; if \(R>1\), a fixed-policy population cannot be stationary.
This is the full-model analogue of \(\nu\bar n=1\), but it accounts for
fertility timing, survival, sex composition, and headship.

The renewal ratio should not be confused with the dominant eigenvalue reported
in the terminal spectral diagnostic. \(R\) governs lifetime replacement and
stationary scale. The approximately 0.96 four-year eigenvalue measures the
speed at which the frozen person block forgets a perturbation. Long-lived
cohorts can make that convergence slow even when lifetime renewal is well below
one.

For the current architecture the scalar \(R=bL\) is enough: births per head is
a scalar and newborn sex composition is fixed. A Perron-root statement becomes
necessary only if fertility and descendant entry types form a genuinely
multi-type next-generation system. Adding that machinery now would make a
simple result look harder than it is.

### 4.3 Long-run equilibrium comparative statics

Let \(z=(\log P,T)\) collect the terminal house price and equal transfer, and
let \(g\) be a policy or preference shifter. Solving the household and
distribution fixed point at \((z,g)\) produces:

- the renewal ratio \(R(z,g)\);
- average occupied housing per head \(\bar h(z,g)\); and
- housing supply \(K(P,g)\).

For a pure rebated property tax with no other outlay, let \(\tau^p\) denote the
model-period tax rate used in the household problem, not the annual statutory
rate. Fiscal balance and housing clearing imply

\[
T=\tau^p P\bar h.
\]

Suppose the policy leaves the survival, headship, and fixed-migration schedules
unchanged, so that \(M_H>0\) is fixed. In the fixed-migration economy, the two
terminal equations can be written as

\[
\begin{aligned}
F_H(z,g)
&=\log M_H-\log[1-R(z,g)]
  +\log\bar h(z,g)-\log K(P,g)=0,\\
F_G(z,g)
&=T-\tau^p(g)P\bar h(z,g)=0.
\end{aligned}
\]

At a regular root,

\[
\boxed{
\frac{dz^*}{dg}
=-\left(\frac{\partial F}{\partial z}\right)^{-1}
  \frac{\partial F}{\partial g}.
}
\]

This is an analytical result even though its derivatives are evaluated
numerically. It decomposes a policy into its effects on fertility/renewal,
housing per head, supply, and the fiscal transfer. The existing terminal solver
already constructs a finite-difference Jacobian in price and transfer, so the
required numerical object is close to what the code computes now.

A transparent special case holds the transfer fixed. Define

\[
\eta_R=-\frac{\partial\log R}{\partial\log P},\quad
\xi_g=\left.\frac{\partial\log R}{\partial g}\right|_P,
\quad
\eta_h=-\frac{\partial\log\bar h}{\partial\log P},\quad
d_g=\left.\frac{\partial\log\bar h}{\partial g}\right|_P,
\]

and let \(\eta_K\) and \(s_g\) be the corresponding supply elasticity and
direct supply shift. With

\[
\lambda=\frac{R}{1-R},
\]

fixed migration implies

\[
\boxed{
\frac{d\log P^*}{dg}
=\frac{\lambda\xi_g+d_g-s_g}
{\eta_K+\eta_h+\lambda\eta_R},
\qquad
\frac{d\log S^*}{dg}
=\lambda\left(
\xi_g-\eta_R\frac{d\log P^*}{dg}
\right).
}
\]

The multiplier \(\lambda\) is the demographic amplification from fixed
migration. As renewal approaches one from below, a small change in births per
head has a large effect on population scale. Housing supply and household
demand absorb part of that effect through the equilibrium price.

For a closed economy, the root is \(R(P^*,T^*,g)=1\). If the transfer is fixed,

\[
\frac{d\log P^*}{dg}=\frac{\xi_g}{\eta_R},
\qquad
\frac{d\log S^*}{dg}
=s_g-d_g+(\eta_K+\eta_h)\frac{\xi_g}{\eta_R}.
\]

With an endogenous rebate, the two-by-two implicit system should be used
instead of suppressing the transfer channel. This matters in the current tax
experiment: the one-date Shapley calculation shows that the rebate can more
than offset the direct fertility effect of the tax.

### 4.4 Household mechanisms that remain analytical

The dynamic policies themselves do not have useful closed forms. Income risk,
saving, borrowing constraints, survival, bequests, and future values make a
closed-form Bellman solution unrealistic. Several local identities are exact,
however.

**Rental access.** The renter block is more tractable than the full Bellman
problem. Conditional on next-period assets \(a'\), let \(X\) be resources left
after subsistence consumption \(\bar c\), subsistence housing \(\bar h\), and
saving have been paid for. On a branch where the consumption floor and rental
cap are slack, the within-period solution is

\[
c-\bar c=\alpha X,
\qquad
h-\bar h=\frac{1-\alpha}{r_t}X.
\]

If this housing choice exceeds \(h_R^{\max}\), total rental housing is instead
set to \(h_R^{\max}\), and consumption absorbs the remaining resources. These
are the formulas used by the active solver. The saving choice and continuation
value remain numerical. Equivalently, conditional on a fixed active set,

\[
\frac{u_h}{u_c}=r_t+\zeta_t^R,
\qquad \zeta_t^R\geq0,
\]

where \(\zeta_t^R\) is the consumption-unit multiplier on the upper rental cap.
The child-space floor shifts \(\bar h\) and therefore reduces the discretionary
space available below that cap; it is not a second upper-cap multiplier. This
is the direct full-model counterpart of the simplified access wedge.

**Discrete housing products.** Owners choose a rung rather than a marginal
quantity. Let \(V_j\) be the value of housing product \(j\), including its
transaction, down-payment, and continuation consequences. With positive
tenure-product taste scale \(\kappa_T\),

\[
\pi_j=\frac{e^{V_j/\kappa_T}}
{\sum_k e^{V_k/\kappa_T}},
\]

and

\[
d\pi_j
=\frac{\pi_j}{\kappa_T}
\left(dV_j-\sum_k\pi_kdV_k\right).
\]

This maps a policy-induced change in alternative values into ownership and
owner-size probabilities. It is preferable to inventing a continuous owner
housing FOC that the model does not have. Feasibility changes and binding debt
constraints still create branch changes, so the formula is local.

**Sequential fertility.** At age \(a\), write the value difference between
attempting and waiting as \(\Delta V_a=V_a^A-V_a^W\). The attempt probability
and birth hazard are

\[
\pi_a^F=\Lambda\!\left(\frac{\Delta V_a}{\kappa(n)}\right),
\qquad
b_a=f_a\pi_a^F.
\]

For a smooth shifter \(g\) that leaves conception risk \(f_a\) and the logit
scale \(\kappa(n)\) fixed,

\[
\boxed{
\frac{\partial b_a}{\partial g}
=\frac{f_a\pi_a^F(1-\pi_a^F)}{\kappa(n)}
  \frac{\partial\Delta V_a}{\partial g}.
}
\]

This is the quantitative version of “housing access changes the price of a
child.” Housing affects fertility only to the extent that it changes the value
of the child state relative to waiting. The formula also explains why a policy
can move ownership but barely move births: the ownership response can occur in
states with small mass at risk for a birth, or it can change tenure without
materially changing child-relevant housing services or \(\Delta V_a\).

The active sequential-fertility code permits a sharper decomposition. Let
\(V_a^B\) be the continuation value if conception succeeds and let \(C_a\) be
the fixed utility cost of the birth, which is nonzero only where the code
specifies one. The attempt alternative is

\[
V_a^A=(1-f_a)V_a^W+f_a(V_a^B-C_a),
\qquad
\Delta V_a=f_a(V_a^B-C_a-V_a^W).
\]

Thus, if \(f_a\), \(C_a\), and \(\kappa(n)\) do not change with \(g\),

\[
\boxed{
\frac{\partial b_a}{\partial g}
=\frac{f_a^2\pi_a^F(1-\pi_a^F)}{\kappa(n)}
 \frac{\partial(V_a^B-V_a^W)}{\partial g}.
}
\]

One factor of \(f_a\) makes a successful birth conditional on an attempt; the
other makes the attempt value less responsive when conception is unlikely.

Combining the identities gives a useful analytical chain:

\[
g
\longrightarrow
\{V_j,\pi_j,h\}
\longrightarrow
\Delta V_a
\longrightarrow
b_a
\longrightarrow
R
\longrightarrow
\{P^*,S^*\}.
\]

Every arrow can be reported as a derivative or decomposition from the solved
model. What the calibration does not currently supply is a causal empirical
target for the housing-cost-to-fertility part of this chain.

### 4.5 Transition analysis

The full transition can be written as a finite residual system without
pretending that the Bellman problem has a closed form. Let \(z_T\) stack the
price and transfer paths. Backward induction maps \(z_T\) into values and
policies; the household and person laws map those policies and the initial
state into distributions; market clearing and fiscal balance produce residuals
\(F_T(z_T,g)\).

The forward part has an exact local accounting equation. If \(\mu_t\) is the
household distribution and \(Q_t\) is the transition operator induced by the
current policies, then

\[
\mu_{t+1}=Q_t\mu_t+e_t,
\qquad
d\mu_{t+1}=Q_t\,d\mu_t+(dQ_t)\mu_t+de_t,
\]

where \(e_t\) contains entrants and any externally imposed population flows.
The second equation separates inherited demographic composition from the
current policy response. It is also the exact reason a one-date preference or
policy change can affect markets for many later dates.

On a branch with stable feasibility and constraint status, a nonsingular path
Jacobian implies

\[
\frac{dz_T}{dg}
=-F_{T,z}^{-1}F_{T,g}.
\]

This is the natural full-model analogue of the finite-transition proposition.
It would yield a local impulse response and show which dates are responsible
for weak convergence. Positive tenure smoothing helps, but debt limits, renter
caps, and interpolation boundaries can still make the derivative branchwise.

The current quantitative evidence does **not** support describing the active
exercise as a completed transition between two verified closed steady states:

- the calibrated 2023 state is a dated transition state, not a steady state;
- the current person-law terminal roots use fixed net migration and have
  renewal below one;
- the H128 reform path passes its market and fiscal gates, but the matched
  baseline path does not; and
- the approximately 0.96 demographic eigenvalue is a frozen-person-law
  persistence measure, not the stable root count of the full equilibrium
  system.

Stationary derivatives at the certified terminal roots remain meaningful even
while the baseline transition algorithm is unresolved. They should be labeled
as terminal comparative statics, not as estimates from a certified transition.

### 4.6 Reallocation and welfare in the quantitative model

The simple reallocation theorem cannot be copied literally for three reasons.
Owner housing is discrete, the model has a seller transaction wedge, and the
active code has no capital-gains tax or stepped-up basis. There are therefore no
quantitative counterparts to \(\ell_t\) or \(q\ell_{t+1}\).

The right full-model object is compensating variation. At fixed prices and a
fixed continuation environment, let \(WTP_Y\) be the largest wealth reduction
that leaves a young household indifferent when access to a larger product is
added. Let \(WTA_O\) be the wealth compensation that leaves an old owner
indifferent when the owner releases or downsizes the corresponding housing.
Both are defined by inverting optimized value functions, so they remain
meaningful with discrete choices and nonlinear utility.

There are two equivalent accounting conventions, but they must not be mixed.
If the value comparisons exclude the real resource cost \(K\) of implementing
the transfer, the gross money-metric surplus is

\[
\mathcal S=WTP_Y-WTA_O-K.
\]

If the actual model transitions used in the value comparisons already include
the seller transaction and moving costs, the net surplus is instead
\(\mathcal S=WTP_Y-WTA_O\). Subtracting \(K\) again would double count it.

If \(\mathcal S>0\), household-specific transfers can make the current
households weakly better off and one strictly better off. With a continuum of
households, a small mass can be moved between adjacent products even though an
individual product is indivisible. This is a local current-household result,
not a global welfare theorem.

The computation would be useful and feasible:

1. recover alternative-specific values for renter and owner products;
2. convert value differences into exact compensating variations by wealth
   inversion rather than dividing by a common marginal utility;
3. compute expected compensating variation within each observed state using
   the model's logit distribution; realized household-specific taste shocks
   are integrated out and should not be presented as observed heterogeneity;
4. rank old states by the compensation required to release a room and young
   states by willingness to pay for it;
5. match released and desired rooms and report the surplus schedule, counting
   seller and moving costs exactly once; and
6. repeat with fertility held fixed and with households allowed to reoptimize
   fertility, clearly labeling the two exercises.

For policy welfare, the least controversial first object is consumption-
equivalent variation for households alive in 2023. It includes their own
future consumption, housing, fertility choices, and bequest motive, but assigns
no utility weight to households not yet born. If the paper later wants a social
welfare comparison across paths with different populations, it must choose a
population criterion explicitly; the model does not supply one.

## 5. Recommended analytical package

### Minimum viable package

1. **One person-law proposition.** State and prove
   \(S^*=M_H/(1-R)\) for fixed migration and \(R=1\) for a nonzero closed
   stationary population.
2. **One housing-to-fertility identity.** Report
   \(\partial b_a/\partial g\) as the logit slope times the change in the
   attempt-versus-wait value gap.
3. **One equilibrium comparative-static system.** Write the two equations in
   \((P,T)\) and the implicit derivative. Report the renewal, housing-demand,
   supply, and fiscal components numerically when those derivatives have been
   computed.

This is enough to show that the simplified theory and quantitative model share
an economic mechanism without pretending they share every institution.

### Valuable technical extension

Construct the compensating-variation reallocation schedule across old and
young states. This would answer whether the calibrated model contains an
economically important intergenerational allocation gap, not merely whether a
borrowing constraint binds somewhere.

### Claims not currently available

- a quantitative capital-gains lock-in effect;
- a certified closed-to-closed quantitative transition;
- an infinite-horizon transition existence or uniqueness theorem;
- a global welfare theorem with endogenous population; or
- an empirically identified housing-cost-to-fertility elasticity.

These are limitations of the present model or evidence, not failures of the
simplified theory.

## 6. Suggested division of labor between the two theory pieces

The simplified OLG should do only two things:

1. show how a housing-access wedge and realization-based taxation can generate
   a Pareto-improving reallocation among current households; and
2. show why one permanent fertility change can move the economy between
   positive steady states while generating price gains only during the
   transition.

The quantitative analysis should then do three different things:

1. measure how housing-product values change the fertility-attempt hazard;
2. translate births per head into renewal and population scale using the exact
   person law; and
3. decompose equilibrium policy effects into prices, transfers, housing use,
   tenure, fertility, and demographic propagation.

Keeping that division explicit solves the apparent reconciliation problem. The
simple model supplies the economic logic and the tax-lock-in illustration. The
quantitative model measures the housing-fertility and demographic mechanisms
it actually contains. It should not be asked to validate a capital-gains tax
that has not been implemented.

## 7. Verification record

- `test_build_simplified_olg_mixed_tenure_theory.py`: 4 tests passed.
- `demographic_transition.tests.test_person_cohort_law`: 6 tests passed.
- `test_run_e5f_perfect_foresight_person_demography.py`: 2 tests passed.
- `test_build_e5f_person_demography_terminal_spectrum.py`: 2 tests passed.
- Independent horizons 28 versus 32 recomputed without changing saved output.
- The protected `latex/JMP_DS_draft/` folder and both simplified-theory source
  files were not modified.
