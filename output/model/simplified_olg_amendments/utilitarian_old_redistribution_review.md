# Balanced redistribution between existing old owners

Analytical verification, September 8, 2026. The result concerns the original
household constraints, timing, and equilibrium accounting in
[the conventional-finance discussion draft](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_conventional_finance.tex).
Only this review was written; no model run or manuscript edit was made.

**Result.** A one-time, balanced cash transfer from richer to poorer existing
old owners can strictly increase utilitarian welfare while preserving every
private borrowing restriction, household optimization, the housing market,
and the entire equilibrium price path. The sufficient regime is that the
selected old owners have slack housing caps and estate floors and different
pre-transfer resources. This is general wealth redistribution between old
households. It proves neither a housing transfer toward young households nor
a distortion attributable specifically to housing finance.

## 1. Policy instrument and theorem

The policy adds household-specific, signed lump-sum transfers to **existing
old owners only** at the intervention date. It is unanticipated before that
date, announced before those old households make their current choices, and
explicitly expires with the current old cohort. Transfer eligibility and
amounts depend on predetermined identities or pre-transfer states, not on
subsequent housing, consumption, or saving choices. The fund is balanced at
the intervention date. The original property-tax rate and equal property-tax
rebate remain unchanged.

This defines the relevant constrained policy class: household budgets,
mortgage and nonnegative-saving restrictions, and private reoptimization are
all retained, while the planner can make balanced individual cash transfers.
If policy instruments were restricted to a uniform property-tax change and
its uniform rebate, this additional redistribution instrument would not be
available; the theorem would not establish implementability in that narrower
class.

**Proposition.** Start from a positive stationary equilibrium. Suppose two
positive-mass sets of existing old owners, \(A\) and \(B\), have slack old
housing caps and estate floors, with their pre-transfer resources uniformly
ordered as
\[
0<\mathop{\rm ess\,sup}_{i\in A}z_i
 <\mathop{\rm ess\,inf}_{j\in B}z_j<\infty.
\tag{1}
\]
Suppose the selected sets can be restricted to positive mass with a uniform
positive housing-cap margin and resources bounded away from zero. Under equal
utilitarian weights within the old cohort, a sufficiently small balanced
transfer from \(B\) to \(A\) strictly raises welfare. Every household then
optimizes subject to its original constraints. Aggregate housing, current
goods use, old financial saving, and estate payments, all young choices, and
all future equilibrium choices and cohort sizes can remain unchanged. The
original price and rebate path continues to clear every market exactly.

No restriction comparing \(\beta\) and \(q\) is needed. Giving each old
household the common lifetime factor \(\beta>0\) instead of weight one simply
multiplies the welfare gain. Arbitrary person-specific welfare weights would
instead require comparing their weighted marginal utilities.

## 2. Exact policy and accounting proof

At the intervention date write
\[
z_i=a_i+P_tH_i+y_i^o+T_t,
\quad p_t=(1+q\tau^p)P_t-qP_{t+1}>0,
\quad K=1+\gamma+\omega_B.
\]
Here \(a_i\) and \(H_i\) are the old household's predetermined net financial
position and housing title. With a cash transfer \(t_i\), its old problem
uses resources \(z_i+t_i\). As long as its housing cap and estate floor are
slack, the unique optimal choices are
\[
c_i^2=\frac{z_i+t_i}{K},\qquad
h_i^2=\frac{\gamma(z_i+t_i)}{Kp_t},\qquad
e_i=\frac{\omega_B(z_i+t_i)}{Kq}.
\tag{2}
\]
Consequently its optimized utility is \(K\log(z_i+t_i)\) plus a constant
common across these old owners at the fixed price path.

Let \(m_A\) and \(m_B\) be the actual household masses of the two sets. For
a total transfer \(\epsilon>0\), give every household in \(A\) the amount
\(\epsilon/m_A\), tax every household in \(B\) the amount
\(\epsilon/m_B\), and leave all other households untouched. Thus
\[
\int t_i\,dM^O(i)=0,
\]
even when the two sets have different masses. Equation (2) implies the exact
aggregate identities
\[
\int\Delta c_i^2\,dM^O=0,\qquad
\int\Delta h_i^2\,dM^O=0,\qquad
\int\Delta e_i\,dM^O=0.
\tag{3}
\]
For example, total housing of the poorer set rises by
\(\gamma\epsilon/(Kp_t)\), and housing of the richer set falls by that
same amount. The transfer reallocates housing **within the old cohort**.

The financial and title accounting also holds individually. The old owner's
financial saving during old age is
\[
a_i^e=q(e_i-P_{t+1}h_i^2),\qquad
\Delta a_i^e=
\frac{1}{K}\left(\omega_B-rac{qP_{t+1}\gamma}{p_t}\right)t_i.
\tag{4}
\]
The coefficient is positive when the estate floor is slack. Positive
post-transfer resources therefore keep \(a_i^e>0\), even for taxed donors.
Summing (4) gives zero change in aggregate financial saving. With
\(A_t=(1+q\tau^p)P_t\), the original cash budget is verified directly:
\[
\Delta c_i^2+\Delta a_i^e+A_t\Delta h_i^2=t_i.
\tag{5}
\]
Predetermined financial claims \(a_i\) and titles \(H_i\) are not forgiven
or altered. Each old owner resizes through the model's original purchase and
sale market; increased purchases by recipients are exactly matched by donor
sales. Existing creditor payments remain in the original budget. No household
borrowing restriction is relaxed and the government incurs no debt.

Choose \(\epsilon\) small enough that every recipient remains below its
housing cap and every donor retains positive resources. The old estate-to-
housing ratio in (2) is constant, so a slack estate floor remains slack after
this scaling. The uniform subset assumptions ensure a single positive
\(\epsilon\) works with a continuum of households.

The old cohort's exact welfare change is
\[
\begin{aligned}
\Delta\mathcal W(\epsilon)
={}&K\int_A\log\left(1+\frac{\epsilon}{m_Az_i}\right)dM^O(i)\\
&+K\int_B\log\left(1-\frac{\epsilon}{m_Bz_i}\right)dM^O(i).
\end{aligned}
\tag{6}
\]
Its derivative at zero is
\[
K\left[
\frac{1}{m_A}\int_A\frac{1}{z_i}dM^O(i)
-\frac{1}{m_B}\int_B\frac{1}{z_i}dM^O(i)
\right]>0
\tag{7}
\]
by the strict resource ordering. Hence a sufficiently small positive transfer
strictly increases utilitarian welfare. Donors lose utility, so this is not a
Pareto improvement.

## 3. Why this is an equilibrium intervention, not a fixed-price conjecture

The original price path is an exact equilibrium continuation for this policy:

- Selected old owners make the optimal choices (2). Other old households face
  unchanged budgets and prices, so their original choices remain optimal.
- Current young households receive no transfer and will not receive or pay
  this expired policy when old. All their current and anticipated prices,
  incomes, rebates, feasible sets, fertility choices, and tenure values are
  unchanged. Their original choices therefore remain optimal without changing
  their choice timing.
- Aggregate old housing is unchanged by (3); young housing is unchanged.
  Housing clears at the original prices. Aggregate owner holdings and thus
  property-tax revenue are unchanged, so the original equal rebate still
  balances the property-tax account. The separate transfer account balances
  by construction.
- Aggregate current goods use, old financial saving, and terminal estate
  payments are unchanged. Thus there is no change in the net external asset
  path or external resource requirement.
- The modified old cohort exits. Current young portfolios and housing titles,
  which generate next period's old distribution, are unchanged. Individual
  estates do change, but the model expressly makes entrant wealth exogenous
  rather than a function of those estates. Future entrant types and every
  cohort size therefore remain unchanged.

This argument exhibits an equilibrium after the intervention; it does not
assume a unique equilibrium selection rule. A previously anticipated recurring
redistribution policy would alter young saving and possibly fertility and
tenure. The proof does not apply to that different experiment. It would also
need revision if individual estates determined future entrant wealth.

## 4. A compact primitive specialization

The generic theorem only needs dispersion in \(z\) among eligible old owners.
A particularly transparent sufficient specialization of the existing model is
\(\boxed{\phi=q}\). Then the stationary owner cash cost equals the housing
service cost, \(L=p\). A predecessor owner with binding mortgage satisfies
\[
qa'+\phi Ph=0
\quad\Longrightarrow\quad a'=-Ph
\quad\Longrightarrow\quad
\boxed{z=a'+Ph+y^o+T=y^o+T.}
\tag{8}
\]
Thus two separated old-income groups directly give two separated old-resource
groups; inherited titles and net mortgage repayment cancel exactly.

Retain conditions (E), (I), and (C) in
[the stationary primitives review](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/conventional_stationary_primitives_review.md).
They ensure stationary existence and a positive-mass owner set with strict
finance and slack young and old housing caps and old estate floor. At
\(\phi=q\), its income condition simplifies to
\[
k=\frac{\beta K}{1+\alpha+\theta},\qquad
qv_{0i}>kw_{0i}+(k-q)_+\bar T.
\tag{9}
\]
Here \(w_0=y^y+b\), \(v_0=y^o\), and \(\bar T\) is the existing primitive
rebate bound. Require its covered set to contain two positive-\(F\)-mass
subsets \(S_A,S_B\) with
\[
\mathop{\rm ess\,sup}_{S_A}v_0
<\mathop{\rm ess\,inf}_{S_B}v_0.
\tag{10}
\]
Finite logistic ownership probabilities give each subset positive owner mass
in every stationary equilibrium. Their stationary predecessors are current
old owners satisfying (8), so adding the common rebate preserves the strict
resource ordering. The same cap and floor bounds provide the eligible regime;
one can restrict to positive-mass subsets with uniform margins if needed.
No \(\beta\ge q\) or lifetime-housing condition (H) is required.

This is a sufficient primitive specialization, not a newly constructed
parameter family or a calibration claim. Strictly constrained predecessors
are useful here only to make old-resource dispersion transparent. The welfare
mechanism is concavity of utility in wealth and can operate without a housing
finance distortion. Young housing does not change anywhere in the theorem.

**Verification.** Equations (2), (4), and (5) were derived from the original
old-owner optimization problem and cash/estate identities. Aggregate
invariance is exact because the selected demands are linear in resources.
The equilibrium argument separately checks young incentives, fiscal balance,
future state propagation, and the model's exogenous entrant-wealth assumption.
No numerical approximation or simulation was used.

## 5. Check of the assembled short proposition

At 16:19 UTC I checked Proposition 2 and its surrounding explanation in
[the assembled utilitarian note](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex:283).
The statement already specifies an unanticipated one-time balanced transfer,
and the proof correctly uses linear old demands, their unchanged aggregates,
and the positive derivative
\(K[1/(z_L+b)-1/(z_H-b)]\). The original young timing is preserved, and the
following paragraph correctly distinguishes wealth redistribution from an
old-to-young housing improvement. The \(\phi=q\) primitive specialization
and the absence of a \(\beta/q\) restriction are correct.

Selecting equal positive submeasures is legitimate in the household continuum
even when the original income groups have different masses; it does not impose
equal total group masses. For a common uniform transfer size, the clearest
wording is to select those submeasures with a uniform positive resource gap,
positive resources, and housing-cap slack. This is an available restriction of
the positive-mass eligible groups, not an additional economic assumption.
Alternatively, balanced matched transfers may have pair-specific small sizes.
Equation (4) above supplies the financial-saving and existing-title accounting
behind the short proof's unchanged-aggregate claim. No substantive blocker was
found in the assembled old-only proposition.
