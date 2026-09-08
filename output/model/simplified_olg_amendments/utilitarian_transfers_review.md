# Transfers with conventional mortgage finance

Analytical implementation review, September 8, 2026. No simulation, numerical
reference equilibrium, model change, or protected-manuscript edit. This note
uses the conventional-finance proposal's budgets and the stationary existence
bounds already established there.

**Finding.** A complete transfer-only equilibrium improvement can be proved.
The clean broad construction needs an additional group of sufficiently wealthy
old owners whose existing housing cap binds. A simpler construction needs two
primitive equalities. Neither construction establishes a universal result for
a uniform old-age tax and young-age grant. Both retain household optimization,
the original mortgage covenant, housing caps, and estate floor.

## 1. Welfare criterion and authority

Fix fertility and tenure at their initial stationary values. Give each current
old household weight one on remaining utility and each current young household
weight one on its lifetime utility, including its private discount factor
\(\beta\). Unaffected future cohorts cancel from the welfare difference.
Thus an improvement here is utilitarian, not Pareto: taxed old households lose.

Allow announced, household-specific lump-sum taxes and grants at dates zero and
one. Their amounts may depend on predetermined identity, age, endowments and
initial tenure, but not on the household's subsequent choices. Government can
save or borrow in the external bond at price \(q\), with its position settled
at date one. Its **present-value** transfer budget must balance exactly. The
construction below with \(\phi\ge q\) only requires government saving and gives
the selected young grants at both dates; no new household loan is introduced.
Existing common property-tax rebates remain unchanged because prices,
population and occupied housing remain unchanged.

Let
\[
p=(1-q+q\tau^p)P,\quad L=(1-\phi+q\tau^p)P,\quad
K=1+\gamma+\omega_B,\quad \ell=L/p.
\]
For an eligible young owner write adult consumption and space as
\(x=c-\chi n\), \(s=h-\kappa n\), and its planned old resources as
\(z=a'+Ph+v\). Assume both age-specific housing caps and its old estate
floor are slack, and its mortgage strictly binds. Its matching current old
counterpart has resources \(z\), since the initial equilibrium is stationary.
Define
\[
m=K/z=1/c^2,\qquad
\lambda=\beta m/q,\qquad
\Lambda=1/x=\lambda+\mu,\qquad
\alpha/s=\lambda p+\mu L,\quad \mu>0.
\tag{1}
\]
Here \(m\) is the old marginal value of cash, \(\Lambda\) the young marginal
value of current cash, and \(\mu\) the multiplier on its current cash limit.
If \(\beta\ge q\), two strict gaps follow:
\[
\Lambda-m=(\beta/q-1)m+\mu>0,\qquad
\alpha/s-pm=(\beta/q-1)pm+\mu L>0.                 \tag{2}
\]
These gaps alone do not implement a transfer equilibrium: the induced housing
demands and government budget must also be reconciled.

## 2. Exact household implementation at unchanged prices

For a small positive adult-space increment \(\varepsilon\), set
\[
s_\varepsilon=s+\varepsilon,\qquad
x_\varepsilon=
\frac{L}{\alpha/(s+\varepsilon)-(p-L)\lambda},\qquad
z_\varepsilon=z.                                      \tag{3}
\]
Give this young household the predetermined grants
\[
G_\varepsilon=x_\varepsilon-x+L\varepsilon,
\qquad
J_\varepsilon=(p-L)\varepsilon/q,                     \tag{4}
\]
respectively now and when old. At small positive \(\varepsilon\), the
denominator stays positive, the cash multiplier remains positive, and the
housing caps retain their original status. Equations (3) satisfy all first
order conditions of the strictly concave household problem, so these are
actual household optima under the announced lump sums.

In particular, net financial wealth changes by
\[
\Delta a'=-\phi P\varepsilon/q.
\]
Consequently \(qa'+\phi Ph=0\) still holds exactly, and next-date resources
change by
\(-\phi P\varepsilon/q+P\varepsilon+J_\varepsilon=0\).
The same old consumption, housing, estate and nonnegative gross financial
assets remain optimal. This verifies the old-age accounting, including the
larger inherited title and mortgage repayment.

An uncapped old owner's housing demand is \(h^2=gz\), where
\[
g=\gamma/(Kp).
\]
Tax the matching current old household \(D_\varepsilon=\varepsilon/g\).
It optimally reduces housing by exactly \(\varepsilon\), as well as reducing
consumption and its estate proportionally. Current housing therefore clears
exactly. The pair's remaining present-value funding requirement is
\[
R_\varepsilon
=G_\varepsilon+qJ_\varepsilon-D_\varepsilon
=x_\varepsilon-x+p\varepsilon-\varepsilon/g.           \tag{5}
\]
This residual cannot simply be omitted.

## 3. A transfer-only equilibrium theorem with explicit restrictions

In addition to the proposal's primitive conditions supplying eligible young
owners, impose
\[
\boxed{\quad \beta\ge q,\qquad
\alpha(1+\omega_B)\le\gamma\min\{\ell,\ell^{-1}\}.\quad}       \tag{6}
\]
Suppose a positive-mass group of other current old owners is strictly at its
owner housing cap and has cash marginal utility strictly below that of every
selected matching old donor. Impose \(\phi\ge q\) if government borrowing or
future taxes on the selected young are undesirable; then \(J_\varepsilon\ge0\).

**Proposition.** A finite, sufficiently small policy of the form (4), funded
by (5) and the matching old taxes, is a competitive-equilibrium utilitarian
improvement. It transfers housing from current old to current young owners,
retains the original financial restrictions, and restores active households'
choices from date one and the inherited state from date two. Fertility and
all cohort masses remain fixed.

**Proof.** Write
\(\rho=\alpha x/s=(\lambda p+\mu L)/(\lambda+\mu)\). Along (3),
\[
A\equiv\frac{dx_\varepsilon}{d\varepsilon}
=\frac{\alpha x_\varepsilon^2}{L s_\varepsilon^2}
=\frac{\rho_\varepsilon^2}{\alpha L}>0.
\]
While the mortgage binds, \(\rho_\varepsilon\) lies between \(p\) and
\(L\). Condition (6) therefore implies
\[
A\ge\frac{\min\{p,L\}^2}{\alpha L}
\ge\frac{(1+\omega_B)p}{\gamma}=1/g-p.
\]
Hence \(R'_\varepsilon\ge0\), so (5) is a nonnegative tax requirement.
Collect it from the wealthy capped old group. A sufficiently small tax leaves
their housing demand at their cap; they optimally reduce consumption and
estates. With multiple recipients, collect the sum of their residuals and
match each recipient with equal measure of its stationary old counterpart.

Let \(m_F\) be the capped group's marginal utility of cash. For common
bundles and one unit of recipient mass, the welfare derivative at zero is
\[
\begin{aligned}
W'(0)
&=\Lambda A+\alpha/s-m/g-m_F(A+p-1/g)\\
&=(\Lambda-m)A+(\alpha/s-pm)
  +(m-m_F)(A+p-1/g)>0.                              \tag{7}
\end{aligned}
Every term is nonnegative and the first two are strictly positive by (2).
The exact utility changes are logarithms: recipient gain
\(\log(x_\varepsilon/x)+\alpha\log[(s+\varepsilon)/s]\), matching old
change \(K\log[1-\varepsilon/(gz)]\), and capped old change
\((1+\omega_B)\log[1-R_\varepsilon/(z_F-pH_O)]\), with the appropriate
mass scaling. These functions are continuously differentiable while their
displayed arguments remain positive. Thus (7) supplies a finite positive
step, restricted also by the strict mortgage and cap margins. This is a local
policy argument around analytically characterized equilibria, not existence
in a neighborhood of a numerical equilibrium.

Date-zero net tax revenue after grants equals \(qJ_\varepsilon\); investing
it in the external bond pays exactly \(J_\varepsilon\) at date one. If this
number is negative, the reverse bond trade is repaid by the date-one tax.
Housing clears at date zero by the matched reductions, and at every later
date because affected young households have exactly their original old
resources and choices. Date-one entrants face unchanged prices, transfers,
fertility and tenure. Their choices regenerate the original date-two state.
Initial old estates change; their utility and dated asset payments are included
in their optimized budgets, and estates do not finance entrants in this model.
This constructs the entire equilibrium path without a transition stability
or convergence assumption. Heterogeneous recipients are handled by integration
on positive-mass subsets with uniform strict margins. QED.

The capped funder condition can itself be imposed in primitives. Use the
proposal's price upper bound \(P_+\), eligible-type resource bound
\(M_{\mathcal S}\), and \(d_p=1-q+q\tau^p\). Assume its old estate condition
\(\omega_Bd_p>q\gamma\). Let
\(\delta=(\phi/q-1)_+\). A positive-mass funder set satisfying
\[
\boxed{\quad
y_f^o>\delta P_+H_O+
\max\left\{
\frac{Kd_pP_+H_O}{\gamma},\;
d_pP_+H_O+\frac{(1+\omega_B)M_{\mathcal S}}{qK}
\right\}
\quad}                                                        \tag{8}
\]
is sufficient. The mortgage covenant implies
\(z_f\ge y_f^o-\delta P_+H_O\). The first bound makes its old housing cap
strictly binding. The second gives
\[
m_F=\frac{1+\omega_B}{z_f-pH_O}
<qK/M_{\mathcal S}\le m_i
\]
for every selected old donor. Logistic tenure gives such funder types a
positive owner mass. These are additional, restrictive implementation
conditions; they must not be silently folded into the original direct-planner
claim. The price and mean-resource bounds in (8) are evaluated on the full
endowment distribution, including the funder group.

## 4. Simpler strict case: only current balanced transfers

There is an exact benchmark requiring no capped funders, future transfers or
government bond trading. Impose
\[
\boxed{\qquad \phi=q,\qquad \alpha(1+\omega_B)=\gamma,\qquad
\beta\ge q.\qquad}                                           \tag{9}
\]
Keep the eligible owners' binding mortgage and slack caps. Now \(L=p\),
\(z=v\), and, putting \(w_n=w-(\chi+p\kappa)n\),
\[
x=w_n/(1+\alpha),\quad
s=\alpha w_n/[p(1+\alpha)],\quad
\frac{\partial h^y}{\partial w}
=\frac{\alpha}{p(1+\alpha)}
=\frac{\gamma}{Kp}=g.
\]
Give a selected young household \(b>0\), taxing its matching current old
counterpart the same amount. Young housing rises by \(gb\), old housing
falls by \(gb\), and planned old resources remain \(v\). These are exact
finite choices. The entire future equilibrium is unchanged.

The exact welfare gain per pair is
\[
\Delta W(b)=(1+\alpha)\log(1+b/w_n)+K\log(1-b/v).
\]
Its derivative is positive throughout
\[
0<b<\frac{(1+\alpha)v-Kw_n}{K+1+\alpha}.                       \tag{10}
\]
The numerator is positive by (2). Also restrict
\(b<q(1+\alpha)v/(\beta K)-w_n\) to preserve strict mortgage binding,
and \(gb<H_O-h\) to preserve the young cap. These bounds explicitly
produce a finite positive transfer. This result is especially simple but the
two equalities in (9) are knife-edge restrictions, not a general theorem.

## 5. Limits relevant to the main note

Without capped funders, the two-date, unchanged-price, unchanged-old-resource
construction requires \(\Delta x+p\Delta h=\Delta h/g\). At \(L=p\)
this is precisely the preference equality in (9). Away from that equality,
the simple construction fails its government budget. That is an obstruction
to this implementation, not proof that every transfer policy fails.

A uniform current old tax/young gift generally changes future old resources
and equilibrium prices. Its welfare and housing directions cannot be inferred
from the fixed-price cash gap alone. The transfer-constrained planner's value
is strictly above the original equilibrium under the theorem because the
displayed policy is feasible; the direction of every optimal transfer or its
global optimum is not characterized.

When private fertility is restored, the coordinated increase in adult space
at fixed old resources raises adult consumption and fertility on the uncapped
binding branch. However, changed births alter next-period entrant mass.
Therefore the finite-tail equilibrium construction above must not be reused
as an endogenous-fertility transition theorem. The fertility response and the
demographic path require a separate statement and accounting.
