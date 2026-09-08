# Transfers, private fertility, and equilibrium paths

Analytical review, September 8, 2026. This note retains the proposed preferences,
positive child goods and space costs, heterogeneous entrants, ordinary mortgage
finance, physical housing caps, and warm-glow estates. It supplements
`simplified_olg_conventional_finance.tex` and the two conventional-finance
reviews. No simulation, equilibrium computation, or convergence argument is
used. The fixed-population welfare comparison holds individual fertility fixed;
the results below restore private fertility choice without assigning a social
value to additional births.

## 1. A cash gift has a stronger conditional result than easier mortgages

Fix a date, tenure, and the relevant price and rebate path. Let \(w\) be cash
available when young, \(v\) income plus transfers when old, \(p>0\) the current
housing-service price, and \(L>0\) the current cash requirement per housing unit.
For owners, \(p=(1+q\tau^p)P_t-qP_{t+1}\) and
\(L=(1-\phi_t+q\tau^p)P_t\); for renters set \(L=p\).
Write \(\mathcal V^m(z)\) for old-age utility maximized over consumption,
housing and estates using resources \(z\), retaining every old-age constraint.
The young problem can be written as:
\[
\max\ \log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta\log n
       +\beta\mathcal V^m(z),\qquad
c+ph+qz=w+qv,\quad c+Lh\le w,\quad h\le H_m.
\]
Thus a gift raises both current cash and lifetime resources, whereas changing
\(L\) changes the price of satisfying the cash restriction.

**Proposition.** A positive current gift, holding \(v,p,L\) and old-age
prices fixed, strictly raises conditional fertility. Conditional young housing
increases until its cap binds. This holds for every \(\chi,\kappa>0\), on
either side of \(L=p\), and with binding old housing caps or estate floors.
The statement extends to finite gifts with the same fixed-price scope.

The proof needs only the concavity of the retained old-age problem. At strict
finance and slack young housing, put \(x=c-\chi n\), \(s=h-\kappa n\),
\(\Delta=L-p\), and \(\Gamma=-\beta\mathcal V^{m\prime\prime}(z)/q^2>0\).
The binding budgets imply \(c=w-Lh\) and \(qz=qv+\Delta h\). Define the
curvatures of the two-choice problem by:
\[
\begin{aligned}
A&=L^2/x^2+\alpha/s^2+\Gamma\Delta^2,\\
B_c&=-L\chi/x^2+\alpha\kappa/s^2,\\
D_n&=\chi^2/x^2+\alpha\kappa^2/s^2+\vartheta/n^2,
\qquad J=AD_n-B_c^2>0.
\end{aligned}
\]
Implicit differentiation gives two strictly positive responses:
\[
h_w=\frac{LD_n+\chi B_c}{x^2J}>0,\qquad
n_w=\frac{LB_c+\chi A}{x^2J}>0,
\]
because their numerators simplify to:
\[
LD_n+\chi B_c=\frac{L\vartheta}{n^2}
 +\frac{\alpha\kappa(\chi+L\kappa)}{s^2}>0,\qquad
LB_c+\chi A=\frac{\alpha(\chi+L\kappa)}{s^2}
 +\chi\Gamma\Delta^2>0.
\]
If young housing binds, \(h=H_m\) and
\(n_w=\chi/(x^2D_n)>0\). If finance is slack and young housing is slack,
the three current adult-consumption, adult-space and fertility expenditures
have fixed expenditure shares; their common expenditure strictly increases
with wealth by strict concavity of old utility. If only young housing binds,
the consumption/fertility first-order conditions likewise give \(n_w>0\).
Old values are continuous, concave and piecewise smooth, and unique household
choices join continuously at all constraint changes. These observations cover
finite gifts and every cap regime.

## 2. Repayment changes the result

Let a marginal current transfer \(dg\) carry a future tax with present value
\(r\,dg\), so \(dw=dg\) and \(q\,dv=-r\,dg\). Holding prices fixed on
the same strictly constrained, young-uncapped branch gives the exact test:
\[
\frac{dn}{dg}
 =\frac{(LB_c+\chi A)/x^2+r\Gamma\Delta B_c}{J}.
\]
A pure gift has \(r=0\); a loan repaid at the bond return has \(r=1\).
When \(\Delta B_c\ge0\), repayment reinforces the conditional fertility
effect. When \(\Delta B_c<0\), fertility rises exactly when:
\[
r<\frac{LB_c+\chi A}{-x^2\Gamma\Delta B_c}.
\]
These are household-allocation tests expressed through observed bundles and
the known old-value curvature, not a universal primitive loan theorem.

The rational counterexample in `conventional_fertility_review.md` makes the
distinction concrete. It has \(p=1,L=1/10,\alpha=\vartheta=\kappa=1\),
\(\chi=1/10\), \(x=1,s=100/91,n=100/101\), and discounted old resources
\(S=10/9\) with old log weight one. At exactly this admissible household:
\[
n_{g,\mathrm{gift}}=\frac{23123000}{210723483}>0,\qquad
n_{g,\mathrm{fair\ loan}}=-\frac{36516490}{210723483}<0.
\]
Small gifts and small loans therefore have opposite fertility effects despite
unchanged positive child costs. A young cap makes the direct fertility effect
positive while finance binds, since repayment then leaves current housing
fixed; crossing into slack finance still requires its own comparison.

A current tax on today's old households need not be a repayment obligation for
today's young recipients. But a permanent age-transfer rule can tax those
recipients next period. For example, grants \(g_t\) per young financed by a
uniform old tax require \(d_t=(Y_t/O_t)g_t\); recipients anticipate
\(d_{t+1}=\nu\bar n_tg_{t+1}\). The corresponding change in \(v\), as well
as ordinary rebates and prices, belongs in their fertility comparison.

**Coordinated transfers holding old resources fixed.** There is another useful
conditional direction. Hold \(z\) fixed and finance strict. When young housing
is uncapped, its first-order condition has
\(\alpha/s=L/x+\lambda(p-L)\), with
\(\lambda=\beta\mathcal V^{m\prime}(z)/q\) fixed. Consequently:
\[
dx=\frac{\alpha x^2}{Ls^2}\,ds>0,\qquad
dn=\frac{n^2}{\vartheta}
       \left(\frac{\chi}{x^2}\,dx+\frac{\alpha\kappa}{s^2}\,ds\right)>0
\quad\text{when }ds>0.
\]
Both current consumption and housing rise. The necessary transfers obey
\(dw=dc+L\,dh>0\) and \(q\,dv=(p-L)dh\). When \(L<p\), this direction
requires a future subsidy as well as a current grant. It is not a pure gift
or a fair loan. With a binding young cap, a current grant at fixed \(z\)
instead raises fertility through consumption alone.

## 3. Tenure choice is a separate restriction

For a common current gift, conditional values satisfy \(W_w^m=1/x_m\).
The logistic tenure rule therefore gives:
\[
\frac{d\bar n_i}{dg}=\pi_i n_{O,w}+(1-\pi_i)n_{R,w}
 +\frac{\pi_i(1-\pi_i)}{\sigma_\xi}
       (n_O-n_R)(1/x_O-1/x_R).
\]
The first two terms are positive; the selection term need not be. A sufficient
condition is \((n_O-n_R)(1/x_O-1/x_R)\ge0\), which can be checked directly
from each type's conditional bundles. Integrating preserves this result for
any fixed heterogeneous entrant distribution.

This qualification is substantive. Take slack caps, strict renter finance at
\(L=p\), and then an owner cash cost \(L=p-\delta\) with small
\(\delta>0\). Write \(E=1+\alpha+\vartheta\),
\(B=\beta(1+\gamma+\omega_B)\), and
\(t_0=Bw/(Eqv)\in(0,1)\). At \(\delta=0\), \(x=w/E\).
The owner's reduction in cash cost increases both adult consumption and
fertility: the first derivatives have
\(x_\delta=x t_0H(p)/E>0\) and
\(n_\delta=\vartheta[x_\delta/(\chi+p\kappa)
 +x\kappa(1-t_0)/(\chi+p\kappa)^2]>0\), where
\(H(p)=\alpha/p+\vartheta\kappa/(\chi+p\kappa)\).
Thus the displayed selection term is negative. Set the taste location to make
\(\pi=1/2\); a sufficiently small positive taste scale makes that term
dominate the conditional gains. A nondegenerate bounded distribution with
\((w_i,v_i)=s_i(w_0,v_0)\) preserves the example: real choices scale with
\(s_i\), conditional value differences and the selection product do not.
Hence even a common fixed-price gift can reduce cohort-average fertility.

An **owner-contingent grant available before purchase** expands only the owner
menu, so ownership rises. If baseline \(n_O\ge n_R\), conditional-owner
monotonicity then gives a finite aggregate fertility increase. One primitive
sufficient baseline ordering, with slack old caps, is the earlier condition:
\[
L\le p,\qquad
\frac{p}{\alpha+\beta(1+\gamma+\omega_B)}
 \le\frac\chi\kappa\le\frac p\alpha.
\]
A grant to households whose tenure is already fixed has no selection term.
Targeting and timing must be specified before either claim is used.

## 4. What can be said along an admissible equilibrium path

There is an exact finite test that includes prices, future taxes, housing caps,
and consumption adjustment. Compare a household's baseline fertility \(n_0\)
with its new privately chosen bundle \((c_1,h_1)\), holding the fertility taste
\(\vartheta\) fixed. If \(c_1>\chi n_0\) and \(h_1>\kappa n_0\), then:
\[
n_1\ge n_0
\quad\Longleftrightarrow\quad
\frac\chi{c_1-\chi n_0}+\frac{\alpha\kappa}{h_1-\kappa n_0}
\le
\frac\chi{c_0-\chi n_0}+\frac{\alpha\kappa}{h_0-\kappa n_0}.
\]
If either feasibility inequality fails, \(n_1<n_0\). The proof is that
\(\vartheta/n-\chi/(c-\chi n)-\alpha\kappa/(h-\kappa n)\) is strictly
decreasing in \(n\). Increasing both total consumption and housing is a
simple sufficient condition, but the exact test allows one to fall. Apply
this to matched types and taste draws, or use the exact finite tenure identity:
\[
\Delta\bar n_i=\pi_1\Delta n_O+(1-\pi_1)\Delta n_R
                 +(\pi_1-\pi_0)(n_{O0}-n_{R0}).
\]
These date-by-date statements require no convergence. They require actual
equilibrium bundles, or independent restrictions that imply their ordering.

The distinction between a dated comparison and a transition theorem is:

| Claim | Required information |
|---|---|
| Fertility at a specified date on two given admissible paths | The preceding bundle/selection test; no convergence |
| Cohort population at a finite date | Cumulative fertility differences; no convergence |
| Positive stationary population limits and their comparison | Convergence to positive stationary equilibria must be assumed or proved |
| Policy determines fertility at every future date | Control of the equilibrium price, tax and tenure responses along the whole path; convergence alone is insufficient |

For paths sharing \(Y_{t_0}\), cohort accounting is exactly:
\[
\frac{Y_t^1}{Y_t^0}=\prod_{s=t_0}^{t-1}
       \frac{\bar n_s^1}{\bar n_s^0},\qquad O_{t+1}^j=Y_t^j.
\]
This can order finite populations even when neither path converges. If both
converge to positive stationary equilibria, both fertility limits equal
\(1/\nu\); their limiting population ratio instead equals the inverse ratio
of stationary young-plus-old housing per cohort. A higher limiting population
does not imply higher fertility at every intervening date.

Each young household directly needs only
\((P_t,P_{t+1},P_{t+2},T_t,T_{t+1})\), its announced age-specific transfers,
and current policy parameters. The second future price enters old services
and the estate restriction. Nevertheless, determining those near-future
prices as equilibrium outcomes can require the entire announced policy and
equilibrium path. Assuming prices become stationary after a chosen date is
not a proof of convergence. Finally, a fixed-fertility welfare construction
that restores baseline states after two dates cannot be carried over without
rechecking it: private fertility changes the next cohort's mass.

Verification: the gift numerators and the rational gift/loan comparison were
checked with exact arithmetic; the cap cases and tenure counterexample follow
from the displayed first-order conditions. No equilibrium or welfare ranking
with endogenous population is claimed.
