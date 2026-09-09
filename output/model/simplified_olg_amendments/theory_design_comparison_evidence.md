# Two-hour comparison of illustrative housing models

September 9, 2026. **Unadopted proposal; the main deck and maintained theory have not changed.**

Reader: `output/pdf/theory_design_comparison.pdf` (eight-page target). Source: `output/model/simplified_olg_amendments/theory_design_comparison.tex`.

This supporting record consolidates the calculations and independent reviews from the 17:54–19:54 UTC design pass. Read the PDF first. The component calculations below were completed sequentially: early statements that fertility or transition is not established describe that component's scope, and are extended by the later components. Nothing here replaces the September 14 seminar deck until Tommaso chooses a specification.

## Lead assessment

The strongest candidate retains two utility ages, the original logarithmic consumption/housing/fertility preferences, a common divisible housing stock, rental and owner upper size limits, 80% mortgage LTV, positive rebated property taxes, and heterogeneous income and cash. It distinguishes cash available at the initial home purchase from later working income. Both current ages consume at the same period end. Young households rent at their housing limit and save; old households choose larger owned homes, can sell and resize, and need not have positive financial estates.

For a finite region of primitive parameters, the stationary price is explicit and every household's tenure and fertility deviations have been checked. Household fertility is a quadratic; heterogeneous replacement fertility is imposed by a unique scalar preference calibration at a specified cohort normalization. This calibration is implicit, and the result does not exclude other stationary regimes globally. The finite existence family is not a numerically solved reference followed by continuity.

A dated planner maximizing the equally weighted utilities of current households chooses all current consumption, housing and tenure. Fertility is fixed first. Individual young continuation opportunities and old net estates are held fixed, and current transfers sum to zero including rental intermediaries. The planner allocates more housing to the young. When it also chooses fertility, an explicit goods-versus-space condition ensures higher mean fertility. A stronger primitive restriction verifies this throughout the displayed family without imposing beta R >= 1.

The same financial and tenure model also has a local perfect-foresight transition with endogenous prices and rebates. A small permanent fertility-taste decline lowers impact fertility and terminal adult-household population. A small permanent rebated property-tax increase introduced at any finite date on that nearby baseline raises fertility on impact and raises the terminal population. Both terminal fertility rates equal replacement. Young housing stays at the rental limit: this tax result does not implement the planner's home upgrades. The tax price falls on impact and rises at the destination in the explicit family. No monotone adjustment or all-date fertility ordering is claimed.

## Assumptions that need an author decision

1. Initial cash includes every receipt already received. Labor income Y genuinely arrives after the housing decision; the household cannot revisit the same period's home purchase after receiving it.
2. Tenure can change between ages. The first illustration omits ownership tastes, while retaining heterogeneous cash, income and fertility. All young rent in the verified example; it does not have a positive realized young mortgage share.
3. Current consumption and property taxes are paid from liquid receipts before title liquidation for BOTH ages. A young mortgage can then be repaid using house-sale proceeds; total wealth after settlement must be nonnegative before the next rebate. The final reader uses the exact two solvency inequalities below, weakening the initial calculation's stronger before-sale financial-asset bound without changing the verified equilibrium. All young rent in that equilibrium.
4. The physical environment has a common divisible housing stock and tenure-specific upper size menus. The planner may reassign tenure; it is not confined to immutable rental and owner inventories.
5. This is a small open economy with outside bond and intermediary finance. Estates provide warm-glow utility and do not endogenously finance the next entrants' cash. The distribution of entrant income and cash is the same exogenous distribution each date. Fertility-based type transmission and an internal bequest distribution would require a different transition state.
6. The finite existence example restricts gross resource dispersion to 1% for the joint-fertility and heterogeneous-transition result and child goods costs to 2% of mean lifetime stage resources at the maximum rental-space fertility. It also requires low purchase cash relative to later resources. These are sufficient analytical restrictions, not an empirical assessment.
7. The planner result is dated utilitarian welfare with fixed future opportunities, not a Pareto improvement, an infinite-horizon social optimum, or tax implementation. The local transition neighborhood is reference-dependent, not a common numerical shock radius for the entire parameter family.

## Alternatives investigated

| Design | Analytical result | Main economic cost / reason not to lead with it |
|---|---|---|
| Original income timing plus a minimum owner size | Explicit fixed-fertility stationary price and age allocation with mobile tenure | The verified family still uses strong old-income and estate demand. The owner minimum changes the housing menu. |
| Add a mature working age | Finite stationary existence with conventional LTV, a full dated planner result, quadratic/cubic household choices | Adds an age and a financial restriction; the useful existence family needs high mature income. A children-only-when-young variant simplifies the fertility calculation. |
| Original two-age logs, weaker estate assumptions | Covers both estate regimes and gives a conditional adult-space welfare comparison | In the checked moderate-taste family old gross homes are much smaller than young homes, so it does not illustrate the intended age pattern. |
| General concave preferences | Increasing differences between child needs and housing support the marginal housing comparison | Monotonicity and concavity alone do not imply greater young housing needs; consumption and aggregate allocation still need their own conditions. |
| Fixed small and large products | Clean swap/sorting and fertility comparisons | Fixing occupied product counts can fix the population and prevent the desired change in terminal population. |
| Real moving or liquidation costs | A direct reason for persistent large old homes | Costs also enter physical planner feasibility, so persistence does not itself imply inefficient retention. |
| Quasilinear consumption | Very short complete household and price formulas | Removes the income effects central to the question; the verified positive-consumption construction adds an unwanted discount restriction. |

These alternatives are not impossibility results. The common-calendar model is recommended for discussion because it delivers the requested combination with fewer additional mechanisms. The earlier three-age full calculation remains separately available in the same output folder; none of these designs has been adopted.

## Verification and navigation

The lead checked the household and settlement accounts, the complete tenure comparison, the finite parameter bounds, the joint fertility comparison, and the transition algebra. Separate agents derived the static, fertility, and transition components. A hostile reviewer checked the independent derivations and qualifications. Verification uses exact algebra and rational inequalities; no numerical equilibrium search or simulated reference proves the main claims. The reader PDF was compiled twice and every page rendered for inspection. The final review status is reported in the independent-review part below.

- Part I: dated accounts, explicit stationary price, global tenure verification, finite primitive family, full fixed-fertility planner.
- Part II: full joint-fertility planner, quadratic private fertility, endogenous heterogeneous stationary closure.
- Part III: endogenous-price transition, stationary endpoints, local stability and shock signs, then exact heterogeneous extension.
- Part IV: alternative repayment calendars, followed by the final common-order settlement specified here and audited in Section 12.4.
- Part V: independent audits, including corrections and scope restrictions.

Cross-sectional means in these calculations are integrals over a fixed type distribution; expectation notation appearing in a derivation is not uncertainty or stochastic perfect foresight.


## Final financing specification in the reader

The initial stationary derivation imposed \(a'_O\ge0\). This is sufficient but unnecessarily strong. The final reader keeps the same consumption-before-title-sale order for both ages and uses, in addition to the origination deposit constraint,

\[
a'_O+P_{t+1}h_y\ge0,\qquad qa'_O+\phi P_t h_y\ge0.
\]

Here \(a'_O\) is the young owner's net financial position just before title liquidation, after subtracting the matured mortgage. The first inequality requires nonnegative total wealth after repayment; the second permits consumption and property taxes to be paid from liquid receipts before the house is sold. These are financial feasibility conditions, not additional conditions on the theorem's parameter family.

Given a candidate satisfying the consolidated budget and these constraints, choose

\[
d=\max\{0,P_t h_y-w-T_t,-qa'_O\},\qquad
k=w+T_t+d-P_t h_y.
\]

The deposit and second solvency conditions imply \(d\le\phi P_t h_y\) and \(k\ge0\). Young liquid income and bond receipts first cover consumption and taxes; their remaining liquid assets equal \(a'_O+d/q\ge0\). Sell the title if needed and repay \(d/q\), leaving \(a'_O+P_{t+1}h_y\ge0\). The next entry rebate arrives after this clearing. Old households take no new loan and use the same consumption-before-sale ordering, so their old financial-estate floor remains unchanged. This avoids granting young households special access to house-sale receipts for current consumption.

A renter replicates every such owner plan with financial saving \(a'_R=a'_O+P_{t+1}h_y\). Its present housing cost is \((1+q\tau_t)P_t-qP_{t+1}\), giving identical current consumption, fertility, and old resources. All owned homes affordable to the young are below the rental cap in the verified family. Enlarging the young-owner feasible set this way therefore leaves every equilibrium, planner, fertility and local-transition result unchanged. Section 12.4 below independently verifies this implementation. The earlier terminal-only relaxation described in Part IV is an alternative, not the final reader's preferred ordering.


---

# Part I. Common calendar and stationary equilibrium

# A common consumption calendar with origination cash

September 9, 2026. Separate, unadopted two-utility-age design. This note
proves an exact fixed-fertility stationary competitive branch, conditional
on the stated primitives. It does not prove endogenous-fertility equilibrium,
a transition, or global equilibrium uniqueness. Only this scratch file is
owned by this task.

## Result and departures from the maintained model

A synchronized end-of-stage consumption calendar admits an exact
heterogeneous equilibrium with positive property tax and fully funded equal
rebates, 80 percent LTV, all young households renting at the rental size cap,
and all old households owning larger homes. No owner minimum size is needed.
The old estate floor binds in the explicit family. The private discount
factor ranges across both sides of the bond discount.

This uses the original two log-utility ages and the original physical size
menus, but makes three substantive changes explicit. Working income received
after closing finances consumption later in the current stage; both ages
consume at stage end. Tenure can change with age, and the first illustration
has no ownership taste shock. Finally, a young mortgage is fully repaid at
the age boundary with nonnegative financial assets before carrying the house
into old age. That last condition is stronger than repayment alone if a
contemporaneous house sale could fund repayment. It is slack for every
household in the constructed branch, since all young rent and save strictly.

Fertility \(n_i\) is fixed in the equilibrium theorem here. Require
\(0<\kappa n_i<H_R\), and, for a stationary population,
\(\nu\mathbb E n_i=1\). A nonempty such schedule requires
\(\nu>\kappa/H_R\). The young and old cohort
masses are normalized to one, with the same stationary endowment
distribution. Hence \(H\) below is total stock per young cohort, and there
are two current adult households per young-cohort unit.

## 1. Calendar, fiscal funding, and individual accounts

At the beginning of stage \(t\), both current ages choose the home occupied
over \(t\) through \(t+1\). All cash already received is available at that
purchase. Current young and old consumption, and the current old estate,
are delivered at \(t+1\). A current young household's own old consumption
is delivered at \(t+2\). Thus all current consumption uses a common delivery
date. The bond price over one stage is \(q\in(0,1)\); the private utility
discount between its two ages is \(\beta>0\).

At stationary price \(P>0\) and tax rate \(\tau>0\), define
\[
A=1+q\tau,\qquad D=A-q>0.
\tag{1}
\]
Owners pay \(\tau Ph\) at stage end. Competitive intermediaries charge
end-of-stage rent \(DPh/q\). They buy at \(P\), collect rent, pay the tax,
and sell at \(P\); their net end receipt is \(P/q\), exactly the required
bond return.

All current households receive the same upfront rebate
\[
T=\frac{q\tau PH}{2}.
\tag{2}
\]
The authority borrows \(2T\) initially and repays \(2T/q=\tau PH\) from
the stage-end property tax. Rebates are therefore an internal fiscal
transfer, not a new resource. The exogenous young cash \(w_i>0\) includes
all other receipts already received; \(Y_i>0\) is genuinely later working
income. Retirement labor income is zero in this illustrative family.

A young owner originates a mortgage \(d\le\phi Ph\) and buys bonds
\(k\ge0\):
\[
k+Ph=w_i+T+d,
\qquad c+a'+\tau Ph=Y_i+k/q-d/q,
\qquad a'\ge0.
\tag{3}
\]
Equivalently,
\[
qc+qa'+APh=w_i+T+qY_i,
\qquad h\le\frac{w_i+T}{(1-\phi)P},\quad a'\ge0.
\tag{4}
\]
The deposit cap and the end-stage financial lower bound are both required.
A young renter satisfies
\[
qc+qa'+DPh=w_i+T+qY_i,
\qquad a'\ge0,\quad h\le H_R.
\tag{5}
\]
Housing is committed before the later income receipt. Permitting a new
purchase after that receipt and before delivering the same young housing
service would reopen the origination constraint and change this model.

An old household begins with resources
\[
Z=a'+\mathbf 1\{\text{previously owner}\}Ph_y+T.
\tag{6}
\]
It can sell its inherited home and choose either tenure. An owner has no
new borrowing and obeys
\[
q c_o+qe+DPh_o=Z,
\qquad e\ge Ph_o,\qquad h_o\le H_O.
\tag{7}
\]
Indeed, buying the old house leaves \(k_o=Z-Ph_o\ge0\), and
\(c_o+a_e+\tau Ph_o=k_o/q\), with financial estate \(a_e\ge0\) and
\(e=a_e+Ph_o\). Equation (7) is exactly this account. Its estate floor and
positive consumption imply \(APh_o<Z\), so the initial no-loan purchase
condition is automatically satisfied. An old renter instead has
\(q c_o+qe+DPh_o=Z\), \(h_o\le H_R\), with no housing estate floor.

Preferences are
\[
\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta\log n
+\beta[\log c_o+\gamma\log h_o+\omega\log e].
\tag{8}
\]
The estate remains warm glow and does not fund the next entrant's \(w_i\).
Goods, bonds, and intermediary finance retain the outside-sector convention
of the maintained small open economy.

## 2. Old-owner solution and a candidate price

Define
\[
K=1+\gamma+\omega,\quad
a=\min\left\{\frac\gamma D,\frac{\gamma+\omega}{A}\right\},
\quad b_e=\max\{\omega,qa\}.
\tag{9}
\]
If its size cap is slack, an old owner chooses
\[
c_o=\frac{Z}{qK},\qquad
h_o=\frac{aZ}{KP},\qquad e=\frac{b_eZ}{qK}.
\tag{10}
\]
The financial-estate floor is slack when \(\omega D\ge q\gamma\) and
binds otherwise. In either case
\[
V_O(Z)=K\log Z+C_O,
\quad C_O=-K\log K-(1+\omega)\log q-\gamma\log P
+\gamma\log a+\omega\log b_e.
\tag{11}
\]
This value relaxes the owner size cap. It is an upper bound on the
owner alternative, not on the renter alternative.

Write \(r=H_R\), \(m=H-r\), and define fixed-fertility net resources
\[
b_i=w_i/q+Y_i-\chi n_i,\quad \bar b=\mathbb E b_i,
\quad B=(1+q)/q,\quad G=1+\beta K.
\tag{12}
\]
If young rent at \(r\) and become uncapped old owners, their candidate is
\[
x_i=c_i-\chi n_i
=\frac{b_i+BT-(D/q)Pr}{G},
\qquad Z_i=\beta Kx_i,
\qquad a_i'=\beta Kx_i-T,
\tag{13}
\]
\[
c_{o,i}=\frac\beta qx_i,
\qquad h_i^o=\frac{a\beta}{P}x_i.
\tag{14}
\]
The old size cap, young rental cap, saving bound, and global tenure choices
must still be checked; they are checked below.

Market and fiscal clearing reduce to the explicit equation
\[
J=Gm+a\beta\left[\frac Dq r-\frac{(1+q)\tau H}{2}\right],
\qquad
\boxed{P=\frac{a\beta\bar b}{J}},\qquad
T=\frac{q\tau PH}{2},
\tag{15}
\]
provided \(J>0\). Heterogeneity survives in particularly simple form:
\[
\boxed{h_i^o=m+\frac JG\left(\frac{b_i}{\bar b}-1\right)}.
\tag{16}
\]
These equations solve prices and rebates jointly. They do not hold the
rebate fixed while solving the price.

## 3. Global tenure verification without an owner minimum

First require
\[
\frac{w_{\max}}P+\frac{q\tau H}{2}<(1-\phi)r.
\tag{17}
\]
Every feasible young owner then occupies strictly less than \(r\). A young
renter can replicate any such owner's \((c,h,Z)\) by setting
\(a'_R=a'_O+Ph\ge0\). Equations (4)–(6) verify the identity exactly.
Both have access to the same future tenure menu. With no ownership taste,
ownership cannot improve on the optimal young renter plan. Once the
strictly optimal renter size is \(r\), the comparison is strict. This
argument would not exclude ownership under an unbounded positive taste
shock; no such shock is included in this illustration.

The current rental cap is optimal within the relaxed future-owner problem
if
\[
\alpha x_i>\frac{DP}{q}(r-\kappa n_i).
\tag{18}
\]
Because (11) is a global upper bound on future ownership, strict concavity
and attainment of the old unconstrained owner size make (13) globally
optimal among all plans ending as owners.

The future-renter alternative must be checked separately. Let
\(K_R=1+\omega\), \(G_R=1+\beta K_R\). Its **joint concave two-age renter
problem**, with both housing caps retained but the saving lower bound
temporarily relaxed, has both caps binding whenever (18) holds and
\(h_i^o>r\). Its solution is
\[
x_i^{RR}=\frac{Gx_i-DPr}{G_R},\qquad
c_o^{RR}=\frac\beta q x_i^{RR},\qquad
e^{RR}=\omega c_o^{RR},\qquad h_y^{RR}=h_o^{RR}=r.
\tag{19}
\]
To verify the active set rather than extrapolate a capped value function,
use \(aD\le\gamma\):
\[
\frac{x_i^{RR}}{x_i}
=\frac{G-a\beta D(r/h_i^o)}{G_R}>1,
\tag{20}
\]
so its young cap condition follows from (18). Its old cap condition
\(\beta\gamma x_i^{RR}\ge DPr\) reduces to
\(h_i^o/r\ge aD/\gamma\), which follows from \(h_i^o>r\).
Below its saving is positive as well. Even without that extra check,
(19) would supply a valid upper bound by relaxing the saving restriction.

The exact value difference between the candidate ending as an owner and
(19) is, writing \(t_i=h_i^o/r\),
\[
\boxed{\Delta_i=
\beta\gamma\log t_i+\beta\omega\log(b_e/\omega)
-G_R\log\left(\frac{G-a\beta D/t_i}{G_R}\right).}
\tag{21}
\]
A strictly positive value proves the global future-tenure choice. No claim
that (11) upper-bounds low-wealth renter utility is needed. A useful purely
primitive sufficient bound is
\[
\Delta_i\ge\beta\left[
\gamma(\log t_i-1)+\frac{aD}{t_i}
+\omega\log(b_e/\omega)\right],
\tag{22}
\]
from \(\log(1+z)\le z\). The bracket is increasing for \(t_i>1\), because
\(aD\le\gamma\).

## 4. A nonempty heterogeneous family with positive taxes

Normalize housing units and the stock per young cohort by
\[
q=\frac12,\quad\phi=\frac45,\quad r=1,\quad H=\frac52,
\quad H_O=2,
\tag{23}
\]
and let
\[
\alpha=\gamma=1,\quad\omega=\frac14,
\quad\beta\in\left[\frac3{10},1\right],
\quad\tau\in\left[\frac1{100},\frac1{50}\right].
\tag{24}
\]
Allow any compact, genuinely heterogeneous positive \((w_i,Y_i,n_i)\)
distribution satisfying fixed \(0<\kappa n_i<1\), \(\nu\mathbb E n_i=1\),
and
\[
\left|\frac{b_i}{\bar b}-1\right|\le\frac1{20},
\qquad \bar b\ge45w_{\max}.
\tag{25}
\]
Here \(\chi,\kappa>0\) and \(\vartheta>0\) are retained. Because fertility
is fixed, \(\vartheta\log n_i\) does not affect the allocation comparisons.
These inequalities concern net stage income relative to liquid resources
at closing. They impose a strong but transparent liquid-wealth restriction;
they do not withhold any income already received. The family is nonempty:
choose a compact heterogeneous \(w_i\), a much larger, mildly heterogeneous
positive \(b_i\), and set \(Y_i=b_i-w_i/q+\chi n_i>0\).

For every member of this family, \(K=9/4\), the old estate floor binds,
and
\[
a=\frac{5}{4(1+\tau/2)},\qquad b_e=qa,\qquad
\frac{505}{804}\le aD\le\frac{255}{404},
\quad\frac{b_e}{\omega}\ge\frac{250}{101}.
\tag{26}
\]
Also
\[
J=\frac32G+a\beta\left(1-\frac{7\tau}{8}\right)>0,
\qquad
\frac JG<\frac32+\frac{250}{201}\frac4{13}<\frac{19}{10}.
\tag{27}
\]
Equation (16) therefore gives the uniform, strict bounds
\[
\frac75<\frac{281}{200}<h_i^o
<\frac{319}{200}<\frac85< H_O.
\tag{28}
\]
All old owners occupy more than young renters, and their physical cap is
slack. For origination,
\[
\frac{J}{a\beta}\le\frac{31033}{4000},\qquad
\frac{w_{\max}}P+\frac{q\tau H}{2}
\le\frac{33283}{180000}<\frac15.
\tag{29}
\]
Thus every feasible young owner is restricted to
\(h<33283/36000<r\). The rebate has been included in this comparison.

The young renter cap is strictly binding: its condition is
\(h_i^o>a\beta D(r-\kappa n_i)/q\), while the right side is at most
\(255/202<7/5\). Young saving is strictly positive since
\[
\frac{Z_i}{P}=\frac{K h_i^o}{a}>
\frac{12663}{5000}>\frac1{80}\ge\frac TP.
\tag{30}
\]
The RR candidate also saves: (20) implies
\(\beta K_Rx_i^{RR}/P>1407/1000\), whereas \(T/P\le1/80\).
Its saving is \(a'_{RR}=\beta K_Rx_i^{RR}+DPr-T>0\).

Finally, (22) is uniformly positive. For \(t\ge7/5\), the bracket is at
least
\[
\log(7/5)-1+\frac{505}{804}\frac57
+\frac14\log(250/101)>0.
\tag{31}
\]
This is an exact inequality, not a numerical value comparison. Using
\(\log z\ge2[u+u^3/3]\), \(u=(z-1)/(z+1)>0\), its left side is bounded
below by
\[
\frac{109}{324}-1+\frac{2525}{5628}
+\frac12\left[\frac{149}{351}
+\frac13\left(\frac{149}{351}\right)^3\right]
=\frac{612657517}{60843676257}>0.
\tag{32}
\]
Consequently the future-owner plan beats the global RR alternative,
the old tenure choice is optimal, and replication rules out current
ownership. All individual, rental-pricing, fiscal, and housing-clearing
equations hold at (15). The stationary old asset distribution is exactly
the preceding young cohort's \(a_i'=\beta Kx_i-T\).

This proves existence of a positive stationary competitive equilibrium
with the stated regime, and gives its price explicitly. The price is
unique within this regime. Other prices and tenure regimes have not been
excluded, so global stationary uniqueness is not claimed.

## 5. What this does and does not give the planner

The synchronized calendar removes the previous alternative's mismatch
between young and old consumption delivery dates. With current young
future resources \(Z_i\), current old estates \(e_i\), and future tax
receipts held fixed, reallocations at the beginning of the stage satisfy
the same current settlement account for either age:
\[
q\,\Delta c_i+DP\,\Delta h_i=\text{current transfer}_i.
\tag{33}
\]
Thus fixed total stage-end goods and current housing imply zero aggregate
transfer. A planner that can relax individual financing can implement the
static resource reallocation before home commitments, with the original
physical size caps retained. The old estate floor is a private financial
restriction and may bind initially; preserving the net estate does not
require freezing old housing. This is the beginning-of-stage planner over
current-stage service and goods delivery, not an intervention after those
homes have already been committed.

The stock must be common, convertible floor space, with no fixed rental and
owner sector inventories. A planner that **freezes realized tenure** cannot
raise young housing here: every young household is already at \(H_R\).
The intended full planner must explicitly be allowed to assign some current
young households to ownership, retaining the physical \(H_R,H_O\) menus.
That permission does not follow merely from relaxing financial inequalities
if the earlier benchmark fixed tenure.

The rental financier must enter the settlement ledger. For a young owner
whose future \(Z\) is fixed, \(\Delta a'=-P\Delta h_y\); for an old owner
with fixed \(e\), \(\Delta a_e=-P\Delta h_o\). Their renter counterparts
have respectively \(\Delta a'=0\) and \(\Delta a_e=0\). Thus households'
total end-stage financial claims change by \(-P\Delta H_O\), where
\(\Delta H_O\) is the change in their owned housing. The rental
intermediary's end-stage financial claims, after rent and tax settlement
and before title liquidation, change by \(-P\Delta H_R\). Indeed its new
end debt is \(-P\Delta H_R/q\), while net rent after tax adds
\((1/q-1)P\Delta H_R\). Since \(\Delta H_O+\Delta H_R=0\), the consolidated
end-claim change is zero. Omitting the intermediary would give incorrect
household-only bond accounting when housing changes tenure.

The quantitative primitive region remains restrictive: all young rent,
there are no ownership tastes, net later income is large relative to closing
cash, and fertility is fixed. A full endogenous-fertility planner, private
fertility equilibrium, global equilibrium uniqueness, and any policy
transition require separate results. Positive tax here funds an actual
rebate in equilibrium; it is not yet a signed tax comparative static.
No realized young household carries a mortgage in this branch. Credit
matters by restricting the owner alternative below the rental ceiling;
the rental size cap is the binding constraint on the chosen allocation.

## 6. Uniform tenure verification for a separate fertility extension

The fixed-fertility price theorem does not itself establish an equilibrium
with endogenous fertility. Its tenure comparison can nevertheless be made
uniform over fertility deviations. Define gross resources
\(g_i=w_i/q+Y_i\) and \(\bar g=\mathbb E g_i\), and replace (25) by
\[
|g_i/\bar g-1|\le1/100,\qquad
\chi r/\kappa\le\bar g/50,\qquad
\bar g\ge46w_{\max}.
\tag{34}
\]
At any actual candidate with \(0<\bar n<r/\kappa\), its mean net resource
\(\bar b=\bar g-\chi\bar n\) lies in \([.98\bar g,\bar g]\) and exceeds
\(45w_{\max}\). For any individual's hypothetical
\(n\in(0,r/\kappa)\),
\[
g_i-\chi n\in[.97\bar g,1.01\bar g],\qquad
\left|\frac{g_i-\chi n}{\bar b}-1\right|
\le\frac3{98}<\frac1{20}.
\tag{35}
\]
Hold the actual candidate \(P=a\beta\bar b/J\) and rebate fixed. Its
conditional future-owner plan for this deviating fertility still satisfies
\[
h_i^o(n)=m+\frac JG\left(\frac{g_i-\chi n}{\bar b}-1\right).
\tag{36}
\]
The denominator is the actual mean, not a recomputed mean for the individual
deviation. Every cap, saving, origination, and global RR comparison above
therefore holds uniformly for every feasible fertility choice. The joint
future-owner problem is strictly concave in adult consumption, adult space,
fertility, and future resources, with linear constraints. Its optimizer
dominates every joint RR plan; young-owner replication also covers all
feasible fertility choices. This verifies the global tenure regime of an
endogenous-fertility candidate satisfying the stated price/resource relation.
The fertility first-order conditions, demographic replacement, and existence
of a candidate satisfying them remain separate equations to solve.

## Verification

All displayed constants in (26)–(32) were checked with exact rational
arithmetic. No numerical household solver, equilibrium run, or root search
was used. The global tenure proof uses the joint renter problem and explicit
value comparison rather than extrapolating a capped continuation formula.


---

# Part II. Joint fertility and endogenous stationary closure

# Common calendar: joint fertility and endogenous stationary closure

Independent mathematical check, September 9. This uses the proposed accounts in common_calendar_model.md. It does not adopt a model into the main note, or establish a policy transition. The main result is a full dated joint-planner fertility increase without imposing \(\beta\ge q\), including an explicit primitive family.

## 1. Reference, planner, and the exact scalar profile

Write \(\theta=\vartheta\). Young utility is
\(\log(c-\chi n)+\alpha\log(h-\kappa n)+\theta\log n\);
old utility is \(\log c_o+\alpha\log h_o+\omega\log e\).
Thus \(\alpha=\gamma\), and both current consumption goods are delivered at the same calendar date.

Take a stationary reference with equal current cohort masses \(N\). Every young household rents at \(r=H_R\); old households own, with mean home \(m>r\) and individual homes at most \(H_O\). Set
\[
x_i=c_i-\chi n_i,\quad n_0=\mathbb E n_i,\quad
\bar x=\mathbb E x_i,\quad d=\beta/q.
\]
The competitive branch implies \(c_{o,i}=d x_i\). Require reference fertility to be privately optimal:
\[
\frac{\theta}{n_i}=\frac{\chi}{x_i}
+\frac{\alpha\kappa}{r-\kappa n_i}. \tag{1}
\]
An equilibrium with arbitrarily fixed fertility does not automatically satisfy (1).

The planner chooses all current \(c,h,n\), gives each currently living household unit utility weight, fixes individual young continuation resources and old estates, and relaxes private finance. Future utility has no direct fertility dependence in this two-age model. Housing is a common divisible stock, there is no ownership taste or fixed sector stock, and tenure may change. Hence the physical owner cap \(H_O\) is available to every planner household. Retaining an individual rental cap instead would be a different problem.

Current resources per young-cohort unit are
\[
\mathcal C=(1+d)\bar x+\chi n_0,\qquad
\mathcal H=r+m<2H_O.
\]
These are physical resources at common delivery dates; no discount factor belongs in their adding-up constraints.

Strict concavity and identical current preferences imply common young fertility \(n\), common young adult consumption and old consumption
\[
X(n)=\frac{\mathcal C-\chi n}{2},
\]
and common young adult space
\[
S(n)=\min\left\{\frac{\mathcal H-\kappa n}{2},\,H_O-\kappa n\right\}.
\]
The young home is \(\kappa n+S(n)\); the old home is \(\mathcal H-\kappa n-S(n)\). The old cap is slack because the young home is weakly larger and total stock is below \(2H_O\).

At the joint optimum every young household also receives a gross home larger than \(r\): its common home exceeds \(\mathcal H/2=(r+m)/2>r\), since \(n_J>0\). This conclusion does not require an individual consumption gain.

The exact fertility derivative of the full optimized resource allocation is
\[
F_J(n)=\frac{\theta}{n}-\frac{\chi}{X(n)}
-\frac{\alpha\kappa}{S(n)}. \tag{2}
\]
It is continuous and strictly decreasing on
\[
0<n<\min\{\mathcal C/\chi,\ \min(\mathcal H,H_O)/\kappa\}.
\]
It tends to \(+\infty\) and \(-\infty\) at the two endpoints. Thus the full joint optimum has a unique interior fertility \(n_J\), and
\[
n_J>n_0\quad\Longleftrightarrow\quad
\frac{\theta}{n_0}>
\frac{2\chi}{(1+d)\bar x}
+\frac{\alpha\kappa}{S_0},\qquad
S_0=\min\left\{\frac{m+r-\kappa n_0}{2},H_O-\kappa n_0\right\}. \tag{3}
\]
This is a global optimum comparison, not a local private response at an assigned bundle.

There is also a quadratic closed form on each physical-cap branch. When the young planner cap is slack, select the feasible root of
\[
\chi\kappa(\theta+2+2\alpha)n^2
-\{(\theta+2)\chi\mathcal H+(\theta+2\alpha)\kappa\mathcal C\}n
+\theta\mathcal C\mathcal H=0.
\]
When it binds, replace the equation by
\[
\chi\kappa(\theta+2+\alpha)n^2
-\{(\theta+2)\chi H_O+(\theta+\alpha)\kappa\mathcal C\}n
+\theta\mathcal C H_O=0.
\]
Feasibility and the applicable cap inequality select the unique root of (2). No utility transformation or descendant welfare term is introduced.

## 2. A cap-valid condition allowing \(\beta<q\)

Let \(s_0=r-\kappa n_0\). The private fertility function defined by (1) is increasing and concave in adult consumption. More generally, in units \(C=x/\chi\), \(R=h/\kappa\), it is
\[
\mathcal N(x,h)=
\frac{R+(\alpha+\theta)C-
\sqrt{[R+(\alpha+\theta)C]^2-4\theta CR}}2.
\]
The square root is a norm: its quadratic form is positive definite, with determinant \(4\alpha\theta>0\). Hence \(\mathcal N\) is jointly concave and homogeneous of degree one. Jensen and (1) give
\[
n_0\le\mathcal N(\bar x,r),\qquad
\frac{\theta}{n_0}\ge\frac{\chi}{\bar x}+\frac{\alpha\kappa}{s_0}. \tag{4}
\]
The first inequality is strict with nondegenerate adult-consumption heterogeneity.

For \(\beta<q\), a sufficient condition for the strict joint fertility gain is therefore
\[
\boxed{\;
\frac{\chi}{\bar x}\frac{q-\beta}{q+\beta}
<\alpha\kappa\left(\frac1{s_0}-\frac1{S_0}\right).
\;} \tag{5}
\]
For a homogeneous reference it is necessary and sufficient. With heterogeneous reference fertility it is sufficient; (3) remains the exact test. When \(\beta\ge q\), (4), \(m>r\), and \(H_O>r\) make (3) automatic.

A simpler sufficient condition is valid even if the planner's young housing cap binds:
\[
\boxed{\;
\frac{\chi r}{\kappa\bar x}\frac{q-\beta}{q+\beta}
<\alpha\frac{m-r}{m+r}\quad(\beta<q).
\;} \tag{6}
\]
Indeed,
\[
\frac1{s_0}-\frac1{S_0}
=\min\left\{\frac1{s_0}-\frac2{m+s_0},
\frac1{s_0}-\frac1{H_O-r+s_0}\right\}
\ge\frac{m-r}{r(m+r)}.
\]
For the second term, use \(s_0<r\) and \(H_O\ge m\); its lower bound
\((H_O-r)/(rH_O)\) is larger than the displayed bound.

Consequently the stronger condition
\[
\boxed{\ \frac{\chi r}{\kappa\bar x}
<\alpha\frac{m-r}{m+r}\ } \tag{7}
\]
works for every \(\beta>0\). This does not say discounting is irrelevant: (6) is sharper and can be much weaker than (7). For illustration, \(m=1.5r,\alpha=1\) makes (7) read \(\chi r/(\kappa\bar x)<.2\). If children occupy \(\kappa n_0=.5r\), this requires their mean goods bill below one tenth of mean adult consumption. This is an illustration, not an empirical calibration.

## 3. Private fertility and a stationary price root

Use gross stage resources
\[
g_i=w_i/q+Y_i,\quad B=(1+q)/q,\quad
A=1+q\tau,\quad D=A-q,\quad
a=\min\{\alpha/D,(\alpha+\omega)/A\},
\]
\[
K=1+\alpha+\omega,\qquad G=1+\beta K.
\]
Conditional on renting at \(r\), saving, and subsequently owning below \(H_O\),
\[
M_i=g_i+BT-DPr/q,\qquad x_i=(M_i-\chi n_i)/G.
\]
The fertility choice maximizes
\(G\log(M_i-\chi n)+\alpha\log(r-\kappa n)+\theta\log n\).
Its unique feasible quadratic root satisfies
\[
\chi\kappa(\theta+G+\alpha)n^2
-\{\kappa M_i(\theta+\alpha)+\chi r(\theta+G)\}n
+\theta M_i r=0. \tag{8}
\]
Call it \(f_\theta(M_i)\). Implicit differentiation gives
\[
0<f_\theta'(M)<1/\chi. \tag{9}
\]
Its derivative is the positive numerator
\(G\chi/(M-\chi n)^2\), divided by
\(\theta/n^2+G\chi^2/(M-\chi n)^2+
\alpha\kappa^2/(r-\kappa n)^2\).
The function is also strictly increasing in \(\theta\).

If \(\nu\) and \(\theta\) are fixed, write \(n_*=\nu^{-1}\),
\(S=\bar g-\chi n_*\), and
\[
\mu_0=\frac{Bq\tau a\beta}{2},\quad
\widetilde G=G-\mu_0,\quad
C_0=D/q-Bq\tau/2=(1-q)(1/q+\tau/2)>0.
\]
Provided \(S>0\) and \(\widetilde G>0\), define
\[
M_i(P)=g_i+\frac{\mu_0 S}{\widetilde G}
-\frac{G}{\widetilde G}C_0Pr. \tag{10}
\]
Then \(\mathbb E f_\theta(M_i(P))=n_*\) has at most one positive feasible root, because every term strictly decreases in \(P\). The mild sufficient condition \(\tau<2/(1-q)\) ensures
\(\mu_0<\beta(\alpha+\omega)\) and \(\widetilde G>1+\beta\).

At a root the eliminated equations reconstruct exactly:
\[
\bar x=\frac{S-C_0Pr}{\widetilde G},\quad
m=\frac{a\beta\bar x}{P},\quad
N=\frac{\bar H}{r+m},\quad
T=\frac{q\tau P(r+m)}2. \tag{11}
\]
Equivalently,
\[
P(N)=\frac{a\beta S}
{G(\bar H/N-r)+(a\beta D/q)r-a\beta Bq\tau\bar H/(2N)}.
\]
The off-root functions in (10) are an elimination device; their tax accounts become jointly consistent at the fertility root.

Here are explicit endpoint tests. Require
\(n_*<\theta r/[\kappa(\alpha+\theta)]\), and define
\[
M_*=\chi n_*+
\frac{G\chi}{\theta/n_*-\alpha\kappa/(r-\kappa n_*)}.
\]
For any \(0<P_L<P_U\) with all \(M_i(P_U)>0\), the primitive inequalities
\[
\min_i M_i(P_L)>M_*,\qquad \max_i M_i(P_U)<M_*
\]
give existence and uniqueness by the intermediate value theorem. Actual household finance, physical caps, and global tenure dominance must also hold at that root or uniformly on the bracket. A scalar fertility root alone is not a general-equilibrium proof.

## 4. A complete nonempty family at fixed cohort normalization

A particularly short completion uses the companion note's normalization \(N=1\), total stock \(H=5/2\), \(r=1\), \(H_O=2\), \(m=3/2\), and
\[
q=\tfrac12,\quad\phi=\tfrac45,\quad
\alpha=\gamma=1,\quad\omega=\tfrac14,\quad
\beta\in[.3,1],\quad\tau\in[.01,.02].
\]
Allow separate heterogeneity in cash and later income subject to
\[
|g_i/\bar g-1|\le.01,\qquad
\chi r/\kappa\le.02\bar g,\qquad
\bar g\ge46w_{\max}. \tag{12}
\]
Define the positive, primitive coefficients
\[
E_r=Dr/q-Bq\tau H/2=1-7\tau/8,\quad
J=Gm+a\beta E_r,\quad k=a\beta E_r/J\in(0,1).
\]
For any proposed mean fertility \(z\in[0,r/\kappa]\),
\[
P(z)=a\beta(\bar g-\chi z)/J,\quad
T(z)=q\tau P(z)H/2,\quad
M_i(z)=g_i-k(\bar g-\chi z).
\]
The map
\[
\mathcal T_\theta(z)=\mathbb E f_\theta(M_i(z))
\]
has derivative strictly below \(k<1\), maps the interval into its interior, and is continuous. It therefore has a unique fixed point \(n_0(\theta)\).

The required goods margin is explicit: \(J/G<2\), hence \(k<1/4\), so
\(M_i(z)\ge(.99-.25)\bar g>\chi r/\kappa\).
As \(\theta\) rises from zero to infinity, the unique fixed point rises continuously and strictly from zero to \(r/\kappa\). Thus, for every specified \(\nu>\kappa/r\), exactly one positive \(\theta\) gives \(n_0=1/\nu\). Selecting \(\theta\) this way is a one-dimensional preference calibration, not a closed-form solution for \(\theta\). At that calibrated value, \(P=a\beta(\bar g-\chi/\nu)/J\) is explicit and every individual fertility follows (8). Alternatively, choosing \(\theta\) first defines the compatible demographic coefficient \(\nu=1/n_0(\theta)\).

The global tenure checks survive endogenous fertility. At the actual candidate let \(S=\bar g-\chi n_0\in[.98\bar g,\bar g]\). For every individual deviation \(n\in(0,r/\kappa)\),
\[
g_i-\chi n\in[.97\bar g,1.01\bar g],\qquad
\left|\frac{g_i-\chi n}{S}-1\right|\le\frac3{98}<.05,\quad
S>45w_{\max}.
\]
Holding the actual price and rebate fixed,
\[
h_{o,i}(n)=m+\frac JG\left(\frac{g_i-\chi n}{S}-1\right).
\]
Consequently the companion note's five-percent-support, saving, housing-cap and fixed-\(n\) RO-versus-RR value certificates hold **uniformly in the deviating fertility**, not merely at the proposed optimum. Any feasible young owner can also be replicated by a young renter at the same \(c,h,n,Z\). This rules out joint deviations in tenure and fertility and completes the competitive branch, using the companion note's verified household inequalities.

Finally, this family meets the discount-independent joint-planner condition (7). Since \(G\le13/4\), \(a\beta<5/4\), and \(E_r<1\),
\[
J<49/8,\qquad
\bar x=\frac{mS}{J}>\frac6{25}\bar g,\qquad
\frac{\chi r}{\kappa\bar x}<\frac1{12}<\frac15
=\frac{m-r}{m+r}.
\]
The full joint planner therefore strictly increases average fertility throughout this family, including \(\beta<q\). The physical cap \(H_O\) remains in the proof.

## Scope and verification

The result establishes a nonempty analytical family of stationary competitive allocations and a full dated joint-planner fertility gain, together with a more general conditional scalar price method. The fixed-\(N\) construction uses a unique scalar calibration of \(\theta\) to satisfy an externally specified replacement coefficient. It does not prove global equilibrium uniqueness across tenure regimes, signed tax effects, a transition initiated by fertility tastes, or that private fertility necessarily rises after the separate fixed-fertility planner allocation. The planner remains a redistributive utilitarian benchmark over the currently living; no Pareto or descendant-welfare claim is made.

All derivations were checked algebraically against the common-calendar accounts. No numerical model solve, root search, simulation, or main-file edit was used.


---

# Part III. Same-model perfect-foresight transition

# Common-calendar model: stationary endpoints and a local perfect-foresight transition

September 9, 2026. Bounded mathematical extension of `common_calendar_model.md` and `common_calendar_fertility.md`. This uses the **same model and accounts**. It does not adopt the candidate into the main theory or claim global equilibrium uniqueness. Sections 1–5 first derive the transition for homogeneous households; Section 7 proves that the same aggregate recursion, stability, and shock signs extend to the actual heterogeneous primitive family.

## 1. Scope and notation

On the maintained branch, every young household rents at the physical rental cap \(r\), saves strictly, and becomes an uncapped old owner. The old financial-estate floor binds. Let

\[
b=\gamma+\omega,\quad K=1+b,\quad G=1+\beta K,
\quad A(\tau)=1+q\tau,\quad a(\tau)=b/A(\tau).
\]

The homogeneous young household has gross stage-end resources
\(B=w/q+Y^{\rm income}\). Its adult consumption is \(x=c-\chi n\). Cohort masses satisfy

\[
Y_{t+1}=\nu n_tY_t,\qquad O_{t+1}=Y_t.
\tag{T1}
\]

The stock \(\bar H\) is fixed. Prices and rebates are dated, and the rebate is

\[
T_t=\frac{q\tau_tP_t\bar H}{Y_t+O_t}.
\tag{T2}
\]

All rent and asset-price terms below use perfect foresight after the shock. The mortgage share remains \(\phi=.8\) in the explicit family, while \(q=.5\); no relation between \(\phi\), \(q\), and \(\beta\) is imposed.

At dated prices, a floor-constrained old owner solves

\[
q c_o+A(\tau_t)P_th_o=Z_t,
\qquad e_t=P_{t+1}h_o,
\]

so

\[
c_{o,t}=\frac{Z_t}{qK},\qquad
h_{o,t}=\frac{a(\tau_t)Z_t}{KP_t}.
\tag{T3}
\]

Its value has resource-dependent part \(K\log Z_t\). The future sale price appears in its additive \(\omega\log P_{t+1}\) term. That term is independent of saving and fertility conditional on this tenure branch; it has not been set equal to the current price. The floor binds when

\[
\omega\big[A(\tau_t)P_t-qP_{t+1}\big]<q\gamma P_{t+1},
\]

and end-of-stage rent is \([A(\tau_t)P_t/q-P_{t+1}]h\). Both inequalities are strict at the reference family and therefore survive sufficiently small local paths.

## 2. Exact transition recursion, including unexpected rebates

The state at the beginning of a date is \((Y_t,O_t,s_t)\), where \(s_t\) is the **inherited financial saving per old household**. The actual old resources are

\[
Z_t=s_t+T_t.
\]

Current market clearing and (T2) give an explicit price:

\[
\boxed{\displaystyle
P_t=\frac{a(\tau_t)O_ts_t}
{K(\bar H-rY_t)-a(\tau_t)O_tq\tau_t\bar H/(Y_t+O_t)}.}
\tag{T4}
\]

The denominator and \(s_t\) must be positive on this saving branch. This is the correct equation at an unexpected tax reform: saving is inherited, while the actual rebate changes. In general one must not set \(Z_t=\beta Kx_{t-1}\) at such an impact date. That equality holds only if the rebate equals what the previous young household anticipated when saving.

The private fertility condition at the binding young rental cap is

\[
\frac{\theta_t}{n_t}=\frac{\chi}{x_t}
+\frac{\alpha\kappa}{r-\kappa n_t}.
\tag{T5}
\]

Equivalently, write

\[
X(n;\theta)=
\frac{\chi n(r-\kappa n)}{\theta r-\kappa(\alpha+\theta)n},
\quad
0<n<\frac{\theta r}{\kappa(\alpha+\theta)}.
\tag{T6}
\]

This inverse schedule is strictly increasing in \(n\), decreasing in \(\theta\), and has elasticity \(nX_n/X>1\). It is the same household first-order condition, not an imposed demographic rule.

Optimal positive saving gives

\[
s_{t+1}+T_{t+1}=\beta Kx_t.
\tag{T7}
\]

The full young resource equation therefore is

\[
Gx_t+\chi n_t
=B+T_t/q+T_{t+1}
-\left[\frac{A(\tau_t)P_t}{q}-P_{t+1}\right]r.
\tag{T8}
\]

The next house price and rebate are determined by the current choices:

\[
P_{t+1}=\frac{a(\tau_{t+1})\beta Y_tx_t}
{\bar H-r\nu n_tY_t},
\quad
T_{t+1}=\frac{q\tau_{t+1}P_{t+1}\bar H}
{Y_t(1+\nu n_t)}.
\tag{T9}
\]

Define \(k_t=rY_t/\bar H\), \(v=\nu n\),

\[
\mathcal M(n;Y_t,\tau_{t+1})
=\frac{\beta a(\tau_{t+1})
[k_t+q\tau_{t+1}/(1+\nu n)]}
{1-k_t\nu n},
\]
\[
L_t=B+T_t/q-A(\tau_t)P_tr/q.
\]

Then a single scalar equation determines current fertility:

\[
\boxed{\displaystyle
[G-\mathcal M(n_t;Y_t,\tau_{t+1})]X(n_t;\theta_t)
+\chi n_t=L_t.}
\tag{T10}
\]

Given its admissible root, reconstruct \(x_t\), \(Y_{t+1}\), \(O_{t+1}\), the prices and rebates in (T9), and \(s_{t+1}=\beta Kx_t-T_{t+1}\). This is an exact finite-dimensional forward equilibrium recursion on the stated branch. Expected resale prices are present in (T8)–(T10). The local result below establishes a unique nearby root and convergence; it does not say (T10) has a unique admissible root globally.

## 3. Closed stationary endpoints with positive property taxes

At a positive stationary population let \(n_*=1/\nu\) and require

\[
n_*<\frac{\theta r}{\kappa(\alpha+\theta)}.
\]

For a homogeneous household, (T6) pins

\[
x_*=\frac{\chi n_*(r-\kappa n_*)}
{\theta r-\kappa(\alpha+\theta)n_*}.
\tag{E1}
\]

Put

\[
c_\tau=\frac{1+q}{2},\quad
C_0=(1-q)(1/q+\tau/2)>0,
\quad \widetilde G=G-c_\tau\tau a(\tau)\beta.
\]

The exact solution is

\[
\boxed{\displaystyle
P_* =\frac{B-\chi n_*-\widetilde Gx_*}{C_0r},\qquad
N_* =\frac{\bar H}{r+a(\tau)\beta x_*/P_*}.}
\tag{E2}
\]

Rebates are \(T_*=q\tau P_*[r+a\beta x_*/P_*]/2\). Thus (E2) solves the fiscal and housing equations jointly. A positive price, positive saving, and the cap/tenure inequalities still define the branch. In the explicit family \(\widetilde G>0\).

### A fertility-taste decline lowers the destination population

For fixed \(\tau\), \(x_{*,\theta}<0\). With \(\widetilde G>0\), a fall in \(\theta\) raises adult consumption at replacement fertility, lowers \(P_*\), and raises the old home. It therefore lowers \(N_*\). This is a change in the population **level**, since fertility at both positive endpoints is still \(1/\nu\).

For completeness, let \(S=B-\chi n_*>0\), \(A=1+q\tau\), and \(b=\gamma+\omega\). The mean old home is

\[
m_*=
\frac{C_0r\beta b x_*}
{AS-AGx_*+c_\tau\tau\beta b x_*},
\]

whose derivative in \(x_*\) is strictly positive:

\[
\frac{\partial m_*}{\partial x_*}
=\frac{C_0r\beta b AS}
{[AS-AGx_*+c_\tau\tau\beta b x_*]^2}>0.
\]

### A permanent tax increase raises the destination population

Hold \(\theta\) fixed, so \(x_*\) and \(n_*\) are fixed. Define

\[
U=B-\chi n_*-Gx_*,\qquad Q=\beta b x_*>0,
\qquad d=1-q.
\]

Then the population has a fractional-linear closed form:

\[
\boxed{\displaystyle
N_*(\tau)=\frac{\bar H}{r}
\frac{(1+q\tau)U+c_\tau\tau Q}
{(1+q\tau)U+(d/q+\tau)Q}.}
\tag{E3}
\]

Differentiating gives

\[
\boxed{\displaystyle
N_{*,\tau}=\frac{\bar H}{r}
\frac{dQ[U+(1+q)Q/q]}
{2[(1+q\tau)U+(d/q+\tau)Q]^2}>0.}
\tag{E4}
\]

The sign follows just from a positive stationary price and \(\tau\ge0\) within this binding-floor branch. Positive price implies

\[
U> -\frac{c_\tau\tau Q}{1+q\tau},
\]

and \(c_\tau\tau/(1+q\tau)<(1+q)/q\). Hence the numerator in (E4) is positive even if \(U<0\). The old mean home falls. The young home remains \(r\), so the young aggregate housing share rises because more young households occupy the stock. This tax result is not the dated planner's per-young housing enlargement.

The purchase-price derivative has its own exact formula:

\[
\boxed{\displaystyle
\frac{P_{*,\tau}}{P_*}
=\frac{c_\tau m_*}{(1+q\tau)C_0r}
-\frac{1-q}{2C_0}.}
\tag{E6}
\]

Its sign is that of \((1+q)m_*-(1-q)(1+q\tau)r\). It is strictly positive in the explicit family \(q=.5,m_*=1.5r,\tau\in[.01,.02]\). Thus the higher stationary population comes with a **higher stationary purchase price** in this family. The old mean home falls because \(a(\tau)/P_*\) falls. There is no universal claim that a property tax lowers the stationary asset price after population adjusts.

### Heterogeneous stationary endpoints

The same endpoint population signs extend to the heterogeneous branch. Write \(g_i=w_i/q+Y_i^{\rm income}\). Every stationary household faces resources \(M_i=g_i+c\), where \(c\) is the common net rebate/rent shift, and has fertility \(f_\theta(M_i)\) from the companion note's strictly concave household problem. Replacement pins the unique shift by

\[
\mathbb E f_\theta(g_i+c_*)=1/\nu.
\tag{E5}
\]

At fixed \(\theta\), (E5) is independent of \(\tau\); therefore every \(n_i\), every \(x_i\), and \(\bar x\) remain unchanged as the tax varies within the branch. Equations (E2)–(E4) apply with \(B=\bar g\) and \(x_*=\bar x\).

Since \(f_M>0\) and \(f_\theta>0\), implicit differentiation of (E5) gives \(c_{*,\theta}<0\). Since
\(G\bar x=\bar g+c_*-\chi/\nu\), mean adult consumption decreases with the taste for children. Thus the taste endpoint sign also survives heterogeneity. These stationary results alone do not establish a transition. Section 7 supplies the additional aggregation and elasticity proof; it never replaces average fertility by representative-household fertility at average resources.

## 4. Local existence and stability in the explicit family

Use the companion branch's reference values

\[
q=1/2,\quad \phi=4/5,\quad r=1,\quad
m_*=3/2,\quad \bar H/N_*=5/2,
\]
\[
\alpha=\gamma=1,\quad \omega=1/4,\quad K=9/4,
\quad \beta\in[.3,1],\quad \tau\in[.01,.02].
\tag{S1}
\]

The representative household's resources and fertility satisfy the common-calendar branch's strict primitive finance/tenure conditions. The homogeneous specification is used first to derive the transition; Section 7 justifies exactly the same calculation for heterogeneous households. No numerical root or simulated equilibrium is used in the following stability argument.

### The scalar current-choice equation is regular

At a stationary point \(k=rN_*/\bar H=2/5\), \(\nu n_*=1\). Write \(\mathcal M_v\) for the derivative of \(\mathcal M\) with respect to \(v=\nu n\) at fixed current young mass. Direct differentiation gives

\[
\mathcal M_*=
\frac{\beta a(k+q\tau/2)}{1-k},\qquad
\mathcal M_v=
\frac{\beta a[k^2+q\tau(3k-1)/4]}{(1-k)^2}.
\]

Using \(a\le5/4\) and \(q\tau\le1/100\),

\[
\mathcal M_*\le\frac{27}{32}\beta,\qquad
0<\mathcal M_v\le\frac{107}{192}\beta,
\]

so

\[
G-\mathcal M_*-\mathcal M_v
\ge1+\frac{163}{192}\beta>0.
\tag{S2}
\]

Because \(nX_n/X>1\), the derivative in \(n\) of the left side of (T10) is strictly positive. The implicit-function theorem gives a unique smooth nearby current root. Current price (T4) also has a strictly positive denominator at (S1), since
\(Km_*-a q\tau(r+m_*)/2>27/8-1/64>0\).

### Complete linearized recursion

Hold primitives constant after a shock. Let \(\ell_t=\delta Y_t/N_*\), \(z_t=\delta x_t/x_*\), and let

\[
e=\frac{x_*}{n_*}\frac{\partial n}{\partial x}\in(0,1),
\qquad h=\frac{r}{m_*}=2/3.
\]

For an anticipated date, proportional price changes satisfy

\[
p_t=\ell_{t-1}+z_{t-1}+h\ell_t,
\qquad \ell_{t+1}=\ell_t+ez_t.
\tag{S3}
\]

This is also a valid local state representation using \(Z_t/(\beta K)\) as the inherited resource coordinate. At an unexpected tax date that coordinate includes the surprise rebate; the truly inherited coordinate remains \(s_t\), as in (T4).

Define positive coefficients

\[
w=\frac{T_*}{x_*}=\frac{q\tau a\beta(1+h)}2,
\quad u=\frac{A P_*r-T_*}{qx_*}
=\frac{a\beta hA-w}{q},
\quad v=\frac{P_*r+T_*}{x_*}=a\beta h+w,
\]
\[
\Lambda=-uh+v(1+h)-\frac{w}{2q}-w,
\]
\[
\mathcal D=G+\frac{\chi n_*e}{x_*}
-v(1+he)+\frac{we}{2}.
\tag{S4}
\]

Linearizing the **full** resource equation (T8) gives

\[
\mathcal D z_t
=\Lambda\ell_t-(u+w/(2q))\ell_{t-1}-u z_{t-1}.
\tag{S5}
\]

Together, (S3) and (S5) give the three-dimensional Jacobian exactly. It has one zero root. Eliminating \(z_{t-1}=(\ell_t-\ell_{t-1})/e\) gives the two remaining roots from

\[
\ell_{t+1}=A_2\ell_t+B_2\ell_{t-1},
\]
\[
A_2=1+\frac{e\Lambda-u}{\mathcal D},\qquad
B_2=\frac{u-e(u+w/(2q))}{\mathcal D}.
\tag{S6}
\]

### Both nonzero roots are strictly inside the unit disk

In (S1), \(vh-w/2=\mathcal M_v\), so (S2) implies

\[
\mathcal D\ge1+163\beta/192.
\]

Also

\[
0<u\le101\beta/60,\qquad 0<w\le\beta/96,
\]

and hence, for the full \(\beta\le1\) interval,

\[
\mathcal D-u\ge1-267\beta/320\ge53/320>0.
\tag{S7}
\]

At these values,

\[
\Lambda=\frac{2a\beta}{9}(1-2\tau)+w>0.
\]

For the quadratic \(z^2-A_2z-B_2\), the three strict unit-disk conditions follow directly:

\[
1-A_2-B_2
=\frac{e}{\mathcal D}
\left[u+\frac{w}{2q}-\Lambda\right]
=\frac{e a\beta h(1+h)C_0}{\mathcal D}>0,
\]
\[
1+A_2-B_2
=\frac{2(\mathcal D-u)+e(\Lambda+u+w/(2q))}{\mathcal D}>0,
\]
\[
1+B_2>1-\frac{w}{2q\mathcal D}>0.
\tag{S8}
\]

These are the elementary quadratic stability inequalities, not an assumed eigenvalue pattern. Together with the zero root, they make the entire smooth forward recursion locally asymptotically stable. A norm with linear contraction factor below one exists; continuity gives a contraction in a sufficiently small neighborhood. Thus small permanent shocks have unique nearby forward paths converging to their corresponding nearby stationary equilibria. The strict private-choice margins keep those paths on the verified household branch, so the forward paths satisfy perfect foresight. No extra jump price is chosen or suppressed.

Negative or complex adjustment modes are not excluded. The result is local convergence, not monotone convergence and not global uniqueness across tenure regimes.

## 5. Shock and policy signs along the same path

### Taste shock

At the old stationary state, an unexpected change in \(\theta\) does not change current \(Y,O,s,\tau\). Equation (T4) therefore leaves the current price and rebate unchanged: \(L_t\) is fixed. In (T10), \(X_\theta<0\), \(G-\mathcal M_*>0\), and the fertility-root derivative is positive. Therefore \(\partial n_t/\partial\theta_t>0\). A small permanent taste decline lowers fertility on impact, lowers the next young mass, and starts the locally convergent transition to the lower endpoint in (E2). Both old and new positive endpoints have replacement fertility.

### Unexpected permanent tax at the old state or nearby baseline states

At an unexpected intervention, current inherited savings are fixed. At an equal-cohort stationary inherited state with old mean home \(m>r\), (T4) becomes

\[
P_t(\tau)=
\frac{b s_t}{Km+q\tau[Km-b(r+m)/2]}.
\]

This includes the surprise old rebate. Define its denominator as \(D_o(\tau)\). Then

\[
L_t(\tau)=B-\frac{P_t(\tau)}q
\left[r+\frac{q\tau(r-m)}2\right],
\]
\[
L_{t,\tau}=
\frac{b s_t(r+m)(Km-br)}{2D_o(\tau)^2}>0,
\tag{I1}
\]

because \(m>r\) and \(K=1+b\). The permanent policy also changes the anticipated next-date coefficient. At fixed current choices,

\[
\mathcal M_\tau
=\frac{\beta b q[1/(1+\nu n)-k_t]}
{(1+q\tau)^2(1-k_t\nu n)}.
\tag{I2}
\]

At the stationary reference \(\nu n=1\) and \(k_t=r/(r+m)<1/2\), this is strictly positive. Differentiating (T10) therefore gives

\[
\frac{\partial n_t}{\partial\tau}
=\frac{L_{t,\tau}+\mathcal M_\tau x_*}
{\partial_n\{(G-\mathcal M)X+\chi n\}}>0.
\tag{I3}
\]

The two channels here are the current cash-resource effect, including the actual rebate to the inherited old, and the anticipated next-date resale/rent/rebate effect. Both are retained.

The current purchase price falls on impact:

\[
\frac{P_{t,\tau}}{P_t}
=-\frac{q[Km-b(r+m)/2]}{D_o(\tau)}<0,
\]

since \(m>r\). Together with (E6), this gives a genuine sign reversal between impact and destination for the explicit family: the unexpected tax lowers the purchase price while inherited cohorts and savings are fixed, but the policy's larger terminal population has a higher purchase price. It does not imply a monotone intervening price path. The actual impact rebate rises, since \(\partial_\tau(\tau P_t)=P_tKm/D_o(\tau)>0\).

These impact inequalities and the stability inequalities are strict. They therefore continue at sufficiently nearby inherited states and nearby endpoints. Start with a sufficiently small permanent taste decline, follow its local perfect-foresight baseline, and at any finite date introduce a sufficiently small unexpected permanent tax increase from exactly the same inherited \((Y,O,s)\). Both paths exist locally and converge. Policy fertility is higher on impact relative to baseline, and its terminal population is higher by (E4). Later fertility gaps need not keep one sign.

This is the requested same-model transition statement, initially derived on the homogeneous local branch. Section 7 extends it to the heterogeneous primitive family under its explicit liquidity/tenure restrictions. It is not a population forecast or a global dynamics theorem.

**Mechanism qualification:** all young households continue renting at \(h_t=r\) on this local branch. Credit limits restrict their owner alternative, and policy does not relax that origination restriction or produce realized young home upgrades. The fertility response operates through goods, saving, and dated rebates and rental costs. The stationary young aggregate housing share rises with the number of young households; their individual home remains fixed. The full dated planner's larger young homes and higher fertility are a separate allocation comparison, not an implementation description of this tax experiment.

## 6. Population meaning and verification

The model counts adult households. At a positive steady state \(Y=O=N\), total adult-household population is \(2N\). A resident-person statement needs the additional fixed headship and dependent-child counting convention. For common inherited young mass,

\[
\frac{Y_T^P}{Y_T^B}
=\prod_{t=t_p}^{T-1}\frac{n_t^P}{n_t^B}.
\]

Both terminal fertility rates equal replacement; the higher policy endpoint records the cumulative transition fertility gap. It does not require fertility permanently above replacement.

Verification: the stationary price and fiscal elimination, the fractional-linear tax derivative, the surprise-rebate price equation, the scalar transition equation, and the three-dimensional linearization were derived directly from the dated accounts. Independent symbolic differentiation verified the tax derivative and surprise-rebate resource derivative; symbolic expansion verified the full resource linearization, the Jacobian's zero-times-quadratic characteristic factorization, and the first stability identity. Exact rational checks verified all displayed stability margins. No numerical equilibrium, characteristic-root search, simulation, figure, or model-code mutation is used.

## 7. Exact heterogeneous extension for the actual primitive family

This extension retains the companion note's nondegenerate endowment distribution:

\[
|g_i/\bar g-1|\le1/100,\quad
\chi r/\kappa\le\bar g/50,\quad
\bar g\ge46w_{\max},
\]

and all parameters in (S1), including \(\beta\in[.3,1]\). The entering type distribution is the same exogenous distribution \(F\) in every cohort, as maintained by the model. It is not selected by parents' fertility, and estates do not endogenously become entrants' cash. Endogenizing that intergenerational selection or inheritance would change the aggregation result.

### Exact aggregation uses a common resource shift, not a representative household

At any nearby dated prices, define the common resource shift

\[
C_t=T_t/q+T_{t+1}
-[A(\tau_t)P_t/q-P_{t+1}]r.
\]

Each young type satisfies

\[
M_{it}=g_i+C_t=Gx_{it}+\chi n_{it},\qquad
n_{it}=f_{\theta_t}(M_{it}),
\tag{H1}
\]

where \(f_\theta\) is the companion note's unique feasible quadratic household solution. Consequently the two mean schedules

\[
\bar n(C,\theta)=\mathbb E f_\theta(g_i+C),\qquad
\bar x(C,\theta)=
\frac{\bar g+C-\chi\bar n(C,\theta)}G
\tag{H2}
\]

are exact. Since \(0<f_M<1/\chi\), \(\bar x_C>0\). Locally invert \(C\mapsto\bar x\) to define the aggregate fertility schedule

\[
\bar n=\mathcal N(\bar x,\theta).
\]

This schedule generally differs from the individual fertility function evaluated at \(\bar x\). Its inverse, denoted \(\mathcal X(\bar n,\theta)\), replaces (T6) in the exact mean version of (T10).

At fixed \(C\), \(\bar n_\theta>0\). At fixed \(\bar x\), (H2) gives the useful exact derivative

\[
\mathcal N_\theta
=\frac{\bar n_\theta|_C}{1-\chi\bar n_C}>0.
\tag{H3}
\]

Hence \(\mathcal X_\theta<0\), just as in the individual inverse schedule.

### The aggregate elasticity is strictly below one, with a large explicit margin

At a private optimum let \(d_i=\partial n_i/\partial x_i>0\). Differentiating the private fertility condition gives

\[
d_i=
\frac{\chi/x_i^2}
{\theta/n_i^2+\alpha\kappa^2/(r-\kappa n_i)^2}.
\]

Since \(\theta>\alpha\kappa n_i/r\),

\[
d_i\le\frac{\chi r}{\alpha\kappa}\frac{n_i}{x_i^2}.
\tag{H4}
\]

Differentiating (H1) with respect to the common shift gives

\[
\rho_i\equiv\frac{\partial x_i}{\partial C}
=\frac1{G+\chi d_i},\qquad
\frac{\partial n_i}{\partial C}=\rho_i d_i.
\]

Thus the aggregate elasticity used in (S3)–(S8) is exactly

\[
e=\frac{\bar x}{\bar n}
\frac{\mathbb E(\rho_i d_i)}{\mathbb E\rho_i}.
\tag{H5}
\]

To bound it entirely by primitives, use the stationary construction's coefficient \(k<1/4\) and \(G\le13/4\). Its individual adult consumption obeys

\[
x_i=\frac{g_i-k(\bar g-\chi\bar n)-\chi n_i}{G}
\ge\frac{(.99-.25-.02)\bar g}{G}
=\frac{18\bar g}{25G},
\]

and \(\bar x\le101\bar g/(100G)\). Let \(x_{\min}=18\bar g/(25G)\). Because \(n_i<r/\kappa\), (H4), with \(\alpha=1\), implies

\[
\chi\max_i d_i
\le\left(\frac{\chi r}{\kappa x_{\min}}\right)^2
\le\left(\frac{13}{144}\right)^2<\frac1{100}.
\tag{H6}
\]

The weights in (H5) lie between \(1/(G+1/100)\) and \(1/G\). Therefore

\[
\begin{split}
e
&\le\frac{\chi r}{\alpha\kappa}
\frac{\bar x}{x_{\min}^2}\frac{G+1/100}{G}\\
&\le\frac{(1/50)(101/100)(13/4+1/100)}{(18/25)^2}
=\frac{16463}{129600}<\frac{13}{100}<1.
\end{split}
\tag{H7}
\]

Every inequality is uniform over the permitted heterogeneous distribution. It follows that \(\bar n\mathcal X_{\bar n}/\bar x=1/e>1\). The strict margin survives in a sufficiently small neighborhood of the reference, and (H3) supplies the required taste derivative there.

### Why only means enter the aggregate dynamic state on this branch

Let \(s_{it}\) be individual inherited saving. At current prices, every old owner's demand is linear in resources:

\[
\bar h_{o,t}=\frac{a(\tau_t)}{KP_t}
(\bar s_t+T_t).
\]

Hence (T4) holds exactly with \(s_t=\bar s_t\). Every young type has
\(s_{i,t+1}+T_{t+1}=\beta Kx_{it}\), so

\[
\bar s_{t+1}=\beta K\bar x_t-T_{t+1}.
\]

The exact **aggregate** state is therefore \((Y_t,O_t,\bar s_t)\), and the aggregate recursion is (T1)–(T10) with

\[
B\mapsto\bar g,\qquad x\mapsto\bar x,\qquad
n\mapsto\bar n,\qquad X\mapsto\mathcal X.
\]

Individual old profiles are still used to recover allocations and check the private regime. They are not literally identical across states sharing a mean. Their mean-zero differences have no aggregate demand effect on this linear branch and disappear when that old cohort exits. Initial profiles must be the constructed stationary profile or uniformly close to it; proximity of the mean alone is insufficient. Every new profile is the smooth type-by-type map from the common price/resource shifts, on the fixed compact endowment support. The nearby price path therefore preserves the companion note's uniform saving, size, estate-floor, and tenure-dominance margins. Thus no unverified infeasible old type is hidden inside aggregation.

### Consequence for the same-model transition theorem

The stability calculation in Section 4 uses only the exact aggregate resource identity, linear old demand, and \(0<e<1\). Equations (H3) and (H7) verify precisely those missing properties. All inequalities (S2)–(S8) therefore apply unchanged to the heterogeneous family. The three-dimensional aggregate recursion has a unique locally stable path, and individual profiles are recovered from (H1) and the old linear demands while retaining the strict branch margins.

The taste-shock and tax-impact proofs in Section 5 likewise apply to **mean fertility**. A small permanent taste decline lowers mean fertility on impact and initiates convergence to the smaller stationary population. A small unexpected permanent tax increase at any finite date along that nearby baseline raises mean fertility on impact and raises the terminal population. It starts from the same full inherited population and savings profile; in particular, the old receive their actual surprise rebate rather than a retroactively changed saving choice. The impact price fall and terminal price rise also survive.

This removes the representative-household limitation for the stated one-percent heterogeneous endowment family. The result remains local to the verified regime; it does not establish all-date fertility ordering, monotone convergence, global equilibrium uniqueness, or a transition under endogenous entrant-type selection. Young housing remains \(r\) and young tenure remains rental throughout this branch.

Additional verification: exact rational arithmetic checked (H6)–(H7), and symbolic differentiation checked the individual elasticity formula and the fixed-mean-consumption taste derivative (H3). The heterogeneous argument is additional to the separately audited homogeneous theorem; its aggregate closure and profile requirements are stated explicitly above.


---

# Part IV. Alternative repayment calendars

# Mortgage repayment with boundary title settlement

September 9, 2026. Bounded accounting appendix to
`common_calendar_model.md`. No original source or earlier scratch memo is
edited by this task.

## Verdict

Yes. The extra requirement of nonnegative financial wealth **before**
liquidating the young home can be replaced by nonnegative wealth **after**
its sale and full mortgage repayment. This gives an exact feasible financing
problem, without debt forgiveness or a new old-age loan. Every expanded
young-owner plan remains replicable by a renter when its origination cap
lies below the rental size ceiling. Consequently the verified all-young-
renter equilibrium and its global tenure proof survive unchanged.

The required timing statement is explicit: stage-end young consumption,
property tax, mortgage repayment, and sale of the young home clear together.
Sale proceeds may finance consumption as well as repay debt. A requirement
that consumption be paid before this sale would instead retain an additional
liquidity inequality. The common calendar permits simultaneous settlement;
it does not make the payment ordering irrelevant without stating it.

## 1. Exact nonstationary accounts

Let \(W_t=w_i+T_t\) include all cash already received at the beginning of
the young stage. The household receives genuinely later working income
\(Y_{i,t}\) at its end. The initial mortgage and bond account is unchanged:
\[
k+P_t h=W_t+d,
\qquad k\ge0,\quad 0\le d\le\phi P_t h.
\tag{1}
\]
The home provides current housing services before being sold at the age
boundary for \(P_{t+1}h\). At that boundary, the household receives its labor
income and matured bonds, fully repays the mortgage including interest,
pays current consumption and property tax, and retains wealth
\[
\boxed{\widetilde a
=Y_{i,t}+\frac{k}{q}+P_{t+1}h
-c-\tau_tP_t h-\frac dq\ge0.}
\tag{2}
\]
There is no mortgage balance after this settlement. Define the accounting
balance before title liquidation by
\[
a'_O\equiv\widetilde a-P_{t+1}h.
\tag{3}
\]
It can be negative; that number is not an unsecured loan carried into old
age. It records the remaining mortgage and other contemporaneous payments
net of nonhousing receipts before the house proceeds are included.

Eliminating the initial bond and loan positions gives
\[
qc+qa'_O+(1+q\tau_t)P_t h=W_t+qY_{i,t},
\qquad a'_O+P_{t+1}h\ge0,
\tag{4}
\]
together with the original deposit restriction
\[
\boxed{(1-\phi)P_t h\le W_t.}
\tag{5}
\]
Equivalently, put
\(u_t=(1+q\tau_t)P_t-qP_{t+1}>0\). The complete reduced constraint is
\[
\boxed{qc+q\widetilde a+u_t h=W_t+qY_{i,t},
\qquad\widetilde a\ge0,}
\tag{6}
\]
with (5), the physical owner size cap, and the utility domains.

These conditions are sufficient as well as necessary. Given a reduced
feasible allocation, choose
\[
d=\max\{0,P_t h-W_t\},\qquad k=W_t+d-P_t h.
\tag{7}
\]
The deposit condition ensures \(d\le\phi P_t h\), while \(k\ge0\).
Equation (6) then gives (2) exactly, so all the mortgage principal, interest,
tax, and consumption bills can be paid from the contemporaneous receipts.
There is no missing interim credit line. At a binding deposit cap,
\(d=\phi P_t h\) and \(k=0\), as expected.

## 2. Old-age entry and exact renter replication

After settlement the old household receives its new stage's upfront rebate
and starts with
\[
Z_{t+1}=\widetilde a+T_{t+1}
=a'_O+P_{t+1}h+T_{t+1}.
\tag{8}
\]
It can rent or buy within the original size menus, using the unchanged
old-age budget and estate floor. Sale and repurchase have no real resource
cost in this frictionless title market. The household can repurchase its
previous home if that home is affordable; the old no-borrowing restriction
still governs its selected bundle. The argument does not guarantee that
retaining the previous size is feasible.

The young renter budget is
\[
qc+qa'_R+u_t h=W_t+qY_{i,t},\qquad a'_R\ge0.
\tag{9}
\]
For every feasible young-owner plan with \(h\le H_R\), set
\[
\boxed{a'_R=\widetilde a=a'_O+P_{t+1}h.}
\tag{10}
\]
Equations (6) and (9) coincide; current consumption, housing, fertility,
and future old resources are identical. The future tenure menu is shared.
Thus, with no current ownership taste, a renter replicates every owner
plan whenever
\[
\frac{W_t}{(1-\phi)P_t}<H_R.
\tag{11}
\]
The same identity holds at nonstationary prices. It uses the future sale
price \(P_{t+1}\), not \(P_t\), in (3), (4), and (10).

## 3. Consequences and boundary cases

The stationary family in `common_calendar_model.md` satisfies (11) strictly.
Its chosen households all rent, save strictly, and later select ownership.
Its renter optimum, old value functions, price/rebate formulas, housing
clearing, and exact RO-versus-RR value comparison are unchanged. Enlarging
the young-owner menu cannot defeat that optimum, because every additional
owner plan is covered by (10). The uniform fertility-deviation comparison
is also unchanged: replication keeps the same \(n\).

The same conclusion applies to a previously verified deterministic path
that remains in this certified regime: its reduced equations and global
tenure argument are unaffected. This appendix proves no new transition and
makes no claim about paths on which the owner origination cap rises above
\(H_R\), or where another maintained strict regime condition fails.

The common-date planner's consolidated account is also unchanged. Fixing
future \(Z\) and the future rebate fixes \(\widetilde a\); (6) then gives
\(q\Delta c+u_t\Delta h\) as the required current transfer. In the
pre-liquidation notation, \(\Delta a'_O=-P_{t+1}\Delta h\). The rental
intermediary and old estate ledgers still have to be included as in the
common-calendar memo. Changing this private repayment restriction supplies
no new aggregate goods or fiscal resources.

The relevant corners are:

- **Zero remaining wealth:** \(\widetilde a=0\) is feasible settlement.
  If \(T_{t+1}=0\) as well, positive old consumption, housing, and estate
  cannot all be financed; the old utility domain then excludes that plan.
- **A price decline or mortgage larger than the sale receipt:** other
  stage-end income may cover the shortfall. If total receipts do not cover
  every payment, (2) fails; no default is being allowed.
- **The next rebate:** (2) imposes solvency before crediting
  \(T_{t+1}\). This preserves the renter's same pre-rebate saving condition.
  If that next rebate is instead spendable in the very same debt settlement,
  both renter and owner solvency permissions must be changed symmetrically;
  imposing only \(Z_{t+1}\ge0\) on owners would invalidate (10)'s renter
  feasibility at negative \(\widetilde a\).
- **Payment ordering:** if consumption, tax, or mortgage repayment must
  precede the sale, the corresponding payment-liquidity restriction returns.
  The relaxation requires simultaneous settlement or costless netting of
  these same-date claims. If consumption and tax precede the sale but the
  mortgage can use its proceeds, then
  \(c+\tau_tP_t h\le Y_{i,t}+k/q\), equivalently
  \(a'_O\ge-d/q\); eliminating the permitted origination debt adds
  \(qa'_O+\phi P_t h\ge0\). If the mortgage must also be repaid before
  title sale, the stronger \(a'_O\ge0\) returns. Neither extra inequality
  is needed under the simultaneous settlement specified in (2).

The substantive assumption that a young owner must retain nonnegative
financial assets before title realization can therefore be removed. The
remaining model still requires the early home commitment, genuinely later
income, and an explicit boundary-settlement convention.

## Retained asymmetry at death and the modeling tradeoff

The simultaneous liquidation permission above applies to a **living young
household entering old age**. It does not extend to the old household's
terminal home sale. In the maintained old problem,
\[
e=a_e+P_{t+1}h_o,\qquad a_e\ge0,
\]
the retained home's proceeds at death are earmarked for the estate. Old
consumption and tax must be financed without those terminal sale proceeds;
the household cannot use them in the same boundary settlement to finance
current consumption. Old households can still sell and resize at the
beginning of their stage, as already allowed.

Thus the two ages share a stage-end goods-delivery date but **do not share
the same access to terminal liquidation receipts**. Young sale proceeds
may fund current consumption; old death-sale proceeds may fund only the
estate. This is an explicit restriction on settlement and pledgeability,
not a consequence of logarithmic warm-glow preferences alone. If old
households were given the same terminal sale-financing permission, the
estate floor would cease to follow from these accounts, and the existing
old-value formulas and theorems would require a new derivation.

Removing the young pre-sale solvency condition is therefore not an
unqualified simplification. It permits sale-funded mortgage repayment but
retains an age-specific estate-earmarking rule. Keeping the original rule
that current consumption, tax, and mortgage payments precede title
liquidation is a cleaner uniform ordering for a short first illustration,
although it explicitly imposes repayment from nonhousing receipts on young
owners. Neither choice changes the verified all-young-renter branch, where
the stronger young condition is already strictly slack. If the relaxed
version is used, the old terminal-earmarking restriction must appear beside
the young settlement permission, rather than being left in an appendix or
described as fully symmetric timing.

## Verification

The forward derivation, constructive reverse implementation (7), and renter
replication (10) were checked directly. No numerical model run is needed
for these accounting identities.


---

# Part V. Independent review

## 8. Updated review: synchronized common-calendar candidate

Independent hostile review, September 9, 2026, of
`common_calendar_model.md`. This is a materially changed candidate, developed
after the comparison above. The following judgment updates the provisional
candidate ranking; it does not repair the earlier asymmetric-calendar model
by relabeling its dates.

**Verdict.** The new construction passes as an analytically nonempty,
fixed-fertility stationary competitive-equilibrium illustration with positive
taxes, equal funded rebates, endogenous price, heterogeneous income and
initial cash, all young renting at the upper rental cap, and all old owning
larger homes. Its common consumption-delivery date also supports the stated
full dated planner settlement. It is now the strongest static candidate of
the alternatives reviewed here: it needs neither an owner minimum nor high
old labor income or a slack estate floor. It remains an explicit model
variant, not a theorem about the maintained lifetime-tenure specification.

### 8.1 Individual and fiscal accounts: pass

Let `w` be cash received before the home commitment, `Y` income received at
the end of the current stage, and `T` the upfront rebate. The young owner
accounts

\[
 k+Ph=w+T+d,\qquad
 c+a'+\tau Ph=Y+(k-d)/q
\]

give exactly

\[
 q(c+a')+(1+q\tau)Ph=w+T+qY,
 \qquad (1-\phi)Ph\le w+T.
\]

The projection requires both the deposit condition and the separately
imposed `a' >= 0`. The latter is stronger than mortgage repayment alone if
the house can be sold at the same boundary to meet that payment. This
departure is stated correctly in the source, and is slack in its realized
all-renter allocation. None of the construction's received income may be
excluded from `w`; the mechanism requires that `Y` actually arrives after
commitment. Letting households wait to buy while obtaining the same entire
current-stage housing service would change this timing friction.

Write `A=1+q tau` and `D=A-q`. The renter's end-stage payment `DPh/q` gives
an intermediary buying and eventually reselling at `P` the bond return
after paying the property tax. There is no double property-tax charge.
The rebate `T=q tau P H/2` is correctly financed by beginning-of-stage
government borrowing of `2T`, repaid from end-stage taxes `tau P H`.
It is included in purchase cash and in the equilibrium price calculation.

For old owners, let `a_e` be the financial estate delivered at the end,
so total estate is `e=a_e+Ph`. The identity

\[
 qc_o+qe+DPh_o=Z,\qquad e\ge Ph_o
\]

is equivalent to buying the house from `Z` and holding the remaining
funds in bonds. In particular, positive consumption and the estate floor
imply `A P h_o < Z`, so a separate initial old purchase constraint does not
invalidate the solution. With `K=1+gamma+omega`, the unconstrained-size
owner solution is correctly

\[
 c_o=\frac{Z}{qK},\qquad
 h_o=\frac{aZ}{KP},\qquad
 a=\min\left\{\frac\gamma D,\frac{\gamma+\omega}{A}\right\},
 \qquad e=\max\{\omega,qa\}\,c_o.
\]

The financial estate is zero when `omega D < q gamma`; this is the strict
regime of the explicit family. The owner value remains `K log Z + constant`
there. This homogeneity is a valid upper bound when the owner size cap is
relaxed, and is not a global upper bound on the renter alternative.

### 8.2 Global tenure comparison: pass

The replication argument uses the correct direction. If every feasible
young owner house is below `H_R`, a renter obtains the same current
`(c,h)` and future old resources by taking `a'_R=a'_O+Ph`. This is feasible
and preserves continuation opportunities because tenure is reselected in
old age. With no ownership taste and a strictly optimal young rental size
`H_R`, every current owner plan is strictly dominated. An unbounded
ownership taste would defeat this exclusion; such a taste is explicitly
absent in this candidate.

The future-renter comparison has been checked from the joint two-age
problem rather than from an extrapolated capped value. Put `r=H_R`, let
`x` be the candidate young adult consumption, set
`G=1+beta K`, `G_R=1+beta(1+omega)`, and let `t=h_o/r>1`.
The joint renter-renter optimum is

\[
 x_{RR}=\frac{Gx-DPr}{G_R},\qquad
 \frac{x_{RR}}x=\frac{G-a\beta D/t}{G_R}>1.
\]

Its young rental cap follows from the candidate's strict cap condition.
Its old rental cap condition is exactly

\[
 \beta\gamma x_{RR}\ge DPr
 \quad\Longleftrightarrow\quad
 t\ge aD/\gamma,
\]

using `G=G_R+beta gamma`. It therefore holds because `aD <= gamma` and
`t>1`. The separately verified positive saving makes this relaxed optimum
attainable. Even without positive saving it would still be an upper bound
on that alternative, provided these cap checks were retained.

The optimized lifetime value difference is exactly

\[
 \Delta=\beta\gamma\log t+
 \beta\omega\log\!\left(\frac{\max\{\omega,qa\}}\omega\right)
 -G_R\log\!\left(\frac{G-a\beta D/t}{G_R}\right).
\]

The lower bound in equation (22) of the source follows from
`log(1+z) <= z`. Its bracket increases for `t>1`. Thus the positive bound
in equations (31)–(32) is genuinely global over the specified family.
It proves the optimal old tenure as well: a better renter continuation at
the candidate's selected `Z` would create a better feasible lifetime
renter-renter plan, contradicting this comparison.

### 8.3 Price and finite primitive family: pass, with one minor clarification

Market and fiscal clearing give the displayed price directly, with

\[
 J=G(H-r)+a\beta\left[\frac Dq r-
                         \frac{(1+q)\tau H}{2}\right],\qquad
 P=\frac{a\beta\bar b}{J},\qquad
 h_i^o=H-r+\frac JG\left(\frac{b_i}{\bar b}-1\right).
\]

No equilibrium price or ownership share is assumed as a primitive.
The current old distribution is the preceding young cohort's savings
distribution. The estate is warm glow and goes outside the entrants'
cash account, as stated; it is not simultaneously counted as their `w`.

The finite bounds in equations (26)–(32) have been independently checked,
including the late expansion to `beta in [3/10,1]`. For example,
`J/G < 3/2+(250/201)(4/13) < 19/10`; this yields
`281/200 < h_i^o < 319/200`, hence `7/5 < h_i^o < 8/5`.
The updated young-cap upper bound `255/202` remains below `7/5`.
The origination bound
`33283/180000 < 1/5` includes the rebate. The rental cap, old owner cap,
positive saving in both alternatives, and the rational lower bound on
the value gap all have strict analytical margins. The price is unique
within the verified regime; other equilibria are not excluded.

For the nonemptiness statement, fixed fertility satisfying both
`kappa n_i < 1` and `nu E n_i=1` requires `nu>kappa`. State this primitive
restriction, or explicitly choose `nu` as part of the family. Conditional
on it, the proposed construction of positive compact heterogeneous `w`,
large mildly heterogeneous `b`, and `Y=b-w/q+chi n` is valid.

### 8.4 Full dated planner: pass with an explicit rental account

The source's equation (33) now has a common economic date: both ages'
current consumption is delivered at the end of this stage. Hold every
current young household's next-stage total old resources `Z` and every
current old household's total estate `e` fixed. Once individual financial
bounds are relaxed, the beginning transfer needed by either age is

\[
 t_i=q\,\Delta c_i+DP\,\Delta h_i.
\]

Hence these transfers sum to zero when current end-stage goods and
occupied housing are conserved. Young future tenure opportunities depend
on `Z`, so changing its current tenure does not change its continuation
value. The planner intervenes before the current home commitments.

The following financier account is necessary when the planner changes
tenure. A young renter preserving `Z` has `Delta a'=0`; an old renter
preserving `e` has zero change in financial estate. For household-owned
title changes, however, the young future financial holding or the old
financial estate changes by `-P Delta h`. Thus household end claims change
by `-P Delta H_O`. A rental intermediary's residual end financial claims,
after rents, taxes and funding service but before title sale, change by
`-P Delta H_R`. Because total floor space is fixed,

\[
 \Delta H_O+\Delta H_R=0,
 \qquad
 \Delta F_{\rm households}+\Delta F_{\rm rental}=0.
\]

Discounted to the beginning date these compensating claim changes are
`-qP Delta H_O` and `-qP Delta H_R`. This is not an assertion that the
rental firm's initial gross purchase loan is only `qP` per unit: that
gross loan is `P`, with the stated rental cash flow covering its funding
service. Government receipts and debt are unchanged when total housing
and the tax rate are fixed. No positive liquid estate of the initially
floor-constrained old is required: the full planner may replace title
with claims and relax private financial restrictions.

The technology must be the stated common divisible floor-space stock,
with no immutable owner/renter sector inventory, and the planner must be
allowed to reassign current tenure. All large homes can then use the
owner menu while its upper cap `H_O` remains in force. An interpretation
with fixed separate physical inventories would require another theorem.

### 8.5 Housing follows; fertility and transition do not yet follow

With common housing weights `alpha=gamma=1`, no ownership taste, positive
child space, equal current age masses, and the common feasible menu, a
full planner cannot assign a larger home to an old household than to a
young household: exchanging their houses while preserving their other
allocations strictly improves housing utility. This argument retains any
selected positive fertility, so it also applies to a dated planner that
chooses fertility jointly. It is stronger than checking only the
competitive marginal-utility gap.

Full use of `H=5/2` is possible and optimal because each of the two
adult households can use up to `H_O=2` and housing utility is increasing.
Age ordering therefore gives planner mean young housing at least
`H/2=5/4`, strictly above the competitive rental size `H_R=1`. This proves
the aggregate young-housing gain, with both consumption and housing
chosen by the planner, rather than inferring it from a single favorable
pair exchange. For fixed fertility, the housing first-order conditions
provide an equivalent check: old size is a common `s`, young size is
`min{H_O,s+kappa n_i}`, and their sum uses `H`.

The exact competitive family has fixed or inherited fertility. Its
`vartheta log n` term is constant in private choice. Consequently it does
not establish private fertility optimality or sign a planner change in
mean fertility. Even the current consumption ordering is conditional:
`c_o/x_y=beta/q`, while the specified beta interval crosses `q=1/2`.
A consumption cushion used in a separate endogenous-fertility result
cannot simply be presumed across this entire region. Changing fertility
also changes future population; the current dated resource comparison
does not by itself establish a feasible stationary destination or a
general-equilibrium transition.

Finally, the realized young households rent and save strictly. The
mechanism is purchase exclusion from cash available at closing, followed
by saving from later income. It does **not** supply a positive mass of
realized young mortgage borrowers. Its most important remaining
limitations are therefore the explicit timing/tenure/taste departures
and the missing endogenous-fertility competitive extension, not an
algebraic failure of the fixed-fertility equilibrium or planner ledger.

### 8.6 Late uniform-tenure lemma: conditional pass

Section 6 of the candidate, added during this review, correctly strengthens
the gross-resource restrictions to verify tenure uniformly over individual
fertility deviations. If `g_i=w_i/q+Y_i`, gross dispersion is at most one
percent, and `chi r/kappa <= mean(g)/50`, then every feasible renter
fertility has net resources in `[.97 mean(g),1.01 mean(g)]`. The actual
candidate mean net resource is in `[.98 mean(g),mean(g)]`, so the ratio of
any deviating individual's net resources to that actual mean is within
`3/98 < 1/20` of one. Also `.98*46 > 45`, as required for purchase
exclusion. Thus all previously checked cap, saving and global tenure
inequalities hold at the candidate's actual price and transfer for every
such deviation. An individual deviation must not recompute the aggregate
mean or equilibrium price; the source correctly holds them fixed.

All owner-feasible current houses are below `r`, so their feasible
fertility is covered by the same range and the replication comparison
remains valid. Strict concavity of the relaxed future-owner problem then
identifies the global private fertility/continuation optimum conditional
on the candidate aggregate state. This lemma removes a possible
off-branch tenure-deviation objection to a future fertility construction.
It does not itself solve the aggregate fertility and replacement
equations or establish an endogenous-fertility equilibrium; the source
correctly separates those tasks.

No numerical model solve, browser search, build, or other-file edit was
used in this review.

## 9. Independent review of joint fertility and stationary closure

September 9, 2026. This section independently checks the newly completed
`common_calendar_fertility.md`, including its full joint-planner profile,
the cap-valid sufficient condition without `beta >= q`, and the
endogenous-fertility stationary construction. It supersedes Section 8's
statement that this fertility extension was still outstanding. All other
model departures and transition limits remain in force.

**Verdict: pass within the stated model and parameter scope.** No fatal
algebraic, cap, contraction-domain, or global-tenure gap was found. The
construction completes an endogenous-fertility competitive branch through
a single fertility-preference calibration and proves a full dated joint
planner increase in both young housing and mean fertility. It does not
establish a policy transition, global equilibrium uniqueness, or a
stationary population with arbitrarily and independently chosen fertility
preferences and demographic replacement.

### 9.1 The full joint-planner profile and cap: pass

Continuation resources and net estates are fixed individually, financial
restrictions are relaxed, and current tenure is assignable under a common
divisible stock. Identical current preferences and strict concavity then
justify equal young fertility and adult consumption, together with common
old consumption. Let `C` and `H` denote current resources per young-cohort
unit. For any candidate fertility `n`, adult consumption is

\[
 X(n)=(C-\chi n)/2.
\]

Equal young and old housing weights imply that old housing and young
adult space coincide until the young gross housing cap binds. Thus

\[
 S(n)=\min\{(H-\kappa n)/2,H_O-\kappa n\}.
\]

When the cap binds, young housing is `H_O` and old housing is the constant
`H-H_O`; the derivative of old housing utility is then zero. When the cap
is slack, both young adult space and old housing move at rate `-kappa/2`.
Both branches consequently give the same optimized fertility derivative

\[
 F_J(n)=\theta/n-\chi/X(n)-\alpha\kappa/S(n).
\]

It is continuous at the cap switch, strictly decreasing, and has opposite
infinite limits at the feasible endpoints stated in the source. The
unique root is the full joint optimum. Both displayed quadratic equations
are the correct branch-specific forms of this derivative; feasibility and
the cap inequality select the applicable root. The old cap is slack
because young housing is weakly larger and total stock is below `2H_O`.

At reference mean fertility `n_0`, total adult-consumption resources equal
`(1+beta/q) mean(x)`. Therefore equation (3) is the exact global comparison
for `n_J>n_0`, including a potentially binding young planner cap. At the
joint optimum the common young home is strictly larger than `H/2>r`.
This uses the complete optimum, not a private fertility response after a
partial allocation experiment.

### 9.2 Jensen and the discount-independent condition: pass

Solving the private fertility first-order condition at adult consumption
`x` and gross housing `h` gives the source's function `N(x,h)`. In units
`C=x/chi`, `R=h/kappa`, its discriminant is

\[
 R^2+2(\alpha-\theta)CR+(\alpha+\theta)^2C^2.
\]

The matrix is positive definite because its determinant is
`4 alpha theta>0`. Its square root is therefore a norm, so subtracting it
from the linear part establishes joint concavity and homogeneity. With
fixed rental housing `r`, it is strictly concave in `x`; Jensen is strict
if adult consumption is nondegenerate. Applying the decreasing private
fertility residual to the Jensen inequality gives equation (4) with the
correct direction.

The consumption cost of redistribution when `beta<q` is exactly

\[
 \frac{\chi}{\bar x}\frac{q-\beta}{q+\beta}.
\]

Subtracting this from the housing-space gain gives condition (5), which
is sufficient under heterogeneity and exact under a homogeneous
reference. The lower bound used for condition (6) is valid on both cap
branches: each reciprocal-space gap decreases as reference adult space
rises to `r`; the second branch also uses the necessary reference fact
`H_O >= m`. Thus

\[
 \frac1{s_0}-\frac1{S_0}
 \ge\frac{m-r}{r(m+r)}.
\]

Condition (7) is consequently sufficient for every positive discount
factor at a reference satisfying it. This means the welfare condition
does not require an ordering of `beta` and `q`; it does not assert
competitive existence for all positive `beta`. The constructed
competitive family retains `beta in [3/10,1]`.

### 9.3 Private fertility and the scalar price method: pass

Differentiating the reduced future-owner objective gives

\[
 \theta/n-G\chi/(M-\chi n)-\alpha\kappa/(r-\kappa n)=0.
\]

Multiplication verifies equation (8). Its feasible root `f_theta(M)` is
strictly increasing in resources, with

\[
 0<f'_\theta(M)=
 \frac{G\chi/(M-\chi n)^2}
 {\theta/n^2+G\chi^2/(M-\chi n)^2+
                  \alpha\kappa^2/(r-\kappa n)^2}<1/\chi.
\]

It is also strictly increasing in `theta`.

For the general fixed-`nu,theta` scalar method, substitution of the funded
rebate into mean young resources gives directly

\[
 \widetilde G\,\bar x=S-C_0Pr.
\]

This reproduces the displayed `M_i(P)`, and every individual term in the
fertility equation decreases strictly with price. The condition
`tau<2/(1-q)` indeed ensures `mu_0<beta(alpha+omega)` and
`tilde G>1+beta`. The inversion level `M_*` and the two endpoint tests have
the right signs and establish existence on a bracket when their stated
inequalities hold. At a root, feasibility of individual consumption
also ensures positive reconstructed mean consumption. As the source
emphasizes, this conditional scalar method alone does not verify tenure
or guarantee that such a bracket exists for arbitrary primitives.

### 9.4 The entire contraction domain and replacement calibration: pass

In the normalized family, `E_r=1-7 tau/8` is positive and

\[
 k=\frac{a\beta E_r}{J}=1-\frac{m}{J/G}.
\]

The verified bound `J/G<2`, with `m=3/2`, implies `0<k<1/4`.
For every proposed mean fertility `z in [0,r/kappa]`, the source correctly
constructs

\[
 M_i(z)=g_i-k(\bar g-\chi z).
\]

In particular, the worst case is the lower endpoint, and the primitive
support restrictions give the uniform strict margin

\[
 M_i(z)>(.99-.25)\bar g=.74\bar g
 >.02\bar g\ge\chi r/\kappa.
\]

Thus the entire map is defined with positive goods remaining even at
the largest possible private fertility. It is not merely defined near
an eventual fixed point. Each finite positive `theta` maps the closed
interval into its interior, and

\[
 |\mathcal T_\theta(z_2)-\mathcal T_\theta(z_1)|
 \le k|z_2-z_1|.
\]

The contraction is uniform in the distribution of types. Compact support
and the derivative bound justify aggregation. It gives a unique
stationary mean fertility for each `theta`, with positive prices and
rebates throughout the construction.

The unique fixed point increases strictly and continuously with `theta`.
For completeness, its endpoint limits use the goods margin, not simply
monotonicity. As `theta` tends to zero, the private first-order condition
bounds `n` above by `theta M/(G chi)`, uniformly over the bounded resource
support. As `theta` tends to infinity, any fixed distance below `r/kappa`
would keep both private marginal child costs bounded while `theta/n`
diverges. Hence convergence to both endpoints is uniform and passes to
the fixed point. Every target `1/nu` strictly inside `(0,r/kappa)` is
therefore attained by exactly one positive `theta`.

The necessary scope sentence is:

> With cohort mass fixed at one, the fertility preference `theta` is
> calibrated to demographic replacement; `theta` and `nu` cannot both be
> specified independently while retaining that normalization.

This is a one-dimensional calibration restriction, rather than an open
parameter box in all independently specified primitives. It does not
select the welfare sign: the tenure regime and welfare inequalities
hold uniformly in `theta`. Alternatively, choosing `theta` first and
setting the compatible `nu` is exactly equivalent, as stated.

### 9.5 Global deviations and the mean-consumption bound: pass

At the fixed point, actual mean net resources `S` lie in
`[.98 mean(g),mean(g)]`. Every individual's feasible fertility deviation
has net resources in `[.97 mean(g),1.01 mean(g)]`, hence lies within
`3/98<.05` of the actual mean after division by `S`. The accompanying
price and rebate must be held fixed during an individual deviation;
the source does this correctly. The already verified fixed-fertility
tenure, rental-cap, old-cap and saving inequalities therefore apply
uniformly to every such deviation. Current-owner replication covers
all of its feasible fertility choices because every owner-feasible
house is below `r`. This rules out joint tenure/fertility deviations,
rather than checking fertility only on an assumed branch.

Finally, the upper-beta expansion is harmless for the new welfare bound:
`G<=13/4`, `a beta<5/4`, and `E_r<1` give `J<49/8`. Consequently

\[
 \bar x=\frac{mS}{J}
 >\frac{(3/2)(49/50)}{49/8}\bar g
 =\frac6{25}\bar g,
 \qquad
 \frac{\chi r}{\kappa\bar x}<\frac1{12}<\frac15.
\]

The strict inequality is valid even when the primitive goods-cost bound
holds at equality. The last fraction is exactly `(m-r)/(m+r)` in this
family. Thus the full dated joint-planner fertility increase is proven
throughout the constructed family, including its `beta<q` portion, with
the owner physical cap retained.

The comparison remains a dated redistributive benchmark over current
households, with their continuation resources and estates preserved.
Its higher fertility is not a proof of a feasible policy transition or
a new stationary population. No numerical model run, root search,
simulation, build, browser, or other-file edit was used for this check.

## 10. Independent audit of the same-model local transition

September 9, 2026. This section audits `common_calendar_transition.md`
directly from the synchronized dated budgets. Its priority is the
three-state recursion and local stability proof, including unexpected
taxation of an inherited state. It is not a heterogeneous-path theorem.

**Verdict: pass as a local homogeneous transition theorem.** The saving
state, current and forward resource equations, three-dimensional
linearization, characteristic polynomial, and strict unit-disk bounds
are correct. The stated taste and permanent-tax impact and endpoint signs
also pass. The late-added price distinction is necessary and correct:
the tax lowers the purchase price on impact but raises it at the
stationary destination in this explicit family. No fatal gap was found
within this scope.

### 10.1 Inherited saving and actual surprise rebates: pass

The true beginning-of-stage state is `(Y,O,s)`, where `s` is financial
saving inherited by each old household from its young rental period.
Current old resources are `Z=s+T`; the current rebate is determined by
the actual tax and the actual current cohort masses. Combining

\[
 P=\frac{aO(s+T)}{K(\bar H-rY)},\qquad
 T=\frac{q\tau P\bar H}{Y+O}
\]

gives equation (T4) exactly. At an unexpected tax change, replacing the
actual `Z` with `beta K x_previous` would be wrong because the preceding
saving decision used the anticipated rebate. The source correctly avoids
this substitution on impact.

The old owner floor uses the next sale price: `e=P_next h_o`.
Substituting this floor into its dated present-value budget leaves
`q c_o+A P h_o=Z`, so its housing choice and resource derivative are
independent of the next sale price. That price remains in the additive
estate-utility term. This is a consequence of the binding floor, not a
stationary-price shortcut. The financial floor and positive rental user
cost have strict margins at the reference and remain valid on sufficiently
small local paths.

### 10.2 The exact forward closure: pass

The young renter's dated cash account is

\[
 x+\chi n+s_{t+1}
 =B+T_t/q-(A_tP_t/q-P_{t+1})r.
\]

Because consumption and financial saving are delivered at the same end
date, its positive-saving first-order condition is
`s_next+T_next=beta K x`, without an additional factor of `q`.
Substitution gives (T8), including both the current rebate and the
anticipated rebate, rent and resale terms.

Next old mass is current young mass, next young mass is `nu n Y`, and
next old resources are `beta K x` under post-shock foresight. Therefore
next market clearing gives

\[
 P_{t+1}=\frac{a_{t+1}\beta Y_t x_t}
                  {\bar H-r\nu n_tY_t}.
\]

The next rebate is then (T9). Dividing `P_next r+T_next` by `x` produces
the coefficient `M` in (T10) exactly. Thus the scalar current-choice
equation retains the full future resource terms. Future young housing is
known to be `r` on this strict branch; it has not been replaced by an
arbitrary expenditure rule.

The scalar root is locally regular. At the reference,
`M_* <= 27 beta/32` and `M_v <= 107 beta/192` are valid bounds, giving
`G-M_*-M_v >= 1+163 beta/192`. Since `n X_n/X>1`, differentiating
`(G-M)X+chi n` gives a strictly positive derivative. The current price
denominator is also strictly positive. Hence the nearby forward map is
smooth. This does not claim global uniqueness of admissible roots.

### 10.3 The full three-state linearization: pass

The source's temporary use of current old resources as a state coordinate
is legitimate. With current masses fixed,

\[
 Z=\frac{s}{1-
 aOq\tau\bar H/[K(\bar H-rY)(Y+O)]},
\]

so `dZ/ds>0` in the verified neighborhood. The transformation from
`(Y,O,s)` to `(Y,O,Z)` is locally invertible. An unexpected tax changes
this coordinate through the actual rebate; it does not change inherited
`s`. Thus no initial financial degree of freedom is silently removed.

Let `ell=delta Y/N`, use `z=delta x/x`, and let `e=(x/n)n_x` lie in
`(0,1)`. Direct differentiation of market clearing gives
`p_t=ell_previous+z_previous+h ell_t`, while the demographic law gives
`ell_next=ell+e z`. The resource linearization before substitutions is

\[
 (G+\chi n e/x)z_t
 =-u p_t+v p_{t+1}
 -\frac{w}{2q}(\ell_t+\ell_{t-1})
 -\frac w2(\ell_{t+1}+\ell_t).
\]

The last two terms are the changes in rebate denominators as cohort
masses change. Substituting the price and demographic identities gives
(S4)–(S5), with exactly the stated `Lambda` and `D`. Omitting either
denominator term would change the stability polynomial.

In coordinates `(ell_t,ell_previous,z_previous)`, the first row of the
Jacobian is the second row plus `e` times the third row. Its zero
eigenvalue is therefore real, rather than an assumed discarded mode.
For the nonzero modes, eliminating
`z_previous=(ell_t-ell_previous)/e` yields (S6) and precisely the quadratic
`lambda^2-A_2 lambda-B_2`. This is the full three-dimensional Jacobian
factored into its zero root and two remaining roots.

### 10.4 The entire stated beta interval is locally stable: pass

At `q=1/2`, `h=2/3`, the identities `v=M_*` and
`v h-w/2=M_v` give

\[
 \mathcal D\ge1+163\beta/192.
\]

The conservative bound `u <= 101 beta/60` is valid (the exact relation
`u=5 beta/3-2w` is even sharper). Consequently

\[
 \mathcal D-u\ge1-267\beta/320\ge53/320
 \quad\text{for }\beta\le1.
\]

The bound `w<=beta/96` and the displayed positive expression for `Lambda`
also check. Each of the three quadratic unit-disk inequalities is then
strict:

\[
 1-A_2-B_2
 =\frac{e a\beta h(1+h)C_0}{\mathcal D}>0,
\]

\[
 1+A_2-B_2
 =\frac{2(\mathcal D-u)+e(\Lambda+u+w/(2q))}
        {\mathcal D}>0,
\]

\[
 1+B_2>1-\frac{w}{2q\mathcal D}>0.
\]

These conditions place both roots strictly inside the unit disk. The
zero root causes no obstruction: local asymptotic stability and the
contraction-norm argument for a smooth forward map do not require the
map to be invertible. Small permanent shocks have a nearby stationary
point and a unique nearby forward path from the inherited state. The
strict branch inequalities and stable neighborhood keep that path on
the verified choice branch and make its prices consistent with foresight.

The local neighborhood may depend on the particular calibrated reference.
For example, `e` can approach zero and a root can approach one. The proof
does not supply a uniform quantitative shock size or convergence rate
over all possible fertility-preference calibrations. Nor does it exclude
negative or complex roots or prove monotone adjustment.

### 10.5 Impact and terminal signs: pass, with the price reversal retained

An unexpected taste change leaves current `(Y,O,s,tau)` unchanged, hence
also current price and rebate. In the scalar equation, `X_theta<0`,
`G-M>0`, and the fertility-root derivative is positive. Thus a small taste
decline lowers fertility on impact and next young mass. At replacement
fertility in the destination, it raises adult consumption, lowers price,
enlarges the old home and lowers stationary population. Fertility at both
positive endpoints remains `1/nu`; later fertility signs are not implied.

For an unexpected permanent tax at a stationary inherited state, the
source's current-price denominator is

\[
 D_o=Km+q\tau[Km-b(r+m)/2].
\]

Differentiating the full current resource term gives exactly (I1), with
positive factor `Km-br=m+b(m-r)`. Differentiating the anticipated
next-date coefficient gives (I2), positive at `nu n=1` because
`r/(r+m)<1/2`. Both effects therefore raise impact fertility in (I3).
They include the surprise rebate to inherited old resources and the
anticipated next-date tax; neither effect has been omitted.

The stationary population derivative (E4) also checks directly. Its
numerator has sign `U+(1+q)Q/q`, which is positive under the stated
positive-price condition. Thus the permanent tax raises the destination
population and lowers the old mean home within this floor-binding branch.

Purchase-price timing must remain explicit. At impact,

\[
 \frac{P_{t,\tau}}{P_t}
 =-\frac{q[Km-b(r+m)/2]}{D_o}<0,
\]

because the bracket equals `m+b(m-r)/2`. At the stationary destination,

\[
 \frac{P_{*,\tau}}{P_*}
 =\frac{(1+q)m}{2(1+q\tau)C_0r}
   -\frac{1-q}{2C_0}>0
\]

in this family. This independently confirms the late-added equation
(E6): the tax lowers the impact price and raises the terminal price.
The actual impact rebate nevertheless rises because
`d(tau P_t)/d tau=P_t Km/D_o>0`. There is no implied monotone price path.

### 10.6 Exact scope for the packet

A concise theorem can claim a unique nearby homogeneous forward path
after sufficiently small permanent shocks, a fertility decline on impact
and lower destination population after a taste decline, and higher impact
fertility and destination population after a subsequent small unexpected
permanent tax increase from the same inherited state. The strict signs
extend to sufficiently nearby states along the first baseline path by
continuity.

The shocks arrive at stage openings before new housing commitments.
Initial old households inherit financial saving because they previously
rented. Young housing remains `r` throughout the competitive paths;
the policy's larger young aggregate housing share comes from population,
and is distinct from the dated planner's per-young housing increase.
Both tax paths converge to replacement fertility, so their terminal
population difference records cumulative transition fertility. No common
sign is proved for all later fertility gaps.

The transition is homogeneous. The separately valid heterogeneous
endpoint formulas do not justify replacing a heterogeneous transition
distribution by its mean. Population here counts adult households; a
resident-person claim needs a further counting convention. There is no
global-regime uniqueness, large-shock, heterogeneous-path, or numerical
population-forecast result. No model run, root search, simulation, build,
browser, or other-file edit was used for this audit.

## 11. Independent audit of the heterogeneous transition extension

September 9, 2026. This checks the newly added Section 7 of
`common_calendar_transition.md`. It supplies a separate review of the
heterogeneous-path argument that was expressly excluded from Section 10.

**Verdict: pass under the stated fixed entrant distribution and uniform
individual-regime restrictions.** The aggregate elasticity bound, inverse
schedule taste derivative, mean-state closure and transfer of the earlier
stability and impact equations are correct. This extends the checked local
transition to the stated heterogeneous family; heterogeneous stationary
endpoints alone would not have provided that extension.

### 11.1 Exact common-shift aggregation: pass

The argument uses a fixed exogenous joint distribution of entrant cash
and later income in each cohort, with common preferences and child costs.
At current and anticipated prices, every young type faces the same
additive shift `C`, so

\[
 M_i=g_i+C=Gx_i+\chi n_i,\qquad n_i=f_\theta(M_i).
\]

Let `n_C=E f_M` and let `n_theta|C=E f_theta` denote derivatives of mean
fertility. The household derivative satisfies `0<f_M<1/chi`, hence

\[
 \bar x_C=(1-\chi n_C)/G>0.
\]

This permits a local inverse from mean adult consumption to the common
shift. At fixed mean adult consumption,

\[
 C_\theta=\frac{\chi(n_\theta|C)}{1-\chi n_C},
 \qquad
 \mathcal N_\theta
 =\frac{n_\theta|C}{1-\chi n_C}>0.
\]

Equation (H3) therefore has the correct denominator and sign. Its inverse
has `mathcal X_theta<0`. This schedule integrates actual household
policies; it does not evaluate an individual policy at mean resources.

### 11.2 The aggregate elasticity bound: pass

Define `d_i=partial n_i/partial x_i` from the private fertility condition
at the common rental cap. Direct differentiation gives (H4), and the
first-order condition implies

\[
 d_i<\frac{\chi r}{\alpha\kappa}\frac{n_i}{x_i^2}.
\]

Because `M_i=Gx_i+chi n_i`, the response of adult consumption to the common
shift is `rho_i=1/(G+chi d_i)`. Thus the aggregate elasticity is exactly

\[
 e=\frac{\bar x}{\bar n}
       \frac{\mathbb E(\rho_i d_i)}{\mathbb E\rho_i}.
\]

The weighting matters; individual elasticities below one would not by
themselves prove that this aggregate elasticity is below one.

The stated primitive support and child-goods bounds, together with
`k<1/4`, give

\[
 x_i>\frac{18\bar g}{25G}=x_{\min},\qquad
 \bar x\le\frac{101\bar g}{100G}.
\]

With `alpha=1`, `n_i<r/kappa`, and `G<=13/4`, this yields

\[
 \chi\max_i d_i
 \le\left(\frac{\chi r}{\kappa x_{\min}}\right)^2
 \le(13/144)^2<1/100.
\]

Consequently `rho_i` is between `1/(G+1/100)` and `1/G`. Applying these
weight bounds and then the bound on `E d_i` gives

\[
 e\le\frac{\chi r}{\alpha\kappa}
      \frac{\bar x}{x_{\min}^2}\frac{G+1/100}{G}
 \le\frac{(1/50)(101/100)(13/4+1/100)}{(18/25)^2}
 =\frac{16463}{129600}<\frac{13}{100}.
\]

The rational constant is correct. In particular, `0<e<1` and
`bar n * mathcal X_n/bar x=1/e>1`. The reference has a strict upper
elasticity margin that persists under small changes in the common shift
and parameters, uniformly over the compact type support. This does not
provide a positive uniform lower bound on `e`, or a uniform convergence
rate across all fertility-preference calibrations.

### 11.3 Why the mean state suffices, and what it does not omit

All current old households have the same old utility coefficients and
are uncapped owners with a binding financial-estate floor. Their housing
demand is exactly linear in resources:

\[
 \bar h_o=\frac{a}{KP}(\bar s+T).
\]

Current aggregate housing and the fiscal equation therefore depend on
the inherited savings distribution only through its mean. Every young
type also obeys `s'_i+T_next=beta K x_i`, so the next mean savings is
`beta K mean(x)-T_next`. Demography depends on actual mean fertility.
The aggregate state is consequently `(Y,O,mean(s))` and its resource
identity has the same form as before.

This is not a claim that arbitrary savings distributions with the same
mean are individually feasible. The initial full savings profile must
be the constructed profile or uniformly close to it, as stated in the
extension. Mean proximity alone could conceal a type hitting its housing
cap, changing tenure, or losing positive resources. The previous uniform
type-level margins rule these changes out under the specified stronger
proximity condition. The estate-floor regime itself has a strict price
inequality even though the financial estate equals zero.

Subsequent individual profiles are continuous functions of the type,
common resource shift, taste parameter and anticipated rebate. Compact
type support and the nearby aggregate path preserve those margins
uniformly. Initial mean-zero old-profile differences have no aggregate
demand effect on the linear branch and do not survive that old cohort's
exit. Their estates do not feed the entrants' cash distribution.

The exogenous entrant distribution is essential. If parental fertility
selected entrant types, or estates determined their initial cash, a new
cohort would generally have a different distribution and this three-state
closure would fail. Both channels are explicitly excluded from the
candidate, rather than silently absorbed into a representative agent.

### 11.4 Stability and impact equations genuinely carry over

The exact mean version of (T10) uses `mean(g)` for gross income and the
aggregate inverse `mathcal X(mean(n),theta)` for adult consumption. The
next-price formula continues to use `beta K mean(x)` because old demand
is linear. The previous scalar regularity proof needs only
`mean(n) * mathcal X_n/mean(x)>1`, now independently established.

The complete resource linearization (S4)–(S5) is unchanged: all rebate
denominators depend on cohort masses, while aggregate child goods enter
as `chi mean(n)`. The demographic derivative is
`ell_next=ell+e z`, with the aggregate elasticity just derived. The
earlier three-dimensional Jacobian, its zero eigenvalue and the two
quadratic roots therefore apply with this `e`. Every unit-disk bound
in (S7)–(S8) uses only `0<e<1`, so it remains valid through `beta=1`.

At fixed mean fertility, the aggregate inverse has a negative taste
derivative, preserving the taste-impact sign. It has no direct tax
dependence: taxes enter the common resource shift, while the function
itself depends on the fixed entrant distribution and preferences. Thus
both terms in the permanent-tax impact derivative remain exactly those
already checked. Current price still uses actual surprise rebates and
the inherited mean saving, with the full inherited profile held fixed
between policy and baseline. The impact-price decline and terminal-price
increase therefore survive as well.

The valid conclusion is a local heterogeneous transition in the stated
one-percent resource-dispersion family, with recoverable individual
allocations and uniform regime feasibility. It retains the earlier
limits: young homes remain at the rental cap, no sign is proved for all
later fertility differences, and convergence need not be monotone.
No claim about endogenous entrant selection, large shocks or global
equilibrium uniqueness is added. No model run, root search, simulation,
build, browser, or other-file edit was used in this check.

## 12. Independent audit of repayment from boundary sale receipts

September 9, 2026. This checks the proposed relaxation of the young
owner's separate nonnegative financial-asset condition, directly from
the dated purchase and settlement accounts. Only this audit section is
edited; the proposed implementation is being documented separately in
`repayment_relaxation.md`.

**Verdict: conditional pass; the age-specific terminal restriction must
also be explicit.** The condition
`a'_O >= 0` is unnecessary for the verified all-young-rental equilibrium,
global tenure comparisons, dated planner results and local transitions.
It can be replaced by `a'_O+P_next h >= 0` if income, bond payoffs, house
sale, consumption, taxes and debt repayment clear at the same boundary.
This netting permission is for the young; retaining the old estate floor
requires a different terminal spendability restriction, detailed in
Section 12.4 below. There is no new borrowing in old age. This is a genuine enlargement of
the off-equilibrium young-owner opportunity set, not a claim that every
other regime of the expanded model has unchanged equilibria.

### 12.1 Exact projection and a constructive implementation

Let `P_t h` be the house purchased before young housing services are
delivered. Initial cash is `w+T_t`, mortgage debt is `d`, and initial bond
holdings are `k`. The closing conditions remain

\[
 k+P_t h=w+T_t+d,\qquad
 k\ge0,\qquad 0\le d\le\phi P_t h.
\]

At the end of the young stage, allow sale of the occupied house for
`P_{t+1}h`. Define the cash left after all boundary obligations by

\[
 b^+=Y_t+k/q+P_{t+1}h-c-\tau_tP_t h-d/q\ge0.
\]

This is the actual nonnegative financial position entering retirement,
before the next rebate. Define only for consolidated accounting

\[
 a'_O=b^+-P_{t+1}h.
\]

This quantity can be negative; it is not an unsecured loan remaining
after settlement. Eliminating the initial bond and mortgage gives

\[
 qc+qa'_O+A_tP_t h=w+T_t+qY_t,
 \qquad a'_O+P_{t+1}h\ge0,
\]

or, equivalently, with `u_t=A_tP_t-qP_{t+1}`,

\[
 qc+qb^++u_t h=w+T_t+qY_t.
\]

The origination condition is unchanged:
`(1-phi)P_t h <= w+T_t`. It cannot be financed by income or sale receipts
that arrive only at the boundary.

These projected conditions are sufficient as well as necessary. Given a
feasible allocation, choose

\[
 d=\max\{0,P_t h-w-T_t\},\qquad
 k=w+T_t+d-P_t h.
\]

The deposit condition gives `d<=phi P_t h`, and `k>=0` by construction.
At the boundary, full liquidation produces enough combined resources to
pay consumption, tax and the full mortgage face value, leaving exactly
`b^+>=0`. Full sale and repurchase, if desired in old age, are a feasible
normalization under the model's costless trade and freely reselected
tenure. Partial sale is unnecessary for this implementation. The old
household then starts with `Z=b^++T_{t+1}` and chooses its old home under
the unchanged no-new-borrowing and estate-floor restrictions.

### 12.2 The settlement timing must be stated, not inferred

The exact terminal-only projection permits sale receipts to fund current
end-stage consumption and property tax as well as mortgage repayment.
That is coherent with the synchronized common calendar, provided these
items clear together after the current housing service has been delivered.
It does not require an interim bridge loan or old-age debt.

If consumption and tax instead had to be paid before any sale receipts
were available, one would also require

\[
 c+\tau_tP_t h\le Y_t+k/q,
 \quad\text{equivalently}\quad a'_O\ge-d/q
\]

for the selected initial mortgage. The terminal collateral condition
alone would then not be an exact projection. If mortgage repayment also
preceded sale, still stronger interim cash conditions would arise.
Accordingly the proposed relaxation should name simultaneous boundary
netting explicitly. It must not quietly impose consumption-before-sale
liquidity after deriving the less restrictive budget.

### 12.3 Renter replication and all previously verified results

A young renter can reproduce any feasible owner allocation at the same
current `c,h,n` by saving

\[
 a'_R=b^+=a'_O+P_{t+1}h\ge0.
\]

Its budget is exactly the owner's projected budget because

\[
 qc+qa'_R+(A_tP_t-qP_{t+1})h
 =qc+qa'_O+A_tP_t h.
\]

The renter and owner enter old age with identical total resources
`b^++T_{t+1}` and face the same future tenure menu. This argument uses
the dated next sale price, not a stationary substitution. It also
preserves fertility and both ages' utility because no ownership taste
or transaction cost is present.

The existing origination bounds keep every affordable owner home
strictly below the rental ceiling `r`. The renter replication therefore
remains physically feasible for every newly admitted owner plan. The
already verified global renter optimum reaches `r` and strictly dominates
all such owner plans; its uniform comparison over fertility deviations
is unchanged. No newly feasible owner plan invalidates the global
tenure argument.

Every household in the constructed competitive reference and local
transition already rents when young and saves strictly. Its realized
policy and every aggregate equilibrium equation are therefore unchanged.
Old no-borrowing and the binding old estate floor are retained. The dated
planner already relaxed individual finance and used the same net title
and claim settlement, so its feasible comparisons are unchanged as well.

This conclusion is specific to the verified all-young-rental branch
and its sufficiently small paths. The expanded private model can have
different choices or equilibria in other regimes. In particular, no
general default or unanticipated-price-solvency theorem for actual
young mortgage borrowers is supplied here. No model run, root search,
simulation, build, browser, or reader/source edit was used in this check.

### 12.4 Late consistency qualification: old estate protection is substantive

The unrestricted young-boundary netting above cannot be described as a
uniform settlement rule for both ages while preserving the existing old
solution. That solution retains

\[
 e=a_e+P_{t+1}h_o,\qquad a_e\ge0.
\]

It therefore prevents housing proceeds from funding the old household's
end-stage consumption. If old households could also liquidate after
receiving the full housing service and use the receipts to consume before
leaving the estate, this floor would disappear. An initial old purchase
cash constraint could remain, but the old policies, values and subsequent
proofs would require a new derivation. No such change has been verified.

Equal consumption delivery dates do not by themselves justify permitting
young consumption to use terminal sale receipts while forbidding old
consumption from doing so. The terminal-only young relaxation is coherent
if the model explicitly makes young title proceeds spendable at retirement
but protects the old occupied house for the estate. That is an
age-specific terminal restriction, not a consequence of merely saying
that old households take no new loans. Old households can still sell or
downsize at the beginning of their final stage; this restriction concerns
the house retained during that final stage.

There is a cleaner intermediate relaxation if a common payment order is
preferred. For both ages, pay consumption and property tax from liquid
receipts first; then liquidate titles, repay any existing young mortgage,
and enter retirement or distribute the old estate. The young must satisfy

\[
 c+\tau_tP_t h\le Y_t+k/q,
 \qquad b^+=a'_O+P_{t+1}h\ge0.
\]

For some admissible initial mortgage, the first condition is equivalent
to `a'_O>=-d/q`. Eliminating that mortgage gives the exact additional
inequality

\[
 qa'_O+\phi P_t h\ge0.
\]

Thus the complete projected owner conditions under this common ordering
are the consolidated budget, the original deposit cap, terminal solvency
`a'_O+P_{t+1}h>=0`, and `qa'_O+phi P_t h>=0`. Sufficiency follows by choosing

\[
 d=\max\{0,P_t h-w-T_t,-qa'_O\},\qquad
 k=w+T_t+d-P_t h.
\]

The two origination inequalities ensure `d<=phi P_t h`. Consumption and
tax can be paid before the sale, and terminal solvency then permits full
mortgage repayment from remaining cash and title proceeds. There is no
loan entering old age. For an old household with no initial loan, the
same consumption-before-liquidation order yields nonnegative financial
estate and the retained housing-estate floor.

This intermediate option strictly weakens `a'_O>=0` while preserving a
common payment order and the current old-age restriction. Its owner set
is contained in the terminal-solvent set already covered by renter
replication, so every previously verified all-young-rental result remains
unchanged. It is the cleaner recommendation if the packet seeks to remove
the unnecessary young asset bound without adding age-specific permission
to finance consumption from end-stage title proceeds. The unrestricted
terminal-only relaxation remains an admissible explicitly asymmetric
alternative, rather than an unqualified simplification.
