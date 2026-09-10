# Pro response — two periods with common positive retirement income

Captured September 9, 2026, after completion in https://chatgpt.com/c/6aa1e33d-d328-83ea-971e-0541c9e3e359. External advice; no automatic adoption into the model or slides. The following is the complete visible response, preserving original TeX and attachment labels.

<!-- BEGIN CAPTURE -->
## A two-period model with common positive retirement income

**Common positive retirement income does not obstruct the result.** The two types can lie on opposite sides of the saving decision: the poorer household would like to borrow against retirement income, while the richer household wants to accumulate additional retirement wealth. Restricted rentals then matter because ownership requires retaining housing equity—not merely paying the housing user cost.

The extra restriction relative to the previous result is substantive: **the richer type must finance its housing equity from wealth accumulated above its retirement income**. Its relevant resources are \(\beta Kx_H-qy^o\), not \(\beta Kx_H\).

The proposition below establishes a mixed-tenure equilibrium, the full dated consumption–housing allocation result, and both requested parental-fertility conclusions. It imposes no ordering between \(\beta\) and \(q\).

### Model

Write \(y=y^o>0\). Young endowments are \(0<w_L<w_H\), with fractions \(f\) and \(1-f\). Both types receive \(y\) when old. Preferences are exactly those in the clarification:

\[
u^y=\log(c^y-\chi n)+\alpha\log(h^y-\kappa n)+\vartheta\log n,
\qquad
u^o=\log c^o+\alpha\log h^o+\omega_B\log e,
\]

with lifetime utility \(u^y+\beta u^o\) and all preference and child-cost coefficients positive.

Maintain the stipulated serviced-interest convention:

\[
c+(1+\tau)Ph+a=w_i,\qquad
z=y+\frac aq+Ph,\qquad a\geq-q\phi Ph,
\tag{1}
\]

where \(0<q,\phi<1\) and \(\tau>0\). Equivalently,

\[
c+ph+qz=w_i+qy,\qquad
q(z-y)\geq\eta Ph,
\tag{2}
\]

with

\[
d=1+\tau-q,\qquad p=dP,\qquad
\eta=q(1-\phi),\qquad \Gamma=\frac{\eta}{d}>0.
\]

Here \(a\) is the net financial position after accounting for interest servicing. The retained-equity requirement is \(z-y\geq(1-\phi)Ph\). Retirement income is unavailable for that equity requirement while young, but is fully available upon retirement.

Renters satisfy \(c+ph+a=w_i\), \(a\geq0\), and \(h\leq r\). In retirement, the same owner budget and borrowing rule apply with current resources \(z\), subsequent payoff \(e\), and no subsequent income. Old households may sell and resize.

The stock \(\bar H\) is common and divisible; ownership has no size minimum, upper bound, taste premium, or adjustment cost. Competitive external finance gives \(p=dP\). For this calculation, **property taxes finance an additively valued public service**, with expenditure \(G=\tau P\bar H\), rather than household rebates. This remains a proposed fiscal closure. Estates do not determine entrants’ endowments.

For equilibrium closure, retain the previous stationary convention: a parent produces \(\nu n\) entering households; entrants receive the fixed type distribution independently of parental type. Both adult cohorts have endogenous mass \(N\). No transition calculation is needed.

## Proposition: housing allocation and parental fertility

Define

\[
K=1+\alpha+\omega_B,\qquad
A=\alpha+\vartheta,\qquad J=1+A,
\]

and the following quantities, all explicit functions of primitives:

\[
x=\frac{w_H+qy}{J+\beta K},\qquad
\ell=\frac{w_L}{J},\qquad
m=\frac{qy}{\beta K},\qquad
B=f\frac yK+(1-f)\frac{\beta x}{q}.
\tag{3}
\]

The interpretation is useful. \(x\) is the richer type’s unrestricted young adult consumption. \(\ell\) is the poorer type’s adult consumption with a freely sized rental and zero saving. At zero saving, \(m\) is the adult-consumption level at which the household is indifferent at the margin about saving.

Suppose \(\ell<m\), and define

\[
k=\frac{1+\Gamma}{1+\Gamma\ell/m}>1,
\qquad b=k\ell<m.
\tag{4}
\]

The quantity \(b\) will be a primitive upper bound on the poorer type’s equilibrium adult consumption. Using it avoids the stronger restriction based on its entire initial endowment.

To define the price interval, let \(n_L(p)\), for \(0<p<w_L/r\), be the unique solution

\[
\frac{\vartheta}{n_L}
=\frac{\chi}{w_L-pr-\chi n_L}
+\frac{\alpha\kappa}{r-\kappa n_L}.
\tag{5}
\]

These are initially **candidate** no-saving renter choices. Set

\[
x_L(p)=w_L-pr-\chi n_L(p),\qquad
s_L(p)=r-\kappa n_L(p),\qquad
M(p)=\frac{\alpha x_L(p)}{s_L(p)}.
\]

Define \(p_-,p_+\) by

\[
M(p_-)=kp_-,
\qquad M(p_+)=p_+.
\tag{6}
\]

Both are elementary primitive cutoffs; their polynomial is given below. They satisfy \(0<p_-<p_+<w_L/r\).

Finally, set

\[
n_H(p)=\frac{\vartheta x}{\chi+\kappa p},
\qquad F(p)=fn_L(p)+(1-f)n_H(p).
\]

**Proposition.** Suppose the following conditions hold:

\[
\boxed{
\omega_B\geq\alpha\Gamma,
\qquad
(\beta K-\Gamma A)x\geq qy.
}
\tag{Finance}
\]

\[
\boxed{
\ell<m,
\qquad
B\geq fb+(1-f)x.
}
\tag{Resources}
\]

\[
\boxed{
F(p_+)<\frac1\nu<F(p_-).
}
\tag{Market}
\]

Then a stationary competitive equilibrium exists with a unique price **within the verified regime**:

\[
F(p^*)=\frac1\nu,\qquad P^*=\frac{p^*}{d}.
\tag{7}
\]

The price is a unique scalar root, not a claimed closed-form expression.

The poorer young rent \(r\), save zero, and strictly prefer this allocation to every ownership alternative. The richer young own a home larger than \(r\). Their choices are

\[
\begin{aligned}
L:\quad&
h_L^y=r,\quad c_L^y=w_L-pr,\quad
n_L=n_L(p),\quad a_L=0,\quad z_L=y;\\
H:\quad&
c_H^y=x+\chi n_H,\quad
h_H^y=\frac{\alpha x}{p}+\kappa n_H,\quad
z_H=\frac{\beta Kx}{q},\\
& a_H=\beta Kx-qy-qPh_H^y.
\end{aligned}
\tag{8}
\]

For each old type,

\[
c_i^o=\frac{z_i}{K},\qquad
h_i^o=\frac{\alpha z_i}{Kp},\qquad
qe_i=\frac{\omega_Bz_i}{K}.
\tag{9}
\]

Housing clears at

\[
N=\frac{\bar H}{D(p)},\qquad
D(p)=fr+(1-f)h_H^y+\frac{\alpha B}{p}.
\tag{10}
\]

At this equilibrium, the full dated, equally weighted planner, initially holding fertility fixed, assigns **more aggregate housing to the young** and **more consumption and housing to every poorer young household**. Those poorer parents subsequently choose higher fertility at their assigned bundles. The joint dated planner choosing fertility as well as consumption and housing also chooses higher mean fertility.

Other equilibrium regimes are not excluded. Old households whose optimal home fits within the rental ceiling may be indifferent about tenure.

## 1. Proof of equilibrium, including all tenure deviations

### Retirement and the richer young

Without financing restrictions, the old optimum is (9). Its ownership constraint is precisely

\[
qe_i\geq\eta Ph_i^o
\quad\Longleftrightarrow\quad
\omega_B\geq\alpha\Gamma.
\]

Thus ownership implements the unrestricted optimum at **every** \(z>0\). No rental choice can improve on that unrestricted optimum. The old value function is consequently

\[
V^o(z)=K\log z+C(p).
\]

The richer young household’s unrestricted lifetime optimum gives (8). Its exact finance test is

\[
\begin{aligned}
q(z_H-y)
&=\beta Kx-qy\\
&\geq \Gamma x
\left[
\alpha+\vartheta\frac{\kappa p}{\chi+\kappa p}
\right].
\end{aligned}
\tag{11}
\]

The second Finance inequality is a price-independent sufficient bound for (11).

**The subtraction of \(qy\) is essential.** The same retirement income that supports lifetime consumption cannot also be counted as privately accumulated equity.

Since this unrestricted optimum is financeable, it dominates all rental and ownership deviations. Notice also that Finance implies

\[
x>m>b.
\tag{12}
\]

### The poorer young: a global comparison

For \(p\in(p_-,p_+)\), define

\[
t=\frac{M(p)}p\in(1,k).
\]

The candidate renter’s budget and fertility condition imply

\[
tw_L-Jx_L=(t-1)(x_L+\chi n_L)>0.
\]

Therefore

\[
x_L<t\ell<k\ell=b<m.
\tag{13}
\]

This establishes zero saving from primitives rather than assuming it.

Let

\[
\xi=\frac{\alpha}{s_L}-\frac p{x_L}>0,
\qquad
\zeta=\frac1{x_L}-\frac1m>0.
\]

Write \(\sigma=1-b/m>0\). The definitions give

\[
k-1=\Gamma\sigma,
\qquad
\zeta\geq\frac{\sigma}{x_L},
\qquad
\xi<\eta P\zeta.
\tag{14}
\]

Consider **any** alternative lifetime plan \((c,h,n,z)\). Concavity of

\[
\log(c-\chi n)+\alpha\log(h-\kappa n)
+\vartheta\log n+\beta K\log z
\]

and the consolidated budget imply

\[
U-U_L\leq \xi(h-r)-q\zeta(z-y).
\tag{15}
\]

The fertility term disappears because the candidate satisfies its fertility first-order condition; the inequality therefore covers different fertility as well as different consumption, housing, and saving.

For another rental plan, \(h\leq r\) and \(z\geq y\), so (15) is nonpositive. For any owner plan,

\[
q(z-y)\geq\eta Ph,
\]

and hence

\[
U_O-U_L
\leq(\xi-\eta P\zeta)h-\xi r<0.
\tag{16}
\]

This proves global tenure preference.

Moreover,

\[
n_L=\frac{\vartheta x_L}{\chi+t\kappa p}
<
\frac{\vartheta x}{\chi+\kappa p}=n_H,
\qquad
s_L<\frac{\alpha x}{p}.
\]

Thus the richer young household’s unrestricted home is indeed larger than \(r\).

### Existence and the justified uniqueness

Differentiating (5) gives \(n_L'(p)<0\), \(x_L'(p)<0\), and \(s_L'(p)>0\). Consequently \(M\) is strictly decreasing, while

\[
M(0)>0,\qquad \lim_{p\uparrow w_L/r}M(p)=0.
\]

The two cutoff roots exist and satisfy \(p_-<p_+\).

Both fertility demands are strictly decreasing, so the Market inequalities give exactly one root of (7) inside this interval. The household proof applies throughout that interval, and (10) supplies a positive market-clearing cohort mass. Old wealth is precisely the wealth generated by the preceding cohort’s choices. This completes the stationary equilibrium construction.

## 2. The full dated consumption–housing planner

The comparison includes **all households currently alive**. It preserves each young household’s \(z_i\), each old household’s net estate \(e_i\), and continuation prices. The authority may relax private finance and reassign tenure within the common stock. This is the dated benchmark specified in the packet, not a stationary lifetime-welfare comparison. oracle_essential_theory_packet

Normalize quantities by cohort mass and write

\[
\bar x=fx_L+(1-f)x,\qquad
\bar s=fs_L+(1-f)\frac{\alpha x}{p},\qquad
\bar h_o=\frac{\alpha B}{p}.
\]

At fixed competitive fertility, current resources are

\[
\mathcal C=\bar x+B+\chi\bar n,\qquad
\mathcal H=\bar s+\bar h_o+\kappa\bar n.
\]

The planner maximizes the sum of current utilities over all four age–type groups, choosing their consumption and housing subject to

\[
\sum_i f_i(c_i^y+c_i^o)=\mathcal C,\qquad
\sum_i f_i(h_i^y+h_i^o)=\mathcal H,
\]

where \(f_L=f\), \(f_H=1-f\). Fixed continuation utility and estate utility are constants in this problem.

The unique real allocation equalizes adult consumption and adult space:

\[
X=\frac{\bar x+B}{2},\qquad
S=\frac{\bar s+\bar h_o}{2},
\tag{17}
\]

\[
c_i^{y,F}=X+\chi n_i,\quad h_i^{y,F}=S+\kappa n_i,
\qquad c_i^{o,F}=X,\quad h_i^{o,F}=S.
\tag{18}
\]

By Resources and (13),

\[
B\geq fb+(1-f)x>\bar x,
\]

while the binding rental ceiling implies

\[
\bar s<\frac{\alpha\bar x}{p}<\frac{\alpha B}{p}=\bar h_o.
\]

Therefore

\[
\boxed{
\bar h_y^F-\bar h_y^E
=S-\bar s
=\frac{\bar h_o-\bar s}{2}>0.
}
\tag{19}
\]

Also \(\bar x>x_L\) and \(\bar s>s_L\), so

\[
X>x_L,\qquad S>s_L.
\]

Every poorer young household receives strictly more total consumption and housing.

**Settlement.** Each household receives current transfer

\[
T_i=\Delta c_i+p\Delta h_i,\qquad \sum_iT_i=0.
\]

Changes in owned titles are offset by

\[
\Delta a_i=-qP\,\Delta h_i^{\mathrm{owned}},
\]

including intermediaries’ corresponding positions. Continuation wealth, estates, public expenditure, and aggregate financial claims remain unchanged. Relaxing finance to a 100% financed balance suffices here: all young continuation targets satisfy \(z_i\geq y\), and estates remain positive.

### What is specifically inefficient?

The Resources condition signs the **equal-weight redistribution**. It is not needed for the underlying efficiency gap.

Throughout the verified tenure regime,

\[
\frac{\alpha x_L}{s_L}>p
=\frac{\alpha c_i^o}{h_i^o}.
\tag{20}
\]

A small transfer of housing \(\epsilon\) from an old household to a poorer family can compensate the old exactly with

\[
\Delta c_o
=c_o\left[
\left(\frac{h_o}{h_o-\epsilon}\right)^\alpha-1
\right]
=p\epsilon+O(\epsilon^2).
\]

Taking these goods from the family leaves it a gain

\[
\left(\frac{\alpha}{s_L}-\frac p{x_L}\right)\epsilon
+O(\epsilon^2)>0.
\]

That is a dated Pareto improvement after relaxing finance, independent of equal welfare weights. Restricted rentals prevent purchasing additional *housing services* directly; the retained-equity requirement prevents reproducing that rental purchase through ownership. Removing either obstruction eliminates this particular marginal housing wedge. The full utilitarian allocation need not choose the compensating transfer.

## 3. The parental-fertility conclusions

### Poorer parents offered their new bundles

Hold the assigned **total** consumption and housing in (18) fixed while the parent rechooses fertility. At its original \(n_L\), the fertility derivative is

\[
\frac{\vartheta}{n_L}-\frac{\chi}{X}
-\frac{\alpha\kappa}{S}
>
\frac{\vartheta}{n_L}-\frac{\chi}{x_L}
-\frac{\alpha\kappa}{s_L}=0.
\]

Strict concavity implies

\[
\boxed{n_L^{\mathrm{new}}>n_L.}
\]

This uses both the consumption gain and the housing gain. There is no small-child-goods-cost assumption.

### The joint dated, parents-only planner

Now let the planner choose fertility too, still preserving parental continuation opportunities and old estates. Identical parental preferences imply common chosen fertility \(n^F\), with

\[
X(n)=\frac{\mathcal C-\chi n}{2},
\qquad
S(n)=\frac{\mathcal H-\kappa n}{2}.
\]

Its unique interior solution satisfies

\[
\frac{\vartheta}{n^F}
=\frac{\chi}{X(n^F)}
+\frac{\alpha\kappa}{S(n^F)}.
\tag{21}
\]

Private fertility obeys

\[
n_i=\vartheta
\left(\frac{\chi}{x_i}+\frac{\alpha\kappa}{s_i}\right)^{-1}.
\]

The weighted harmonic mean is concave. Hence

\[
\frac{\vartheta}{\bar n}
\geq\frac{\chi}{\bar x}+\frac{\alpha\kappa}{\bar s}
>
\frac{\chi}{X(\bar n)}+\frac{\alpha\kappa}{S(\bar n)}.
\]

The planner’s fertility derivative is strictly decreasing, establishing

\[
\boxed{n^F>\bar n.}
\tag{22}
\]

The objective values existing parents; no independent welfare weight on unborn people has been added. This is a dated parental optimum, not a new stationary population.

**Outside the sufficient resource condition**, a housing gain alone is insufficient. For an affected parent, the exact fertility test is

\[
\chi\left(\frac1{x_L}-\frac1X\right)
+\alpha\kappa\left(\frac1{s_L}-\frac1S\right)>0.
\tag{23}
\]

If consumption falls, the space benefit must exceed the goods loss in this expression. For the joint planner, the exact mean-fertility test is

\[
\frac{\vartheta}{\bar n}>
\frac{\chi}{X(\bar n)}+\frac{\alpha\kappa}{S(\bar n)}.
\]

Neither follows merely from a welfare improvement.

## 4. The conditions have nonempty interior—and identifiable costs

The key financial and distributional restrictions are not mutually contradictory.

Fix preference and financing parameters satisfying

\[
\omega_B>\alpha\Gamma,\qquad
\beta K>\Gamma A.
\]

For any \(y>0\), write \(\rho=\beta/q\) and \(m=qy/(\beta K)\). Consider the open endowment ranges

\[
0<w_L<
\frac{Jm\,\min\{1,\rho\}}
{1+\Gamma(1-\rho)_+},
\tag{24}
\]

\[
w_H>
\frac{qy(J+\Gamma A)}{\beta K-\Gamma A}.
\tag{25}
\]

They imply \(\ell<m\), \(b<y/K\), the strict richer-type finance condition, and \(w_H>w_L\).

Next choose any

\[
\frac{(1-\rho)_+x}
{\,y/K-b+(1-\rho)_+x\,}<f<1.
\tag{26}
\]

The denominator exceeds the numerator, so this interval is nonempty, and Resources holds strictly.

For arbitrary positive \(\chi,\kappa,r\), the two price cutoffs are distinct and \(F\) is strictly decreasing. Therefore the permissible replacement level has the nonempty open interval

\[
F(p_+)<1/\nu<F(p_-).
\]

Any positive housing stock then determines a positive \(N\). Thus the sufficient region has nonempty interior in primitive parameter space; it is not tied to a reference calibration or an equality-calibrated fertility preference.

The restrictions nevertheless have economic content:

The old financing inequality ensures that their unrestricted housing–estate mix is feasible. The richer-type inequality is a conservative, price-independent version of the exact test (11). Neither is necessary for all possible equilibrium regimes.

The distributional restriction is what permits the full planner to move **both goods and space** toward young households. When \(\beta<q\), it reads

\[
f\left(\frac yK-b\right)
\geq(1-f)\left(1-\frac\beta q\right)x.
\tag{27}
\]

The poorer households’ rise in available resources at retirement must offset the richer households’ declining consumption profile. A very large rich-young endowment or a small poorer share can violate it. That restriction should remain visible; borrowing constraints alone do not determine the direction of an equally weighted age redistribution.

## 5. Brief relation to Coven–Golder–Gupta–Ndiaye

Their analytical model treats housing as a dividend-paying asset, excludes housing services from its stated utility, and has the old sell the stock; its capitalization lemma and intergenerational-welfare proposition are different from the joint consumption–housing allocation comparison here (pp. 15–18, equations (4)–(16)). [AbdouEcon+1](https://abdouecon.github.io/research/papers/Property_Tax.pdf)

Their quantitative model includes recurring earnings and pensions, an owner-size minimum, segmented tenure stocks, LTV/PTI constraints, adjustment costs, and flexible repayment after origination (pp. 20–23). [AbdouEcon+1](https://abdouecon.github.io/research/papers/Property_Tax.pdf) Our illustration shares the role of finance in tenure access but retains the stipulated common stock, rental ceiling, and serviced-interest constraint. Its additional result is the dated allocation theorem and parental-fertility implication—not a property-tax implementation theorem.

**The common retirement payment changes the thresholds, not the feasibility of the mechanism.** What is established is one mixed-tenure equilibrium regime with a genuine financial/rental efficiency gap and a separately stated resource condition delivering the full planner and fertility signs.

### Explicit cutoff calculation

For direct verification, the price satisfying \(M(p)=tp\) is the unique positive root of

\[
\kappa tr(t+\alpha+\vartheta)p^2
+\Big[
\chi r\{\alpha+t(1+\vartheta)\}
-\kappa t(\alpha+\vartheta)w_L
\Big]p
-\alpha\chi w_L=0.
\]

Use \(t=k\) for \(p_-\) and \(t=1\) for \(p_+\).

Given \(p\), let \(v=w_L-pr\). The poorer type’s fertility is the smaller feasible root of

\[
\chi\kappa J n^2
-\big[(\vartheta+\alpha)\kappa v
       +(\vartheta+1)\chi r\big]n
+\vartheta rv=0.
\]

Thus household demands and the regime cutoffs are elementary; only the equilibrium price in (7) requires the stated monotone scalar root.
<!-- END CAPTURE -->
