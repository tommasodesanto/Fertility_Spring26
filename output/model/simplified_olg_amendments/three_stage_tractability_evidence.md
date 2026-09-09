# Three-stage tractability test: checked supporting calculations

September 9, 2026. This preserves a bounded model experiment, not a replacement for the current theory or the September 14 deck. The reader summary is `output/pdf/three_stage_tractability_test.pdf`. Four independent scopes cover household algebra, finance, the full dated planner, and stationary existence. A second independent check of the stationary construction passed. No numerical model runs or numerical-reference-point existence proofs were used.

The established result is analytical existence of a deliberately strong heterogeneous stationary family at any fixed positive LTV below one, a small-positive-tax extension, and a full fixed-fertility planner housing gain. Further finite resource bounds below yield a positive conditional private fertility response. No global uniqueness, policy implementation, joint-fertility planner theorem, or equilibrium transition theorem is claimed. The full equilibrium price has not been solved in closed form.


---

# Three-stage tractability test — not an adopted model

Time limit: 17:22–18:07 UTC, September 9. User requests testing the proposed extra working-age stage before any rewrite. Deliver explicit choices, an economically interpretable regime condition and welfare/fertility implication if feasible; otherwise precise algebraic obstruction. No numerical-reference-point existence proof. Current main deliverable is the September 14 deck, not another full note.

Reference source: latex/JMP_DS_suggestions/simplified_olg_consolidated_theory.tex; existing reviewer assessment: output/model/simplified_olg_amendments/README.md. Keep original notation q, beta, alpha, gamma, vartheta, chi, kappa, phi, tau^p. Introduce superscript m only for mature workers. Do not edit the author-controlled manuscript.

## Primary candidate, to test honestly

Three equal-length economic stages: young family y, mature working family m, old household o. q in (0,1) is the one-stage bond price; beta discounts one stage. Retain income/entry-wealth heterogeneity with compact positive support (w=y_y+b, y_m, y_o), independent logistic owner taste, fixed lifetime tenure and caps H_R<H_O, positive property tax tau rebated equally at every stage, free housing resizing. All current income at each date is available to pay bills and purchase housing.

Young utility remains u_y=log(c_y-chi n)+alpha log(h_y-kappa n)+vartheta log n. The primary mature-family utility is u_m=log(c_m-chi n)+alpha log(h_m-kappa n): children still use resources, fertility was chosen when young. Old utility unchanged log c_o+gamma log h_o+omega_B log e. Lifetime objective u_y+beta u_m+beta^2 u_o, plus ownership taste. This child timing is deliberate; if simplifying to needs only when young is analytically consequential, identify it explicitly as a DIFFERENT, unadopted variant.

At stationarity let p=(1-q+q tau)P, D=1-q+q tau. Budgets below use the same date-end rent/tax timing as the maintained model.

Young renter: c_y+q a_m+p h_y=w+T, a_m>=0.
Young owner: c_y+q a_m+(1+q tau)P h_y=w+T, q a_m+phi P h_y>=0.

Mature renter: c_m+q a_o+p h_m=a_m+y_m+T, a_o>=0.
Mature owner: c_m+q a_o+(1+q tau)P h_m=a_m+P h_y+y_m+T, a_o>=0.

Old renter: c_o+q e+p h_o=a_o+y_o+T.
Old owner: c_o+q e+p h_o=a_o+P h_m+y_o+T, e>=P h_o.

All housing obeys retained tenure cap. Log domains positive. The young mortgage is repaid on entry to mature work; mature net financial saving cannot be negative, so retirement house equity is not encumbered by a new mortgage. This is an explicit institutional restriction of the test, not an innocuous relabeling. Do not hide it or claim realistic amortization from nothing. Check whether allowing mature refinance changes any claim if you can do so briefly.

At a positive SS, each stage has mass N and mean fertility1/nu; housing N(hbar_y+hbar_m+hbar_o)=Hstock; fiscal3NT=q tau P Hstock. Entry endowments remain exogenous; estates do not generate them.

First derive a candidate regime with young owners financially constrained, mature households saving, old financial-estate floors slack, and selected housing caps slack. CONDITIONAL formulas are only step1; say what inequalities verify them and whether rental caps can be active. Retain taxes throughout or clearly label a tractability specialization. Do not introduce phi=q automatically. In particular test whether origination phi=.8 can be compatible analytically, without pretending one numerical equilibrium proves the result.

The fixed-fertility dated planner includes all THREE current ages, chooses all current consumption/housing, preserves tenure, inherited claims, future real opportunities and net estates, and may relax private financial lower bounds. Each current adult household receives weight1 on remaining utility. Preserve original settlement logic. A y-o-only transfer holding m fixed may be a diagnostic direction, but is not the full planner optimum. For joint fertility, current m fertility is predetermined; changing young n changes their future mature-family needs and cannot be assumed costless or ignored in a fixed-future-opportunities comparison.

Desired standard: keep the original mechanism (young family needs, borrowing, small rentals, older housing) while seeking a short explicit benchmark. Do not broaden to general multi-age existence, a new tax transition polynomial, calibration, numerical solver or a different welfare criterion. Lead combines work; mathematical tasks are read-only.


---

# Conditional three-stage household derivation

Working algebra only, September 9. This does not establish equilibrium existence or adopt the candidate model. Child needs remain present in both young and mature stages. The branch has binding young-owner finance, strictly positive mature financial saving, and slack housing caps and old estate floor.

Let
\[
W=w+T,\quad M=y_m+T,\quad V=y_o+T,\quad U=M+qV,
\]
\[
D=1-q+q\tau,\quad p=DP,\quad
L=(1-\phi+q\tau)P,\quad
\delta=(1-\phi/q)P=(L-p)/q.
\]
Write \(x=c_y-\chi n,\ s=h_y-\kappa n\), and
\[
t=\chi+\kappa p,\quad t_L=\chi+\kappa L,\quad
d=t-\kappa\delta,\quad
K=1+\gamma+\omega_B,\quad G=1+\alpha+\beta K,
\quad C=\beta G,\quad E=1+\alpha+\vartheta.
\]
Here \(C\) is a continuation-utility coefficient, not aggregate consumption.

## Mature continuation, with children

Binding young finance gives \(a_m=-(\phi/q)Ph_y\) and \(c_y+Lh_y=W\). The mature surplus available after child needs is
\[
B=U+\delta h_y-tn=U+\delta s-dn.
\]
On the stated mature/old branch,
\[
x_m=B/G,\quad s_m=\alpha x_m/p,\quad
c_o=\beta x_m/q,\quad h_o=\gamma\beta x_m/(qp),
\quad e=\omega_B\beta x_m/q^2.
\]
The young reduced objective is
\[
\log x+\alpha\log s+\vartheta\log n+C\log B+\text{constant},
\qquad x+Ls+t_Ln=W.
\]

## One scalar cubic

Define \(r=\beta x/x_m=Cx/B\). The first-order conditions give
\[
s=\frac{\alpha x}{L-\delta r},\qquad
n=\frac{\vartheta x}{t_L+dr}.
\]
Put \(R=U/W\) and
\[
F(r)=1+\frac{\alpha L}{L-\delta r}
+\frac{\vartheta t_L}{t_L+dr}.
\]
Then
\[
(1+Rr)F(r)=E+C,\qquad
x=\frac{W+Ur}{E+C}.
\]
Equivalently, the polynomial is
\[
(1+Rr)\big[(L-\delta r)(t_L+dr)
+\alpha L(t_L+dr)+\vartheta t_L(L-\delta r)\big]
-(E+C)(L-\delta r)(t_L+dr)=0.
\]
It is generically cubic, with leading coefficient \(-R\delta d\).
There is no generic quadratic reduction. Cardano is available but unnecessary
for sign comparisons. Special loci, including \(\phi=q\), reduce the degree.
At \(\phi=q\), the reduced quadratic is
\[
(1+\alpha)Rr^2+[ER-\vartheta-C]r-C=0.
\]

The exact young mortgage multiplier is
\[
\mu_y=\frac{1-r/q}{x}.
\]
Thus strict borrowing requires \(0<r<q\). Both denominators above are
automatically positive on this interval: their endpoint values are
\(L,p\) and \(t_L,(1+q)t\).

## Unique root and resource tests without radicals

Define
\[
\mathcal R(r)=\frac{(E+C)/F(r)-1}{r}.
\]
It is strictly decreasing wherever both denominators are positive. Indeed,
with \(z_1=L/(L-\delta r)\), \(z_2=t_L/(t_L+dr)\),
\[
F+rF'=1+\alpha z_1^2+\vartheta z_2^2,\qquad
F^2\le E(F+rF')<(E+C)(F+rF').
\]
Hence \(\mathcal R'<0\). Moreover \(\mathcal R(r)\to\infty\) as
\(r\downarrow0\). The unique root lies below \(q\) exactly when
\[
U/W>\mathcal R(q).
\]
For any threshold \(0<r_* \le q\),
\[
r<r_*\quad\Longleftrightarrow\quad U/W>\mathcal R(r_*).
\]
These statements solve the relaxed continuation branch; the remaining
constraints must also be checked.

## Exact regime checks and age comparisons

Mature financial saving is
\[
a_o=\left(\frac{\beta K}{q}-\frac{\alpha}{D}\right)x_m
-V-\kappa Pn>0.
\]
Its coefficient must therefore be positive. Old estate slackness is
\(\omega_BD>q\gamma\). Check all three owner homes below \(H_O\):
\[
s+\kappa n<H_O,\qquad
\alpha x_m/p+\kappa n<H_O,\qquad
\gamma\beta x_m/(qp)<H_O.
\]
Together with \(0<r<q\), these conditions verify the original budgets and
all relevant complementary-slackness conditions. Relaxing mature borrowing
would leave this optimum unchanged, because \(a_o>0\) already.

The age comparisons are especially short in \(k=x_m/x=\beta/r\):
\[
c_o>x\iff k>q/\beta,
\qquad
h_o>s\iff
k>\frac{\alpha qp}{\gamma\beta L}+\frac{\beta\delta}{L}.
\]
These are conditional resource comparisons, not equilibrium conclusions.

## Fixed fertility gives a quadratic

For fixed \(n\), let \(A_y=W-t_Ln\), \(A_m=U-dn\).
The young problem is
\[
\max_s\ \log(A_y-Ls)+\alpha\log s+C\log(A_m+\delta s).
\]
Its solution is the unique feasible root of
\[
\alpha A_yA_m+
[\delta A_y(\alpha+C)-(1+\alpha)LA_m]s
-\delta L(1+\alpha+C)s^2=0.
\]
Feasibility means \(s>0,\ A_y-Ls>0,\ A_m+\delta s>0\).
This quadratic is valid at general \(\phi\), including \(\delta<0\).
At \(\delta=0\), \(s=\alpha A_y/[(1+\alpha)L]\).
Do not choose a quadratic sign without checking that feasible interval.

## Full three-age fixed-fertility housing comparison

With all planner caps slack, an aggregate young gain requires
\(\alpha(\bar s_m+\bar h_o)>(\alpha+\gamma)\bar s_y\).
A sufficient pointwise condition is \(r<r_H\), where
\[
A_H=\beta(\alpha+\gamma\beta/q),\quad
\ell=L/p,\quad g=\delta/p=(\ell-1)/q,\quad
r_H=\frac{A_H\ell}{\alpha+\gamma+A_Hg}.
\]
For \(\beta\le q\), the denominator is positive and \(r_H\le q\), with
equality only at \(\beta=q\). This threshold is price-independent.
The exact conditional income test is \(U/W>\mathcal R(r_H)\).

A conservative bound eliminates the remaining price dependence. For
\(\epsilon=\kappa p/(\chi+\kappa p)\in(0,1)\),
\[
\frac{t_L}{t_L+dr}
=\frac{1+(\ell-1)\epsilon}
{1+r+(\ell-1)(1-r/q)\epsilon}.
\]
This is monotone in \(\epsilon\), so set
\[
m(r)=\min\left\{\frac1{1+r},
\frac{\ell}{\ell+r(1-g)}\right\},\quad
F_{\min}(r)=1+\frac{\alpha\ell}{\ell-gr}+\vartheta m(r).
\]
Then
\[
\mathcal R(r)\le
\mathcal R_{\max}(r)=\frac{(E+C)/F_{\min}(r)-1}{r}.
\]
Thus \(U/W>\mathcal R_{\max}(r_H)\) is sufficient uniformly in price.
It still requires the saving, estate, and competitive/planner cap checks.
Since \(U/W=[y_m+qy_o+(1+q)T]/(w+T)\), a pure pretax-endowment
certificate may additionally need a bound on the endogenous rebate.

No active-rental-cap formula, stationary existence result, or joint-fertility
planner claim is established here.


---

# Three-stage finance audit

September 9, 2026. Bounded analytical check of `benchmark.md`; no model run,
paper/slide edit, or change to the proposed budgets. The household-root
calculation was coordinated with `general_preferences`. The results below
are conditional household results, not a general-equilibrium existence
theorem.

## Verdict

The repayment accounting is correct. An intervening earning period can
support an 80 percent origination LTV with strictly constrained young
owners, strictly positive mature financial saving, and old housing above
young housing. It does not require forgiving the initial mortgage.

The proposed mature restriction \(a_o\ge0\) is a substantive additional
credit restriction: mature workers cannot carry negative net financial
wealth into retirement, even though they still work and own collateral.
However, it is strictly slack in the certificate below. Allowing mature
refinancing would therefore leave those conditional choices unchanged.

High mature income alone does not guarantee mature saving. A necessary
preference/price condition is
\[
\boxed{\beta K D>\alpha q},\qquad
K=1+\gamma+\omega_B,\quad D=1-q+q\tau.
\]
That condition must be stated and checked. The old estate condition
\(\omega_BD>q\gamma\) is separate.

## 1. What is actually repaid

Write resources inclusive of the common rebate as
\(W=w+T\), \(M=y_m+T\), \(V=y_o+T\). Define
\[
L=(1-\phi+q\tau)P,\quad p=DP,\quad
\delta=(1-\phi/q)P=(L-p)/q.
\]
When the young mortgage binds,
\[
a_m=-\phi Ph_y/q,
\qquad Z_m=M+a_m+Ph_y=M+\delta h_y.
\tag{1}
\]
The face repayment is the original principal \(\phi Ph_y\) divided by
\(q\). The housing title is worth \(Ph_y\) on entry to maturity. If
\(\phi>q\), the difference is negative. For example, at \(\phi=.8\),
\(q=.5\), the mature household must cover \(.6Ph_y\) in addition to its
current spending; its initial mortgage repayment is \(1.6Ph_y\).
Equation (1) retains that obligation in full.

Thus mature earnings can fill the original liquidity gap before retirement.
The added stage does not manufacture equity. It changes the income timing:
repayment can be covered by peak working-age earnings rather than by the
retirement-income component alone. A restriction comparing only \(V\) with
\(W\) consequently need not imply the original theorem's small-LTV bound.

The mature constraint is a restriction on *net* financial wealth. With
identically priced loans and bonds, gross borrowing offset by financial
assets is not separately identified. A nonnegative net position admits a
no-mortgage representation; it should not be described as identified gross
amortization without that qualification.

## 2. Exact continuation and financial-regime tests

Retain mature child costs. Put
\[
t=\chi+\kappa p,\quad A=1+\alpha,\quad G=A+\beta K,
\qquad B=M+qV+\delta h_y-tn.
\]
Combining the mature and discounted old budgets leaves \(B\) for mature
adult bundles and old consumption, housing, and estates. When mature saving
and all relevant size/old-estate restrictions are slack, the exact choices
are
\[
x_m=B/G,\quad s_m=\alpha x_m/p,\quad h_m=s_m+\kappa n,
\]
\[
c_o=\beta x_m/q,\quad h_o=\gamma\beta x_m/(qp),\quad
e=\omega_B\beta x_m/q^2.
\tag{2}
\]
The continuation value is \(G\log B\) plus terms independent of current
\(h_y,n\). In particular, mature children cost \(tn\); dropping this term
would change the proposed model.

Define
\[
\Psi=\beta K/q-\alpha/D.
\]
The exact mature financial position is
\[
\boxed{a_o=\Psi x_m-V-\kappa Pn.}                               \tag{3}
\]
Because \(V,n>0\), \(\Psi>0\) is necessary for this saving regime.
If it fails, arbitrarily high mature income does not validate the
unconstrained continuation: higher income also raises desired mature
housing and leaves its purchase too costly relative to desired old-age
resources. The nonnegative-financial-wealth constraint would bind instead.

When \(\Psi>0\), the exact checks are
\[
\begin{array}{ll}
\text{young mortgage strictly binding:}&x_m>\beta x_y/q,\\
\text{mature financial saving positive:}&\Psi x_m>V+\kappa Pn,\\
\text{old housing above young:}&\gamma\beta x_m/(qp)>h_y.
\end{array}                                                       \tag{4}
\]
The first follows from the original young mortgage multiplier
\(\mu_y=1/x_y-\beta/(qx_m)\). These conditions are mutually compatible.
Once (3) is strictly positive, replacing \(a_o\ge0\) with a mature
collateral borrowing limit leaves the solution unchanged: it already
solves the unconstrained continuation and satisfies the weaker restriction.

There is no corresponding implication that old housing exceeds mature
housing. Indeed \(h_o/s_m=\beta\gamma/(q\alpha)\). With equal adult
housing weights and \(\beta/q\le1\), old housing is below mature total
housing because mature children still require space.

## 3. An exact income-ratio certificate

The household-root worker verified the following reduction, including the
future mature child costs. Write
\[
r=\beta x_y/x_m,\quad t_L=\chi+\kappa L,\quad d=t-\kappa\delta,
\quad E=1+\alpha+\vartheta,
\]
\[
F(r)=1+\frac{\alpha L}{L-\delta r}
       +\frac{\vartheta t_L}{t_L+dr},\qquad
\mathcal R(r)=\frac{(E+\beta G)/F(r)-1}{r}.
\tag{5}
\]
The young solution has
\[
x_y=W/F(r),\quad s_y=\alpha x_y/(L-\delta r),\quad
n=\vartheta x_y/(t_L+dr),\qquad
(M+qV)/W=\mathcal R(r).
\tag{6}
\]
For \(\phi>q\), both denominators are positive for \(r>0\), and
\(\mathcal R\) is strictly decreasing. This follows from
\[
F+rF'=1+\alpha\left(\frac L{L-\delta r}\right)^2
 +\vartheta\left(\frac{t_L}{t_L+dr}\right)^2,
\qquad F^2\le E(F+rF').
\]
Consequently bounds on \(r\) translate into explicit inequalities on
mature versus young resources, rather than an assumed multiplier sign.

For a general conservative certificate with \(\phi>q\), set
\(d_c=1-\phi+q\tau\), \(\ell=L/p=d_c/D\), \(\rho_o=V/W\). Since
\(F(r)\le E\), it suffices to have
\[
r<\min\left\{q,
\frac{\beta\Psi}{E\rho_o+\vartheta/d_c},
\frac{\gamma\beta^2\ell}{q(\alpha+\vartheta)}\right\}.
\tag{7}
\]
The exact income test is \((M+qV)/W>\mathcal R(\bar r)\), where
\(\bar r\) is the right side of (7). It can be conservative; the following
family verifies the exact tests with a less demanding earnings profile.

## 4. An 80 percent LTV family with mature resources between 2.5 and 2.8 times young resources

Take
\[
q=1/2,\quad \phi=4/5,\quad \tau=1/10,\quad
\alpha=\gamma=\vartheta=\kappa=1,\quad\omega_B=2,
\]
\[
49/100\le\beta\le1/2,\qquad
z:=\chi/(\kappa p)\ge10,\qquad 19/20\le V/W\le21/20.
\tag{8}
\]
Here \(D=11/20\), \(d_c=1/4\), \(\ell=5/11\), and
\(\delta/p=-12/11\). Impose the open income-ratio interval
\[
\boxed{\quad
\mathcal R(2/5)<\frac{M+qV}{W}<\mathcal R(39/100).
\quad}                                                          \tag{9}
\]
It selects \(39/100<r<2/5\) by strict monotonicity. The interval can be
imposed type by type with heterogeneous \(W,V/W,M/W\); it does not require
a degenerate or proportional-income distribution.

At \(r=2/5\),
\[
\frac{L-\delta r}{p}=49/55,\qquad
\frac{t_L+dr}{\kappa p}=(77z+71)/55,
\]
\[
F(2/5)=1+25/49+(55z+25)/(77z+71).
\]
The following bounds verify the financial and age restrictions throughout
the interval (9), rather than at a numerical solution:

* \(\mu_y x_y=1-r/q>1/5\).
* The old-versus-young housing comparison, after multiplying by \(r\),
  is hardest at the largest \(r\). At that endpoint,
  \[
  \frac{p h_o}{x_y}\ge\frac{2401}{2000}
  >\frac{55}{49}+\frac{55}{841}
  \ge\frac{p h_y}{x_y}.
  \]
  Thus \(h_o>h_y\) throughout (9).
* \(F(r)<9/4\) throughout the interval. Moreover
  \[
  \frac{\beta\Psi}{r}-\frac{\kappa Pn}{x_y}
  \ge\frac{14161}{5500}-\frac{100}{841}>0.
  \]
  To see this lower bound for smaller \(r\), multiply the expression by
  \(r\): the subtracted term is proportional to \(r/(t_L+dr)\), which
  increases with \(r\). Equation (3) therefore gives
  \[
  \frac{a_o}{W}>
  \frac{4}{9}\left(\frac{14161}{5500}-\frac{100}{841}\right)
  -\frac{21}{20}
  =\frac{575543}{13876500}>0.
  \]

Exact endpoint bounds in (9) also give
\[
\frac52<\frac MW<\frac{14}{5}.                                 \tag{10}
\]
For example, the lower bound exceeds \(5/2\) by at least
\(1487/54500\), and the upper bound is below \(14/5\) by at least
\(33632551/1263192840\). These are rational interval calculations, not
numerical model evaluations.

The old estate floor is strictly slack because
\(\omega_BD=11/10>q\gamma=1/2\). Given a positive trial price, finite
owner caps such as \(H_O>2W_{\max}/p\) make all three owner housing
choices slack in this family. Rental caps and rental policies remain a
separate check.

## Interpretation and unresolved scope

This supplies a conditional high-LTV financial certificate without an
extreme tenfold mature-income assumption. It does not establish empirical
plausibility: \(W\) includes entry wealth, so mature income of 2.5–2.8
times \(W\) can imply a still larger ratio to young labor earnings. The
child-cost restriction \(\chi\ge10\kappa p\) is also substantive: goods
costs dominate child-space expenditure in this demonstration. The estate
weight and discount factor must satisfy the saving restriction; they cannot
be ignored by asserting that mature earnings are high.

Conditions (8)–(10) concern resources inclusive of the endogenous rebate
and the equilibrium-relative child-cost ratio. The general-equilibrium
worker must establish price/rebate closure and verify these inequalities
at the resulting stationary equilibrium before they can become a primitive
equilibrium theorem. Finite caps must also be chosen consistently with
those price bounds. No property-tax transition, Pareto claim, or full
three-age planner/fertility result is established here.

The financial implication is nevertheless precise: interim earnings can
repay the entire young mortgage, finance the mature home, and leave positive
financial saving. The strict-saving certificate makes the added mature
credit ban locally inessential. It avoids the original theorem's tight
origination-LTV bound by changing the earnings timeline, not by forgiving
debt or assuming that old households cannot resize.


---

# Three-stage dated planner and fertility audit

September 9, 2026. Bounded analytical review of `benchmark.md`; no adoption, model run, source edit, or equilibrium-existence claim. Equal remaining-utility weights apply to all three currently living ages.

## 1. Full fixed-fertility optimum

At the stationary reference, let \(Q\) be the common distribution of type and retained tenure, \(H_i\) its cap, and \(n_i\) the fertility of each matched young/mature type. The current mature cohort's fertility was chosen previously; stationarity supplies the matching distribution, not a new mature fertility choice. Define
\[
x_i^y=c_i^y-\chi n_i,\quad x_i^m=c_i^m-\chi n_i,
\qquad s_i^y=h_i^y-\kappa n_i,\quad s_i^m=h_i^m-\kappa n_i.
\]
Bars are means under \(Q\), and \(\bar n=1/\nu\). Current goods and housing per cohort mass are
\[
\mathcal C=\bar x^y+\bar x^m+\bar c^o+2\chi\bar n,
\qquad
\mathcal H=\bar s^y+\bar s^m+\bar h^o+2\kappa\bar n=\bar H/N.
\]
With fertility, future real allocations and old net estates fixed, the planner maximizes
\[
\int[\log x_i^y+\alpha\log s_i^y+
\log x_i^m+\alpha\log s_i^m+
\log c_i^o+\gamma\log h_i^o]\,dQ.
\]
The young continuation utility has weights \(\beta,\beta^2\), and current mature continuation utility has weight \(\beta\); these terms are constant here. Dropping them does not change the unit weights on current young, mature and old utility.

The exact solution is
\[
X=\frac{\bar x^y+\bar x^m+\bar c^o}{3},\qquad
c_i^{y,F}=c_i^{m,F}=\chi n_i+X,\quad c_i^{o,F}=X,
\tag{1}
\]
\[
h_i^{y,F}=h_i^{m,F}=\min\{H_i,\kappa n_i+t\},\qquad
h_i^{o,F}=\min\{H_i,(\gamma/\alpha)t\},
\tag{2}
\]
where \(t=\alpha/\lambda_H\) clears all **three** housing demands. Strict concavity gives the unique allocation.

If planner caps are slack, put \(S=\bar s^y+\bar s^m+\bar h^o\). Then
\[
t=\frac{\alpha S}{2\alpha+\gamma},\qquad
\boxed{\bar h^{y,F}>\bar h^{y,eq}
\iff \alpha(\bar s^m+\bar h^o)>(\alpha+\gamma)\bar s^y.}
\tag{3}
\]
Individual young housing rises exactly when \(s_i^y<t\). Thus old housing exceeding young housing is insufficient on its own: the mature cohort's adult-space claim must enter the comparison. At \(\alpha=\gamma\), (3) says that mean mature-plus-old adult space exceeds twice mean young adult space.

Two useful capped sufficient versions are available.

* If \(\alpha\ge\gamma\), \(n_i>0\), and \(\mathcal H<3\int H_i\,dQ\), then (2) gives
  \[
  \bar h^{y,F}=\bar h^{m,F}>\mathcal H/3.
  \]
  Consequently \(\bar h^{m,eq}+\bar h^{o,eq}\ge2\bar h^{y,eq}\) suffices for a strict young gain. Old planner housing cannot be capped everywhere unless all three ages are capped everywhere, contradicting unused capacity. Where old housing is uncapped, each family home is strictly larger.
* Without requiring \(\alpha\ge\gamma\), retain the strict resource comparison in (3) and impose
  \[
  H_i-\kappa n_i\ge\bar s^y\quad\forall i,
  \qquad Q\{H_i-\kappa n_i>\bar s^y\}>0.
  \tag{4}
  \]
  At \(t=\bar s^y\), each family receives mean housing \(\kappa\bar n+\bar s^y\), while old housing is at most \((\gamma/\alpha)\bar s^y\). Condition (3) leaves excess housing, so clearing requires \(t>\bar s^y\). Condition (4) makes the resulting young aggregate increase strict.

These are allocation conditions. A separate household/equilibrium argument must derive the required reference comparisons. A strict allocation change implies a strict full-planner utility-sum gain; it does not imply every household gains utility.

There is a short link to the regime being tested by the household agent. Suppose young, mature and old housing and the old estate floor are slack, and mature saving is interior. Write \(b=\beta/q\), \(\ell=(1-\phi+q\tau)/D>0\), and let \(\mu_y\) be the young cash multiplier. The exact conditions give
\[
s_m=\frac{\alpha x_m}{p},\quad c_o=bx_m,\quad
h_o=\frac{\gamma b x_m}{p},\qquad
\frac\alpha{s_y}=\frac{pb}{x_m}+\ell p\mu_y.
\]
Consequently the typewise comparison underlying (3) is
\[
(\alpha+\gamma b)(b+\ell\mu_y x_m)>\alpha+\gamma.
\]
At \(\beta=q\), a strictly positive young finance multiplier makes this strict for any positive \(\ell\); no \(\phi<q\) restriction is needed for this allocation implication. If all other types have weak inequalities and a positive owner mass has strict finance, aggregation gives (3). Thus \(\phi=.8\) is not excluded by the planner argument. Establishing that this household regime exists at that LTV is a separate task; binding mature/rental caps require the capped comparisons above rather than these uncapped equalities.

## 2. Financial settlement includes the mature cohort

At date \(t\), write \(p_t=(1+q\tau_t)P_t-qP_{t+1}\). For current young owners set
\(\Delta a_m=-P_{t+1}\Delta h_y\); for current mature owners set
\(\Delta a_o=-P_{t+1}\Delta h_m\). These adjustments preserve the respective future entering resources, including the value of housing titles. Because fertility is fixed, each future decision problem and its real optimum are unchanged.

For an old owner, financial estate saving satisfies
\(e=a^e/q+P_{t+1}h_o\). Set \(\Delta a^e=-qP_{t+1}\Delta h_o\) to keep its net estate fixed. Renters retain their future financial positions. Each of the three current ages then requires the transfer
\[
g_i^a=\Delta c_i^a+p_t\Delta h_i^a.
\]
The transfers sum to zero by the two current resource constraints. As in the existing settlement, intermediaries offset the future value of changes in the **entire** rental stock with bond positions, so owner and intermediary financial adjustments cancel in aggregate.

Young mortgage limits, mature \(a_o\ge0\), and old financial-estate floors can all obstruct these offsets privately. Their relaxation is part of this planner benchmark. This is not implementation by the property tax. Allowing mature refinancing in the competitive model changes reference allocations and any mechanism argument, but not formulas (1)–(4) conditional on a reference allocation.

## 3. Private fertility remains a short conditional comparison

In the primary candidate, the correct interior lifetime fertility condition is
\[
\boxed{\frac{\vartheta}{n_i}
=\frac\chi{x_i^y}+\frac{\alpha\kappa}{s_i^y}
+\beta\left(\frac\chi{x_i^m}+\frac{\alpha\kappa}{s_i^m}\right).}
\tag{5}
\]
The mature terms follow from the continuation envelope, including when mature choices reoptimize. The two-age condition containing only current child costs is not valid here.

Nevertheless, give a parent a new current gross bundle while preserving the same mature continuation problem \(V_m(A_m,n)\), retained tenure and future prices. Its continuation derivative is the same in both comparisons at the original \(n_i\). Strict concavity therefore gives the exact finite test
\[
\boxed{n_i^{S}>n_i\iff
\chi\left(\frac1{x_i^y}-\frac1{x_i^{y,F}}\right)
+\alpha\kappa\left(\frac1{s_i^y}-\frac1{s_i^{y,F}}\right)>0.}
\tag{6}
\]
Here \(x_i^{y,F}=X\), and all residuals are evaluated at the original fertility, which remains feasible. If goods fall by \(d\ge0\) and housing rises by \(\Delta h\ge0\), the weak condition is equivalently
\[
d\le\frac{\alpha\kappa(x_i^y)^2\Delta h}
{\chi s_i^y(s_i^y+\Delta h)+\alpha\kappa x_i^y\Delta h}.
\]
Thus the finite sign can remain short even though the baseline fertility equation changes. The future terms affect the initial choice and response magnitude; their cancellation in (6) does not remove them from lifetime utility. Private reoptimization after the transfer can change future housing demand, so (6) alone does not establish a future market-clearing path.

## 4. Joint fertility needs an explicit future convention

The fixed-fertility planner is complete. Reusing its objective with \(\vartheta\log n\) added, while treating every continuation term as constant, is invalid: young fertility changes that household's future mature needs.

One coherent conditional version fixes promised future **gross** mature bundles \((C_{mi},H_{mi})\) and adds
\[
\beta\log(C_{mi}-\chi n_i)
+\beta\alpha\log(H_{mi}-\kappa n_i)
\tag{7}
\]
to the joint objective, retaining the corresponding positive domains. Future aggregate goods and housing then remain unchanged; extra child needs reduce the parent's adult resources within its promised bundle. This is a commitment convention, and generally prevents unconstrained future reoptimization. It has not been adopted.

Under this convention, the joint young fertility condition is
\[
\frac{\vartheta}{n_i}
=\chi\lambda_C+\kappa(\lambda_H+\eta_i^y)
+\beta\left[\frac\chi{C_{mi}-\chi n_i}
+\frac{\alpha\kappa}{H_{mi}-\kappa n_i}\right],
\tag{8}
\]
where \(\eta_i^y\) is the young housing-cap multiplier. Even with slack caps, heterogeneous future promises generally produce heterogeneous joint fertility.

Alternatively, preserving future **adult** consumption and space requires additional future goods \(\chi\Delta n_i\) and housing \(\kappa\Delta n_i\). Their funding and physical reallocation must be specified. Keeping only each future entering financial resource fixed and using \(V_m(A_m,n)\) is a well-defined conditional individual comparison, but does not automatically preserve future aggregate resource feasibility after all households reoptimize.

Therefore there is no unconditional joint-fertility implication from (3). Once young fertility changes, current mature fertility remains predetermined, so the equality \(h_y^F=h_m^F\) and the one-third argument cannot simply be reused for the joint optimum. Removing child needs from the mature stage eliminates (7), but is a substantively different, unadopted model.


---

# Three-stage stationary-equilibrium audit

Independent bounded analytical review, September 9, 2026. The candidate in `benchmark.md` is a test, not an adopted model. No numerical runs, builds, or source edits were used. This note establishes a low-dimensional clearing reduction and an analytically nonempty, deliberately strong benchmark. It does not establish general closed-form equilibrium, a transition theorem, or a joint-fertility welfare theorem.

## 1. Exact clearing dimension with the original tenure choice

For each candidate stationary rental price and rebate, `(p,T)`, set `P=p/D`, where `D=1-q+q tau`. Solve each type's two conditional lifetime problems, including all young and mature child needs and the actual caps and financial constraints. Let `V_i^d`, `n_i^d`, and `h_{ji}^d` denote the resulting tenure-conditional value, young fertility, and stage-j housing. Use the original independent logistic taste to obtain

\[
 \pi_i(p,T)=\operatorname{logistic}
 \left((V_i^O(p,T)-V_i^R(p,T)+\bar\xi)/\sigma_\xi\right),
 \qquad \sigma_\xi>0.
\]

Define total housing per stationary cohort and mean fertility by

\[
 \mathcal H(p,T)=\int\left[\pi_i\sum_jh_{ji}^O+
 (1-\pi_i)\sum_jh_{ji}^R\right]dF(i),
 \quad
 \mathcal N(p,T)=\int[\pi_i n_i^O+(1-\pi_i)n_i^R]dF(i).
\]

The complete stationary clearing system is

\[
 \nu\mathcal N(p,T)=1,\qquad
 T=\frac{q\tau p}{3D}\mathcal H(p,T),\qquad
 N=\frac{\bar H}{\mathcal H(p,T)}.
\]

Thus there are **two scalar clearing unknowns**. At zero tax there is **one scalar price equation**, since `T=0`. This is exact for arbitrary compact income-wealth heterogeneity and endogenous type-specific tenure probabilities. Caps and regime changes complicate conditional demand but do not add stationary clearing unknowns.

This reduction is not a general explicit solution for price. Even the uncapped borrowing-owner branch generically requires a cubic household root; arbitrary heterogeneity and logistic mixing then produce integrals of those roots and values. There is no justified claim that these aggregate equations become a polynomial or have closed-form roots for an arbitrary distribution. Population is endogenous through the stated stock equation; treating it as independently fixed would change this reduction.

## 2. The mature-saving test that high income cannot bypass

Write `K=1+gamma+omega_B`, `J=1+alpha+beta K`, `W_i=w_i+T`, `U_i=y_mi+T+q(y_oi+T)`,

\[
 L=(1-\phi+q\tau)P,\quad
 \delta=(1-\phi/q)P,\quad t=\chi+\kappa p.
\]

On the branch with binding young owner finance and slack mature/old owner caps and old estate floor,

\[
 R_i=U_i+\delta h_{yi}-tn_i,\quad x_{mi}=R_i/J,
 \quad h_{mi}=\kappa n_i+\alpha x_{mi}/p,
 \quad z_{oi}=\beta Kx_{mi}/q.
\]

The mature owner's actual net financial saving is therefore

\[
 a_{oi}^O=x_{mi}\left(\frac{\beta K}{q}-\frac{\alpha}{D}\right)
          -P\kappa n_i-(y_{oi}+T).
\]

With positive children and old income, **`beta K D>q alpha` is necessary** for positive saving on this uncapped owner branch. Raising mature income cannot rescue the branch if that coefficient is nonpositive. Given the positive coefficient, a sufficiently large mature-income ratio does deliver strict saving. The old owner floor is slack exactly when `omega_B D>q gamma` on this branch. Both restrictions must be stated separately from young borrowing.

## 3. Active renter caps are compatible with mature saving

Suppose the young renter has `a_m=0` and slack housing cap, while mature and old renters occupy `H_R`. Set `K_c=1+omega_B`, `J_c=1+beta K_c`. Direct substitution in the frozen budgets and the mature saving Euler equation gives

\[
 x_m=\frac{y_m+T+q(y_o+T)-\chi n-(1+q)pH_R}{J_c},
 \quad c_o=\frac{\beta x_m}{q},\quad
 e=\frac{\omega_B\beta x_m}{q^2},
\]
\[
 a_o=\frac{\beta K_c x_m}{q}+pH_R-(y_o+T).
\]

The cap inequalities are

\[
 \frac{\alpha}{H_R-\kappa n}\ge\frac{p}{x_m},\qquad
 \frac{\gamma}{H_R}\ge\frac{p}{c_o}.
\]

At fixed bounded young income, fertility, old income, price and rental cap, sufficiently high mature income makes both cap inequalities and `a_o>0` strict. There is no inconsistency between binding mature rental housing and positive saving.

The young renter's fertility condition is, however,

\[
 \frac{\vartheta}{n}
 =\frac{\chi+\kappa p}{x_y}
  +\frac{\beta\chi}{x_m}
  +\frac{\beta\alpha\kappa}{H_R-\kappa n},
 \quad (1+\alpha)x_y+(\chi+\kappa p)n=w+T.
\]

The last term survives as mature income becomes large. Dropping it would remove the candidate's mature child-space needs. Consequently the uncapped high-income fertility limit is not the exact limit for a capped renter. The following construction controls this term rather than omitting it.

## 4. A finite, primitive, nonempty benchmark with active rental caps

This subsection first sets `tau=0`; positive tax is handled separately below. It accommodates **every fixed `phi` in `(0,1)`, including `.8`**, without imposing `phi=q`. It retains the original logistic tenure rule, finite tenure caps, and a compact heterogeneous distribution. The sufficient conditions are intentionally stronger than necessary and are not a calibration claim.

Choose positive preferences and needs, `q` in `(0,1)`, and `alpha>=gamma`. Choose the estate weight so that

\[
 \omega_B(1-q)>q\gamma,\qquad
 \beta K(1-q)>q\alpha.
\]

Let young endowments have any compact distribution on `[w_min,w_max]`, with `w_min>0`, and write `mu_w=E[w]`. Let old-income ratios lie in a compact positive interval with `y_o/w<=o_bar`. Parameterize mature income by

\[
 u_i=\frac{y_{mi}+q y_{oi}}{w_i}\in[u_{\min},u_{\max}],
 \qquad u_{\min}>q\bar o.
\]

The joint distribution within this compact support may be arbitrary, including arbitrary dependence between these ratios and young endowments. Thus the construction does not require identical household types or constant tenure shares.

Choose `0<chi<nu vartheta mu_w/(2E)`, where `E=1+alpha+vartheta`. Define price coefficients

\[
 d_0=1-q,\quad \ell_R=1,\quad \ell_O=\frac{1-\phi}{d_0},\quad
 \ell_{\min}=\min(\ell_R,\ell_O),\quad
 \ell_{\max}=\max(\ell_R,\ell_O),
\]
\[
 \varepsilon_O=\frac{\ell_O-1}{q},\quad
 \varepsilon_{\max}=|\varepsilon_O|.
\]

The young owner cash cost is `L=ell_O p` and the mature equity term is `delta=epsilon_O p`. The following **primitive price bracket** has positive endpoints:

\[
 p_- =\frac{\nu\vartheta\mu_w/(2E)-\chi}{\kappa\ell_{\max}},
 \qquad
 p_+ =\frac{2\nu\vartheta\mu_w/E-\chi}{\kappa\ell_{\min}}.
\]

### 4.1 Owner bounds, uniform over the bracket

Set

\[
 C_e=\varepsilon_{\max}/\ell_{\min},\quad
 t_+=\chi+\kappa p_+,\quad C_R=C_e+t_+/\chi.
\]

Young cash feasibility gives `h_y<w/(ell_min p)` and `n<w/chi`, hence

\[
 (u_{\min}-C_R)w\le R_i\le(u_{\max}+C_e)w,
 \qquad r_i:=\beta x_{yi}/x_{mi}
 \le\frac{\beta J}{u_{\min}-C_R}.
\]

For an explicit fertility bound, put

\[
 K_\rho=\varepsilon_{\max}/\ell_{\min},\quad
 C_n=t_++\kappa\varepsilon_{\max}p_+,\quad
 B_+=\chi+\kappa\ell_{\max}p_+,
\]
\[
 C_*=(1+\alpha)C_n+2\alpha K_\rho(B_++C_n),
\]
\[
 \bar r_O=\min\left\{q/2,1,\frac1{2\max(K_\rho,1)},
                         \frac{\nu\vartheta\mu_w}{4C_*}\right\}.
\]

Choose `u_min>C_R+beta J/bar r_O`. To check the bound, write `rho=L-delta r` and use the household first-order conditions to obtain

\[
 \frac{n_i^O}{w_i}
 =\frac{\vartheta}{\mathcal D_i},\qquad
 \mathcal D_i=[\chi+\kappa\rho+tr](1+\alpha L/\rho)
                  +\vartheta(\chi+\kappa L).
\]

For `r<=bar r_O`, `rho>=L/2` and

\[
 |\mathcal D_i-E(\chi+\kappa L)|\le C_*r
 <\nu\vartheta\mu_w/4.
\]

It follows pointwise that

\[
 n_i^O(p_-)>\frac{4w_i}{3\nu\mu_w},\qquad
 n_i^O(p_+)<\frac{4w_i}{7\nu\mu_w}.
\]

Strict owner saving is ensured by the additional finite threshold

\[
 u_{\min}>C_R+\frac{J}{G_O}
       \left(\frac{p_+\kappa}{d_0\chi}+\bar o\right),
 \quad G_O=\frac{\beta K}{q}-\frac{\alpha}{d_0}>0.
\]

To make both later owner houses larger than the young house, impose

\[
 u_{\min}>C_R+J\max\left\{\frac1{\alpha\ell_{\min}},
                              \frac{q}{\gamma\beta\ell_{\min}}\right\}.
\]

These are sufficient inequalities; all are finite for any fixed permitted `phi`.

### 4.2 Renter cap and income bounds

Choose the finite rental cap before choosing the final large mature-income bound:

\[
 H_R>\max\left\{\frac{w_{\max}}{p_-},
 \frac{\kappa w_{\max}}\chi+
 \frac{8\beta\alpha\kappa w_{\max}}{\nu\vartheta\mu_w}\right\}.
\]

This leaves all young renter housing below the cap and bounds the persistent mature child-space cost. Set

\[
 \bar r_R=\min\left\{q/2,
              \frac{\nu\vartheta\mu_w}{8(1+\alpha)\chi}\right\},
 \quad
 X_{\min}=\frac{u_{\min}w_{\min}-w_{\max}-(1+q)p_+H_R}{J_c}.
\]

Increase `u_min` until

\[
 X_{\min}>\max\left\{
 \frac{\beta w_{\max}}{\bar r_R},\quad
 \frac{p_+H_R}{\alpha},\quad
 \frac{q p_+H_R}{\gamma\beta},\quad
 \frac{q\bar o w_{\max}}{\beta K_c}\right\}.
\]

Then the capped-renter formulas in Section 3 have `x_m>=X_min`; young saving is zero with a strictly positive finance multiplier, both later caps have strictly positive multipliers, and mature saving is strictly positive, uniformly in prices and types.

For clarity, the exact renter fertility identity is

\[
 \frac{n_i^R}{w_i}=\frac{\vartheta}{
 E t+(1+\alpha)\chi r_i+
 (1+\alpha)\beta\alpha\kappa x_{yi}/(H_R-\kappa n_i)}.
\]

Since `(1+alpha)x_y<=w`, the last two denominator terms together are less than `nu vartheta mu_w/4`. Consequently

\[
 n_i^R(p_-)>\frac{4w_i}{3\nu\mu_w},\qquad
 n_i^R(p_+)<\frac{w_i}{2\nu\mu_w}.
\]

Choose any finite `u_max>u_min`. Finally choose a finite owner cap strictly above `H_R` and above

\[
 w_{\max}\max\left\{
 \frac1{\ell_{\min}p_-},\quad
 \frac\kappa\chi+\frac{\alpha(u_{\max}+C_e)}{Jp_-},\quad
 \frac{\gamma\beta(u_{\max}+C_e)}{qJp_-}\right\}.
\]

This verifies slack owner caps on the entire bracket. All parameter choices are finite, ordered primitive restrictions; no equilibrium price or tenure share was selected to make the proof work.

### 4.3 Actual equilibrium and full-planner housing direction

Both tenure branches have pointwise replacement brackets proportional to `w_i`. Therefore arbitrary endogenous **type-specific** logistic mixing preserves the bounds:

\[
 \mathcal N(p_-,0)>4/(3\nu),\qquad
 \mathcal N(p_+,0)<4/(7\nu).
\]

The conditional objective is strictly concave in the continuous choices and all displayed regime inequalities are strict. Its policy and value functions are continuous on the compact price bracket. The logistic probabilities and aggregate fertility are continuous as well. The intermediate value theorem therefore supplies an actual stationary price in `(p_-,p_+)`, followed by the positive cohort mass `N=Hstock/mathcal H`. At every finite primitive choice the logistic rule gives positive renter and owner probabilities; the renter caps bind for a positive mass of actual households. This does **not** supply a uniform lower bound on the rental share as mature income grows.

The market allocation has `h_m>h_y` and `h_o>h_y` pointwise for both tenures, so aggregate young housing is strictly below one third of the stock. Under the candidate's specified fixed-fertility full-planner settlement and common current-resource constraints, the housing solution for each stationary type-tenure pair is

\[
 h_y^F=h_m^F=\min\{H_d,\kappa n+\alpha/\lambda_H\},\qquad
 h_o^F=\min\{H_d,\gamma/\lambda_H\}.
\]

Since `alpha>=gamma` and `n>0`, each family house is at least as large as its old counterpart. Some cap capacity is unused because all market young houses are strictly below their caps and the housing stock is unchanged. Hence aggregate young planner housing is strictly above one third of the stock. **The full fixed-fertility planner therefore gives more housing to the young in this benchmark.** This argument includes the mature cohort and retains its inherited children and housing cap. It uses the stated common-resource planner closure; it does not independently replace the required title, rental-intermediary, and transfer ledger.

## 5. Positive taxes: an analytical local extension, not arbitrary-tax existence

The construction above uses zero tax only as the starting point. Its strict household regime inequalities and pointwise fertility brackets hold uniformly on a compact support and price interval. Conditional optima, values, and logistic weights are continuous in `(p,T,tau)` there. Therefore some strictly positive `T_bar` and `tau_bar` preserve every margin throughout `[p_-,p_+] x [0,T_bar] x [0,tau_bar]`.

Let `g(p,T)=mathcal N(p,T)-1/nu`. On this rectangle `g(p_-,T)>0` and `g(p_+,T)<0`. Since `mathcal H<=3H_O` and `D>=d_0`, reduce `tau_bar` if necessary so that

\[
 \frac{q\tau p_+H_O}{d_0}<T_{\rm bar}.
\]

For any such positive tax, the continuous map

\[
 (p,T)\longmapsto
 \left(\operatorname{proj}_{[p_-,p_+]}[p+\eta g(p,T)],
       \frac{q\tau p}{3D}\mathcal H(p,T)\right),\qquad \eta>0,
\]

maps the rectangle into itself. Brouwer gives a fixed point. The strict fertility signs exclude either price boundary, so the fixed point satisfies both exact clearing equations, with positive tax and rebate. All regime and fixed-fertility housing inequalities survive. This proves an analytically nonempty positive-tax family around the explicitly constructed primitive benchmark. The allowable tax interval is established by continuity, not an explicit sharp tax bound. It is not a claim for every property-tax rate or a numerical-neighborhood proof.

There is also a general explicit rebate bound, independently checked from the frozen budgets. Both tenures satisfy the same discounted lifetime identity:

\[
 c_y+qc_m+q^2c_o+q^3e+p(h_y+qh_m+q^2h_o)
 =w+qy_m+q^2y_o+(1+q+q^2)T.
\]

Because weighted housing is at least `q^2` times total housing, aggregate fiscal clearing implies, for `0<tau<3q/(1+2q)`,

\[
 T<\frac{\tau\,\mathbb E[w+qy_m+q^2y_o]}
 {(1-q)[3q-\tau(1+2q)]}.
\]

This is an a priori bound on any equilibrium, not by itself an existence or regime-verification theorem. It provides an explicit compact rebate range when additional analysis requires one.

## 6. Scope and remaining limits

- The two-equation reduction is general; the nonempty regime is a sufficient family with high mature relative to young and old income, strong enough estate demand for unencumbered mature owner saving, and large owner housing capacity relative to rental capacity. Arbitrary income supports need not satisfy it.
- The finite inequalities give a genuine analytical existence argument with endogenous price and tenure. Their conservatism and potentially large mature-income ratios are economic limitations, not hidden numerical calibration.
- No uniqueness, monotone equilibrium comparative statics, or stability result follows from this scalar/two-dimensional existence argument.
- The current full-planner housing gain is a fixed-fertility result. With young fertility allowed to change, future mature child goods and housing costs must enter the feasible intervention and welfare derivative. The capped-renter term derived above shows explicitly why this additional step is substantive.
- Allowing mature refinancing changes the candidate's household constraints. Where the present optimum already has strict positive mature financial saving, relaxing its lower bound locally leaves that optimum unchanged; globally and in other regimes it can change policies. The current construction does not justify suppressing that institutional assumption.



---

# Lead review and final strengthening

## Independent verification of the equilibrium construction

A second reviewer checked the owner fertility-denominator constant, both fertility boundary signs, the renter's persistent mature child-space term, saving and housing-cap inequalities, noncircular ordering of finite parameters, IVT, positive-tax Brouwer rectangle, and the rebate bound. The proof constructs the constrained-young candidate first and then verifies every financing/cap inequality. Positivity of the continuation surplus throughout the young cash simplex justifies this construction. The conclusion is existence within the verified bracket, not a statement about all equilibria globally.

## A finite strengthening for each young household and private fertility

Let \(M_x=\sup_i x_i^y\) and \(M_s=\sup_i s_i^y\), where the supremum includes actual tenure types. Put
\[
\Phi(z)=\int\min\{H_i,\kappa n_i+z\}\,dQ,
\qquad \Psi(z)=\int\min\{H_i,(\gamma/\alpha)z\}\,dQ.
\]
The sufficient resource tests are
\[
\bar x^y+\bar x^m+\bar c^o>3M_x,
\qquad \bar H/N>2\Phi(M_s)+\Psi(M_s).
\]
They imply planner adult consumption \(X>M_x\) and its housing root \(z>M_s\). Since every competitive young home in the constructed regime is strictly below its own cap,
\[
\Delta h_i^y=\min\{H_i-h_i^y,z-s_i^y\}>0
\]
for every young type. It is unnecessary to require each cap's residual capacity to exceed the maximum adult space of other households. A more conservative housing condition is
\[
\bar s^y+\bar s^m+\bar h^o>(2+\gamma/\alpha)M_s.
\]

These conditions can be attained through finite primitive strengthening of the zero-tax equilibrium construction. In that construction's notation define
\[
B_y=\frac{w_{\max}}{\ell_{\min}p_-},\qquad
b_\beta=\beta/q,\qquad
L_O=\frac{(u_{\min}-C_R)w_{\min}}{J},\qquad L_R=X_{\min}.
\]
Young adult consumption is at most \(w_{\max}\), young adult space is below \(B_y\), mature adult consumption is at least \(\min\{L_O,L_R\}\), and old consumption equals \(b_\beta x^m\). The finite bounds
\[
(1+b_\beta)\min\{L_O,L_R\}>3w_{\max},
\]
\[
2H_R-\frac{\kappa w_{\max}}\chi>(2+\gamma/\alpha)B_y,
\qquad
\frac{(\alpha+\gamma b_\beta)L_O}{p_+}>(2+\gamma/\alpha)B_y
\]
therefore suffice, independently of the endogenous tenure shares. The first housing bound controls mature-plus-old renter adult space, the second controls the corresponding owner space. First enlarge \(H_R\), then the mature-income lower bound, then \(H_O\). This preserves the ordering of primitive choices and all earlier conditions. Increasing mature income also preserves the mature/old rental caps. All strict inequalities survive the established small-positive-tax extension.

The fixed-fertility planner consequently gives every young parent strictly more current gross consumption and housing. At the same continuation problem, the exact private fertility derivative is positive. This is a conditional individual choice result: the future aggregate goods and housing feasibility of all parents taking that response has not been established. In particular it is not an implemented transfer policy or a general-equilibrium fertility path.

## Main unresolved choices

The candidate keeps child needs in both family ages, retains tenure for life, and introduces a mature net financial lower bound. These remain unadopted model choices. The high-LTV analytical existence family uses young rental caps that are slack and mature/old rental caps that bind. Its required income ratios can be much larger than in the separate conditional 2.5--2.8 income-ratio example. The latter example is not a stationary-equilibrium proof. The proof of a positive-tax interval does not supply its numerical or sharp primitive endpoint. The existing two-age transition theorem cannot simply be attached to this three-age benchmark.
