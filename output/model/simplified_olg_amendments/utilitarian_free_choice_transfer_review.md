# A balanced transfer improvement with fertility and tenure freely chosen

Analytical review, September 8, 2026. This is a new construction for the
conventional-finance model, not an extension of the fixed-fertility proof by
assertion. No simulation or numerical reference equilibrium is used.

**Result.** At \(\phi=q\), four selected groups support a finite utilitarian
improvement with original household optimization, mortgage restrictions,
housing caps and estate restrictions. Young recipients buy more housing and
have more children; a second young group reduces fertility by exactly the
same amount. Housing moves from old to young, while aggregate fertility and
every future cohort mass are unchanged. The government balances at the
intervention date. No transition convergence assumption is needed.

## 1. Scope, information and stationary baseline

Use equal weights on current young lifetime utility and current old remaining
utility. These weights count young private old-age utility with \(\beta\).
Future entrants have the same mass and exogenous type distribution under both
paths, so their utility differences vanish. Estates do not determine entrant
wealth, and parental identity does not determine entrants' type distribution.

The government can identify selected households by predetermined identity,
endowments and the realized ownership taste \(\xi\). Each selected young
household has a uniform strictly positive initial ownership-value margin.
Its announced cash transfer applies regardless of its actual tenure choice.
Thus it still optimizes over both tenures; a sufficiently small intervention
preserves optimal ownership. This is stronger information than anonymous
endowment-only transfers. Eligibility must not be interpreted as a subsidy
conditional on subsequently owning. The original logistic taste distribution
is retained, although a taste-dependent transfer generally prevents using the
usual endowment-only logit formula without conditioning on eligibility.

Start at a stationary equilibrium with
\[
\phi=q,\quad p=(1-q+q\tau^p)P=L,\quad
E=1+\alpha+\vartheta,\quad K=1+\gamma+\omega_B.
\]
Let \(w\) denote young income, liquid wealth and the common rebate, and
\(v\) old income plus that rebate. A strictly finance-bound young owner
has \(z=v\), independently of its chosen current housing. Define
\[
a_n=\frac{\vartheta}{E(\chi+p\kappa)},\qquad
a_h=\frac{\alpha/p+\vartheta\kappa/(\chi+p\kappa)}{E},\qquad
g=\frac{\gamma}{Kp}.
\tag{1}
\]
An uncapped young owner chooses \(n=a_nw\), \(h=a_hw\), and adult
consumption \(x=w/E\). An uncapped old owner with a slack estate floor
chooses housing \(gv\).

Select equal positive masses of four groups:

- **Young A:** finance strictly binds; young housing is uncapped. Its old
  resources are \(v_A\), with \(gv_A<H\), where \(H\) is the owner cap.
- **Young B:** finance strictly binds and young housing is strictly at \(H\).
  Denote its privately chosen fertility and adult consumption by
  \(n_B,x_B\), and its cash by \(w_B\).
- **Old A:** stationary counterparts of Young A, with resources \(v_A\),
  uncapped housing and a slack estate floor.
- **Old B:** owners strictly at \(H\), with a slack estate floor. They may
  be stationary counterparts of Young B, as in the construction below.

Old tenure is retained as stipulated by the original model. Equal selected
masses can be drawn from different-sized groups because the logistic taste
distribution is atomless and assigns positive mass to strict owners.

## 2. Compact sufficient conditions and the actual transfer policy

For a capped, finance-bound young household, fertility is the unique solution
of
\[
\frac{\vartheta}{n_B}
=\frac{\chi}{w_B-pH-\chi n_B}
 +\frac{\alpha\kappa}{H-\kappa n_B},\qquad
x_B=w_B-pH-\chi n_B>0.
\tag{2}
\]
Define the positive derivative and the fertility-offset tax ratio
\[
b_B=\frac{\chi/x_B^2}
 {\vartheta/n_B^2+\chi^2/x_B^2
                  +\alpha\kappa^2/(H-\kappa n_B)^2},\qquad
q_B=\frac{a_n}{b_B},\qquad d=\frac{a_h}{g}.
\tag{3}
\]
These expressions contain no borrowing or housing multipliers.

**Proposition.** Suppose the strict group conditions above hold and
\[
\boxed{\quad
0<q_B<1,\qquad q_B+d>1,\qquad
\frac{E}{w_A}>\frac{q_B}{x_B}+\frac{Kd}{v_A}.
\quad}                                                        \tag{4}
\]
Then there is a finite strictly positive, date-zero, balanced cash-transfer
policy that raises utilitarian welfare with fertility and tenure privately
chosen. It raises aggregate young housing and lowers aggregate old housing.

Give each Young A household \(G>0\). Choose a tax \(Q(G)>0\) on each
Young B household to satisfy exactly
\[
n_B(w_B-Q(G))=n_B(w_B)-a_nG.                                  \tag{5}
\]
Since \(b_B>0\), this determines a unique small tax with
\(Q(0)=0\), \(Q'(0)=q_B\). Tax each Old A household
\(D(G)=dG\), and give each Old B household
\[
S(G)=Q(G)+D(G)-G.                                            \tag{6}
\]
By (4), \(S(G)>0\) for small positive \(G\), and
\(G+S=Q+D\) balances the government budget exactly, without bonds or later
taxes. All four amounts are fixed in advance; none depends on an individual's
actual choice.

**Household optimality and market clearing.** Small transfers preserve Young
A's uncapped, finance-bound regime. Its changes are
\(\Delta n_A=a_nG\), \(\Delta h_A=a_hG\). Young B still chooses its cap,
and (5) uses its actual fertility optimum after the tax. Old A optimally
reduces housing by \(gD=a_hG\), while Old B optimally stays at its cap
after a positive grant. Thus housing clears at the original price, and the
two fertility changes cancel exactly. The total housing stock and cohort
masses are unchanged, so the original property-tax rebate also balances.

For either affected young owner,
\[
a'=-Ph,\qquad qa'+\phi Ph=0,\qquad a'+Ph+v=v.
\tag{7}
\]
Consequently its inherited mortgage and title exactly offset when it becomes
old. Its old resources and all old consumption, housing, financial saving and
estate choices remain at their original optima. The date-one state may contain
different titles and net assets for A, but these have the same net resources;
both components are explicitly accounted for in (7).

Next-date young and old cohort masses are unchanged by (5). Entrants still
draw from the original exogenous distribution, face the original prices and
rebates, and receive no policy transfer. Their choices regenerate the original
date-two state. Initial old estates change, but their utility and asset
payments are included in their budgets and do not fund entering wealth.
Completed fertility and parentage change for selected households, not the
future payoff-relevant demographic or financial state.

The consolidated goods account also closes. Per unit of selected group mass,
young total nondurable spending changes by
\(\Delta C_Y=(1-pa_h)G-Q\). The two initial old budgets give
\(\Delta C_O+q\Delta E_O=Q-G+pa_hG\). Therefore
\[
\Delta C_Y+\Delta C_O+q\Delta E_O=0.                          \tag{7a}
\]
All other future goods and estate choices are unchanged. Equation (7a)
includes the initial old estates paid at date one; neither those payments nor
existing mortgage claims are omitted from the external resource account.

**Welfare.** Let \(m_{OB}>0\) be Old B's cash marginal utility. The exact
derivative, per unit of selected group mass, is
\[
\mathcal W'(0)
=\frac{E}{w_A}-\frac{q_B}{x_B}-\frac{Kd}{v_A}
 +m_{OB}(q_B+d-1)>0.                                        \tag{8}
\]
The last term is nonnegative and the preceding strict inequality is (4).
Young A's exact lifetime gain is \(E\log(1+G/w_A)\); its old utility is
unchanged. Young B's loss is its optimized current-utility change, including
its chosen fertility reduction. Old A's exact loss is
\(K\log(1-dG/v_A)\). When Old B has resources \(v_B\), its exact gain is
\((1+\omega_B)\log[1+S/(v_B-pH)]\).
These continuously differentiable expressions, together with strict branch,
solvency and ownership margins, give a finite positive welfare gain. This
constructs an equilibrium path; it does not infer its existence from an
unaccounted price response or a presumed convergent tail. QED.

## 3. An explicit admissible range for Young B

The first two inequalities in (4) are compatible with a strictly binding cap.
Put
\[
r=\frac{\chi}{p\kappa},\quad a=\alpha+\vartheta,
\quad t=\frac{\kappa n_B}{H},\qquad
f(t)=\frac{Et}{a}-\frac{\alpha}{a^2}
 +\frac{\alpha\vartheta}{a^2(\vartheta-at)}.
\]
Solving (2) for wealth gives the explicit inverse
\[
w_B=pH+\frac{\chi H}{\kappa}f(t),\qquad
f'(t)=\frac Ea+
 \frac{\alpha\vartheta}{a(\vartheta-at)^2}>0,\quad
f''(t)=\frac{2\alpha\vartheta}{(\vartheta-at)^3}>0.
\tag{9}
\]
The cap threshold is
\(w_*=H/a_h\), corresponding to
\(t_*=\vartheta/[\alpha(r+1)+\vartheta]\). Its tax ratio simplifies to
\[
q_*=1-pa_h\left(1-\frac1{\alpha r}\right),\qquad
q_*+d-1=a_h\left(\frac1g-p+\frac p{\alpha r}\right)>0.
\tag{10}
\]
If \(\alpha r>1\), then \(0<q_*<1\). More strongly, define
\[
t_\dagger=\frac{\vartheta}{a}
 \left[1-\sqrt{\frac{\alpha r}{E(\alpha r+a)}}\right].
\tag{11}
\]
One has \(t_\dagger>t_*\) exactly when \(\alpha r>1\).
Choosing any \(t_B\in(t_*,t_\dagger)\) in (9) makes the cap strict and
ensures (4)'s first two inequalities. This is an explicit parameter interval,
not an unspecified neighborhood of a numerical example.

It also supplies an exact finite implementation of (5):
\[
Q(G)=\frac{\chi H}{\kappa}
 \left[f(t_B)-f\left(t_B-\frac{\kappa a_nG}{H}\right)\right].
\tag{12}
\]
For \(0<G<H(t_B-t_*)/(\kappa a_n)\), Young B remains strictly capped.
Moreover \(Q'(G)\) decreases from \(q_B\) toward \(q_*\), so
\(0<Q'<1\) and \(S'>0\) throughout this finite interval. Also retain
\(G<H/a_h-w_A\), \(G<qEv_A/(\beta K)-w_A\), and \(dG<v_A\), plus
the selected households' ownership-value margins.

## 4. Analytical nonemptiness: a constructed family of primitive economies

The following parameterization verifies that the conditions can hold in an
actual stationary equilibrium. It chooses primitives; it is not a universal
theorem for already-fixed empirical endowments, prices or child requirements.

Choose finite \(\beta>0\), \(q\in(0,1)\), \(p,H>0\), positive preference
and child-cost parameters with \(\chi/(p\kappa)>1/\alpha\), and
\(\omega_B(1-q)>q\gamma\). Choose \(0<\eta<1\) and set
\(v_A=\eta H/g\), so Old A is uncapped. Choose \(t_B\) strictly inside
(11), giving \(w_B,n_B,x_B,q_B\) explicitly. Choose
\[
v_B>\max\left\{H/g,\;pH+\frac{\beta(1+\omega_B)x_B}{q}\right\}.
\tag{13}
\]
Then Old B is strictly capped and Young B's finance restriction strictly
binds. Finally choose positive \(w_A\) below all three explicit bounds
\[
\boxed{\quad
w_A<\min\left\{
\frac H{a_h},\;
\frac{qEv_A}{\beta K},\;
\frac E{q_B/x_B+Kd/v_A}
\right\}.
\quad}                                                       \tag{14}
\]
Every bound is strictly positive, so the interval is nonempty. These choices
establish all owner-branch and welfare conditions without a numerical anchor.

Give the two endowment types positive probabilities and retain any logistic
taste distribution with positive scale. Compute both tenure menus at the
auxiliary total resources \((w_i,v_i)\). Their finite values imply positive
masses of strict owners in both types. Selecting equal masses from sufficiently
high taste draws gives the required uniform ownership margins. Other
households keep their original optimizing choices.

For a positive property tax, write \(d_p=1-q+q\tau^p\), set \(P=p/d_p\),
and let \(\bar h_{\rm life}\le2H\) be average young-plus-old housing under
these auxiliary menus. Set
\[
N=\frac{\bar H}{\bar h_{\rm life}},\qquad
T=\frac{q\tau^pP\bar h_{\rm life}}2.
\tag{15}
\]
Choose \(\tau^p>0\) small enough that
\(T<\min_i\{w_i,v_i\}\); the explicit bound
\(\tau^p<(1-q)\min_i\{w_i,v_i\}/(2qpH)\) suffices.
Define primitive young income plus liquid wealth as \(w_i-T>0\), split
between two positive components, and primitive old income as \(v_i-T>0\).
All real menus remain unchanged: \(p,w_i,v_i\) are unchanged and the
estate floor stays slack under the stated preference inequality. Equation
(15) closes housing and the common rebate exactly.

If \(\nu\) is given, start with provisional \((\chi_0,\kappa_0)\), compute
their positive mean fertility \(\bar n_0\), and multiply both child
requirements by \(s=\nu\bar n_0\). Each household's fertility becomes
\(n_0/s\), leaving goods and space used by children, all other real choices,
and tenure-value differences unchanged. Both tenure values receive the same
\(-\vartheta\log s\) shift. Thus mean fertility becomes \(1/\nu\).
The ratio \(r\), the admissible cap interval, transfers and welfare
inequalities are preserved. The initial old distribution is the one generated
by these young choices, completing stationarity.

This proves a full-choice transfer improvement in a nonempty analytically
constructed family, with heterogeneous incomes, positive ownership and renting,
and positive property taxes. It makes no claim that the same restrictions hold
for a previously calibrated economy. Aggregate fertility is deliberately held
constant through voluntary offsetting responses; the proposition is not a
policy-induced population-growth result.

Verification: ten exact symbolic checks passed for the capped optimum,
inverse derivatives, cap threshold, tax-ratio interval, rebate, fiscal budget
and restored old resources. The consolidated goods identity was checked
directly from the dated household budgets. No model run was required.
