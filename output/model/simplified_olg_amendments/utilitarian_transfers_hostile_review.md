# Hostile audit: utilitarian transfers with conventional finance

September 8, 2026. Sources: `utilitarian_transfers_review.md` (all sections) and
`latex/JMP_DS_suggestions/simplified_olg_conventional_finance.tex` (full file).
Analytical verification only; no simulations or numerical reference equilibrium.

**Verdict:** equations (3)–(8) implement the claimed finite improvement in the
economy conditional on previously committed fertility and tenure. The primitive
conditions are nonempty, including with positive property taxes and
\(\phi>q\). No financing, estate, or title-accounting blocker was found. An
unqualified claim of equilibrium in the original model, with fertility freely
reoptimized at intervention, is false. The qualifications and missing explicit
arguments below should accompany the theorem.

## 1. Exact household implementation — PASS

Write \(A_P=(1+q\tau^p)P\). Under predetermined transfers \(G,J\), the original
young budget and mortgage are
\[
c+qa'+A_Ph=w+G,\qquad qa'+\phi Ph\ge0,
\]
and old resources are \(z=a'+Ph+v+J\). The proposed changes satisfy
\[
\Delta a'=-\phi P\varepsilon/q,\quad
\Delta c+q\Delta a'+A_P\Delta h=\Delta x+L\varepsilon=G,
\]
\[
q\Delta a'+\phi P\Delta h=0,\qquad
\Delta a'+P\Delta h+J=0.
\]
Thus both cash accounting and the existing mortgage covenant hold exactly.
The mortgage principal increases by \(\phi P\varepsilon\), with zero young
gross bond assets, and is fully repaid next date. The larger inherited title
and \(J\) exactly offset its larger repayment.

With fixed \(n\), eliminate \(a'\). The lifetime and current-cash constraints
are affine in \((x,s,z)\); their multipliers obey
\(\lambda=\beta K/(qz)\), \(\mu=1/x-\lambda\). Equation (3) gives exactly
\[
\alpha/(s+\varepsilon)=\lambda p+\mu_\varepsilon L.
\]
For sufficiently small positive \(\varepsilon\), \(\mu_\varepsilon>0\).
Indeed \(\varepsilon<\alpha/(p\lambda)-s\) preserves this margin. The
denominator in (3) then remains positive. The young cap also remains slack.

Global optimality does not rely on a local second-order test. In the full
variables \((x,s,c^2,h^2,e)\), utility is a positive weighted sum of logarithms,
strictly concave on its positive domain. All budgets, caps, and the estate
floor are affine restrictions. The displayed allocation satisfies the full
KKT conditions and is therefore the unique conditional optimum. Old choices
are unchanged because inclusive old resources and prices are unchanged; their
gross financial saving \(a^e=q(e-Ph^2)>0\) is unchanged too.

The transfers must be fixed amounts attached to predetermined identities.
An equation used to calculate the schedule from the baseline is permissible;
recalculating grants from realized housing would change the household FOCs.

## 2. Residual funding and welfare — PASS

For \(\phi\ge q\), \(0<L\le p\), so condition (6) reduces to
\(\alpha(1+\omega_B)\le\gamma L/p\). On the binding branch,
\[
A_\varepsilon=\frac{\rho_\varepsilon^2}{\alpha L}
\ge\frac L\alpha\ge\frac{(1+\omega_B)p}{\gamma}=1/g-p.
\]
Consequently \(R_\varepsilon\ge0\), starting from \(R_0=0\). The matching
old tax \(D=\varepsilon/g\) induces exactly the required housing reduction,
with consumption and estate proportional to resources. Its estate floor
remains slack under proportional scaling.

The welfare derivative (7) is correct. In particular
\(\Lambda-m=(\beta/q-1)m+\mu>0\) and
\(\alpha/s-pm=(\beta/q-1)pm+\mu L>0\). Low marginal utility among capped
funders makes the remaining term nonnegative. Continuity supplies a strictly
positive finite policy step. This uses the stated welfare weights: weight one
on each current old household's remaining utility and on each current young
household's lifetime utility, the latter containing \(\beta\).

## 3. Heterogeneous masses — PASS, make integration explicit

Stationarity supplies identical type-and-owner measures among current young
and current old. Select the same submeasure in both cohorts; a continuous
ownership-taste draw permits subdivision even when endowments have atoms.
Each selected young owner is matched by type with an old owner having its
baseline \(z_i\). No equality between total recipient and funder masses is
needed.

For selected measure \(\mathcal A\) and funder mass \(M_F>0\), let
\(\mathcal R(\varepsilon)=\int_{\mathcal A}R_i(\varepsilon)\,d\mu\) and
tax each funder \(\mathcal R/M_F\). Write
\(\bar m_F=M_F^{-1}\int_Fm_f\,d\mu\). Then
\[
W'(0)=\int_{\mathcal A}\left[
(\Lambda_i-m_i)A_i+(\alpha/s_i-pm_i)
+(m_i-\bar m_F)(A_i+p-1/g)\right]d\mu>0.
\]
Equation (8) ensures \(\bar m_F<m_i\). Restrict groups to positive-measure
level sets with uniform positive mortgage, cap, resource, and funder margins.
Their smooth allocations and bounded derivatives justify differentiation and
a common finite step. A tiny funder mass only reduces the admissible policy
size; it does not invalidate existence.

## 4. Funder bound and its nonemptiness — PASS

The covenant implies
\[
z_f\ge y_f^o-(\phi/q-1)_+P_+H_O.
\]
At \(z_f>KpH_O/\gamma\), old housing is strictly capped and
\(m_f=(1+\omega_B)/(z_f-pH_O)\). The maintained
\(\omega_Bd_p>q\gamma\) makes the estate floor slack at the cap threshold
and everywhere above it. Meanwhile \(z_i<M_{\mathcal S}/q\), hence
\(m_i>qK/M_{\mathcal S}\). Both parts of (8) are therefore sufficient.

Its dependence on the full-distribution mean is not circularly impossible.
Here is a wholly analytical compatibility construction. Fix \(q\in(0,1)\)
and \(0<\tau^p<2\). Choose \(\phi>q\) sufficiently close to \(q\) that
\(\ell>1/(1+q)\). Choose \(\beta\ge q\), \(\gamma>0\),
\(\omega_B>q\gamma/d_p\), and
\(0<\alpha<\gamma\ell/(1+\omega_B)\). Increasing \(\vartheta\) makes
\[
k=\frac{1+\alpha+\vartheta+\beta K}
{1+\ell(\alpha+\vartheta)}-1\longrightarrow\ell^{-1}-1<q.
\]
Fix \(\nu,\kappa>0\), \(H_R>\kappa/\nu\), and \(w_0=W>0\). Choose
\(\vartheta\) sufficiently large that \(k<q\) and
\(\vartheta\nu>\alpha\kappa/(H_R-\kappa/\nu)\); then choose \(\chi>0\)
sufficiently small to satisfy the full existence inequality. Choose ordinary old income
\(V_S\) with \(qV_S>kW\). The eligible-income restriction now has no rebate
term because \(k<q\).

For any fixed small \(c>0\), add funders with old income \(V>V_S+c\) and
mass \(\eta=c/(V-V_S)\), replacing the same mass of ordinary types. Then
\[
\bar M_0=W+qV_S+qc
\]
is independent of \(V\). Consequently \(\bar T,P_+,M_{\mathcal S}\) are
independent of \(V\); \(P_-\) is also unchanged. Choose a finite \(H_O\)
satisfying the eligible-cap bound. The entire right side of (8) is now a
finite constant, so sufficiently large finite \(V\) satisfies (8), with
strictly positive funder mass. Every resulting distribution is bounded.
Logistic tenure supplies positive owner mass in both groups. No equilibrium
price, multiplier, or numerical reference point was assumed.

## 5. Property taxes, estates, titles, and finite tail — PASS

At date zero, \(D+R-G=qJ\). Government buys that amount of external bonds;
their date-one proceeds pay exactly \(J\), leaving zero program assets and
liabilities. With \(\phi\ge q\), this is saving, not borrowing. Aggregate
occupied housing and both cohort masses remain fixed, so the original
property-tax revenue and common rebate are unchanged. Rental demand and
intermediary pricing are unchanged.

Estate payments must also enter external resource accounting. For one matched
pair, write \(\Delta x=x_\varepsilon-x\). Across its matching old donor and
funders, the total consumption reduction is
\[
D/K+R/(1+\omega_B)=\Delta x/(1+\omega_B),
\]
because \(D=Kp\varepsilon/\gamma\) and
\(R=\Delta x-(1+\omega_B)p\varepsilon/\gamma\). Therefore aggregate current
consumption and next-date current-old estate payments change by
\[
\Delta C_0=\frac{\omega_B\Delta x}{1+\omega_B},\qquad
\Delta E_1=-\frac{\omega_B\Delta x}{q(1+\omega_B)},\qquad
\Delta C_0+q\Delta E_1=0.
\]
The same identities integrate across unequal masses. Extra consumption is
financed by smaller estates, fully counted in the old households' utility;
there is no missing external subsidy.

At date one, affected young households inherit \(\varepsilon\) more title;
current old households leave \(\varepsilon\) less title for death sales.
The former sell their extra title when choosing their unchanged old housing,
exactly offsetting the latter reduction in death sales. Their old portfolios
and estates are their original ones. Date-one entrants face the original
prices and rebates; estates do not determine entrant wealth. From date two
the inherited household state is identical to baseline. All new mortgages
and government trades settle at date one, and no terminal refinancing or
convergence assumption is required.

## 6. Knife-edge policy — PASS

Under (9), \(L=p\), \(J=R=0\), and
\(\alpha/[p(1+\alpha)]=\gamma/(Kp)=g\). A grant \(b\) and matched old tax
\(b\) therefore clear housing exactly. The displayed welfare expression,
positive derivative interval (10), mortgage-binding bound, and young-cap
bound are correct. Their upper bounds are strictly positive initially.
Interval (10) already implies \(b<v\), preserving donor resources; the
estate floor remains slack by proportional scaling.

## 7. Scope and wording — FIX; unrestricted-fertility claim BLOCKED

Call the result an equilibrium **conditional on committed fertility and
tenure**, or specify an unexpected intervention after those commitments and
before purchase. It is not an equilibrium with unrestricted date-zero
fertility choice. At the baseline optimum,
\(\vartheta/n=\chi/x+\alpha\kappa/s\). Both \(x\) and \(s\) strictly
increase under (3), so at the proposed unchanged \(n\),
\[
\vartheta/n-\chi/x_\varepsilon-\alpha\kappa/s_\varepsilon>0.
\]
A small fertility increase holding total consumption, housing, and financial
choices fixed is feasible and improves utility. Thus the original fertility
FOC fails. Restoring fertility requires a separate equilibrium and demographic
argument. The review correctly warns against reusing its finite-tail proof
for that extension.
