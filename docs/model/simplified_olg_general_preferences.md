# General preferences: housing allocation and fertility

Working mathematical memo, September 8, 2026. This is a separate exploration,
not a change to the household specification or a policy-transition theorem.
The maintained planner is the dated benchmark recorded in
`simplified_olg_utilitarian_work.md`: it reallocates **all current consumption
and housing**, first fixing individual fertility and tenure. It preserves
future real continuation opportunities and net estates, relaxes young and old
financing restrictions while honoring obligations, and retains physical caps.
The narrower planner discussion in the existing LaTeX note is superseded for
this exercise.

## 1. Primitives and the planner

Let young utility be \(U^y(x,s,n)\), where \(x=c-\chi n>0\) and
\(s=h-\kappa n>0\), with \(\chi,\kappa>0\). Let old utility be
\(U^o(c_o,h_o,e)\). Both are increasing and jointly concave in their arguments;
assume differentiability and interior goods choices for the formulas below.
Strict concavity gives uniqueness where imposed. The induced gross-bundle
utility \(u^y(c,h,n)=U^y(c-\chi n,h-\kappa n,n)\) remains concave, but
\(u^y_n=U^y_n-\chi U^y_x-\kappa U^y_s\) need not be positive.
All comparisons use the specified **cardinal** utility scales.

Each age has mass \(N\) and the same probability law \(Q\) of paired types
and retained tenures. Write \(H_i\in\{H_R,H_O\}\), \(H_R<H_O\), for a
household's physical cap. Suppressing constant continuation and taste terms,
the planner maximizes
\[
N\int[U^y(c_i^y-\chi n_i,h_i^y-\kappa n_i,n_i)
+U^o(c_i^o,h_i^o,e_i)]\,dQ
\]
subject to
\[
N\int(c_i^y+c_i^o)dQ=C,\qquad
N\int(h_i^y+h_i^o)dQ=\bar H,\qquad h_i^a\le H_i.
\]
Fertility \(n_i>0\), estates \(e_i>0\), \(C\), and \(\bar H\) are fixed at
the competitive reference. The owner estate floor is a financing restriction,
so it is relaxed here; the promised estate itself is unchanged.

For resource multipliers \(\lambda_C,\lambda_H\) and cap multipliers
\(\eta_i^{a*}\ge0\), optimality requires
\[
U^y_x=U^o_c=\lambda_C,\qquad
U^y_s=\lambda_H+\eta_i^{y*},\qquad
U^o_h=\lambda_H+\eta_i^{o*}.
\]
Together with complementary slackness, these conditions are sufficient
under concavity and feasibility. With
nonseparable utility, consumption changes housing marginal utilities, so
solving housing alone generally does not solve this planner.

## 2. The competitive wedge needs no logarithms

At stationarity, \(q=1/R_f\), \(p=(1-q+q\tau^p)P\),
\(L_R=p\), and \(L_O=(1-\phi+q\tau^p)P\). Here \(p\) is the cost
of housing services, and \(L_d\) is housing's coefficient in the current
financing constraint. Define current resources \(w=y^y+b+T\) and old income
plus rebate \(v=y^o+T\); \(z\) denotes resources in the old-age budget,
after repayment of any mortgage from youth. Conditional young choices satisfy
\[
c+ph+qz=w+qv,\qquad c+L_dh\le w.
\]
Let \(\Lambda_i,\mu_i\) be their budget and financing multipliers, and
\(\eta_i^y\) their housing-cap multiplier. Write \(V_d(z)\) for optimized
old utility and \(m_i=V_d'(z_i)>0\). The envelope theorem and young
first-order conditions give
\[
q\Lambda_i=\beta m_i,\quad
U^y_x=\Lambda_i+\mu_i,\quad
U^y_s=p\Lambda_i+L_d\mu_i+\eta_i^y.
\]
These identities hold even when old choices meet their cap or estate floor.

The old budget is \(c_o+ph_o+qe=z\). Put \(\rho_i\ge0\) on the owner's
constraint \(e-Ph_o\ge0\), with \(\rho_i=0\) for renters, and let
\(\eta_i^o\ge0\) be the old housing-cap multiplier. Then
\[
U^o_c=m_i,\qquad U^o_e=qm_i-\rho_i,\qquad
U^o_h=pm_i+P\rho_i+\eta_i^o.
\]
Stationarity matches this future old household with a current old counterpart.
Their **direct current housing-utility gap** is exactly
\[
\boxed{U^y_s-U^o_h
=\left(\frac\beta q-1\right)pm_i
+L_d\mu_i+\eta_i^y-P\rho_i-\eta_i^o.}
\]
The estate floor raises old direct housing marginal utility. Omitting its
negative contribution to the age gap silently excludes a regular regime.
No restriction on \(\beta R_f\) was used. Its threshold in earlier sufficient
proofs therefore does not originate in logarithmic utility.

A positive gap supports a small housing transfer to an **uncapped** young
recipient from its old counterpart, with current consumption unchanged.
This proves an improving direction, not the full optimum's aggregate
direction. A positive young cap multiplier is not transferable room.

## 3. What concavity cannot establish

Increasing concavity alone orders neither age's housing marginal utility nor
its aggregate housing gain. Both \(\beta/q-1\) and the old financial wedge
can oppose the young financing wedge. Even with fixed cardinal scales,
concavity imposes no cross-age ordering of marginal utilities.

The normalization issue is substantive: replace \(U^o\) by \(A U^o\) and
\(\beta\) by \(\beta/A\), for any \(A>0\). All competitive choices and
young lifetime utilities remain unchanged, while the equally weighted
current-household planner assigns \(A\) times the previous weight to old
utility. This gives an admissible family under unrestricted preferences and
patience with identical market allocations but different planner objectives.
It is not an innocuous normalization of a fixed social criterion.
For uncapped log housing this gives
\[
H_Y^*(A)=C_n+\frac{\alpha}{\alpha+A\gamma}(\bar H-C_n),
\qquad C_n=\kappa N\bar n.
\]
With sufficiently large finite caps, varying \(A\) moves this allocation
across any interior reference young housing total.

## 4. A useful nonlogarithmic housing theorem

Consider the more structured class
\[
U^y=f(x)+\alpha g(s)+v(n),\qquad
U^o=f_o(c_o)+\gamma g(h_o)+b(e),
\]
with \(\alpha,\gamma>0\), increasing concave components, and strictly
concave \(g\). Housing
separates from consumption at fixed fertility. For an interior adult-space
solution, define \(r_y=(g')^{-1}(\lambda_H/\alpha)\) and
\(r_o=(g')^{-1}(\lambda_H/\gamma)\). Then
\[
h_i^{y*}=\min\{H_i,\kappa n_i+r_y\},\qquad
h_i^{o*}=\min\{H_i,r_o\}.
\]
Assume interior positive adult space exists. Extend the inverse to
\(+\infty\) below the derivative range, so a household demanding beyond its
cap is correctly clipped. The standard conditions \(g'(0+)=\infty\) and
\(g'(\infty)=0\) avoid inverse-range qualifications.

If \(\alpha=\gamma\) and planner caps are slack, adult space is identical
across both ages **for every such \(g\)**:
\[
r=\frac{\bar H/N-\kappa\bar n}{2},\qquad
H_Y^*=\frac{\bar H+\kappa N\bar n}{2},\qquad
\bar n=\int n_i\,dQ.
\]
Thus \(H_Y^*>H_Y^{eq}\) exactly when old mean reference housing exceeds
young mean reference adult space. An individual young household receives
more housing exactly when \(s_i^{eq}<r\).

More generally, \(\alpha\ge\gamma\) implies \(r_y\ge r_o\), hence
\(h_i^{y*}\ge h_i^{o*}\), strictly wherever the old counterpart is uncapped.
If
\[
\kappa N\bar n<\bar H<2N\int H_i\,dQ,
\]
not everyone can be capped, giving \(H_Y^*>\bar H/2\). Consequently the
reference ordering \(H_Y^{eq}\le H_O^{eq}\) is sufficient for a strict
aggregate young gain, allowing binding planner caps. It remains a condition
on the reference, not an implication of borrowing restrictions alone.

For the same class, the paired ordering \(h_i^{o,eq}\ge h_i^{y,eq}\)
and \(n_i>0\) directly imply
\[
\alpha g'(s_i^{eq})>\gamma g'(h_i^{o,eq})
\quad(\alpha\ge\gamma).
\]
This local child-space result still requires an uncapped young recipient.

## 5. Mapping and remaining gaps

The existing specification sets \(f=f_o=g=\log\),
\(v(n)=\vartheta\log n\), and \(b(e)=\omega_B\log e\).
The wedge and common-housing-utility results survive; logarithms additionally
give constant expenditure shares and explicit old-estate regime boundaries.
General primitives delivering the competitive age ordering, stationary
existence, and policy implementation are not established here. Neither an
aggregate housing gain nor a utilitarian gain implies every young household
receives more housing.

## 6. Fertility without logarithms

This section is the lead's derivation. Hold a parent's gross goods and housing
bundle, tenure and continuation opportunities fixed when defining its fertility
choice, and assume \(U^y\) is twice continuously differentiable. Write
\[
\widetilde U(c,h,n)=U^y(c-\chi n,h-\kappa n,n).
\]
Increasing utility in children is imposed on the primitive net-resource
utility, before paying their goods and space costs. It need not hold for
\(\widetilde U\) when gross goods and housing are held fixed.

At an interior choice the fertility condition is
\[
F:=U_n^y-\chi U_x^y-\kappa U_s^y=0.
\]
At a regular optimum, \(\widetilde U_{nn}<0\). Implicit differentiation gives
\[
dn=-\frac{A\,dc+B\,dh}{\widetilde U_{nn}},\qquad
A=U_{nx}^y-\chi U_{xx}^y-\kappa U_{sx}^y,\quad
B=U_{ns}^y-\chi U_{xs}^y-\kappa U_{ss}^y.
\]
Thus \(A\ge0\) and \(B\ge0\) are exactly the local conditions for fertility
to respond weakly positively to each resource separately. Concavity supplies
the denominator's weak sign, but a strictly negative second derivative is
needed for this derivative formula. Increasingness and joint concavity alone
do not determine the signs of \(A\) and \(B\).

For example, take any increasing, strictly concave functions \(f,g,v\) with
strictly negative second derivatives and set
\[
U^y(x,s,n)=f(x)+g(s+a n)+v(n),\qquad a>\kappa.
\]
This utility is increasing and jointly strictly concave. Yet
\(\widetilde U_{nh}=(a-\kappa)g''<0\), so fertility falls with additional
housing at every regular interior optimum. Children and space enter the
same utility component as substitutes. An interior example can be obtained
with square-root functions; no logarithmic functional form is needed for
this counterexample.

### A useful generalization of the existing preferences

The additive class
\[
U^y(x,s,n)=f(x)+\alpha g(s)+v(n)
\]
with increasing functions and \(f'',g'',v''<0\) gives
\[
dn=
\frac{-\chi f''(x)\,dc-\alpha\kappa g''(s)\,dh}
{-v''(n)-\chi^2f''(x)-\alpha\kappa^2g''(s)}.
\]
Both individual resource effects are strictly positive. When consumption
falls and housing rises, the numerator supplies the precise local tradeoff.
This result requires additive separability and curvature, not logarithms.
The current specification is the special case
\(f(x)=\log x\), \(g(s)=\log s\), and \(v(n)=\vartheta\log n\).
It reproduces the existing conditional fertility differential exactly.

### Joint choice by the planner

When the dated planner chooses fertility too, its young-household conditions
in net resources are
\[
U_x^y=\lambda_C,\qquad
U_s^y=\lambda_H+\eta_i^{y*},\qquad
U_n^y=\chi\lambda_C+\kappa(\lambda_H+\eta_i^{y*}),
\]
where \(\lambda_C,\lambda_H\) are the resource multipliers and \(\eta_i^{y*}\)
is the multiplier on \(s_i+\kappa n_i\le H_{d_i}\). This values children
through currently living parents only. The conditional derivative above does
not by itself establish greater average fertility at this joint optimum.
The logarithmic proof's explicit allocations and aggregation inequalities
require a new argument under general preferences. A policy transition also
requires the actual funded equilibrium changes in bundles and distributions.

Lead verification checked the implicit derivatives, the substitute example
and the specialization back to the current logarithmic formula algebraically.
No simulation, model revision or paper/slides edit was performed.
