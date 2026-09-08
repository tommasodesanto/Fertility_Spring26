# Independent fertility and transfer audit

September 8, 2026. Read-only review of `simplified_olg_conventional_finance.tex`,
`utilitarian_fertility_path_review.md`, and `utilitarian_transfers_review.md`.
No equilibrium simulation, source edit, or convergence argument. This note is
the sole new artifact.

**Verdict: conditional fertility results pass; the endogenous-fertility
implementation needs an explicit qualification and corrected transfers.**
The existing fixed-fertility funding and welfare restrictions do not
automatically survive restored fertility, even at the intervention date.

## 1. Current gifts: pass, with precise kink wording

At fixed tenure and prices, a current gift strictly raises fertility for all
positive child costs and every young/old cap and estate-floor regime. On a
strict-finance, uncapped-young branch, the review's Hessian calculation gives
\[
n_w=\frac{LB_c+\chi A}{x^2(AD_n-B_c^2)}>0,
\quad LB_c+\chi A=
\frac{\alpha(\chi+L\kappa)}{s^2}+\chi\Gamma(L-p)^2>0.
\]
Here \(A=L^2/x^2+\alpha/s^2+\Gamma(L-p)^2\),
\(B_c=-L\chi/x^2+\alpha\kappa/s^2\), and
\(D_n=\chi^2/x^2+\alpha\kappa^2/s^2+\vartheta/n^2\).
The housing numerator is also strictly positive. With both young finance and
housing binding, \(n_w=\chi/(x^2D_n)>0\).

The omitted slack-finance/young-cap calculation is valid. Set
\[
k_c=\frac{\chi}{x^2D_n},\qquad
\Gamma=-\beta\mathcal V''(z)/q^2>0.
\]
Here \(k_c\) is fertility's response to total consumption at fixed housing.
The first-order conditions imply
\[
c_w=\frac{\Gamma}{\Gamma+(1-\chi k_c)/x^2}>0,
\qquad n_w=k_c c_w>0,
\]
since \(\chi k_c<1\). When young housing is also slack, write
\(E=1+\alpha+\vartheta\). Current expenditures are proportional to \(x\),
and \(x_w=\Gamma/(x^{-2}+E\Gamma)>0\).

Old optimized utility is strictly concave, not merely concave. Its smooth
regimes have the form \(A_o\log(z-z_0)+C\), with
\(A_o\in\{K,1+\omega_B,1\}\). The last case retains both the owner cap
and estate floor. The marginal value \(\mathcal V'=1/c^2\) is continuous
at regime boundaries; its derivative can jump. Thus there is no old
marginal-utility discontinuity invalidating strict gift monotonicity.
**Fix the wording:** use positive one-sided responses at kinks and strict
finite increases, rather than an everywhere-defined derivative. Replace
“needs only concavity” by “uses strict concavity of the retained log old-age
problem.” General concavity alone cannot establish strictness.

The reported rational example is exact: gift response
\(23123000/210723483>0\), fair-loan response
\(-36516490/210723483<0\). Repayment can reverse the gift result.

## 2. Finite bundle test: pass

For given total consumption and housing define
\[
F(n;c,h)=\vartheta/n-\chi/(c-\chi n)
                   -\alpha\kappa/(h-\kappa n).
\]
Its derivative in \(n\) is strictly negative. Every private optimum has
\(F(n;c,h)=0\): the finance constraints and physical cap contain total
\(c,h\), so holding them fixed permits this fertility first-order condition.
For unchanged \(\vartheta\) and baseline \(n_0\), feasibility of \(n_0\)
at the new bundle therefore gives exactly
\[
n_1\ge n_0\iff
\frac{\chi}{c_1-\chi n_0}+\frac{\alpha\kappa}{h_1-\kappa n_0}
\le
\frac{\chi}{c_0-\chi n_0}+\frac{\alpha\kappa}{h_0-\kappa n_0}.
\]
If \(c_1\le\chi n_0\) or \(h_1\le\kappa n_0\), positivity of the
new adult bundle directly implies \(n_1<n_0\), including equality cases.
This tests existing counterfactual bundles; it does not manufacture them or
hold consumption fixed during a mortgage reform.

## 3. Tenure: pass

For a common gift, \(W_w^m=1/x_m\), so the stated selection term is
\(\pi(1-\pi)(n_O-n_R)(1/x_O-1/x_R)/\sigma_\xi\). It has no universal sign.
At strict finance with \(L=p\), let \(x=w/E\) and
\(t_0=Bx/(qv)\in(0,1)\). An owner cash-cost reduction \(\delta=p-L>0\)
has
\[
x_\delta=xt_0H(p)/E>0,\qquad
n_\delta=\vartheta\left[
\frac{x_\delta}{\chi+p\kappa}
+\frac{x\kappa(1-t_0)}{(\chi+p\kappa)^2}\right]>0.
\]
Hence the selection term is negative for sufficiently small positive
\(\delta\). Centering ownership at one half and choosing sufficiently small
positive taste dispersion makes it dominate both positive conditional
responses. Bounded proportional endowment heterogeneity preserves the example:
allocations scale, tenure-value differences and the selection product do not.

An owner-contingent grant raises ownership and owner fertility. Baseline
\(n_O\ge n_R\) is sufficient for a finite increase, since
\[
\Delta\bar n=\pi_1\Delta n_O+
(\pi_1-\pi_0)(n_{O0}-n_{R0}).
\]
The stated primitive sufficient ordering also passes, with slack old caps:
\(L\le p\) and
\(p/(\alpha+B)\le\chi/\kappa\le p/\alpha\). At \(L=p\), the upper
bound makes relaxing the young cap fertility-increasing; the lower bound
makes subsequently lowering \(L\) fertility-increasing. Retained old estate
floors change constants, not the common uncapped coefficient \(K\).

## 4. Actual transfer construction: amend, do not import its theorem

At the construction's unchanged stationary prices, its fixed-old-resource
direction extends to the privately chosen fertility of an initially selected,
tenure-fixed household. Let
\(\lambda=\beta\mathcal V'(z)/q\), and for small \(\varepsilon>0\) set
\[
s_\varepsilon=s+\varepsilon,\qquad
x_\varepsilon=\frac{L}{\alpha/s_\varepsilon-(p-L)\lambda},\qquad
n_\varepsilon=\frac{\vartheta}
 {\chi/x_\varepsilon+\alpha\kappa/s_\varepsilon}.
\]
Both \(x\) and \(n\) strictly increase. Put
\(\Delta c=\Delta x+\chi\Delta n\) and
\(\Delta h=\varepsilon+\kappa\Delta n\). The correct grants are
\[
G=\Delta c+L\Delta h,
\qquad qJ=(p-L)\Delta h.
\]
Announce these as predetermined amounts, not payments contingent on actual
births. They preserve the binding covenant through
\(\Delta a'=-\phi P\Delta h/q\), and preserve old resources through
\(\Delta a'+P\Delta h+J=0\). The fertility and housing first-order
conditions then establish actual conditional optima. With a binding young
cap, fixed \(z\) instead gives \(\Delta h=0\), \(J=0\), and fertility
increases through consumption.

The matching current old tax must now be \(D=\Delta h/g\), not
\(\varepsilon/g\), where \(g=\gamma/(Kp)\) is that uncapped old donor's
housing response to resources. Its remaining funding requirement is
\[
R=\Delta x+(p-1/g)\varepsilon+
           [\chi+\kappa(p-1/g)]\Delta n.
\]
**This is a substantive restriction.** In the existing simple case
\(L=p\), \(\alpha(1+\omega_B)=\gamma\), one has
\(dx/d\varepsilon=1/g-p=p/\alpha\), hence
\[
R'=(\chi-\kappa p/\alpha)n'.
\]
It is negative whenever \(\chi/\kappa<p/\alpha\), although the original
fixed-fertility residual condition holds at equality. The old funding proof
therefore fails to establish a nonnegative residual tax. Its welfare
decomposition, now using total \(c',h'\), is
\[
W'=(\Lambda-m)c'+(\alpha/s-pm)h'+(m-m_F)R'.
\]
The final term need not be nonnegative. This is an obstruction to reusing
that proof, not a global impossibility result for endogenous-fertility transfers.

## 5. Dated claim versus equilibrium continuation

Corrected grants can establish the conditional recipient response and, subject
to a feasible residual allocation, date-zero housing/fiscal accounting.
They do not construct the announced equilibrium continuation. More births
change \(Y_1=\nu\bar n_0Y_0\); unchanged positive per-household choices
then add housing demand. Common property-tax rebates also change with cohort
mass. Future compensation must satisfy household optimization, housing
clearing, fiscal present-value balance, and demographic accounting, including
its feedback into today's promised resources. It cannot simply be asserted.

Given two admissible equilibrium paths, the finite bundle test and population
product require no convergence. A positive stationary-limit comparison does
require such limits to exist. Thus the main note may link fertility to the
corrected recipient construction above, but must preserve both qualifications:
its fixed-fertility welfare proof needs rechecking, and its unchanged-price
finite tail is not an endogenous-fertility transition theorem.
