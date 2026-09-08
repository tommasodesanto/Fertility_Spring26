# Assembled-note fertility review

September 8, 2026. Reviewed
`latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex`, including
Proposition 3, the complete fertility/path section, and Appendix
`Private fertility and the limits of the policy comparison`.
Read only; no build, simulation, or source edit.

**Verdict: mathematical pass subject to one explicit scope correction and
one type-level precision correction.** The restored-fertility funding
obstruction is now accurately stated. No equation needs changing.

## Required corrections

1. **Restrict the repayment derivative to uncapped young housing.**
   At lines 768–769, replace “the binding-finance calculation instead gives”
   with “on the branch with binding finance and uncapped young housing, the
   calculation instead gives.” The displayed derivative in terms of
   \(A,B_c,D_n,J_n\) differentiates two unconstrained choices \((h,n)\).
   It therefore does not apply when the young housing cap binds. In that
   regime \(h=H\), \(c=w-LH\), and
   \[
   \frac{dn}{db}=\frac{\chi}{x^2D_n}>0
   \]
   for a current unit transfer, regardless of the repayment rate while the
   binding-finance and binding-housing regime persists. This does not weaken
   the claim that repayment can reverse fertility on another branch.

2. **State the owner/renter ordering for each affected endowment type.**
   At lines 402–404, use: “A grant available only to buyers raises ownership
   at fixed prices; if \(n_{i0}^O\ge n_{i0}^R\) for every affected endowment
   type, both contributions are nonnegative.” The preceding finite identity
   is type-specific. An ordering of observed aggregate owner and renter
   fertility is insufficient because reform-induced tenure probabilities
   reweight types differently. The intended conditional ordering gives
   \[
   \Delta\bar n_i
   =\pi_{i1}^O\Delta n_i^O+
   (\pi_{i1}^O-\pi_{i0}^O)(n_{i0}^O-n_{i0}^R)\ge0.
   \]
   Conditional renter choices are unchanged under an owner-only grant.

## Claims that pass

**Current gifts and all cap regimes.** The constrained Hessian signs are
correct. With slack finance and binding young housing, the new formula
\[
c_w=\frac{\Gamma}{\Gamma+(1-\chi k_c)/x^2},\qquad n_w=k_c c_w
\]
correctly proves strict fertility growth because
\(0<\chi k_c<1\) and \(\Gamma>0\). The retained old log problem has
continuous marginal utility and strictly negative one-sided curvature in
every old cap/estate regime. The assembled text now handles old-regime kinks
with one-sided responses and finite monotonicity, avoiding an unjustified
everywhere-differentiable claim. Young housing grows strictly while uncapped.

**Finite comparison.** The fertility first-order condition uses total
consumption and housing, and remains valid under all the financial and
physical constraints. Its residual is strictly decreasing in fertility.
Consequently the displayed necessary-and-sufficient bundle test has the
correct inequality direction. The appendix correctly states that
\(c_1\le\chi n_0\) or \(h_1\le\kappa n_0\) implies \(n_1<n_0\),
including equality cases. The unchanged preference weight is explicit.

**Tenure and repayment.** The exact finite tenure identity and common-gift
selection derivative are correct. The common-gift selection term can be
negative; conditional gift monotonicity alone does not sign total fertility.
The repayment derivative is algebraically correct on its uncapped branch,
and its additional term can reverse the gift response.

**Connection to the actual transfer construction.** Holding old resources
fixed makes \(\lambda\) constant. Increasing adult space in the displayed
allocation raises adult consumption, and the fertility first-order condition
then gives the stated positive \(dn\). The corrected grants use
\(\Delta c=\Delta x+\chi\Delta n\) and
\(\Delta h=\varepsilon+\kappa\Delta n\). Their present values and the
matching old housing reduction produce exactly the displayed residual.
For complete explicitness one can calculate the intended fertility as
\[
n_\varepsilon=
\frac{\vartheta}{\chi/x_\varepsilon+
                  \alpha\kappa/s_\varepsilon},
\]
so the grants are predetermined numbers rather than payments contingent on
realized births. This extra display is optional: it already follows from
the stated fertility condition.

**Funding obstruction.** Under \(L=p\) and
\(\alpha(1+\omega_B)=\gamma\), the assembled appendix correctly obtains
\(R'=(\chi-\kappa p/\alpha)n'\). The original condition can hold at
equality while \(R'<0\). Thus the original proof of a nonnegative residual
tax fails, and its welfare decomposition loses a signed term. The text
correctly presents this as a failure to extend that construction's theorem,
rather than claiming that every free-fertility transfer policy is impossible.

**Dated and limiting paths.** The population product and old-cohort law
correctly order finite-date populations when fertility is ordered along
existing paths. They require no convergence. Positive stationary endpoints
imply replacement fertility in both economies and the stated inverse
housing-use ratio for adult population. The assembled note explicitly
requires an equilibrium path with fiscal feedback and declines to infer an
endogenous-fertility equilibrium theorem from the fixed-fertility finite
tail. Its distinction between a conditional recipient response, a complete
equilibrium continuation, and a stationary-limit comparison is sound.

Optional notation clarification: before the appendix proof, the reduction
can be explicitly dated by setting \(p=u_t\),
\(L=(1-\phi_t+q\tau^p)P_t\),
\(v=y_i^o+T_{t+1}\) including announced transfers, and
\(z=a'+P_{t+1}h+v\). The same proof then applies verbatim to a fixed
nonstationary price path. No additional mathematical restriction is needed.
