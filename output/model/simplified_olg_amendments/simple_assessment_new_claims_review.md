# September 5: independent review of the new simple-model claims

Read-only reviewer report on renter replication, the constructed equilibrium
family and slack-rental example, and stationary stock scaling. The lead
accepted its clarifications about the varying primitives across the family
and fixed adult/child counting conventions. The lead also corrected the scaling display to hold the normalized old-state
distribution G fixed, as defined in the original specification’s
\(\eqref{eq:amend_finite_memory}\); only the corresponding unnormalized
cohort measure would scale with population. The economic scaling conclusion
is unchanged. The report is otherwise preserved as returned. Equation labels remain useful; line numbers refer to the reviewed
checkpoint. The separate negative-fertility original-problem optimization was
added by the lead and is not part of this reviewer’s claimed coverage.

## 1. Renter replication

**Proved, conditional on the stated feasibility conditions.** From the original budgets [proposal, \(\eqref{eq:amend_young_renter}\)–\(\eqref{eq:amend_young_owner}\), lines 99–122] and user-cost identity \(\eqref{eq:amend_usercost}\), set
\[
a_R'=a_O'+qP_{t+1}h-\phi P_th.
\]
Then
\[
a_R'+qr_th
=a_O'+[(1-\phi)+q\tau^p]P_th,
\]
so the young renter exactly reproduces the owner’s \(c,h,n\) and budget. At \(t+1\),
\[
R_fa_R'=R_fa_O'-R_f\phi P_th+P_{t+1}h,
\]
which is exactly the owner’s old-age financial wealth plus title-sale value. Hence the renter can reproduce the same \(c^2,h^2,e\) in \(\eqref{eq:amend_old_renter}\), provided the rental menu contains both \(h\) and \(h^2\).

The calculation is nonstationary; it uses \(P_{t+1}\), not \(P_t=P_{t+1}\). It requires: no capital-gains tax; \(a_R'\geq0\); rental caps allowing the replicated current and old housing; and no ownership taste payoff. The last condition is indispensable: otherwise matching material allocations does not match lifetime value.

If \(MV^Y>u_t\) and there is slack above the replicated current rental size, hold \(n,a_R'\), and continuation choices fixed, set \(h\mapsto h+\epsilon\), \(x\mapsto x-u_t\epsilon\). The derivative is
\[
\frac{\alpha}{s}-\frac{u_t}{x}
=\frac{MV^Y-u_t}{x}>0.
\]
Thus a small rental deviation strictly improves utility. Assessment equation \(\eqref{eq:simple_replication}\) and its surrounding language are correct [assessment lines 156–170].

I also ran the read-only verifier’s replication check. Its fixed-price witness has budget errors below \(7\times10^{-16}\) and positive gains; this is supporting evidence, not the proof [JSON lines 133–137, 405–409; verifier lines 204–223].

## 2. Appendix A mixed-tenure constructions

**First construction: proved as a parametric family of equilibria, not a tax comparative static.** With \(x_O=h_O=1,n_O=.75\),
\[
\vartheta=.75\left(.15+\frac{.2}{.625}\right)=.3525,\quad
\rho=1+.4(1+.3+.4)=1.68.
\]
For \(u=(1+\tau^p)/2\), the old-owner FOCs give
\[
(c^2,h^2,e)=\left(.8,\frac{.24}{u},.64\right),
\]
and \(MV^O=.3(.8)/(.24/u)=u\), while \(MV^Y=.4/.625=.64>u\) for \(0\leq\tau^p<.28\). Retention and estate inequalities are strict: \(h^2\leq.48<1\) and \(e=.64>h^2\). The purchase-cash restriction binds exactly: \((1-\phi)Ph=.2=b\). Saving is
\[
a_O'=.98-\tfrac12T>0.
\]

The renter equation \(\eqref{eq:simple_witness_renter}\) has exactly one root in \((0,.5)\): its left side decreases from infinity and its right side increases to infinity. Its displayed resource expression implies \(x_R>1.30\). Therefore both rental-cap multipliers are strictly positive: at the young cap \(\alpha x_R/(h_R-\kappa n_R)>u\), and at the old cap \(\gamma c_R^2/h_R^2>u\). These are not merely feasible cap choices.

\[
\pi^O=\frac{.5-n_R}{.75-n_R}\in(0,1)
\]
makes mean fertility \(1/\nu=.5\); \(Y=O=2/d_h\) clears housing, and \(T=q\tau^p d_h/2\) satisfies \(2YT=q\tau^p P\bar H\). The logistic-location choice implements \(\pi^O\). The dated constraints include the original purchase, mortgage repayment, retention, estate, saving, and cap inequalities [verifier lines 89–114]. Concavity of the log objective and convexity of this feasible set make the verified KKT solutions global optima; read-only independent SLSQP checks match choices within \(8.1\times10^{-7}\) across the tested witnesses [verifier lines 117–164].

**Second construction: proved.** At \(\tau^p=0\), \(n_O=.49\),
\[
\vartheta=.49\left(.15+\frac{.2}{.755}\right)
=\frac{61397}{302000},
\]
and \(y+b=2.2535=w\). The unrestricted renter solution is
\[
n_R=\frac{276716279}{551645600}=.50161966>.5>.49,\qquad
\pi^O=.139389713>0,
\]
with \(h_R=1.04036834\), \(h_R^2=.47373511\), both below \(1.5\). The owner wedge remains strict:
\[
MV_O^Y=.52980132>.5=MV_O^O.
\]
Thus a positive owner mass exists with slack rental caps. Removing the redundant \(1.5\) cap changes nothing in this *second* economy.

The note correctly prevents the false causal reading: the two columns differ in \(\tau^p,\vartheta,y,n_O\), and equilibrium composition, not merely the rental cap [assessment lines 455–475]. The only qualification is terminological: “for every \(\tau^p\)” means a family with \(y,T,\bar\xi,\pi^O\) adjusted by \(\tau^p\), not one economy held fixed while varying the tax.

## 3. Stationary stock scaling and transition language

**Proved under the actual closure.** Given a positive stationary equilibrium, scale
\[
(\bar H,Y,O,G)\mapsto(a\bar H,aY,aO,G),
\]
while retaining all per-household distributions, choices, \(P,T,q\), and tenure shares. Individual problems are unchanged because the entrant distribution and world bond price are fixed. Housing clearing and the rebate budget both multiply by \(a\):
\[
a(Y\bar h^Y+O\bar h^O)=a\bar H,\qquad
a(Y+O)T=q\tau^pP(a\bar H).
\]
The demographic recursion is homogeneous, so stationary replacement and \(Y=O\) remain valid. This proves assessment lines 266–275 and \(\eqref{eq:simple_population_stock}\). It fails if stock expansion changes entrant income, supplier/land income, fiscal transfers beyond the ordinary rebate, domestic asset-market clearing, or \(q\).

The figure’s qualification is adequate: it explicitly calls fertility prescribed and says it is neither a solved GE transition nor a speed prediction [assessment lines 279–285; verifier lines 256–266]. It should not be reported as an equilibrium transition result.

One small precision improvement is advisable for the person-count sentence (assessment lines 287–290): state that adult counts, child conversion, and children’s timing/residence convention are fixed across steady states. Under those conditions person counts scale with household counts; otherwise the ranking need not carry over.