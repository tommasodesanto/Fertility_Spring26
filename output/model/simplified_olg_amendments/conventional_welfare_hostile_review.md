# Hostile review of the conventional-finance welfare argument

Independent bounded review, September 8, 2026. Sources read in full:
docs/model/simplified_olg_conventional_finance_proposal.md and
output/model/simplified_olg_amendments/conventional_efficiency_review.md.
Mandatory project context was checked. No model code, paper source, or other
review file was changed.

## Verdict

**PASS for the fixed-fertility household algebra. FIX before presenting the
displayed welfare proof as covering continuous endowment heterogeneity.
FIX the omitted individual financial implementation. BLOCK any claim that
these calculations alone establish a fully primitive equilibrium theorem or
characterize the global first best.**

I found no counterexample to the conditional lifetime Pareto construction once
the two repairs below are made. Both repairs stay within its existing
allocation and financing benchmark. Neither needs equal group masses, an
income atom, a numerical equilibrium, or changes to the private mortgage.

The result establishes that **one feasible Pareto improvement increases
current young housing and decreases current old housing**, with individual
fertility and tenure fixed. It does not establish a current housing
marginal-value advantage for young owners. That advantage has the opposite
sign in the stated \(q<\phi\) regime. It also does not characterize aggregate
young housing at every efficient allocation.

## 1. Fixed-fertility algebra: pass, with an exact scope distinction

Maintain the source notation
\[
 B=\beta(1+\gamma+\omega_B),\quad
 p=(1-q+q\tau^p)P,\quad L=(1-\phi+q\tau^p)P,\quad
 \delta=p-L>0.
\]
The proof requires positive utility weights, \(p>0\), \(L>0\), a genuinely
positive cash multiplier, and the stated strict young-cap, old-cap, and
old-estate-floor margins. The floor condition
\[
 \omega_Bp>q\gamma P
\]
is correct and independent of old wealth at stationary prices.

Fix the actual chosen \(n\), and use the source's adult cash \(a\), lifetime
adult resources \(M\), multiplier ratio \(r\in(0,1)\), adult consumption \(j\),
and \(\rho=p-r\delta\). The exact identity is
\[
 LM-pa
 =j\left[\frac{BL}{1-r}-\delta\right].
\]
Together with the independently checked housing identity this gives
\[
 s_F-s_C
 =\frac{\alpha r}{p(1+\alpha+B)\rho}(LM-pa).
\]
Therefore
\[
 h_F>h_C
 \iff LM>pa
 \iff Lqv>\delta(w-\chi n).
\]
This derivation does not evaluate the restricted objective at a potentially
infeasible \(s_F\); it also covers that boundary case directly.

The displayed primitive restriction
\[
 \boxed{\beta(1+\gamma+\omega_B)(1-\phi+q\tau^p)\ge\phi-q}
\]
is sufficient for this strict housing sign at every strictly constrained
allocation in the reduced regime. At equality the sign is still strict,
because \(r>0\). It is a sharp uniform threshold over unrestricted admissible
cash-multiplier ratios: if \(BL<\delta\), choosing
\(0<r<1-BL/\delta\) reverses the sign. It is **not** a necessary condition
for a particular household; the exact household test is \(LM>pa\).

The old-resource sign also passes:
\[
 S_C/M>B/(1+\alpha+B)=S_F/M.
\]
Under the slack-floor and slack-cap old solution, old consumption, housing,
and estate spending are all proportional to \(S\); all decline when young
finance is relaxed. A finite step toward the relaxed allocation preserves
positive log arguments, respects the lifetime budget, and strictly improves
utility. The source's explicit young-cap step bound suffices; old housing
falls, so its cap remains slack. The consumption fee
\(\eta=x_A(1-\exp(-g/2))\) leaves exactly \(g/2>0\) of the utility gain.

As a check on scope, the proposal's counterexample with actual \(n=1\) has
\(p=1,L=1/10,B=1,\alpha=31/40,M=10/3,a=11/10\).
Its **fixed-\(n\)** relaxed young housing is
\[
 h_F=1+\frac{310}{333}=\frac{643}{333}<2=h_C.
\]
Thus a binding ordinary mortgage alone does not deliver the desired
direction, even after fertility is held fixed. This is a household
counterexample at given prices, not an equilibrium counterexample.

## 2. Continuous heterogeneity: replace common bundles by uniform bounds

The source's Section 2 assumes positive masses sharing a common old bundle
and then multiplies a common selected-household increment by the selected
mass. A positive continuous endowment density generally gives zero mass to
an exact income type or exact bundle. Its finite-type statement therefore
does not by itself cover the requested heterogeneous economy.

The following is a direct repair. Let the household measures be atomless;
a continuous endowment density suffices. Current old donors and next-period
old recipients may have arbitrary, unequal positive masses \(m_D,m_R\).
Reserve a positive-measure recipient group before selecting households.
Every selected household \(i\) has its own finite step, gain \(g_i>0\), fee
\(\eta_i>0\), extra young housing \(H_i>0\), and reduction in old housing
\(K_i>0\). Its individual \(n_i\) is unchanged.

Restrict the available selected set to a positive-measure subset on which
\[
 \eta_i\ge\underline\eta>0,\qquad
 H_i\le\overline H<\infty,\qquad K_i\le\overline K<\infty.
\]
Such a subset exists whenever the eligible group has positive measure and
these quantities are finite and strictly positive: take countable level
sets. This is a measurable-subset argument, not an appeal to an unspecified
neighborhood of a numerical equilibrium.

Similarly restrict donor and recipient groups to positive-measure subsets
with baseline housing bounded below and consumption bounded above:
\[
 h_i\ge\underline h_D>0,\ c_i\le\overline c_D
 \quad(i\in D),\qquad
 h_i\ge\underline h_R>0,\ c_i\le\overline c_R
 \quad(i\in R).
\]
For recipients require the uniform margin
\[
 \min\{h_O^{\max}-h_i,\ e_i/P-h_i\}\ge\underline d_R>0.
\]
Again these uniform bounds follow on suitable positive-measure subsets from
the corresponding pointwise strict conditions. No globally bounded income
support is required. The source's slack old first-order conditions give
\(\gamma c_i/h_i=p\) individually, despite heterogeneous bundles.

Choose an atomless subset \(S_\epsilon\) of available selected households of
mass \(\epsilon>0\), and put
\[
 H_\epsilon=\int_{S_\epsilon}H_i\,d\mu,\quad
 K_\epsilon=\int_{S_\epsilon}K_i\,d\mu,\quad
 E_\epsilon=\int_{S_\epsilon}\eta_i\,d\mu.
\]
Take \(s_D=H_\epsilon/m_D\) housing from each donor and give
\(s_R=K_\epsilon/m_R\) housing to each recipient next period. Apply the
source's exact compensation formulas separately to each household:
\[
 D_i(s)=c_i[(1-s/h_i)^{-\gamma}-1],\qquad
 R_i(s)=c_i[1-(1+s/h_i)^{-\gamma}].
\]
These preserve each household's old utility with its own estate fixed.

For \(s_D\le\underline h_D/2\), valid uniform second-derivative bounds are
\[
 \overline A_D=
 \frac{\overline c_D\gamma(\gamma+1)2^{\gamma+2}}{\underline h_D^2},
 \qquad
 \overline A_R=
 \frac{\overline c_R\gamma(\gamma+1)}{\underline h_R^2}.
\]
Hence the total compensation cost above the linear service value obeys
\[
\begin{aligned}
 &\int_D D_i(s_D)d\mu-pH_\epsilon
 +q\left[pK_\epsilon-\int_R R_i(s_R)d\mu\right]\\
 &\quad\le
 \frac{\overline A_D}{2m_D}H_\epsilon^2+
 q\frac{\overline A_R}{2m_R}K_\epsilon^2
 \le C\epsilon^2,\\
 &C=\frac{\overline A_D\overline H^2}{2m_D}
    +q\frac{\overline A_R\overline K^2}{2m_R}.
\end{aligned}
\]
The selected fees satisfy \(E_\epsilon\ge\underline\eta\epsilon\).
If \(m_S>0\) is the available selected mass, an explicit sufficient bound is
\[
 0<\epsilon<
 \min\left\{
 m_S,\frac{\underline\eta}{2C},
 \frac{m_D\underline h_D}{2\overline H},
 \frac{m_R\underline d_R}{2\overline K}
 \right\}.
\]
It guarantees every physical and estate margin and a strictly negative
present-value resource cost, bounded above by
\(-\underline\eta\epsilon/2\).

This proves the heterogeneous extension conditional on positive-measure
eligible groups. Neither positive continuous density alone nor a single
strict type proves that those groups exist in equilibrium. A primitive
coverage result must put the relevant restrictions on a positive-measure
endowment set and establish the supporting price bounds. Full-support
ownership tastes then give a positive owner share where both tenure values
are finite. No equality of cohort masses or of donor/recipient masses is used.

## 3. Individual finance: an exact missing implementation

The source's aggregate goods budget is correct, but an aggregate
present-value surplus alone does not prove individual cash feasibility.
The following entries make its asserted implementation verifiable.

Let \(A=(1+q\tau^p)P=p+qP\). For selected household \(i\), let
\(\Delta z_i<0\) be its change in old total resources before the fee.
Preserve its original private mortgage and its original young bond position,
so its original net financial claim \(a'_{C,i}\) is unchanged.
Assign the following extra current and old transfers:
\[
 f_{0,i}=\Delta c_i-\eta_i+A\Delta h_i,\qquad
 f_{1,i}=\Delta z_i-P\Delta h_i<0.
\]
The first is net current support; the second is an enforceable old-age
payment to the program. They satisfy
\[
 f_{0,i}+qf_{1,i}=-\eta_i.
\]
Equivalently, lend the selected household
\[
 b_i^{\mathrm{bridge}}=q(P\Delta h_i-\Delta z_i)>0
\]
when young, require repayment \(b_i^{\mathrm{bridge}}/q\) on entering old
age, and collect the current fee \(\eta_i\). The loan covers precisely
\(\Delta c_i+A\Delta h_i\); the fee is financed by the specified reduction
in current consumption. This bridge crosses the young private financing
limit, which is exactly the benchmark's permission. It does not alter the
pre-existing mortgage or introduce borrowing initiated in old age.

On entering old age, this household has
\[
 a'_{C,i}+P(h_{C,i}+\Delta h_i)+v_i+f_{1,i}
 =z_{C,i}+\Delta z_i>0.
\]
Its new old bond saving is
\[
 a^e_{A,i}=q(e_{A,i}-Pk_{A,i})>0
\]
by the preserved strict estate floor. Thus its old consumption, housing
purchase, bond saving, and old-age payment fit the actual cash budget without
new old borrowing.

For donor \(i\), whose housing falls by \(s_D\), preserve its incoming
financial claims and original estate. Its old bond saving increases by
\(qPs_D\). The required extra current transfer is
\[
 f_{D,i}=D_i(s_D)-ps_D\ge0.
\]
For recipient \(i\), whose old housing increases by \(s_R\), old bond saving
decreases by \(qPs_R\), remaining positive by the uniform estate margin.
Its required extra old transfer is
\[
 f_{R,i}=ps_R-R_i(s_R)\ge0.
\]
These formulas follow from the cash budget with purchase price \(A\);
they are not assumed service-price trades. They show explicitly how title
payments and bond replacement generate the service-price terms.

The program's discounted budget is
\[
 \int_{S_\epsilon}(f_{0,i}+qf_{1,i})d\mu
 +\int_D f_{D,i}d\mu+q\int_R f_{R,i}d\mu
 \le C\epsilon^2-\underline\eta\epsilon<0.
\]
Its bridge finance and transfer commitments therefore have complete
repayment funding. The proposed benchmark must explicitly permit these
household-specific intertemporal arrangements. If only a uniform LTV change
or transfers financed by a common rebate were allowed, this proof would not
establish implementability under that narrower instrument set.

## 4. Goods, titles, estates, and future households: pass after clarification

With the heterogeneous integrals substituted, the dated additional goods uses
are exactly
\[
\begin{aligned}
 G_0&=\int_{S_\epsilon}(\Delta c_i-\eta_i)d\mu
          +\int_D D_i(s_D)d\mu,\\
 G_1&=\int_{S_\epsilon}\Delta c_i^2d\mu
          -\int_R R_i(s_R)d\mu,\\
 G_2&=\int_{S_\epsilon}\Delta e_i\,d\mu.
\end{aligned}
\]
Each selected lifetime identity and both housing balances imply
\[
 G_0+qG_1+q^2G_2\le C\epsilon^2-\underline\eta\epsilon<0.
\]
The consolidated incremental external balance at \(t+2\) is consequently
\(G_0/q^2+G_1/q+G_2<0\). A negative balance is a terminal surplus; it can
be disposed of or returned without lowering any utility. There is no
unpaid external debt or indefinite rollover. The program's own bridge
account can in fact settle at \(t+1\); selected old estate portfolios finish
at \(t+2\). These are distinct accounts.

At \(t\), extra selected young titles \(H_\epsilon\) exactly equal donor
title reductions. Donors replace \(PH_\epsilon\) of terminal title value
with bonds, keeping their estates fixed. At \(t+1\), the smaller title
supply from those donors' death sales is offset by the extra titles inherited
by the selected old households. Total current old occupancy is unchanged
because selected reductions \(K_\epsilon\) equal recipient additions.
At \(t+2\), selected and recipient death-title changes cancel. Future
entrants therefore face the same aggregate title supply and can execute
exactly the original purchases at the original prices.

Owner and renter occupancy totals separately remain fixed at every date,
so the rental intermediary, tax receipts, and common rebates can remain
unchanged. This does not require separate physical stocks by tenure.

The initial old retain their incoming claims and exactly the same estates;
donors receive exact utility compensation. Current young recipients have
unchanged young utility and unchanged old utility. Every selected young
household gains strictly after the fee. Everyone else can keep the baseline.
Each selected \(n_i\), every nonselected fertility choice, and all future
households' allocations remain fixed. Future births and cohort sizes are
therefore unchanged relative to the baseline.

This last conclusion uses the proposal's explicit convention that estates
are warm-glow spending, not entrants' inherited wealth. It would not carry
over unchanged to a model in which selected estates fund descendants or
give utility to additional named beneficiaries. No such change should be
silently imported.

The phrase “all other dated external claims are held fixed” should refer to
claims already inherited or contracted at the intervention date.
Counterfactual future bond purchases cannot all be fixed: donor and
recipient bond replacement is essential to the proof.

## 5. Remaining blockers and claims that must stay separate

1. **Primitive equilibrium coverage remains outside this proof.** The
   price-independent condition \(BL\ge\delta\) does not establish a positive
   constrained owner mass, cap margins, or stationary market clearing for
   an arbitrary housing stock and continuous endowment distribution.
   The source correctly acknowledges this obligation. It must be discharged
   by the separate analytical equilibrium construction before calling the
   combined result a fully primitive equilibrium theorem.
2. **No global first-best characterization has been proved.** A feasible
   Pareto improvement refutes Pareto efficiency of the baseline under the
   stated benchmark. It does not characterize Pareto weights, all optimal
   allocations, the global direction of housing reallocation, endogenous
   tenure at the optimum, or an efficient population.
3. **No current MRS-gap theorem is available here.** At the constrained
   young owner \(u_h^y/u_c^y=p-r\delta<p\), while an interior old owner has
   \(u_h^o/u_c^o=p\). The welfare gain comes from a complete lifetime
   adjustment that brings consumption and housing forward and lowers later
   consumption, housing, and estates. Calling this a larger current
   consumption-equivalent marginal value of young housing would be false.
4. **Fertility stays fixed in this theorem.** The exact comparison is valid
   conditional on the actual equilibrium \(n_i\), but the preferred fertility
   after a financing intervention is a separate comparative-static
   question. The Pareto construction establishes neither a birth increase
   nor the welfare effect of changing population.

## Verification

Independently reduced the two budgets and first-order conditions; checked
the fixed-\(n\) housing identity symbolically with an exact zero residual;
factored the derivative at the relaxed comparison as a separate cross-check;
checked the old-resource sign; and derived every individual transfer and
the aggregate terminal budget directly. The displayed counterexample uses
exact fractions from the proposal. No model solve, numerical equilibrium,
heavy computation, or literature claim was used.
