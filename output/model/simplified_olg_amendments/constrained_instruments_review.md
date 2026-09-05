# A finite market implementation with explicit extra instruments

**PROVED UNDER STATED ASSUMPTIONS.** An advance against current income, together with two household-specific housing incentives, can implement a Pareto improvement toward young owners at the original equilibrium prices. Households choose their allocations optimally. Both dated government budgets balance exactly. The intervention ends when the affected young household becomes old; its last changed estate expenditure occurs one date later.

This is a fallback result for an **expanded instrument set**, not a theorem about one-time unconditional cash gifts. It does not establish that this package is the smallest enlargement among all possible equilibrium mechanisms. Its advantage is a complete finite construction that includes private optimality, all title holders, future market clearing, and fiscal and external accounts.

The parallel investigation concerns a different enlargement: committed individualized cash transfers and later taxes. A result for that sequence should be called constrained Pareto inefficiency **relative to that committed multi-date transfer class**. It does not establish the same result for one-time gifts. Compulsory later taxes can provide government-enforceable claims on resources that a household cannot pledge privately, even when the nominal mortgage share is unchanged.

## 1. Benchmark and permissions

Take a complete competitive equilibrium of the proposal, including ultimate ownership of all initial housing titles. Keep fertility and every cohort mass fixed. Households optimize their remaining consumption, housing, saving, estate, and tenure choices conditional on those fertility levels. This is not a claim that an unrestricted endogenous fertility first-order condition remains satisfied after reform.

Choose a young owner \(i\) and an old owner \(j\), of equal positive mass for notation. The young household strictly prefers owning at its actual taste draw; its down-payment constraint has a positive multiplier, saving is positive, the physical owner cap is slack, and its subsequent old-age retention and estate constraints are slack. The old donor has slack retention and estate-composition constraints. These are eligibility conditions, not statements about every household or equilibrium.

For a continuum, use equal-mass matched groups with a uniform positive value gap, tenure margin, and physical slack. A zero-mass type does not imply positive aggregate reallocation.

Keep the single fixed housing stock, rental cap, owner cap, prohibition on additional old-age purchases, estate-composition constraints, original property tax, and original common rebate. Add only:

1. A government advance before closing, recovered from the recipient's current income after income arrives **at the same date**. Repayment is enforceable. There is no within-date interest in the model and no public loan remains next date.
2. A tax on the selected old donor's retained housing, paired with a predetermined lump sum.
3. A subsidy to the selected recipient's old-age housing one period later, paired with a predetermined lump-sum tax.

The government may target identified households and their baseline bundles; private information and incentive-compatible reporting are not modeled. No quantity is mandated. Housing incentives depend on actual housing; the lump sums are constants when households compare choices. A rebate tied to actual housing would cancel the marginal incentive and invalidate the proof. Government commitment is required for the announced next-date policy.

## 2. Exact allocation and welfare

Write \(u_t=q r_t\) and \(u_{t+1}=q r_{t+1}\) for unchanged user costs. Preserve the author's adult bundles \(x_i=c_i-\chi n_i\) and \(s_i=h_i-\kappa n_i\). Define \(B=1+\omega_B\) and \(J=1+\beta B\). The eligible baseline owners satisfy:

\[
q e_j=\omega_Bc_j^2,\qquad
\frac{\gamma c_j^2}{h_j^2}=u_t,\qquad
c_i^2=\frac{\beta x_i}{q},\qquad
q e_i=\omega_Bc_i^2,\qquad
\frac{\gamma c_i^2}{h_i^2}=u_{t+1}.
\]

These follow from the original old consumption-estate conditions, housing conditions, and young saving condition under the stated slackness.

Move \(\epsilon>0\) units of current housing from \(j\) to \(i\). Let \(C(\epsilon)\geq0\) be its real cost in date-\(t\) goods, with \(C(0)=0\) and finite \(C'(0+)=K\geq0\). Compensate the donor by scaling its consumption and estate together:

\[
z_j=\left(\frac{h_j^2}{h_j^2-\epsilon}\right)^{\gamma/B},
\qquad D(\epsilon)=Bc_j^2(z_j-1).
\]

The donor's new bundle is \((z_jc_j^2,h_j^2-\epsilon,z_je_j)\). Its utility change is exactly zero:
\(B\log z_j+\gamma\log(1-\epsilon/h_j^2)=0\).
The quantity \(D\) is its extra consumption-plus-estate expenditure in date-\(t\) value.

The recipient finances this expenditure and the real cost by the common scaling:

\[
z_i=1-\frac{D(\epsilon)+C(\epsilon)}{Jx_i}.
\]

Its target allocation is:

\[
c_i^*=\chi n_i+z_ix_i,\quad h_i^*=h_i+\epsilon,\quad
c_i^{2*}=z_ic_i^2,\quad h_i^{2*}=h_i^2,\quad e_i^*=z_ie_i.
\]

All other occupancy, at every date, is unchanged. Small enough \(\epsilon\) preserves positivity, the current owner cap, positive saving, and the future estate constraint \(z_ie_i>P_{t+2}h_i^2\).

The recipient's exact lifetime gain is:

\[
\Delta U_i=J\log z_i+\alpha\log(1+\epsilon/s_i).
\]

Since \(D'(0)=\gamma c_j^2/h_j^2=u_t\), its right derivative is:

\[
\Delta U_i'(0+)=
\frac{\alpha x_i/s_i-u_t-K}{x_i}.
\]

**Proposition.** Under the stated eligibility conditions and instrument permissions, if \(\alpha x_i/s_i-u_t>K\), the construction below is an equilibrium of the reformed fixed-fertility economy and a Pareto improvement for sufficiently small positive \(\epsilon\). Young housing increases, old housing decreases, and all other households and outside claimants retain their utilities or contractual returns.

This construction changes the affected households' chosen warm-glow estate spending and includes both changes in utility and dated resources. The source does not give these estates to separate heirs as entrant \(b_i\). If each estate amount is instead independently protected as an entitlement, this construction does not meet that stronger requirement.

## 3. Policies and private optimality

### Old donor at date \(t\)

Set the per-unit retained-housing tax and fixed lump sum to:

\[
\theta_{j,t}=
\frac{\gamma z_jc_j^2}{h_j^2-\epsilon}-u_t>0,\qquad
L_{j,t}=D(\epsilon)-u_t\epsilon+
\theta_{j,t}(h_j^2-\epsilon).
\]

Its modified original budget is:

\[
c^2+qe+(u_t+\theta_{j,t})h^2
=a_j+P_tH_j+T_t+L_{j,t}.
\]

The target satisfies this budget. Its consumption-estate ratio remains optimal; its housing condition is
\(\gamma z_jc_j^2/(h_j^2-\epsilon)=u_t+\theta_{j,t}\).
The objective is strictly concave and constraints linear, so these conditions give its global optimum. Its retention and estate-composition bounds relax.

The government's net payment to the donor is:

\[
N_{j,t}=L_{j,t}-\theta_{j,t}(h_j^2-\epsilon)
=D(\epsilon)-u_t\epsilon.
\]

This is positive for \(\epsilon>0\), because \(D\) is strictly convex with \(D(0)=0\) and \(D'(0)=u_t\). It is second order for small transfers. The compensated marginal incentive causes less housing; a cash gift alone would instead increase an interior old owner's demand at unchanged prices.

### Recipient when old at \(t+1\)

Set a negative housing tax and accompanying fixed lump-sum tax:

\[
\theta_{i,t+1}=(z_i-1)u_{t+1}<0,\qquad
L_{i,t+1}=\theta_{i,t+1}h_i^2<0.
\]

Here \(h_i^2\) in the lump-sum formula is the **baseline quantity**, a predetermined number. The subsidy is paid on actual housing, but \(L_{i,t+1}\) does not change with that choice. The effective price is \(z_i u_{t+1}>0\), which makes the unchanged old housing optimal:

\[
\frac{\gamma z_i c_i^2}{h_i^2}=u_{t+1}+\theta_{i,t+1}.
\]

The net government payment at the equilibrium quantity is zero:

\[
L_{i,t+1}-\theta_{i,t+1}h_i^2=0.
\]

### Recipient at date \(t\)

Advance the extra down payment:

\[
\ell_i=(1-\phi_t)P_t\epsilon.
\]

The nominal mortgage share remains \(\phi_t\), but the closing cash becomes \(b_i+\ell_i\), placing the cap exactly at \(h_i+\epsilon\). The net current transfer is:

\[
N_{i,t}=-D(\epsilon)-C(\epsilon)+u_t\epsilon
=-N_{j,t}-C(\epsilon)\leq0.
\]

Pay \(\ell_i\) before closing and collect \(\ell_i-N_{i,t}\) after income arrives. Positive current income supports this recovery for small enough \(\epsilon\). This explicitly makes current income pledgeable to government; it is not an unconditional gift.

The recipient chooses saving \(a_i'+\Delta a_i'\), where:

\[
\Delta a_i'
=\beta Bx_i(z_i-1)+(\phi_tP_t-qP_{t+1})\epsilon.
\]

Positive baseline saving makes the new value nonnegative for small transfers. Its original current budget holds because:

\[
(z_i-1)x_i+\Delta a_i'
+[(1-\phi_t)+q\tau^p]P_t\epsilon
=Jx_i(z_i-1)+u_t\epsilon=N_{i,t}.
\]

Its next-date financial wealth plus title liquidation resources change by:

\[
\Delta(a_{i,t+1}+P_{t+1}H_i)
=\frac{\Delta a_i'}q-\frac{\phi_tP_t\epsilon}q+P_{t+1}\epsilon
=Bc_i^2(z_i-1).
\]

Therefore its future target also satisfies the original old-owner budget with the announced tax and subsidy. The saving condition remains
\(z_i c_i^2=\beta z_i x_i/q\), and the future consumption-estate and housing conditions hold.

Its current financing cap binds with a positive multiplier, since for small enough \(\epsilon\):

\[
\frac{\alpha z_i x_i}{s_i+\epsilon}>u_t.
\]

Conditional on fixed fertility and owning, the full two-budget problem has a concave log objective and linear constraints. These first-order conditions and the cap multiplier therefore establish global conditional optimality, not just budget feasibility. The strict baseline ownership-taste margin preserves the optimal tenure for a sufficiently small policy. No other household's price or policy changes, so every other household remains at its original optimum.

## 4. Markets, title holders, and finance

**Current housing.** Demand rises by \(\epsilon\) for the young recipient and falls by \(\epsilon\) for the old donor. The single stock clears. Both remain owners, so rental holdings and tenure segmentation are unchanged.

**Future housing.** At \(t+1\), the recipient carries \(h_i+\epsilon\), keeps its original \(h_i^2\), and sells the extra title when the donor's estate would have sold it. Every subsequent occupancy is unchanged. The future subsidy is essential to make that occupancy optimal despite lower goods consumption.

**All initial title owners.** Every house price remains unchanged. Passive owners of residual initial titles receive identical sale values, and rental intermediaries receive identical rental and resale flows. Thus all their contractual returns and welfare remain unchanged, regardless of which outside claimant owns those titles. Their ownership must still appear in the complete initial balance sheet; their price exposure is zero here.

**Original property-tax budget.** The price path, stock, property tax, and common rebates stay fixed, so the original fiscal identity still holds.

**New date-\(t\) budget.** The recipient's current net tax finances the donor's net compensation and real cost:

\[
N_{i,t}+L_{j,t}
-\theta_{j,t}(h_j^2-\epsilon)+C(\epsilon)=0.
\]

A lender may provide the government temporary cash before closing and receive the principal back after current income and the new levies arrive. There is no new interperiod government debt. The instrument requires this within-date financing and enforceable recovery, not free resources.

**New date-\(t+1\) budget.** The recipient's lump-sum tax exactly finances its housing subsidy at the equilibrium quantity. No new government payments occur at \(t+2\) or later. The announced lumps remain constant when households compare deviations; fiscal balance is evaluated at equilibrium choices.

**External dated resources.** Count estate spending when it is actually delivered, and do not count domestic title payments as goods costs. The only changes in aggregate goods use are:

\[
\begin{aligned}
\Delta G_t&=(z_i-1)x_i+(z_j-1)c_j^2+C(\epsilon),\\
\Delta G_{t+1}&=(z_j-1)e_j+(z_i-1)c_i^2,\\
\Delta G_{t+2}&=(z_i-1)e_i.
\end{aligned}
\]

They satisfy:

\[
\Delta G_t+q\Delta G_{t+1}+q^2\Delta G_{t+2}
=Jx_i(z_i-1)+Bc_j^2(z_j-1)+C(\epsilon)=0.
\]

External bonds finance these dated differences at the original world terms. Original household budgets hold and every lender is repaid. No outside resources are donated; the incremental net external position returns to zero after the last changed estate expenditure.

The recipient's different old state exists only at \(t+1\). The next young cohort's income, entrant wealth, choices, and mass remain unchanged, because estates do not feed entrant wealth in this model. From \(t+2\), the living-household distribution, housing allocation, prices, and policies are the original continuation. Only the already-accounted estate spending at \(t+2\) differs. After that expenditure, the full flow allocation and net external position return to baseline.

## 5. Minimality within this route and remaining limits

Each added control has a distinct role **within this fixed-price finite-support route**:

- An interior old owner with an unchanged housing price and only a lump-sum wealth change scales housing, consumption, and estate together. Reducing its housing through that channel makes it worse off. A donor housing wedge or another substitution instrument is needed to reduce its housing at fixed prices while preserving utility.
- Without a future marginal incentive, unchanged old housing would force recipient old consumption unchanged by \(h_i^2=\gamma c_i^2/u_{t+1}\). Its estate ratio and saving condition would then force its other nonhousing utility arguments unchanged. The future subsidy permits lower goods spending while holding occupancy fixed. Otherwise one must allow future markets to adjust or use a different marginal instrument.
- The recipient's net current transfer is nonpositive, but its binding closing constraint needs positive extra cash. Separating a preclosing advance from after-income collection supplies that missing liquidity. A gift paired with a compulsory later tax is also a fiscal credit arrangement, not the original one-time gift class.

These arguments do not exclude a simpler policy that changes prices and many future allocations. The parallel committed-cash-sequence construction is such an alternative; its transfer dates, enforceability, and commitment should be stated as part of its benchmark.

The theorem also requires an eligible positive-mass pair and a value gap above cost. It uses personalized instruments, not a uniform tax and common rebate. It protects household lifetime welfare under warm-glow estates, not separate heirs in a different inheritance model. A positive fixed moving cost needs the exact finite inequality
\(J\log z_i+\alpha\log(1+\epsilon/s_i)>0\) at a feasible transfer. Local existence of this reformed equilibrium does not imply uniqueness of all equilibria of the policy.

## 6. Three checks

**Derivation.** Joint donor consumption-estate compensation preserves the original consumption-estate condition. The recipient's common nonhousing scaling preserves its saving and estate conditions. The two housing wedges supply the remaining housing conditions. Exact welfare and all budget equalities then establish the result.

**Original constraints and accounts.** The checks separately cover the young's two budgets, both old-owner problems, mortgage repayment, nonnegative saving, the modified closing cap, original retention and estate composition, the single housing stock, rental intermediaries, outside initial title owners, property rebates, both dated government budgets, and all three dates of real goods. The only changed permissions are those listed in Section 1.

**Adversarial and independent calculation.** An in-memory Python check directly optimized the recipient's original two-budget owner problem with SLSQP, retaining the future retention and estate inequalities. Twelve constructed-household cases covered financed shares \(0.4,0.8\), transfer sizes \(10^{-3},10^{-5}\), and marginal costs \(0,0.05,0.3\). The last cost exceeds the constructed baseline value gap \(0.25\) and correctly produces a welfare loss.

The largest original-budget, fiscal, or external-account residual was \(4.44\times10^{-16}\). The largest difference between six independently optimized choices and the analytical target was \(1.10\times10^{-7}\). At \(\epsilon=10^{-5}\), young welfare slopes were \(0.24998669,0.19998651,-0.05001465\), converging to \(0.25,0.20,-0.05\). Seller utility was unchanged to numerical precision.

Inputs were \(q=0.5\), \(P_t=P_{t+1}=P_{t+2}=1\), \(\tau^p=0.1\), \(\alpha=\beta=0.4\), \(\gamma=0.3\), \(\omega_B=0.4\), \(n_i=0.5\), \(\kappa=0.5\), \(\chi=0.15\), \(x_i=1\), \(s_i=0.5\), \(h_j^2=1\), and \(C(\epsilon)=K\epsilon+0.1\epsilon^2\). Other bundles were derived from the displayed original conditions and budgets. An initial checker rejected a financed-share variant whose saving was accidentally held fixed during baseline reconstruction; saving was then derived from the original old budget before every case. No theorem formula changed to pass that gate.

These are checks on constructed conditional households, not a numerical full-economy equilibrium or calibration. The equilibrium theorem instead starts from any complete baseline equilibrium containing an eligible pair.

## 7. Decision supported by this pass

There is a proved finite market implementation after explicitly enlarging the instrument set. If the author prefers transfers alone, the parallel multi-date cash-sequence theorem should be assessed under its own commitment and enforcement assumptions. If that class fails, the specific missing margin can be supplied by the compensated housing incentives and current-income advance above. These are different constrained planner problems and should not share an unqualified theorem label.

Only this report was written. No source, model code, calibration, ledger, memory, or other worker's output was edited.
