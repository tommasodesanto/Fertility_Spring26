# Mortgage timing in the two-age housing model

Independent review, September 7–8, 2026. This is a proposed-model audit, not an adopted specification or an equilibrium inefficiency theorem. The protected manuscript and model code were not edited.

The comparison sources are [the conventional-finance proposal](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/simplified_olg_conventional_finance_proposal.md) and [the original analytical model](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/oracle_analytic_efficiency_model.tex).

## Main finding

The condition \(qP_{t+1}\geq\phi P_t\) is tied to **capitalizing every mortgage interest payment until old age**. A loan of \(d\) then becomes a liability of \(d/q\). With \(q=1/2\), a mortgage initially financing 80% of a constant-price house requires a retirement payment equal to 160% of its purchase price. This is the financing contract in the current conventional-finance proposal; it is not an amortizing mortgage.

Servicing interest during working age changes the relevant terminal liability. An interest-only mortgage leaves principal \(d\) due at retirement; a partially amortizing mortgage leaves less. Correct accounting restores the sufficient housing-wedge condition at constant prices. But the earlier payments must be paid from young resources. The new restriction is a genuine restriction on when debt can be serviced, not a relabeling that creates equity for free.

There are two viable simple variants:

1. **Smallest algebraic change:** keep all young resources available initially and service interest during young age. Only the interpretation of the mortgage and one net-asset bound change. The initial down payment remains 20%; the household must also fund the working-age interest bill. This is an interest-only mortgage, not a fully amortizing fixed-payment mortgage.
2. **Ordinary amortization with a genuine closing constraint:** distinguish the paycheck already received from earnings arriving during the remaining working years. With level income, consumption, and mortgage payments, all intermediate financial constraints reduce exactly to a closing check and a retirement check. There are still only two utility ages and one housing/fertility decision. Closing cash is \(b+y_{\mathrm{current}}\), not \(b\) alone.

Neither variant proves that the required constraint pattern occurs in a market-clearing equilibrium. Neither makes the mortgage constraint bind for arbitrarily patient households while holding all endowments fixed.

## 1. A general repayment identity

Let \(d\) denote the principal advanced at purchase, with
\[
0\leq d\leq\phi P_t h.
\]
Keep the existing property-tax reserve and write \(\zeta_t=q\tau^p\), so purchasing and reserving the tax costs \((1+\zeta_t)P_t h\).

Within the young age, let \(Q_j\) be the discount factor to payment date \(j\), with \(Q_J=q\) at retirement. The mortgage makes scheduled payments \(m_jd\) during young age and leaves principal \(\rho d\) after those payments, immediately before the old-age problem begins. Payments classified as young-age payments must precede access to old income. Competitive lending at the saving rate requires
\[
\boxed{\quad A+q\rho=1,\qquad A=\sum_{j=1}^{J}Q_jm_j.\quad}
\tag{1}
\]
Here \(A\) is the present value of payments made while young per dollar originally borrowed; \(\rho\) is the remaining principal per dollar borrowed. They are not independent parameters.

Three contracts illustrate the distinction:

| Contract | Young payment PV \(A\) | Remaining balance \(\rho d\) |
|---|---:|---:|
| All interest capitalized until old age | \(0\) | \(d/q\) |
| Interest serviced during young age | \(1-q\) | \(d\) |
| Fully repaid before old income arrives | \(1\) | \(0\) |

A normal amortization schedule has declining principal after each regular payment. The interest and principal portions must be distinguished: interest service is a financing cost; principal repayment builds equity. [CFPB, mortgage amortization](https://www.consumerfinance.gov/ask-cfpb/how-does-paying-down-a-mortgage-work-en-1943/).

For a level-payment loan with gross payment-period rate \(R\), original term \(N\), and \(J\leq N\) payments before old age,
\[
q=R^{-J},\qquad
m=\frac{R-1}{1-R^{-N}},\qquad
\rho_J=\frac{1-R^{-(N-J)}}{1-R^{-N}},\qquad
A_J=1-q\rho_J.
\tag{2}
\]
Thus \(0\leq\rho_J\leq1\), including zero when the mortgage is fully repaid. These formulas derive from \(D_j=RD_{j-1}-md\), with \(D_0=d\). No house-price appreciation is used.

## 2. Exact aggregation when all young resources are initially available

Let \(w=b+y^y+T_t\) be the proposal's young resources, all available at purchase. Let \(k'\geq0\) be gross financial assets entering old age, after making the scheduled young-age payments. Define the author's next-age **net** financial wealth by
\[
a'=k'-\rho d.
\]
Using present-value nondurable expenditure \(c\), the gross budget is
\[
c+qk'+(1+\zeta_t)P_t h+Ad=w+d.
\tag{3}
\]
Substitute \(k'=a'+\rho d\) and use (1). The mortgage terms cancel from the resource identity, giving
\[
\boxed{\quad c+qa'+(1+\zeta_t)P_t h=w,\qquad
qa'+q\rho\phi P_t h\geq0.\quad}
\tag{4}
\]
The budget is unchanged; the feasible net-asset floor changes. Keeping the old floor \(qa'+\phi P_t h\geq0\) after adding amortization is generally incorrect.

This equivalence goes both ways. For \(\rho>0\), any allocation satisfying (4) can be implemented with
\[
d=\max\{0,-a'/\rho\}\leq\phi P_t h,\qquad k'=a'+\rho d\geq0.
\]
The household can reserve \(Ad\) for future coupons, reserve the nondurable expenditure, and invest \(qk'\) for old age. All reserves are nonnegative. With all \(w\) initially available, the budget pays for every reserve and the house. For \(\rho=0\), (4) requires \(a'\geq0\), and the same allocation can be implemented without a mortgage. There is no separately available unsecured loan in this construction.

Let \(v=y^o+T_{t+1}\), \(z=a'+P_{t+1}h+v\), and retain the no-arbitrage service price
\[
p=(1+\zeta_t)P_t-qP_{t+1}.
\]
Equation (4) becomes
\[
c+ph+qz=w+qv,\qquad c+L_\rho h\leq w,
\quad L_\rho=(1+\zeta_t-q\rho\phi)P_t.
\tag{5}
\]
The coefficient \(L_\rho\) is the **total cash burden over young age per unit of housing**, including the initial equity payment, the young-age repayment PV, and the tax reserve:
\[
L_\rho P_t^{-1}=1-\phi+A\phi+\zeta_t.
\tag{6}
\]
It is not the down-payment ratio. At \(P=100\), \(\phi=.8\), \(q=.5\), and zero tax:

| Contract | Actual initial down payment | PV of young payments | Aggregate young cash burden |
|---|---:|---:|---:|
| Capitalized-interest balloon | 20 | 0 | 20 |
| Interest serviced, principal repaid when old | 20 | 40 | 60 |
| Fully amortized while young | 20 | 80 | 100 |

The interest-servicing case does not impose a 60% initial equity requirement. It does require a household whose entire young endowment is already available to cover both 20 of initial equity and 40 of interest in present value. Calling all 60 a down payment would be misleading. If \(y^y\) means only an annual paycheck rather than the young period's endowment, this budget needs an explicit income-timing interpretation before it is used quantitatively.

### The minimal interest-servicing variant

Set \(\rho=1\), so \(a'=k'-d\). The exact constraint is
\[
\boxed{\quad a'+\phi P_t h\geq0.\quad}
\tag{7}
\]
It differs from the current proposal's \(qa'+\phi P_t h\geq0\). The new bound is supported by actual young-age coupons, not a lower interest rate. For example, continuous coupons at rate \(r\) over young-age length \(T\), with \(q=e^{-rT}\), have PV
\[
\int_0^T e^{-rs}rd\,ds=(1-q)d.
\]
Principal \(d\) remains due on entering old age. Coupons paid strictly before old income arrives cannot be financed by that income. If all interest instead becomes due when old income arrives, and there is no earlier servicing or reserve requirement, the liability is again \(d/q\); simply calling its components interest and principal changes nothing.

### The full-amortization limitation

If \(\rho=0\), (4) and (5) are independent of \(\phi\). With all young income pooled at purchase, a fully repaid mortgage only rearranges gross payments within young age. It cannot transfer resources from old income to young spending. Housing may still be restricted by the inability to borrow against old income, but there is no remaining origination-LTV comparative static in this aggregated feasible set. A claimed \(\phi\) effect would need another explicitly specified timing restriction.

## 3. Actual income timing, without a quantitative life-cycle model

If earnings arrive during young age, their present value cannot all be used at closing without borrowing against future wages. An aggregate budget alone then loses a real financial restriction.

Let \(W_j=b+\sum_{s=0}^{j}Q_sy_s\), including actual rebates when received; \(C_j\) is cumulative discounted nondurable expenditure, including child goods; and \(A_jd\) is cumulative discounted mortgage payments. If the existing tax reserve is posted at purchase, no unsecured borrowing and no refinancing require, at every payment date,
\[
C_j+(1+\zeta_t)P_t h+A_jd\leq W_j+d.
\tag{8}
\]
The closing restriction uses \(W_0=b+y_0\), where \(y_0\) is the paycheck already received. It never excludes that paycheck. A restriction using \(b\) alone is valid only if the closing actually precedes income receipt, or if \(b\) is redefined to include that income. Neither convention may be imposed silently.

The terminal version of (8) gives (4) with total young income PV in \(w\). Earlier versions need not follow from it: two income paths can have the same PV and different cash at closing. A single two-period PV-income budget cannot generally distinguish those paths.

### An exact two-checkpoint construction

There is a particularly simple way to retain ordinary amortization. Keep just one housing/fertility choice and two utility ages. At closing the household receives \(y_0\), so liquid cash is
\[
w_0=b+y_0.
\]
During the remaining young age, it receives a constant income payment \(y_w\), purchases a constant nondurable flow \(C\), and makes the constant mortgage payment \(md\) at dates \(1,\ldots,J\). Define
\[
\Lambda_j=\sum_{s=1}^j R^{-s},\quad \Lambda=\Lambda_J,
\quad c=\Lambda C,\quad w=w_0+\Lambda y_w.
\]
This defines \(c\) as the PV of the consumption composite entering the existing young utility; it adds no new within-young consumption choices. The same fixed flow schedule can include the expenditure on children. The house provides the existing young housing bundle throughout.

Let \(\ell_j\geq0\) be liquid financial wealth after date \(j\)'s payments. Then
\[
\ell_0=w_0+d-(1+\zeta_t)P_t h,
\qquad
R^{-j}\ell_j=\ell_0+(y_w-C-md)\Lambda_j.
\tag{9}
\]
The right side is a convex combination of its initial and terminal values. Consequently **all** intermediate \(\ell_j\geq0\) restrictions are equivalent to \(\ell_0\geq0\) and \(\ell_J\geq0\). With mortgage principal \(d\leq\phi P_th\), this reduces exactly to
\[
\boxed{
\begin{aligned}
c+qa'+(1+\zeta_t)P_t h&=w_0+\Lambda y_w,\\
a'+\rho\phi P_t h&\geq0,\\
(1-\phi+\zeta_t)P_t h&\leq w_0=b+y_0.
\end{aligned}}
\tag{10}
\]
To verify sufficiency, choose \(d=\phi P_th\), so the first inequality in (9) follows from the closing check and the terminal one from \(\ell_J=a'+\rho d\geq0\). Equation (9) then verifies every intermediate date. Gross borrowing and saving can coexist because their returns are equal, as in the existing model.

The old endowment \(v\) is added only after these working-age cash constraints. It remains known when housing is chosen and can have the same heterogeneity as in the proposal. It has not been set to zero or required to finance an arbitrary share of the loan.

This construction adds a distinction between current cash and subsequent young earnings. It does not require a new independent income parameter: one may specify \(y_w=y_0\), representing a known level earnings stream, if economically appropriate. That equality is an explicit income-profile choice, not an accounting identity. A growing or irregular income path generally requires checking its intermediate constraints rather than assuming (9).

If consumption must also be purchased at the instant of closing, add its actual initial payment to the final inequality in (10). Its exclusion here is justified by closing preceding the subsequent flow purchases, while the current paycheck has already arrived. The timing is therefore different from the original model's exclusion of all current income.

For a fully amortized mortgage, \(\rho=0\), (10) still contains \(\phi\) through the genuine closing condition. This is the cleanest way to preserve the 20% initial payment mechanism with full amortization. It is a transparent timing extension to the proposal's assumption that all young income is already available at purchase.

## 4. The housing sign at constant prices

For (5), the condition in the proposal becomes
\[
\boxed{\quad L_\rho\geq p
\iff P_{t+1}\geq\rho\phi P_t.\quad}
\tag{11}
\]
With constant prices and an interest-servicing or amortizing mortgage, \(0\leq\rho\leq1\) and \(\phi<1\), the inequality is strict for every \(q\in(0,1)\). It does not require appreciation or an artificially high generational bond price.

By contrast, capitalized interest has \(\rho=1/q\), so (11) is exactly the old \(q\geq\phi\) condition. Mortgage design, rather than a failure of present-value accounting, explains the difference.

These marginal-value statements use the proposal's old-owner resizing problem, so continuation utility depends on total resources \(z\). If the original restriction \(h^2\leq h\) is retained and binds, the continuation value also depends separately on the inherited home. That restriction cannot be removed by correcting mortgage accounting.

Let \(\lambda>0\) and \(\mu\geq0\) be the lifetime-budget and aggregate-young-cash multipliers. With a slack young housing cap, the young household's consumption-valued marginal benefit of housing is
\[
MV_Y\equiv\frac{\partial u^y/\partial h}{\partial u^y/\partial c}
=\frac{\lambda p+\mu L_\rho}{\lambda+\mu}
=p+\frac{\mu}{\lambda+\mu}(L_\rho-p).
\tag{12}
\]
Thus strict financial restriction and (11) generate \(MV_Y>p\), irrespective of the numerical value of \(\beta\). In the uncapped homogeneous regime, the proposal's proof also yields \(h<h^*\), relative to removal of that financial restriction. When \(\rho>0\) is fixed, its positive local credit-response proof applies after replacing the effective financed share by \(q\rho\phi\): \(\partial h/\partial\phi>0\). These are conditional-owner, fixed-price results.

Under (10), let \(\eta\geq0\) be the additional closing multiplier and \(L_0=(1-\phi+\zeta_t)P_t\). Because the closing payment excludes subsequently purchased nondurables, (12) becomes
\[
MV_Y-p=\frac{\mu(L_\rho-p)+\eta L_0}{\lambda+\mu}.
\tag{13}
\]
Both financial channels can therefore produce a positive housing wedge at constant prices. Full amortization removes the terminal \(\phi\) channel but preserves the closing channel.

An old owner's housing marginal value equals \(p\) only when its estate floor and size cap are slack. Otherwise those restrictions can make the old marginal value exceed \(p\). In the homogeneous old-age problem the estate floor is strictly slack when
\[
\omega_Bp_{t+1}>q\gamma P_{t+2}.
\tag{14}
\]
The old housing choice must also be below its cap. A young wedge above \(p\) alone does not prove a young-old MRS gap if these old restrictions bind. Likewise, with the original inherited-home restriction binding, \(\beta V_H\) enters the young housing first-order condition and subtracts \(\beta V_H/(\lambda+\mu)\) from the right side of (12). An equilibrium theorem must state and establish these branches, not infer them from loan accounting.

### What high patience does and does not permit

Retain the proposal's \(B=\beta(1+\gamma+\omega_B)\), \(D=1+\alpha+\vartheta+B\), and slack relevant housing caps. Strict financial restriction in (5) is characterized by
\[
\frac{w+qv}{w}>
\frac{D}
{1+\alpha L_\rho/p+
\vartheta(\chi+L_\rho\kappa)/(\chi+p\kappa)}.
\tag{15}
\]
The \(\kappa\) in (15) is the existing child-space requirement; \(\rho\) is the mortgage residual, a different object.

No upper bound on \(\beta\) is needed for the conditional sign in (12). But increasing \(\beta\) while fixing \(w,v\) eventually makes the unconstrained household devote very little to young consumption and housing, so the aggregate-young cash restriction becomes slack. For any given finite \(\beta\), a sufficiently high future-resource/current-resource ratio can satisfy (15), subject to the stated cap checks. This is an explicit income-timing restriction, not a universal claim about high-patience households. The closing construction allows small initial liquidity alongside substantial working-age income; its binding test compares the desired house with \(w_0/L_0\).

## 5. Changes that are not equivalent

- **Reducing the final balance without paying earlier installments:** violates (1), unless an interest subsidy or lender loss is introduced. Changing the lower bound is valid only with the corresponding contract interpretation in (3).
- **Charging young repayments twice:** once \(a'=k'-\rho d\) is used, the reduced budget is (4). Subtracting \(Ad\) again from its right side double-counts the repayments already embodied in the tighter asset floor and the gross-to-net substitution.
- **Making interest an old-age deduction while claiming it was serviced when young:** if old resources are \(v+P_{t+1}h-d-I\), with \(I=(1/q-1)d\), the total obligation remains \(d/q\). The old constraint and sign obstruction return.
- **Using a short-period \(q\) with a generational income/housing horizon:** changes the timing units. All asset prices, service costs, income aggregation, and preferences need a common horizon.
- **Restricting repayment to future collateral instead of origination LTV:** a zero-payment balloon with \(d/q\leq\phi P_{t+1}h\) genuinely changes the initial maximum advance to \(q\phi P_{t+1}h\). At constant prices, \(q=.5\), and \(\phi=.8\), this permits only 40% initial borrowing. It cannot be advertised as the same 80% origination mortgage.
- **Allowing costless repeated cash-out borrowing during young age:** changes the prescribed-amortization feasible set. Payments can potentially be refinanced; (4) and (10) describe a fixed origination schedule with no such refinancing. The existing once-young mortgage can be retained, but this restriction should be explicit.
- **Using initial assets separately:** is economically legitimate when an actual earlier payment must occur before subsequent earnings arrive. Once current income has arrived, its amount belongs in the early cash check. Equation (10) gives that chronology directly.

## 6. What the cited papers support

The August 1, 2026 version of Coven, Golder, Gupta, and Ndiaye uses a two-period illustration in Section 3.1 with young labor income, borrowing \(B_t\leq\lambda P_tH\), and repayment \(RB_t\) when old. It therefore has capitalized-interest one-period debt. Housing is an asset in that illustration; the old sell it, and direct housing utility is omitted. Its annual-looking numerical example uses \(R=1.05\), not a generational \(R=2\). It does not establish the present model's young-old service allocation or warm-glow theorem. The richer model uses amortization to evaluate origination payment-to-income limits and then leaves repayment speed unconstrained. [August 2026 paper, Sections 3.1 and 3.3](https://abdouecon.github.io/research/papers/Property_Tax.pdf#page=16).

The January 31, 2025 version's Section 2 instead gives old households utility from bequeathed wealth alone. Its initial down-payment equation is expressly described in footnote 2 as a restriction on the initial state, rather than an additional purchase constraint. Current income appears in that equation. That setup neither preserves the current old consumption/housing/warm-glow problem nor justifies excluding an already received paycheck from closing cash. Its borrowing bound also discounts the future depreciated collateral value. These differences must be translated explicitly rather than imported as an equivalent 80% initial mortgage. [January 2025 paper, Section 2](https://www.bwl.uni-mannheim.de/media/Lehrstuehle/bwl/Area_Finance/Finance_Area_Seminar/FSS_2025/Arpit_Paper.pdf#page=8).

Chambers, Garriga, and Schlagenhauf separately study down-payment requirements, payment schedules, and amortization schedules in a model with liquidity-constrained households. That separation supports treating these as distinct economic features, rather than using an arbitrary terminal balance to obtain a desired sign. Their model is a richer life-cycle exercise; it is not needed to implement the two-checkpoint algebra above. [Mortgage Contracts and Housing Tenure Decisions, primary working paper](https://fraser.stlouisfed.org/files/docs/publications/frbsl_wp/2007-040.pdf).

## 7. Recommendation and remaining decision

The interest-servicing variant is the smallest coherent amendment if the author accepts the proposal's initial availability of the entire young endowment. Describe it as an interest-only mortgage with interest paid during young age, keep the true 80% origination share, and explain that the aggregate cash requirement includes interest. It preserves old income, both size caps, warm-glow utility, and net-asset notation; the old resizing decision remains separate.

If the economic point is specifically that a family has enough earnings to service an ordinary mortgage but insufficient cash to close, the exact two-checkpoint construction is closer to that mechanism. It accommodates full amortization and uses current income in the down payment. Its additional ingredient is the date at which the remaining young earnings arrive, not extra housing choices or a quantitative model.

The next specification decision is therefore concrete: **use the interest-servicing contract with the existing initially available young endowment, or use a level-payment mortgage with current cash distinguished from subsequent working-age earnings?** After that decision, the lead still needs explicit primitive conditions for the equilibrium branch pattern and a separately defined feasible welfare comparison. This review does not certify either.

Validation: checked the original and proposed budgets, verified (1) exactly for the three limiting contracts, checked the level-payment recurrence against (2), and checked every discounted intermediate balance against the endpoint interpolation in (9). An illustrative 20-payment young horizon with a 30-payment loan and \(q=.5\) gives \(\rho=.4530818393\), \(A=.7734590803\), and aggregate cash coefficient \(1-q\rho\phi=.8187672643\) at \(\phi=.8\). This is an accounting example, not a selected calibration or equilibrium. No model solve was run.
