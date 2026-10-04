# Credit, housing costs and fertility: decision memo

October 4, 2026. Reference: the author-adopted post-interest, soft-financing chain 13, with the retained wealth target. No preferences, earnings, transfers, calibration targets, production model or manuscript were changed. Results below are diagnostic, pending the full acceptance limitation stated at the end.

**There is a useful unchanged-model result: housing costs affect fertility substantially, while mortgage relief mainly changes tenure.** That gives us a paper direction worth developing. It does not yet validate a mortgage-led fertility claim or establish that correcting the income gradient will restore one. The most informative new finding is that credit can have a positive value and bind within a housing option while providing little additional value to having a child rather than waiting.

**1. Consumption, saving and purchase choices share the same budget.**

Tommaso's clarification is correct. For each housing option, the solver optimizes consumption and ending assets; housing choice uses those optimized conditional values, and fertility anticipates them. The computational order does not commit consumption or saving before the purchase decision. Numerical interpolation approximates off-grid optimized values; birth outcomes and nested taste shocks retain their specified timing.

Let \(b\) be beginning financial assets, \(R\) their gross return, \(S\) sale proceeds, \(y\) period income, \(Q=ph\) the purchase price, \(K\) carrying costs and \(b'\) ending assets. A new purchaser satisfies

\[
Rb+S+y=Q+K+c+b',\qquad b'\ge-\phi Q.
\]

The purchase screen \(Rb+S+y\ge(1-\phi)Q\) therefore establishes only a necessary condition. Consumption and saving still compete for resources. Screen passage is not evidence of slack financing.

One actual age-22 household illustrates the difference. It starts renting with zero assets and four-year after-tax income 2.061. A four-room home costs 3.104; its down payment is .621 and carrying costs .304. Even though it passes the screen, its maximum consumption at the borrowing limit is only **1.137**. Its conditional four-room plan hits that limit both while waiting and after a successful birth.

Credit has a positive shadow value in both plans. But the household chooses four-room ownership with probability **4.763% while waiting versus .0123% after success**; the successful-birth branch almost always rents. Raising financing to 95% raises the optimized lifetime value of waiting by .00778 utility units and of success by .00607. Thus the relative value of having the child falls by .00171. All 48 saved endpoint budgets and saving optimality inequalities in the small illustration packet pass.

This extends beyond one household. Across all 249 supported age-22 renter states, covering **21.5% of baseline first-birth responsiveness**, chosen ownership hits a debt floor in 10.7% of waiting plans and 7.7% of successful-birth plans. Holding future choices and current menus fixed, current borrowing relief improves waiting more on average. The child-relative effect is negative in the low- and middle-income groups and positive in the high group. The permanent 80%→95% policy also lowers the average child-relative value in this slice. This is a local constraint calculation and a separate permanent-policy comparison, not a decomposition of the full policy's causes. Every positive-probability owner option was checked; zero-probability options do not establish what opening a new menu would do. [Exact controls, coverage and source checks](diagnostics/joint_budget/README.md).

The origination census is useful only with this qualification. Nearly all first-birth susceptibility passes a parent-sized purchase screen, but ending debt limits remain relevant. Removing income from that inequality without reoptimizing is a screen comparison, not an alternative model or a reason to undo yesterday's adopted timing. We retain the current joint budget. [Census](census/README.md).

**2. A weak mortgage effect coexists with a live housing-cost channel.**

We ran two small, matched diagnostics against the saved baseline. The first changes only the financed share; the second raises the imposed house price by 10%, which also raises rents and changes collateral and incumbent housing wealth. They are normalized lifecycle cross-sections, not market-clearing policy equilibria or dated transitions.

| Diagnostic | Change in explicit births | Changed policies at baseline household states | Changed exposure to those policies | Ownership, ages 18–29 |
|---|---:|---:|---:|---:|
| Financing 80% → 95% | −0.250% | +0.156% | −0.406% | +10.112 percentage points |
| House price and rent +10%; financing unchanged | −6.434% | −7.417% | +0.983% | −3.818 percentage points |

The accounting identity is \(\Delta B=\sum_s G_0(s)\Delta h(s)+\sum_s\Delta G(s)h_1(s)\), where \(G\) is prebirth household exposure and \(h\) the realized birth hazard. A symmetric allocation gives the same qualitative result. The negative credit exposure term includes tenure, assets and birth histories; it does not prove that debt alone reduces fertility. A positive policy term is a birth count, not a welfare result.

Credit barely changes young households' rooms (−.548%). The price diagnostic reduces them by 8.605%. The common-state first-birth policy response to price is −14.912%; the resulting first-birth flow falls only 5.015% because exposure also changes. The price test establishes housing-price sensitivity, not that the room floor alone causes it. The two shocks differ in economic content and are not welfare rankings. [Full diagnostics, controls and graphs](diagnostics/README.md).

The mechanism to develop is therefore **child-relative value**, rather than ownership or the total value of a loan. A family chooses a birth when its optimized value after success sufficiently exceeds its value while waiting. A financing policy can improve both options, or improve waiting more. Rental housing supplies family space even when ownership is financially constrained. The rental alternative is also central to weak housing responses in [Kaplan, Mitman and Violante](https://gregkaplan.me/s/kaplan_mitman_violante_sep2019.pdf), though that paper does not model fertility. Equal saving and borrowing rates do not make credit constraints valueless.

**3. The income gradient is a real concern; its causal role remains unresolved.**

The pasted CPS/model numbers reproduce. Under equal population thirds and approximate interview-age matching, young children-ever-born means are .596/.448/.226 in CPS versus .124/.508/.973 in the model. The opposed gradients survive family-unit checks. Older CPS fertility is approximately flat, while the model remains strongly increasing. But CPS measures current family money income in 2024; the model has current Markov earnings and is calibrated to 2004/2006 fertility. We have not obtained a matched prebirth-resource validation. [Measurement audit](measurement/AUDIT.md).

The additive child reward and CRRA-2 material utility create a force toward richer households having more children: under homothetic proportional costs, the material utility loss declines with resources while the reward stays fixed. Housing necessities and dynamic constraints break exact homogeneity, so this is not a theorem about the complete model. A proportional earnings deduction alone retains that scaling force; realistic time costs or labor-supply mechanisms cannot be dismissed generally. [Jones, Schoonbroodt and Tertilt](https://www.nber.org/system/files/chapters/c8406/c8406.pdf).

The earlier “benefit world” does **not** isolate the gradient: it adds an unfinanced child-dependent cash floor and substantially lowers the child-utility parameter. Its larger credit response cannot be attributed to correcting sorting alone. Moving births toward constrained households would also be insufficient unless credit improves their successful-birth value relative to waiting. [Recovered original experiments](evidence/RECOVERY.md).

Our next empirical requirement is successful first-birth hazards by **prebirth household resources and liquid wealth**, with matched ages and interview exposure. The current PSID history panel lacks the necessary resource sidecar. A legacy all-zero hazard file drops its birth events and cannot fill that gap. We identified the exact fields and risk-set construction needed, rather than treating that file as evidence. [Bounded feasibility check](measurement/prebirth/README.md).

**Recommendation.** Keep the adopted model while validating this incidence margin. Develop the conditional finding that financing changes tenure much more than fertility, whereas housing costs materially change fertility. If a matched prebirth-resource gradient rejects the model, revise the child-cost/resource block using measured parental costs or benefits, retain current identifying targets and add the relevant joint moment. A financing revision needs independent evidence that the full household credit constraint is misrepresented; a per-child space revision needs measured incremental room demand. Neither a mortgage spread nor a transfer should be selected because it makes credit effects larger. No particular revision is yet identified by these diagnostics.

The baseline also misses its existing young-fertility stock target (.534 versus .810). That limits quantitative policy confidence independently of the gradient. Its verified loss is 13.771; the [complete 14-row target table](evidence/baseline/target_fit.csv) and [31-row parameter/restriction table](evidence/baseline/parameters.csv) accompany this memo. Ten scored moments and ten free parameters do not by themselves establish informative rank or search convergence.

The analytical reallocation theorem holds fertility fixed and assumes bridge finance and tailored transfers. It does not promise a positive aggregate birth response. Its cash-before-income convention and the quantitative interval budget must be labelled distinctly; neither model's result automatically validates the other's welfare claim. Stationary policy prices solve birth renewal, and household population scale clears housing along a supply curve with \(H_0\) fixed. Prices can differ across policies. The tax population arithmetic is recovered, but supplier tax incidence, fiscal acceptance and the transition remain separate unresolved objects; it is not a rescue of the credit mechanism.

**Verification limit.** Both fresh cases saved complete inputs/arrays and 17 standard graphs; reached-source identity, probability, distribution, aggregate-policy and exact birth-flow reconstruction checks pass. Full canonical dated-budget/purchase audits are unexecuted: frozen authentication requires a preexisting deleted archive initializer. The local joint-budget checks above pass but do not replace those full audits. These are inspectable diagnostics, not production policy certification. [Exact blocker and receipts](audit/README.md). Fable 5.1 and Astra max completed three bounded debate rounds; the lead rejected universal credit irrelevance, sorting-only attribution and welfare interpretations of birth accounting.
