# Structural review (ChatGPT Pro), September 15, 2026

Verbatim copy of ChatGPT Pro's answer to
`docs/prompts/MASTER_PROMPT_model_review_20260915.md`, pasted by the author
on September 15, 2026. ChatGPT worked from the prompt text, the September 14
slides and the ledger text; it did not read the repository. Inline source
chips from the ChatGPT interface (site names such as "users.ssc.wisc.edu",
"Associazione Economica Americana", "Pasted text.txt") are reproduced as
pasted. The consolidated reading of this review and the Fable review is in
`structural_model_review_consolidated_20260915.md`.

---

Structural review

I would retain the core housing–fertility mechanism, but not proceed to a new baseline calibration until the demographic accounting, financial closure, and empirical population are settled. More preference parameters are not the immediate priority. The mechanism in your brief can operate without them, but several current assumptions have stronger implications than their descriptions acknowledge.

My confidence is highest in the analytical and accounting corrections below. Recommendations about additional economic structure are more conditional: neither a good fertility-level fit nor a successful housing-around-birth fit establishes that the model's responses to housing costs are correct.

## P1. Utility, equivalence scales, and the linear child benefit

(a) Stakes. The specification determines what children cost, why households want them, and how these trade-offs change with resources. Those are separate questions. An equivalence scale does not, by itself, provide a theory of fertility preferences.

(b) Defensibility. CRRA over consumption and housing, a family-space requirement, and a reduced-form benefit of parenting are defensible building blocks. But the literature citation needs qualification. Scholz, Seshadri, and Khitatrakun use an equivalence scale with exogenous family structure and an objective containing \(e\,U(c/e)\), not simply \(U(c/e)\). Your specification imports the scale's shape, not their complete objective. When family size is chosen, that distinction affects behavior and welfare. Likewise, the quantity–quality interaction in Becker–Lewis is not a derivation of a linear \(\psi m\) term. users.ssc.wisc.edu

(c) Minimal fix. Retain the specification provisionally, but explicitly describe it as parental utility over equivalized living standards plus the benefit of children at home. Do not multiply by the equivalence scale merely to match a citation: that would be a substantive change. Also acknowledge that linear child benefits supply no diminishing marginal benefit. In the displayed specification, with \(\sigma=2\) and equivalence-scale exponent \(0.7\), flow utility is actually convex in additional children at fixed consumption and housing among households already parenting: the equivalence-scale burden increases less than proportionally. This does not prove corner lifetime fertility, but it places considerable responsibility for stopping on dynamics, other constraints, and shocks. The displayed parameters and functional form establish this implication.

(d) Discriminating observable. Subsequent-birth probabilities by existing number of children, resources, and age; child-related consumption changes; and fertility responses to income changes. A completed-fertility mean cannot distinguish these mechanisms.

(e) Recommendation. Keep, provisionally. Do not introduce concave child benefits or a different cardinal scaling until a specific behavioral failure justifies the change.

## P2. Does the model actually make children increase the marginal value of space?

(a) Stakes. Three different properties must not be conflated: children increasing the marginal utility of housing; children increasing housing demand relative to other consumption; and lower housing costs increasing optimized fertility.

(b) Defensibility. The prompt's concern about outer CRRA is algebraically reversed. For the non-child-benefit component,
\[
u_s=C^{-\sigma}C_s,
\qquad
u_{sm}
=
C^{-\sigma}
\left(C_{sm}-\sigma\frac{C_sC_m}{C}\right).
\]
Consequently, if \(C_m<0\), \(C_s>0\), and \(C_{sm}>0\), CRRA strengthens rather than overturns the positive cross-partial.

However, your specific aggregator does not satisfy that assumed \(C_{sm}>0\) for additional children among existing parents. With
\[
C=\frac{c^{\alpha_0}(s-h_P)^{1-\alpha_0}}{e(m)},
\qquad m>0,
\]
the parenthood requirement is constant, so
\[
C_{sm}=-\frac{e'(m)}{e(m)}C_s<0,
\qquad
u_{sm}=(\sigma-1)\frac{e'(m)}{e(m)}u_s.
\]
Thus \(u_{sm}>0\) at \(\sigma=2\), despite \(C_{sm}<0\). The current presentation source still states the restriction on \(C\), whereas the displayed quantitative utility gives the opposite sign for this margin. See latex/september_14_presentation.tex:127–145, 414–428. The implementation's equivalence-scale multiplier \(e^{\sigma-1}\) is consistent with the displayed utility; I checked it in code/model/intergen_eqscale_seq_optimized/solver.py:2241–2280.

(c) Minimal fix. State the property on \(u\), not \(C\), and use finite differences because children are discrete. More importantly,
\[
\frac{u_s}{u_c}
=
\frac{1-\alpha_0}{\alpha_0}\frac{c}{s-h_P},
\qquad m\geq1.
\]
Additional children do not directly shift the intratemporal housing–consumption trade-off among parents. They increase both marginal utilities through the equivalence scale. The first child has a distinct space effect; subsequent children need not. Larger families can still occupy larger homes through income, saving, and selection, but that is a different mechanism. Retain parenthood-only space needs if that is the intended claim. Add a per-child space component only if the paper intends and empirically supports a separate space requirement for each additional child.

(d) Discriminating observable. Housing versus nonhousing expenditure changes at first and subsequent births, conditional on resources—not merely an unconditional room gap between differently sized families.

(e) Recommendation. Change the shape restriction and economic interpretation. Positive \(u_{sm}\) alone does not establish the sign of the optimized fertility response to prices.

## P3. The first-birth cost \(\xi\)

(a) Stakes. This parameter separates the decision to become a parent from decisions about additional children. It can also conceal a missing cost of parenthood if interpreted too broadly.

(b) Defensibility. Its role differs from taste dispersion. Writing
\[
\Delta V^H=V^H(n+1,m+1)-V^H(n,m),
\]
your attempt equation implies
\[
\log\frac{p^A}{1-p^A}
=
\frac{\pi_a}{\kappa(n)}
\left[\Delta V^H-\xi\mathbf 1\{n=0\}\right].
\]
Holding continuation values fixed, \(\xi\) shifts the deterministic first-birth index; \(\kappa_1\) rescales that index and its response to costs. They are not interchangeable. But their separate identification does not follow from assigning one moment to each parameter. The ledger explicitly says the identification mapping has not been verified.

(c) Minimal fix. Keep \(\xi\) as an explicit fixed utility cost of entering parenthood, not an unspecified goods expenditure. Compare its implications with a zero-\(\xi\) restriction while preserving the same initial completed-fertility target through the existing \(\psi_0\) root. Removing it might preserve the childlessness mean; there is no analytical guarantee that it preserves the joint age–wealth pattern of first births.

(d) Discriminating observable. First-birth hazards across age and liquid wealth, together with completed childlessness. Similar childlessness rates can coexist with very different responsiveness to affordability.

(e) Recommendation. Keep, unless the restricted model reproduces those joint patterns without it. Calling it a childlessness parameter is descriptively fair; that alone is not an objection.

## P4. Children ever born versus children at home

(a) Stakes. Collapsing \(n\) and \(m\) would confuse reproductive history with current household needs.

(b) Defensibility. Keeping both is necessary in the current model. A household with \(m=0\) may never have had children or may be an empty nest. Those households face different first-birth costs, taste-shock scales, remaining fertility opportunities, and completed-fertility statistics. Removing child-dependent bequests does not remove these distinctions. These roles follow directly from the fertility problem in the brief.

(c) Minimal fix. Preserve both states, but store only feasible combinations, such as \(0\leq m\leq n\). Outside fertile ages, some continuation computations can potentially be shared across histories with identical remaining preferences and constraints.

(d) Discriminating observable. Empty-nest rates and completed births by age. A one-state model cannot generally reproduce both without additional history information elsewhere.

(e) Recommendation. Keep both. This is useful state information, not redundant notation.

## P5. Bequests and old households' large homes

(a) Stakes. Bequests influence old-age saving and the housing available to younger households. But a motive to preserve wealth is not automatically a motive to preserve a particular house.

(b) Defensibility. A nonhomothetic warm-glow bequest motive is a reasonable reduced form. De Nardi's intergenerational-link framework provides a rationale for distinguishing intentional bequests from accidental estates and for studying their effects on wealth accumulation. It does not establish that households specifically prefer to bequeath housing rather than equally valuable financial assets. Federal Reserve Bank of Chicago

(c) Minimal fix. Define \(w^e\) as net bequeathable wealth, with a consistent death-date house price, mortgage deduction, and treatment of any actual liquidation costs. Keep the bequest motive independent of \(n\). Multiplying \(B(w^e)\) by the number of children is not innocuous: for \(\sigma>1\), the displayed unnormalized \(B\) is negative, so child-number scaling can introduce a fertility incentive driven by its level normalization. Genuine child-specific altruism requires a deliberate specification of recipients and estate division.

(d) Discriminating observable. Old-age downsizing and portfolio composition conditional on total wealth, alongside estate sizes. Aggregate bequest flows do not identify why wealth remains in housing.

(e) Recommendation. Keep the net-estate warm-glow motive; do not restore mechanical child-number scaling. Attribute old-house retention to the actual housing frictions and preferences that produce it.

## F1. Within-period fertility and housing timing

(a) Stakes. Your timing gives households the option to adjust housing after observing whether a birth occurs. That reduces the cost of committing to children relative to a model requiring housing to be chosen beforehand.

(b) Defensibility. This is a coherent information structure, not merely a solution technique. Simultaneous choice is equivalent only if it permits the same birth-contingent housing plans and information. Committing to a house before birth realization is a genuinely different economic problem. The displayed sequence clearly grants post-birth adjustment.

(c) Minimal fix. Keep the sequence initially, but define what a four-year period represents: a household may plan a birth and adjust housing during that interval. Do not interpret the sequence as evidence that real households wait until delivery to move. A commitment-before-birth variant should be a targeted robustness exercise, not an automatic replacement.

(d) Discriminating observable. Moves, ownership transitions, and room changes before versus after first birth, including households that move before unsuccessful or postponed fertility plans where those are observable.

(e) Recommendation. Keep, with an explicit informational interpretation. Change it only if the pre-birth commitment margin proves consequential.

## F2. Fertility and tenure taste shocks

(a) Stakes. Economic preference heterogeneity and numerical smoothing produce similar-looking choice probabilities but have different welfare meanings.

(b) Defensibility. Fertility shocks can represent unobserved motivations for having a child. Housing shocks can also represent genuine unobserved preferences. But the ledger records \(\kappa_H=0.005\) as a search-bound solution, not an interior estimate, and says its proposed auxiliary validation exercise was not implemented. That does not prove nonidentification, but it does not establish a structural interpretation either.

(c) Minimal fix. Remove \(\kappa_H\) from the estimated parameter vector for now and treat it as a fixed numerical approximation, with a stated limiting target as smoothing vanishes. Check choice-menu dependence: for \(K\) equal-valued alternatives, a log-sum value adds \(\kappa_H\log K\). Refining a housing grid must not create an economically meaningful new taste benefit. Also specify the extreme-value location normalization. Because fertility shock scales differ by reproductive history, omitted scale-dependent constants need not be harmless for earlier fertility choices.

(d) Discriminating observable. Conditional tenure transition probabilities and persistence, plus stability of economic outcomes across smaller smoothing values and alternative housing grids. Prediction accuracy alone is not a complete identification argument.

(e) Recommendation. Change the treatment of housing smoothing. Retain fertility shocks as economic primitives, with explicit utility units and normalization.

## F3. Stochastic child maturation

(a) Stakes. The child clock determines housing needs, the duration of parenting benefits, the spacing of births, and the timing of new adult households.

(b) Defensibility. Deterministic adult aging alongside stochastic departure from dependency is not inherently inconsistent. The problem is interpreting a memoryless departure process as biological maturation. With per-period departure probability \(\mu\),
\[
\Pr(D=4k)=(1-\mu)^{k-1}\mu,
\qquad
E[D]=\frac4\mu.
\]
There is no minimum duration and no maximum duration. The ledger confirms that newborns enter the first maturation draw and identifies very late dependent-child observations as a consequence of this approximation.

(c) Minimal fix. My preferred change is a bounded stochastic dependency duration, represented by a short child-age-stage or birth-history state. Stochasticity can remain; memorylessness need not. This explicitly reopens the recorded decision to set bounded durations aside. If you retain the current process instead, call it departure from dependency and stop treating its exits as literal age-18 maturation without an additional consistency bridge.

(d) Discriminating observable. Dependents by parent age, time since birth, and sibling age configuration—not just the average duration of dependency.

(e) Recommendation. Change the duration structure for a model making long-run demographic claims. Independence across siblings is secondary to the missing age dependence.

## F4. Dependents when a parent household exits

(a) Stakes. This is both a resource-accounting issue and a population-accounting issue. Surviving dependents require a residence and resources even when their original household disappears.

(b) Defensibility. An important source qualification: the reported roughly 1.2% flow comes from an earlier checkpoint, not a fresh measurement of the retained baseline. Moreover, the ledger says dependents disappear from the household distribution while a separate recorded-birth queue can preserve their future adult entry. Thus the established problem is not simply "the population loses these people"; it is a missing link between their ongoing care, housing, and future entry.

(c) Minimal fix. If the current dependency process remains, prefer an aggregate care pool over automatic joint death or household-level reassignment. The pool must include goods and housing needs, funding, and a consistent exit rule. It must not generate a second adult-entry stream alongside the birth queue. Joint death changes the mortality model; it is not just bookkeeping. Reassignment adds household-state and fertility effects. If F3 eliminates dependency before any positive parental-death hazard, establish that invariant instead: a care pool may then be unnecessary on the model's support.

(d) Discriminating observable. Conservation of dependent persons, care-pool stocks and resources, and subsequent adult entry, conditional on parental exit.

(e) Recommendation. Change the accounting. Use a care pool only where the chosen dependency and mortality specification requires it.

## F5. Exogenous fecundity and endogenous attempts

(a) Stakes. This separation determines the distinction between delayed desired births and births that never occur.

(b) Defensibility. Separating an age-dependent success technology from a fertility decision is coherent. But the success probability must represent success conditional on the modeled attempt, not an observed age-specific birth rate that already contains choices. The model's attempt/success distinction is explicit in the brief.

(c) Minimal fix. Define the four-year attempt precisely, and translate biological information to that horizon consistently. Separately examine the restriction of at most one birth per four-year period: it imposes spacing and can limit catch-up after delayed first births. Do not attribute all resulting losses in completed fertility to biology. With utility shocks attached to attempts but no explicit attempt cost, also avoid interpreting estimated attempt probabilities too literally when success probabilities are very low.

(d) Discriminating observable. Intended versus realized fertility, age-specific unsuccessful attempts where available, birth intervals, and completed fertility conditional on age at first birth.

(e) Recommendation. Keep the biological/choice split. Validate the time aggregation and birth-count restriction before making strong timing claims.

## H1. The housing-finance constraint set

(a) Stakes. The mechanism requires a distinction between lifetime affordability and resources available to enter ownership. It does not require every mortgage-market institution.

(b) Defensibility. Down-payment restrictions, collateralized borrowing, uninsurable income risk, and transaction costs are well-established housing-model ingredients. But the current description contains two substantive ambiguities. First, freely rechoosing signed financial assets subject to a house-value borrowing floor allows additional borrowing against an unchanged house; that is not literally "no refinancing." Second, the displayed down-payment inequality does not switch off for a nonmover: its right-hand side becomes zero, leaving a net-equity restriction. The source confirms that housing trades deliberately precede current earnings; the gross-return factors are therefore not, by themselves, an accounting error. kamilasommer.net

(c) Minimal fix. State the down-payment condition as a logical restriction on the relevant purchase transactions, and specify separately the restrictions on incumbent debt and new borrowing. Describe the signed-asset formulation as a reduced-form collateral credit technology unless gross mortgage debt and liquid assets are separately tracked. Keep \(\phi=0.80\). Do not add rate risk, amortization, and HELOCs simultaneously. The more immediate concern is the strength of excluding four years of current earnings from purchase liquidity.

(d) Discriminating observable. Ownership entry conditional on liquid assets and expected earnings; equity withdrawal by incumbents; and the incidence of binding purchase versus ongoing debt restrictions.

(e) Recommendation. Change the formal financial description and purchase-condition logic. Preserve the parsimonious core rather than adding a full mortgage industry.

## H2. Underwater rollover

(a) Stakes. A price decline can either reduce a borrower's equity or trigger forced repayment. Those are very different transmission mechanisms.

(b) Defensibility. Grandfathering existing debt has an economic rationale. A model that reapplies an origination collateral requirement to all outstanding debt can create artificial forced deleveraging. Kaplan, Mitman, and Violante explicitly distinguish origination restrictions from restrictions on existing mortgages. The particular age taper described in your ledger, however, requires its own justification. Giovanni L. (Gianluca) Violante

(c) Minimal fix. Keep the principle that a fall in the collateral price does not automatically call existing debt. Replace—or explicitly rationalize—the age taper with a contract rule distinguishing scheduled outstanding debt from new borrowing. Inherited negative equity is not the same as permission to originate unsecured credit. Any sale shortfall or terminal unpaid debt must have a specified creditor and settlement rule.

(d) Discriminating observable. Debt repayment, moving, and refinancing conditional on negative equity and age, especially around price declines.

(e) Recommendation. Change the rollover specification and expose it in the model. Do not simply delete it from the code to match an incomplete presentation.

## H3. The rental size cap

(a) Stakes. This is one of the assumptions carrying the main mechanism. A hard cap makes ownership necessary for sufficiently large housing consumption by construction.

(b) Defensibility. A restricted rental menu can be a useful approximation. It is not the same claim as "large rental homes do not exist." In the publicly accessible December 2024 version of Greaney, Parkhomenko, and Van Nieuwerburgh, the rental upper bound is calibrated using housing-size distribution information. That supports modeling segmentation; it does not establish a universal six-room cap for your population. Associazione Economica Americana

(c) Minimal fix. Treat the cap as an empirical approximation, not an externally established institutional fact. Discipline the full joint size–tenure distribution, and include a targeted alternative with scarce but nonzero large-rental availability—through size-specific supply or a steep size-dependent rental premium. Keep one pooled market. The purpose is to test whether the result requires literally zero access or survives expensive and limited access.

(d) Discriminating observable. Rental shares and rents across the housing-size distribution, particularly the upper tail relevant to families, with comparable definitions of rooms and quality.

(e) Recommendation. Change the empirical foundation and require a soft-tail robustness comparison. A hard-cap baseline can remain, but not as an unquestioned fact.

## H4. Housing supply and forward-looking rents

(a) Stakes. Instantaneous supply can absorb demand changes that would otherwise move prices and redistribute an inherited durable stock.

(b) Defensibility. The rent equation is forward-looking, not static: it includes the next-period asset price. Under a deterministic continuation and competitive landlords it is coherent. The more serious issue is using an instantly reversible stock-supply schedule throughout a long transition. Also, Baum-Snow and Han's headline estimates distinguish floor-space and housing-unit supply; neither automatically supplies your selected aggregate elasticity or its adjustment horizon. Their reported average urban elasticities are approximately 0.5 for floor space and 0.3 for units, not a direct derivation of 0.63 for this model. RCNi Company Limited

(c) Minimal fix. Introduce an inherited aggregate stock and a construction margin, for example
\[
H_{t+1}=(1-\delta_{\mathrm{net}})H_t+I_t,
\qquad I_t\geq0,
\]
with construction costs or an investment-supply schedule. Define \(\delta_{\mathrm{net}}\) consistently with paid maintenance: do not count the same physical depreciation twice. Retain conditional deterministic user cost,
\[
r_t=(R_b+\delta_H+\tau_t^p)P_t-E_tP_{t+1}.
\]
With successive surprises, use the forecast known at that date, not prices generated by later unanticipated news. Alternatively, if \(H^s\) means occupied services rather than physical stock, explicitly account for the inherited stock and vacancies.

(d) Discriminating observable. Construction, vacancies, rents, and prices following demand changes at four-year and longer horizons.

(e) Recommendation. Change supply dynamics; keep conditional forward-looking rental pricing unless aggregate risk is central to the intended claims.

## H5. Who owns the rental stock?

(a) Stakes. Rental income, property-tax payments, capital gains, and losses must accrue somewhere. Otherwise the resource and welfare accounts are incomplete.

(b) Defensibility. An outside competitive landlord sector is an admissible closure when the financing return is exogenous. Domestic household landlords are another admissible closure, but not an equivalent one: Sommer–Sullivan–Verbrugge explicitly model household rental-property decisions, while Kaplan–Mitman–Violante discuss the distinction between household landlords and a rental company. kamilasommer.net

(c) Minimal fix. Specify an outside rental company, its financing, stock ownership, tax payments, and cash flows. State that resident-household welfare excludes nonresident owners unless their welfare is separately counted. Zero expected excess return does not mean that rent payments or unexpected capital losses disappear. If the model is instead financially closed, allocate landlord claims and returns to domestic owners. Do not silently distribute them through the property-tax rebate.

(d) Discriminating observable. Rental-stock ownership and rental-income incidence. Internally, every rent payment and asset revaluation should have a counterpart.

(e) Recommendation. Change the closure documentation and accounts. A separate optimizing household-landlord state is not automatically necessary.

## H6. The owner housing-service premium

(a) Stakes. \(\chi\) can represent nonfinancial benefits of ownership, but it can also compensate for a misspecified ownership budget or rental menu.

(b) Defensibility. An owner-service premium is a legitimate reduced form; Greaney and coauthors use such a distinction. It should not be described as a missing financial motive when ownership already provides an asset position and collateral services. Associazione Economica Americana

(c) Minimal fix. Keep \(\chi\), but separate physical space from services. If the family requirement is physical, then applying utility to \(\chi(h-\bar h(m))\) is not equivalent to applying it to \(\chi h-\bar h(m)\): the latter lowers the owner's physical minimum to \(\bar h(m)/\chi\). State which interpretation is intended and verify it consistently across code and exposition. The displayed specification uses services inside the requirement.

(d) Discriminating observable. Ownership choices conditional on dwelling characteristics and financial user cost, together with tenure-specific housing adjustment at birth.

(e) Recommendation. Keep, after resolving the physical-space interpretation. Do not use \(\chi\) as an unrestricted ownership-fit residual.

## D1. Household, person, dependent, and estate accounting

(a) Stakes. Without consistent units, fertility can affect future housing demand through a conversion convention rather than an economic mechanism.

(b) Defensibility. A unitary household need not literally contain two adults. Conversely, using a two-adult-normalized equivalence scale does not establish a two-person demographic unit. The ledger confirms that person and household aggregates can satisfy their separate checks while the dependent-child link remains unresolved.

(c) Minimal fix. Use an unnormalized household measure \(G_t\), with household mass \(\mathcal H_t=\int dG_t\). Define total dependents as
\[
D_t=\int m\,dG_t+Q_t,
\]
where \(Q_t\) contains dependents outside ordinary households. With consistent within-period timing and no migration, require
\[
D_{t+1}=D_t+B_t-M_t-\Delta_t^D,
\]
\[
A_{t+1}=A_t+M_t-\Delta_t^A,
\]
where \(B_t\) is births, \(M_t\) is movement into the adult population, and \(\Delta_t^D,\Delta_t^A\) are deaths. Household formation is a separate mapping from adult persons to \(\mathcal H_t\); it is not automatically one household per maturing child. Orphan transfers are internal movements, not person losses. Physical housing clears using \(h+h^R\), plus any care housing—not \(\chi h+h^R\). Estates must reconcile net assets, recipients, and any losses exactly once.

(d) Discriminating observable. Reproduction of person stocks, household stocks, household size, dependent counts, and estate flows under the same units and transition operator.

(e) Recommendation. Change. This is a prerequisite for population and intergenerational-allocation claims, not an optional extension.

## D2. Historical and counterfactual population closure

(a) Stakes. The counterfactual requires additional births to become additional future adults and households through a consistent delayed mechanism.

(b) Defensibility. The uploaded prompt's mixed historical-headship description is not the latest retained-branch description. The September 15 M26 correction identifies a retained branch with endogenous household propagation, a recorded-birth queue, future household entry based on births divided by 2.1, and no age-mass rescaling. It explicitly describes 2.1 as an external replacement normalization—not modeled child mortality. This is a ledger-based branch identification, not my independent verification of the running transition code.

(c) Minimal fix. Choose one demographic interpretation. A literal person model needs survival, adult entry, and household formation to agree. An abstract reproductive-household model can use a conversion normalization, but must not present resulting household counts as literal persons. If observed historical household masses are imposed in another branch, describe that history as conditional on those masses and identify the implied formation/rescaling residual. For the counterfactual, preserve the same inherited 2023 distribution and pre-2023 birth pipeline, then allow policy-induced births to affect entry only at the specified later dates.

(d) Discriminating observable. Cohort-specific adult entry, household formation, age distributions, and the lag between a birth change and a household-count change.

(e) Recommendation. Change the demographic contract. Matching initial completed fertility to 2.1 does not itself establish replacement under the model's survival and formation rules.

## D3. PAYGO pensions and the property-tax rebate

(a) Stakes. The experiment changes households' holding costs, housing wealth, and transfers. Later demographic changes also affect pensions.

(b) Defensibility. Fixed payroll taxes with pensions adjusting to balance the budget are coherent:
\[
\varpi_t N_t^{\mathrm{ret}}
=
\tau^w wL_t.
\]
Equal property-tax rebates are also coherent, but define a particular redistribution experiment. They are not an accounting-neutral implementation detail. The brief specifies an equal-rebate baseline and counterfactual with pensions balanced separately.

(c) Minimal fix. State the tax base and recipients explicitly. If both owner-occupied and rental properties are taxed at the rate entering user cost, full rebate requires
\[
T_t\mathcal H_t
=
\tau_t^pP_t\left(H_t^O+H_t^R\right),
\]
absent other uses of that revenue. Any care-pool funding deducted from these receipts changes the full-rebate experiment and must be identified. Report the separate roles of capitalization, net tax payments, rebates, and pension adjustment rather than inferring incidence from age alone.

(d) Discriminating observable. Model-consistent net transfers and housing expenditure by age, wealth, and tenure; retirement-income replacement rates; and the timing of demographic fiscal effects.

(e) Recommendation. Keep these closures as the stated experiment. Make their distributional and intergenerational consequences explicit.

## D4. What economy is being modeled?

(a) Stakes. Combining different geographic populations can make parameters reconcile differences in samples rather than economic behavior.

(b) Defensibility. A pooled model is perfectly acceptable. Calling it national while using selected-metropolitan housing moments is not justified by pooling alone. The active ledger identifies the housing targets as selected-metro moments and the fertility and PSID moments as national, and leaves the scope decision unresolved.

(c) Minimal fix. I recommend a national pooled benchmark, given the intended national fertility interpretation. Remeasure the housing moments nationally and reconcile prices, supply, earnings, and demographic units with that population. The alternative is a pooled selected-metro economy with correspondingly scoped fertility and demographic evidence. That remains one market, not a spatial model, but raises additional migration and sample-support questions.

(d) Discriminating observable. The same target definitions measured nationally and in the selected metros, with population weights and uncertainty shown.

(e) Recommendation. Change the target population contract before recalibration. Do not solve this by relabeling the existing estimates.

## E1. Permanent types plus persistent earnings risk

(a) Stakes. Permanent earnings differences determine lifetime resources; persistent shocks determine how current circumstances forecast future resources. Both can matter for a young household's liquidity constraint.

(b) Defensibility. Permanent heterogeneity is not inherently contrary to your mechanism. Consider
\[
\log y_{i,a}=\mu_a+\theta_i+x_{i,a},
\qquad
x_{i,a+1}=\rho x_{i,a}+\varepsilon_{i,a+1}.
\]
Under independence and a stationary persistent component,
\[
\operatorname{Cov}(\log y_{i,a}-\mu_a,\,
\log y_{i,a+k}-\mu_{a+k})
=
\sigma_\theta^2+\rho^k\sigma_x^2.
\]
The permanent component captures covariance that does not decay. Double counting arises if the persistent process was already fitted to dispersion containing that component. The ledger says the two components currently come from different empirical sources; it does not establish either compatibility or duplication.

(c) Minimal fix. Reconcile the income definition, population, taxes, age effects, variances, and autocovariances before deleting states. For annual-to-four-year conversion, \(\rho_4=\rho_{\mathrm{ann}}^4\), while innovation variance must accumulate across the intervening innovations. A fair persistent-only alternative must be refitted to the same income evidence; merely deleting permanent variance compares different earnings distributions, not just different persistence structures.

(d) Discriminating observable. Earnings autocovariances at multiple horizons, cohort dispersion, forecastability, and liquid wealth relative to predicted lifetime earnings. A short panel may not sharply distinguish a permanent component from a highly persistent one.

(e) Recommendation. Keep permanent types provisionally. Remove them if a comparably disciplined persistent-only process reproduces the relevant evidence—not because less heterogeneity makes the liquidity story cleaner.

## E2. Earnings uncertainty and fertility

(a) Stakes. At the same current wealth, households facing different future earnings distributions may make different fertility decisions.

(b) Defensibility. Your dynamic problem already permits this: the earnings state changes continuation values and the option value of waiting. No additional "uncertainty penalty" is required. Sommer's lifecycle fertility model provides a direct precedent for earnings risk affecting fertility timing and completed fertility, although its findings do not constitute a universal sign theorem. Federal Reserve

(c) Minimal fix. Retain risk in the continuation problem and make its timing explicit. Compare states or environments with similar expected resources but different uncertainty. Do not assert that every mean-preserving increase in risk lowers fertility: the relevant object is its effect on the difference between having a child and waiting, not on either value separately.

(d) Discriminating observable. Fertility responses to changes in earnings uncertainty conditional on expected income and wealth, including whether postponed births are subsequently recovered.

(e) Recommendation. Keep the existing channel. Do not add a second direct uncertainty disutility without independent evidence.

## X1. Which assumptions must be decided jointly?

(a) Stakes. Several apparent single-parameter changes alter multiple mechanisms. Equivalence-scale normalization affects \(\psi,\xi,\kappa_1,\kappa_C\); dependency duration affects parenting benefits, old-age space needs, and entrants; rental segmentation interacts with credit and \(\chi\); landlord ownership interacts with taxes and welfare.

(b) Defensibility. These interactions are features of the model, not reasons to abandon it. They do rule out interpreting independent "moment-to-parameter" labels as established identification or adding separate mechanism effects as though they had no interactions.

(c) Minimal fix. Decide units and demographic closure first; then preference and credit interpretations; then earnings and geographic scope; then identify the revised parameter vector. Keep the initial \(\psi_0\) root as an initial calibration restriction, but do not reapply it to counterfactual environments to restore fertility to 2.1. That would remove part of the response being studied. Also coordinate F3 and F4: changing the child-duration process may eliminate the particular parental-exit problem rather than requiring two separate additions.

(d) Discriminating observable. Joint responses of fertility, ownership, housing size, liquid wealth, and dependency by age. Later comparisons should separate fixed-parameter mechanism changes from recalibrated-model comparisons.

(e) Recommendation. Change the decision and comparison sequence. The active ledger already correctly warns against attributing weak responses to either calibration or economics without separating these experiments.

## X2. The smallest defensible change set

(a) Stakes. The objective should be a model whose restrictions can be understood and challenged—not the largest model compatible with the topic.

(b) Defensibility. The existing architecture does not require endogenous marriage, bargaining, adjustable-rate mortgages, spatial equilibrium, or a fully stochastic aggregate economy before it can contribute. Those are potential extensions, not repairs to every current weakness.

(c) Minimal fix. My minimum package is: correct the preference interpretation and shock normalization; reconcile dependents, persons, and entrants; formalize purchase credit and incumbent debt; specify landlord, estate, and fiscal accounts; choose one empirical population and a compatible earnings process; and give the housing stock a defensible transition interpretation. Retain the separate child states, the fertility attempt/success structure, the core borrowing mechanism, and the basic warm-glow bequest motive. Whether bounded dependency replaces the present process is the largest substantive household-state decision in this package.

(d) Discriminating observable. Before final re-estimation, require one agreed mapping from every important model object to its empirical counterpart. Subsequently, distinguish fit to levels from validation of responses. Housing adjustment around childbirth does not by itself validate the reverse causal effect of housing costs on childbirth.

(e) Recommendation. Change these specific components and definitions. Treat alternative child-benefit curvature, reversed within-period timing, and richer mortgage instruments as conditional extensions rather than automatic prerequisites.

## X3. Additional structural issues

(a) Stakes. The largest additional issue is welfare with endogenous fertility. A second is whether omitted nonhousing costs of children are being absorbed into housing and fertility preferences.

(b) Defensibility. A unitary household model can isolate housing mechanisms, but its fertility shocks and fixed costs may summarize bargaining, time costs, and earnings consequences that respond differently to resources. Doepke–Kindermann's bargaining framework is a reminder that the allocation of child costs within a couple can affect fertility, not a requirement that you reproduce their model. Couillard's housing–fertility framework likewise demonstrates that multiple housing attributes can matter; it does not validate this model's particular decomposition. Associazione Economica Americana

(c) Minimal fix. Define welfare first for an explicitly identified population—such as households alive in 2023 and their remaining utilities—and report future-cohort outcomes separately. More births are not automatically a welfare gain. Any welfare criterion that sums over different numbers of future households requires explicit weights and a utility-level convention. Also ensure that the mathematical feasible set enforces the stated either-rent-or-own choice, that annual flows and four-year utility/discounting conventions match, and that the one-birth-per-period restriction is acknowledged. Finally, consider an externally disciplined child goods/time-cost schedule as a targeted alternative if the model attributes implausibly much of the resource response to housing.

(d) Discriminating observable. Birth-related earnings and nonhousing expenditure changes; birth intervals; welfare changes by incumbent age and tenure; and sensitivity to explicitly different population-welfare criteria. An annual discount factor at its cap is a diagnostic of model fit or identification, not by itself proof of mathematical invalidity. The boundary estimate is documented in the brief. Pasted text.txt

(e) Recommendation. Change the welfare and measurement contract. Do not compensate for omitted mechanisms by tuning parameters until the housing experiment produces a preferred fertility response.

## Decision table

"Blocks recalibration" means it should be settled before estimating the next defensible baseline. It does not mean that ongoing numerical work must stop. Costs refer to the recommended change or audit, not to every possible extension.

| Label | Recommendation | Cost | Blocks recalibration? |
|---|---|---|---|
| P1 — Utility and child benefits | Keep | Low | Yes: cardinal interpretation must be chosen |
| P2 — Complementarity | Change | Low | No: the analytical correction is available now |
| P3 — First-birth cost | Keep | Low | No |
| P4 — Two child states | Keep | Low | No |
| P5 — Bequest object | Keep | Low | Yes: net-estate accounting must be settled |
| F1 — Within-period timing | Keep | Low | No |
| F2 — Taste shocks and smoothing | Change | Low | Yes |
| F3 — Dependency duration | Change | Medium | Yes |
| F4 — Dependents after parental exit | Change | Medium | Yes, jointly with F3 |
| F5 — Fecundity and attempts | Keep | Medium | Yes: horizon and birth-count mapping |
| H1 — Financial constraints | Change | Medium | Yes |
| H2 — Underwater rollover | Change | Medium | Yes |
| H3 — Rental segmentation | Change | Medium | Yes: empirical support and baseline menu |
| H4 — Supply and user cost | Change | Medium | Yes |
| H5 — Landlord closure | Change | Medium | Yes |
| H6 — Owner premium | Keep | Low | Yes: physical-space interpretation |
| D1 — Accounting units | Change | Medium | Yes |
| D2 — Population closure | Change | High | Yes |
| D3 — Fiscal closure | Keep | Low | Yes: complete tax base and recipients |
| D4 — Geographic scope | Change | Medium | Yes |
| E1 — Permanent earnings types | Keep | Medium | Yes: reconcile income-process evidence |
| E2 — Earnings uncertainty | Keep | Low | No |
| X1 — Interactions and sequencing | Change | Low | Yes |
| X2 — Minimum repair package | Change | Medium | Yes |
| X3 — Welfare and omitted margins | Change | Medium | Partly: welfare comparisons require it |

The central distinction is between preserving the mechanism and preserving every implementation choice. The first is warranted. The second is not: especially where the same model currently uses one account of dependent children, another account of future entrants, and a preference interpretation that is stronger than the displayed functional form.
