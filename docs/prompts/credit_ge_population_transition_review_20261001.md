# Mortgage borrowing, fertility and equilibrium population

## Request to ChatGPT Pro

Please act as a mathematically careful quantitative macroeconomics coauthor. I want to understand a surprising result in a lifecycle housing–fertility model and decide which transition experiment and diagnostics to run next. This document is intended to be self-contained: it describes the current experimental model, the financing interventions, the stationary equilibrium calculation, the results, and the questions. You do not have repository access. Source paths in the final section are an audit trail, not attachments you should pretend to have read.

The motivating mechanism is that down-payment constraints make family-sized housing difficult for young households to acquire, potentially delaying or reducing fertility. We can compute stationary general equilibria. In the completed mortgage-only experiments, relaxing the purchase financing constraint raises ownership but slightly **reduces** the equilibrium population scale. Prices decline. Both endpoints satisfy demographic replacement. We want an explanation of this result, not an argument that it is impossible merely because fertility is at replacement in both steady states.

Please distinguish a valid stationary comparative static, a causal decomposition of household behavior, and the existence or stability of a transition connecting the endpoints. Do not infer a transition from two stationary solutions. Do not propose changing preferences, targets, or equilibrium closure simply to produce a desired sign. If the present mechanism does not support the intended positive effect, say so and identify the evidence that would discriminate between economic failure, closure dependence, and numerical error.

The central questions are:

1. Why can easier mortgage borrowing lower stationary equilibrium prices and population while raising ownership?
2. What mathematical conditions determine the sign, and which derivatives or household-level objects would we need to measure?
3. Could births rise initially but population eventually settle at a lower level? What would have to happen along a genuine equilibrium transition?
4. Is the population scale economically pinned under the stated closure, or could its interpretation be compromised by the treatment of entry, estates, fiscal transfers, or housing supply?
5. What is the smallest decisive transition and mechanism exercise to do next?

## 1. Reference identity and boundaries of the evidence

The experiments below use the **same frozen experimental calibration**, selected floor-search chain 7, case 0173. Its reported original-weight loss is 31.28400725566496, its owner asset price is \(q_0=0.719168368828958\), and its population scale is \(N_0=0.9222615666500361\). This is a numerical reference, not a production adoption or a claim of optimizer convergence. The complete target and parameter tables appear below.

This reference differs from a later candidate with loss 30.99909262396404 and from a subsequently launched calibration that normalizes baseline population to one and internally calibrates the supply level. **No results from those later calibrations are substituted into the comparisons in this document.** All policy comparisons here preserve the reference's preferences, earnings, entrant law, target definitions and weights, and supply parameters unless explicitly stated.

There are three distinct credit experiments:

- **Purchase-only mortgage relaxation:** raise the financed share at origination from 0.80 to 0.90 or 1.00. Keep the existing-owner collateral rule at 0.80 and renter unsecured borrowing at zero.
- **Purchase and existing-owner relaxation:** raise both financed shares together. Keep renter unsecured borrowing at zero.
- **Broad credit diagnostic:** replace renter and owner borrowing restrictions with finite-grid continuation-solvency support and remove effective purchase gates. This is a bundled experiment and is not an LTV-only intervention.

A fourth, subsequently completed diagnostic changes **only renter unsecured borrowing** at the fixed reference price. It has no new GE root. It must not be used as if it were the mortgage-only GE result.

## 2. Economic environment and household states

The implemented experiment is one pooled housing market. It is not a spatial equilibrium with migration across cities, although implementation arrays retain a singleton location dimension. Time is discrete; a period represents four years. Households enter adult life at age 18, age deterministically between decision dates, and face age-specific survival. Fertility opportunities are age dependent. Housing size is measured in physical rooms.

Write the inherited household state as

\[
x=(a,b,h^-,z,n,m).
\]

Here \(a\) is adult age; \(b\) is signed net financial assets; \(h^-\) is inherited owner housing, with zero denoting a renter; \(z\) is the current income state; \(n\) is children ever born; and \(m\) is children currently at home. Net financial assets can be negative for owners because mortgage borrowing is represented through this same asset position. The model does not separately track gross mortgages, deposits, and home-equity credit lines. An owner with negative \(b\) cannot therefore automatically be classified as using unsecured credit.

Children ever born and children at home are separate states. The first-birth utility cost and the remaining fertility opportunities depend on reproductive history, while current utility and housing needs depend on children at home. The upper children-ever-born category is a top-coded three-or-more group. The demographic renewal calculation contains an explicit top-code accounting adjustment; raw entries into the top category and adjusted births must not be silently interchanged.

The numerical comparison uses 120 wealth nodes and nine income states. Its wealth support extends from −12 to 3000 in model units. The native distribution and household budget checks are retained. This resolution is an experimental computational grid, not a grid-convergence certificate.

Renters choose housing up to a six-room cap. The owner menu is \(\{2,4,6,8,10\}\) rooms. Owners may keep, sell, or change their product. Consequently the rental and owner choice sets differ. Large housing requirements may induce tenure changes, but this does not establish that a marginal relaxation of purchase liquidity disproportionately benefits parents.

The 17 adult decision ages are 18, 22, …, 82. Fertile decision ages are 18 through 42; retirement begins at 66. Survival to the following age is one through age 62, followed by approximately 0.9391263, 0.9184976, 0.8849522 and 0.8300468 at ages 66, 70, 74 and 78. Age 82 is terminal. Each resident child independently leaves dependency with probability \(2/9\) per period in this experiment. The implied mean dependency duration is 18 years, but this is not deterministic departure at child age 18 and does not track each child's actual age. Adult-entry timing is handled by a separate birth-vintage queue described below.

### Earnings and fiscal objects

Working income is an age profile multiplied by a nine-state Markov income process. The income states are approximately \(0.103468,0.171312,0.283642,0.469626,0.777561,1.287408,2.131562,3.529229,5.843348\); their entry weights are \((1,8,28,56,70,56,28,8,1)/256\). Transitions follow the pinned nine-by-nine matrix. This document does not attach the full matrix, so do not infer its entries from the stationary weights or substitute the old five-state-plus-permanent-types process.

The after-tax four-year working-income profile before multiplication by the income state is approximately 2.65083 at ages 18–22, 3.46647 at 26–30, 4.07820 at 34–42, 4.01703 at 46–54, and 3.81312 at 58–62. Retirement income is the PAYGO pension rather than that working-income process. The baseline payroll tax is 0.08028070961950022 and the reported four-year pension is 0.917784047463731. The stationary PAYGO calculation balances pension payments against payroll receipts on the solved distribution. On a transition it must use the actual dated age/income distribution.

The gross four-year bond return is \(R=1.02^4=1.08243216\). Annual depreciation is 0.01416143718381309, mapped into four-year depreciation 0.05545379079326218. Annual property tax is 0.010598360773872594, with the experiment's four-year convention 0.042393443095490375. The stationary user-cost coefficient is their sum with \(R-1\), approximately 0.1802793939. There is no active means-tested transfer floor, child earnings penalty, birth purchase grant, parental down-payment waiver, payment-to-income screen, rental wedge or owner-size cost in this packet. These absent mechanisms must not be supplied implicitly to explain the result.

### Within-period timing

The household starts with inherited financial assets, tenure/product, income and family states. During fertile ages it chooses whether to attempt a birth. Birth success is realized before current housing and saving choices. The household then chooses tenure/housing, consumption and next-period net financial assets conditional on the realized family state. The child-related housing requirement therefore activates in the birth period itself. Continuation values account for income risk, survival, and children leaving home.

This timing allows birth-contingent housing adjustments. It is not a model in which the household must irreversibly buy a large house before learning whether the attempted birth succeeds. A change to that timing would be an additional economic intervention.

## 3. Preferences and the fertility decision

The active experiment uses a physical housing floor when children are at home:

\[
\bar h(m)=h_P\mathbf 1\{m>0\},\qquad h_P=2.3.
\]

There is no additional per-child room floor. The floor is removed when no children remain at home. It is not multiplied by the household equivalence scale. For renters, discretionary housing services are \(h-\bar h(m)\); for owners the owner premium is applied after subtracting the physical floor, giving \(\chi[h-\bar h(m)]\).

Let \(\omega\) equal one for renters and \(\chi\) for owners. Material living standards are

\[
C(c,h,m,\omega)
=\frac{c^{\alpha}\{\omega[h-\bar h(m)]\}^{1-\alpha}}{e(m)},
\qquad
e(m)=\left(\frac{2+0.7m}{2}\right)^{0.7}.
\]

On the feasible consumption/housing domain, flow utility is

\[
u(c,h,m,\omega)
=\frac{C(c,h,m,\omega)^{1-\sigma}}{1-\sigma}
+\psi_c m^{1-\nu}.
\]

The child term is zero at \(m=0\). In this reference, \(\sigma=2\), \(\alpha=0.733\), \(\chi=1.093187485240681\), \(\psi_c=0.17156192800028292\), and \(\nu=0.09689406390932934\). The direct child benefit is therefore mildly concave. The coefficient reported as `psi_child` multiplies \(m^{1-\nu}\); it is not the coefficient of \(m^{1-\nu}/(1-\nu)\). The parameter table includes the corresponding alternative coefficient separately.

The consumption share is constant across family states. The historical child-dependent expenditure-share and compensating-utility-factor specification is **not active** here. In particular, do not import that earlier specification's \(A(m)\) term into this experiment.

The household discounts four-year continuation utility using \(\beta=\beta_{\rm annual}^4\), with \(\beta_{\rm annual}=0.9670100471250203\). The warm-glow function is normalized to zero at zero estate:

\[
B(W)=\theta_0\frac{(\theta_1+\max\{W,0\})^{1-\sigma}-\theta_1^{1-\sigma}}{1-\sigma},
\quad \theta_0=0.07862721952114816,\quad\theta_1=0.008193084126995582.
\]

**A valuation distinction requiring review:** the utility calculation in this experiment uses \(W=b'+qh\), valuing the house before selling costs. Death-solvency and the estate-funding ledger use net sale proceeds. Do not silently replace gross housing in donor utility with net housing and call it the same model. The bequest child-loading parameter is zero, so this motive is not multiplied by the number of children. There is no estate tax in this packet.

Let \(V^H_a(x;n,m)\) denote the optimized value after current fertility is known and housing/saving choices remain. During fertile ages, write \(\pi_a\) for the exogenous probability that an attempted birth succeeds. The two deterministic choice values are

\[
W^0_a=V^H_a(n,m),
\]
\[
W^A_a=(1-\pi_a)V^H_a(n,m)
+\pi_a\left[V^H_a(n+1,m+1)-\xi\mathbf1\{n=0\}\right].
\]

The first-birth cost \(\xi=0.4009349064519595\) is a utility cost, not a cash down payment. With extreme-value fertility shocks, the attempt probability has the logit representation

\[
p^A_a(x)=\Lambda\!\left(
\frac{\pi_a[\Delta V^H_a(x)-\xi\mathbf1\{n=0\}]}{\kappa(n)}
\right),\qquad f_a(x)=\pi_a p^A_a(x),
\]

where \(\Delta V^H_a=V^H_a(n+1,m+1)-V^H_a(n,m)\), \(\Lambda\) is the logistic function, and \(\kappa(n)\) differs for first and subsequent births. The two scales are 0.13056430591307258 and 0.364420274899206. Tenure/housing choice also has a taste-shock scale, 0.012988772333020718. No new fertility taste parameters are fitted during the credit comparisons.

At a fixed inherited state and price, with the fecundity schedule and taste scales unchanged, the sign of the fertility response depends on how the policy changes the **difference** between the optimized child and waiting values. A higher value of homeownership or a positive utility complementarity between children and housing does not by itself sign that difference.

### A mathematical distinction we want reviewed

For an existing parent, treating \(m\) as locally continuous only to inspect the equivalence-scale channel,

\[
\frac{u_h}{u_c}=\frac{1-\alpha}{\alpha}\frac{c}{h-h_P},
\qquad
u_{hm}=(\sigma-1)\frac{e'(m)}{e(m)}u_h.
\]

Thus a positive cross-partial at \(\sigma=2\) does not imply that additional children alter the intratemporal housing–consumption ratio among parents. The first child introduces a discrete floor shift and requires a finite-difference analysis. Please check these statements and explain which object actually determines the mortgage-policy fertility response after reoptimization.

## 4. Financing, budgets and estate solvency

Use \(q\) for the price per physical owner room, \(r\) for rent per room, \(R\) for the four-year gross financial return, and \(s\) for the selling-cost fraction. The stored variable `psi` in some code is the selling cost; it must not be confused with the child-benefit parameter \(\psi_c\). Here \(s=0.06\).

Separate the purchase financed share \(\phi_B\) from the existing-owner collateral share \(\phi_S\). Both equal 0.80 in the reference. Raising a financed share lowers the required down payment: 80% financing means a 20% down payment, not an 80% down payment.

For a purchase, define post-transaction financial wealth before current cash income as

\[
x=b+(1-s)qh^- - qh,
\]

with sale proceeds included only when a house is actually sold; a same-house stayer does not sell and repurchase each period. The current native purchase rule allows current income to contribute to the liquidity test:

\[
x+\frac{y}{R}\ge-\phi_B qh.
\]

For a renter buying a house, this can be written \(b+y/R\ge(1-\phi_B)qh\). This is materially different from the older presentation's restriction that excluded current earnings. The experiment uses the current rule.

The buyer's subsequent asset choice is subject to the purchase mortgage floor, and to estate solvency when death is possible. For a surviving same-house owner, the native incumbent rule allows inherited debt to be carried forward without requiring it to be reduced immediately to the current collateral limit:

\[
b'\ge\min\{b,-\phi_S qh\}.
\]

At a date with positive death probability, net liquidated estate wealth must additionally be nonnegative:

\[
b'+(1-s)qh\ge0.
\]

The operative same-house floor is then the maximum of the incumbent floor and the estate floor. A buyer likewise cannot use a 100% purchase financing rule to leave an unfunded negative estate after selling costs. The 100% experiments include a reviewed correction enforcing that buyer estate restriction; the earlier failed 100% attempt is not treated as a valid result.

Conditional on a transaction-adjusted financial state \(\widetilde b\), the period budgets are

\[
c+rh+b'=R\widetilde b+y\quad\text{for renters},
\]
\[
c+(\delta+\tau_H)qh+b'=R\widetilde b+y\quad\text{for owners}.
\]

For a same-house owner, \(\widetilde b=b\); a purchaser subtracts the purchase price and adds any net proceeds from selling the inherited house. A selling owner who becomes a renter adds the net proceeds without subtracting a new purchase. In the fixed-zero-unsecured-credit regime, the owner-to-renter transaction also requires nonnegative net financial wealth after sale. This transaction restriction should not be conflated with the next-period renter saving floor. Current-income timing in the purchase gate and the gross-return treatment in the subsequent budget are part of the implemented convention; please check their economic interpretation rather than replacing them with a textbook mortgage budget.

Schematically, the post-birth housing value optimizes current utility plus \(\beta\) times survival-weighted expected continuation and death-weighted warm-glow utility, over these budgets and financing restrictions. Expected continuation integrates the income matrix and children-at-home departure kernel; a same-house owner uses the incumbent saving branch, while a purchaser uses the buyer branch. Housing taste shocks smooth discrete product/tenure choices. A fertility decision compares these optimized housing values after accounting for birth success, not merely the instantaneous utility of a larger house.

Renters in the reference and all mortgage-only interventions satisfy \(b'\ge0\). This is held fixed when \(\phi_B\) or \(\phi_S\) changes. The signed-asset specification is a reduced-form collateral credit technology, not a separately modeled mortgage industry with amortization schedules, mortgage spreads and distinct HELOC contracts.

The broad-credit diagnostic instead removes the effective purchase gates and uses a continuation-feasibility lower bound for renter and owner assets. The lower bound is the first point on the existing grid at which every positive-probability continuation branch remains feasible. At positive death risk a renter has no housing collateral and must have a nonnegative liquid estate. Feasibility on this finite grid is not a proof of the true natural debt limit. The isolated renter-only diagnostic changes just the renter floor, preserving \(\phi_B=\phi_S=0.80\).

## 5. Distribution, entry and stationary equilibrium

Distinguish a normalized household distribution \(\mu\), whose total mass is one, from the aggregate household scale \(N\). The unnormalized distribution is \(G=N\mu\). Do not interpret \(N\) automatically as a headcount of persons; household size is endogenous and a person-population measure requires its own accounting.

The household law of motion propagates current fertility, tenure/product choices, saving, income transitions, adult survival and children-at-home transitions. Adult entrant types and wealth are taken from the frozen entry law in these comparisons. The current entry construction is nonnegative and mean preserving relative to the original five wealth-to-income bin means: negative bin means are set to zero and positive bin means are scaled down to preserve the original weighted mean. This is an explicit experimental reference assumption, not raw individual survey wealth or a policy-dependent inheritance distribution.

The untransformed wealth-to-annual-income bin means are approximately −2.22253, −0.05264, 0.104245, 0.351573 and 3.10345, with approximately one-fifth weight each. Entrants start as nonowners; the wealth law is projected onto the fixed grid preserving its stated mean. The reference mean entrant financial wealth is about 0.186519679. The type distribution is fixed, but its positive financial endowments are not described as free resources: the stationary ledger funds them from available positive net estates and assigns any excess to a residual sink. There is no direct transfer to existing adults and no own-parent-to-own-child inheritance mapping.

This estate-funding ledger is explicitly provisional. It does not certify creditor counterparties or physical housing settlement, and it has no completed transition implementation. In the 100%/100% GE receipt, available estates per normalized four-year cohort are about 0.10854755, positive entrant funding 0.01151450, and the residual sink 0.09703305. Negative-estate liabilities remain recorded, with their creditor rule unresolved. In that particular receipt the net negative estate amount is approximately numerical zero, but that does not establish a complete financial-market closure. Please evaluate whether any conclusion requires more than the ledger proves.

Write \(\mu(q,\phi_B,\phi_S)\) for the stationary normalized household distribution obtained under the resulting policies and the frozen entry law. At a trial price, the model computes adjusted births \(\mathcal B\), the adult entry rate \(\mathcal E\), and mean housing demand \(\bar h\). The stationary renewal root is

\[
F(q,\phi_B,\phi_S)
=\frac{\mathcal B(q,\phi_B,\phi_S)}{2.1\,\mathcal E(q,\phi_B,\phi_S)}-1=0.
\]

The 2.1 conversion is part of the demographic renewal contract, including the model's household/entry convention and birth accounting. It is not evidence that credit is irrelevant. Completed fertility near replacement at both stationary roots is expected under this closure; the equilibrium population **level** can nevertheless differ. The child-utility coefficient is held fixed during this root search: the routine adjusts price, not child preferences.

The top-code adjustment has the form \(\mathcal B=B_{\rm raw}+(w_{3+}-3)B_3\), with \(w_{3+}\simeq3.6023594\) and \(B_3\) the top-bin entry flow used by the native accounting. Dated adult entry splits each adjusted birth vintage equally between 16- and 20-year delays, then converts children to entrant households by division by 2.1. In a four-year time grid,

\[
E_t=\frac{\mathcal B_{t-4}+\mathcal B_{t-5}}{2\times2.1}.
\]

There is no outside entry in the reported stationary closure and retention is one. The selected roots pass a scaled forward-distribution check with constant birth prehistory and these queues. That is a stationarity check, not a test that a policy transition from a different history reaches the selected root. Resident-child departure and adult-entry timing remain distinct model objects.

Stationary rent is linked to price through the fixed user-cost coefficient \(u\), \(r=uq\). The supply rule is

\[
H^s(q)=H_0\left(\frac{uq}{\bar r}\right)^{\eta},\qquad\eta=0.63.
\]

Housing is measured in physical rooms. Given the root price and normalized demand, population scale clears housing:

\[
N\bar h(q,\phi_B,\phi_S)=H^s(q),\qquad
N=\frac{H^s(q)}{\bar h(q,\phi_B,\phi_S)}.
\]

Supply has no separate construction-adjustment lag in this stationary equation. Rental and owner demand are included in physical housing demand; the owner utility premium does not multiply physical room quantities. Payroll/pension and estate/budget checks are separate from housing and renewal clearing. A small normalized-cohort housing residual at \(N=1\) is not the same thing as the final residual after population scaling.

### Calibration normalization versus a policy experiment

A later calibration normalizes reference population to one and profiles the supply level from reference demand:

\[
H_0=\frac{\bar h(q_0)}{(uq_0/\bar r)^{\eta}}.
\]

That is a baseline normalization/calibration step. In a policy comparison, \(H_0\) must then be held fixed. Reprofiling it after every credit intervention would absorb the population response and change the question. None of the numerical results below is claimed to come from that later normalized calibration.

## 6. Completed mortgage-only stationary GE results

All rows below retain the same renter constraint \(b'\ge0\), preferences, entry law, supply parameters and target system. The baseline is \(\phi_B=\phi_S=0.80\). Prices and population are solved again for each intervention.

| Purchase / existing-owner financed share | Owner price | Population scale | Ownership, all ages | Rooms per household |
|---|---:|---:|---:|---:|
| 0.80 / 0.80, baseline | 0.71916837 | 0.92226157 | 0.6593736 | 5.9771254 |
| 0.90 / 0.80 | 0.71770027 | 0.91912454 | 0.6678165 | 5.9898095 |
| 1.00 / 0.80 | 0.71617704 | 0.91751919 | 0.6804383 | 5.9922636 |
| 0.90 / 0.90 | 0.71751282 | 0.91138002 | 0.6896946 | 6.0397144 |
| 1.00 / 1.00 | 0.71226701 | 0.90601584 | 0.7326132 | 6.0474517 |

Relative to baseline, the purchase-only 90% and 100% changes reduce population scale by approximately 0.34% and 0.51%; changing both shares to 100% reduces it by about 1.76%. These are stationary GE comparative statics. They are not fixed-price effects and are not a transition simulation. Price falls, mean rooms increase, and the supply/demand identity gives a lower population scale. That identity checks the result but does not explain why the renewal-compatible price falls.

Selected roots passed household, purchase, estate and fiscal gates. Renewal residuals are within \(10^{-6}\), and housing clears after scaling. The same-point repeat matches all 14 target rows and 31 reported estimates within \(2\times10^{-12}\); each selected packet contains the standard 17 diagnostic plots. Do not silently upgrade this mortgage repeat statement to exact image-hash equality for every mortgage arm. The exact image-hash verification reported later applies to the separately identified broad-credit and renter-only packets.

## 7. Fixed-price diagnostics and what they do not establish

At the reference price, applying mortgage policies to the same inherited baseline distribution gives the following births per normalized household per model period:

| Purchase / existing-owner financed share | Births per household at common inherited states |
|---|---:|
| 0.80 / 0.80 | 0.115253846 |
| 0.90 / 0.80 | 0.115237704 |
| 1.00 / 0.80 | 0.115221425 |
| 0.90 / 0.90 | 0.115266648 |
| 1.00 / 1.00 | 0.115181114 |

The effects are small and are not uniformly positive. Ownership rises, but that does not imply a larger relative value of the child branch. A prior branch-value diagnostic found that purchase relaxation can benefit the waiting/no-birth option more than the birth-success option. This is a hypothesis about the relevant optimized value difference, not a universal theorem that credit lowers fertility.

In the **broad bundled credit** experiment, the reference-price impact on the same inherited states raises births from 0.115253845622 to 0.118598501878. Under the treatment's own stationary normalized distribution at that same price, births instead equal 0.109746607782 and completed fertility equals approximately 2.002188. Its separately solved GE has \(q=0.6550572019731878\), \(N=0.8199246807842695\), ownership 0.748113 and mean rooms 6.339064. The population change is about −11.10%. This bundled GE passed a native repeat, full 14/31 tables and all 17 plot hashes, but only finite-grid occupied support is checked.

In the **renter-only** experiment, purchase and incumbent-owner financing both stay at 0.80 and price stays at \(q_0\). The sole economic change replaces the renter's zero unsecured floor with finite-grid continuation support:

| Object | Control | Renter-only treatment |
|---|---:|---:|
| Births per household, same inherited states | 0.115253845622 | 0.118830411249 |
| Births per household, own stationary normalized distribution | 0.115253845622 | 0.110663632459 |
| Completed fertility, own stationary distribution | 2.099999968 | 2.019078615 |
| Ownership, own stationary distribution | 0.659373600 | 0.628540293 |

The common-state birth effect is +3.103%; the stationary-distribution fixed-price effect is −3.983%. The treatment repeat has identical 14/31 tables and all 17 plot hashes. There is **no renter-only GE result** in this packet. It shows that renter credit alone can generate the fixed-price reversal in this reference; it does not identify the mortgage-only GE mechanism or the transition.

### Distribution accounting

Let \(f_k(x)\) be births under policy regime \(k\), and \(\mu_k\) its distribution. The three-corner accounting is

\[
\int f_1\,d\mu_1-\int f_0\,d\mu_0
=\underbrace{\int(f_1-f_0)\,d\mu_0}_{\text{policy effect at inherited states}}
+\underbrace{\int f_1\,d(\mu_1-\mu_0)}_{\text{distribution component under policy 1}}.
\]

This is an accounting identity with an explicit ordering. The second component is not automatically a causal saving or debt channel. Fertility history, assets and tenure are jointly endogenous, and the reverse ordering need not assign the same contributions.

For the broad-credit comparison, a saved-policy decomposition grouped households by age, income, inherited tenure/product, children ever born and children at home. Of the total −0.008851894096 distribution contribution to births, −0.008259773500 is associated with conditional asset-distribution changes within these groups and −0.000592120596 with group-mass changes. This is about 93.3% in the within-group asset term. Ages 22–30 account for about 90.7% of the first-birth distribution term. Among never-parent renters at age 22, mean inherited financial assets change from about 0.3425 to −0.0067 and the negative-asset share becomes about 54.16%.

These facts locate the distributional shift; they do not show that all of it is caused by deliberate early borrowing rather than changed fertility and tenure selection. That decomposition used the original broad-credit PRE distribution; the later mortgage/renter-only local packets use a separately reconstructed PRE array. Their content hashes differ, even though reported baseline aggregates agree numerically. No elementwise identity across those two PRE arrays is claimed. Each within-experiment comparison holds its own PRE array fixed.

## 8. Mathematical mechanism to derive and test

The following is a proposed local analysis, not an already estimated comparative-static derivative. Let \(\phi\) denote precisely one financing parameter, with the other rules held fixed. If all normalized fiscal and distribution objects are solved internally, and the reduced renewal equation is smooth with \(F_q\ne0\), then

\[
\frac{dq^*}{d\phi}=-\frac{F_\phi}{F_q}.
\]

If fertility/renewal falls with price, \(F_q<0\), then the sign of the GE price response is the sign of the fixed-price **stationary** renewal effect \(F_\phi\). A negative same-price stationary effect would require a lower equilibrium price to restore renewal. The common-inherited-state impact is not \(F_\phi\), because \(F_\phi\) includes the stationary distribution's response.

With the fixed supply function,

\[
\frac{d\log N^*}{d\phi}
=\eta\frac{d\log q^*}{d\phi}
-\frac{d\log\bar h(q^*(\phi),\phi)}{d\phi},
\]

or, separating housing's partial derivatives,

\[
\frac{d\log N^*}{d\phi}
=\left(\frac{\eta}{q^*}-\frac{\bar h_q}{\bar h}\right)
\frac{dq^*}{d\phi}-\frac{\bar h_\phi}{\bar h}.
\]

For two computed endpoints, the finite log identity is exact when the supply primitives remain fixed:

\[
\log\frac{N_1}{N_0}
=\eta\log\frac{q_1}{q_0}
-\log\frac{\bar h_1}{\bar h_0}.
\]

For the 100%/100% endpoint, the supply-price term is approximately −0.006075 and the per-household housing-demand term is approximately −0.011697. They sum to about −0.017772 log points, corresponding to the 1.76% lower level. This decomposition is descriptive; the two terms are jointly determined equilibrium outcomes, not independent causal interventions.

Please verify when this reduced system is legitimate. If population scale affects normalized household incentives, entry types, rebates or pensions, or if a fiscal variable is jointly endogenous, write the full implicit system and identify the additional derivatives instead of silently assuming triangularity. Discrete housing choices and grid kinks may limit differentiability even with taste shocks; propose finite-difference checks where appropriate.

We need more than the arithmetic identity. Explain what makes the stationary renewal response negative, how liquidity versus lifetime resources and family timing enter the optimized fertility value gap, and what identifies those channels. We have not run tighter-than-baseline mortgage constraints, so the finite relaxations above do not establish a global theorem that tightening constraints raises population.

## 9. The transition we actually want

The proposed policy experiment is a permanent, unexpected increase in the **purchase** financed share from 0.80 to 0.90, holding the incumbent-owner rule at 0.80 and renter unsecured borrowing at zero. Start from the complete baseline stationary distribution and inherited demographic queues. Preserve the calibrated supply function, preferences, survival/fecundity schedules, earnings, entry-type distribution and fiscal rules. Existing households retain their actual assets, house product and family history at announcement. Do not initialize them at the new stationary distribution.

This is a proposed next experiment, not a completed path. It is intentionally a modest mortgage-only change before considering 100% financing or changes to incumbent collateral rules. The document does not authorize extra economic assumptions to obtain a path.

For a dated equilibrium we need household policies based on the anticipated price/rent/fiscal path, a forward distribution law with endogenous births and delayed adult entry, housing clearing in levels, and the appropriate fiscal and estate accounting every date. The full state includes inherited cohorts and the birth-to-adult-entry pipeline, not only a normalized household histogram. Initial population scale is fixed by the inherited baseline; \(N_t\) should subsequently follow demographic accounting rather than being freely reset to \(H^s_t/\bar h_t\) each date. The latter is an equilibrium condition the price path must satisfy, not permission to overwrite the population law.

Please resolve the asset-price/rent issue explicitly. A stationary rent map \(r=uq\) is valid for the stationary experiments; it must not automatically be used on a transition if the landlord no-arbitrage condition includes expected capital gains. Similarly distinguish instantaneous flow supply, an installed housing stock with depreciation, and construction adjustment. A stationary supply curve alone does not uniquely specify construction dynamics. Identify which transition conclusions depend on those additional choices.

Report births per household, total births, cohort completed fertility, adult household population, person population where defined, entry flows, ownership and purchase rates, rooms per household, total housing, house prices, rents, net assets/debt, binding credit constraints and welfare separately. A fall in a rate, a fall in total births, and a fall in the eventual population level are different claims.

We particularly want to know whether credit could raise fertility for the initially constrained cohorts but lower fertility later through wealth/family composition, and whether such a path can lead to the lower stationary household scale we computed. Do not assume it does: assess stability, multiple stationary roots, the role of delayed births/entry, and terminal-horizon errors. An endpoint with replacement fertility does not itself establish dynamic reachability.

## 10. Requested answer and practical research plan

Please provide:

1. **A verdict in plain economics.** Is the observed stationary response internally coherent? What is established and what is only conjectured? State whether the result threatens the intended mortgage–fertility mechanism, and why.
2. **A mathematical derivation.** Write the household fertility-value comparison and the stationary equilibrium Jacobian; derive transparent sign conditions for prices and population. Check the log decomposition and its assumptions. Separate direct policy, endogenous distribution, price and housing-demand channels.
3. **A minimal illustrative model.** If useful, give a two- or three-age example in which easing borrowing helps initially but lowers later fertility. Explicitly label it illustrative and show which assumptions are necessary. Do not present it as a proof of our numerical mechanism.
4. **A transition specification.** List the inherited state, equilibrium unknowns, forward and backward equations, demographic lags, fiscal and asset-price/rent requirements, and terminal conditions. Identify missing closure choices rather than inventing them.
5. **A short ranked diagnostic sequence.** Prefer measurements from saved policies and distributions before new solves. Then specify the smallest finite differences and one mortgage-only transition that would discriminate among the leading explanations. For each, give the measured object, held-fixed objects and what outcomes would reject the hypothesis.
6. **Numerical and economic falsification checks.** Include debt and estate feasibility, unchanged entry law, no hidden H0 reprofiling, choice-set/solver comparison, full mass conservation, birth/entry accounting, multiple roots, horizon extension and exact repeat. Do not treat solver success as economic validation.
7. **Interpretation for the paper.** Explain what we could honestly claim if the lower-population result survives, and what we cannot claim. Do not recommend choosing a closure, target weight or preference change merely to reverse the sign. Cite literature only when you know it supports the claim, and distinguish a general intuition from an established theorem.

The author is technically trained in economics but does not want code jargon to substitute for the mechanism. Define each object, use a worked example when helpful, and be explicit about whose choices and which prices are held fixed. Do not use “parity”; say children ever born, number of children or fertility.

## Appendix A. Complete reference target table

| Moment | Role | Target | Model | Gap | Weight | Loss Contribution |
| --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | normalization | 2.1 | 2.099999967793076 | -3.220692423866467e-08 | — | — |
| cps_childlessness | scored | 0.19827875100684264 | 0.20253298392882108 | 0.004254232921978435 | 35532.3042455214 | 0.6430813285978316 |
| cps_exactly_one | scored | 0.21365532522014702 | 0.21063334432983166 | -0.0030219808903153567 | 26952.820824310795 | 0.2461430919206547 |
| nchs_mean_age | scored | 25.976263860992496 | 25.956378033847628 | -0.01988582714486853 | 139.82806784479274 | 0.05529446706909022 |
| nchs_share30 | validation | 0.2492780130410667 | 0.22931081311513798 | -0.019967199925928714 | 0.0 | 0.0 |
| wealth_earnings | scored | 6.92658379107299 | 6.488652892252826 | -0.4379308988201638 | 7.595098472533724 | 1.4566143563186387 |
| bequest_wealth | scored | 0.007291023472616158 | 0.006659276336892266 | -0.0006317471357238915 | 5165289.256198346 | 2.061489894087505 |
| old_dispersion | validation | 3.51593508651872 | 2.975984443818792 | -0.5399506426999281 | 0.0 | 0.0 |
| mean_rooms | scored | 5.729434240102641 | 5.977125401038308 | 0.24769116093566712 | 128.02070205233477 | 7.854186724098859 |
| ownership_30_55 | scored | 0.6762604168538028 | 0.6569654337652767 | -0.019294983088526063 | 2339.3623724673616 | 0.8709361249670908 |
| first_birth_rooms | scored | 1.465 | 1.2231930952117365 | -0.2418069047882636 | 137.5652749002964 | 8.043521301678819 |
| family_rooms | validation | 0.38509964969278165 | 0.3004909671601643 | -0.08460868253261733 | 0.0 | 0.0 |
| recent_parent_ownership | scored | 0.12760836356692162 | 0.11837410900018575 | -0.009234254566735878 | 27055.822957508266 | 2.3070894548319165 |
| early_fertility | scored | 0.8095276384290021 | 0.5312175503600474 | -0.27831008806895474 | 100.0 | 7.745650512094934 |

## Appendix B. Complete reported parameter table

There are 31 reported fields, not 31 independently estimated parameters. The historical calibration searched ten coordinates; the experiments fix all reference parameters except the explicitly varied financing rule. Bounds and near-bound flags below are retained report metadata, not newly imposed policy bounds. Purchase LTV overrides are recorded separately in the experiment receipts and cannot be recovered from the single generic `financed_share` row alone.

| Parameter | Estimate | Lower | Upper | Near Bound |
| --- | --- | --- | --- | --- |
| H0 | 6.293507689200028 | 0.2 | 80.0 | False |
| beta_annual | 0.9670100471250203 | 0.94 | 0.99 | False |
| chi | 1.093187485240681 | 0.1 | 5.0 | False |
| first_birth_fixed_cost | 0.4009349064519595 | 0.0 | 8.0 | False |
| kappa_fert | 0.13056430591307258 | 0.02 | 50.0 | True |
| kappa_fert_continuation | 0.364420274899206 | 0.02 | 50.0 | True |
| theta0 | 0.07862721952114816 | 0.0 | 8.0 | True |
| delta_alpha_jump | 0.0 | — | — | — |
| child_benefit_curvature | 0.09689406390932934 | 0.0 | 0.8 | False |
| tenure_choice_kappa | 0.012988772333020718 | 0.001 | 0.1 | False |
| psi_child | 0.17156192800028292 | 0.01 | 0.5 | False |
| child_benefit_CRRA_coefficient | 0.15493859558421574 | — | — | — |
| theta1 | 0.008193084126995582 | — | — | — |
| sigma | 2.0 | — | — | — |
| alpha_cons | 0.733 | — | — | — |
| delta_alpha | 0.0 | — | — | — |
| h_P | 2.3 | 0.1 | 2.3 | True |
| utility_reference_rent | 0.11046592704873838 | — | — | — |
| q_annual | 0.020000000000000018 | — | — | — |
| financed_share | 0.8 | — | — | — |
| housing_supply_elasticity | 0.63 | — | — | — |
| payroll_tax | 0.08028070961950022 | — | — | — |
| pension_period | 0.917784047463731 | — | — | — |
| annual_depreciation | 0.01416143718381309 | — | — | — |
| period_depreciation | 0.05545379079326218 | — | — | — |
| annual_property_tax | 0.010598360773872594 | — | — | — |
| period_property_tax | 0.042393443095490375 | — | — | — |
| selling_cost | 0.06 | — | — | — |
| rental_cap | 6.0 | — | — | — |
| wealth_grid_nodes | 120.0 | — | — | — |
| income_states | 9.0 | — | — | — |

## Appendix C. Source map and limitations

This document was prepared on October 1, 2026 from the current experiment receipts and model sources. It is a research-review handoff, not a replacement for the author-owned manuscript. The September 15 model-review prompt was used only as a format precedent: its old utility, earnings, finance, entry and calibration statements are not presumed current.

The source map below uses repository-relative paths for portability. The local repository root is `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

- `code/model/refactor_lab/engine/household.py`
  - SHA-256: `2a34f5f28c0c63ca3d24d9e759cd1b5b92aadaddb5fe04795e65abbc0b448082`.
- `code/model/refactor_lab/engine/kernels.py`
  - SHA-256: `fe7d43afdac234af71edd09e1260666e21c309dda83de4adc6a3abff35c3d5d4`.
- `code/model/refactor_lab/engine/shared.py`
  - SHA-256: `1fbe36f75befe5fed67c490c39a65663f4c8896b078851693fb4bca533d07773`.
- `code/model/refactor_lab/engine/child_preferences.py`
  - SHA-256: `887763b066bdd64efb2f23b97fac977abfa9630d1e194adb2148da24b2ff7123`.
- `code/model/refactor_lab/engine/adult_entry.py`
  - SHA-256: `e752178ad90d64c77575b70d0718a94d6f3b04463b03637d9af12e23f57bce3f`.
- `output/model/fixed_reference_economics_20260928/utility_floor_round2_v1/inputs.py`
  - SHA-256: `0306ba05089e7c0b16b21380690f2b6fe749efb7585d20db00f3dd2c1bad4285`.
- `output/model/fixed_reference_economics_20260928/credit_ge_v1/run_ge.py`
  - SHA-256: `59fa315ef6d50aaa8be2abcfb29f01cf406f3756e363a0bd088e473a424311b8`.
- `output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/purchase_ltv_v1/local_run/RESULTS.md`
  - SHA-256: `308c826b3baeb5ebfe401ab6ffce1865b30f1c5ef0f20e09aeb6707d0d5889f7`.
- `output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/collected/RESULTS.md`
  - SHA-256: `0c9f602ab0c24dac7a3bbc02c2bc5cdbfae31a63c65f376009665fc31b01a114`.
- `output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/diagnosis/renter_only_credit_v1/RESULTS.md`
  - SHA-256: `17d3035a9d701dd297a5c3a7fc3340c49a766f7898a6c4485643c94568731865`.
- `output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/diagnosis/wealth_composition_saved/README.md`
  - SHA-256: `9d2a4c4c4c38d967a6bc2af56922c5d49e714e412e8699dba3d7e8d59347e7e1`.
- `output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/diagnosis/wealth_composition_saved/lifecycle/README.md`
  - SHA-256: `62d49d32d0e71246baaf1e668bac3da62e22448365c1bdc23b384e93dc160486`.
- `output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/diagnosis/renter_only_credit_v1/results_20261001T223416Z/baseline_d0/target_fit.csv`
  - SHA-256: `80edf254f850d3ed2ad55f2c702316a2ca9dc378e56e23c639ef2bcc2bffd3e2`.
- `output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/diagnosis/renter_only_credit_v1/results_20261001T223416Z/baseline_d0/parameters.csv`
  - SHA-256: `ae092fafab15985f773b61999a29d4a8d505495623e615440cfff1f107c75406`.

Known limitations: no mortgage transition is completed in this packet; the experimental calibration is not adopted and optimizer convergence is not certified; full natural borrowing support and grid convergence are not certified; an independently active child-preference transition project is not evidence about the mortgage intervention; full empirical estimator provenance for every target is outside this packet, so do not infer identification merely from the target list. The preserved historical housing-market plot can display an unscaled-cohort residual in its title even when the population-scaled market clears; use closure receipts for the equilibrium residual.
