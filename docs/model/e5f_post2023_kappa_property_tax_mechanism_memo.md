# Why the property-tax reform moves housing much more than fertility

## Executive answer

The current result is not “nothing happens.” The housing equilibrium moves a great deal. Raising the rebated annual property-tax rate from 1% to 2% lowers the asset price by roughly 3–5%, raises the rent-equivalent user cost by roughly 17–27%, lowers aggregate ownership by roughly 2–3 percentage points, and lowers housing services per household head by roughly 6–8%. Annual births per household head nevertheless rise by only 0.16–0.60% along the transition.

That combination gives a sharper diagnosis than “the fertility mechanism is weak.” The reform currently delivers the wrong *housing first stage* for the Coven et al. narrative. It contracts housing services and aggregate ownership; it does not generate their large capitalization-driven expansion of housing access for young constrained buyers. The small positive birth response is therefore likely a net income effect from the equal rebate, after higher housing user costs partially offset it, rather than a large “cheaper homes permit children” response.

There is also a calibration issue. The model contains a housing cost of children, but the active fit underpredicts the observed increase in rooms after a first birth by 0.34 rooms (0.380 in the model versus 0.720 in the data) and overpredicts average occupied rooms by 0.93 rooms (6.710 versus 5.780). Prospective parents therefore begin with too much housing and add too little when the first child arrives. The calibrated first-child housing floor is only 0.515 rooms, while young parents in the baseline consume roughly 5 rooms as renters and 8–9 rooms as owners. No positive-mass young-parent group is numerically at the child-space floor in any inspected year.

The correct conclusion at this stage is:

1. The model has a theoretical housing–fertility link, but the counterfactual elasticity of fertility with respect to housing costs is not identified tightly by the current target set.
2. The present positive-tenure-taste run is an informative sensitivity, not a recalibrated E5f model. It borrows `tenure_choice_kappa = 0.0100175` from the older M5 calibration while holding the E5f parameters fixed.
3. Re-estimating that tenure-taste scale jointly may change incidence and smooth the response, but it cannot by itself create the missing Coven first stage or a strong fertility transmission.
4. Before adding a new family-specific tenure preference, we should determine whether the missing result comes from capitalization, financial constraints, age reallocation, or the child-housing block. The decomposition in this packet is designed to make exactly that distinction.

## 1. What was run

The experiment compares two deterministic perfect-foresight paths with the same reconstructed 2023 household and person distributions:

- a rebated 1% annual property tax;
- a rebated 2% annual property tax.

Both paths use 128 four-year dates (2023–2531), the same terminal demographic closure, the same housing-supply elasticity of 1.75, and a positive post-2023 tenure-choice logit scale of 0.0100175. The 1% path has maximum market and fiscal residuals of $1.4942\times 10^{-4}$ and $1.4131\times 10^{-5}$. The 2% path has maximum market and fiscal residuals of $1.1567\times 10^{-4}$ and $1.6511\times 10^{-5}$. Both pass the official $2\times 10^{-4}$ market and $2.5\times 10^{-5}$ fiscal gates, history reproduction, distribution feasibility, probability, accounting, and H128 terminal checks.

The 2% scheduler allocation reached its one-hour limit after the complete two-evaluation packet, plots, and summary had already been written. The numerical driver’s batch step exited successfully; the Slurm allocation was marked `TIMEOUT` because the outer allocation remained alive through the wall-time boundary. This is a scheduler-cleanup artifact, not a failed equilibrium evaluation. The immutable output hashes and all scientific gates were audited independently.

The positive tenure-taste scale deserves careful labeling. E5f itself was estimated with deterministic tenure choice, `tenure_choice_kappa = 0`. The value 0.0100175 was jointly estimated in the earlier M5 target system. Here it is imposed from 2023 onward on the E5f household block without re-estimating the other E5f parameters. These paths therefore answer, “What does the E5f policy experiment look like after adding a plausible amount of tenure smoothing?” They do not yet answer, “What is the fully recalibrated positive-taste E5f counterfactual?”

## 2. The aggregate result

All entries below are the 2% path relative to the 1% path.

| Year | Asset price | Rent-equivalent price | Ownership (pp) | Housing/head | Annual births/head | Resident persons | Equal rebate |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2023 | -3.196% | +27.063% | -2.299 | -5.525% | +0.158% | +0.000% | +82.906% |
| 2035 | -4.574% | +19.494% | -2.864 | -7.864% | +0.321% | +0.021% | +75.838% |
| 2051 | -4.653% | +18.191% | -2.605 | -8.004% | +0.415% | +0.068% | +75.428% |
| 2079 | -4.344% | +17.314% | -2.230 | -7.572% | +0.605% | +0.180% | +76.824% |
| 2103 | -4.224% | +17.772% | -2.167 | -7.515% | +0.535% | +0.319% | +77.156% |
| 2355 | -4.163% | +18.930% | -2.629 | -7.568% | +0.302% | +0.476% | +77.164% |
| 2531 | -4.164% | +18.944% | -2.633 | -7.538% | +0.307% | +0.434% | +77.226% |

![Paired transition effects](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/paired_transition_effects.png)

The birth response peaks at 0.605% in 2079. The largest ownership decline is 2.936 percentage points in 2031. The largest asset-price decline is 4.690% in 2043.

A worked numeric example resolves the earlier unit confusion. In 2079, annual births per household head rise from 0.0152279 to 0.0153200. The level change is therefore 0.0000921 annual births per head, and $0.0000921/0.0152279=0.00605$, or 0.605%. The number near 0.000006 discussed during the solver work was a fiscal-equilibrium residual, not a fertility effect. It measured numerical convergence and has no direct demographic interpretation.

Three observations matter.

First, the capitalization channel is present but modest. A 3–5% asset-price decline is economically meaningful, yet much smaller than the 12.6% decline in the current Coven et al. California counterfactual.

Second, the relevant flow price rises sharply. Doubling the property tax lowers the capitalized asset value but raises the annual rent-equivalent cost of housing. Households react by consuming substantially less housing. Thus “houses are cheaper to buy” and “housing is cheaper to use” are not the same statement in this experiment.

Third, the equal rebate rises by roughly 75–83%. The government returns the additional property-tax revenue equally across household heads. That is a broad positive income effect, not a transfer targeted to young constrained households or prospective parents. A small positive fertility response can therefore coexist with a large contraction in housing services.

## 3. Exact fertility decomposition

For each selected year, let $s$ index the complete at-risk state: liquid wealth, origin tenure, age, income, parity, and children currently at home. Let $w_s^j$ be state mass per household head and let $p_s^j$ be expected child units produced in the four-year period, where $j=0$ is the 1% path and $j=1$ is the 2% path. Annual births per household head are

\[
F^j = \frac{1}{4}\sum_s w_s^j p_s^j.
\]

The exact symmetric two-factor decomposition is

\[
F^1-F^0
= \underbrace{\frac{1}{8}\sum_s (w_s^0+w_s^1)(p_s^1-p_s^0)}_{\text{within-state fertility-policy effect}}
+ \underbrace{\frac{1}{8}\sum_s (p_s^0+p_s^1)(w_s^1-w_s^0)}_{\text{distribution/composition effect}}.
\]

The first term asks whether the *same type of household* becomes more or less likely to have a child. The second asks whether the reform moves mass toward states that already have high or low fertility. This is an accounting decomposition, not a claim that the two components are separately structural causal channels.

| Year | Total birth effect | Within-state fertility policy | Distribution/composition |
|---:|---:|---:|---:|
| 2023 | +0.158% | +0.158% | -0.000% |
| 2035 | +0.321% | +0.232% | +0.089% |
| 2051 | +0.415% | +0.244% | +0.170% |
| 2079 | +0.605% | +0.258% | +0.347% |
| 2103 | +0.535% | +0.227% | +0.308% |
| 2355 | +0.302% | +0.175% | +0.126% |
| 2531 | +0.307% | +0.175% | +0.132% |

![Age decomposition of the birth response](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/fertility_age_decomposition.png)

The decomposition resolves the main quantitative-versus-theory question. The same-household behavioral response is positive, so the housing–fertility link is not absent or wrong-signed. It is, however, small: holding the normalized distribution fixed, the reform raises annual births per head by only 0.16–0.26% across the inspected dates. At the 2079 peak, the within-state term is 0.258%, while 0.347%, or 57% of the total response, comes from changes in the distribution of households across wealth, tenure, age, parity, and child states. The initial 2023 composition term is zero to numerical precision because both policies begin from the same reconstructed distribution.

The age decomposition is also informative. In 2079, ages 25–34 contribute 0.338 percentage point of the 0.605% aggregate response, ages 35–44 contribute 0.176 point, and households under 25 contribute 0.091 point. These are contributions to the aggregate birth rate, not group-specific percentage changes, and they add exactly to the total.

Tenure reveals strong offsetting incidence. In 2079, renter-origin states contribute +0.699% of baseline aggregate births, while owner-origin states contribute -0.094%. At later dates, much larger positive renter-composition terms offset large negative owner-composition terms. The wealth-quartile decomposition has similar offsetting movements. This warns against describing the net positive birth effect as one clean housing-price channel: it is a general-equilibrium mixture of transfers, tenure changes, wealth transitions, and birth behavior.

Most importantly, the reform does not deliver more child-relevant space. By 2079, ownership among childless 25–34-year-olds rises 2.379 percentage points, which superficially resembles Coven et al.’s young-buyer channel. Yet their mean housing services fall 7.437%. Among 25–34-year-olds with dependent children, ownership falls 1.081 percentage points and housing services fall 6.200%; among 35–44-year-old parents, ownership falls 2.943 points and housing falls 6.590%. The policy can move some young childless households from renting to owning while simultaneously shrinking the amount of housing they consume. It changes tenure, but it does not relax the family-space margin.

![Young-household ownership and housing-service incidence](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/young_household_incidence.png)

The fertility cutoff is not densely populated either. No positive-mass at-risk state lies within 0.05 value units of birth indifference. Only 2.3–2.6% of at-risk mass lies within 0.25 units, and about 13–14% lies within 0.50 units. Combined with zero numerical child-floor binding, this explains why even large price, rent, and tenure changes move relatively few birth probabilities by much.

![Child-space floor binding is zero in every inspected young-parent group](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/child_space_floor_binding.png)
## 4. What Coven et al. establish

The current August 2026 version of Joshua Coven, Sebastian Golder, Arpit Gupta, and Abdoulaye Ndiaye’s [“Property Taxes and Housing Allocation Under Financial Constraints”](https://abdouecon.github.io/research/papers/Property_Tax.pdf) is a paper about housing allocation across generations. It does not contain endogenous fertility, and the text does not estimate or simulate a fertility effect. “Young families” in that paper means young households facing housing-market constraints, not households choosing whether to have a child.

Their mechanism has two complementary pieces:

1. **Embedded leverage.** A higher property tax is capitalized into a lower purchase price. With an 80% loan-to-value ceiling, a 12.6% price reduction lowers the required down payment dollar-for-dollar times the 20% equity requirement. This disproportionately helps low-wealth young buyers.
2. **Old-owner holding cost.** The same higher recurring tax raises payment burdens for incumbent owners. Their payment-to-income constraint, moving costs, capital-gains treatment, and bequest motive generate age-dependent lock-in and downsizing.

The model is built to make those channels quantitatively important. It includes a 20% down-payment requirement, a 36% payment-to-income limit, a minimum owner-occupied house size, transaction costs, mortgages, segmented owner and rental housing stocks, and a California housing-supply elasticity of 0.232. In the main counterfactual, the property tax rises from 0.8% to 2.0%, house prices fall 12.6%, aggregate ownership falls 2.1 percentage points, ownership among ages 25–34 rises 0.6 percentage points, ownership among ages 35–44 is roughly unchanged, and ownership among ages 65–74 and 75+ falls 5.0 and 7.6 percentage points.

The aggregate ownership response is therefore not the headline. The headline is the *age reallocation of owned housing*. Their empirical evidence is also about that allocation: high-tax areas have more owner-occupied bedrooms held by younger households, less crowding among people aged 20–40, younger migration into high-tax areas, and older downsizing after reassessment shocks. Their two-region extension strengthens young ownership because migration adds another adjustment margin.

An especially useful decomposition in their paper holds the tax system fixed and imposes the counterfactual price and rent paths. That experiment raises aggregate ownership by 14.3 percentage points. Yet the ownership increase is larger for older households, so capitalization alone does not generate the desired generational reallocation. The full result needs both price relief for constrained buyers and higher recurring holding costs for older owners.

This distinction is central for our paper. Coven et al. provide a housing-allocation first stage and a theory of which ages receive the housing. They do not establish that the reallocation changes fertility. Our paper’s potential contribution is to add that missing transmission and show when it is large, small, or offset.

## 5. Why our first stage is different

The comparison below is diagnostic rather than like-for-like: Coven et al. report stationary California age groups, while our numbers come from a national-style transition sensitivity. Its purpose is to locate the missing mechanism, not rank the models.

| Object | Coven et al. | Current E5f sensitivity |
|---|---:|---:|
| Property-tax change | 0.8% to 2.0% | 1.0% to 2.0% |
| Housing-supply elasticity | 0.232 | 1.75 diagnostic value |
| Asset-price response | -12.6% | -3.2% initially; -4.7% peak |
| Aggregate ownership | -2.06 pp | -2.30 pp initially; -2.94 pp peak |
| Payment-to-income limit | 36%, active | Inactive |
| Owner/rental stocks | Segmented | One housing-service market |
| Spatial migration | Two-region extension | Absent |
| Young/family incidence | Ownership +0.62 pp at ages 25–34 | Childless ages 25–34 eventually own more, but parents own less and every young group consumes less housing |

### 5.1 Capitalization is about one-third as large

Our diagnostic supply elasticity is 1.75, compared with 0.232 in Coven et al.’s California calibration. More elastic supply means that the housing stock contracts more and the asset price capitalizes less. Our price decline is 3–5%, versus 12.6% in their main experiment. The amount of down-payment relief is correspondingly smaller.

For a starter owner rung of two housing units and a 20% down payment, a 4% decline from an asset price around 0.62 reduces the required down payment by only

\[
0.20\times 2\times (0.04\times 0.62) \approx 0.010
\]

model income units. That is unlikely to move many households across the purchase threshold unless substantial mass is already concentrated very close to it.

### 5.2 The model has no active payment-to-income constraint

The current E5f block has the 20% down-payment requirement but sets `use_pti_constraint = False`. A higher property tax therefore cannot mechanically make an old owner fail a mortgage-payment-plus-tax-to-income test. Coven et al. use precisely that cash-flow constraint to generate downsizing among older owners. Our old-owner response must instead come through preferences, user costs, and the budget constraint, which is a weaker and less sharply age-targeted mechanism.

### 5.3 Owner and rental housing are too substitutable for the family mechanism

The model lets households obtain housing services as renters. Ownership has a general utility premium, but there is no child-specific benefit of tenure such as residential stability, school access, control over space, or protection from displacement. Consequently, a tax-induced change in ownership need not change the housing services available to prospective parents. Coven et al. separately clear owner and rental stocks; our one-market closure cannot reproduce all of that reallocation.

### 5.4 There is no spatial sorting margin in this run

The current path has one housing market. Coven et al.’s evidence and two-region extension emphasize migration and spatial sorting. If high-tax places free up large family units but prospective parents can neither move toward them nor select into a different local school/amenity bundle, an important family-relevant margin is absent.

## 6. Why the housing–fertility block produces a small effect

### 6.1 Housing and fertility are linked through a difference in indirect values

The model does not say that any movement in aggregate housing must generate a large fertility response. An eligible household chooses whether to attempt a birth and then optimizes consumption, housing, tenure, and saving within each fertility branch. The birth decision depends on the difference between the optimized “attempt” and “wait” values. A higher housing price affects both branch values. It changes fertility only through the *difference* in their housing exposures.

This is an envelope-theorem point with economic content. If the attempt and wait branches choose similar interior housing bundles, a common change in housing cost lowers both values by almost the same amount and largely cancels from the birth value gap. The response becomes first-order when a child creates a sufficiently large extra housing requirement, when one branch is at a housing-access constraint and the other is not, or when tenure delivers a child-specific service. The measured absence of floor binding and the small first-child housing shift therefore explain why a large aggregate housing response need not translate mechanically into births in the current architecture.

### 6.2 Birth-to-housing is not the same derivative as housing-to-birth

The calibrated moment “rooms after a first birth” asks how housing changes conditional on a birth. The policy counterfactual asks how the probability of a birth changes when housing prices, rents, ownership, wealth, and transfers all move. Those are different total derivatives. A model can match the first and still have a small second if households can substitute renting, reduce nonessential housing, or absorb the price change through consumption and saving.

More importantly, the current fit does not match the first derivative well. The data target is a 0.720-room increase after a first birth; the model produces 0.380. The model therefore under-delivers the empirical housing adjustment at exactly the family-formation margin we want to use.

### 6.3 Prospective parents have too much housing slack

The model’s mean occupied rooms are 6.710 against a target of 5.780. In the accepted 1% path, young parents consume about 5.0–5.3 housing units when renting and 8.0–8.9 when owning. The calibrated first-child subsistence shift is

\[
\bar h_{1} = 0.19945 + 0.31527 = 0.51471.
\]

No positive-mass young-parent group is at the numerical child-space floor at any of the seven inspected dates. The Stone–Geary shift still costs resources even away from the corner, so it is not literally irrelevant. But it is small relative to the housing bundle already consumed. A policy can remove several percent of surplus housing without making the child floor decisive.

### 6.4 The current target system does not identify the policy elasticity

E5f uses 11 free transition parameters and 12 dated targets. It targets total fertility, childlessness, age at first birth, the share of first births after age 30, the housing response to a first birth, room differences by parity, ownership, the family ownership gap, wealth, bequests, and average rooms. It does not contain a moment measuring fertility responses to exogenous housing-cost, house-price, down-payment, property-tax, or rent variation.

The housing-cost semi-elasticity of fertility is therefore a model-implied cross-derivative, not a directly disciplined object. The architecture gives it a sign under some local conditions; the current data do not pin down its magnitude.

### 6.5 The equal rebate pushes in the opposite direction from housing costs

An earlier, deterministic-tenure 2023 partial-equilibrium Shapley exercise is not directly comparable to the current positive-taste perfect-foresight path, but it is informative about signs. Holding the 2023 pre-fertility distribution fixed, the higher tax-rate component lowered births per adult by 1.086%, the asset-price component raised births by 0.099%, and the equal-rebate component raised births by 1.326%. The three components summed to a small positive response.

The current dynamic path shows the same qualitative ingredients: higher rent-equivalent costs and lower housing consumption, modest price relief, and a much larger equal transfer. The small positive fertility result should therefore not be described as evidence that cheaper housing strongly raises births. It is the net of offsetting housing-cost and income channels.

## 7. Quantitative problem or theory problem?

The answer depends on which row of the mechanism chain fails.

| Diagnostic | Interpretation |
|---|---|
| Young ownership and young housing rise, but within-state birth probabilities barely move | The housing first stage works; the fertility block or its calibration is quantitatively too weak. |
| Aggregate ownership moves, but young/parent housing services do not rise | The model reallocates tenure or financial claims, not child-relevant space. This is mainly a Coven-first-stage problem. |
| Within-state birth probabilities move, but distribution effects offset them | The mechanism exists at the household level; aggregation and endogenous sorting hide it. |
| Lower-wealth young renters remain far from owner entry after capitalization | The price decline is too small or financial constraints are incompletely represented. |
| Child-space costs are small relative to baseline housing and the first-birth room response remains underfit | The model’s housing–fertility complementarity is weak at the relevant states. |

The aggregate paths already establish two facts. Housing moves substantially, so the null is not caused by a frozen housing block. But the model’s response is a contraction in housing services, not Coven et al.’s reallocation toward young owned space. The state decomposition will decide whether there is an additional fertility-block failure after conditioning on that first stage.

## 8. What a fully recalibrated positive tenure taste can and cannot do

Bringing back a positive tenure taste is sensible. Deterministic tenure choice created discontinuous perfect-foresight mappings, while the positive scale produces smooth, convergent paths. Jointly re-estimating it inside E5f could move the mass of households near tenure margins and change age-specific incidence.

It is not, however, a free mechanism parameter for fertility. The older M5 estimate came from a different 15-moment system and a different parameter vector. Transplanting only its tenure scale into E5f is not internally calibrated. More importantly, the scale smooths the owner-versus-renter choice; it does not add a payment-to-income constraint, create spatial migration, lower the supply elasticity, make ownership child-specific, or identify the fertility response to housing costs.

My prior is therefore that full positive-taste recalibration can change the effect meaningfully but is unlikely by itself to turn a 0.2–0.6% birth response into a several-percent response. A several-percent response requires either a much stronger first stage for young family housing or a stronger, empirically disciplined transmission from housing access to birth timing.

## 9. Recommended sequence

### Step 1: finish the current decomposition

Use the paired exact state tables to report, by age, wealth, origin tenure, parity, and child state:

- the within-state change in birth probabilities;
- the distribution/composition contribution;
- young-parent ownership and housing-service changes;
- mass near the tenure and fertility margins;
- child-space floor incidence.

This is the cheapest way to locate the broken link.

### Step 2: run first-stage sensitivities, not a broad recalibration

Run a small predeclared set of one-at-a-time structural diagnostics:

1. lower the supply elasticity from 1.75 toward the California value 0.232;
2. introduce the already-coded payment-to-income constraint with an explicitly calibrated limit;
3. test segmented owner/rental supply or, at minimum, a reduced-form wedge that reproduces the owner-versus-renter price response;
4. measure the age profile of ownership and owned housing, not only the aggregate ownership rate.

The goal is to reproduce the Coven first stage before asking fertility to amplify it.

### Step 3: repair identification of the family transmission

Re-estimate the child-housing block so the first-birth room response is not missed by 0.34 rooms and average rooms are not high by 0.93. Add a moment that identifies the reverse link: a fertility or birth-timing response to plausibly exogenous housing costs, property taxes, house prices, rents, down-payment relief, or housing supply shocks. Without such a moment, a strong policy fertility effect would be a modeling assumption rather than a calibrated result.

### Step 4: then recalibrate positive tenure taste jointly

Estimate `tenure_choice_kappa` jointly with the E5f parameters under the revised target contract. Report its bound, identifying moments, and covariance/substitutability with the ownership premium, child-space parameters, and wealth process. Only then treat the positive-taste policy path as the production counterfactual.

### Step 5: add new theory only if the first three steps demand it

If young households gain owned housing services and the empirically disciplined fertility block still does not respond, then a missing child-specific value of ownership becomes a serious theory candidate. Possible objects include residential stability, school access, landlord restrictions, moving risk during pregnancy/early childhood, and the ability to customize space. Those mechanisms should enter only with corresponding empirical moments; otherwise they can mechanically manufacture the desired answer.

## 10. Bottom line for the paper

The Coven et al. mechanism is real and useful, but its current quantitative result is an age reallocation of housing, not a fertility result. Our model can make a strong contribution by showing when that reallocation reaches prospective parents and changes family formation. The current run says that this transmission is small under the present closure.

If that finding survives the exact state decomposition and the first-stage sensitivities, the smart explanation is not “housing does not matter for fertility.” It is:

> Property taxes reorganize the price, financing, tenure, and age allocation of housing. Fertility responds only when the households gaining child-relevant space are also near a birth margin. In the present calibration, the reform raises recurring housing costs, contracts housing services, rebates revenue broadly, and provides too little capitalization relief to young constrained buyers. The housing and fertility margins therefore do not overlap enough to generate a large aggregate birth response.

That is a coherent economic result. But because the current calibration underfits first-birth housing adjustment and does not directly identify the housing-cost elasticity of fertility, it is not yet the final result we should defend.

## Appendix A. Complete active target fit

The active E5f transition loss is 36.0992. The table reports every dated target used in that objective. “Std. gap” is the gap in the target’s uncertainty units; the loss contribution is the weighted squared gap.

| Moment | Target | Model | Gap | Weight | Std. gap | Loss contribution |
|---|---:|---:|---:|---:|---:|---:|
| Total fertility rate | 1.9180 | 1.9183 | +0.0003 | 1,425.74 | +0.009 | 0.000 |
| Childless rate | 0.1880 | 0.1859 | -0.0021 | 17,180.74 | -0.280 | 0.078 |
| Mean age at first birth | 26.0446 | 26.3273 | +0.2827 | 44.44 | +1.885 | 3.552 |
| Share of first births at age 30+ | 0.2603 | 0.2421 | -0.0182 | 10,000.00 | -1.821 | 3.317 |
| Rooms gained after first birth | 0.7202 | 0.3798 | -0.3405 | 137.57 | -3.993 | 15.946 |
| Rooms: parity 3+ minus parity 1–2, ages 30–55 | 0.3677 | 0.3593 | -0.0084 | 2,958.51 | -0.459 | 0.211 |
| Family ownership gap | 0.1677 | 0.1769 | +0.0092 | 14,229.59 | +1.098 | 1.206 |
| Ownership rate | 0.5755 | 0.5461 | -0.0294 | 1,207.85 | -1.022 | 1.045 |
| Mean occupied rooms, ages 18–85 | 5.7800 | 6.7101 | +0.9301 | 11.97 | +3.218 | 10.359 |
| Wealth / annual gross labor earnings | 6.8731 | 6.7664 | -0.1067 | 6.29 | -0.268 | 0.072 |
| Annual bequest flow / wealth | 0.0088 | 0.0086 | -0.0002 | 5,165,289.26 | -0.402 | 0.162 |
| Old-age wealth p90/p50, ages 76–84 | 3.4481 | 3.3963 | -0.0518 | 56.96 | -0.391 | 0.153 |

The two housing-size rows account for 26.304 of the 36.099 loss points. This is why the weak family-housing first stage is not a subtle residual concern: it is already the dominant visible calibration tension.

## Appendix B. Complete active free-parameter table

The active transition system has 11 free parameters and 12 dated targets. Only `theta1` is flagged near a search bound. The positive tenure-choice scale used in the policy sensitivity is not in this table because it was not estimated in E5f.

| Parameter | Estimate | Lower bound | Upper bound | Transform | Near bound? |
|---|---:|---:|---:|---|---|
| $\beta_{annual}$ | 0.994959 | 0.94 | 0.9995 | discount | No |
| `kappa_fert` | 2.267631 | 0.02 | 50.0 | log | No |
| `kappa_fert_continuation` | 1.765427 | 0.02 | 50.0 | log | No |
| `chi` | 1.033945 | 0.1 | 5.0 | log | No |
| `H0` | 16.381189 | 0.2 | 80.0 | log | No |
| `theta0` | 0.563219 | 0.0 | 8.0 | soft-zero | No |
| `theta1` | 0.099707 | 0.02 | 16.0 | log | **Yes** |
| `hbar_child_rooms` | 0.315267 | 0.1 | 1.8 | log | No |
| `first_birth_fixed_cost` | 4.545167 | 0.0 | 8.0 | soft-zero | No |
| `hbar_first_child_jump` | 0.199446 | 0.0 | 0.5 | soft-zero | No |
| `psi_child_change_2023` | -0.350000 | -1.5 | 0.2 | asinh | No |

For completeness, the derived or externally normalized fertility intercepts are `psi_child_2007 = 0.280341` and `psi_child_2023 = -0.069659`.

## Reproducible artifacts

The audited packet is [`output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/`](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/). Its main files are:

- [`summary.json`](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/summary.json), with input hashes and all path gates;
- [`fertility_decomposition_summary.csv`](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/fertility_decomposition_summary.csv) and [`fertility_decomposition_by_group.csv`](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/fertility_decomposition_by_group.csv);
- [`cross_section_incidence.csv`](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/cross_section_incidence.csv);
- the complete [`target-fit table`](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/active_calibration_target_fit.csv) and [`parameter table`](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/active_calibration_parameters.csv);
- [paired transition effects](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/paired_transition_effects.png), [age decomposition](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/fertility_age_decomposition.png), [young-household incidence](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/young_household_incidence.png), and [child-space floor incidence](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/child_space_floor_binding.png);
- [`artifact_hashes.json`](../../output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/artifact_hashes.json), which fingerprints every reported table and plot.

The complete 1% and 2% state diagnostics are immutable inputs. The collector validates path, source, terminal, and state-table hashes; equilibrium, history, accounting, feasibility, and probability gates; exact state-grid equality; and exact decomposition add-up before writing any reported output.
