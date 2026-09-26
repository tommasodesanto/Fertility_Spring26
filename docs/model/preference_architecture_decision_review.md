# Decision review: household needs, housing, and the benefit of children

September 25, 2026. Independent review for the author; no specification is adopted by this note. Small text reads, primary-paper checks, and arithmetic only: no model import, solve, empirical estimation, or job launch.

## Decision

**For the next working calibration, retain the equivalence scale, the parenthood housing requirement, and a separate linear flow benefit while children live at home. The first challenger should change only the benefit to a mildly concave function of children at home.** I would not restore two subsistence jumps, remove all child dependence from housing demand, put fertility inside the consumption composite, or add freely estimated benefit curvature to the present shock specification today.

This recommendation is based on the mechanisms and the evidence, not the amount of work already invested. The scale supplies an additional-child cost without requiring a fixed amount of nonhousing spending. The housing requirement supplies a distinct, nonhomothetic space-demand channel. A separate child benefit makes its assumptions visible. Linearity is a defensible restriction over the model's small child-count support; the current calibration has not established that it is the source of the important fit failures. The claim that intermediate family sizes are produced only by shocks is currently unverified.

**The strongest reservation is the housing requirement's hard lower support, not its coexistence with a scale.** Its interpretation must be a reduced-form family housing requirement in model services, not a measured universal bedroom standard. A material, reliably measured population of comparable parents below the implied threshold would reverse my ranking in favor of a smooth crowding cost. Absence of bunching at the threshold would not: Stone–Geary utility tends to minus infinity there and does not generally predict an atom at the minimum.

**Ranked runner-up:** the same material-needs specification with \(v(m)=m^{0.86}\) in place of \(m\), retaining both existing decision-shock scales and re-normalizing \(\psi\). The exponent is one mild, literature-anchored sensitivity restriction, not an externally identified estimate for this model. A bounded comparison is specified below. I prefer this to an immediate \(\log(1+m)\) baseline: the closest published estimates do not justify that particular amount of curvature.

## 1. What each proposed cost device buys

Let \(m\) denote children currently at home, \(n\) children ever born, \(c\) total nonhousing consumption, and \(s\) housing services. The maintained utility is

\[
u(c,s,m)=U\!\left(\frac{c^{\alpha}[s-h_P\mathbf1\{m>0\}]^{1-\alpha}}{e(m)}\right)+\psi m,
\quad U(x)=-x^{-1},\quad
e(m)=\left(\frac{2+0.7m}{2}\right)^{0.7},\quad\alpha=0.733.
\]

The [adapter](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_parenthood_utility.py:18) fixes the nonhousing intercept and housing slope at zero. The [solver](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2435) applies the scale, housing requirement, and child benefit to **current dependents**. Thus the requirement vanishes when the last child leaves and returns if another child subsequently arrives. It is not an irreversible payment charged only at the first-ever birth. The separate first-birth utility cost is a different object.

### Equivalence scale versus two jumps

I interpret “two jumps” as a parenthood step in required nonhousing spending and a parenthood step in required housing. That matches the economic question, but it is not a literal description of every historical parameter: the July predecessor also had childless intercepts and per-child slopes. The [July 18 memo](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/equivalence_scales_sequential_fertility_readings_20260718.tex) and [July 20 specification](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/eqscale_seq_specification_note_20260720.tex) document removing that committed bundle after its affordability constraint obstructed the intended income risk. They do not refute the author's recollection that housing requirements subsequently returned to strengthen the housing response; the exact later motivation is incompletely documented.

The economic distinction can be stated without solving the dynamic model. For an interior renter, let \(E_0\) be the expenditure required for a childless household to attain a given **material** living standard at the same rent \(r\). The alternative cost systems imply:

| Device | Expenditure needed for that material standard | Economic restriction |
|---|---|---|
| Scale only | \(e(m)E_0\) | Additional needs rise with the living standard; no positive subsistence bill is introduced. |
| Two pure parenthood jumps, no scale | \(E_0+c_P\mathbf1\{m>0\}+rh_P\mathbf1\{m>0\}\) | A fixed setup cost in each good; no additional material requirement for a second or third child while children remain at home. |
| Current scale plus housing requirement | \(e(m)E_0+rh_P\mathbf1\{m>0\}\) | Progressive household needs plus a disproportionately costly housing commitment for low-resource parents. |

These are compensated **material-needs** comparisons, holding the material standard fixed and excluding the pleasure of children. They are not indifference scales for total utility including \(\psi m\).

The gain from the scale is therefore substantial. It permits families to reduce nonhousing consumption when resources fall, while still making another dependent child costly. Replacing it by two pure parenthood steps would make the second and third child share all the modeled material commitments, an unnecessarily strong economies-of-scale assumption. Adding a nonhousing per-child floor would fix that omission, but would impose a new absolute expenditure commitment and require a measured magnitude. It is not a free simplification.

For example, with a childless material budget of 100 and a parenthood housing requirement costing 10, the current first-child requirement is approximately 133.4: 123.4 for the scaled discretionary bundle plus 10 for housing. A nonhousing jump of 20 and housing jump of 10 instead imply 130, but the same absolute 30 cost at a childless budget of 50 or 200. The scale's cost adjusts with the living standard. This is illustrative arithmetic, not an estimate.

There is a genuine calibration qualification: the imported scale is conventionally a **total-needs** scale. Once a separate housing requirement is added, the model's total needs ratio exceeds it. Calling the scale “nonhousing costs” does not remove that issue. Under Cobb–Douglas, scaling only \(c\) by \(e_c(m)\) simply produces an aggregate scale \(e_c(m)^\alpha\); choosing \(e_c=e^{1/\alpha}\) reproduces the current utility exactly. Any different coefficients are an economic change, not a new channel created by notation.

My recommendation is to describe \(e\) honestly as an imposed scale for the discretionary material bundle and validate the **combined** expenditure pattern. Do not add a second scale parameter that fertility tastes could offset. Ordinary conditional demand identifies allocations, not the household's utility ranking across family sizes: that distinction is the point of [Pollak and Wales (1979), pp. 216–221](https://news.fbc.keio.ac.jp/~hayami/pdf/consumption/PollakWales1979.pdf).

### Why add the housing requirement to a scale?

For a renter with fixed expenditure \(X=c+rs\), the scale alone cancels from the within-period allocation: housing is \((1-\alpha)X/r\). Children can still change housing through saving, total expenditure, tenure, and selection, but the scale alone does not tilt the conditional housing allocation. With the parenthood requirement,

\[
s=\alpha h_P\mathbf1\{m>0\}+(1-\alpha)X/r.
\]

It adds an absolute housing-demand intercept, and its housing-budget-share effect is largest at low expenditure. This is exactly the useful asymmetry between adaptable nonhousing spending and a larger housing commitment. The same feature causes the hard feasibility risk: parents must finance the requirement even after bad income shocks. Removing the nonhousing floor did not remove every potential affordability problem.

[Dustmann, Fitzenberger, and Zimmermann (2022), Eq. (1a), p. 1712](https://www.rfberlin.com/wp-content/uploads/2026/02/ueab097.pdf), provide a close housing precedent: \(n^{-\phi}(h-n^\phi\bar h)^a x^{1-a}\), with a household-size scale, a housing minimum, and no nonhousing minimum. Their German expenditure analysis motivates housing as a necessity. Family composition is exogenous, the minimum scales with household size, and fertility responses are explicitly outside their scope. This supports the combination's coherence, not our parenthood-only shape, U.S. parameter, or cardinal fertility utility.

A child-dependent housing share instead raises desired housing **in proportion to expenditure**. It is attractive if richer families' birth-related housing expansion scales with their budgets and poor parents readily accommodate children in small homes. It does not deliver the same absolute commitment. Moreover, changing Cobb–Douglas exponents across child states changes utility levels and the composite's units. An endogenous-fertility version needs an explicit common reference-price/expenditure normalization; matching budget shares alone does not supply it.

A smooth crowding penalty permits below-standard housing at a finite utility cost. It is the best replacement if the hard-support restriction fails. [Szabó's December 2022 working draft, p. 19](https://benceszabo.github.io/papers_currentversion/bence_szabo_fertility_housing_current.pdf), illustrates a crowding mechanism with two housing sizes, desired-family-size heterogeneity, and multiplicative utility. That is a possible mechanism, not a ready-made U.S. calibration. In particular, multiplying our negative CRRA term by a decreasing crowding factor would give the wrong sign; any transfer requires a fresh utility-level check. [Van Doornik et al., BCB Working Paper 612, pp. 13–14](https://www.bcb.gov.br/content/publicacoes/WorkingPaperSeries/WP612.pdf), instead impose \(n/h\leq s\) and explicitly mention continuous crowding costs as an alternative. Their Brazilian credit-lottery evidence supports a space mechanism, not our particular floor.

The present parenthood-only shape remains restrictive. It allows shared child space but gives later children no direct additional housing intercept. The family-room gap is consequently an informative check. I would not add a housing slope and benefit curvature together: both can reduce the incentive for additional children, while the present room observer does not yet cleanly identify the former.

## 2. Where the benefit of children should enter

**Keep it additively separate from material utility.** This makes the amount and curvature of the child benefit distinct from risk aversion and housing substitution. Linear \(\psi m\) says that, conditional on the period and a child's remaining time at home, another child contributes the same flow benefit. It does not say that every birth has the same lifetime value: age, remaining opportunities, resources, the housing requirement, and future child departures all matter.

Concave \(\psi v(m)\) supplies declining marginal enjoyment of additional **simultaneous dependents**. That is economically plausible and the best next alternative. It also gives a reason to space births to reduce overlap. It does not merely change completed-family-size preferences. Replacing \(m\) by \(n\) would continue paying benefits after children leave and change the duration of rewards. A one-time birth reward would also change the architecture unless its state- and age-dependent present value were matched. Neither timing change is recommended here.

The primary literature supports concavity as a modeling choice, but not a universal amount. [Sommer (2016), pp. 30–34 and Table 2](https://www.kamilasommer.net/Fertility.pdf), uses separable consumption utility plus \(\zeta(nq)^{1-\kappa}/(1-\kappa)\), where \(n\) means children at home and \(q\) is endogenous child quality. The estimated \(\kappa=0.14\) implies mild curvature, and goods, time, minimum quality, and fertility timing are jointly disciplined. Its quality minimum is not a housing floor. [Baudin, de la Croix, and Gobbi (2015), Eq. (1), p. 1860 and Table 3, p. 1867](https://perso.uclouvain.be/david.delacroix/pdfpubli/aer15.pdf), use \(\log(n+\nu)\), but estimate \(\nu=9.362\), making curvature mild over small families. That is quite different from simply selecting \(\log(1+m)\). Neither estimate transfers directly to our flow benefit and cost system.

Conversely, [Doepke and Kindermann (2019), quantitative-model preference section](https://www.nber.org/system/files/working_papers/w22072/w22072.pdf), use a linear birth reward with heterogeneous partner preferences and explicit child costs. Their result shows that linearity is not intrinsically defective; their partner heterogeneity, bargaining, and birth timing cannot be equated to our decision-specific logit shocks.

Putting children inside a curved composite is a larger and less identified change. If the bundle is multiplied by a child index \(b(m)\), then \(U(Qb(m)/e(m))\) is just a new effective scale \(e(m)/b(m)\). With Cobb–Douglas, this does not create an independently identified benefit channel. If the inside term is instead \(Q/e+\psi v(m)\), children become substitutes for material consumption inside the same risk aggregator, changing how the value of a child varies with consumption and risk. Those are substantive restrictions on complementarity, substitution, and insurance. They should not be selected merely to obtain curvature.

An outside family weight, \(e(m)U(Q/e(m))\), is also substantive. [Scholz, Seshadri, and Khitatrakun (2006), pp. 615 and 619](https://users.ssc.wisc.edu/~aseshadr/Publications/optimality.pdf), use that aggregation with exogenous family composition. We borrow their scale's shape, not their complete objective. With endogenous fertility, a family-dependent multiplier interacts with the utility-level normalization. Adding it just to turn declining cost increments into increasing ones would hide the fertility assumption in material utility.

### What the calibrated model presently establishes

At fixed material composite, the scale increments decline: 0.234, 0.216, and 0.203. That arithmetic is correct. The inference that total birth costs necessarily decline, deterministic choices have only endpoint family sizes, or intermediate families arise only from shocks is not. The [current incentives audit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/utility_fertility_rationale_review/current_incentives_audit.md) checks that the solver compares optimized birth and wait values, including housing, saving, risk, maturation, and future opportunities. No compact saved output establishes the mass-weighted distribution of those gaps.

Similarly, concavity alone does not guarantee an interior choice or a correct fertility–income gradient. With additive child benefits, even a homothetic material scale can make children more attractive as resources rise because the material-utility penalty falls. A concave benefit does not itself supply a wage-dependent opportunity cost. These are reasons to inspect the actual choices, not to redesign unrelated earnings or labor margins in this exercise.

The latest selected low-tax experiment matches normalized completed fertility at 2.100 but has childlessness 0.183 versus 0.198, exactly-one-child share **among mothers** 0.231 versus 0.214, and first-birth mean age 26.342 versus 25.976. These are real misses; they are not a visible collapse to only childless and maximum-size families. Mean rooms is too high while the birth-related room response and larger-family room gap are too low. That pattern does not point uniquely to insufficient benefit curvature.

The selected \(\psi=0.135\), \(h_P=1.890\), first-birth utility cost 0.507, and first-birth/continuation shock scales 0.209/0.482 must be interpreted together. Raw shock scales cannot be called large or small without the optimized gaps in the same utility units. The compact receipt omits some branch switches, so a run-specific logit-gap extraction must first confirm them. \(\psi\) is derived from the author-chosen 2.1 replacement normalization; its successful normalization is not empirical validation of fertility tastes.

The housing floor likewise is not the event-study coefficient. At these estimates the fixed-expenditure interior-renter formula gives \(\alpha h_P=1.385\) services, but the actual model observer gives 0.789 rooms against 1.465. Tenure, expenditure, anticipation, selection, and the observer all intervene. Owners' service premium also makes the physical threshold differ from the service threshold. Do not set \(h_P\) equal to the empirical room response.

The earlier same-objective comparison selected a floor case with loss 464.852 and a share case with loss 851.271; the latter nevertheless had a slightly larger birth-room response. It supplies evidence against a confident claim that shares already fit better. It does not establish the best attainable fit of either architecture: searches were bounded, some cases failed, and the principal arms lack selected repetitions. Those losses cannot be compared to today's 280.411 under a different target contract. Complete historical [fit tables](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/utility_overnight_v1/final_readout/selected_target_fits.csv), [parameter tables](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/utility_overnight_v1/final_readout/selected_parameters.csv), and [limitations](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/utility_overnight_v1/final_lead_review.json) are retained.

## 3. Calibration strategy and the smallest useful comparison

The current selected contract has eight searched coordinates and twelve weighted rows, plus the separately imposed fertility normalization. The material scale, consumption share, and risk aversion are externally fixed restrictions. Only three searched coordinates belong directly to the fertility block: first-birth utility cost and the two decision-shock scales. Parameter count alone does not establish local identification.

Keep those restrictions for the comparison. Do not freely estimate both a new curvature coefficient and a more flexible family scale against the same few fertility moments. At a common completed-fertility age, with the 3+ representative fixed, the mean, childlessness, and exactly-one share already determine the remaining two count bins. A 3+ row is then not independent identifying variation. Here the CPS stock window and the 2.1 completion normalization differ; their difference should be exposed as a late-fertility check, not exploited as an additional clean curvature moment.

For eventual empirical discipline, use the existing data roles precisely:

| Object | Available evidence and immediate use | What is still needed to identify it more directly |
|---|---|---|
| Combined scale plus housing requirement | CEX childless-renter share fixes \(\alpha\); existing consumption builders permit family comparisons. | Matched nonhousing and housing expenditure by resources, two-adult composition, and current dependents. This tests combined needs; conditional spending shares alone do not identify welfare-scale levels. |
| Absolute housing requirement versus smooth demand | Existing ACS rooms and family-room contrast; PSID birth-event evidence. | A matched lower-tail distribution of housing services/physical rooms by current dependents, income, rent, and tenure; plus room responses by pre-birth resources. Zero below-floor support is the restriction to confront. |
| Benefit curvature versus continuation shocks | CPS stocks and NCHS first-birth timing are already available. | Birth-order-specific age/hazard or spacing variation, with current children at home distinguished from ever-born children. Curvature shifts systematic higher-child-count incentives; shock scale changes sensitivity around those incentives. Existing first-birth timing alone is weak for this separation. |

**One proposed overnight comparison: linear versus mild concavity, with a small matched fertility-block refit.** Use \(v_0(m)=m\) and \(v_1(m)=m^{0.86}\), both zero at zero and one at one. The concave values at one, two, and three dependents are 1.000, 1.815, and 2.572. This preserves the duration, state definition, first-child benefit units, and number of free search coordinates. The number 0.86 is an explicit experimental restriction motivated by Sommer's mild curvature; it is not an imported U.S. estimate for our specification.

Freeze the same nonpreference launch contract in both arms, including the adopted pension rule for a new calibration. Declare that common fiscal update separately from the single experimental utility change. Keep the housing floor, equivalence scale, income, entry wealth, survival, bequests, housing market, targets, weights, and all numerical gates identical across the arms. Re-normalize \(\psi\) to 2.1 in each case and compare only the new arms to one another.

Use two identical anchor evaluations per arm to smoke-test the exact loop; then the same six predeclared positive/negative coordinate moves in first-birth cost and the two shock scales; then two exact repetitions of each arm's selected point. That is **at most twenty normalized objectives**, not twenty stationary solves. Keep the five other structural coordinates fixed. Select using the unchanged full objective and display every row; do not silently drop housing or wealth targets because this is a fertility comparison. This is a local architecture diagnostic, not a complete re-estimation or a convergence claim.

The recent full normalized objective took about ten minutes at the selected point, while the broader recent timing reference was about fourteen minutes. Four one-thread workers suggest roughly one to two hours of numerical work before overhead; use a four-hour global cap with one hour reserved for repetitions and export. Before launch, the controller must also pin the inherited normalization iteration limit and thus the maximum stationary-solve count; stop on the existing numerical failures rather than enlarging the budget. Coordinate resources with other Torch work. No computation is launched by this note.

The decisive output is the complete fit and parameter comparison plus two compact additions from each solved checkpoint: (i) mass-weighted optimized birth-minus-wait gaps divided by the relevant shock scale, by age and children already born/currently at home; (ii) age-cell first-birth flows and the child-count distribution at 40–44 versus completed ages. Recovering gaps from binary probabilities is valid only after confirming the branch switches; zeros and ones are censored. An argmax of saved values is not a zero-shock equilibrium because continuation values already include future shocks.

**Decision rule:** favor the concave challenger if its matched local refit improves the full economic fit, particularly count and timing restrictions, without shifting the 2.1 normalization into a worse late-birth tail or failing numerical checks. If it merely trades the same stock fit for a different shock scale, retain linearity and report that curvature is not selected by these data. If the two remain competitive, the next discriminating moment is a conditional higher-order birth hazard/spacing contrast, not another aggregate fertility mean. A failure of the housing-support restriction would instead change the housing block and warrants a separate smooth-crowding comparison; it is not a reason to bundle several preference changes tonight.

## Evidence attached to the current numerical claims

The September 25 low-tax selection below is **experimental**, uses the frozen 8.751% tax, and predates implementation of the adopted pension rule. It passed its saved gates and exact repetitions but is not a converged calibration. Its first-birth observer remains an unmatched stationary proxy for the empirical panel estimator. The [full-precision target table](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/payroll_tax_review/overnight_pair_20260925_preparation/final_results/oasi_087510/target_fit.csv), [all 25 parameter/restriction rows](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/payroll_tax_review/overnight_pair_20260925_preparation/final_results/oasi_087510/parameters.csv), and [receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/payroll_tax_review/overnight_pair_20260925_preparation/final_results/oasi_087510/receipt.json) are authoritative for this point. Values below are rounded to three decimals; gaps are model minus target and weights apply to the unrounded source units.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Completed-fertility normalization | 2.100 | 2.100 | \(3.007\times10^{-7}\) | — | — |
| Childlessness, ages 40–44 | 0.198 | 0.183 | −0.015 | 35,532.304 | 8.108 |
| Exactly one child among mothers, ages 40–44 | 0.214 | 0.231 | 0.017 | 26,952.821 | 8.236 |
| Mean first-birth age, model-cell convention | 25.976 | 26.342 | 0.366 | 139.828 | 18.755 |
| First births at 30+ | 0.249 | 0.247 | −0.002 | 13,866.065 | 0.063 |
| Wealth / gross earnings | 6.927 | 6.045 | −0.882 | 7.595 | 5.902 |
| Annual bequests / wealth | \(7.291\times10^{-3}\) | \(6.904\times10^{-3}\) | \(-3.867\times10^{-4}\) | 5,165,289.256 | 0.772 |
| Old-age dispersion | 3.516 | 2.958 | −0.558 | 10.616 | 3.306 |
| Mean rooms | 5.608 | 6.544 | 0.936 | 128.021 | 112.085 |
| Ownership, ages 30–55 | 0.676 | 0.757 | 0.081 | 2,339.362 | 15.161 |
| First-birth room response | 1.465 | 0.789 | −0.676 | 137.565 | 62.892 |
| Larger-family room contrast | 0.385 | 0.220 | −0.165 | 280.528 | 7.629 |
| Recent-parent ownership gap | 0.128 | 0.090 | −0.037 | 27,055.823 | 37.501 |
| **Total weighted loss** | | | | | **280.411** |

| Searched parameter | Estimate | Search interval | Bound assessment |
|---|---:|---:|---|
| Housing supply scale \(H_0\) | 8.646 | [0.200, 80.000] | Interior |
| Annual discount factor | 0.962 | [0.940, 0.990] | Interior |
| Owner service premium \(\chi\) | 1.127 | [0.100, 5.000] | Interior |
| First-birth utility cost | 0.507 | [0.000, 8.000] | Interior |
| Parenthood housing requirement \(h_P\) | 1.890 | [0.100, 2.300] | Interior |
| First-birth shock scale | 0.209 | [0.020, 50.000] | Broad 1%-of-range lower-proximity flag; not at bound |
| Continuation shock scale | 0.482 | [0.020, 50.000] | Same flag; not at bound |
| Bequest scale \(\theta_0\) | 0.088 | [0.000, 8.000] | Interior under the saved screen |

The fixed bequest shift is \(\theta_1=0.008\); \(\psi=0.135\) is derived rather than searched. The linked parameter table lists all remaining external, adopted, endogenous, and diagnostic settings. Neither the spread of the raw shock bounds nor their proximity flags establishes that shocks dominate fertility.
