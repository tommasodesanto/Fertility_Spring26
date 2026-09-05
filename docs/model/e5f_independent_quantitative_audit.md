---
title: "Housing, tenure, and fertility"
subtitle: "Independent quantitative audit and reconciliation"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026; original review dated 4 September"
lang: en-US
---

# Discussion guide

**The model is worth keeping. The current quantitative package is not yet ready for circulation as identified policy estimates.** The selected calibration is reproducible, and the supply-versus-credit contrast is informative. The first task is to establish that the relevant moments, numerical choices, and demographic objects support the claims the paper will make.

This document brings together the complete fourteen-section audit and the subsequent reconciliation with concurrent calibration work. It preserves the full target-fit table, all estimated parameters and bounds, equation-to-code findings, source receipts, policy results, historical validation, and proposed priorities. Section 15 preserves the September 4 reconciliation. **Sections 16-18 contain the overnight verification and supersede earlier pending-work descriptions.** The initial review was read-only; the overnight work includes a verified reconstruction, numerical experiments, and a separately labeled rental-cap economic sensitivity.

## Start our discussion with these five issues

| Issue | What is established | What we need to decide or verify |
|---|---|---|
| **1. Scope of the paper** | The active model is a national, one-location housing-access model. | Lead with conditional household mechanisms; decide separately whether resident-population, spatial, tax-basis, or welfare claims belong in a future extension. |
| **2. First-birth housing response** | Model 0.436 rooms versus target 0.720. Exact repeats preserve the miss. | Separate measurement, numerical, and economic explanations before changing targets or weights. |
| **3. Numerical and housing-opportunity sensitivity** | Global saving and the denser wealth grid preserve the birth contrast. Larger rentals leave supply births near 1.29% but reverse its ownership response. | Use the completed conditional checks; validate wealth support, history, and the tenure-specific opportunity set. |
| **4. Demographic handoff** | Historical household totals and age shares are imposed; the post-2023 queue changes entry sharply. | Reconcile the 2023-2027 handoff and distinguish household birth rates, total births, household heads, and resident persons. |
| **5. Identification and tenure evidence** | The reported young fit is 31.0984% versus 34.1166%; aligned-age diagnostics fit less well. The bequest direction is weak. | Choose age/geography/time definitions and confront the large old-age ownership miss before further calibration. |

**Reading route.** Start here and with the completed overnight evidence in Sections 16-18, then read Sections 1 and 13. Use Sections 3-8 for implementation, targets, identification, and demographic details; Sections 9-12 contain policy interpretation and validation. Section 14 states the proposed research direction; Section 15 preserves the earlier reconciliation.

\clearpage

**Meaning of the judgments.** A verified implementation fact is different from an empirical assumption, a passed numerical gate, or an unperformed verification. The overnight audit now shows that reporting floors affect negligible occupied mass and that rare non-global saving choices exist. Their relevance for the tiny credit-fertility effect is assessed separately in Section 16; a local value gain is not a bound on an equilibrium policy effect.

## Status at document preparation

The production reference remains September 4 **task_010**, with loss **30.4829667**, eleven free parameters, twelve production moments, supply elasticity fixed at **0.63**, and tenure-choice dispersion fixed at **0.005**. Two exact repetitions reproduce its parameter and target-fit records.

The separate young-ownership profile adds a thirteenth diagnostic row with a synthetic scale equal to 5% of its target. It is default-off and has not been promoted. The current reported young-ownership gap is **3.0182 percentage points**; the older near-zero ownership failure is historical.

**Updated work snapshot, 5 September 2026.** All 23 cases in Torch panel 16963855 and collector 16964275 have completed. Ridge-plan job 16964564 rejected a flow-versus-stock measurement comparison. Independent read-only receipts explain that rejection and reconstruct both Jacobians; both fail the declared relative-rank threshold. The panel perturbs normalized transformed coordinates by 0.02, not parameters by 2%. The H128 baseline remains unaccepted. No calibration proposal has been launched or promoted.

The seven-case smaller-step follow-up, nine-case saving comparison, four-case nested-grid check, and three-case rental-cap sensitivity are complete. The numerical checks preserve the initial supply-versus-credit contrast. Allowing eight-room rentals leaves the supply birth response near 1.29% but reverses its ownership response. This economic variant is uncalibrated. All planned work is collected; remaining target, identification, and model choices await author review.

**Scope of authoring.** The author-owned manuscript under latex/JMP_DS_draft remains read-only. This is an audit and discussion document, not a revised manuscript. The overnight reconstruction now supplies the standard seventeen-figure diagnostic packet for the actual selected 2023 state; supplemental validation figures are labeled separately.

\clearpage

\tableofcontents

\clearpage

# Full independent audit

## 1. Overall assessment

**I would retain this model, but I would not yet sign off on its quantitative policy magnitudes for circulation.** It contains a coherent and useful mechanism: children increase the value of space; renting and owning provide different housing opportunities; liquid wealth and financing constraints affect access to those opportunities. The September 4 results demonstrate that mechanism at the maintained discretization.

My sign-off would be narrower:

| Object | Judgment |
|---|---|
| Household implementation | Substantial central machinery survives inspection, including sequential births, independent child maturation, transaction accounting, and the financed-share convention. The overnight checks find rare local saving misses and negligible reporting-floor excess; their impact on the tested initial policy effects is negligible. The 120-to-239-node conditional impact check preserves the supply-versus-credit contrast; wealth-support and full-history sensitivity remain unverified. |
| Calibration | **A reproducible local calibration, not a convincingly identified structural estimate.** The score is **30.4829667**. Important housing moments remain wrong, one family-ownership analogue is materially mismatched, and numerical derivatives are unstable. |
| Historical transition | A conditional simulation with imposed household totals and age distributions. It is **not an explanation of those demographic paths**, and it misses the housing-price collapse and subsequent recovery pattern. |
| Finite policy experiments | Useful conditional mechanism comparisons. Their accounting and market gates pass, but demographic closure, external restrictions, and numerical uncertainty limit their interpretation. |
| Perfect-foresight comparison | **Not accepted.** The H128 baseline fails market and fiscal gates. A converged reform cannot rescue an unconverged comparison. |
| Job-market-paper suitability | Worth developing and potentially worth presenting. **Material work is required before circulation**, rather than simply adding limitations to the existing tables. |

Three decisions matter most:

1. **Choose the contribution:** a national housing-access-and-fertility model, or a quantitative account of spatial reallocation and tax lock-in. The current implementation supports the former.
2. **Choose the demographic claim:** conditional household responses from the fitted 2023 state, or resident-population dynamics. The existing birth-to-household normalization does not support the latter.
3. **Make identification and numerical validation precede further searches.** Align age-bounded and family-ownership targets, use the completed numerical checks, and validate the policy-relevant fertility response before trying to reduce the score again.

The initial September 4 review used saved records and bounded arithmetic. The subsequent authorized overnight work reconstructed the dated state and completed separate cluster experiments, without changing production. This document incorporates these stages; Sections 16-18 record exactly what was run and what remains outside the completed checks.

## 2. Source-of-truth map

| Classification | Object | Governing evidence |
|---|---|---|
| **Active and paper-relevant** | One-market lifecycle model with sequential fertility, independent children-at-home counts, permanent and persistent income heterogeneity, housing floors, and warm-glow bequests | [Active package](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2344) and [September 4 status](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md:61) |
| **Active calibration** | September 4 `task_010`, selected from the completed 23-point coordinate panel | [Selected record](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010/summary.json:2) and [complete fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/report/best_target_fit.csv:1) |
| **Active, conditional mechanism evidence** | Five-case 2023-2063 temporary-equilibrium policy packet | [Policy comparison](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_post2023_policy_mechanisms_eta063_kappa005_task010_production_20260904a/report/summary.json:1) |
| **Active but diagnostic** | Equal-rebate impact calculation and eight-cell tax/price/transfer decomposition | [Rebate packet](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_post2023_rebated_tax_impact_eta063_kappa005_task010_20260904a/summary.json:1), [Shapley packet](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_rebated_tax_shapley_diagnosis_eta063_kappa005_task010_20260904a/summary.json:1) |
| **Experimental and unpromoted** | Person-cohort law, household-head bridge, associated terminal roots, and perfect-foresight policy paths | [Person extension](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_perfect_foresight_person_demography.py:169) |
| **Failed numerical comparison** | H128 baseline-reform pair using a preceding calibration | [Explicitly unaccepted comparison](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_pf_person_policy_comparison_20260903b_eta063_kappa005_h128_runtime_diagnostic/summary.json:2) |
| **Historical or superseded** | M4/M5, one-shot fertility, earlier E-series calibrations, August historical-validation figures, September 3 calibration tables | Their findings require revalidation against the September 4 contract. Their losses and policy magnitudes are not interchangeable. |
| **Author manuscript: not integrated** | Quantitative sections and appendix in the protected manuscript | The six specified files contain section/label declarations and `% Author text.`; for example [quantitative model placeholder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_draft/sections/03_quantitative_model.tex:1). |

Four source-of-truth issues require explicit disposition.

**First, parameter labels have been corrected in the live status note.** `chi` is an owner housing-service premium; `H0` is a housing-supply intercept; `theta0` and `theta1` govern bequest utility. The original status table mislabeled these objects as a generic housing weight, household housing scale, and housing-preference intercept and slope. The correction changes documentation only; the equations, estimates, bounds, and production selection are unchanged.

**Second, neither supporting model document describes the complete live system.** The [paper draft](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/intergenerational_housing_fertility_paper_draft.tex:352) still specifies one-shot family-size categories and shared child aging. Its income process, entry equilibrium, calibration, and policy tables also contain superseded objects. The older [dynamic model note](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/dynamic_intergenerational_housing_fertility_model.tex:341) is not a substitute for a current specification.

**Third, the decision ledger is not numerically current despite its September 4 date.** It still contains older losses, rank findings, supply elasticity, and closed-root conclusions. Those remain historical evidence, not descriptions of `task_010`. See [ledger points 5-10](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/ACTIVE_DECISION_LEDGER.md:76).

**Fourth, saved-code provenance and the working tree are not identical.** The original inspection matched **49 of the 50 recorded scientific-file hashes** to then-current files; optional target-profile machinery explained the driver mismatch. The overnight reconstruction instead uses a frozen Git snapshot matching **all fifty production source hashes**, and reproduces all seventy-five historical entries checked. This supports reproduction of the **saved contract**, not blanket reproduction from whatever the current checkout defaults happen to be. Concurrent working-tree changes remain separate from the audited snapshot.

## 3. Quantitative model map

The household is an adult decision unit, interpreted economically as a family with a potential mother. That interpretation is not automatically equivalent to a Census household head or a resident person.

The state is

\[
x=(b,H,i,j,z,n,m),
\]

where \(b\) is liquid wealth, including secured debt; \(H\) is inherited owned housing, with zero denoting renting; \(i\) is location; \(j\) is age; \(z\) combines permanent and persistent earnings components; \(n\) is children ever born, top-coded at \(3+\); and \(m\) counts children currently at home. The active location count is \(i=1\).

| Within-period stage | Economic operation |
|---|---|
| Beginning | Inherit liquid wealth, owned house, income state, parity, and children-at-home count. |
| Fertility | During ages beginning at 18 through 42, choose whether to attempt a birth. Success depends on age-specific fecundity. |
| Successful birth | Move \((n,m)\rightarrow(n+1,m+1)\); pay the first-birth utility cost if \(n=0\). The new child affects current housing needs and policy eligibility. |
| Housing | Choose renting or an owner size, incorporating purchase costs and net sale proceeds. Location choice is degenerate in the active one-location model. |
| Consumption and saving | Choose liquid saving and consumption; renters also choose continuous housing services. |
| End of period | Apply mortality, income transitions, and independent child maturation; preserve completed parity. |
| Next date | Advance survivors and insert entrant households. Historical age reweighting or the chosen post-2023 demographic law determines their masses. |

One period is four years. There are 17 age cells, starting at \(18,22,\ldots,82\), with retirement at 66 and forced terminal exit after the last cell. Survival is certain before retirement; the four subsequent survival probabilities are approximately \(0.9391,0.9185,0.8850,0.8300\).

The active flow utility can be written as

\[
u_t(c,s,m)=
\frac{e(m)^{\sigma-1}
       [c^\alpha s^{1-\alpha}]^{1-\sigma}}
     {1-\sigma}
+\psi_t m,
\qquad
e(m)=\left(\frac{2+0.7m}{2}\right)^{0.7}.
\]

Usable housing is rented space net of the child floor, or owned space net of that floor multiplied by the owner premium. At the selected point,

\[
\bar h(m)=0.364931\,\mathbf 1\{m\ge1\}+0.282210\,m.
\]

Thus the first child adds approximately **0.647 rooms to the floor**, while each additional modeled child adds approximately **0.282 rooms**. These are minimum-needs parameters, not predictions for realized room changes.

Renters choose up to six rooms. Owners choose among \(2,4,6,8,10\) rooms. The mortgage block represents collateralized negative liquid wealth, not a separate amortizing loan.

The principal classifications are:

| Object | Classification |
|---|---|
| Eleven parameters in Section 7 | Estimated jointly |
| Consumption share, income process, entry wealth, survival, financing conventions, housing menu, equivalence scale | Externally measured or maintained restrictions, with differing evidentiary strength |
| Supply elasticity \(0.63\), tenure dispersion \(0.005\) | Fixed restrictions; not estimates from this calibration |
| Old fertility-preference intercept | Normalized to old completed fertility \(2.1\) |
| 2007-2023 preference change | Estimated; linear calendar path imposed |
| Historical household totals and head-age shares | Empirically imposed |
| Initial housing-supply normalization | Accounting normalization |
| Births divided by \(2.1\) | External replacement normalization |
| Post-2023 zero outside entry and full retention | Maintained closure |
| Person-cohort extension and H128 | Experimental |

The resulting economy clears a physical housing market with externally supplied wages and the safe return. It does not contain a complete domestic capital market, construction resource account, landlord balance sheet, or inheritance-recipient equilibrium.

## 4. Equation-to-code audit

“Correct” below refers to the specified object, not to its empirical adequacy.

| Economic object | Written specification | Implementation | Verdict | Evidence |
|---|---|---|---|---|
| Consumption and housing utility | Supporting draft uses Stone-Geary consumption and housing needs | Active equivalence scale, zero consumption subsistence intercept, child housing floor, fixed consumption share | **Inconsistent across sources; coherent active formula** | [Preference construction](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2344), [active floor profile](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/e5f_floor_profile.py:31) |
| Renter housing choice | Continuous rented space with upper bound | Analytical consumption/housing allocation conditional on saving, with renter cap | **Correct conditional allocation; reporting-floor qualification below** | [Renter kernel](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/kernels.py:760) |
| Owner services | Premium on usable owned space | \(\chi[H-\bar h(m)]\); physical demand remains \(H\) | **Correct as implemented** | [Owner kernel](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/kernels.py:957) |
| Purchase timing | Draft pays bond return before housing transactions | Code evaluates the branch at post-transaction wealth and then applies the gross bond return | **Inconsistent across sources**: \(Rb+T_H\) differs from \(R(b+T_H)\) | [Draft budget timing](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/intergenerational_housing_fertility_paper_draft.tex:402), [branch resources](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2509) |
| Down payment | Financed share \(\phi\); cash requirement \((1-\phi)PH\) | Purchases screened using that financed-share convention, including sale proceeds for movers | **Correct as implemented** | [Tenure transitions](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2647) |
| Debt constraint | Draft summarizes \(b'\ge-\phi PH\) for owners | Collateral floor plus an age-tapered rollover allowance for existing debt exceeding that floor | **Correct only under the fuller coded debt rule; draft incomplete** | [Debt helpers](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:110), [kernel floor](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/kernels.py:976) |
| Moving costs | Sale friction and possible location cost | Six-percent sale wedge; location utility cost is inactive with one location | **Correct, with limited scope** | [Household housing choice](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2640) |
| Sequential fertility | Supporting draft still chooses completed family size once | First and continuation attempt logits; successful birth changes both parity and at-home count; first-birth cost paid on success | **Correct active implementation; written specification superseded** | [Fertility block](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2720) |
| Child maturation | Older text uses a shared clock | \(m'\mid m\sim\mathrm{Binomial}(m,7/9)\), conditional on the post-birth count | **Correct as implemented**; memoryless approximation | [Transition matrix](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/parameters.py:902) |
| Family indexing | Parity and dependent children are separate | Seven-dimensional state; flattened family index uses parity fastest | **Correct on the inspected active path** | [State layout](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/state_layout.py:1), [fallback decoding](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2590) |
| Saving optimization | Global Bellman maximization | Golden-section search on the feasible interval; no active previous-policy narrowing | **Rare non-global choices verified; policy relevance assessed in Section 16** | [Search controls](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2476) |
| Bequests | Different documents describe different estate conventions | Normalized warm glow in nonnegative \(b'+PH\), gross of sale costs; \(\theta_n=0\) | **Coherent conditional specification**; not donor-recipient closure | [Bequest function](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:7170) |
| Income | Older draft has five low-risk states | Five persistent states crossed with three permanent levels; explicit transition matrix | **Correct discretization of maintained process**, not a direct estimate of all household earnings risk | [Annual-to-period conversion](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/local_panel.py:1049), [permanent component](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/e6b_profile.py:1) |
| Distribution propagation | Choices followed by saving, aging, mortality, and entry | Pre-fertility, post-fertility, and realized housing distributions are separated | **Verified strength** on the active path | [Forward distribution](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:4411), [current-choice realization](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:3420) |
| Housing market | Draft supply responds to rental user cost | Dated supply is \(K_0(P/P_0)^\eta\), responding to asset price | **Different economic specification under tax changes** | [Dated supply rule](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_dynamic_population_transition.py:105) |
| Family ownership target | Recent parents minus households without resident children | Any dependent-child parent minus never-parent | **Material measurement mismatch** | [Data definition](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/audit_ownership_targets.R:64), [model definition](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:6382) |
| Historical population | Supporting draft has endogenous entry choice | Age-specific household masses reweighted to observations | **Accounting bridge, not entry equilibrium** | [Reweighting](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_open_population_transition.py:567) |
| Birth-to-entry lag | Twenty-year delay | Four-slot queue, whose popped flow enters the next distribution | **Correct indexing for twenty years** | [Queue](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_open_population_transition.py:351), [next-state insertion](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_open_population_transition.py:1205) |
| Fiscal closure | Draft requires taxes, pensions, and rebates to satisfy a government rule | Rebate experiments balance the property-tax account; pensions remain separately prescribed | **Property-tax closure verified; consolidated fiscal closure incomplete** | [Pension rule](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/parameters.py:741), [pension residual](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:5970) |

The transaction-timing discrepancy is not proof that the code lacks an economic interpretation. Beginning-of-period trading can justify applying interest to post-trade wealth. But the manuscript must state that convention, because financing and the opportunity cost of housing depend on it.

## 5. Major findings

| Severity | Finding and claim status | Evidence | Consequence | Required resolution |
|---|---|---|---|---|
| **Fatal before use** | Treating the H128 pair as an accepted policy comparison would be **false** | [Comparison rejection](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_pf_person_policy_comparison_20260903b_eta063_kappa005_h128_runtime_diagnostic/summary.json:2) | Baseline-reform differences are not certified equilibrium effects | Keep the pair diagnostic until both paths satisfy unchanged gates |
| **Fatal before use** | Resident-population, spatial-reallocation, or welfare claims from the maintained packet are **unsupported or absent by construction** | [Population units](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_post2023_policy_mechanisms.py:160), [policy scope](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_post2023_policy_mechanisms.py:1) | These claims would exceed the model’s objects | Narrow claims or change and validate the required architecture |
| **Major** | Family-ownership target is **inconsistent across data and model** | Data and model references in Section 4 | Apparently good fit does not validate the intended parenthood-tenure margin | Align the empirical comparison or implement an event-conditioned model analogue; adopt a new contract |
| **Major** | Historical-to-renewal handoff creates a **demonstrated entry discontinuity** | [Baseline path](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_post2023_policy_mechanisms_eta063_kappa005_task010_production_20260904a/baseline/policy_path.csv:2) | Subsequent composition and price dynamics inherit an artificial household-formation change | Reconcile young-head entry with the historical bridge; quantify closure sensitivity |
| **Major** | Rare non-global saving choices and reporting-floor defects are **verified** | [Golden search](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2476), [post-search floors](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/kernels.py:890) | Repetition does not quantify the policy effect of those defects | Use the occupied-state receipts and conditional global-saving comparison in Section 16 and nested-grid check in Section 17; retain their scope limits |
| **Major** | Local identification is **numerically fragile** | Raw coordinate records and [Jacobian construction](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/analyze_e5f_transition_calibration_panel.py:186) | Central-difference rank overstates what is established about parameter separation | Multiple step sizes, tighter numerical checks, and profiles in weak directions |
| **Major** | Central housing fit failures are **reproducibly demonstrated** | [Full target table](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/report/best_target_fit.csv:1) | Counterfactuals operate through a margin the calibration materially misses | Diagnose measurement, numerical behavior, and mechanism before altering weights |
| **Major** | Household heads, potential mothers, and persons are **not fully reconciled** | [Historical bridge](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_open_population_transition.py:567), [person alignment](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_perfect_foresight_person_demography.py:169) | National fertility and population interpretations require additional assumptions | Explicit unit correspondence and validation against female exposures/headship |
| **Major** | Historical housing-price dynamics fail under the maintained shocks | [Current dated path](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/report/best_transition_path.csv:2) | The dated exercise cannot explain the housing cycle | Present it as conditional calibration, or identify additional historical shocks |
| **Major** | Fiscal closure is **conditional on an externally supported rest of government** | [Fixed pension formula](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/parameters.py:741) | “Fiscally closed” is too broad even for the rebate experiments | State or solve the residual pension/general-government account |
| **Major** | External supply and tenure restrictions lack sufficient policy validation | [Fixed restrictions](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md:66) | Magnitudes depend on objects not estimated by the target system | Source the exact supply aggregation; justify tenure dispersion and assess sensitivity |
| **Closed overnight** | The seventeen standard graphs now exist for the selected dated state | [Dated diagnostic packet](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/diagnostics_dated_units/summary.json) | Shapes and boundary exposure can now be inspected | Keep this packet and the labeled supplemental validation; continue discretization checks |
| **Verified strength** | Exact repeats, target-score reconstruction, and common-state comparisons survive | Selected and repeated records; [collector gates](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/collect_e5f_transition_calibration.py:547) | Results are traceable and reproducible under the saved contract | Preserve these controls |
| **Verified strength** | Sequential birth destinations, independent maturation, and down-payment convention survive | Section 4 references | No evidence that these previously problematic core objects remain broken | Retain the implementation and focused tests |

I did not establish a fatal active-path coding error that invalidates every September 4 result. The distinction is between **verified mismatches**, **unverified numerical properties**, and **claims the model cannot support**.

## 6. Solver and numerical verdict

**Bellman solution.** Finite-lifetime backward induction is appropriate. It does not require an infinite-horizon contraction argument. The active Markov path uses unrestricted feasible saving intervals; an older comment about narrowing around previous policies does not describe the active calls.

One verified issue is global maximization. Golden-section search is reliable for a unimodal objective. Discrete owner choices, borrowing thresholds, and interpolated continuation values can generate nonconcavity. The production code does not establish unimodality or compare the result against a global search over relevant interpolation intervals. A tolerance of \(10^{-3}\) in saving is not a bound on Bellman error or policy-effect error. The separate overnight experiment now makes that comparison and finds negligible changes to the initial supply and credit responses; this does not certify other grids or parameter vectors.

**Reported consumption and housing.** After optimizing, the kernels floor consumption at \(0.04\) in four-year units and excess housing at \(0.01\) rooms. The optimized objective uses much smaller positivity thresholds. If those reporting floors bind, reported consumption and housing need not satisfy the budget associated with the selected saving choice.

This is an algebraic inconsistency conditional on binding floors. The overnight audit now measures it at `task_010`: excess expenditure affects only \(2.91\times10^{-9}\) of household mass, with aggregate excess \(6.19\times10^{-12}\) in model resource units. It is quantitatively negligible in this dated state, while remaining a reporting defect. Section 16 records its scope; existing dead-state and distribution gates are different checks.

**Distribution and timing.** The inspected propagation separates fertility, housing realization, saving, mortality, and child maturation. Current physical housing demand uses realized tenure rather than assigning conditional renter housing to owners. Wealth accounting separates inherited balance sheets from current choices. These are important strengths.

**Interpolation and boundaries.** Wealth interpolation clips outside the grid. The active grid has 120 wealth points over a finite support, and owner housing has five rungs. The overnight extraction finds zero mass at the lower wealth node and 0.781% at the upper node. No occupied adjacent-wealth comparison violates value monotonicity at \(10^{-7}\). The remaining checks include:

- stability under a higher wealth ceiling and a refined historical path; Section 17 completes a conditional 120-to-239-node density check;
- stable ownership at down-payment thresholds;
- stable wealth quantiles;
- robustness to a denser owner menu or renter-cap changes.

A denser housing menu changes both approximation quality and, with independent option tastes, the economic choice set. It needs more care than simply calling it a numerical grid refinement.

**Equilibrium.** The scalar equilibrium searches use brackets and actual excess-demand residuals. The selected dated calibration has maximum relative market residual approximately \(2.85\times10^{-5}\); the five policy paths range from \(3.70\times10^{-5}\) to \(4.78\times10^{-5}\), below the declared \(2\times10^{-4}\) gate. Mass errors are at floating-point scale.

Those results demonstrate clearing **at the maintained discretization**. They do not prove equilibrium uniqueness, global household optimality, or grid-independent policy magnitudes.

**Repeatability.** Both exact repeats reproduce the selected parameter record and target fit. That is strong reproducibility evidence. It is not an independent algorithmic validation.

The smallest decisive numerical package is:

1. Occupied-state consumption/housing budget residuals.
2. A global one-dimensional saving oracle on representative occupied states, emphasizing borrowing and tenure boundaries.
3. Fixed-parameter verification at a denser wealth grid, followed by matched baseline-policy comparisons.
4. Value monotonicity, probability adding-up, boundary mass, and distribution diagnostics from that same solution.

The overnight work completes the dated budget, saving-oracle, value/probability, and boundary-exposure checks. It also generates the standard diagnostic packet. Section 17 also completes matched current baseline-policy comparisons on a denser wealth grid. Reconstructed history and re-estimated target fits on that grid remain unperformed; Section 16 documents the separate saving-method experiment.

## 7. Calibration and identification verdict

### Complete target fit

The score is

\[
L(\theta)=\sum_{k=1}^{12}w_k
\big[m_k(\theta)-m_k^{\mathrm{data}}\big]^2
=30.4829667077.
\]

I independently reconstructed it from the saved rows.

| Moment | Target | Model | Model minus target | Weight | Contribution |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.922907 | +0.004907 | 1,425.739 | 0.034334 |
| Completed childlessness | 0.188000 | 0.189407 | +0.001407 | 17,180.744 | 0.034013 |
| Mean age at first birth | 26.044627 | 26.256032 | +0.211404 | 44.444 | 1.986301 |
| First births at age 30+ | 0.260327 | 0.237579 | -0.022748 | 10,000.000 | 5.174788 |
| First-birth rooms response | 0.720246 | 0.436418 | -0.283828 | 137.565 | 11.082056 |
| Rooms: 3+ minus 1-2 resident children | 0.367700 | 0.403441 | +0.035741 | 2,958.515 | 3.779341 |
| Family ownership gap | 0.167662 | 0.161699 | -0.005963 | 14,229.591 | 0.505992 |
| Ownership, ages 30-55 | 0.575472 | 0.544488 | -0.030984 | 1,207.846 | 1.159555 |
| Mean occupied rooms, ages 18-85 | 5.779970 | 6.317291 | +0.537321 | 11.973 | 3.456812 |
| Aggregate wealth / annual gross labor earnings | 6.873100 | 6.932652 | +0.059552 | 6.288 | 0.022299 |
| Annual bequest flow / aggregate wealth | 0.008800 | 0.008433 | -0.000367 | 5,165,289.256 | 0.697133 |
| Living-old wealth/income p90 / p50 | 3.448111 | 3.236511 | -0.211600 | 56.960 | 2.550344 |
| **Total** | | | | | **30.482967** |

The model produces only **60.6% of the first-birth rooms target**, a miss of **3.33 reported standard errors**. Mean rooms are **9.3% too high**; prime-age ownership is **3.10 percentage points too low**; the late-first-birth share is **2.27 percentage points too low**. These are economically relevant failures, not decorative imperfections.

### Target provenance and measurement

The [complete provenance ledger](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/e5_target_provenance.csv:1) records six measured rows, five provisional rows, and one external normalization. “Measured” does not mean every weight is a sampling variance.

| Rows | Authoritative source, sample, estimator, and uncertainty | Model analogue and assessment |
|---|---|---|
| Completed fertility; childlessness | [CPS builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/cps_fertility/build_cps_fertility_targets.py:1). June 2024, national women 40-44, valid `PTSF1` 0-5, positive supplement weights; \(N=3,278\). Weighted mean and zero share; no fixed effects; 1,000 person bootstrap draws. SEs \(0.026484,0.007629\). | Model age-42 completed-parity distribution. Properly completed fertility, not period TFR. Cohort/age-cell alignment and the survey’s five-child top code remain approximations. |
| First-birth mean age; share 30+ | [NCHS builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/nchs_natality_timing/build_first_birth_timing.R:1). National 1987-2023 files, maternal cohorts 1979-1984, 10,578,871 first births. Count-weighted, boundary-collapsed four-year bins, midpoint labels 20-44; no fixed effects or clustering. | Synthetic cohort hazards reaching age 42 in 2023. Declared scales \(0.15\) years and \(0.01\) are cohort-window stability choices, **not sampling SEs**. The cohort construction is distinct from current-period timing. |
| First-birth rooms response | [PSID estimator](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do:1) and [receipt finalizer](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/finalize_first_birth_rooms_target_receipt.py:1). Women 18+, reference person/spouse, valid dwelling and weights, one woman per household-year; 49,457 household-years, 4,112 women. Sun-Abraham contrast \(+3-(-1)\), person/year effects and age/education controls; person clustering; SE \(0.085260\). | Same model risk set branched into first birth versus childless control, advanced four years using dated policies. An informative conditional comparison, but not the same selection process as empirical treated versus confirmed-childless women. Earlier empirical leads reject a flat prepath; **provisional, not clean causal evidence**. |
| Rooms: 3+ minus 1-2 | [ACS receipt builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/moment_standard_errors/build_active_acs_room_target_receipt.R:1). Pooled 2012-2023, 42 matched metros, heads 30-55 with own children under 18; \(N=918,595\). Household-weighted mean difference, no fixed effects. | Must use resident children, not completed parity. Active scale \(0.018385\) is synthetic. Audited metro-bootstrap SE \(0.069003\) is **not used**. |
| Family ownership gap; ownership rate | [ACS ownership builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/audit_ownership_targets.R:59). Pooled 2012-2023 matched center/periphery metros, heads 30-55, standard structures and positive rooms; ownership sample \(N=1,806,067\). Household-weighted levels/difference, no fixed effects; synthetic scales \(0.008383,0.028774\). | Overall rate is a reasonable conditional analogue. Family gap is materially mismatched: empirical recent parents/no resident children versus model dependent-child parents/never-parents. |
| Mean rooms | [Housing builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/build_intergen_one_market_housing_targets.R:1), checked by the ACS receipt builder. Pooled 2012-2023, 42 metros, heads 18-85; \(N=3,986,589\). Household-weighted literal rooms; no fixed effects. | Realized physical rooms in the model. Synthetic scale \(0.288999\), versus audited metro-bootstrap SE \(0.113293\). Pooled metropolitan target is not a national 2023 level. |
| Wealth / earnings | [PSID wealth builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/audit_aggregate_wealth_earnings_ratio.R:1). Waves 2005-2019; wealth ages 18-85, RP/spouse earnings 18-65; 49,550 family-years, 10,432 people. Ratio of weighted sums, no fixed effects; 999 reference-person bootstrap draws; SE \(0.398836\), rounded to \(0.3988\). | Aggregate net wealth divided by annualized gross labor earnings, excluding pensions from that denominator. Broadly appropriate, conditional on the model earnings and portfolio concepts. |
| Bequests / wealth | Historical literature normalization \(0.0088\); no project microdata builder or measured SE. Synthetic scale \(0.00044\). | Model death-weighted estates, including forced terminal estates, divided by wealth and annualized. **External normalization**, not a freshly measured target. |
| Old wealth/income p90/p50 | [PSID old-wealth builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/audit_intergen_bequest_family_size_targets.R:1). Waves 1984-2019, living reference persons 76-84, observed completed children, income \(>\$1,000\); 4,015 family-years, 1,685 people. Weighted quantile ratio, no fixed effects; 499 reference-person bootstrap draws; SE \(0.132477\), rounded to \(0.1325\). | A living-household statistic, not estates at death. Model retirement income is homogeneous, whereas observed family income is not; consequently the ratio does not isolate the same wealth dispersion in both economies. |

Age windows require explicit interval treatment. Rounding requested ages to four-year state indices does not exactly reproduce samples such as 30-55 or 76-84. This matters especially for retirement and fertility-window boundaries.

### Complete estimated-parameter table

Identification entries describe joint economic discipline, not one-to-one instruments.

| Parameter | Estimate | Bound | Status / near bound | Main moments | Assessment |
|---|---:|---|---|---|---|
| Annual discount factor \(\beta_a\) | 0.995117 | [0.94, 0.9995] | Estimated; collector flag no | Wealth/earnings; bequests; tenure | High patience; jointly substitutable with wealth and housing-price mechanisms |
| First-birth dispersion \(\kappa_1\) | 2.168173 | [0.02, 50] | Estimated; no | Childlessness; first-birth timing; fertility | Joint with first-birth cost and preference change |
| Continuation dispersion \(\kappa_C\) | 1.736471 | [0.02, 50] | Estimated; no | Completed fertility conditional on parenthood; housing/family composition | No direct continuation-hazard target |
| Owner-service premium \(\chi\) | 1.043472 | [0.1, 5] | Estimated; no | Ownership; family gap; rooms | Affects utility services, not physical rooms; derivative unstable |
| Supply intercept \(H_0\) | 14.562959 | [0.2, 80] | Estimated; no | Mean rooms; price-mediated wealth and tenure | Joint demand/supply normalization; not household preference |
| Bequest strength \(\theta_0\) | 0.528428 | [0, 8] | Estimated; no | Bequest flow; wealth/earnings; old wealth | Conditional on mortality, terminal estate treatment, and income process |
| Bequest wealth shifter \(\theta_1\) | 0.107249 | [0.02, 16] | Estimated; **collector flags lower proximity** | Old wealth dispersion; bequest flow | Particularly unstable numerical derivative; no clean separate identification |
| Per-child room floor \(\bar h_n\) | 0.282210 | [0.1, 1.8] | Estimated; no | Family-size rooms gap; first-birth rooms; mean rooms | Strongly entangled with housing supply and wealth |
| First-birth utility cost \(K_1\) | 4.559138 | [0, 8] | Estimated; no | Childlessness and timing | Helps separate first from later births, but is not an observed monetary cost |
| First-child room jump \(\bar h_J\) | 0.364931 | [0, 0.5] | Estimated; no | First-birth rooms versus additional-child rooms | Dominates the weakest central-Jacobian direction |
| Preference change \(\Delta\psi\) | -0.328714 | [-1.5, 0.2] | Estimated; no | Dated cohort fertility and timing | Conditional on linear path, old normalization, and imposed demographics |

The [parameter record](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010/parameter_table.csv:1) also contains:

| Fixed or derived object | Value/restriction | Interpretation |
|---|---|---|
| Old child-preference level | \(0.288017\) | Derived to match old completed fertility \(2.1\) |
| 2023 child-preference level | \(-0.040697\) | Old level plus estimated change |
| Tenure dispersion | \(0.005\) | Fixed; not estimated |
| Dated supply elasticity | \(0.63\) | Fixed; old stationary initialization still uses \(1.75\) |
| Consumption share | \(0.733\) | Rounded CEX estimate \(0.732605\), SE \(0.009801\) |
| Risk aversion | \(2\) | Maintained |
| Child-dependent bequest scale | \(\theta_n=0\) | Child-invariant warm glow |
| Child expenditure-share tilts | Both zero | Replaced by housing-floor specification |
| Equivalence scale | Adult base 2; child weight 0.7; exponent 0.7 | Maintained functional form |
| Baseline financed share | \(0.80\) | Twenty-percent down payment for a new purchase |
| Safe return / depreciation / property tax | Annual 2% / 1.1% / 1% | Period conversion imposed |
| Sale cost | 6% | Fixed |
| Child maturation | \(2/9\) each four-year period | Independent geometric duration |
| Persistent income | \(\rho_a=0.9136\); innovation variance \(0.0426\), scaled by tax progressivity | Literature restriction |
| Permanent log-income variance | Approximately \(0.393053\), SE \(0.033255\) | Project PSID estimate; three discrete permanent levels |
| Entry wealth | PSID young childless-renter distribution | Fixed distribution, with feasibility censoring |
| Housing menu and renter cap | Owners 2,4,6,8,10 rooms; renters at most 6 | Maintained segmentation |

A near-bound flag is not evidence that a constraint binds. In particular, \(\theta_1\) is near the bottom of its **raw** range but is not at the lower endpoint of its log coordinate.

### Weights and identification

There are twelve rows for eleven free parameters. That passes a count check. It does **not** establish eleven informative independent empirical restrictions.

I rebuilt the Jacobian from the saved September 4 coordinate panel. It differentiates standardized residuals with respect to transformed unit coordinates, using radius \(0.02\). It is evaluated at the panel anchor, **not at the selected `task_010`**.

**September 5 update.** The separate 23-case panel centered on `task_010` is now complete. Section 16 reconstructs its original twelve-row and augmented thirteen-row Jacobians: both have algebraic rank eleven and relative rank ten at \(10^{-3}\); their weakest direction is the bequest wealth shifter \(\theta_1\). The table below is retained as a **historical derivative diagnostic at the earlier anchor**, not as the selected-point rank result. The change in the weakest direction across nearby anchors is itself a reason to avoid treating one local matrix as a global identification certificate.

| Derivative construction | Condition number | Rank at relative \(10^{-3}\) |
|---|---:|---:|
| Central differences | 117.75 | 11 |
| Forward differences | 7,706.77 | 10 |
| Backward differences | 99.07 | 11 |

At a relative \(10^{-2}\) threshold, central rank falls to ten. Approximately **95.7% of the squared loading of the weakest central direction** belongs to the first-child room jump.

Forward/backward disagreement is particularly large for:

- \(\theta_1\): 0.840;
- \(\chi\): 0.664;
- annual discounting: 0.402.

These ratios compare the norm of the derivative difference with the sum of the two derivative norms. They demonstrate sensitivity to derivative direction; they do not by themselves distinguish numerical kinks from economic nonlinearities.

The strongest central column association is between the child-room floor and supply intercept, with cosine approximately \(-0.796\). Fertility dispersion, the first-birth cost, and the preference change also have meaningful cross-loadings.

**My identification verdict is therefore weak and numerically fragile local separation, not proven global underidentification and not credible precision.** The first-child jump has not yet earned a strong separate empirical interpretation. Nor do the available rows establish that the two fertility scales independently identify responses to housing costs.

Synthetic weights matter materially. For the additional-child rooms gap, the synthetic scale gives roughly **14 times** the weight implied by the available metro-bootstrap SE. For mean rooms, the synthetic scale gives roughly **one-sixth** of that alternative weight. The choice therefore prioritizes one housing discrepancy over another.

This is a valid explicitly weighted minimum-distance calibration criterion. It is **not an efficient-SMM statistic with a conventional chi-squared interpretation**. No joint covariance matrix, parameter confidence set, or policy uncertainty interval has been established.

The absence of a causal housing-cost-to-fertility target does not by itself prove structural underidentification. It means that the policy derivatives rely on the maintained model and restrictions, and have not been validated against independent causal variation. The local rank and finite-step evidence are separate reasons for caution about parameter precision.

The selected point is the **best valid point found in a local search sequence**. The completed panel, transparent failure handling, and exact repeats support that statement. They do not establish a global optimum. The collector’s exclusion of failed cases through eligibility gates is preferable to assigning convenient artificial losses.

### External restrictions

Several external inputs have real provenance:

- The CEX receipt reports \(\alpha=0.732605\) from 4,928 interviews and 2,452 consumer units, with consumer-unit clustering. It uses childless cash renters aged 30-55, one- or two-person households, expenditure/rent trimming, and a cross-fitted local room-price measure. The [saved estimate](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/cex_child_cost/output/les_2019_2023/les_parameter_estimates.csv:2) supports rounding to \(0.733\). Its transfer from a childless-renter expenditure regression to this household utility system remains a maintained structural mapping.
- The persistent-income numbers match Table IV of [Flodén and Lindé (2001)](https://martinfloden.net/files/Floden%20M%20-%20Linde%20J%20-%20RED%202001.pdf). They are hourly-wage estimates from an older PSID sample, not a current estimate of total family earnings.
- The \(0.181\) tax-progressivity parameter matches [Heathcote, Storesletten, and Violante](https://www.jonathanheathcote.com/hsv_taxation_final.pdf). Scaling log innovations by \(1-\tau\) is correct for their log-linear mapping. Combining that restriction with separately estimated permanent family-income heterogeneity still requires validation of the resulting income distribution.
- The \(0.0088\) bequest ratio appears in Table 2 of [De Nardi and Yang (2014)](https://users.nber.org/~denardim/research/De-Nardi_Yang_EER.pdf), attributed there to Gale and Scholz. That establishes a borrowed historical normalization, not a project-specific contemporary estimate.

I did **not** establish a primary empirical receipt deriving the exact national elasticity \(0.63\). The published abstract of [Baum-Snow and Han](https://www.journals.uchicago.edu/doi/abs/10.1086/728110) reports an average urban floor-space elasticity of \(0.5\), which is not itself the project’s \(0.63\). An aggregation or other source may justify the restriction, but that derivation must be documented.

Likewise, \(0.005\) is presently a maintained tenure-dispersion restriction. Its small size does not make it innocuous: it governs substitution across a discrete housing menu and helps determine the smoothness of equilibrium demand.

## 8. Transition and demographic-closure verdict

### What the dated bridge establishes

The 2007 initial state is not an observed joint distribution of household wealth, tenure, fertility, and income.

The construction:

1. Solves an old stationary household economy.
2. Chooses its fertility-preference level to deliver completed fertility \(2.1\).
3. Preserves conditional economic-state distributions within age cells.
4. Reweights age-cell masses to observed 2007 household-head shares.
5. Normalizes housing supply to clear at the retained old price.
6. Imposes observed household totals and age shares at the subsequent historical dates.

The zero initial housing residual is consequently partly **by construction**. At the selected point, the initial supply stock is scaled by approximately \(1.01044\).

There are also two supply elasticities in this construction: \(1.75\) in the old stationary initialization and \(0.63\) in the dated supply rule. The saved metadata disclose this. Describing the entire construction as simply “the model with elasticity \(0.63\)” conceals that distinction.

The estimated preference change is linear over 2007, 2011, 2015, 2019, and 2023. At each date, households solve as if that date’s preferences and prices remain in place for their remaining lifetime. This is a **temporary-equilibrium expectations rule**, not perfect foresight about the estimated historical path.

The cohort-timing construction appropriately splices age-specific hazards rather than raw counts, preventing imposed changes in age-cell mass from mechanically changing the synthetic cohort’s childlessness. That is a useful measurement repair. It does not turn the imposed household path into an endogenous demographic result.

### What the queue establishes

The twenty-year delay is correctly implemented. For example, a birth recorded in 2023 is appended to the queue; it is released after four subsequent queue advances and enters the next household distribution in **2043**.

The economic conversion is much weaker. Dividing births by \(2.1\) makes replacement fertility reproduce a normalized entrant flow. It is not an estimated mapping from children through survival, partnering, independent residence, and household-head formation.

The maintained transition also uses two different child clocks:

- independent geometric maturation governs household utility and housing needs;
- the deterministic birth-vintage queue governs future household entry.

They can coexist as approximations, but they are not one internally unified childhood-to-household life history.

The handoff is quantitatively important. In the saved baseline:

\[
E_{2023}=0.014337,\qquad E_{2027}=0.063908,
\]

so the youngest household cohort increases by a factor of **4.46** when historical reweighting ends. This is not caused by post-2023 policy births. It comes from the inherited queue and the change in population rule.

Common-state policy differences remain comparable under that rule. Their subsequent composition is nevertheless conditional on this artificial entry pattern.

### Closed and open renewal

For a stationary closed household economy, the relevant renewal condition is that delayed entrants generated by births equal the entrants needed to sustain the stationary age distribution. A normalized completed-fertility number alone does not establish that condition in a reweighted, nonstationary head population.

Historical failures to find a positive closed root establish failure within the investigated parameter/price domain and closure. They do not prove nonexistence everywhere, and their numerical bounds cannot be transferred to `task_010` without checking the parameter contract.

An open closure with fixed outside inflow and retention can be a transparent sensitivity. Its old-state calibration is an accounting identity, not identified migration behavior.

### Person-cohort extension

The extension explicitly tracks annual persons by sex and age, survival, births, migration, and heads. That is a real improvement over adding household mass to children-at-home counts.

However:

- migration is inferred from an external demographic path;
- headship follows an externally supplied, aligned age-sex profile;
- formation and dissolution reconcile head stocks;
- household economic states are reweighted within age cells.

New heads do not arise from an estimated decision about leaving home, partnering, or housing access. Nor are their assets generated by a closed inheritance or household-splitting account.

Its stationary demographic calculation is coherent **conditional on fixed primitives and a fixed births-per-head response**. It is an externally anchored demographic layer, not a fully endogenous population equilibrium.

### Finite horizon versus H128

| Object | Proper interpretation |
|---|---|
| 2007-2023 exercise | Dated household calibration with an imposed demographic bridge |
| 2023-2063 maintained packet | Finite temporary-equilibrium mechanism simulation |
| Old normalized stationary solution | Initialization object |
| Closed-root searches | Existence diagnostics within their searched domain |
| Open-entry cases | Closure sensitivities |
| Person-cohort terminal roots | Conditional stationary endpoints under external demographic primitives |
| H128 paths | Attempted perfect-foresight transitions to those endpoints |

A finite mechanism exercise does **not** require a terminal steady state for its stated finite-date claims.

H128 does require the relevant path and endpoint conditions. Its baseline maximum market residual is \(3.4855\times10^{-4}\), above \(2\times10^{-4}\); its fiscal residual is \(3.6914\times10^{-5}\), above \(2.5\times10^{-5}\). The reform passes both. The [continuation snapshots](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_h128_baseline_stage02_snapshots_20260903/continuation_diagnostic.json:1) show that seven further evaluations under fixed damping did not contract.

The demographic spectral radius \(0.96102\) per four-year period implies a roughly 69.7-year half-life **for the frozen demographic operator**. It does not certify convergence of the coupled nonlinear household, price, and fiscal system.

## 9. Policy and mechanism verdict

### Exact contracts

The five maintained cases share the same reconstructed **pre-fertility 2023 distribution and birth queue**, and reproduce the fitted history. They run through 2063, hold the 2023 preference level fixed, set outside entry to zero and retention to one, and solve temporary equilibria date by date.

| Case | Changed primitive | Eligibility and duration | Fiscal and supply treatment |
|---|---|---|---|
| Baseline | None | Maintained settings throughout | One-percent property tax; revenue discarded; baseline supply curve |
| Supply expansion | Supply intercept multiplied by \(1.2\) | Entire market, permanently from 2023 | No construction subsidy cost or explicit construction resource account |
| Dependent-child credit | Financed share rises from \(0.80\) to \(0.95\) | Households with \(m>0\), including a successful current birth; all eligible owner rungs; applies while eligible | No purchase grant; no default or lender-loss account |
| Combined | Both preceding changes | Same eligibility | Same treatment |
| Unrebated tax | Annual property tax rises from 1% to 2% | All taxed housing, permanently from 2023 | Revenue discarded; supply responds to equilibrium asset price |
| Equal-rebate impact | One-percent or two-percent tax with equal transfer | All modeled adult household units, impact calculation | Jointly solves housing clearing and property-tax revenue = transfer outlays |
| Shapley decomposition | Eight combinations of endpoint tax, asset price, and transfer | Fixed initial distribution | Endpoints are equilibria; hybrid cells generally are not |

For the credit policy, the closing requirement is correct. An eligible buyer of a house worth \(PH=4\) needs \(0.8\) units of liquid wealth under 80% financing and \(0.2\) under 95% financing. Existing excess debt is additionally governed by the rollover rule described above.

The policy is not a first-time-buyer intervention. It has no past-ownership or prior-receipt eligibility state.

### Maintained policy results

These are **percentage changes in births per adult household**, not automatically percentage changes in total births. Ownership changes are percentage points; housing use is physical rooms per adult household.

| Policy | Births/household 2023 | 2043 | 2063 | Ownership 2023 | Ownership 2063 | Rooms/household 2023 | 2063 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Supply +20% | +1.291% | +1.704% | +2.050% | +1.784 pp | +2.671 pp | +5.754% | +8.603% |
| Dependent-child LTV 95% | +0.018% | +0.045% | +0.076% | +0.745 pp | +1.677 pp | +0.044% | +0.020% |
| Combined | +1.314% | +1.766% | +2.155% | +2.963 pp | +4.230 pp | +5.897% | +8.720% |
| Tax 1%→2%, no rebate | -1.100% | -1.296% | -1.623% | -2.301 pp | +0.722 pp | -4.951% | -6.187% |

The [underlying comparison](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_post2023_policy_mechanisms_eta063_kappa005_task010_production_20260904a/report/policy_effects_relative_to_baseline.csv:1) also establishes:

- Supply expansion lowers the impact asset price and user cost by **18.179%**.
- Credit raises them slightly, by **0.060%**.
- The unrebated tax lowers the asset price by **7.746%** but **raises user cost by 14.522%**.

Asset prices and housing costs must therefore remain separate in the tax interpretation.

By 2063, total-birth effects differ from births-per-household effects because household masses differ. Reconstructing total births gives approximately **+2.459%** for supply, **+0.083%** for credit, **+2.575%** combined, and **-1.931%** for the unrebated tax.

The positive 2063 aggregate ownership effect under the tax cannot be interpreted as a positive individual treatment effect. At that date, dependent-child ownership remains slightly below baseline, and the population composition has changed.

### Equal rebate and decomposition

The impact packet contains three distinct equilibria: the unrebated one-percent status quo, rebated one-percent tax, and rebated two-percent tax.

Relative to the **unrebated status quo**, births rise approximately:

- **2.246%** when the existing one-percent tax is rebated;
- **2.756%** under the rebated two-percent tax.

The reform between the two **rebated** equilibria raises births by **0.498%**. Those comparisons must not be conflated.

| Shapley contribution | Births | Ownership | Housing use |
|---|---:|---:|---:|
| Tax rate | -1.323% | -2.471 pp | -6.851% |
| Asset price | +0.291% | +0.671 pp | +1.598% |
| Equal transfer | +1.531% | +0.831 pp | +2.156% |
| **Total** | **+0.498%** | **-0.969 pp** | **-3.097%** |

The decomposition averages marginal contributions across all six orderings and adds exactly. See [implementation](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_rebated_tax_shapley_diagnosis.py:92).

The property-tax budget residuals pass the impact gates: approximately \(7.07\times10^{-6}\) and \(3.30\times10^{-6}\) in absolute value for the rebated one- and two-percent cases. This verifies the **property-tax account**, not all government activity.

The hybrid Shapley cells recompute household choices at assigned tax/price/transfer combinations. They need not clear markets or budgets individually. The exercise is a model incidence decomposition, not three independently identified causal channels.

### Economic interpretation

The small credit effect is **not evidence of a broken fertility transition**. Ownership can change with little change in physical housing or the value of having a child. There is no independent utility payoff from the legal status “owner,” although the owner-service premium means tenure does affect the utility value of a given amount of space.

The credit result is still insufficiently disciplined as a finding about actual fertility policy. Its impact fertility response, **0.018%**, is very small relative to unresolved numerical, parameter, and closure uncertainty. No uncertainty envelope establishes its precision.

The supply and unrebated-tax signs are plausible. Their magnitudes remain predictions of this functional form and calibration. Ratios of birth changes to housing-use changes are policy-specific arc ratios, not structural causal elasticities of fertility with respect to housing services.

Discarded revenue is useful for isolating a gross tax burden. It is a poor benchmark for a realistic fiscal reform unless the omitted expenditure or transfer destination is stated. Conversely, the positive rebated-tax birth response is principally a **transfer result**.

**No verified welfare result exists.** Births, ownership, transfers, and clearing conditions do not rank welfare. A welfare analysis would additionally require a feasible consolidated resource/fiscal account, treatment of affected cohorts and outsiders, and a normative treatment of endogenous population.

## 10. Validation and diagnostics

### Current historical implications

The August historical-validation **figures** are stale. However, their empirical input receipts remain usable: I verified the recorded hashes of the ACS, NCHS, and historical price inputs, then combined those unchanged data with the saved September 4 model rows.

For first-birth timing, I used the existing comparison operator: model age-specific first-birth rates multiplied by national female-age exposure shares. I did not compare raw household-head-weighted birth counts directly with maternal-age data.

| Historical object | Data | Current model | Assessment |
|---|---|---|---|
| Mean first-birth age, 2007→2023 | 25.840→28.090 | 26.913→27.913 | Later endpoint reasonably close; postponement understated |
| Share first births 30+, 2007→2023 | 0.2359→0.3766 | 0.2784→0.3401 | Initially too late, subsequently insufficient postponement |
| Ownership ages 18-85, 2007→2023 | 0.6722→0.6517 | 0.6279→0.6356 | Initial level too low; net change has opposite sign |
| Rooms, 2011→2023 | 5.8683→5.9140 | 6.5513→6.3173 | Too high throughout; change has opposite sign |
| Real house-price index, 2011; 2007=1 | 0.7133 | 1.0403 | Misses the collapse |
| Real house-price index, 2023; 2007=1 | 1.1562 | 1.1371 | Endpoint close despite wrong intervening path |

Evidence: [current dated records](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/report/best_transition_path.csv:1), [current age-specific fertility](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/report/best_dated_period_fertility_by_age.csv:1), [unchanged empirical receipts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_historical_validation_jump11_polish_r2_20260818/classification_and_provenance.json:1), and [comparison operator](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/build_e5f_transition_historical_validation.py:890).

These are genuinely untargeted **time-path implications**, but not wholly independent evidence: they reuse some survey families and outcomes appearing in calibration. Household totals and age shares are not holdouts at all.

The current diagnostic called period TFR falls from approximately \(2.100\) to \(1.429\). Its denominator is model household mass by age. It should not be presented as a directly validated national female TFR without the corresponding unit mapping.

### Other validation

**September 4 reconciliation.** The current young-ownership diagnostic is now available: 31.0984% against an ACS household-head target of 34.1166%, a reported gap of -3.0182 percentage points. This replaces any impression that young ownership is unmeasured or near zero in the selected calibration. The selected dated lifecycle packet is now complete. Section 16 adds annual-age measurement sensitivities and the large older-age ownership discrepancy; geography and reference-period qualifications remain.

| Dimension | Current evidentiary status |
|---|---|
| Lifecycle ownership and housing by tenure/parenthood | Standard dated packet and age-level holdouts available; parenthood-specific validation still needed |
| Parity and continuation fertility | Literal states exist, but completed mean/childlessness and first-birth timing do not validate the full parity or spacing distribution |
| Wealth accumulation/decumulation | Selected aggregate and old-tail targets fit reasonably; realistic lifecycle portfolio paths remain unestablished |
| Bequests | Flow normalization targeted; no verified inheritance-recipient distribution or donor-recipient accounting |
| Mobility | One location prevents spatial mobility validation; within-market sales/downsizing need separate assessment |
| Income distribution | Literature and project inputs exist; the combined discretized earnings distribution needs direct validation |
| Borrowing/renter-cap/tenure thresholds | Occupied-state budgets and grid/menu exposure measured; transaction-interpolation sensitivity remains |
| Value monotonicity and saving optimality | Selected baseline passes monotonicity; global saving removes two local-credit flags and barely changes tested impact aggregates |
| Grid and menu robustness | No current verification sufficient to bound reported policy errors |

### Figures

The September 4 calibration packet contains target-fit and pilot plots. The pilot has a single point and overlapping horizontal-axis labels; it is not an informative search or policy diagnostic.

The five-policy comparison is readable and correctly separates user cost and household mass. It also reveals a sharp common decline and rebound in births per household. That shape needs interpretation alongside age composition and the queue handoff. The near-overlapping population paths conceal small policy differences unless the CSV is consulted.

The Shapley figure is useful for its intended impact decomposition. The separate overnight packet now supplies dated lifecycle, distribution, and equilibrium-residual diagnostics.

A minimum stable packet should contain:

1. Lifecycle ownership, rooms, fertility/parity, liquid wealth, and total wealth, with data where available.
2. Consumption, saving, rented space, tenure probabilities, and fertility probabilities over wealth, with occupied-mass overlays.
3. Separate childless, parent, high-parity, youngest, oldest, poorest, and richest views.
4. Owner-rung shares, renter-cap mass, and borrowing/closing-threshold mass.
5. Asset prices, user costs, physical demand/supply, and market and fiscal residuals.
6. Historical and policy paths with separate birth rates, total births, household heads, and any person counts.

The existing diagnostic code covers much of this. The overnight reconstruction now supplies its seventeen standard figures for the actual selected 2023 state, with dated housing quantities correctly expressed per household. Supplemental age validation and boundary receipts complement that set. The existing dense income legends and modal-housing convention are documented. Standard fertility-by-age curves are childless first-birth propensities, and the income-state completed-children panel uses parity capped at three; neither is the top-bin-adjusted aggregate birth or completed-fertility measure. Conditional renter consumption is also distinct from consumption realized after tenure choice. Read the output README before interpreting these panels as aggregates.

### Resolution of previous audits

After reconstructing the primary system, I read the specified earlier evaluations.

- Older claims that fertility was one-shot, children matured together, income risk was the low-risk numerical construction, or the CPS targets lacked builders are **superseded for the active profile**.
- Concerns about the mortgage shortcut, wealth/tenure validation, mixed weights, grids, and population units survive, with the qualifications above.
- The September 3 audit’s numerical tables are superseded. Its service-versus-tenure distinction broadly survives, but its parameter labels are wrong and its horizon evidence certifies only a frozen demographic block.
- The September 1 hostile audit correctly warns that a logit should be analyzed through response slopes, not merely mass “near indifference.” Its suggestion that an equal rebate is fertility-neutral by construction is **false**: equal transfers can change the parent-versus-childless value difference through different marginal utilities and future choices.
- That audit’s extrapolated fertility ranges and proposed large-effect thresholds are not established bounds for the current model.
- The August recommendation to restore the circulated M5 is a historical specification recommendation, not authority to override the September 4 maintained model.

The earlier PDF is evidence of what was assessed then: [Earlier audit PDF](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/e5f_housing_tenure_fertility_overnight_numerical_audit_20260903_v3.pdf).

## 11. Economic assessment

The interesting mechanism is **access to useful space over the family lifecycle**. Saving and collateral connect housing choices to fertility; children change housing needs; sale costs and bequests affect the release of housing later in life.

Sequential fertility is useful because it distinguishes entering parenthood from having another child and permits timing responses. It also introduces identifying demands. Four-year periods, at most one modeled birth per period, no explicit birth spacing, and memoryless child maturation are consequential approximations.

The separate \(n\) and \(m\) states are productive: completed fertility survives after children leave, while current housing needs disappear. The \(3+\) top code is less satisfactory for policy. The extra \(0.602\) children used in completed-fertility and renewal accounting do not receive their own housing-cost, maturation, or fertility-choice histories. High-parity policy effects inherit that fixed correction.

Tenure is economically useful even without a separate fertility payoff: it changes access to larger homes, collateral, housing services, and transaction costs. But the current menu and renter cap are strong structural restrictions. They need validation through tenure-specific room distributions and transitions, not just an aggregate ownership rate.

Saving is central to the question. Bequests may also matter if the paper emphasizes older households retaining space. Their current specification is less securely justified as an essential two-parameter block. Retention, downsizing, late-life liquid wealth, and estate distribution evidence should determine how much bequest complexity to retain.

The quantitative model connects cleanly to analytical arguments about space needs, financing at purchase, and housing-price feedback. It does not quantify capital-gains taxation, acquisition tax basis, stepped-up basis, cross-jurisdiction sorting, or a planner’s reallocation of titles and pledgeable funds.

The related [Coven-Golder-Gupta-Ndiaye work](https://arpit-gupta-t2ui.squarespace.com/papers2) studies intergenerational housing allocation under financial constraints. Its mechanisms cannot be reduced to spatial relocation alone; capitalization and tenure reallocation also matter. The local primary text confirms additional location, capital-gains, ownership-duration, and financing structure absent here. Their paper is a useful comparison, not a mechanism the current model reproduces in full.

I trust the **conditional qualitative contrast** more than the magnitudes: supply changes move space and fertility more than this credit relaxation does. I do not yet trust the estimated strength of housing-cost effects, the near-zero credit coefficient as a precise policy conclusion, or long-run population magnitudes.

## 12. Minimum viable quantitative package

I would retain:

- one market, explicitly presented as national aggregation;
- sequential first and continuation births;
- distinct completed parity and children-at-home counts;
- independent maturation, disclosed as an approximation;
- income risk and empirically anchored entry wealth;
- liquid saving, owner/renter choice, collateral constraints, and sale costs;
- housing floors and an upward-sloping supply curve;
- the fitted 2023 distribution as a conditional starting point.

I would make **impact and short-horizon household mechanisms** the first circulation target. I would not require a successful H128 transition before presenting that narrower contribution.

For parameters:

- Continue estimating patience, the owner premium, housing-supply scale, broad child housing needs, and fertility preferences jointly with suitable moments.
- Retain both fertility dispersions only if completed parity composition and continuation/timing evidence separate them sufficiently.
- Give the first-child room jump a specific decision: retain it with independent information separating first-child from additional-child needs, or fix it explicitly and report the sensitivity. Do not retain it merely because it is already in the search.
- Keep the bequest wealth shifter free only with stable derivative/profile evidence and suitable old-wealth validation. Otherwise impose a disclosed external restriction and preserve wealth/retention evidence as validation.
- Do not free supply elasticity or tenure dispersion without adding independent information about supply responses and tenure transitions.

The minimum target system should retain the completed-fertility, childlessness, timing, first-birth housing, additional-child housing, tenure, wealth, and bequest information, **after aligning their definitions**. Every proposed removal or reweighting must state which parameter block loses information.

Specific identifying additions or replacements are:

| Weak block | Required information |
|---|---|
| First-child jump versus per-child floor | Consistently measured first-birth response and additional-child room responses |
| Owner premium versus housing access | Tenure-specific rooms, ownership entry/exit, and liquid wealth around purchase |
| First versus continuation fertility | Completed parity shares and coherently defined continuation/timing moments |
| Bequest strength versus wealth shifter | Late-life portfolio composition, decumulation, retention/downsizing, and distributional moments |
| Housing-cost-to-fertility transmission | Credible exogenous cost variation, or an explicit external elasticity validation with matched units |

A causal cost-to-fertility moment is not a mathematical prerequisite for every structural calibration. Here it is particularly valuable because the intended policy response is weakly validated by the existing level and reverse-direction housing moments.

For presentation:

- **Main text:** supply expansion and the contrasting credit experiment, once numerical and identifying checks pass.
- **Appendix:** unrebated tax mechanism, rebated impact comparison, Shapley decomposition, external-parameter and closure sensitivities.
- **Park:** H128 policy magnitudes, resident-population forecasts from the queue model, capital-gains or spatial-policy interpretations, and welfare rankings.

## 13. Priority ladder

| Priority | Next concrete action | Evidence produced | Author decision? |
|---|---|---|---|
| **Fatal before use** | Exclude unaccepted H128 differences and unsupported population/spatial/welfare claims | A claim set consistent with actual objects | Yes, on scope |
| **Completed for the selected state and tested impacts** | Occupied-state budgets, global saving, and matched 120/239-node supply/credit comparisons | The numerical checks preserve the supply-versus-credit contrast; see Sections 16-17 | No production change adopted |
| **Important before circulating** | Align age-bounded and family-ownership data and model definitions | Comparable moments and an explicitly approved target fingerprint | Yes, target-contract choice |
| **Important before circulating** | Reconcile the 2023-2027 entrant jump | Explicit formation/entry accounting and matched closure sensitivities | Yes |
| **Important before circulating** | Validate wealth support, historical discretization, and tenure-specific housing opportunities | Extend the completed conditional grid check and rental-cap sensitivity with empirical discipline | Yes for an economic menu change |
| **Important before circulating** | Extend the completed three-column step check when defining an identified estimation strategy; profile weak directions | Stable local information and policy variation along nearly equivalent fits | Yes on identifying evidence or restrictions |
| **Important before circulating** | Explain the first-birth rooms and mean-room tension before changing weights | Diagnosis separating measurement, numerics, and missing mechanism | Yes if changing model or target weights |
| **Important before circulating** | Source the exact \(0.63\) restriction and justify \(0.005\); evaluate their effects | External provenance and policy sensitivity | Yes on maintained restrictions |
| **Important before circulating** | State and audit pension/general-government financing | Complete meaning of fiscal closure | Yes |
| **Completed diagnostic packet; validation remains active** | Retain the seventeen dated-state figures and supplemental age validation | Inspectable policies, distributions, residuals, and holdouts | No new economic choice for regeneration |
| **Useful improvement** | Validate period fertility with female exposure and parity-specific flows | Evidence about demographic speed beyond completed fertility | Yes on measurement convention |
| **Useful improvement** | Validate young wealth, ownership transitions, renter-cap mass, and old downsizing | Evidence that the access/retention mechanism operates in plausible states | No initially |
| **Useful improvement** | Evaluate policy effects over credible near-fit parameter sets | Magnitude ranges reflecting calibration uncertainty | Yes on admissible-set definition |
| **Optional** | Develop explicit household formation, inheritance transmission, mortgage contracts, or multiple locations | A richer model for broader claims | Yes; substantial scope expansion |
| **Optional** | Resume H128 with a materially different solver strategy | Accepted coupled path, if successful | Yes; not another unchanged continuation |

## 14. Final recommendation

If this were my job-market paper, I would stop broadening the quantitative model for now. I would keep the sequential housing-access mechanism and make the near-term contribution narrower: **which housing interventions change usable space for families, and how does that affect fertility within a disciplined household model?**

I would first resolve the age and family-ownership mappings, the household-entry handoff, and empirical housing-opportunity validation. The completed saving-method and wealth-grid comparisons support retaining the qualitative supply-versus-credit contrast. I would then evaluate policy variation across credible near-fit parameter values. The present score and exact repeats are a sound platform for that work, but they are not the endpoint.

I would lead with supply and credit only after those checks. I would present the rebated-tax result as incidence driven substantially by the transfer rule. I would leave spatial reallocation, capital-gains lock-in, long-run resident population, and welfare outside the quantitative claims until the model actually contains and validates those objects.

**The model is worth keeping. The present quantitative package is not yet ready for circulation as a set of identified policy estimates.**

\clearpage

## 15. Reconciliation with concurrent calibration work

This reconciliation checked the live files, the young-ownership smoke, the completed sensitivity cases, and Torch’s queue. **The handoff largely agrees with my audit. The young-ownership work adds useful evidence and changes that part of the assessment; it does not overturn the broader findings.**

The audit already used `task_010`, loss **30.48297**, its two exact reproductions, and the updated policy results in items 1-8. The complete [twelve-target fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010/target_fit_long.csv) and [parameter estimates, bounds, and bound flags](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010/parameter_table.csv) remain the production reference.

The six amendments below are incorporated into the assessment and supersede older descriptions of young-ownership evidence and pending work.

1. **Young ownership is no longer an unmeasured concern or a near-zero failure.**\
   The saved diagnostic records **31.0984%**, against the ACS target of **34.1166%**, a shortfall of **3.0182 percentage points**. I verified that the smoke preserves the production parameter vector and all twelve original target-fit rows exactly. The historical record also supports the move from a hard target to validation during the July redesign.

   This materially improves the assessment of the current tenure fit. The older near-zero results should appear only as historical findings.

2. **The new comparison is useful, but “cleanly matched” needs qualification.**\
   The ACS target uses household heads in the project’s matched metropolitan housing sample, pooled over **2012-2023**. It is not a national 2023 ownership rate. Moreover, the model statistic includes whole age nodes **26, 30, and 34**, without prorating the four-year cells to the literal ACS interval 25-34. This follows from the [age-index mapping](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:6102) and [ownership calculation](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:6351).

   Therefore, **3.0182 percentage points is the verified reported gap**, with age, geography, and reference-period qualifications. The household-head definition is an improvement; exact measurement alignment remains to be established.

3. **The additional target is a separate diagnostic contract.**\
   I verified the default-off profile, separate fingerprint, unchanged original rows, and synthetic scale equal to 5% of the target. At the unchanged parameters, the additional row contributes **3.13063**, making the thirteen-row diagnostic objective **33.61359**. This is an accounting addition, not deterioration in the original fit. The [complete thirteen-row table](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_youngown_overid_smoke_20260904a/task_001/target_fit_long.csv) documents it.

   The scale is not an empirical standard error, and “overidentifying” does not by itself establish stronger identification or a valid statistical specification test.

4. **The winner-centered sensitivity analysis is now underway.**\
   Torch confirms **14/23 completed**, with cases 15-23 pending on priority. Collector `16964275` waits for the panel; Jacobian/ridge-plan job `16964564` waits for the collector. The latter prepares diagnostics and candidate proposals; its script does not launch another calibration.

   One correction to item 16: these are **±0.02 changes in normalized, transformed coordinates**, not ±2% changes in parameter values. For example, the positive steps raise \(\chi\) about **8.14%** and \(H_0\) about **12.73%**. The partial trade-offs in item 18 are supported, but should not be described as responses to tiny 2% perturbations.

   The review will now mark the selected-point Jacobian as **in progress**. Once complete, it should compare the original twelve rows with the augmented thirteen rows at the same center, including forward/backward derivative agreement.

5. **The policy wording needs two distinctions.**\
   The quoted 2063 effects concern **births per model household**, not aggregate births. Recalculating from the saved paths gives:

   | Policy | Births per household, 2063 | Total births, 2063 |
   |---|---:|---:|
   | Supply +20% | +2.050% | +2.459% |
   | Dependent-child LTV95 | +0.076% | +0.083% |
   | Unrebated tax increase | -1.623% | -1.931% |

   Also, the **+0.498%** rebate result compares **rebated 2% taxation with rebated 1% taxation**. It is not the effect relative to the unrebated status quo.

6. **The mechanism conclusion should remain conditional.**\
   Item 9 accurately summarizes the supply and credit experiments: supply moves housing consumption and fertility, while credit mostly moves ownership. But housing services are not the only fertility channel-the rebate operates through household resources and can increase births even while aggregate housing use falls.

   The parameter labels also need correction: \(\chi\) is the **owner housing-service premium**, \(H_0\) the **housing-supply intercept**, and \(\theta_0,\theta_1\) govern **bequest preferences**.

The remaining audit concerns still stand: the first-birth rooms shortfall, family-ownership target mismatch, numerical verification gaps, demographic handoff, fiscal scope, and fragile identification of policy fertility responses. Passing the existing gates does not resolve those separate questions.

**This reconciliation supports keeping calibration frozen until the full panel and Jacobian are assessed.** The reconciliation was read-only. The present PDF contains both the full audit and these amendments.

\clearpage

## 16. Overnight verification, 4-5 September 2026

The author subsequently authorized up to twelve hours of further verification and exploratory quantification, including use of Torch. This section supersedes the earlier pending-work descriptions. The production `task_010`, its twelve target rows, weights, eleven free parameters, and external restrictions remain unchanged. The protected manuscript remains read-only. The [morning discussion brief](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_overnight_morning_review.md) contains the compact readout and complete production and unpromoted-panel tables.

### What is now established numerically

The independent reconstruction uses a separate source snapshot that matches all **fifty source hashes** recorded for the production calibration. It reproduces **all seventy-five checked historical entries exactly**, covering fifteen fields at each date from 2007 to 2023. Torch job `16981272` completes this exact-loop smoke in 6 minutes 43 seconds, using about 1.2 GB. A hash-pinned 2023 state checkpoint now permits inspection without another historical reconstruction. The standard seventeen-figure diagnostic packet has been generated from this dated state, rather than from a substituted stationary solution.

The occupied-state budget audit finds that the post-optimization consumption/housing reporting floors create excess expenditure in only **\(2.91\times10^{-9}\) of household mass**. The mass-weighted excess is **\(6.19\times10^{-12}\)** in model resource units. A maximum occupied-cell discrepancy of 0.0245 occurs in a negligible tail. This resolves the earlier concern in a qualified way: a reporting-floor defect exists, but it is not quantitatively important in this particular 2023 solution. It does not certify all parameter vectors, policies, or grid resolutions.

The independent saving oracle maximizes the same objective interval by interval under the saved piecewise-linear continuation value. Within an interpolation interval, utility is concave and the continuation value is affine; endpoints, rental-cap kinks, and feasible first-order roots exhaust the candidates. Independent bounded smooth optimizers verify the interval calculations. The full audit, Torch job `16981404`, uses **50,000 household-weighted draws**, a fixed seed, and additional separately weighted boundary states. It inspects **30,541 distinct states**. Nine sample draws, **0.018%**, have a value gain above \(10^{-3}\) and a saving change exceeding 0.01. The largest saving change is 1.123 model units, and the largest value gain is 0.01089. These are rare non-global choices, not a general failure of saving monotonicity or a measured welfare loss.

No occupied adjacent-wealth comparison exhibits a value decline larger than \(10^{-7}\). All saved fertility, continuation-birth, tenure, and location probabilities are finite and lie in \([0,1]\). The [numerical receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/numerical_full/summary.json), state tables, and unchanged standard figures support these statements.

The boundary receipts make the remaining grid question concrete. The wealth grid has 120 nodes over \([-12,30]\): 0.781% of current household mass is at its top node, compared with 0.0351% within fertile ages. Original saving policies exceed the upper node for 0.142% of household mass, but the largest occupied excess is only 0.00128 model units and its weighted aggregate is \(4.22\times10^{-7}\). Housing-menu exposure is more common: **13.22% of renters** occupy the six-room cap, rising to **20.40% among renters with dependent children**; **33.49% of owners** choose the largest offered home, ten rooms. Choosing the largest discrete option does not prove that a larger home would be chosen if offered. It does make housing-menu sensitivity an economically relevant check. The [boundary receipts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/boundary_exposure_summary.json) record the exact distributions and thresholds.

### Complete global-saving policy comparison

Torch job **16985060** completed the nine-case comparison in **14 minutes 6 seconds**, after corrected exact-loop smoke **16984897** passed. The global method maximizes over all saving interpolation segments at each backward step. Parameters, reported-policy floors, current resource equations, the inherited 2023 pre-choice distribution, and historical population remain fixed. Six cases separately clear the current housing market: two saving methods times baseline, supply +20%, and dependent-child LTV 95%. Three further global-method cases use the corresponding local-method policy and baseline prices.

The requested market residual is **\(10^{-5}\)**, tighter than the original experiment; the maximum achieved residual is **\(7.05\times10^{-6}\)**. The original acceptance gate remains \(2\times10^{-4}\). Each policy change below uses its own method's baseline. Ownership covers all adult households; births and rooms are measured per household.

| Impact policy | Saving method | Births, % | Rooms, % | Ownership, pp |
|---|---|---:|---:|---:|
| Supply +20% | Production local search | 1.290982 | 5.753224 | 1.783382 |
| Supply +20% | Global saving | 1.290873 | 5.753410 | 1.783440 |
| Dependent-child LTV 95% | Production local search | 0.017845 | 0.041904 | 0.744055 |
| Dependent-child LTV 95% | Global saving | 0.017850 | 0.041861 | 0.744088 |

Global saving changes the supply birth response by **-0.000109 percentage points**, and the credit response by **+0.0000046 percentage points**. The latter is **0.026% of the small credit effect**. Tightening market clearing moves the original credit estimate from approximately 0.01836% to 0.01785%, a larger change than replacing the optimizer. These results support the qualitative contrast while discouraging excessive precision for the tiny credit coefficient. Separate clearing reaches identical recorded prices under the two methods; the additional fixed-price comparisons therefore coincide at this numerical resolution. This is not proof that the exact equilibrium price correction is zero.

![Supplemental comparison of saving methods. The panels use different vertical scales. Both methods use the same calibrated parameters and inherited 2023 population.](output/model/e5f_overnight_independent_verification_20260905a/global_saving_receipts/supplemental_global_saving.pdf){width=90%}

All nine cases retain valid choice probabilities and negligible reported-budget excess: the maximum affected mass share remains \(2.91\times10^{-9}\), and the maximum weighted resource gap is \(6.31\times10^{-12}\). The local credit case has two adjacent-wealth value declines, at indebted renter states with three dependent children and ages 30 and 38. Their lower-node mass is \(1.08\times10^{-7}\) of households and the maximum decline is 0.00805. Global saving removes both; all other cases have no occupied decline at the \(10^{-7}\) screen. At the original common baseline price, global saving weakly improves the value function across **2,692,140 states above the production value cutoff**, with no decline exceeding \(10^{-7}\). Value comparisons across different prices are not dominance tests. None of these utility differences is a welfare calculation.

The housing-cap exposure rises under supply expansion: approximately **26.7% of renters**, and **40.1% of renters with dependent children**, reach the six-room rental maximum. This gives housing-menu sensitivity a concrete priority. Adding owner options under option-specific tastes changes the choice set and its inclusive value; it should not be described as pure numerical grid refinement.

**A diagnostic-adapter error was caught and repaired before interpretation.** The first checkpoint-based smoke, job 16983266, failed the exact local reproduction gate because the new audit driver had not reinstalled three sequential calendar operators after loading saved objects. Production reconstruction used the correct initialization and was unaffected. A zero-solve replay restored exact equality of all saved distributions and aggregate quantities. The corrected driver asserts the operator identities before every case, and its replacement smoke passed before the full panel was submitted. The failed run is retained as evidence, not counted as a successful comparison. No production kernel or numerical gate was changed.

The [complete numerical readout](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/global_saving_receipts/summary.json) verifies the nine case checkpoint hashes, operator contract, requested clearing precision, and seventeen standard figures per case. Its `all_case_quantities.csv` and `policy_effects.csv` record every comparison. These are **conditional impact experiments**, not a new calibration, a revised 2063 transition, an H128 solution, or confidence intervals. The common historical distribution was deliberately held fixed; re-solving the historical path with global saving is outside this experiment.

### The completed panel and the failed follow-up

All **23/23** young-ownership cases and their collector completed successfully. The follow-up ridge planner rejected a comparison between a birth-flow ratio and reconstructed completed fertility. Across all cases, the exact queue, outside-entry, retention, and renewal identities close within **\(2.22\times10^{-16}\)**. The flow-versus-stock measurement difference is at most **\(2.62\times10^{-9}\)**, inside the calibration driver's declared **\(2\times10^{-8}\)** measurement gate. The ridge planner instead applies **\(2\times10^{-10}\)** to that comparison, rejecting eighteen cases. The planner and production gates were left unchanged; the [independent receipts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/panel_receipts/summary.json) preserve both the original passes and the stricter planner rejections.

The independently reconstructed weighted Jacobian has algebraic rank eleven in both target systems, but relative rank ten at the declared \(10^{-3}\) threshold. The smallest singular value changes only from **0.560326** with twelve rows to **0.560335** with thirteen. The weakest direction places **99.08%** of its squared loading on the bequest wealth shifter \(\theta_1\). The extra ownership row mostly strengthens sensitivity to the owner-service premium \(\chi\). A larger condition number after adding that row does not mean information was destroyed: it means the already strongest direction became stronger while the weakest barely moved.

This is **weak local identification**, not proof of structural underidentification. Deleting the \(\theta_1\) column leaves relative rank ten out of ten at the \(10^{-3}\) threshold (condition number 118), but does not supply an economic justification for fixing that parameter. An external restriction or new identifying evidence remains an author decision. No ridge proposal was submitted.

### What the smaller-step panel establishes

The new seven-case subset completed in Torch smoke **16982612** and array **16985059**. It uses the original thirteen-row diagnostic profile, seed 2026090505, and exactly the anchor plus/minus three coordinates: annual patience, owner premium, and bequest wealth shifter. All source, target, historical, market, accounting, feasibility, and population contracts pass. The anchor reproduces the original point. Each step is **0.005 in the normalized transformed coordinate**, one quarter of the original 0.02; this is three checked columns, not a new eleven-column Jacobian or a rank test.

For each column, the cosine measures agreement between the forward and backward weighted derivatives. The final column measures the change in the entire central derivative vector relative to the norm of its larger-step estimate.

| Parameter | Side cosine, 0.02 | Side cosine, 0.005 | Central-vector change |
|---|---:|---:|---:|
| Annual patience \(\beta\) | 0.0849 | 0.9590 | 35.2% |
| Owner premium \(\chi\) | 0.4714 | 0.9493 | 52.8% |
| Bequest shifter \(\theta_1\) | 0.9717 | 0.8767 | 17.6% |

The first two columns become much more internally consistent at the smaller step, but their central estimates remain step dependent. The bequest-shifter column is small at both radii, with weighted norms 2.495 and 2.482; a small column alone is not the same as a weak conditional singular direction. The full larger-step matrix and the smaller-step column check together motivate caution about local numerical identification.

An approximately 2% increase in the owner premium raises reported young ownership from **31.10% to 38.97%**, while reducing the first-birth rooms response from **0.4364 to 0.4041**. Thus the housing-response trade-off survives the narrower perturbation. Greater patience has a much smaller favorable effect on first-birth rooms, and the bequest shifter barely moves that response in these one-coordinate comparisons. These observations do not prove that a jointly varied parameter vector cannot fit better.

![Supplemental step-size comparison. Dashed lines are the corresponding target values; connected points are actual evaluations at both radii, with no missing derivative columns imputed.](output/model/e5f_overnight_independent_verification_20260905a/smallstep_receipts/supplemental_step_stability.pdf){width=95%}

The [seven-case receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/smallstep_receipts/summary.json) pins every task summary and records the full step comparison. Its `all_target_fits.csv` gives every target, model value, gap, weight, and contribution; `all_parameters.csv` gives every estimate, bound, and restriction. The complete corresponding 23-case tables remain in `panel_receipts/`. The best inherited twelve-row score, **29.3859** at original panel task 003, is unpromoted and fully tabulated in the morning brief. No additional optimization was launched.

### Measurement changes the tenure assessment

The production young-ownership statistic is reproduced at **31.0984%**. But the model selects entire nodes 26, 30, and 34, whereas the ACS estimator averages literal ages 25-34. Holding the model solution fixed and weighting both sides by the same annual-age ACS household weights yields:

| Window and diagnostic model measurement | Model | ACS | Gap |
|---|---:|---:|---:|
| Ages 25-34, constant rates over \([a_j,a_j+4)\) | 26.9792% | 34.1166% | -7.1374 pp |
| Ages 25-34, linear interpolation between nodes | 30.2169% | 34.1166% | -3.8997 pp |
| Ages 80-84, constant rates over model cells | 99.2320% | 74.9713% | +24.2607 pp |

These are diagnostic discretization choices, not newly approved targets. The model does not identify annual within-cell choices. The original three-percentage-point young-ownership comparison remains a correct reproduction of its reported operator, but is not a fully aligned age-window estimate. The old-age ownership miss is large under either interpolation convention. The ACS sample remains metropolitan household heads, pooled over 2012-2023; reference period and geography qualifications still apply. Exact annual-age weights, estimates, and their model mappings appear in the [measurement receipts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/panel_receipts/ownership_annual_age_receipts.csv).

This strengthens the case for direct late-life tenure, portfolio, and housing-retention validation. A rising cross-sectional age profile alone does not prove that individual older households never downsize; the current result is a failure to match the age-specific ownership level.

### Entry accounting and priorities

The saved paths verify that current age-18 mass rises from **0.0143371 in 2023 to 0.0639084 in 2027**, a factor of **4.4575**. The identity is exactly \(E_{2027}=M_{2023}+\rho_{2023}B_{2023}=0+1\times0.0639084\). The due queue contains historical births, while the 2023 entrant cohort comes from the imposed historical age distribution. This is a historical-to-forward handoff, not a 2023 policy fertility effect. Policy-induced births affect entry only twenty years later under the present contract.

The priority ordering is therefore more precise than in the initial review: establish one ownership measurement convention; confront the older-age tenure miss and weak bequest block; then use the additional grid and rental-cap evidence in Sections 17-18 to specify the remaining transaction, wealth-support, and housing-opportunity checks. The saving-method comparison is complete and has negligible effects on the two initial policy responses tested. Smaller steps clarify local derivatives without supplying a new identification certificate. Keep the first-birth rooms miss central. Keep household formation, resident-population forecasts, and the fiscal meaning of alternative rebates as distinct specification decisions. The conditional supply-versus-credit interpretation survives this numerical test; its empirical identification and long-horizon interpretation remain separate questions.

\clearpage

## 17. Age alignment and nested wealth grids

These additional checks keep the production selection, target system, weights,
and economic parameters fixed. The continuation design is recorded in the
overnight run plan; it separates saved-array measurement checks from a bounded
nested-wealth-grid experiment.

### The young-ownership discrepancy is mainly the age window

The national ACS age builder assigns each annual age to a four-year cell using
integer floor division. The historical bridge uses closed integer bins
\([18+4j,21+4j]\), preserving the distribution of all other states conditional on
age. Thus nodes 26, 30, and 34 receive age-bin masses for **26-29, 30-33, and
34-37**. In contrast, `age_to_index` rounds the requested age to the nearest
label. Applying that operator to the bounds 25 and 34 selects these three whole
nodes. The apparent 25-34 ownership comparison therefore combines different age
windows under the historical bridge's implemented convention.

The table isolates window selection and weighting in the unchanged 2023 state.
The original model rate is reconstructed from the saved lifecycle masses and
rates, not copied from a hard-coded scalar. ACS rows retain the same metropolitan
household-head sample, housing restrictions, and pooled 2012-2023 period.

| Diagnostic measurement | Model | ACS | Gap, pp |
|---|---:|---:|---:|
| Stored comparison: three model nodes versus ACS 25-34 | 31.0984% | 34.1166% | -3.0182 |
| ACS 26-37, original model age weights | 31.0984% | 39.5101% | -8.4117 |
| ACS 26-37, common ACS age weights on both sides | 31.1798% | 39.5101% | -8.3303 |
| ACS 25-34, common weights and constant within-cell rates | 26.9792% | 34.1166% | -7.1374 |
| ACS 25-34, common weights and linear rates between labels | 30.2169% | 34.1166% | -3.8997 |

In an ordered arithmetic decomposition, changing the ACS window from 25-34 to
26-37 contributes **-5.39346 pp** to the original gap. Standardizing the model to
ACS age weights contributes **+0.08141 pp**. Together these move the gap from
-3.01822 to -8.33027 pp. This establishes that the particular discrepancy is
mainly about the selected age window, not differing age weights. It does not
establish a new preferred target: the model has no identified annual within-cell
ownership path, and geography and period alignment remain separate choices.

The review should extend to the existing age-bounded targets. The active prime
ownership statistic selects nodes 30 through 54, which the bridge treats as full
cells covering **30-57**, while its ACS target uses literal ages 30-55. The saved
rate reproduces at **54.4488%**, versus **57.5472%** in ACS. Applying common ACS
weights to literal ages 30-55 and holding rates constant within model cells gives
**53.4737%**, a **-4.0735 pp** gap. Matching the full 30-57 cells with common ACS
weights instead gives model **54.9512%**, ACS **58.5175%**, and a **-3.5663 pp**
gap. These remain diagnostic comparisons; the complete production fit and its
loss contribution stay unchanged. The living-old wealth/income dispersion
operator uses a direct label test, \(76\leq a_j\leq84\), selecting nodes 78 and
82 (bridge cells 78-85), rather than the requested literal 76-84 interval.
Consequently, bequest-block identification should also be assessed after an
explicit age-measurement contract is agreed. This source observation does not
establish how much that wealth-quantile target would change.

Fertility measurement uses a different reporting convention. Timing statistics
explicitly label flows at cell midpoints, \(20,24,\ldots,44\). Completed-fertility
stocks select node 42 through `age_to_index`, although their docstring calls it
a cohort centered at 42. The executable age-bin map and the written descriptions
should be reconciled before changing empirical target rows. No model operator,
target builder, or production statistic was modified in this check.

Source evidence is in
[the national age builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/build_national_householder_age_path.py:154),
[the historical age bridge](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_open_population_transition.py:567),
[the age-index operator](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:6102),
and [fertility measurement](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_transition_calibration.py:1256).
The [continuation receipts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/continuation_receipts/summary.json)
pin the original lifecycle table, annual-age weights, and decomposition.

### The visible first-birth spike has no current household exposure

The unchanged standard renter-policy plots show a narrow first-birth attempt
probability spike under local saving at age 42, wealth **-4.72093**, and income
state 7. Replacing local saving by the global interval method changes its attempt
probability from **1 to 0**. The inherited pre-choice mass at that exact state is
**zero**, for baseline, supply expansion, and family credit. Its contribution to
the current first-birth difference is therefore exactly zero.

A broader screen considers fertile childless states above the production value
cutoff where the two methods' attempt probabilities differ by more than 0.1. There are **19 baseline,
20 supply, and 19 credit grid states** in that screen; their total inherited mass
is also exactly zero. These are counts of numerical state cells, not sampled
households or estimated frequencies. The comparison uses identical saved prices
and pre-choice measures for each local/global policy pair. It weights probability
differences by age-specific fecundity and childless mass, and independently
checks the implied first-birth difference against changes in post-fertility
childless stocks.

That independent stock replay agrees within **\(6.36\times10^{-17}\)**. Across
occupied fertile childless states (mass above \(10^{-12}\)), the largest attempt
probability differences are **0.000642** at baseline, **0.000870** under supply,
and **0.000643** under family credit. These are probability units, not percentage
changes in aggregate births.

This resolves why those particular plots can change visibly while aggregate
births barely move. It does not certify off-distribution behavior, the historical
path under another method, or other parameter vectors. The separate two local
credit value declines occur at occupied states with three dependent children;
they are a different finding and retain their small positive exposure reported
in Section 16. The standard plots have been retained intact.

\clearpage

### Nested wealth-grid check

The finer grid adds each interval midpoint, raising the node count from
**120 to 239**. Every original node and its inherited 2023 pre-choice mass stay
in place; new nodes initially have zero mass. Both grids use global saving,
unchanged parameters and housing choices, wealth endpoints **-12 and 30**, and
separately cleared current markets at a requested residual of \(10^{-5}\).

Smoke **16986180** exactly reproduces the 120-node global baseline and solves
the full 239-node baseline market in 16 minutes. Job **16986566** completes both
fine-grid policies in about 29 minutes. All three fine markets pass: maximum
residual **\(2.29\times10^{-6}\)**, valid population/probabilities/feasibility,
and no occupied wealth-value declines above \(10^{-7}\). Reporting-floor excess
remains negligible and separately measured. All four checkpoints and their
seventeen-graph packets are verified.

Each effect below is measured against its own grid's baseline. Ownership covers
all households; births and rooms are per household.

| Policy | Wealth nodes | Births, % | Rooms, % | Ownership, pp |
|---|---:|---:|---:|---:|
| Supply +20% | 120 | 1.290873 | 5.753410 | 1.783440 |
| Supply +20% | 239 | 1.288387 | 5.759105 | 1.769706 |
| Dependent-child LTV 95% | 120 | 0.017850 | 0.041861 | 0.744088 |
| Dependent-child LTV 95% | 239 | 0.018494 | 0.041947 | 0.747284 |

![Supplemental resolution comparison. The two panels have different vertical scales. Both use global saving and the same inherited household measure.](output/model/e5f_overnight_independent_verification_20260905a/continuation_receipts/supplemental_nested_grid.pdf){width=90%}

The supply birth response changes by **-0.002486 pp**, or **-0.193%** of its
coarse-grid value; credit changes by **+0.000644 pp**, or **+3.61%**. The contrast
survives. Baseline births per household move from **0.084776385 to 0.084791202**
(**+0.01748%**); ownership moves by **+0.00443 pp**. The baseline birth-level
change is similar in size to the credit effect, making matched baselines essential.

The predeclared refinement triggers were 2% of the supply effect, 5% of the
credit effect, 1% in baseline births, or 0.1 pp in baseline ownership. None is
crossed: the branch stops after four scenarios, with **no 477-node launch**.
These are disclosed materiality thresholds, not a convergence certificate or
statistical test. The
[stop receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/nested_grid/stop_decision.json)
records every calculation.

This check leaves the wealth bounds, housing choices, and historical path fixed.
It does not re-estimate the model or the matched **2019-2023 first-birth rooms
target**, revise the 2063 path, or compute welfare. Full levels, effects, and
checkpoint hashes are in `continuation_receipts/summary.json` and its CSV tables.

\clearpage

## 18. Rental opportunities: a separate economic sensitivity

### Why test the six-room restriction?

The rental cap is an upper bound on continuously chosen rented space. The
production model sets it to six rooms, while owners choose among 2, 4, 6, 8,
and 10 rooms. This is an economic restriction on housing access. Refining wealth
interpolation leaves it in place; allowing larger rentals changes it.

The existing ACS/MMS housing-stock packet supplies a useful descriptive check.
Independent aggregation of its metro-location-size-tenure cells reproduces a
**9.1661%** share of renter households in the **seven-or-more-room** bin,
compared with **48.3463%** of owners. Thus its room distribution does not support
a literal universal six-room rental limit. It does not identify a replacement
upper bound, an elasticity, or a new target weight.

The packet pools **2012-2023** metropolitan household heads aged 18 or above,
with positive household weights and rooms, GQ codes 1 or 2, and a matched MMS
location. Middle PUMAs are assigned to the center. Its **4,103,889** records
reproduce when the mutually exclusive cells are summed. This is a saved-cell
arithmetic check and builder inspection, not a new raw-IPUMS extraction; the
sample also lacks the active ownership target's standard-structure filter.
No empirical standard error was computed from these saved cells.
The [builder](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/analyze_family_size_supply_menu.R:157)
and [stock cells](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_family_size_supply/acs_family_size_supply_cells.csv)
are recorded by hash in `continuation_receipts/summary.json`.

### Maintained experiment and baseline change

The diagnostic raises only the continuous rental limit from **six to eight
rooms**. Eight is an explicit diagnostic value, not an empirical estimate or
production default. All eleven fitted parameters, supply elasticity 0.63,
tenure dispersion 0.005, owner choices and their taste distribution, the 239-node
wealth grid, and the exact inherited 2023 pre-choice measure stay fixed.
Global saving and separately cleared current markets are used throughout.
The underlying history is not re-estimated under the new opportunity set.

Baseline smoke **16987040** passes the complete market loop in sixteen minutes,
with residual **\(1.56\times10^{-7}\)**, unchanged population, valid probabilities,
no occupied wealth-value declines above the declared screen, and all seventeen
standard graphs verified. The sole public-parameter change is independently
checked across the saved parameter records. Reporting-floor budget excess remains
negligible and separately reported.

| Baseline quantity | Six rooms | Eight rooms | Change |
|---|---:|---:|---:|
| Births per household | 0.08479120 | 0.08477440 | -0.01982% |
| Ownership, all households | 63.5671% | 61.6049% | -1.9622 pp |
| Rooms per household | 6.316906 | 6.341289 | +0.38599% |
| Housing asset price | 0.693013 | 0.697262 | +0.61304% |
| Renters at their own limit | 13.2451% | 6.6445% | Different limits |
| Parent renters at their own limit | 20.3493% | 6.8552% | Different limits |

Changing rental access moves tenure substantially at this inherited state,
while baseline births barely move. Aggregate rooms rise slightly, from an
already high baseline. This is not a demonstrated solution to the original
joint-fit problem. Prices and tenure composition also adjust; the baseline
birth change cannot be read as a pure housing-services derivative.

\clearpage

### Matched initial policy effects

Job **16987657** completes the two remaining cases in **33 minutes**. All three
cap-eight markets pass, with maximum residual **\(5.23\times10^{-6}\)**.
Every checkpoint and seventeen-graph packet is verified. Each effect below is
measured against its own rental-limit baseline; both limits use global saving
and 239 wealth nodes. Ownership covers all households; births and rooms are
per household in the current four-year period.

| Policy | Rental limit | Births, % | Rooms, % | Ownership, pp |
|---|---:|---:|---:|---:|
| Supply +20% | 6 rooms | 1.288387 | 5.759105 | 1.769706 |
| Supply +20% | 8 rooms | 1.292605 | 6.076940 | -0.232000 |
| Dependent-child LTV 95% | 6 rooms | 0.018494 | 0.041947 | 0.747284 |
| Dependent-child LTV 95% | 8 rooms | 0.013773 | 0.010498 | 0.626012 |

![Supplemental economic sensitivity. Every policy uses its own rental-limit baseline. Births and ownership use different units; panels have separate vertical scales.](output/model/e5f_overnight_independent_verification_20260905a/continuation_receipts/supplemental_rental_cap.pdf){width=90%}

**Supply's birth response barely changes, while its ownership response reverses
sign.** The birth effect moves from 1.2884% to 1.2926%; ownership changes from
**+1.7697 pp to -0.2320 pp**. Similar aggregate fertility responses can therefore
coexist with opposite tenure responses in this model. Tenure predictions need
their own housing-opportunity validation.

Credit's birth effect falls from **0.018494% to 0.013773%**, a **25.5%** reduction
in an already small coefficient. Its ownership response falls from **0.7473 pp
to 0.6260 pp**. The qualitative supply-versus-credit contrast survives, without
certifying the exact credit coefficient or identifying the correct rental limit.

Inherited input distributions and grids agree **exactly**. Existing feasibility
preparation moves **\(8.31\times10^{-31}\)** units of household mass under supply.
The processed-array identity rejection is retained; saved-policy replay is exact.
[The continuation receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/continuation_receipts/summary.json)
records both objects, hashes, and outcomes. Production gates remain unchanged.

This experiment does not recompute the matched **2019-2023 first-birth rooms
response**, SMM fit, history, 2063 transitions, or welfare. Value differences
across cleared prices are not a welfare test. Cap eight remains unpromoted;
its role is to identify housing opportunities that need empirical validation.
