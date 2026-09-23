# Second pass: corrected conclusions on the calibration targets

**Access correction.** I wrongly reported that the daily memory files do not exist. `memory/daily/2026-09-22.md` reads fine at its direct path; my file search did not follow the memory symlink. That was an access limitation on my side, not a missing file.

Scope of this pass: I read code and sources only. No model runs, no data writes, and nothing that duplicates the four Luna workers.

## (a) Claim, evidence, status

Evidence labels E1–E13 are listed with absolute paths at the end.

| # | First-report claim | What the code and sources show | Status |
|---|---|---|---|
| 1 | 2.1 and the two CPS rows cannot all hold. B's share of women with 3+ children is ≈0.66. Part of the CPS-row loss is mechanically forced. | The two sides measure different women. The 2.1 normalization applies to top-coded completed fertility, $p_1+2p_2+T\,p_{3+}$ with $T=3.602359422009$, over all households in age cells 46 and older, weighted by population mass (E1, E2). The fertility utility scale ψ is solved until that statistic equals 2.1 within 5e-4 (E3). The CPS rows instead come from a 40–44 projection that assumes births fall uniformly within each four-year cell, applied to the model's household reproductive-member mass, which is "not certified female exposure" (E4). The frozen contract already records the same-population algebra with weighted CPS data: 2.1 would require a 3+ share of 0.4168 against an observed 0.2863. It labels this an "algebraic conditional diagnostic only" (E5). | **Withdrawn**: the impossibility claim, B's 3+ share, and the forced loss. What survives is conditional algebra that the contract already contains. |
| 2 | Unweighted counts: 3+ share ≈0.29–0.30, capped mean within the 3+ group ≈3.50. | The contract has weighted values: 0.2863 and 3.4827933 (E5). Check: at those weighted shares, the observed capped mean of 1.8566 is reproduced; with $T=3.602$ the same shares give ≈1.891. | **Withdrawn**; use the contract values. |
| 3 | 3.602 is the 2024 vintage. | The June 2024 CPS builder defines it as $E[\min(N,5)\mid N\ge3]$ for women 40–44 (weight PWSSWGT; cohorts roughly 1979–84) (E6). The 2004/06 equivalent is 3.4828 (E5), and the design report already flags the mismatch. | **Verified.** Which vintage to use is an author decision. |
| 4 | Role of 2.1. | It is an author-selected stationary benchmark normalization. It is unscored, pins ψ, and also enters the population renewal law ("divides topcode-adjusted births by 2.1"). It is not a measured period TFR, and the author "has not approved treating the old completed fertility observation as literally equal to the 2.1 normalization" (E5, E7). | **Verified.** It is not an empirical target and not a policy normalization. |
| 5 | B's entry law departs from July, presented as a new finding. | B's "inherited heterogeneous marginal with diagnostic income-rank coupling" is a declared experimental design. Harmonizing the income denominator is already listed as pending (E8, E9). Daily memory records that inherited negative entry debt failed the feasibility gate under the new income support earlier that day (E13). July's 18–24 sample stands as reviewed. | **Corrected.** It is not new and not arbitrary. Quantifying the conversion gap remains useful. |
| 6 | Pass rules for entry: 1 SE, 20%, Spearman 0.3. | These were my own thresholds, and they conflate sampling uncertainty with economic sensitivity. | **Withdrawn.** |
| 7 | Timing wedge between beginning-of-period wealth and the PSID interview date ≈ +1.7%. De Nardi–Yang's 6.90 on a gross basis ≈ 5.3–5.5. PSID's omission of employer DC accounts biases the target low. | The 1.7% came from a simplified stationary identity that ignores housing transactions, prices, transaction costs, entrant normalization and survival timing. The 5.3–5.5 figure used an assumed tax rate. The DC point has no stated contract for which assets the model's b covers, so no direction follows. | **Withdrawn.** The DC omission stays as a coverage question (Boston Fed working paper). |
| 8 | Asset income in family income compresses the data's old-age p90/p50. | This holds only under the stated transform $x/(1+rx)$ with a common return r. Heterogeneous returns and income sources remove the sign. | **Downgraded** to conditional algebra. |
| 9 | Survey-reported bequest flows ≈0.5% of wealth, making the 5% tolerance ≈10× too tight. | 0.5% used an assumed net worth / disposable income multiple. The FEDS Note reports only ≈3% of disposable personal income. | **Withdrawn**, including the multiple. |
| 10 | 0.0088 comes from Gale–Scholz (primary not seen). | Now verified in the working-paper version. Table 4 (p. 13): bequests $105.00B per year = 0.88% of 1986 SCF net worth ($11,976B). Method (p. 12): mortality-imputed from 1986 net worth excluding pensions, less trusts. Children receive 25% of the estate if the head is married and 75% if the head is single or both spouses die. Life insurance ($7.84B) is separate. | **Verified (working paper).** The concept is narrower than "all estates". The published JEP (1994) table was not re-checked. |
| 11 | The first-birth room observer spans 0–8 years after the birth, mean ≈4. | The code fixes the comparison in model dates, not calendar years (see F1 below). | **Withdrawn.** |
| 12 | The range of alternative coefficients (0.43–0.81) bounds the target; a gap over 0.2 rooms proves economics; the parenthood housing floor causes B's 1.34. | 0.426502 is a different estimand (see F2). Alternative specifications are not draws around the matched-branch estimand. No attribution exercise exists. | **Withdrawn.** |
| 13 | Five coordinates are at or near bounds. | The flag is $\min(v-lo,\,hi-v)\le 0.01\,(hi-lo)$, measured in levels (E11). On [0.02, 50] the threshold is 0.4998, so κ_fert = 0.293 and κ_fert continuation = 0.454 are flagged although interior on their log search scale. θ₀ (threshold 0.08) and θ₁ (0.16) are flagged the same way. Only h_P = 2.3 sits exactly at its bound. | **Corrected.** A flag is not a binding optimum. |
| 14 | Kaplan–Mitman–Violante and De Nardi–French–Jones are precedents for counting current income toward the down payment. | The Kaplan–Mitman–Violante article shows a loan-to-value cap at origination (eq. 5, p. 3294) and a payment-to-income cap on current income (eq. 6, p. 3295). The buyer's full budget constraint is in their Appendix A, which I could not verify. That is a different object from immediately available down-payment cash. De Nardi–French–Jones has no housing purchase. | **Withdrawn** as precedent. Our Y/R eligibility is reconciled to the budget, but four-year income known at the decision remains an economic timing assumption. No general-equilibrium direction is claimed. |

**Corrected narrow findings**

- **F1. What the first-birth room observer computes** (E12).
  - At model date t, for fertile cells covering ages 18–42, treated mass is fecundity × attempt probability × the childless "settled" mass in the pre-fertility distribution. The control is identical mass left childless (lines 1201–1212).
  - Both branches then advance one period: date-t location, tenure, saving and housing choices, plus survival, income and child-maturation transitions (lines 1228–1273).
  - At t+1 the treated branch goes through the fertility step, so second births are allowed. The control is held childless with no fertility step (lines 1353–1382).
  - Both realize the t+1 cross-section under the t+1 policies, and the outcome is uncapped mean realized housing (lines 1383–1424).
  - The estimand is therefore $\Delta=E[h_{t+1}\mid\text{first birth at }t]-E[h_{t+1}\mid\text{same state, no birth at }t\text{ or }t+1]$. Because the pre-states are identical, this equals a difference in changes from the housing held entering t.
  - The data counterfactual is different: never-treated women with person and year fixed effects and cohort weights, identified under parallel trends.
  - The code contains no within-period birth date. Mapping model dates to event years is therefore an explicit convention that someone must choose.
- **F2. The two empirical numbers are different estimands** (E10).
  - **0.7202462624** is from `sa_rooms_first_birth_household_aligned_v1`: β(+3) − β(−1), with −2 (and −6) omitted, person and year fixed effects, age and education controls, IW weights, one woman per household-year, single-family-unit dwellings, confirmed never-treated controls. It uses 49,457 household-years and 4,112 women, and already includes the rooms-timing correction ("current code moves that value forward", line 53). Its SE of 0.08526 is computed from the covariance of the two coefficients.
  - **0.426502** is a binned contrast: +3/+4 against −3/−2, verified timing assignment, the 2019 cohort as sole control, only cohorts observed in all six windows. It uses 149,402 observations and 14,453 clusters.
  - On September 12 the author "reaffirmed −2, not −1, as the main empirical comparison" (line 433), while the scored contrast subtracts β(−1). Whether the author meant the reference year or the scored contrast is unresolved.
- **F3. The bequest external restriction (0.0088) is now pinned down.** It is the children's share of mortality-imputed net-worth estates, on the 1985 population aged 25+. The model's numerator is the full non-negative household estate at death. These are not the same convention.

## (b) Decisions that genuinely need the author

**1. Fertility level and top-bin convention.**
- *Default:* keep 2.1 unscored and T = 3.602 frozen tonight. Add visible diagnostics: the model's distribution of women over 0/1/2/3+ children at the 40–44 observer and at 46+, next to the weighted CPS distribution.
- *Alternative:* use the early-vintage T = 3.4828 with 2.1 kept, or run a robustness benchmark normalized to the same-window CPS level (1.8566 capped). Either is a different benchmark, not a correction.
- *What it identifies:* the normalization fixes ψ. κ_fert, κ_fert continuation and the first-birth cost face the CPS shares and the timing rows. T also enters the renewal law, so population closure moves.
- *Evidence missing:* the model's distributions at both ages; the sensitivity of ψ to T (needs normalized solves, not tonight).

**2. Entry-wealth conversion.**
- *Default:* keep B's declared level-preserving, rank-coupled marginal for the frozen comparison.
- *Alternative:* the July convention (each ratio × annual gross income at the entrant's state under the new process), or a joint wealth–earnings distribution re-estimated on the same 18–24 sample.
- *What it affects:* young ownership, rooms and first-birth timing through affordability; parameters β, χ, first-birth cost and h_P. No target changes.
- *Evidence missing:* implied entry wealth/income by rank under B; same-sample earnings versus family income; whether the lowest entry-wealth group (quintile mean −2.22 × income) is feasible under the new income support.

**3. First-birth room estimand and observer dating.**
- *Default:* keep 0.7202 and the matched-branch observer tonight.
- *Alternatives:* score β(+3) against the −2 reference; measure the model outcome at the birth date t rather than t+1; or apply the event-study estimator to simulated model panels, which gives the model the same counterfactual as the data.
- *What it identifies:* h_P (at its bound), the first-birth cost and the per-child rooms parameter (fixed at zero in B). It interacts with the family-rooms and recent-parent rows.
- *Evidence missing:* β(+3), β(−1) and their covariance from the v1 regression; the model split of the response between t and t+1 and by second births.

**4. Housing geography and the household/age unit.**
- *Default:* 42-metro ACS rows (the maintained approximation).
- *Alternative:* national ACS rows (rooms 5.6079, ownership 0.6763, recent-parent gap 0.1276), matching the national PSID, CPS and NCHS rows.
- *What it affects:* H₀, χ and the family-housing rows. IPUMS MET2013 assigns households by the PUMA's population majority and suppresses metros with at least 15% match error, so "42 metros" is itself a PUMA-based approximation. The age concept differs by row: householder in the ACS, women in CPS/NCHS and the first-birth panel, reference person in PSID wealth.
- *Evidence missing:* which population the price and supply block represents; national bootstrap standard errors.

**5. Wealth and bequest rows: mapping choices.**
- *Default:* keep 0.0088 (signed off July 24) with its verified definition recorded, and keep the old-age ratio row.
- *Alternative:* align the model's estate numerator to the children's-share convention, or state the synthetic 5% scale explicitly as a tolerance. For old age, add a p90/p50 of wealth levels as a robustness row. Replacing either row requires a θ₀/θ₁ Jacobian check first.
- *What it identifies:* θ₀, θ₁, and β through the wealth rows.
- *Evidence missing:* the published JEP table; model estates by age at death; the old-wealth receipt.

**6. Within-period income and purchase timing.**
- *Default:* keep Y/R eligibility with $b'\ge-\phi pH$ for the frozen comparison.
- *Alternative:* allow only beginning-of-period cash plus sale proceeds for the down payment, as a declared robustness specification.
- *What it affects:* how much the $(1-\phi)$ down-payment threshold binds for young renters, and how much entry wealth matters. No equilibrium sign is established.
- *Evidence missing:* a parameter-fixed controlled comparison.

## (c) Questions the incoming receipts can settle

**CPS/NCHS receipt**
1. Weighted CPS 2004+2006 shares with 0, 1, 2 and 3+ children, and the capped mean within the 3+ group. Is it 3.4828?
2. B's distribution at the 40–44 observer and at 46+. How large is the difference, and does the 40–44 mean evaluated at T = 3.602 fall short of 2.1 by more than births after 44 account for?
3. The NCHS first-birth mean computed from single-year counts, with and without +0.5, against the published 25.0 for 2006. How much comes from collapsing births under 18 (a 0.0773 share) into the first model cell?

**Entry-units receipt**
4. Implied entry wealth / B's annual gross income at 18, by persistent state and pooled: mean, median and quintile means, set against July's 0.2589, 0.0985 and the five quintile means.
5. In the same 1,835 family-years: the ratio of summed reference person + spouse gross earnings to summed family income, the quantiles of wealth/earnings, and the count of zero or low earners.
6. The descriptive weighted rank correlation between wealth and earnings, with no pass/fail threshold.

**Old-wealth receipt**
7. p90/p50 of net-worth levels at ages 76–84 in 2003/05, with bootstrap SEs. The current ratio without the children-history filter. Does the PSID-SHELF file carry any income component that excludes asset income?
8. The model's levels-based p90/p50 of b + pH over the same age overlap.

**First-birth receipt**
9. From the v1 regression: β(+3), β(−1), β(+3) − β(−1), and their covariance.
10. For the model: Δ decomposed by origin age; the treated–control difference at date t; the share of treated households with a second birth at t+1; and Δ excluding those second births.

## (d) Verified primary citations

- **De Nardi & Yang (2014)**, "Bequests and heterogeneity in retirement wealth," *European Economic Review* 72: 182–196. Table 2 and targeting text, p. 187; demographic timing, p. 185; Table 1, p. 186. [PDF](https://users.nber.org/~denardim/research/De-Nardi_Yang_EER.pdf)
- **Gale & Scholz**, IRP Discussion Paper 1019-93 (September 1993), the working-paper version of *JEP* 8(4) (1994): 145–160. Table 4, p. 13; method, p. 12; Table 6, p. 18. The published table was not checked. [PDF](https://www.irp.wisc.edu/publications/dps/pdfs/dp101993.pdf)
- **Kaplan, Mitman & Violante (2020)**, *JPE* 128(9): 3285–3345. Eq. (5), p. 3294; eq. (6), p. 3295; household decisions, p. 3297; Table 1, p. 3303; Table 2, p. 3304; entry wealth, p. 3305. Appendix A not verified. [PDF](https://violante.economics.princeton.edu/sites/g/files/toruqf5621/files/documents/kaplan-et-al-2020-the-housing-boom-and-bust-model-meets-evidence.pdf)
- **Sommer & Sullivan (2018)**, *AER* 108(2): 241–274. Targets, pp. 256–257; Table 5, p. 258. [PDF](https://kamilasommer.net/Taxes.pdf)
- **De Nardi, French & Jones (2010)**, *JPE* 118(1): 39–75. Eq. (18), p. 47; p. 48. [PDF](https://users.nber.org/~denardim/research/De_Nardi_French_Jones_JPE_2010.pdf)
- **Sommer**, FEDS Working Paper 2014-32. Table 2, p. 20; Table 4, p. 21. The published *JME* 83 (2016): 27–38 tables were **unavailable** (paywall). [PDF](https://www.federalreserve.gov/pubs/feds/2014/201432/201432pap.pdf)
- **Reports and documentation (not papers):** NCHS Data Brief 21 ([link](https://www.cdc.gov/nchs/products/databriefs/db21.htm)); Feiveson & Sabelhaus FEDS Note, 2018 ([link](https://www.federalreserve.gov/econres/notes/feds-notes/how-does-intergenerational-wealth-transmission-affect-wealth-concentration-20180601.html)); Cooper, Dynan & Rhodenhiser, Boston Fed WP 19-6, abstract page only ([link](https://www.bostonfed.org/publications/research-department-working-paper/2019/measuring-household-wealth-in-the-panel-study-of-income-dynamics-the-role-of-retirement-assets.aspx)); IPUMS [NCHILD](https://usa.ipums.org/usa-action/variables/NCHILD), [ELDCH](https://usa.ipums.org/usa-action/variables/ELDCH) (counts step-, adopted and biological children of any age; code 99 = none present), [MET2013](https://usa.ipums.org/usa-action/variables/MET2013).

**Evidence locations**
- E1: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/earnings_entry_battery_20260922_v1/source/code/model/intergen_eqscale_seq/solver.py:5294-5298`; `.../intergen_eqscale_seq/parameters.py:33-34`
- E2: `.../tmp/earnings_entry_battery_20260922_v1/source/code/model/intergen_eqscale_seq_optimized/calibration.py:1281-1302`
- E3 and E12: `.../tmp/earnings_entry_battery_20260922_v1/source/code/model/tools/run_e5f_transition_calibration.py:765-806, 1170-1432`
- E4: `.../output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/target_fit.csv`
- E5: `.../tmp/earnings_entry_battery_20260922_v1/inputs/objective_source_files/fertility_provenance_file_sha256/fertility_target_contract.json:283-305`
- E6: `.../tmp/paper_baseline_sep14/code/data/cps_fertility/README.md:1-50`
- E7: `.../output/model/e5f_matched_pf_20260909a/OVERNIGHT_PLAN.md:8,22,76,90`
- E8: `.../earnings_entry_battery_v1/README.md:13-15`
- E9: `.../code/model/tools/e5f_earnings_wealth_contract.py:16-79`
- E10: `.../code/data/psid_followup_mar2026/output/first_birth_correction_review/README.md:50-60, 396-406, 433-443, 546-565`
- E11: `.../tmp/earnings_entry_battery_20260922_v1/inputs/validator/panel_validator.py:59`; `.../inputs/scorer/score_initial.py:190,280`
- E13: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/memory/daily/2026-09-22.md:17-21`

All running calibration inputs remain frozen. This pass is not paper adoption.
