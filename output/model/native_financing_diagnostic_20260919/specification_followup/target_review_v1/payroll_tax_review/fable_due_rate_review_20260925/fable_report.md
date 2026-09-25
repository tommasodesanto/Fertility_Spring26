# Provenance of DUE's 17.9% tax rate — final answer (2026-09-25, Fable review)

Scope: primary-source check of the tax rate in *Dynamic Urban Economics* (Greaney, Parkhomenko, Van Nieuwerburgh; "DUE") and of the OECD source it cites, compared with this project's proposed 8.751% OASI proxy. No model runs, no project-file changes other than `preliminary.md` and this file. Facts, arithmetic, and inference are labeled.

## Lead answer: where 17.9% came from

**Fact (three DUE versions compared).** The 17.9% is the **2012–2016 average of the OECD Taxing Wages "Average income tax rate"** for a single worker without children at 100% of the average wage. Both earlier public versions say so explicitly: the January 2025 AEA/ASSA program PDF (printed p. 26) and the February 2025 NBER Working Paper 33512 (printed p. 28) contain the sentence "We use the payroll tax rate reported by the OECD for the U.S. The average tax rate for the period 2012–2016 was 17.9%." with footnote 29: "See https://stats.oecd.org/index.aspx (variable 'Average income tax rate' located in 'Public Sector, Taxation and Market Regulation'/'Taxation'/'Taxing Wages'; accessed on September 14, 2023)." Those versions calibrate to 2012–2016 ACS data, and the February 2025 Table 1 (printed p. 25) prints the row `τ_z | payroll tax | 0.179 | OECD`.

**Fact (Dec. 5, 2025 version, local copy `tmp/pdfs/DUE_latest_checked.pdf`).** The calibration window moved to 2005–2009. The prose on printed p. 31 now reads "We use the payroll tax rate reported by the OECD for the U.S. The average tax rate for the period 2005–2009 was 15.6%." with footnote 33: "Accessed from https://data-explorer.oecd.org on May 15, 2025. Reference area: United States; Measure: Income tax; Household type: single person, no children; Earnings of the principal: 100% of average wage." Footnote 33 attaches to the prose sentence only. Table 1 (printed p. 27, "Externally calibrated" panel) still prints the identical row `τ_z | payroll tax | 0.179 | OECD`, with no footnote, no data year, and no value in the target column. The table note does not mention taxes.

**Conclusion (inference from the facts above): this is a version inconsistency, not two tax concepts.** The Table 1 entry 0.179 is the stale 2012–2016 number carried unchanged from the earlier versions; the prose was updated to the 2005–2009 window and the table was not. The paper does not say which value the Dec. 2025 model was actually solved with.

**Does the cited OECD source verify the printed numbers?** Partly.
- 17.9%: the current OECD series reproduces it approximately. The 2012–2016 mean of "Average income tax rate" (share of gross wage earnings) is **17.98%** in the vintage retrieved today, which rounds to 18.0%; the 0.08 pp gap versus the printed 17.9% is consistent with a Sept-2023 data vintage or truncation. Not exactly reproduced, but the identity of the object is verified by DUE's own footnote 29.
- 15.6%: the current series does **not** reproduce it as a gross-wage-based rate. The 2005–2009 mean of "Average income tax rate" is **16.95%**. However, dividing each year's income tax by labour costs (gross wage plus employer social security contributions), which is the OECD tax-wedge decomposition convention, gives a 2005–2009 mean of **15.58%**, matching the printed 15.6%. Footnote 33's measure label "Income tax" (rather than "Average income tax rate") is consistent with the Data Explorer's tax-wedge decomposition dataflow, whose components are shares of labour costs. This is my arithmetic plus inference; direct retrieval of the decomposition dataflow failed (HTTP 404/422/403 on six URL forms), so the labour-cost reading of 15.6% is strongly supported but not independently retrieved.

**Tax concept (fact, from OECD metadata).** The OECD measure is the **personal income tax** (central plus sub-central), unit "Percentage of gross wage earnings", household "Single person, no children", earnings "100% of average wage". It **excludes** employee social security contributions (a separate OECD series, 7.65% of gross wage for the US in every year 2013–2025) and employer contributions (a separate series, 8.1–8.9%). DUE's label "payroll tax" is therefore a misnomer relative to the source: the number is an average personal income tax rate, not a Social Security or Medicare contribution rate.

## Comparison table

| Object | Value | Definition, tax base, denominator | Period / vintage | Source location | What it finances in that model |
|---|---:|---|---|---|---|
| DUE Table 1, τ_z | 17.9% | OECD "Average income tax rate": personal income tax (central + sub-central) for a single childless worker at 100% AW, as % of gross wage earnings; excludes employee/employer SSC. Stale from the 2012–2016 calibration. | 2012–2016 mean; OECD.Stat accessed 2023-09-14 (per Feb 2025 fn 29) | Dec 2025 version Table 1 p. 27 (source column "OECD", no footnote); NBER Feb 2025 Table 1 p. 25, prose p. 28 and fn 29; AEA Jan 2025 p. 26 and fn 29 | Flat tax on working-age earnings only (eq. 2.1); revenue "go[es] toward pension transfers and wasteful government spending" (pp. 10–11). No budget constraint in the equilibrium definition (p. 20). |
| DUE prose, p. 31 | 15.6% | "Income tax", same household/earnings spec. Reproduced only as income tax ÷ labour costs (gross wage + employer SSC), i.e. the tax-wedge decomposition component; not reproduced as % of gross wage (16.95%). | 2005–2009 mean; data-explorer.oecd.org accessed 2025-05-15 | Dec 2025 version p. 31, fn 33 | Same as above. |
| OECD DF_TW_COMP, AV_ITR, USA, S_C0, AW100 | 2005–09 mean 16.95%; 2012–16 mean 17.98%; 2023: 16.68%; 2024: 16.75%; 2025: 16.66% | Average income tax rate, % of gross wage earnings (unit PT_WG_EARN_G) | Annual 2000–2025; retrieved 2026-09-25 from sdmx.oecd.org (dataflow v2.1) | SDMX URL in appendix | Not program-specific; a general personal income tax measure. |
| Same series converted to labour-cost base (my arithmetic) | 2005–09 mean 15.58%; 2024: 15.50% | AV_ITR / (1 + employer SSC rate), i.e. income tax as % of labour costs | 2005–2009 | Arithmetic below | Same. |
| OECD employee SSC, same spec | 7.65% (2000–2010 and 2013–2025; 5.65% in 2011–12) | Employee social security contributions, % of gross wage (OASDI 6.2 + HI 1.45) | Annual | Same dataflow, AV_R_EMPEE_SSC | Statutory OASDI+HI employee side; program-wide, not OASI-only. |
| This project's proposal | 8.7510174% | 2007 combined employee+employer OASI statutory rate 10.6% × (taxable covered earnings / uncapped covered earnings); flat rate applied to uncapped household gross working earnings; no cap, no employer wedge; excludes pension and asset income | 2007, SSA 2025 statistical vintage | `paygo_rate_review_20260924.md`; SSA 2025 Supplement Tables 4.B1–4.B2 | All receipts fund the common PAYGO pension (pension adjusts to exhaust receipts). |

## Arithmetic (all inputs are retrieved OECD/SSA values; date convention: OECD tax year = calendar year)

OECD "Average income tax rate", USA, single no children, 100% AW, % of gross wage (retrieved 2026-09-25):

| Year | AV_ITR | Employer SSC (% gross wage) | AV_ITR / (1 + employer SSC) |
|---|---:|---:|---:|
| 2005 | 16.621145 | 8.834222 | 15.2720 |
| 2006 | 16.760309 | 8.843590 | 15.3985 |
| 2007 | 17.191223 | 8.795160 | 15.8015 |
| 2008 | 17.674972 | 8.758899 | 16.2515 |
| 2009 | 16.522381 | 8.758477 | 15.1918 |
| 2012 | 17.987980 | 8.844236 | — |
| 2013 | 17.803432 | 8.748946 | — |
| 2014 | 17.950062 | 8.837450 | — |
| 2015 | 17.981947 | 8.466867 | — |
| 2016 | 18.172957 | 8.415232 | — |

- 2005–2009 mean of AV_ITR = 84.77003 / 5 = **16.954%** (does not match 15.6%).
- 2005–2009 mean of the labour-cost-based column = 77.9153 / 5 = **15.583%** (matches the printed 15.6%).
- 2012–2016 mean of AV_ITR = 89.896378 / 5 = **17.979%** (printed 17.9%; rounds to 18.0% in this vintage).
- Identity check, 2005: (16.621145 + 7.65 + 8.834222) / 1.08834222 = 30.418%, equal to the retrieved average tax wedge 30.418159%, confirming the labour-cost denominator convention.
- Project proxy: 5,268,200 / 6,381,306 = 0.8255677; × 0.106 = **0.0875102** (dimensionless share of gross working earnings; SSA amounts in millions of current dollars, 2007 rows, 2025 vintage).

## DUE fiscal structure (Dec 2025 version, exact locations)

- Section 2.1.3, pp. 10–11: "Retirees receive a pension benefit 𝔟. Earnings are subject to a payroll tax of rate τ_z. Payroll tax revenues go toward pension transfers and wasteful government spending."
- Equation (2.1), p. 11: after-tax labor income is the discounted integral of (1 − τ_z) e^{z(ζ,a)} w_js for a < A_ret, and 𝔟 otherwise. So τ_z applies to working-age labor earnings only; pensions and asset income are untaxed; there is no cap and no employer wedge.
- Equation (2.4), p. 13: property tax τ_h on owner-occupied housing; equation (2.27), p. 20: τ_h also enters REIT rent pricing. τ_h = 0.0071 (Table 1). These are the only two taxes in the model.
- Equilibrium definition, Section 2.4, p. 20: six conditions (household optimization, density consistency, firm optimization, labor-market clearing, floorspace clearing, REIT rent equation). **No government budget constraint appears.** τ_z and 𝔟 are not linked by any displayed equation; the residual between revenue and pensions is absorbed by "wasteful government spending" by assumption.
- Pension level: 𝔟 has no row in Table 1 and no calibration sentence in Section 3.2 (pp. 26–32). Its numerical calibration was not located in the pages read (pp. 9–32) and remains unestablished.

Implication: DUE does not have a pension-financing closure. All revenue is not pension funding; the paper does not state the pension share.

## Comparison with the project's 8.751%

- Concept: DUE's number is an average **personal income tax** rate on a representative single worker; ours is a **program contribution** (OASI) rate scaled to the effective taxable base. They are different tax objects; neither is an approximation of the other. The OECD's own Social Security component for the same worker is the 7.65% employee series plus the employer series, and both are OASDI+HI, not OASI.
- Base: DUE applies its rate to individual efficiency earnings in a single-worker model; our model applies the rate to uncapped household gross working earnings. Both exclude pensions and asset income and have no cap or employer wedge, so the base treatment is structurally similar; the difference is what the rate is meant to measure.
- Use of revenue: DUE sends revenue to an exogenous pension plus a wasteful residual with no budget equation; our model sends all receipts to a common PAYGO pension that adjusts to exhaust them. Importing 17.9% (or 15.6%) into our closure would fund pensions with personal-income-tax-scale revenue, roughly double an OASI-scale contribution, which is why the earlier local comparison shows pension/mean-earnings ratios of about 51% versus 25%.

## Recommendation for the paper

1. Neither 17.9% nor 15.6% can calibrate a **pension-only** levy. Both are OECD personal income tax averages (different periods; the 15.6% additionally appears to be on a labour-cost base). Describing either as a payroll or pension contribution would misstate the source.
2. Using 17.9% would require adopting DUE's structure: a flat general labor-income tax, an exogenous pension, and an explicit residual for non-pension (wasteful) spending, with no PAYGO budget equation. That is a different fiscal model, not a re-parameterization.
3. If the paper keeps a pension-only PAYGO levy, the OASI-based proxy is the right kind of object; the unresolved mismatch is institutional: the model's uncapped household earnings versus SSA taxable covered earnings, and attributing both statutory sides to households. Those caveats are already recorded in `paygo_rate_review_20260924.md` and are not resolved by DUE.
4. If DUE is cited for τ_z at all, cite the OECD series directly (dataflow DSD_TAX_WAGES_COMP@DF_TW_COMP, USA, S_C0, AW100, measure AV_ITR), name the window, and label it personal income tax. Do not cite DUE's Table 1 value as a payroll rate.

**Strongest unresolved caveat.** The 15.6% is reproduced only under the labour-cost denominator reading; the OECD decomposition dataflow that would confirm the "Income tax" measure label and unit could not be retrieved (404/422/403 on six URL forms). A one-call check of `DSD_TAX_WAGES_DECOMP@DF_TW_DECOMP` for USA/S_C0/AW100 2005–2009 (for example through the Data Explorer UI) would close it. Also unverified: which value (0.179 or 0.156) the Dec 2025 model was actually solved with, and DUE's pension level 𝔟.

## Primary-source appendix

- DUE, Dec. 5, 2025 public version (local copy `tmp/pdfs/DUE_latest_checked.pdf`, SHA-256 in `lead_review.json`; author page https://sites.google.com/view/briangreaney/research, Drive link https://drive.google.com/file/d/1RLT_B_MBqpr2NrYDPjFcB-53kATJj5ly/view?usp=drive_link, listed as "American Economic Review, Revise and Resubmit (2nd Round)"): Section 2.1.3 pp. 10–11; eq. (2.1) p. 11; eq. (2.4) p. 13; Section 2.4 p. 20; Table 1 p. 27; "Taxes" paragraph pp. 30–31 and fn 33 p. 31.
- DUE, NBER Working Paper 33512, February 2025: https://www.nber.org/system/files/working_papers/w33512/w33512.pdf — Table 1 printed p. 25 (row `τ_z payroll tax 0.179 OECD`); printed p. 28 ("2012–2016 was 17.9%"); fn 29 (OECD.Stat, "Average income tax rate", accessed Sept 14, 2023).
- DUE, AEA/ASSA 2025 program version (January 2025): https://www.aeaweb.org/conference/2025/program/paper/R5DzakN3 — printed p. 26, same sentence and fn 29.
- OECD SDMX, comparative indicators, USA, single no children, 100% AW, all measures from 2000: https://sdmx.oecd.org/public/rest/data/OECD.CTP.TPS,DSD_TAX_WAGES_COMP@DF_TW_COMP,/USA...S_C0.AW100._Z.A?startPeriod=2000&dimensionAtObservation=AllDimensions&format=csvfilewithlabels (measures AV_ITR "Average income tax rate", AV_R_EMPEE_SSC, AV_R_EMPER_SSC, AV_TW, NPATR; unit PT_WG_EARN_G "Percentage of gross wage earnings", PT_COS_LB for the wedge). Income-tax-only subset: same URL with key `USA.AV_ITR..S_C0.AW100._Z.A`. Raw rows were re-read from the downloaded CSV, not only from the fetch summary.
- OECD Taxing Wages 2025 country note, United States (tax wedge formula with labour-cost denominator; 2024 wedge 30.1%; employee net average tax rate 24.4%): https://www.oecd.org/content/dam/oecd/en/publications/reports/2025/04/taxing-wages-2025-country-notes_16d47563/united-states_e44f6362/4a2a60b9-en.pdf
- SSA 2025 Supplement Tables 4.B1–4.B2 (2007 rows used for the 8.751% proxy): https://www.ssa.gov/policy/docs/statcomps/supplement/2025/4b.html#table4.b1
- Failed retrievals (stated for completeness): OECD decomposition dataflow `DSD_TAX_WAGES_DECOMP@DF_TW_DECOMP` (six URL forms: 404/422/403); OECD Taxing Wages 2026 US HTML page (403); the Parkhomenko-site PDF (over size limit); OECD DSD codelists (truncated responses).
