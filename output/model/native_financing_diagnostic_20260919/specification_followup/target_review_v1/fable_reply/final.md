**1. Bottom line for tonight.** Yes, run the floor-versus-child-share comparison tonight, but strictly as a within-design contrast under one frozen, labeled convention: single AR(1) with the fitted four-year coefficients, common age profile, entry at 18 with a stationary income draw, cell B's inherited 18–24 wealth marginal with rank coupling, current-income purchase accounting, payment-to-income screen off, and the frozen 13-row target and weight system. The deliverable is two complete target tables plus the 17 diagnostics side by side. Nothing in the battery certifies B as a calibration: the floor sits at its 2.3 upper bound, and both bequest parameters sit near lower bounds.

Selected B fit, the frozen system the comparison would inherit:

| Row | Target | Model | Loss |
|---|---|---|---|
| Normalization, completed fertility | 2.1 | 2.1000 | unscored |
| Childless women 40–44 | 0.1983 | 0.1689 | 30.7 |
| Exactly one child among mothers | 0.2137 | 0.2553 | 46.7 |
| Mean first-birth age | 25.98 | 27.30 | 243.6 |
| First births at 30+ | 0.2493 | 0.3020 | 38.5 |
| Wealth / annual gross earnings | 6.146 | 6.059 | 0.06 |
| Bequest flow / wealth | 0.0088 | 0.00633 | 31.5 |
| Old p90/p50 wealth-income, 76–84 | 3.516 | 3.728 | 0.48 |
| Mean rooms, cap 9 | 5.561 | 6.245 | 60.0 |
| Ownership 30–55 | 0.648 | 0.453 | 89.0 |
| First-birth room response | 0.720 | 1.341 | 53.1 |
| Rooms, 3+ vs 1–2 children | 0.347 | 0.306 | 0.48 |
| Recent-parent ownership gap | 0.163 | 0.118 | 54.9 |
| Total | | | 648.9 |

Searched coordinates: beta 0.986, chi 1.009, H0 7.58, first-birth cost 0.315, kappa_fert 0.293, kappa_fert_cont 0.454, theta0 0.040 (lower bound 0), theta1 0.092 (lower bound 0.02), h_P 2.3 (at upper bound). Full 17-row parameter table: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/parameters_actual_bounds.csv`.

**2. Household timeline.**

Current implementation, from the receipts and solver comments. A period is four years. Age 18 labels the cell covering ages 18–21. The state entering a cell is net financial wealth b carried from the previous period, tenure and housing, children ever born and at home, and the persistent income state z. For owners b includes mortgage debt, so b can be negative. The cell's z is realized at the start, so the whole four-year after-payroll-tax flow is known when the household picks consumption, housing, tenure, fertility, and next-period wealth. Only future z and survival after 66 are uncertain. Purchase requires cash available before purchase to cover (1−φ)pH and end-of-period wealth not below −φpH. The period's income is part of that cash. Entrants draw z from the stationary AR(1) distribution and b from the inherited five-node marginal, rank-coupled to z. I did not re-read the compiled tenure kernel tonight. The cash composition is taken from the battery README and the solver docstring.

Proposed convention to freeze, no code change tonight. Age 18 is household formation at the start of the 18–21 cell. b_18 is the pre-income, pre-choice stock. Period income is known; future income and survival are not. Keep full-period resources for the down payment as the stated convention and record beginning-of-period-cash-only as a later sensitivity. Two reasons: both arms share it, and it is generous to ownership while the model under-predicts ownership by 20 points. Timing is therefore not the mechanism behind the housing misses. Tightening it would move the wrong way. The stationary draw at 18 makes earnings dispersion flat over age instead of rising. Disclosed limitation, not tonight's problem.

18–24 wealth. Acceptable as a frozen diagnostic initial condition. Misleading as a paper initial condition. Acceptable because the population concept is right: the model's 18-year-olds are independent households, since there is no living-with-parents state, so conditioning on PSID childless renter heads is correct. The denominator error is second order for this group. Young childless renters' family income is mostly their own labor earnings plus transfers, so ratios are understated by roughly the transfer share, not by a factor. Misleading for adoption for three reasons. Rank coupling imposes perfect wealth-income rank correlation at entry. It is imposed, not measured, and it decides who can buy early. The bottom node is debt of 2.2 years of income for one fifth of entrants, and whether renters may roll unsecured debt across periods is not documented in the receipts I read. Heads aged 18–21 are a selected minority.

Minimal construction. Rerun the existing PSID entry builder on the same 18–24 childless renter sample, with the denominator set to the same head-plus-spouse annual gross labor earnings concept used for the wealth target and the income process. Tabulate the joint distribution: wealth-ratio quintile nodes within earnings terciles, mapped to model income ranks. That closes denominator and coupling in one pass with no transfer, floor, or fallback.

**3. Wealth targets.**

Wealth to earnings 6.146. Useful: pins beta given the income process, with De Nardi–Yang precedent. Mismatch: PSID NETWORTHR includes business equity and vehicles the model lacks, and the model counts b plus owner pH. The earnings variable's coverage must be confirmed: if it is head-only rather than head-plus-spouse, the ratio is biased up by roughly the spouse share. The window matters: 2003/05 gives 6.15, 2005/07 gives 6.93, more than two standard errors apart, because house prices moved. Correction: match the earnings concept to the income process, strip business and vehicles if the builder allows, and declare 2003/05 as the pre-announcement economy. Not tonight; before adoption. Affects beta, with H0 and chi cross-loading through pH.

Old dispersion 3.516. Intended for theta1, weakly. Mismatch: PSID family income includes Social Security, pensions, asset income and transfers; the model denominator is pension only. Asset income is the main wedge at p90, so the model ratio should exceed the data for the same wealth, and B shows exactly that. Correction: change the data side, not the model. Either recompute with Social Security plus pension income as denominator, or move to the denominator-free p90/p50 of net worth levels at 76–84. Both keep the same identifying block. Not tonight; before adoption, since theta1 sits near its bound and was the least-identified direction on September 7.

Bequest flow 0.0088. Useful: the only level moment for theta0. Mismatch: external 1986 US aggregate flow against model estates of decedents 66+ with forced terminal death and nonnegativity. With theta0 near 0.04, model bequests are almost entirely accidental. The 5 percent synthetic scale yields weight 5.2 million; the 0.0025 gap costs 31.5 loss units, the same order as the childlessness miss. Correction: keep as an external restriction with explicit author acceptance and a stated cross-study tolerance instead of 5 percent. Not tonight; before adoption.

Denominator rule: harmonize each target to the model's own income concept at that age. Gross labor earnings is right for the working-age aggregate and for entrants. It is wrong for the 76–84 row, where labor earnings are near zero.

**4. Source check on 0.0088.** Verified at the intermediate. De Nardi and Yang, NBER WP 20058, Table 2, p. 12: "Bequest-wealth ratio, Data 0.0088, Model 0.0089", attributed in the text to Gale and Scholz (1994). The same table sources wealth-to-after-tax-earnings 6.90 to Hendricks (2007a) and the p90 bequest 4.53 to Hurd and Smith (2002), single decedents. Their model has five-year periods, entry at 20, no death before 65, and no inter vivos transfers (p. 7 and footnote 1). Not verified: the Gale–Scholz primary. Three fetches failed, two AEA 403s and an unreachable mirror. Numerator coverage, bequests versus inter vivos gifts, gross versus net of estate tax, spousal transfers, and the 1986 wealth base remain unverified by me. One mapping point is established: De Nardi–Yang's 6.90 is wealth to after-tax earnings, while the project's 6.146 is to gross earnings. Do not compare them.

**5. Actions.**

Tonight:
1. Write the Section 2 convention into the run contract and launch the two arms under identical earnings, grids, entry law, targets and weights. Pre-commit the decision rule: full tables and diagnostics, not the scalar. Prefer 15 income states if the objective time fits; otherwise 7 states plus a 15-state check of both selected points, since 7 versus 15 moved first-birth age half a year and timing is the largest loss share.
2. Report h_P's bound status. If the floor arm sits at 2.3 again, the comparison is bound-versus-free, not floor-versus-share.

Tomorrow:
3. Rebuild the entry joint distribution as in Section 2.
4. Confirm the earnings variable's coverage, decide the old-age denominator, and get the author's bequest-tolerance decision.
5. Settle the three housing-by-children rows. Those are where the arms will differ most, and their empirical definitions are open, so tonight's loss ranking can flip when they are re-measured.

Disagreements with the lead. Do not harmonize the entry denominator by rescaling tonight; only the joint rebuild fixes denominator and coupling together, so run B as-is and label it. Retaining frozen targets experimentally is not enough; pre-commit a non-loss decision rule because three deciding rows are undefined. Timing is not the cause of the housing misses; do not spend tonight on it.

Sources: [De Nardi and Yang, NBER WP 20058](https://www.nber.org/system/files/working_papers/w20058/w20058.pdf); [Gale and Scholz 1994, AEA page, not retrievable](https://www.aeaweb.org/articles?id=10.1257%2Fjep.8.4.145).
