# October 1 housing supply level check

Read-only empirical calculations; no model imports, solves, simulations, recalibration or project edits. This is a level diagnostic, not a production target change.

Verdict: the adopted reference passes a useful preliminary national2007 quantity/cost plausibility check. Its implied rent is2.59% below the AHS contract-rent-per-room stock ratio and purchase price0.31% above the AHS self-reported-value-per-room stock ratio. The October1 floor experimental point instead has rent11.31% below and value8.67% below. These are point comparisons without joint statistical inference.

| Cost per physical room | AHS2007 | Adopted block0506 | Experimental floor chain7 0173 |
|---|---:|---:|---:|
| Monthly contract rent (dollars) | 168.22 | 163.86 | 149.20 |
| Purchase value (dollars) | 43495.49 | 43629.03 | 39723.80 |

The room target is5.729434240102641, literal national AHS2007 occupied rooms for heads18--85, public topcode21. AHS room statistics exactly reproduce the prior target record:37793 records, weighted107194201.15842237 households, weighted614162126.4575154 rooms.

Primary empirical unit prices are ratio of weighted annual contract rent (or owner VALUE) to weighted rooms, within the corresponding tenure sample. These are proxies for the model common price, not hedonic quality-adjusted prices or imputed owner rent. Rent is RENT times FRENT (annual payments); RENT=1 income-dependent records and FRENT=53 topcoded frequency are excluded. Missing/nonpositive rent/value are excluded; no arbitrary winsorization. Excluding reported public/subsidized/voucher/controlled units changes rent to172.45 dollars per room/month.

Money conversion: the model wage is1 and its profile has equal-working-age-cell average1 (12four-year cells,18--65). PSID2007 reference persons, IW>0, nonnegative RP/spouse gross EARNINDRRC produce4981 observations. IW-weighted means within those12cells, equally averaged, give77963.06267809066 real2022 dollars. CPI-U207.342/292.655 gives55235.7463286145 dollars in2007 purchasing power per model monetary unit. Period rent is divided by4 to make annual rent; purchase price is a stock and is not divided by4.

Using the observed room quantity and rent together, the externally anchored coefficient is H0=6.06484 in the retained rbar=.16 form, or A=19.24082 in the absorbed form. Current values are6.29351 and19.96628: the external anchor is3.63% lower. Excluding reported subsidies/control gives H0=5.97067,5.13% lower. These are alternatives computed with elasticity.63 unchanged, not recalibrated solutions.

Closure judgment: in the experimental birth-renewal-price/endogenous-population model, H0 scales population and cannot be identified by per-household moments. Recommend an externally measured stock/cost/population anchor for H0 there, with mean rooms retained to discipline household demand and observed cost checked against the renewal-determined price. Changing H0 alone cannot correct the experimental low price. With normalized fixed population and price clearing housing, joint estimation of H0 from room quantities is legitimate; the new cost diagnostic supports its present level approximately.

Limits: PSID interview2007 labor income refers to a tax year (primarily2006), whereas AHS expenses refer to interview/current payments and AHS incomes past12months. Purchasing-power conversion is explicit, but exact production timing/deflator construction needs harmonization. PSID across-working-head instead of equalagecell averaging gives59339.10 nominal dollars and a diagnostic H0=6.34491; that is a weighting sensitivity, not the same model money normalization. Owner value is self-reported; owner/renter quality differs; utilities may be bundled with contract rent. No claims of statistical equivalence, elasticity validation, adopted target changes or equilibrium robustness.

Sources:
- https://www2.census.gov/programs-surveys/ahs/2007/
- https://www.census.gov/data-tools/demo/data/uccb/ahsdict/minicodebooks/ahs_mini_2007National.pdf
- https://www.bls.gov/regions/northeast/data/consumerpriceindex_us_table.htm
- Local retained PSID shelf, EARNINDRRC label: RP/SP combined earnings (realUSD2022), tax year.

All raw files, calculate.py, psid2007_earnings_by_age.csv, ahs_price_quantities.json and anchor_comparison.json are in this temporary packet. Model references and full calibration tables remain in the Fertility project.

Full saved calibration tables (no selected-row fit comparison):
- [output/model/fertility_identification_20260928/resume_v1/selected_export/primary/target_fit.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/resume_v1/selected_export/primary/target_fit.csv)
- [output/model/fertility_identification_20260928/resume_v1/selected_export/primary/parameters.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/resume_v1/selected_export/primary/parameters.csv)
- [output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/target_fit.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/target_fit.csv)
- [output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/parameters.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/parameters.csv)

Public AHS ZIP files remain in /tmp/fertility_ahs2007_supply_check; source URLs and hash are in the JSON record. The PSID source remains in its existing shelf location.
