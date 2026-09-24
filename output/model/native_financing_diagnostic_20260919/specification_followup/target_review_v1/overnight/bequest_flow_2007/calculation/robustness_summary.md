# 2007 bequest-flow robustness

The adopted target remains **0.7291023472616158%** ($674.4908bn / $92.5096tn; trust-excluded positive estate base over aggregate signed net worth). Rows below change one saved scenario dimension at a time unless identified otherwise. Dollar totals are 2022 USD. Percentage-point differences are relative to the adopted ratio.

| Scenario | Annual flow ($bn) | Denominator ($tn) | Ratio | Change (pp) | Interpretation |
|---|---:|---:|---:|---:|---|
| adopted_baseline | 674.490 | 92.5096 | 0.729102% | +0.000000 | Adopted child-directed construction |
| spouse_death_equal_split_child_flow | 612.823 | 92.5096 | 0.662442% | -0.066660 | Same child-directed estimand; spouse-death/wealth-allocation alternative |
| cohabiting_partner_treated_as_married | 646.181 | 92.5096 | 0.698502% | -0.030600 | Same child-directed estimand; partner-recipient mapping alternative |
| topcoded_age_maps_to_qx99 | 679.803 | 92.5096 | 0.734846% | +0.005744 | Same child-directed estimand; 95+ mortality mapping |
| topcoded_age_maps_to_qx100plus | 723.974 | 92.5096 | 0.782593% | +0.053491 | Same child-directed estimand; terminal mortality mapping |
| trusts_included_in_estate_base | 695.099 | 92.5096 | 0.751380% | +0.022278 | Child-directed measurement check |
| trusts_excluded_from_denominator | 674.490 | 91.0732 | 0.740602% | +0.011500 | Child-directed denominator measurement check |
| whole_estate_not_child_directed | 1,385.032 | 92.5096 | 1.497177% | +0.768075 | Different estimand: gross all-estate transfer |

The quantified same-child-directed scenarios retaining trust-excluded estate wealth and aggregate-net-worth denominator span **0.662442%–0.734846%** when q95/q99, spouse allocation, and partner mapping alternatives are considered. This is a scenario envelope, not a confidence interval. The terminal q100+ mapping raises the result to 0.782593%, an extreme top-code sensitivity because NCHS sets terminal qx to 1.0. Among quantified alternatives, independent spouse deaths plus equal estate splitting move the estimate down 0.066660 percentage points; partner-as-married moves it down 0.030600 points. The q95-to-q99 map moves it up 0.005744 points.

The all-estate row is a different estimand: removing the child shares raises the spouse-death proxy to 1.497177%; do not call it a robustness estimate of the adopted child-directed target.

The published source documents the 25%/75% child shares, trust-excluded bequest base, and aggregate-net-worth denominator. It does not give the household estate ownership rule at single-spouse versus joint death. The equal-split alternative assumes independent spouse deaths and is not source-prescribed. SCF cohabiting partners receive 75% under the baseline legal-marriage mapping; treating them as married (25%) is a measurement alternative. Public ages are top-coded at 95; q95 and q99 are reported, with q100+ isolated as an extreme mechanical mapping.

Lifetime-child eligibility is not quantified: `KIDS` counts current rostered children, not all living offspring or lifetime child status. The adopted calculation mechanically applies the child shares to all PEUs; no valid alternative eligibility measure exists in the saved inputs. No statistical uncertainty interval is implied.

Exact scenario numerators, denominators, formulas and source assumptions are in `robustness_summary.json` and `robustness_summary.csv`; the underlying saved 72-scenario grid and adoption receipt remain unchanged.
