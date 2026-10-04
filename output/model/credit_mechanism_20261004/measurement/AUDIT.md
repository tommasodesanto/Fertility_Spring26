# Fertility by income: measurement audit

**Decision.** The recovered numbers are correct under their stated code, but the comparison is not a matched fertility–income moment. The model’s sharply positive gradient survives reasonable grouping, age and family-unit checks. The 2024 CPS young-age gradient is sharply negative; the all-women completed-fertility gradient is approximately flat, with its small negative sign uncertain. This is a validation concern, not an identified causal explanation of weak credit effects or an adopted calibration target.

Verified October 4, 2026. One cached-only tabulation; no solve, recalibration, target edit, download or historical CPS rebuild. [Reproduction script](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/credit_mechanism_20261004/measurement/recompute.py), [complete tables and identities](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/credit_mechanism_20261004/measurement/tables.json).

## Exact provenance

The authoritative builder for the pasted comparison is [the preserved scratch script](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/credit_mechanism_20261004/evidence/preserved/fertility_by_income.py), with [its original result](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/credit_mechanism_20261004/evidence/preserved/fertility_by_income.json). It selects national June 2024 women (`PESEX=2`), valid `PTSF1` 0–5, positive `PWSSWGT`, ages 24–26 or 40–44, and estimates weighted means of \(\min(N,3)\), where \(N\) is children ever born. No fixed effects, regression, clustering or uncertainty appeared in the original result. All twelve pasted CPS/model means reproduce within \(2.3\times10^{-16}\) using the retained fixed-price chain-13 cache.

CPS groups are **fixed brackets**, not terciles: `HEFAMINC` 1–10 below $40,000; 11–14 $40,000–99,999; 15–16 $100,000+. Their older-woman population shares are 18.2%, 34.3%, 47.5%; younger shares are 26.7%, 39.8%, 33.5%. `HEFAMINC` measures combined family money income over the previous twelve months, including earnings, business/rental income, pensions, investment income and transfers. `PTSF1=5` denotes five-plus births; `PWSSWGT` is the basic CPS final person weight. The canonical README’s “supplement weight” wording is imprecise. [Census documentation, pp. 3-3, 6-4/5, 7-1](https://www2.census.gov/programs-surveys/cps/techdocs/cpsjun24.pdf).

The model groups **current persistent Markov earnings states** \(z\), not permanent types (`permanent_income_levels_enabled=false`). Cumulative-mass midpoints assign states 1–4 / 5 / 6–9 to three groups with shares 36.3% / 27.3% / 36.3%. At one working age, gross earnings and \(z\) have identical ranks; that does not equate their economic content with family money income. The scratch script solved at the chain-13 price but saved only summaries, without its own pinned inputs/arrays. The retained phi=.8 cache reproduces those summaries exactly; this verifies the tabulation, not every original execution detail.

`g_beginning_distribution` contains current birth flows before tenure transactions: the forward map moves births into the current child-count cell before saving that array. Consequently, j=1 is the **end** of [22,26), and j=6 is the **end** of [42,46), rather than an average interview-age stock. Both script sides cap children at three; the model’s 3+ weight 3.602 is deliberately excluded from this comparison.

## Bounded sensitivity checks

Each triple is low / middle / high mean children ever born capped at three.

| Comparison | CPS | Model | High minus low: CPS / model |
|---|---|---|---|
| Original young | .635 / .447 / .226 | .079 / .529 / 1.096 | −.409 / +1.017 |
| Equal-weight thirds; uniform ages 24–26 | .596 / .448 / .226 | .124 / .508 / .973 | −.369 / +.848 |
| Original older | 1.771 / 1.765 / 1.710 | 1.477 / 1.918 / 2.220 | −.061 / +.743 |
| Equal-weight thirds; uniform ages 40–44 | 1.784 / 1.702 / 1.733 | 1.318 / 1.782 / 2.104 | −.051 / +.786 |

Equal-weight thirds fractionally split tied CPS income categories or model states, applying the same fraction to every observation within a tie. Model age exposure uniformly mixes exact reconstructed prebirth and postbirth stocks within each four-year cell. Native flow and forward-map reconstruction errors are below \(4\times10^{-17}\). These are transparent diagnostic approximations, not a new empirical observer contract.

The original CPS younger sample has 1,673 women; the older sample 3,278. Approximate household-cluster sandwich SEs for the fixed-bracket high-minus-low contrasts are .058 young and .060 older; equal-third SEs are .052 and .050. The older fixed-bracket approximate 95% interval is [−.178,.057]. These condition on original weights and omit CPS PSU/strata/replicate-weight corrections; they cannot establish survey-design uncertainty or SMM weights.

| CPS selection sensitivity, fixed brackets | Cap-three means | High minus low |
|---|---|---:|
| Ages 24–26, family reference person/spouse | 1.235 / .896 / .535 | −.700 |
| Ages 40–44, family reference person/spouse | 2.081 / 2.036 / 1.873 | −.208 |
| Ages 40–44, all women, public cap-five | 1.999 / 1.987 / 1.838 | −.161 |

The reference/spouse restriction removes many young women living in another family and changes fertility levels materially. It is a selection sensitivity, not proof of exact household equivalence. Older all-women childlessness falls from .200 to .183 with income, while reference/spouse childlessness rises from .080 to .112: unit selection can reverse this particular margin. Capping at three removes more births from lower-income older women and flattens their count gradient.

## What remains unidentified

Current family income can reflect marriage, coresidence, the partner’s earnings and fertility-related labor supply. Thus the CPS gradient does not identify the effect of exogenous permanent earnings or the model’s child-cost primitive. Selection and fertility–earnings endogeneity prevent using its numerical gap as a causal explanation of the credit result. The robust finding is the model’s positive sorting versus the negative young-age count gradient across all tested CPS selections; the older all-women sign is not statistically established.

The live fertility calibration uses June **2004/2006**, not 2024. A narrow check of the retained 35/34-MB June partitions verifies their authoritative hashes but finds `HHINCOME` and `FTOTVAL` blank for all 15,999 selected women; the extract lacks `HEFAMINC`. Therefore the baseline-vintage income gradient remains **unavailable**, not zero or flat. [Support receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/credit_mechanism_20261004/measurement/historical_income_support.json). No giant source rebuild was attempted. Matching the baseline vintage and a household earnings definition is the next empirical requirement before changing a target or preference specification.
