# Saved-case mechanism figures

Three 11 × 4 inch figures for the mechanism slides, built without a model solve. Each has a PDF, PNG, and plotted CSV in this folder. Regenerate with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  code/model/.venv/bin/python \
  output/model/credit_mechanism_20261004/mechanism_slides/responses/build_mechanism_figures.py
```

The script checks the saved case kind, matched baseline source hashes, price ratio, within-cell flow identity, branch-value identity, and renter-cap field audits. Exact array hashes and plotted row counts are in `plot_receipt.json`. It reads only existing diagnostic CSVs and receipts.

## First-birth response to housing cost

`first_birth_price110_by_beginning_tenure.pdf` plots the four-year first-birth probability by age and **beginning** tenure. At each age and tenure, the baseline prebirth mass (G_0(s)) weights both the baseline probability (h_0(s)=\pi_j a_0(s)) and the price-case probability (h_1(s)=\pi_j a_1(s)). The source is the actual `price110` matched-case `diagnostics/price110_decomposition/common_state_groups.csv`, whose receipt identifies baseline price 0.7760569760 and imposed price 0.8536626736. Thus the plotted difference is a common-state policy response, not a change in the population composition. Across plotted first-birth states, (\sum_sG_0(s)[h_1(s)-h_0(s)]=-0.0074670508), matching the saved first-birth policy term. Beginning renters show a larger percentage-point decline than beginning owners at most ages.

The imposed price change also raises rents and purchase costs and changes collateral and incumbent housing wealth. This is a fixed-price household diagnostic, not a general-equilibrium price change. It does not isolate the parent room requirement.

## Credit and the value of a successful first birth

`credit_phi095_wait_success_gains_debt_states.pdf` uses `diagnostics/phi095_decomposition/first_birth_income_wealth_summary.csv`, the receipt explicitly labelled `pair_changes=["phi"]`. The comparison holds the price fixed and changes the financed share from 0.8 to 0.95. For households with beginning financial assets (b\leq0), the bars are changes in optimized lifetime utility on the wait and successful-first-birth branches. They are baseline-prebirth-exposure-weighted means **conditional on exact logit-inversion support**, pooled across fertile ages and beginning tenure. The low/middle/high cuts in annual gross current earnings are 0.338391 and 0.732667. These are utility units, not income equivalents or welfare comparisons.

For the lowest earnings group, the wait branch gains 0.01836 utility and the success branch 0.01430: the successful-birth-minus-wait gap falls by 0.00406. The corresponding gap rises by 0.000654 for the middle and 0.001419 for the high group. Logit support is 92.509% of low-group exposure and 100% in the other two displayed groups; all exposure remains in the underlying birth-flow accounting. The chart establishes heterogeneous relative value responses in this asset range. It does not attribute the full stationary fertility response to these cells or show that expanding credit raises aggregate births.

## Renter space and tenure

`renter_room_cap_births_ownership.pdf` reads `price_decomposition/lead_followup/cap/results.csv`. With all other recorded inputs unchanged, reducing the maximum renter unit from six to four rooms changes explicit births by −2.144% and ages-18–29 ownership by +16.359 percentage points. Raising the cap from six to twenty rooms changes births by +0.367% and ownership by −11.269 percentage points. The figure has separate axes for percent birth changes and ownership percentage-point changes. Both arms have completed fixed-price status and an audit reporting only `hR_max` changed. The age-18–29 ownership statistic is a cross-sectional rate; births are stationary household flows. This is a cap sensitivity, not a direct estimate of household willingness to pay for rooms.

## Related price-channel evidence

The newer `price_decomposition/lead_followup/README.md` finds an order-free decomposition of an approximately equivalent price change: the money group contributes −2.973 percentage points to births, the parent room requirement −3.780 points, and other room inputs +0.268 points, summing to −6.485%. The actual matched price case changes births by −6.434%; the 0.051-point discrepancy is an available grid/equivalence error measure, and no finer-grid check was run. The parent requirement (h_P=2.594) is close to its **2.6 search bound**. It is a different object from the six-room renter cap. These qualifications make the actual matched price case the primary price figure; the decomposition is supporting evidence, not an exact causal allocation.

The source receipts also distinguish these household diagnostics from an accepted equilibrium policy comparison. No dated-budget/purchase audit or general-equilibrium acceptance is claimed here.
