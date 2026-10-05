# Household mechanism slides

Five completed frames are inserted immediately before `Policy` in the continuing
`latex/JMP_slides/JMP_slides.tex`. Their shared editable source is
`latex/JMP_slides/mechanisms/frames.tex`; the standalone five-frame source is
`latex/JMP_slides/mechanisms/mechanism_excerpt.tex`. The reader outputs are
`output/pdf/JMP_mechanisms.pdf` and `output/pdf/JMP_slides.pdf`.

## Reference and scope

All five frames use the **same post-interest, soft-financing chain-13 point,
without Estate A**, at fixed baseline price 0.7760569760205563 (except the
explicit price experiment). This is the working old-wealth-target reference,
not a new calibration. Baseline arrays SHA-256:
`04a1d60ad526402fb0c2a52c853b77795b1b1d04488bb222684f0f1ed0d43100`.
The saved case is
`output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/credit_mechanism/credit_relaxation/phi_080`.

There were **zero new model solves** for this slide packet. The experiments
below were already saved and are reused without recalibration or adoption.
Neither frozen September 14 sources/results nor the author manuscript/mock
were changed. The pre-existing calibration, historical transition and policy
frames retain their earlier sources; adding these frames does not reconcile
those historical numerical results with chain 13.

The retained new 2023 state comes instead from the **one-birth Estate-A
experiment at chain-13 coordinates**, with a permanent 2007 preference shock.
It is not the subsequently recalibrated Estate-A candidate. It is not used in
these figures. Its terminal checks remain false; see
`output/model/transition_readiness_v1/current_baseline_20261003/retained_one_shock_v1/README.md`.

The [complete working-reference target fit](../evidence/baseline/target_fit.csv)
and [parameter/restriction table](../evidence/baseline/parameters.csv) remain
available. Early fertility is below its target, and the model's increasing
fertility-income gradient is not validated against matched prebirth resources.
These weaknesses limit quantitative policy interpretation; the figures explain
the behavior of this parameter point rather than certify a paper benchmark.

## Evidence map

| Frame | Figure and source | Economic comparison |
|---|---|---|
| Household Choices at First Birth | [household](household/README.md), `household_policies.pdf` | Initially childless renters, ages 22–25, fifth of nine earnings states; success versus waiting at identical beginning wealth. Both branches average actual tenure choices and native post-transaction controls. |
| The Parenthood Space Requirement | [floor](floor/README.md), `main_floor_state.pdf` | Same household states; only the parenthood floor changes from 2.593759507 to 2.781445432 rooms (+7.236%). Original case: `../price_decomposition/cases/arm_floor_only`. |
| Housing Costs and Fertility | [responses](responses/README.md), `first_birth_price110_by_beginning_tenure.pdf` | House price and implied rent both +10%; first-birth probabilities averaged using baseline prebirth exposure within each age/beginning-tenure group. |
| Mortgage Finance and Fertility | [responses](responses/README.md), `credit_phi095_wait_success_gains_debt_states.pdf` | Financed share 0.8→0.95; optimized branch-value gains at the same states with nonpositive beginning financial assets, grouped by current earnings. |
| Rental Space and Homeownership | [responses](responses/README.md), `renter_room_cap_births_ownership.pdf` | Only maximum rental rooms change (4, 6, 20); fixed-price lifecycle birth flows and young ownership, relative to cap 6. |

The floor and cap cases hold earnings, entrant wealth/income distributions,
timing, transfers, other preferences, house prices/rents, finance and targets
fixed. The higher-floor case exceeds the calibration search upper bound 2.6;
it is an explicit diagnostic perturbation, not an admissible calibration
candidate. The floor is a parenthood-only jump, not a per-child requirement.
The rental cap is six physical rooms in the reference, distinct from that
search bound. The clean price comparison changes prices, collateral values,
and incumbent housing wealth together; it does not identify rents alone.

## Definitions and checks

- The illustrative household has four-year after-tax income 2.061181307.
  Figure wealth units divide by this household's current annual gross earnings,
  using the same baseline denominator across branches and floor cases.
  The middle earnings **state** is not asserted to be the population median.
- Its 28 positive-exposure wealth nodes are audited; the 12 plotted nodes
  contain 99.668% of that age/earnings-state's prebirth exposure. The selected
  age/earnings state is 6.295% of fertile childless-renter exposure. This is an
  illustration, not an aggregate policy function.
- Success means the post-outcome state with one child ever born and one at home;
  waiting has neither. Constant child maturation makes these the correct
  saved branches. A failed first-birth attempt shares the no-child housing
  menu. Saving is ending net financial wealth; it includes the effect of
  purchasing a house and should not be read as total wealth accumulation.
- Housing is **physical rooms**, including the tenure mixture. Conditional
  renter room demand is never labeled realized housing. Native two-node
  transaction interpolation and separate stayer controls are used. Positive
  choice budget residuals are at most 5.4e-15; no infeasible positive exposure
  is removed. Baseline floor-extraction controls match the independent
  household extraction to 3e-14 or better.
- Price curves keep baseline prebirth weights fixed. In the matched price
  case, the first-birth policy term reconstructs exactly from the plotted
  renter/owner age cells. This is a household-probability comparison, not a
  national birth forecast or a market-clearing transition.
- Credit values are cardinal lifetime utility, not consumption-equivalent
  welfare. Positive-choice logit inversion covers 92.509% of the low-earnings,
  nonpositive-assets group and all of the corresponding middle/high groups;
  values condition on that support. The omitted value support is retained in
  all birth-flow accounting. Credit is valuable in both branches; the
  child-relative gain is heterogeneous. The slide's ownership/birth totals
  describe the full lifecycle cross-section, not just the graphed asset group.
- Cap panels intentionally use separate units: percent changes in birth flows
  and percentage-point changes in ownership, ages 18–29. Twenty rooms is an
  expanded menu, not a proof that every rental constraint is absent.

Subfolder receipts pin source files and record the detailed numerical checks.
The original diagnostics' incomplete canonical dated-budget/purchase audits
are not reclassified as passed by these local extraction checks.

## Established and outstanding

The numerical comparisons establish a parenthood space cost, a housing-price
response, and substantial substitution between rental and owned space. At this
point, increasing mortgage finance changes tenure much more than births.
The branch comparisons explain why a positive value of borrowing does not
guarantee a positive first-birth incentive. They do not establish welfare,
optimal policy, a universal credit null, or that correcting the income gradient
would restore a large mortgage-fertility response.

**Unconstrained credit is outstanding.** Financed share 1 still restricts debt
to collateral and leaves renter and terminal debt rules. A proper benchmark
requires an explicit repayment-at-death contract. Recommended definition:
remove collateral/origination and renter limits while retaining a natural
solvency bound under the maintained risk structure. Actuarially insured debt
would change the mortality-credit market too. No arbitrary negative grid
endpoint was substituted for this choice.

**Fixed physical housing stock is outstanding.** The retained transition fixes
the supply coefficient H0, not housing quantity. A valid comparison keeps the
same actual 2007 distribution, owned houses and both queues, the same preference
shock, fiscal and entry rules, and freezes actual H2007. Prices must clear the
dated housing market. Vacancy, depreciation replacement and the terminal
continuation require an explicit consistent contract. Failure of a stationary
endpoint does not prove a dated path infeasible. Starting from inherited 2023
instead would identify prospective adjustment, not the same 2007 experiment.
No synthetic fixed-stock result or unconstrained-credit placeholder is in the
main slides.

## Reproduction

From the repository root:

```sh
bash output/model/credit_mechanism_20261004/mechanism_slides/regenerate.sh
```

This reads saved data, regenerates figures, copies the five selected PDFs into
the continuing deck's mechanism assets, and compiles both PDFs twice from the
existing LaTeX source. It performs no model solves. Large saved source arrays
remain local; the plotting scripts, tidy CSVs, receipts and slide PDFs are the
compact reproducibility packet.
