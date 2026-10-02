# Independent assessment of purchase financing and the calibration plateau — October 1, 2026

Read-only evidence packet behind [ASSESSMENT.md](ASSESSMENT.md). Nothing here solves the
model, changes a target, or adopts a specification. Every number is arithmetic on outputs
that already existed.

## Facts (source-linked)

| Object | Verified fact | Evidence |
|---|---|---|
| Points used | (a) 1,277 passed, **unverified** search candidates of `normalized_calibration_v2` (best sampled loss 30.0778, chain 2 / `0053_nm`); (b) saved arrays of the earlier **winner31** point (chain 7 / `0173_nm`, loss 31.284, price 0.719168, `h_P` = 2.3), which is not the latest calibration | `../normalized_calibration_v2/deployment/calibration_plateau_diagnosis_v1/`; `../utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/purchase_ltv_v1/local_run/retry5/results/` |
| Units | One model unit = equal-working-age mean annual gross earnings ($55,236 in 2007 dollars); one period = four years; `R = 1.02^4` | `../README.md`; `../room_unit_audit_v1/price_level_readout.md` |
| Executed engine | `code/model/refactor_lab/engine/` (byte-identical to the pinned `small_credit_lab` copy); purchase rules are native flag-gated code, not runtime-rewritten text | `../normalized_calibration_v2/source_pins.json` |
| Buyer budget | `c + (delta+tau_H) Q + b' = R (b - Q) + y`, `b' >= -phi Q`; screen `b + y/R >= (1-phi) Q` | `code/model/refactor_lab/engine/household.py:296-299, 423-437, 894-901`; `kernels.py:403-411, 885-913` |
| Screen redundancy | The screen is implied by the floor and positive consumption; slack is at least `0.1513 Q` | algebra in `ASSESSMENT.md` §1; `02_purchase_tests_saved_arrays.out.txt` |
| Origination LTV in the model (winner31) | Purchases by never-parent renters aged 18–42: 71.6% close above 80% LTV, 53.6% above 90%, 31.0% at 100%; 7.4% end the purchase period at the 80% limit | `02_purchase_tests_saved_arrays.out.txt` |
| Alternative purchase tests (winner31) | Share of those purchases that would fail: stock-only test 71.6%; one-year-income test 8.7%; implemented screen 0% | same |
| Financing arms (winner31, fixed price) | Purchase LTV 80% → 100%: `early_fertility` 0.5312 → 0.5291, `first_birth_rooms` 1.2232 → 1.2295, `mean_rooms` 5.977 → 5.976–6.031; `recent_parent_ownership` 0.118 → 0.005 | `04_financing_arms_target_table.out.txt` |
| Early-fertility ceiling (winner31) | Observer reproduced (0.5313); with every 18–21 mother having a second birth at 22–25 the observer is 0.674; target 0.8095 | `03_early_fertility_bound_and_selection.out.txt` |
| Who has first births (winner31) | First-birth probability per period at ages 18–25 is 0.000–0.011 for income states ≤ 0.47 (40% of childless at risk) and 0.24–0.78 above | same |
| Rental cap (winner31) | 56.6% of new parents rent after the birth; 70.6% of them are at the six-room cap. Ownership of the same households is 43.4% with the birth and 46.6% without | same |
| Search cloud | Regression Jacobian from 1,277 candidates (out-of-sample correlation of predicted and actual loss 0.98); loss gradient at the best point is not zero; `h_P` at its bound with the loss falling outward; singular values per 1% parameter change 5.64 … 0.0001 | `01_search_cloud_jacobian.out.txt` |
| Papers | Greaney, Parkhomenko and Van Nieuwerburgh (Feb. 16, 2025) eqs. 2.2–2.4 and Sommer, Sullivan and Verbrugge (2013) eqs. 3, 6–12 and §3.3 read directly from PDF | `/Users/tommasodesanto/Downloads/Dynamic_Urban_Economics.pdf`; `../utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/diagnosis/purchase_timing_review_v1/opus55/input/references/SSV2013.pdf` |

## Files

- `run_all.sh` regenerates the four `.out.txt` files from the repository root
  (`sh output/model/fixed_reference_economics_20260928/independent_assessment_20261001/run_all.sh`).
- `01_search_cloud_jacobian.py`: local Jacobian by regression on the search candidates; bounded
  linear least squares; Levenberg–Marquardt path; singular values; weight illustration.
- `02_purchase_tests_saved_arrays.py`: origination and end-of-period LTV of renter-buyers;
  budget-identity check of the array reading; pass rates of alternative purchase tests.
- `03_early_fertility_bound_and_selection.py`: ceiling of the age-25 children-ever-born
  observer; first births by income state; rental-cap incidence among new parents.
- `04_financing_arms_target_table.py`: the 14 target rows across the saved financing arms.

## Limits

- The search candidates are exploratory and not post-checked. The regression Jacobian is
  local to a neighbourhood a few thousandths wide in each parameter; extrapolations beyond
  it are labelled as such.
- The array tabulations are at the winner31 point, whose ten free parameters differ from the
  latest sampled best by at most a few percent. They are not results at the latest point,
  whose arrays are not retained locally.
- The purchase tabulations cover never-parent renters (including the first-birth branch),
  the group for which pre-fertility mass, attempt probability and branch tenure choice can
  be combined exactly from saved arrays.
- No alternative purchase rule was solved. Pass rates are computed on the baseline
  distribution and say how much a rule would bind on impact, not what the new equilibrium is.
