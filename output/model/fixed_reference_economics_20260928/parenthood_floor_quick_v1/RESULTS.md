# Fixed-parameter utility diagnostics: terminal readout

**No new arm obtained a birth-renewal equilibrium in the tested price interval.** All three second-attempt jobs (Torch18902151), lasting 2m11s, 2m19s and 2m15s respectively and each completing three native lifecycle solves, failed with `Renewal root unbracketed on [.85,1.15] qref` after three successful prescribed-price evaluations. No selected equilibrium, independent equilibrium repeat, target-fit table, parameter-fit table or standard 17-plot packet was generated for any new arm. This is a bracket failure; it does not establish absence of equilibrium, household infeasibility, or poor calibrated fit.

These are fixed-parameter mechanism experiments, not estimation. All arms use the authenticated nonnegative-wealth pilot-selected guesses, 120 assets × 9 income states, zero renter borrowing, 2% annual interest, fixed child-benefit scale and housing supply scale, and the unchanged target/weight contract. Price must clear actual birth renewal; population clears physical housing supply. Compensation is disabled in all three. The floor arm additionally fixes the child expenditure share and imposes a 1.890-room parenthood-only physical requirement; no_A retains the inherited first-child expenditure-share change with no floor; constant_alpha fixes the child expenditure share with no floor. Nonlinear child benefits and corrected owner/mortality rules are retained. No production adoption follows.

The renewal residual is birth-implied entrant flow divided by the required entrant flow, minus one. Thus a negative number means births cannot replenish the required entrants at that prescribed price. Birth-implied entry equals adjusted births per normalized household times the fixed conversion factor 0.476. These are closure diagnostics at prescribed prices, not target-fit substitutes.

| Arm | Point | Price | Renewal residual | Adjusted births / household | Birth-implied entry | Required entry |
|---|---|---:|---:|---:|---:|---:|
| floor | lower_085 | 0.671 | -0.220 | 0.101 | 0.048 | 0.062 |
| floor | qref | 0.790 | -0.302 | 0.090 | 0.043 | 0.062 |
| floor | upper_115 | 0.908 | -0.377 | 0.081 | 0.038 | 0.062 |
| no_A | lower_085 | 0.671 | 0.410 | 0.183 | 0.087 | 0.062 |
| no_A | qref | 0.790 | 0.390 | 0.180 | 0.086 | 0.062 |
| no_A | upper_115 | 0.908 | 0.368 | 0.177 | 0.084 | 0.062 |
| constant_alpha | lower_085 | 0.671 | 0.141 | 0.148 | 0.070 | 0.062 |
| constant_alpha | qref | 0.790 | 0.115 | 0.145 | 0.069 | 0.062 |
| constant_alpha | upper_115 | 0.908 | 0.092 | 0.142 | 0.067 | 0.062 |
| authenticated_control | selected | 0.799 | -7.914e-10 | 0.130 | 0.062 | 0.062 |

All nine prescribed-price points have zero recorded absolute housing residual and zero occupied-renter floor violation mass. The floor arm has negative renewal residuals throughout the interval; no_A and constant_alpha have positive residuals throughout. Over these three tested prices, raising price reduces the renewal residual in every arm. The observed direction therefore suggests testing lower prices for the floor arm and higher prices for no_A and constant_alpha; this is a direction for a new check, not a guarantee of a root. A separately authorized broader price-bracket check would be needed to locate a renewal root before evaluating equilibrium fit or estimating parameters.

[Full-precision price/flow diagnostics](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/parenthood_floor_quick_v1/collected/price_endpoint_comparison.csv) · [Full 14-target comparison, with new equilibrium moments explicitly unavailable](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/parenthood_floor_quick_v1/collected/target_comparison.csv) · [All 31 observed parameter inputs/restrictions/bounds per arm](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/parenthood_floor_quick_v1/collected/parameter_comparison.csv)

**Any CSVs under `preflight/ge_loop_mock/` are synthetic mock-loop artifacts and are not economic results.** The real new runs produced prescribed-price closure/gate JSONs only. Their target/model/gap/loss cells remain unavailable in the comparison CSV. No parameter was estimated in these new runs. Historical near-bound flags for both fertility dispersion parameters are carried from the fixed control and do not imply a new estimate. The reported actual-input dictionary, rather than its inherited `free` label, establishes the floor and constant-share overrides.

## Reused authenticated control

The control is the previously verified nonnegative-wealth pilot-selected equilibrium, loss 19.697. Its complete14-target and31-parameter tables and17PNG hashes are authenticated; both independent control repeats have byte-identical CSV/closure files. It was not rerun in this diagnostic.

[Control target fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/parenthood_floor_quick_v1/reference_control/target_fit.csv) · [Control parameters/bounds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/parenthood_floor_quick_v1/reference_control/parameters.csv) · [Existing control standard plots](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/nonnegative_mean/011_gn_0.5/phase_b_ge/selected_root/standard_diagnostics)

![Previously verified control fertility profile](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/collected/nonnegative_mean/011_gn_0.5/phase_b_ge/selected_root/standard_diagnostics/fertility_by_age.png)

## Failure and collection evidence

First attempt Torch18901961 failed before any lifecycle solve because repeated in-process native authentication requires a fresh interpreter. Its failure receipts, tracebacks and successful zero-solve preflights are preserved in collected/attempt1. The bounded orchestration correction retained the original deadline; second-attempt failures and three price-point closures per arm are preserved in collected/attempt2. No further retry or bracket expansion was performed.

All60 retained first-attempt and114 second-attempt downloaded files have independently matching remoteSHA256 hashes. No large model arrays or source snapshots are retained. New standard PNGs are unavailable because reporting occurs only after a selected renewal root. Existing twelve calibration searches are outside this collection.
