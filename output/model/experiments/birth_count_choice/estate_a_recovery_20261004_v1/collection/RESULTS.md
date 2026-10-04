# Estate-A recovery result, October 4, 2026

Fifteen isolated chains ran under the experimental Estate-A/new-wealth contract with unchanged economics, 14 reported moments (ten scored), weights, and ten free-parameter bounds. Eleven passed the full fresh native selected-point and exact-repeat gate; four count-three chains failed with `native GE acceptance failed: uncomputed_bounded_budget`. The latter have only provisional saved best cases. No optimization-convergence certificate, grid certification, rank/identification certification, or scientific adoption is claimed.

Best **binary** chain 1: verified loss **21.275413361** after 137 calls. Best **count-three** chain 3: verified loss **78.859933021** after 126 calls. Both searches stopped at the 800-second prelaunch guard and completed a fresh native child with exact-repeat diagnostics. A larger household choice set does not mathematically nest the old aggregate fit, so these losses alone are not a causal interpretation.

The table contains every target row. Blank weights and loss contributions are normalization rows; the exact strings are in [target_fit_joint.csv](target_fit_joint.csv).

| Moment | Role | Target | Weight | Binary model | Binary gap | Binary loss | Count-three model | Count-three gap | Count-three loss |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| initial_normalization | normalization | 2.1 | — | 2.1000001 | 1.3241635e-07 | — | 2.1 | 3.1514386e-08 | — |
| cps_childlessness | scored | 0.19827875 | 35532.304 | 0.19937083 | 0.0010920786 | 0.042377097 | 0.21938057 | 0.021101818 | 15.822064 |
| cps_exactly_one | scored | 0.21365533 | 26952.821 | 0.21983342 | 0.006178092 | 1.0287574 | 0.22074977 | 0.0070944488 | 1.3565679 |
| nchs_mean_age | scored | 25.976264 | 139.82807 | 25.973636 | -0.0026278925 | 0.0009656273 | 26.032172 | 0.055907929 | 0.43706011 |
| nchs_share30 | validation | 0.24927801 | 0 | 0.23830065 | -0.010977364 | 0 | 0.24412359 | -0.0051544243 | 0 |
| wealth_earnings | scored | 4.4583871 | 7.5950985 | 5.0291101 | 0.57072296 | 2.4739112 | 6.8866059 | 2.4282188 | 44.782573 |
| bequest_wealth | scored | 0.0072910235 | 5165289.3 | 0.0064609532 | -0.00083007031 | 3.5589706 | 0.0065242896 | -0.00076673384 | 3.0365743 |
| old_dispersion | validation | 3.5159351 | 0 | 3.2495809 | -0.26635424 | 0 | 2.8818535 | -0.63408162 | 0 |
| mean_rooms | scored | 5.7294342 | 128.0207 | 5.6513691 | -0.078065099 | 0.78017859 | 5.8069762 | 0.077541987 | 0.76975773 |
| ownership_30_55 | scored | 0.67626042 | 2339.3624 | 0.6794061 | 0.0031456813 | 0.023148718 | 0.65512818 | -0.021132241 | 1.0446929 |
| first_birth_rooms | scored | 1.465 | 137.56527 | 1.2462912 | -0.21870877 | 6.5802323 | 1.3421469 | -0.12285309 | 2.0762564 |
| family_rooms | validation | 0.38509965 | 0 | 0.30252636 | -0.082573288 | 0 | 0.21802565 | -0.167074 | 0 |
| recent_parent_ownership | scored | 0.12760836 | 27055.823 | 0.12865866 | 0.0010502954 | 0.029845833 | 0.12524263 | -0.0023657304 | 0.15142279 |
| early_fertility | scored | 0.80952764 | 100 | 0.54958484 | -0.2599428 | 6.757026 | 0.50321139 | -0.30631625 | 9.3829647 |

The next table contains all 31 reported parameters. The ten searched coordinates are exactly those in the recovery start-plan bounds; the other rows are fixed, derived, or retained inputs. `Near` is the native report flag. Full-precision values are in [parameters_joint.csv](parameters_joint.csv).

| Parameter | Status | Binary estimate | Binary bounds | Binary near | Count-three estimate | Count-three bounds | Count-three near |
|---|---|---:|---:|:---:|---:|---:|:---:|
| H0 | derived housing supply coefficient at N0=1; reference bounds advisory | 6.1245553 | [0.2, 80] | False | 6.3354776 | [0.2, 80] | False |
| beta_annual | estimated (free; recovery search bounds) | 0.9444438 | [0.93, 0.99] | False | 0.96756659 | [0.93, 0.99] | False |
| chi | estimated (free; recovery search bounds) | 1.088197 | [0.1, 5] | False | 1.0392911 | [0.1, 5] | False |
| first_birth_fixed_cost | estimated (free; recovery search bounds) | 0.28270677 | [0, 8] | False | 0.36201349 | [0, 8] | False |
| kappa_fert | estimated (free; recovery search bounds) | 0.11612031 | [0.02, 50] | True | 0.082600907 | [0.02, 50] | True |
| kappa_fert_continuation | estimated (free; recovery search bounds) | 0.41642923 | [0.02, 50] | True | 0.52860295 | [0.02, 50] | False |
| theta0 | estimated (free; recovery search bounds) | 0.14389227 | [0, 8] | False | 0.1187683 | [0, 8] | False |
| delta_alpha_jump | fixed zero under experimental utility contract | 0 | — | — | 0 | — | — |
| child_benefit_curvature | estimated (free; recovery search bounds) | 0.10313884 | [0, 0.8] | False | 0.067833317 | [0, 0.8] | False |
| tenure_choice_kappa | estimated (free; recovery search bounds) | 0.015769047 | [0.001, 0.1] | False | 0.016223367 | [0.001, 0.1] | False |
| psi_child | estimated (free; recovery search bounds) | 0.17788163 | [0.01, 0.5] | False | 0.10503548 | [0.01, 0.5] | False |
| child_benefit_CRRA_coefficient | derived from supplied benefit and curvature | 0.15953512 | — | — | 0.097910576 | — | — |
| theta1 | fixed external restriction | 0.0081930841 | — | — | 0.0081930841 | — | — |
| sigma | fixed | 2 | — | — | 2 | — | — |
| alpha_cons | fixed CEX childless expenditure share | 0.733 | — | — | 0.733 | — | — |
| delta_alpha | fixed zero under experimental utility contract | 0 | — | — | 0 | — | — |
| h_P | estimated (free; recovery search bounds) | 2.5426884 | [0.1, 2.6] | False | 2.5461512 | [0.1, 2.6] | False |
| utility_reference_rent | retained inactive normalization; compensation off | 0.11046593 | — | — | 0.11046593 | — | — |
| q_annual | author-retained 2% annual real rate | 0.02 | — | — | 0.02 | — | — |
| financed_share | inherited credit contract | 0.8 | — | — | 0.8 | — | — |
| housing_supply_elasticity | fixed provisional external mapping | 0.63 | — | — | 0.63 | — | — |
| payroll_tax | derived from adopted pension ratio | 0.08028071 | — | — | 0.08028071 | — | — |
| pension_period | balanced PAYGO | 0.91778405 | — | — | 0.91778405 | — | — |
| annual_depreciation | adopted | 0.014161437 | — | — | 0.014161437 | — | — |
| period_depreciation | compounded | 0.055453791 | — | — | 0.055453791 | — | — |
| annual_property_tax | adopted | 0.010598361 | — | — | 0.010598361 | — | — |
| period_property_tax | linear period convention | 0.042393443 | — | — | 0.042393443 | — | — |
| selling_cost | retained | 0.06 | — | — | 0.06 | — | — |
| rental_cap | retained provisional | 6 | — | — | 6 | — | — |
| wealth_grid_nodes | retained exact grid | 120 | — | — | 120 | — | — |
| income_states | retained B15 | 9 | — | — | 9 | — | — |

All chain outcomes are in [chain_outcomes.csv](chain_outcomes.csv). The four failed count-three chains (0, 1, 4, 6) retain provisional best losses 152.924831, 550.605309, 105.082178, and 484.766329, respectively; these values did not pass a final native postcheck and are excluded from winner selection.

Selected actual 14-row CSVs and 31-row parameter CSVs match their receipts. The selected and exact-repeat target rows match, and the 17 selected PNG hashes per arm match both the repeat PNGs and native receipt. The child search SHA and selected residuals also match. See [audit_receipt.json](audit_receipt.json) and the downloaded [binary](binary/selected_root/standard_diagnostics/fertility_by_age.png) and [count-three](count3/selected_root/standard_diagnostics/fertility_by_age.png) diagnostics. The remaining standard diagnostic PNGs are in adjacent directories.

Source: Slurm array 19141024; pinned recovery inventory `974a14e243da6a2ad0572bb9825b47ab349828f9144cbbad69e74b40a4408b22`; target fingerprint `c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70`; weight fingerprint `f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4`. The full native receipts are in `receipts/`, and the strict [final collection](final_collection.json) retains all 15 task outcomes. The original parent optimizer failed; its saved checkpoints were used solely as provisional new search starts. The recovery made a numerical time-budget correction and changed the starts; it did not change earnings, initial wealth, transfers, floors, preferences, targets, weights, or housing/estate equations.
