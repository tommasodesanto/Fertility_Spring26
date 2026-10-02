# Fixed-price credit diagnostic readout

Strict wealth-only purchase origination, with the financed share changed from 80% to 90% for both buyers and owner-stayers. The price, housing-supply coefficient, initial population, and ten calibrated coordinates are fixed. This is a partial-equilibrium diagnostic; birth replacement and housing residuals are measured, not cleared. The 80% replay matches the completed strict-80 target fit; the 90% repeat matches the native tables, closure and 17 plots exactly.

Scored loss: 217.206885943 at 80%; 100.508199287 at 90%; change -116.698686656.

Birth replacement residual: -3.61973095941e-09 to 7.92690752152e-05. Housing-market residual: -8.881784197e-16 to 0.0748797491518.

## Full target comparison

| moment | role | target | strict80_model | strict90_model | strict80_gap | strict90_gap | weight | strict80_loss_contribution | strict90_loss_contribution | model_relative_change_percent |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | normalization | 2.1 | 2.099999992398618 | 2.1001664650580056 | -7.601382190358663e-09 | 0.00016646505800554934 |  |  |  | 0.007927269523348655 |
| cps_childlessness | scored | 0.19827875100684264 | 0.2031372818716155 | 0.20314024526905705 | 0.004858530864772864 | 0.004861494262214411 | 35532.3042455214 | 0.8387514889430407 | 0.8397749720086967 | 0.001458815149165727 |
| cps_exactly_one | scored | 0.21365532522014702 | 0.21059017196374813 | 0.21054026774639237 | -0.0030651532563988892 | -0.003115057473754651 | 26952.820824310795 | 0.2532261849848665 | 0.26153893569922837 | -0.023697315449437262 |
| nchs_mean_age | scored | 25.976263860992496 | 26.00913106187919 | 26.002359803972496 | 0.03286720088669526 | 0.02609594297999962 | 139.82806784479274 | 0.1510496749694374 | 0.09522266810705661 | -0.02603415658364715 |
| nchs_share30 | validation | 0.2492780130410667 | 0.23180034086546736 | 0.23209191131374518 | -0.01747767217559934 | -0.017186101727321518 | 0.0 | 0.0 | 0.0 | 0.12578516804125178 |
| wealth_earnings | scored | 6.92658379107299 | 6.55855642338012 | 6.520670766027885 | -0.3680273676928705 | -0.4059130250451055 | 7.595098472533724 | 1.0287116064302901 | 1.2514093155949568 | -0.5776523812035712 |
| bequest_wealth | scored | 0.007291023472616158 | 0.006832509632245122 | 0.006661290707220777 | -0.0004585138403710356 | -0.0006297327653953808 | 5165289.256198346 | 1.0859242862179517 | 2.048364441180339 | -2.50594487589597 |
| old_dispersion | validation | 3.51593508651872 | 2.902032662787911 | 2.9835538476086745 | -0.613902423730809 | -0.5323812389100455 | 0.0 | 0.0 | 0.0 | 2.809106384848479 |
| mean_rooms | scored | 5.729434240102641 | 5.916861125965659 | 5.991740875117472 | 0.18742688586301792 | 0.26230663501483065 | 128.02070205233477 | 4.497218444704822 | 8.808435058884278 | 1.2655316316823644 |
| ownership_30_55 | scored | 0.6762604168538028 | 0.6016840916994495 | 0.6446341320339892 | -0.0745763251543533 | -0.031626284819813555 | 2339.3623724673616 | 13.01066391274162 | 2.3398814571025177 | 7.138304124549459 |
| first_birth_rooms | scored | 1.465 | 1.0322110554439563 | 1.1731404989123089 | -0.4327889445560438 | -0.2918595010876912 | 137.5652749002964 | 25.766838595999705 | 11.718080896076694 | 13.65316160150392 |
| family_rooms | validation | 0.38509964969278165 | 0.31824396860574655 | 0.3250236085681504 | -0.0668556810870351 | -0.060076041124631274 | 0.0 | 0.0 | 0.0 | 2.1303278714459215 |
| recent_parent_ownership | scored | 0.12760836356692162 | 0.050093476866603814 | 0.07851620668396919 | -0.07751488670031781 | -0.04909215688295243 | 27055.822957508266 | 162.56647228335314 | 65.2056119734634 | 56.739383239565406 |
| early_fertility | scored | 0.8095276384290021 | 0.5265430193314089 | 0.5277497193425142 | -0.2829846190975932 | -0.28177791908648786 | 100.0 | 8.008029464580991 | 7.939879568471129 | 0.22917405925114745 |

## All reported parameters

| parameter | strict80_estimate | strict90_estimate | lower | upper | strict80_near_bound | strict90_near_bound | strict80_status | strict90_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H0 | 6.778473404808042 | 6.778473404808042 | 0.2 | 80.0 | False | False | fixed strict-80 physical supply coefficient for diagnostic | fixed strict-80 physical supply coefficient for diagnostic |
| beta_annual | 0.9670843817494936 | 0.9670843817494936 | 0.94 | 0.99 | False | False | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| chi | 1.093429912990833 | 1.093429912990833 | 0.1 | 5.0 | False | False | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| first_birth_fixed_cost | 0.38540954564501306 | 0.38540954564501306 | 0.0 | 8.0 | False | False | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| kappa_fert | 0.126367646227124 | 0.126367646227124 | 0.02 | 50.0 | True | True | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| kappa_fert_continuation | 0.3511943513636155 | 0.3511943513636155 | 0.02 | 50.0 | True | True | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| theta0 | 0.10416394115934847 | 0.10416394115934847 | 0.0 | 8.0 | False | False | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| delta_alpha_jump | 0.0 | 0.0 | 0.0 | 0.25 | False | False | free in evening DUE calibration | free in evening DUE calibration |
| child_benefit_curvature | 0.09889818563569602 | 0.09889818563569602 | 0.0 | 0.8 | False | False | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| tenure_choice_kappa | 0.013862580239215504 | 0.013862580239215504 | 0.001 | 0.1 | False | False | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| psi_child | 0.17141586769606804 | 0.17141586769606804 | 0.01 | 0.5 | False | False | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| child_benefit_CRRA_coefficient | 0.15446314939175837 | 0.15446314939175837 |  |  |  |  | derived from fixed benefit and proposed curvature | derived from fixed benefit and proposed curvature |
| theta1 | 0.008193084126995582 | 0.008193084126995582 |  |  |  |  | fixed external restriction | fixed external restriction |
| sigma | 2.0 | 2.0 |  |  |  |  | fixed | fixed |
| alpha_cons | 0.733 | 0.733 |  |  |  |  | fixed CEX childless expenditure share | fixed CEX childless expenditure share |
| delta_alpha | 0.0 | 0.0 |  |  |  |  | fixed zero later-child loading | fixed zero later-child loading |
| h_P | 2.2992335366442824 | 2.2992335366442824 | 0.1 | 2.3 | True | True | fixed at strict-80 reference estimate for this experiment | fixed at strict-80 reference estimate for this experiment |
| utility_reference_rent | 0.11046592704873838 | 0.11046592704873838 |  |  |  |  | fixed substantive utility normalization | fixed substantive utility normalization |
| q_annual | 0.020000000000000018 | 0.020000000000000018 |  |  |  |  | author-retained 2% annual real rate | author-retained 2% annual real rate |
| financed_share | 0.8 | 0.9 |  |  |  |  | externally fixed policy input for this experiment | externally fixed policy input for this experiment |
| housing_supply_elasticity | 0.63 | 0.63 |  |  |  |  | fixed provisional external mapping | fixed provisional external mapping |
| payroll_tax | 0.08028070961950022 | 0.08028070961950022 |  |  |  |  | derived from adopted pension ratio | derived from adopted pension ratio |
| pension_period | 0.917784047463731 | 0.917784047463731 |  |  |  |  | balanced PAYGO | balanced PAYGO |
| annual_depreciation | 0.01416143718381309 | 0.01416143718381309 |  |  |  |  | adopted | adopted |
| period_depreciation | 0.05545379079326218 | 0.05545379079326218 |  |  |  |  | compounded | compounded |
| annual_property_tax | 0.010598360773872594 | 0.010598360773872594 |  |  |  |  | adopted | adopted |
| period_property_tax | 0.042393443095490375 | 0.042393443095490375 |  |  |  |  | linear period convention | linear period convention |
| selling_cost | 0.06 | 0.06 |  |  |  |  | retained | retained |
| rental_cap | 6.0 | 6.0 |  |  |  |  | retained provisional | retained provisional |
| wealth_grid_nodes | 120.0 | 120.0 |  |  |  |  | retained exact grid | retained exact grid |
| income_states | 9.0 | 9.0 |  |  |  |  | retained B15 | retained B15 |
