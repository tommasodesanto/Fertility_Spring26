# Matched soft timing calibration readout

Collection status: **complete**; {'no_admissible_candidate': 2, 'verified': 46}. No calibration is adopted by this collector.

## original
Lowest verified native loss: **18.4453051941** (chain 15, v3).
Full CSVs: `original_target_fit.csv` and `original_parameters.csv`.

### Target fit

| moment | role | target | model | gap | weight | loss_contribution |
| --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | normalization | 2.1 | 2.0999999985729323 | -1.4270677972660906e-09 |  |  |
| cps_childlessness | scored | 0.19827875100684264 | 0.200261507466865 | 0.0019827564600223557 | 35532.3042455214 | 0.13968897131071656 |
| cps_exactly_one | scored | 0.21365532522014702 | 0.2148021536591023 | 0.0011468284389552774 | 26952.820824310795 | 0.03544876686505552 |
| nchs_mean_age | scored | 25.976263860992496 | 25.972822848693905 | -0.003441012298591062 | 139.82806784479274 | 0.0016556434154984963 |
| nchs_share30 | validation | 0.2492780130410667 | 0.23234928838554628 | -0.016928724655520422 | 0.0 | 0.0 |
| wealth_earnings | scored | 6.92658379107299 | 6.581195706298087 | -0.34538808477490335 | 7.595098472533724 | 0.9060415436254782 |
| bequest_wealth | scored | 0.007291023472616158 | 0.00718250772463619 | -0.0001085157479799679 | 5165289.256198346 | 0.06082472913043338 |
| old_dispersion | validation | 3.51593508651872 | 2.8792360061960385 | -0.6366990803226815 | 0.0 | 0.0 |
| mean_rooms | scored | 5.729434240102641 | 5.945572774974916 | 0.2161385348722744 | 128.02070205233477 | 5.980597995169983 |
| ownership_30_55 | scored | 0.6762604168538028 | 0.6603880980410459 | -0.015872318812756858 | 2339.3623724673616 | 0.5893567426895043 |
| first_birth_rooms | scored | 1.465 | 1.351172101322728 | -0.11382789867727205 | 137.5652749002964 | 1.7824044493356326 |
| family_rooms | validation | 0.38509964969278165 | 0.3200979931492336 | -0.06500165654354806 | 0.0 | 0.0 |
| recent_parent_ownership | scored | 0.12760836356692162 | 0.12065829000860762 | -0.006950073558314007 | 27055.822957508266 | 1.306891552063457 |
| early_fertility | scored | 0.8095276384290021 | 0.5330788221814269 | -0.2764488162475752 | 100.0 | 7.64239480046856 |

### Parameters and restrictions

| parameter | estimate | lower | upper | near_bound | status |
| --- | --- | --- | --- | --- | --- |
| H0 | 6.79390271408785 | 0.2 | 80.0 | False | derived calibrated housing supply coefficient at N0=1 |
| beta_annual | 0.9670336401689759 | 0.94 | 0.99 | False | free in experimental utility calibration |
| chi | 1.0905005359076558 | 0.1 | 5.0 | False | free in experimental utility calibration |
| first_birth_fixed_cost | 0.28606371352720844 | 0.0 | 8.0 | False | free in experimental utility calibration |
| kappa_fert | 0.10777020144717818 | 0.02 | 50.0 | True | free in experimental utility calibration |
| kappa_fert_continuation | 0.3560939268151674 | 0.02 | 50.0 | True | free in experimental utility calibration |
| theta0 | 0.15831562642314595 | 0.0 | 8.0 | False | free in experimental utility calibration |
| delta_alpha_jump | 0.0 |  |  |  | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.08407417663994812 | 0.0 | 0.8 | False | free in experimental utility calibration |
| tenure_choice_kappa | 0.016927769397298127 | 0.001 | 0.1 | False | free in experimental utility calibration |
| psi_child | 0.17361208520421256 | 0.01 | 0.5 | False | free in experimental utility calibration |
| child_benefit_CRRA_coefficient | 0.15901579208592387 |  |  |  | derived from fixed benefit and proposed curvature |
| theta1 | 0.008193084126995582 |  |  |  | fixed external restriction |
| sigma | 2.0 |  |  |  | fixed |
| alpha_cons | 0.733 |  |  |  | fixed CEX childless expenditure share |
| delta_alpha | 0.0 |  |  |  | fixed zero under experimental utility contract |
| h_P | 2.6 | 0.1 | 2.6 | True | free in experimental utility calibration |
| utility_reference_rent | 0.11046592704873838 |  |  |  | retained inactive normalization; compensation off |
| q_annual | 0.020000000000000018 |  |  |  | author-retained 2% annual real rate |
| financed_share | 0.8 |  |  |  | inherited credit contract |
| housing_supply_elasticity | 0.63 |  |  |  | fixed provisional external mapping |
| payroll_tax | 0.08028070961950022 |  |  |  | derived from adopted pension ratio |
| pension_period | 0.917784047463731 |  |  |  | balanced PAYGO |
| annual_depreciation | 0.01416143718381309 |  |  |  | adopted |
| period_depreciation | 0.05545379079326218 |  |  |  | compounded |
| annual_property_tax | 0.010598360773872594 |  |  |  | adopted |
| period_property_tax | 0.042393443095490375 |  |  |  | linear period convention |
| selling_cost | 0.06 |  |  |  | retained |
| rental_cap | 6.0 |  |  |  | retained provisional |
| wealth_grid_nodes | 120.0 |  |  |  | retained exact grid |
| income_states | 9.0 |  |  |  | retained B15 |

## alternative
Lowest verified native loss: **13.7711314635** (chain 13, v3).
Full CSVs: `alternative_target_fit.csv` and `alternative_parameters.csv`.

### Target fit

| moment | role | target | model | gap | weight | loss_contribution |
| --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | normalization | 2.1 | 2.0999999998568617 | -1.4313839002966233e-10 |  |  |
| cps_childlessness | scored | 0.19827875100684264 | 0.19928030307036224 | 0.0010015520635195951 | 35532.3042455214 | 0.03564268662570389 |
| cps_exactly_one | scored | 0.21365532522014702 | 0.21567067757428574 | 0.0020153523541387164 | 26952.820824310795 | 0.10947279293768183 |
| nchs_mean_age | scored | 25.976263860992496 | 25.96765158349927 | -0.008612277493227793 | 139.82806784479274 | 0.010371232871325497 |
| nchs_share30 | validation | 0.2492780130410667 | 0.23247004842704286 | -0.01680796461402384 | 0.0 | 0.0 |
| wealth_earnings | scored | 6.92658379107299 | 6.740519719992495 | -0.18606407108049527 | 7.595098472533724 | 0.26294108286804535 |
| bequest_wealth | scored | 0.007291023472616158 | 0.00683518502095964 | -0.0004558384516565178 | 5165289.256198346 | 1.0732887087221668 |
| old_dispersion | validation | 3.51593508651872 | 2.9036902986504103 | -0.6122447878683097 | 0.0 | 0.0 |
| mean_rooms | scored | 5.729434240102641 | 5.886397429645046 | 0.15696318954240507 | 128.02070205233477 | 3.1541027331613143 |
| ownership_30_55 | scored | 0.6762604168538028 | 0.6713911717422071 | -0.004869245111595699 | 2339.3623724673616 | 0.05546522435834509 |
| first_birth_rooms | scored | 1.465 | 1.3775251968185076 | -0.08747480318149248 | 137.5652749002964 | 1.0526276370214847 |
| family_rooms | validation | 0.38509964969278165 | 0.29906067130487646 | -0.0860389783879052 | 0.0 | 0.0 |
| recent_parent_ownership | scored | 0.12760836356692162 | 0.12369235439971915 | -0.003916009167202472 | 27055.822957508266 | 0.41490450272300256 |
| early_fertility | scored | 0.8095276384290021 | 0.5338046821453732 | -0.27572295628362886 | 100.0 | 7.602314862178392 |

### Parameters and restrictions

| parameter | estimate | lower | upper | near_bound | status |
| --- | --- | --- | --- | --- | --- |
| H0 | 6.40569359569417 | 0.2 | 80.0 | False | derived calibrated housing supply coefficient at N0=1 |
| beta_annual | 0.9663191380998087 | 0.94 | 0.99 | False | free in experimental utility calibration |
| chi | 1.0500762402240174 | 0.1 | 5.0 | False | free in experimental utility calibration |
| first_birth_fixed_cost | 0.3045994545418478 | 0.0 | 8.0 | False | free in experimental utility calibration |
| kappa_fert | 0.11652185618155607 | 0.02 | 50.0 | True | free in experimental utility calibration |
| kappa_fert_continuation | 0.40070847699255835 | 0.02 | 50.0 | True | free in experimental utility calibration |
| theta0 | 0.10097400014250629 | 0.0 | 8.0 | False | free in experimental utility calibration |
| delta_alpha_jump | 0.0 |  |  |  | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.0629608477756522 | 0.0 | 0.8 | False | free in experimental utility calibration |
| tenure_choice_kappa | 0.014123854930856623 | 0.001 | 0.1 | False | free in experimental utility calibration |
| psi_child | 0.17892072066041628 | 0.01 | 0.5 | False | free in experimental utility calibration |
| child_benefit_CRRA_coefficient | 0.1676557204030058 |  |  |  | derived from fixed benefit and proposed curvature |
| theta1 | 0.008193084126995582 |  |  |  | fixed external restriction |
| sigma | 2.0 |  |  |  | fixed |
| alpha_cons | 0.733 |  |  |  | fixed CEX childless expenditure share |
| delta_alpha | 0.0 |  |  |  | fixed zero under experimental utility contract |
| h_P | 2.593759507364224 | 0.1 | 2.6 | True | free in experimental utility calibration |
| utility_reference_rent | 0.11046592704873838 |  |  |  | retained inactive normalization; compensation off |
| q_annual | 0.020000000000000018 |  |  |  | author-retained 2% annual real rate |
| financed_share | 0.8 |  |  |  | inherited credit contract |
| housing_supply_elasticity | 0.63 |  |  |  | fixed provisional external mapping |
| payroll_tax | 0.08028070961950022 |  |  |  | derived from adopted pension ratio |
| pension_period | 0.917784047463731 |  |  |  | balanced PAYGO |
| annual_depreciation | 0.01416143718381309 |  |  |  | adopted |
| period_depreciation | 0.05545379079326218 |  |  |  | compounded |
| annual_property_tax | 0.010598360773872594 |  |  |  | adopted |
| period_property_tax | 0.042393443095490375 |  |  |  | linear period convention |
| selling_cost | 0.06 |  |  |  | retained |
| rental_cap | 6.0 |  |  |  | retained provisional |
| wealth_grid_nodes | 120.0 |  |  |  | retained exact grid |
| income_states | 9.0 |  |  |  | retained B15 |

Search optimization convergence is not certified by a selected-point native check.
