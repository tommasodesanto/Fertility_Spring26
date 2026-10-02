# Strict purchase origination readout

Experimental fixed-coordinate comparison with normalized v2 chain 2, case `0064_nm`. Only purchase origination eligibility excludes current income. The ten calibrated coordinates and target/weight fingerprints are unchanged. Initial population remains \(N_0=1\); the fertility-price root normalizes completed fertility to 2.1. The equilibrium price may change and \(H_0\) is derived separately from housing demand, so neither is held fixed.

Scored loss: baseline 29.9696787536; strict 217.206885943; strict minus baseline +187.237207189.

## Full target fit

| moment | role | target | baseline_model | strict_model | baseline_gap | strict_gap | weight | baseline_loss_contribution | strict_loss_contribution |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | normalization | 2.1 | 2.099999969287473 | 2.099999992398618 | -3.071252718811479e-08 | -7.601382190358663e-09 |  |  |  |
| cps_childlessness | scored | 0.19827875100684264 | 0.20345029271216752 | 0.2031372818716155 | 0.0051715417053248836 | 0.004858530864772864 | 35532.3042455214 | 0.9503059201463746 | 0.8387514889430407 |
| cps_exactly_one | scored | 0.21365532522014702 | 0.21027063732321577 | 0.21059017196374813 | -0.003384687896931249 | -0.0030651532563988892 | 26952.820824310795 | 0.3087745383817932 | 0.2532261849848665 |
| nchs_mean_age | scored | 25.976263860992496 | 25.98292181303565 | 26.00913106187919 | 0.006657952043152449 | 0.03286720088669526 | 139.82806784479274 | 0.006198344092724217 | 0.1510496749694374 |
| nchs_share30 | validation | 0.2492780130410667 | 0.23098517457792037 | 0.23180034086546736 | -0.018292838463146333 | -0.01747767217559934 | 0.0 | 0.0 | 0.0 |
| wealth_earnings | scored | 6.92658379107299 | 6.523137846278569 | 6.55855642338012 | -0.403445944794421 | -0.3680273676928705 | 7.595098472533724 | 1.2362437759076668 | 1.0287116064302901 |
| bequest_wealth | scored | 0.007291023472616158 | 0.006852110588063577 | 0.006832509632245122 | -0.000438912884552581 | -0.0004585138403710356 | 5165289.256198346 | 0.9950646705902234 | 1.0859242862179517 |
| old_dispersion | validation | 3.51593508651872 | 2.9401999570243076 | 2.902032662787911 | -0.5757351294944124 | -0.613902423730809 | 0.0 | 0.0 | 0.0 |
| mean_rooms | scored | 5.729434240102641 | 5.985245018284809 | 5.916861125965659 | 0.2558107781821679 | 0.18742688586301792 | 128.02070205233477 | 8.377566466768986 | 4.497218444704822 |
| ownership_30_55 | scored | 0.6762604168538028 | 0.6593314627576178 | 0.6016840916994495 | -0.016928954096185 | -0.0745763251543533 | 2339.3623724673616 | 0.6704366617429866 | 13.01066391274162 |
| first_birth_rooms | scored | 1.465 | 1.225240668315453 | 1.0322110554439563 | -0.239759331684547 | -0.4327889445560438 | 137.5652749002964 | 7.907876152780069 | 25.766838595999705 |
| family_rooms | validation | 0.38509964969278165 | 0.31027003500843353 | 0.31824396860574655 | -0.07482961468434812 | -0.0668556810870351 | 0.0 | 0.0 | 0.0 |
| recent_parent_ownership | scored | 0.12760836356692162 | 0.11975154873235194 | 0.050093476866603814 | -0.007856814834569681 | -0.07751488670031781 | 27055.822957508266 | 1.6701434877591281 | 162.56647228335314 |
| early_fertility | scored | 0.8095276384290021 | 0.5294014394503301 | 0.5265430193314089 | -0.28012619897867197 | -0.2829846190975932 | 100.0 | 7.847068735423853 | 8.008029464580991 |

## All reported parameters

| parameter | baseline_estimate | strict_estimate | lower | upper | baseline_near_bound | strict_near_bound | baseline_status | strict_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H0 | 6.866751098158112 | 6.778473404808042 | 0.2 | 80.0 | False | False | derived calibrated housing supply coefficient at N0=1 | derived calibrated housing supply coefficient at N0=1 |
| beta_annual | 0.9670843817494936 | 0.9670843817494936 | 0.94 | 0.99 | False | False | free in experimental utility calibration | fixed at reference estimate for this experiment |
| chi | 1.093429912990833 | 1.093429912990833 | 0.1 | 5.0 | False | False | free in experimental utility calibration | fixed at reference estimate for this experiment |
| first_birth_fixed_cost | 0.38540954564501306 | 0.38540954564501306 | 0.0 | 8.0 | False | False | free in experimental utility calibration | fixed at reference estimate for this experiment |
| kappa_fert | 0.126367646227124 | 0.126367646227124 | 0.02 | 50.0 | True | True | free in experimental utility calibration | fixed at reference estimate for this experiment |
| kappa_fert_continuation | 0.3511943513636155 | 0.3511943513636155 | 0.02 | 50.0 | True | True | free in experimental utility calibration | fixed at reference estimate for this experiment |
| theta0 | 0.10416394115934847 | 0.10416394115934847 | 0.0 | 8.0 | False | False | free in experimental utility calibration | fixed at reference estimate for this experiment |
| delta_alpha_jump | 0.0 | 0.0 |  |  |  |  | fixed zero under experimental utility contract | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.09889818563569602 | 0.09889818563569602 | 0.0 | 0.8 | False | False | free in experimental utility calibration | fixed at reference estimate for this experiment |
| tenure_choice_kappa | 0.013862580239215504 | 0.013862580239215504 | 0.001 | 0.1 | False | False | free in experimental utility calibration | fixed at reference estimate for this experiment |
| psi_child | 0.17141586769606804 | 0.17141586769606804 | 0.01 | 0.5 | False | False | free in experimental utility calibration | fixed at reference estimate for this experiment |
| child_benefit_CRRA_coefficient | 0.15446314939175837 | 0.15446314939175837 |  |  |  |  | derived from fixed benefit and proposed curvature | derived from fixed benefit and proposed curvature |
| theta1 | 0.008193084126995582 | 0.008193084126995582 |  |  |  |  | fixed external restriction | fixed external restriction |
| sigma | 2.0 | 2.0 |  |  |  |  | fixed | fixed |
| alpha_cons | 0.733 | 0.733 |  |  |  |  | fixed CEX childless expenditure share | fixed CEX childless expenditure share |
| delta_alpha | 0.0 | 0.0 |  |  |  |  | fixed zero under experimental utility contract | fixed zero under experimental utility contract |
| h_P | 2.2992335366442824 | 2.2992335366442824 | 0.1 | 2.3 | True | True | free in experimental utility calibration | fixed at reference estimate for this experiment |
| utility_reference_rent | 0.11046592704873838 | 0.11046592704873838 |  |  |  |  | retained inactive normalization; compensation off | retained inactive normalization; compensation off |
| q_annual | 0.020000000000000018 | 0.020000000000000018 |  |  |  |  | author-retained 2% annual real rate | author-retained 2% annual real rate |
| financed_share | 0.8 | 0.8 |  |  |  |  | inherited credit contract | inherited credit contract |
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
