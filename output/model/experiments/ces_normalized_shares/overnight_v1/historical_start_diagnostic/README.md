# Historical fitted-vector starting-point test

Source: September 27 primary final de_0093, all eleven transferable coordinates. Current experimental utility, targets, timing, grid and fixed inputs retained. No reference-rent correction.

The original numerical price caps did not bracket birth renewal. A fixed-price diagnostic at 10.110 found the opposite sign. A separate diagnostic doubled the numerical bracketing origin only, preserving the identical economic-input fingerprint. The stationary root and exact repeat passed; the adapter rejected derived housing supply coefficient H0=0.198580 below its retained 0.2 calibration lower bound. No bound was relaxed. The four existing searches remain unchanged.

Diagnostic loss: 4417.699161134. This is not an accepted calibration candidate. Price: 7.670590398. All 17 plots and both complete tables repeat exactly.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| initial_normalization | 2.1 | 2.0999999967570244 | -3.2429756657847975e-09 |  | 0.0 |
| cps_childlessness | 0.19827875100684264 | 0.20186855800915585 | 0.0035898070023132056 | 35532.3042455214 | 0.45789465372507976 |
| cps_exactly_one | 0.21365532522014702 | 0.21159897400838648 | -0.0020563512117605376 | 26952.820824310795 | 0.113972167331766 |
| nchs_mean_age | 25.976263860992496 | 26.13830013317497 | 0.16203627218247263 | 139.82806784479274 | 3.6712912821046038 |
| nchs_share30 | 0.2492780130410667 | 0.23656577348123806 | -0.012712239559828642 | 0.0 | 0.0 |
| wealth_earnings | 6.92658379107299 | 5.737049077695018 | -1.189534713377972 | 7.595098472533724 | 10.747009914675166 |
| bequest_wealth | 0.007291023472616158 | 0.007033089001481633 | -0.0002579344711345251 | 5165289.256198346 | 0.34364768284838393 |
| old_dispersion | 3.51593508651872 | 3.3994452149791954 | -0.11648987153952461 | 0.0 | 0.0 |
| mean_rooms | 5.729434240102641 | 0.7727316144235331 | -4.9567026256791085 | 128.02070205233477 | 3145.3279443576553 |
| ownership_30_55 | 0.6762604168538028 | 0.09285968341728557 | -0.5834007334365172 | 2339.3623724673616 | 796.2169922901752 |
| first_birth_rooms | 1.465 | 0.3195144341083669 | -1.1454855658916332 | 137.5652749002964 | 180.5045121027938 |
| family_rooms | 0.38509964969278165 | 0.09642698225932556 | -0.2886726674334561 | 280.52808370152104 | 23.37694072140004 |
| recent_parent_ownership | 0.12760836356692162 | 0.031699895744653664 | -0.09590846782226796 | 27055.822957508266 | 248.87120720189523 |
| early_fertility | 0.8095276384290021 | 0.5254898106761208 | -0.28403782775288133 | 100.0 | 8.06774875945755 |

| Parameter | Value | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| H0 | 0.19858026882344829 | 0.2 | 80.0 | True | derived housing supply coefficient at N0=1; reference bounds advisory |
| beta_annual | 0.9602456707427426 | 0.94 | 0.99 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.0889178301004045 | 0.1 | 5.0 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 0.6595894538306366 | 0.0 | 8.0 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.19482764369932495 | 0.02 | 50.0 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.3474664229788561 | 0.02 | 50.0 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.11358885656635645 | 0.0 | 8.0 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.14198430659647016 | 0.0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_curvature | 0.07130705056729134 | 0.0 | 0.8 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.005 | 0.001 | 0.1 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.14416490417555738 | 0.01 | 0.5 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.1338849300634822 |  |  |  | derived from supplied benefit and curvature |
| theta1 | 0.008193084126995582 |  |  |  | fixed external restriction |
| sigma | 2.0 |  |  |  | fixed |
| alpha_cons | 0.733 |  |  |  | fixed alpha0=.733 in normalized CES-limit share experiment |
| delta_alpha | 0.0 | 0.0 | 0.25 | True | supplied primitive; reference calibration bounds advisory |
| h_P | 0.0 | 0.0 | 0.0 | True | fixed zero; housing floor removed in this experiment |
| utility_reference_rent | 0.11046592704873838 |  |  |  | inactive legacy input; unused by normalized CES-limit utility |
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

Reproduction packet: `point.json`, `check.py`/`launch.sh` (original caps, job 19136046), `upper_price.py`/`upper_launch.sh` (one higher fixed price, job 19136629), and `wide_check.py`/`wide_launch.sh` (wider numerical caps only, job 19136864). Each Slurm task used one CPU, 24 GiB, a 30-minute wall cap, at most 1,500 seconds for GE and 32 lifecycle solves. Total observed work: eleven original-cap trials, one higher-price solve, and eight wider-cap lifecycle solves. The fixed-price packet helper is the unchanged `code/cluster/ces_normalized_shares_calibration/diagnose_starts.py`. Original parameter table: [September 27 primary final](../../../../supervised_calibration_20260927/primary_final/parameters.csv).
