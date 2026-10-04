# CES normalized-share collection

Status: **collected_no_adoption**. This collector neither adopts nor publishes a calibration.

## Chain 0

Postcheck: `full_native_postcheck_passed`. Root report: `/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v4/results/smoke_chain_0/run/native_postcheck/selected_root/phase_b_ge/selected_root`. Internal repeat: `/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v4/results/smoke_chain_0/run/native_postcheck/selected_root/phase_b_ge/selected_repeat_final`.

### Full target fit

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| initial_normalization | 2.1 | 2.0999993753504897 | -6.246495103390259e-07 |  | 0.0 |
| cps_childlessness | 0.19827875100684264 | 0.3245948224305346 | 0.12631607142369197 | 35532.3042455214 | 566.9445599092364 |
| cps_exactly_one | 0.21365532522014702 | 0.05523808330378385 | -0.15841724191636317 | 26952.820824310795 | 676.4085988261188 |
| nchs_mean_age | 25.976263860992496 | 22.564931248295224 | -3.4113326126972723 | 139.82806784479274 | 1627.2058200325907 |
| nchs_share30 | 0.2492780130410667 | 0.033946270427197986 | -0.21533174261386873 | 0.0 | 0.0 |
| wealth_earnings | 6.92658379107299 | 7.049017057194719 | 0.12243326612172911 | 7.595098472533724 | 0.1138498019352048 |
| bequest_wealth | 0.007291023472616158 | 0.007312026315797355 | 2.1002843181197606e-05 | 5165289.256198346 | 0.002278509409576336 |
| old_dispersion | 3.51593508651872 | 2.8891320808999343 | -0.6268030056187857 | 0.0 | 0.0 |
| mean_rooms | 5.729434240102641 | 4.5952760659884575 | -1.1341581741141837 | 128.02070205233477 | 164.6749191360439 |
| ownership_30_55 | 0.6762604168538028 | 0.612080580280774 | -0.06417983657302884 | 2339.3623724673616 | 9.63595390814985 |
| first_birth_rooms | 1.465 | 1.2003573569974555 | -0.2646426430025446 | 137.5652749002964 | 9.634484243308407 |
| family_rooms | 0.38509964969278165 | 0.5918820523944097 | 0.20678240270162807 | 280.52808370152104 | 11.995089689737888 |
| recent_parent_ownership | 0.12760836356692162 | -0.05668328117906818 | -0.1842916447459898 | 27055.822957508266 | 918.9080167372239 |
| early_fertility | 0.8095276384290021 | 0.6970569293383831 | -0.11247070909061896 | 100.0 | 1.2649660403346636 |

### Parameters

| Parameter | Estimate | Lower | Upper | Status | Near bound |
|---|---:|---:|---:|---|---|
| H0 | 4.119280677333723 | 0.2 | 80.0 | derived housing supply coefficient at N0=1; reference bounds advisory | False |
| beta_annual | 0.9663191380998087 | 0.94 | 0.99 | supplied primitive; reference calibration bounds advisory | False |
| chi | 1.0500762402240174 | 0.1 | 5.0 | supplied primitive; reference calibration bounds advisory | False |
| first_birth_fixed_cost | 1.9 | 0.0 | 8.0 | supplied primitive; reference calibration bounds advisory | False |
| kappa_fert | 0.11652185618155607 | 0.02 | 50.0 | supplied primitive; reference calibration bounds advisory | True |
| kappa_fert_continuation | 0.40070847699255835 | 0.02 | 50.0 | supplied primitive; reference calibration bounds advisory | True |
| theta0 | 0.10097400014250629 | 0.0 | 8.0 | supplied primitive; reference calibration bounds advisory | False |
| delta_alpha_jump | 0.07780442689806688 | 0.0 | 0.25 | supplied primitive; reference calibration bounds advisory | False |
| child_benefit_curvature | 0.0629608477756522 | 0.0 | 0.8 | supplied primitive; reference calibration bounds advisory | False |
| tenure_choice_kappa | 0.014123854930856623 | 0.001 | 0.1 | supplied primitive; reference calibration bounds advisory | False |
| psi_child | 0.17892072066041628 | 0.01 | 0.5 | supplied primitive; reference calibration bounds advisory | False |
| child_benefit_CRRA_coefficient | 0.1676557204030058 |  |  | derived from supplied benefit and curvature |  |
| theta1 | 0.008193084126995582 |  |  | fixed external restriction |  |
| sigma | 2.0 |  |  | fixed |  |
| alpha_cons | 0.733 |  |  | fixed alpha0=.733 in normalized CES-limit share experiment |  |
| delta_alpha | 0.03897536437154123 | 0.0 | 0.25 | supplied primitive; reference calibration bounds advisory | False |
| h_P | 0.0 | 0.0 | 0.0 | fixed zero; housing floor removed in this experiment | True |
| utility_reference_rent | 0.11046592704873838 |  |  | inactive legacy input; unused by normalized CES-limit utility |  |
| q_annual | 0.020000000000000018 |  |  | author-retained 2% annual real rate |  |
| financed_share | 0.8 |  |  | inherited credit contract |  |
| housing_supply_elasticity | 0.63 |  |  | fixed provisional external mapping |  |
| payroll_tax | 0.08028070961950022 |  |  | derived from adopted pension ratio |  |
| pension_period | 0.917784047463731 |  |  | balanced PAYGO |  |
| annual_depreciation | 0.01416143718381309 |  |  | adopted |  |
| period_depreciation | 0.05545379079326218 |  |  | compounded |  |
| annual_property_tax | 0.010598360773872594 |  |  | adopted |  |
| period_property_tax | 0.042393443095490375 |  |  | linear period convention |  |
| selling_cost | 0.06 |  |  | retained |  |
| rental_cap | 6.0 |  |  | retained provisional |  |
| wealth_grid_nodes | 120.0 |  |  | retained exact grid |  |
| income_states | 9.0 |  |  | retained B15 |  |
