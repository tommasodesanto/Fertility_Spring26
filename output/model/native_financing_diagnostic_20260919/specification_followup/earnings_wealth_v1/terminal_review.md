# Earnings and wealth: terminal overnight review

**September 22, 08:15 EDT. The earnings candidate is implemented and has passed the full numerical smoke, but the requested longer calibration did not finish.** Search stopped at 07:50 EDT on a proposal timeout. Nothing was adopted, no target or gate was relaxed, and no restart was launched.

## What is established

The candidate uses directly estimated four-year PSID persistent AR(1) plus iid log earnings, without a fixed permanent type. The reviewed annual and multi-year literature supports this process family and direct estimation at the model frequency; it does not supply a drop-in four-year U.S. parameter set. [Literature review](literature/multiyear_review.md), [timing review](literature/timing_review.md), [complete data validation](data_validation/README.md).

The four-year estimates are persistence 0.77614467, persistent innovation SD 0.43743545, and iid SD 0.16499061. Persistent stationary variance is 0.48126268. The 199-draw person bootstrap gives persistence percentiles [0.72114,0.83477], persistent-variance percentiles [0.43019,0.52492], and iid-variance percentiles [approximately 0,0.06398]. These are bootstrap uncertainty summaries, not proof of specification validity. Three covariance moments estimate three parameters; exact point fit is mechanical. The corrected same-observer simulation audit is saved with the data evidence.

The implementation has 45 income states and a 160-node wealth grid: all 120 original knots retained, plus 40 upper-tail knots through 3000. Zero assets at age18 and a stationary initial income distribution remain explicit external candidate assumptions. Four-year income is known at the decision. Current income enters purchase eligibility under the reviewed accounting; it is not added to transaction wealth a second time. The frozen September14 reference is preserved.

Two exact starting-point repetitions and a nearby annual-discount-factor test passed the unchanged smoke gates. The smoke completed 18 stationary solves. Its starting loss is 1160.3761384370134; the single-repeat beta 0.985 test gives 1040.4468628084285. [Full smoke target/parameter tables and both original 17-figure sets](smoke_v5/README.md).

## Why the longer run stopped

Eight proposals were attempted; seven were scored. Proposal 006 completed six stationary solves taking 424–478 seconds each and exceeded its 3,100-second native limit during the seventh. Its final native heartbeat records seven started solves, while the saved solve table records six completed solves. It has no complete objective.

The inner timeout arrived at the controller as traceback text. The frozen classifier recognizes typed outer timeouts but treats this traceback as a fatal contract error. The recorded label is therefore `fatal_contract_error`; the observed cause is a runtime timeout, with no evidence of a target/source fingerprint mismatch. The controller stopped before final selection verification. The failure and frozen code are preserved. [Lead failure/count review](staging/v5_terminal_lead_review.json), [collected terminal evidence](search_readout/README.md).

Search totals: 49 stationary solves started,48 completed, 1 incomplete. V5 including smoke: 67 started, 66 completed, 1 incomplete, 10 scored objectives. All attempts including earlier failures and the deliberately stopped runtime pilot: 102 started, 98 completed, 4 incomplete, 12 scored objectives. These counts distinguish a completed stationary solve from a complete normalized/scored candidate.

**The planned original-selection versus both-repetition checks of price, value functions, household distributions, all native moments, fertility normalization, and full 13/17 tables were unrun.** There is no final verified selected calibration. No reporting override can supply that missing evidence.

## What the completed proposals show

The best valid search proposal is005, loss 1087.2816148857435, evaluated once. It improves the repeated anchor by 6.3%, but is worse than the single-repeat smoke test 1040.4468628084285. The frozen search excluded that smoke test from its candidate pool; this limitation is preserved rather than rewritten after the run.

The improvement trades fit across moments. First-birth age moves from 28.21 at the anchor to 26.42 against 25.98, and the 30+ first-birth share nearly matches. But childlessness falls to 11.25% against 19.83%; mean rooms rise to 6.94 against 5.56; ownership ages 30–55 is 50.87% against 64.83%. The first-birth housing response still exceeds its target. These seven valid proposals neither establish an optimum nor show that the remaining targets are unreachable.

### Complete target table: best valid search proposal 005

The unchanged objective contains 12 scored rows and one separately imposed fertility normalization. Weights are the saved objective weights; this diagnostic is not certified SMM. The recent-parent observer compares current births from previously empty-dependent homes with current empty homes, including former parents. It is not the obsolete lifetime-childless comparison.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Fertility normalization (unscored) | 2.1 | 2.10014054 | 0.000140539608 | — | — |
| Childlessness | 0.198278751 | 0.112474749 | -0.0858040021 | 35532.3042 | 261.600435 |
| Exactly one child | 0.213655325 | 0.295433276 | 0.0817779513 | 26952.8208 | 180.250582 |
| Mean first-birth age | 25.9762639 | 26.4236886 | 0.44742472 | 139.828068 | 27.9920243 |
| Share of first births at 30+ | 0.249278013 | 0.250909283 | 0.00163126948 | 13866.0654 | 0.0368981563 |
| Wealth / annual labor earnings | 6.14586139 | 6.02159224 | -0.124269155 | 7.59509847 | 0.11728976 |
| Annual bequests / wealth | 0.0088 | 0.00436250343 | -0.00443749657 | 5165289.26 | 101.711652 |
| Old-age wealth dispersion | 3.51593509 | 4.43515002 | 0.919214934 | 10.6163615 | 8.97035937 |
| Mean occupied rooms | 5.56109738 | 6.93629067 | 1.3751933 | 128.020702 | 242.107197 |
| Ownership ages 30–55 | 0.648334034 | 0.508724682 | -0.139609352 | 2339.36237 | 45.595977 |
| First-birth rooms response | 0.720246262 | 1.02193473 | 0.301688469 | 137.565275 | 12.5206317 |
| Family rooms response | 0.347066932 | 0.286763258 | -0.0603036737 | 280.528084 | 1.02014965 |
| Recent-parent ownership difference | 0.162895509 | 0.0757739013 | -0.0871216078 | 27055.823 | 205.358419 |

### Complete parameter review: best valid search proposal 005

The actual search bound for annual beta is 0.99, not the scorer's generic0.9995. Near-bound means within1% of the raw bound width; this differs from distance in transformed search coordinates. The form and units are checked against the [reviewed literature rubric](literature/parameter_validation_rubric.md). Utility coefficient levels are not portable literature estimates.

| Parameter | Value | Actual bounds or restriction | Near bound | Units and review |
|---|---:|---|---|---|
| `beta_annual` | 0.99 | [0.94, 0.99] | Yes | Annual; $\beta_4=0.96059601$. At upper bound. |
| `kappa_fert` | 0.350501962 | [0.02, 50] | Yes | First-birth taste-shock scale; model-specific utility units. |
| `kappa_fert_continuation` | 0.676803076 | [0.02, 50] | No | Subsequent-birth taste-shock scale; no required ordering versus first-birth scale. |
| `chi` | 0.9233887 | [0.1, 5] | No | Owner housing-service multiplier; model-specific. Value below one is permitted. |
| `H0` | 10.5161337 | [0.2, 80] | No | Housing supply intercept in model quantity units; interpret with rents and rooms. |
| `theta0` | 0 | [0, 8] | Yes | Bequest-utility scale is zero at lower bound. Intentional bequest utility switches off. |
| `theta1` | 0.11962943 | [0.02, 16] | Yes | Wealth shift inside bequest utility; inactive at theta0=0, hence unidentified here. |
| `first_birth_fixed_cost` | 0.123550335 | [0, 8] | No | One-time first-birth utility cost; not dollars. |
| `h_P` | 1.91460312 | [0.1, 2.3] | No | 1.9146 model room units whenever resident children are present, independent of their number. |
| `hbar_child_rooms` | 0 | zero restriction | — | Zero additional per-child housing-floor slope; imposed restriction. |
| `psi_child` | 0.0784853802 | normalized to 2.1 | — | Derived to meet the unchanged 2.1 fertility normalization; not a free estimate. |
| `payroll_tax` | 0.179 | externally fixed | — | Fixed rate 0.179 on the maintained payroll base. |
| `pension_period` | 2.04636139 | budget derived | — | Four-year benefit from the maintained fiscal balance rule. |
| `housing_supply_elasticity` | 0.63 | externally fixed | — | Fixed elasticity for the maintained housing-supply object. |
| `tenure_choice_kappa` | 0.005 | externally fixed | — | Fixed tenure taste-shock scale in model utility units. |
| `alpha_cons` | 0.733 | externally fixed | — | Fixed consumption exponent in the consumption/housing composite. |
| `sigma` | 2 | externally fixed | — | Fixed utility curvature; not an estimated coefficient. |

At zero bequest scale, the wealth-shift parameter drops out of the implemented bequest utility. Its reported value is not identified at this point. Nonzero accidental bequest flows can still arise. Annual beta is at its upper bound. The housing floor is interior, but the rooms and ownership fit does not validate it. Nine free coordinates and12 scored moments do not by themselves establish local identification.

## Visual review and remaining numerical limits

All 17 original case005 figures were visually inspected and their original files retained. [Individual figure links, full-precision tables, and fingerprint manifest](search_readout/README.md).

Fertility flows have a life-cycle hump; ownership rises during working life and drops sharply in the last age cell. Liquid wealth accumulates then decumulates, with a steep terminal decline. Housing use is heavily concentrated at the largest owner rung; the saved renter 25–45 cap share is 55.97%. Housing profiles show a retirement discontinuity. These are descriptive observations requiring economic and numerical review, not proof of an impossible model or a policy mechanism. The housing-market relative residual 7.52e-6 passes the existing numerical gate.

The policy plots retain a high-wealth housing decline from10 to 6 near the expanded upper boundary. The160-node domain correction is not convergence evidence. Low-wealth regions and boundary mass cannot be read precisely on a wealth axis extending to 3000, and 45-state legends overlap or clip. The stable diagnostic set was preserved rather than replaced. Its income-state profiles are cross-sectional conditional profiles, not paths of fixed permanent types.

The allocation-output correction fixed budget reporting at positive current surplus. It did not remove finite infeasibility values from continuation interpolation, which may distort choices in tiny-mass states. That limitation, wealth-grid density/domain sensitivity, and the incomplete best-point repeats still prevent numerical certification.

## Recommended next step, not launched

Keep the persistent-plus-iid direct-period earnings construction as the candidate to assess; this timeout supplies no economic reason to abandon it. Do not adopt the joint parameter vector or interpret policies from it.

Before another long search, prepare and independently review a separate operational correction that preserves typed inner-timeout information, retains a usable terminal best-candidate report on failure, and explicitly accounts for normalization cases requiring more than six stationary solves. Test timeout/failure paths with stub evaluations before spending model time. Any changed duration or selection pool must be declared in a new frozen plan; the current evidence stays unchanged.

Then prioritize the existing unresolved continuation/interpolation and wealth-grid diagnostics and exact repetitions of a declared candidate. The eventual calibration should address the observed timing/childlessness/rooms trade-off under unchanged targets, rather than dropping moments or inferring structural impossibility from seven proposals. Initial wealth at 18, initial earnings risk, coarse-period information, and retirement wealth-denominator alignment remain choices or validation gaps to close before paper adoption.

The source, empirical work, smoke checks, and partial search are recoverable. The calibrated, fully defended earnings-and-wealth block requested for the paper remains unfinished.
