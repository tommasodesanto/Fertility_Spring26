# Optimized two-birth diagnostic

The author approved one fully solved two-birth experiment on September 28.
Reference: **2007 stationary reference — block0506, September 28 verified
export**. The preceding fixed-policy replay is in `../two_births_v1/`; its
0.698 age-25 count and 2.362 lifetime count are not results for this model.

## Experimental economic contract

After one successful birth during a four-year period, a household with room
below the existing maximum of three children may choose one additional
attempt. There is no additional opportunity after waiting or failure, and
at most two new births may occur in the period. Housing, consumption and
saving respond to the resulting number of children.

The new conditional wait/try decision uses the existing later-birth choice
scale. Its expected maximized value enters the preceding birth decision.
Thus the experiment adds a choice/shock opportunity as well as relaxing
spacing; it is not the previous fixed-policy mechanical experiment. The
second conception draw is conditionally independent and uses the existing
age-specific probability. These within-period assumptions are experimental,
not adopted. They introduce no new estimated parameter.

The ten reference calibration coordinates stay fixed. Child benefit is
renormalized to completed fertility 2.1 and demographic renewal is enforced.
Housing prices and the inherited pension balance are solved normally. Income,
entry wealth/income, housing/fiscal primitives, targets, weights, bounds and
scientific tolerances are unchanged. The empirical/model within-cell age
projection remains linear in pre/post stocks and does not date the two births
separately. First-birth housing and recent-parent observers must count the
correct families, including families whose first success leads to two children.

## Isolation and verification

No shared source file is edited. Patch builders generate full source copies
inside the Torch job's output, with original/effective hashes and diffs. Only
that worker process installs changed function bodies and the extended owned
policy cache. Original ancestry authentication remains distinct from the
additional effective-source manifest. Checkpoints stay on Torch outside Git.

Before the experimental evaluation, synthetic Bellman/flow/observer/cache
tests and a full reference Bellman/distribution replay with the option disabled
must pass. The replay uses the saved reference price and child benefit and
compares every core value/policy/population array. The normalized experimental
evaluation then retains all market, probability, household budget, transaction,
estate, value, pension and renewal gates. It reports all 14 target rows, all
31 parameter rows and the 17 standard plots. No automatic promotion follows.

## Bounded execution plan

One Torch smoke job: synthetic tests plus one fixed-price reference solve,
20-minute internal cap. One separately authorized diagnostic evaluation after
that smoke: at most 23 stationary solves, 70-minute internal total cap,
one model worker. The reference normalization used six stationary solves and
971 seconds; the changed model's time and normalization difficulty are unknown.
The hard cap replaces any assumption of a speed gain. There is no calibration
search, continuation of an expired job, transition or automatic retry.

The supervisor writes a heartbeat every 15 seconds and native normalization
writes each stationary solve. Unknown failures, scientific gate failures or
deadlines stop the experiment and preserve the failure. A passed point still
requires review of all targets and plots before any decision about further
recalibration. Source changes after the smoke require a new version and smoke.

All implementation changes received lead review and an independent integration
review. The final case also audits extra-probability menu sums, occupied dead
menus, and equality plus independent ownership of solution/policy/parameter
caches. The standard fertility probability plots retain their original meaning:
they show the outer attempt, not the added conditional opportunity.

Smoke **18753550** stopped before any model solve: 14 of 15 synthetic tests
passed; a later canonical import bypassed the already-patched private recent-
parent observer module. `failed_smoke_v1/` preserves the driver and exact log.
The installer now loads canonical modules first and patches all loaded aliases.
Review also identified the live fertility observer's frozen source path; its
bytes match the authenticated current definition (SHA `ca120b5bddc7cf6a8ad48237c821320432c90aeb737b6905ad156d6883474208`).
The installer now patches that exact live function and records every source
path/module name. No target, definition or tolerance was changed for this repair.

Queued smoke **18753837** was cancelled before starting to include both alias
fixes in one test. Revised smoke **18753907** runs the same 15 synthetic tests
and one full flag-off reference replay. At that stage no optimized evaluation
had launched.

**Revised smoke passed:** job 18753907 completed in 105.1s, all 15 tests and
one full reference-price solve. All 12 value/policy/population arrays match
exactly (maximum difference zero). Receipt and effective-source diffs are in
`smoke_v3/`; the receipt SHA is
`c998e63b3655e5c83320a77c689ee871bdb31e6473e3653fa16ef9cf55f60498`.

After review of that receipt and every source pin, the single experimental
evaluation was submitted as **Torch 18754308**. `launch.json`
records the changes, budget and exact source pins. No additional calibration,
transition or automatic promotion follows. Latest/best summaries on Torch
explicitly report no completed experimental point until one passes.

The job failed after 28m20s, following ten stationary evaluations and native
export. Completed fertility reached its normalization target (2.100022 against
2.1), but the additional `audit_extra_cache` check rejected probability-menu
sums or reached invalid menus. The native export is provisional and the full
experiment did not pass. `audit_saved_cache.py` inspects the saved arrays on
Torch, splitting zero menus from nonzero normalization discrepancies, with
no model solve or change to the failed point. Job 18757242 owns this inspection.
Preserve the failed source and gates; no automatic retry or promotion.

## Saved-point diagnosis and provisional complete tables

Torch 18757242 inspected the failed checkpoint without solving the model.
All three caches have the correct shape, finite probabilities in [0,1], exact
equality and independent ownership. Twenty-two menus fail the 1e-12 sum check;
maximum error is 9.575e-08. Their reached mass is exactly zero, as is reached
mass at zero menus. Largest error among materially reached menus is 1.332e-14.
Native probability calculation `exp(a-logsumexp(a))` loses precision at large
common value levels. A shifted-exponential calculation is being tested in
`test_probability_stability.py`; the original failed sources and threshold
remain unchanged. Native export is not a passed full experiment.

Torch 18757604 passed that regression and all eight existing solver tests in
three seconds, with zero model solves. Shifted-exponential probabilities reduce
the synthetic error from 2.845e-08 to 2.220e-16. Inclusive values are bitwise
unchanged; the dead-menu mask and 1e-12 tolerance are unchanged. The generated
candidate source remains on Torch, pinned by `probability_repair_v1/receipt.json`.
The correction has not been installed in a model run. A full corrected-point
verification remains necessary before accepting or interpreting the experiment.

All ten reference coordinates were held fixed; only child benefit was adjusted
to its completed-fertility normalization target. Zero coordinates were searched
in this experiment. The raw parameter status column retains historical source
labels because the worker stopped before its final reporting rewrite. The
explicit statuses below describe this experiment. All 17 native plots are
retained under `run_v1/worker/case/standard_diagnostics/`; their presence does
not override the failed additional check.

These are provisional outputs, not an adopted calibration or a verified
economic conclusion. Loss 394.225; the comparison uses the unchanged original
target system and weights.

### Targeted rows

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Childlessness, ages 40–44 | 0.198 | 0.242 | 0.043 | 35532.304 | 66.937 |
| Exactly one child among mothers, ages 40–44 | 0.214 | 0.190 | -0.024 | 26952.821 | 15.092 |
| Mean age at first birth | 25.976 | 27.399 | 1.423 | 139.828 | 283.043 |
| Wealth / annual earnings | 6.927 | 6.314 | -0.613 | 7.595 | 2.850 |
| Annual bequests / wealth | 0.007 | 0.007 | -2.244e-04 | 5165289.256 | 0.260 |
| Mean rooms | 5.729 | 5.825 | 0.096 | 128.021 | 1.173 |
| Ownership, ages 30–55 | 0.676 | 0.655 | -0.021 | 2339.362 | 1.032 |
| First-birth rooms response | 1.465 | 1.729 | 0.264 | 137.565 | 9.623 |
| Recent-parent ownership gap | 0.128 | 0.142 | 0.014 | 27055.823 | 5.625 |
| Children by age 25 | 0.810 | 0.516 | -0.293 | 100.000 | 8.590 |

### Untargeted checks

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| First births at age 30+ | 0.249 | 0.317 | 0.067 | 0.000 | 0.000 |
| Older wealth/income, p90/p50 | 3.516 | 3.070 | -0.445 | 0.000 | 0.000 |
| Rooms gap: 3+ versus 1–2 resident children | 0.385 | 0.286 | -0.100 | 0.000 | 0.000 |

### Normalization target

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Completed fertility (normalization) | 2.100 | 2.100 | 2.223e-05 | — | — |

### Parameter estimates and restrictions

| Parameter | Estimate | Lower | Upper | Near bound | Status in experiment |
|---|---:|---:|---:|---|---|
| H0 | 6.294 | 0.200 | 80.000 | False | Held at reference estimate |
| beta_annual | 0.963 | 0.940 | 0.990 | False | Held at reference estimate |
| chi | 1.094 | 0.100 | 5.000 | False | Held at reference estimate |
| first_birth_fixed_cost | 0.621 | 0.000 | 8.000 | False | Held at reference estimate |
| kappa_fert | 0.176 | 0.020 | 50.000 | True | Held at reference estimate |
| kappa_fert_continuation | 0.332 | 0.020 | 50.000 | True | Held at reference estimate |
| theta0 | 0.125 | 0.000 | 8.000 | False | Held at reference estimate |
| delta_alpha_jump | 0.135 | 0.000 | 0.250 | False | Held at reference estimate |
| child_benefit_curvature | 0.061 | 0.000 | 0.800 | False | Held at reference estimate |
| tenure_choice_kappa | 0.012 | 0.001 | 0.100 | False | Held at reference estimate |
| psi_child | 0.090 | — | — | — | normalized to completed fertility 2.1 |
| child_benefit_CRRA_coefficient | 0.085 | — | — | — | derived from normalized one-child benefit |
| theta1 | 0.008 | — | — | — | fixed external restriction |
| sigma | 2.000 | — | — | — | fixed |
| alpha_cons | 0.733 | — | — | — | fixed CEX childless expenditure share |
| delta_alpha | 0.000 | — | — | — | fixed zero later-child loading |
| h_P | 0.000 | — | — | — | no housing floor |
| utility_reference_rent | 0.110 | — | — | — | fixed substantive utility normalization |
| q_annual | 0.020 | — | — | — | author-retained 2% annual real rate |
| financed_share | 0.800 | — | — | — | inherited credit contract |
| housing_supply_elasticity | 0.630 | — | — | — | fixed provisional external mapping |
| payroll_tax | 0.080 | — | — | — | derived from adopted pension ratio |
| pension_period | 0.918 | — | — | — | balanced PAYGO |
| annual_depreciation | 0.014 | — | — | — | adopted |
| period_depreciation | 0.055 | — | — | — | compounded |
| annual_property_tax | 0.011 | — | — | — | adopted |
| period_property_tax | 0.042 | — | — | — | linear period convention |
| selling_cost | 0.060 | — | — | — | retained |
| rental_cap | 6.000 | — | — | — | retained provisional |
| wealth_grid_nodes | 160.000 | — | — | — | retained exact grid |
| income_states | 15.000 | — | — | — | retained B15 |
