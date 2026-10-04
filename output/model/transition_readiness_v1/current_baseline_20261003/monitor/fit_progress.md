# Completed transition trial — not a fitted shock

Local worker 7 completed both diagnostic horizons (24 and 32 dates), with convergence, exact replay and horizon checks passing. Its final fertility gap is above the required 0.005 tolerance. The scalar search is continuing; this is not an accepted estimate or a production policy result.

| Birth years | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| 2008–2011 | 1.975 | 2.071 | +0.097 | 0 | 0.000 |
| 2012–2015 | 1.861 | 2.074 | +0.213 | 0 | 0.000 |
| 2016–2019 | 1.755 | 2.075 | +0.319 | 0 | 0.000 |
| 2020–2023 | 1.646 | 2.075 | +0.430 | 1 | 0.185 |

| Parameter | Tested value | Lower bound | Upper bound | Near bound |
|---|---:|---:|---:|---|
| psi_child | 0.175342 | 0.001789 | 0.357841 | No |

This tested permanent level is 0.98 times the unchanged baseline. It is a numerical trial, not the final estimated parameter. Earlier windows carry zero fitting weight and remain validation rows. Full precision and the complete gate receipt are in [fit_progress.json](fit_progress.json).

Baseline parameters remain fixed: [all 31 parameters](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/experiments/birth_count_choice/estate_a_v1/single/cases/20261003T212605039706Z_a739edc3/parameters.csv).

Evidence: [completed candidate](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/local_workers_v7/worker_7/run/candidate_0001/complete.json) and [pinned plan](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/plans_v7/fit_plan.json). Terminal failures and outstanding estate/policy closure remain disclosed in those receipts; this diagnostic does not establish 104/128-period certification.
