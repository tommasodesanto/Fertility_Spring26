# Perfect-foresight result reconciliation — September 9, 2026

Bounded read-only audit of canonical status, saved result receipts, and the frozen September 7 worktree. No model solve, numerical-source edit, target change, or endorsement of a fallback.

## What was actually achieved

**The perfect-foresight work was implemented and produced a fully audited paired H128 result. It was not merely a plan.** The September 1 result used the coherent person-demography law and positive tenure taste dispersion κ = 0.010017488787185433, imposed after 2023 while older E5f parameters were held fixed. Both rebated 1% and rebated 2% property-tax paths passed the unchanged market (2e-4), fiscal (2.5e-5), history, distribution, probability, accounting, and terminal-state gates. This was explicitly a sensitivity, not a production recalibration or the current simple-fertility-nest calibration.

Canonical evidence: `CALIBRATION_STATUS.md:2479–2501`. Saved audited mechanism receipt: `output/model/e5f_post2023_kappa_m5_property_tax_mechanism_decomposition_20260901/summary.json`, status `complete_audited_h128_property_tax_mechanism_decomposition`. Its README records the collector and inputs. Input folders:

- `output/model/e5f_pf_person_policy_rebated-tax1-baseline_20260831e_post2023_kappa_m5_h128_stage03_convergence/`
- `output/model/e5f_pf_person_policy_rebated-tax2-reform_20260831f_post2023_kappa_m5_h128_stage04_convergence/`

Both input summaries retain the conservative label `complete_unpromoted_person_demography_policy_path`; passing numerical gates is distinct from promotion. Baseline market/fiscal maxima: 1.4941657914915922e-4 / 1.413057030391629e-5. Reform: 1.156722366450562e-4 / 1.6510816389980754e-5. The baseline saved terminal audit explicitly reports `all_checks_pass: true`, terminal year 2535, and passing person, head, distribution, price, rent, transfer, and fertility-intercept checks. The recorded source provenance pins the calendar Bellman solver, PF solver, selected August 18 report, and stationary source.

## The latest attempt before the new nest

The September 3 pair used the then-new η = 0.63 and κ = 0.005 joint calibration. The reform passed H128; the baseline did not. The baseline best market/fiscal maxima were 3.4854883947494e-4 / 3.6913691039766605e-5, both above tolerance, concentrated at calendar 2359. One authorized continuation reproduced the seed but six updates did not contract: the last market/fiscal maxima were 7.659065437571115e-4 / 7.335046616901986e-5. All feasibility and accounting checks passing does not cure that equilibrium failure.

Evidence: `CALIBRATION_STATUS.md:2440–2477`; saved comparison `output/model/e5f_pf_person_policy_comparison_20260903b_eta063_kappa005_h128_runtime_diagnostic/summary.json`, status `complete_unaccepted_person_demography_policy_comparison`, accepted false, baseline path_converged false, reform true, both terminal roots converged. Seven continuation snapshots are preserved under `output/model/e5f_h128_baseline_stage02_snapshots_20260903/`. September 4 explicitly retained the unresolved baseline solver problem and did not certify a new H128 path (`CALIBRATION_STATUS.md:2329–2331`). The diagnosed remaining problem was the coupled price/transfer fixed point, not insufficient horizon; H128 retains approximately 0.62% of the slow demographic mode.

## Three distinct transition systems

1. **August 26 birth-queue PF:** calendar foresight existed, but the old birth-to-household queue demographic closure was superseded. H256 paths were not accepted: baseline market/fiscal residuals failed; reform fiscal failed; terminal population gaps remained large. Evidence: `CALIBRATION_STATUS.md:2878` onward and `output/model/e5f_perfect_foresight_rebated_tax_h160_h256_exact_early_diagnostic_20260826q3/`.
2. **August 26–September 3 coherent-person PF:** annual age/sex person accounting, survival, migration and headship map persons to household-head distributions, alongside calendar-time household decisions and endogenous prices/transfers. This is the implementation behind the verified September 1 pair and later September 3 convergence issue.
3. **Current simple-fertility-nest 2023–2063 diagnostic:** a per-date temporary-equilibrium calculation with the older household-entry birth queue. It does not use the preceding full calendar PF stack. The driver itself states this in `tmp/e5f_fertility_nest_compute_20260907a/code/model/tools/run_e5f_simple_fertility_tax_transition.py:28`, and advances through `baseline.advance_from_evaluation` at line240. No evidence in this bounded record establishes author approval to replace PF with this diagnostic. It must not be presented as completing the earlier PF objective.

## Existing implementation and concrete adaptation gap

All source anchors below are relative to the frozen worktree `tmp/e5f_fertility_nest_compute_20260907a/`:

- `code/model/tools/run_e5f_perfect_foresight_person_demography.py:693` calls `pf.backward_value_path` over prices, rents, transfers and fertility intercepts; lines740–768 solve/reproduce each date with the next-calendar-date continuation value; lines770–812 forward the household and person blocks. This is an implemented backward/forward transition chain.
- `code/model/tools/run_e5f_perfect_foresight_transition.py:233–251` constructs rents using the next asset price (`carrying_factor * current_price - next_price`), so expected capital gains enter. `solve_date_policy` at416–433 supplies `continuation_V`; `backward_value_path` at436–477 propagates calendar continuation values.
- The current Bellman solver already accepts the calendar continuation argument (`code/model/intergen_eqscale_seq_optimized/solver.py:2350`) and contains the simple-fertility-nest branch. This infrastructure need not be recreated from scratch.
- **Verified integration blocker:** the older PF `policy_from_objects` at380–413 builds a `calendar.PolicyBundle` without its `joint_choice`. Current joint-nested forward accounting requires that exact owned object: `run_dynamic_population_transition.py:93–106` raises if it is absent. Thus simply switching the new nest on in the old PF wrapper is not a valid ready-to-run adaptation.

The remaining work is to connect and preserve each dated joint-choice object, initialize the selected new calibration with an explicitly reconciled person/headship contract, regenerate terminal endpoints under its parameters and supply rule, and pass an exact backward/forward reproduction smoke before solving/certifying the coupled full path. The old converged paths are useful historical evidence and possible controlled seeds, not certification of the new nest. This read-only audit does not establish that the joint-choice omission is the only integration issue.

## Interpretation for the author

The honest status is: **we did build and validate perfect foresight in the earlier model; the latest pre-nest baseline had a documented convergence problem; the recent new-nest policy exercise temporarily used a different, simpler transition calculation without carrying that PF work through.** Calling the current exercise the completed final transition would overstate it. The next technical decision should concern adapting the existing PF machinery and resolving its coupled convergence issue, with the diagnostic results kept clearly separate.
