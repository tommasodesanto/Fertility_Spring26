# Successive unexpected permanent preference changes

Prepared September12,2026. Author requested code preparation while separately checking an empirical target. **No calibration, endpoint solve, transition solve, cluster launch, or target change is authorized by this preparation.** The announced-history implementation is untouched.

## Economic experiment

At each four-year historical decision date, the current preference parameter is observed and believed permanent. Solve a complete expected equilibrium path of housing prices and balanced PAYGO pensions conditional on that parameter. Implement only its first period; the next surprise arrives with the resulting distribution inherited. Future realized preference levels never enter the current household solve. Preferences are constant within each forecast, but prices and pensions generally are not.

This is a sequence of deterministic conditional forecasts with expectation revisions, not a rational-expectations aggregate-risk model. The retained observed household totals/age marginals through2023 and fixed2023 person-cohort anchor remain known conditioning inputs. Only information about fertility preferences changes.

Liquid wealth/debt and the owned housing product carry separately. New house prices revalue housing through current-price transaction maps; do not add a capital-gain shift to the liquid-wealth grid. Existing occupied-state feasibility checks remain active. Each vintage retains its expected next price and rent; realized subsequent prices must not be used to reconstruct an earlier rent.

## Code and scope

- `code/model/tools/e5f_successive_surprises.py`: experimental adapter and controller.
- `code/model/tools/test_e5f_successive_surprises.py`: local tests using small artificial household responses and the existing price/pension root.
- Runtime dependencies are the approved utility/PAYGO snapshot in `tmp/e5f_matched_pf/code/model`; production model modules are not overwritten.

`solve_surprise` takes an `InheritedState(year, households)`, a scalar current `psi`, the original approved initial-state reference, a verified terminal solved at that SAME psi, frozen demographic primitives, explicit price/pension guesses and numerical budgets. It calls the existing backward/forward household routines and joint fiscal/housing root. A validation-only synthetic preference line reuses the old wrapper's initial/terminal primitive checks; that line is never used in a household solve. All solved forecast preferences are constant.

`evaluate_forecast` joins the remaining historical dates to the existing person-demography tail from2023. Reindexing the historical bridge preserves its original2007 normalization. No distribution is reset to the stationary benchmark when another shock arrives.

`first_period_state` uses the accepted forecast's price at date1 and value at date1 as the one-period boundary. The original stationary terminal is NOT used for this replay. The replay must reproduce the first-period value and outcome/fiscal entries. Historical state carries both raw and adjusted birth queues. The2019-to2023 step hands the historical distribution to the retained2023 person anchor; the2023-to2027 step carries annual person/head cohorts.

`run_sequence` accepts explicit `(year, psi)` pairs, a `solve_episode(inherited=..., psi=...)` callback, and a persistence callback. The solver callback receives no future shock vector. `persist_episode` saves each full expected path as a CSV, its root receipt, the single implemented observation, and an exactly reloaded next-state checkpoint. Nonfinite scores from rejected trials are saved as labeled strings. Output directories cannot be overwritten.

Current surprise dates are2007,2011,2015,2019,2023. A terminal2027 state can be saved; additional surprises after2023 are explicitly unsupported by this adapter because the underlying person evaluator has a2023 calendar anchor. Four shocks at2007–2019 correspond to the existing four observed birth windows2008–2011 through2020–2023. A surprise at2023 concerns2024–2027 births and cannot be fitted to the preceding window.

## Calling contract once the target discussion is closed

The experiment is an API, not an automatic launcher. The run owner must:

1. Pin the approved initial checkpoint, all scientific sources, empirical inputs, complete target/weight fingerprint, and numerical controls. Do not select between the unrestricted and capped-beta candidates implicitly.
2. Build the original approved initial-state reference once. Keep its original supply curve, replacement2.1 and birth conversion1/2.1.
3. For each trial psi, solve or load an independently verified balanced terminal with identical structural parameters, supply, demography and this psi. The endpoint from the final historical shock is not valid for earlier vintages unless psi is identical.
4. Pass explicit dated price/pension guesses, terminal household audits and a monotonic deadline to `solve_surprise`. Wrap the process in the existing watchdog. Every mapping and best-so-far receipt must be saved with the callback.
5. Retain all forecast vintages and the stable17 diagnostic graphs through the existing observer/report writer. Generate diagnostic graphs at the accepted period before deleting its evaluation; do not infer a stationary state from a path checkpoint.
6. Run the native exact-loop smoke, then horizon extensions, before fitting shocks or reporting numerical effects. The local tests below are not a model-run certificate.

Each accepted diagnostic stage requires the unchanged housing2e-4, fiscal1e-6, replay2e-10 and household gates, plus the existing terminal distance checks and an explicit pension-tail tolerance no looser than1%. Failure persists the failed receipt and does not advance the realized state. `horizon_verified` and `production_eligible` remain false: endpoint proximity alone does not establish insensitivity to a longer horizon.

## Sizing before any launch

For S surprises, N dates per forecast and K maximum root mappings, the path budget is at most2SNK dated Bellman calls, plus up to2S first-period replay calls and separate terminal solves. Four surprises at N=28,K=8 mean1792 dated Bellman calls plus at most8 replay calls. The prior28-date mapping took roughly12–17 minutes: exhausting that budget implies roughly6.4–9.1 hours for one four-shock chain, excluding endpoints and queue time. Early convergence can lower this; a fresh native mapping must establish actual timing. Shock fitting multiplies this by the number of trial forecasts. Future forecast vintages depend on accepted previous states; independent terminal calculations and candidate shocks conditional on the same inherited state can be parallelized.

No empirical objective or optimizer is chosen here while the target discussion is open. A full stochastic process is also outside this preparation.

## Verification

From the project root:

```sh
PYTHONPATH=code/model/tools:tmp/e5f_matched_pf/code/model/tools:tmp/e5f_matched_pf/code/model /opt/anaconda3/bin/python -B -m unittest test_e5f_successive_surprises -q
```

Checks cover all five restart dates, exactly one2023 transition, continuation-value indexing, constant expected preferences, future-surprise exclusion, original demographic scale, both birth queues, existing joint-root updates and fresh final selection, household-audit rejection, unchanged2.1 normalization, rejection of a terminal for the wrong psi, failure stopping, one-period replay boundaries, accepted-state advancement, checkpoint reload, and refusal to overwrite outputs. An independent read-only review checked these interfaces and liquid-wealth accounting against the active source. Native model execution is intentionally deferred.
