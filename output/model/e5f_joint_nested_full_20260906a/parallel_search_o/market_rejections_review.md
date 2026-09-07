Read-only diagnosis complete. The three rejections certify only that the dated market-clearing routine could not meet the unchanged \(2\times10^{-4}\) residual gate after its fallback search; they do not establish market nonexistence or a coding defect.

| Case | Last completed date | Failing period | Residual evidence | Retry behavior |
|---|---:|---:|---|---|
| 023 | \(t=00\) | \(t=01\) | Final fallback residual \(4.361\times10^{-3}\), 21.8× gate | Strict \(5\times10^{-5}\), 18-iteration attempt failed; fallback uses \(2\times10^{-4}\), 60 bisections and also failed. Initial-attempt residual was not retained. |
| 024 | \(t=00\) | \(t=01\) | Strict-attempt residual \(3.226\times10^{-3}\); fallback \(1.286\times10^{-3}\), 6.43× gate | Same 18-then-60 retry. Fallback improved the recorded residual but remained outside the unchanged gate. |
| 027 | \(t=02\) | \(t=03\) | Final fallback residual \(2.696\times10^{-3}\), 13.5× gate | Same 18-then-60 retry. Initial-attempt residual was not retained. |

Evidence: the task receipts record completed periods 1, 1, and 3 respectively, and the logs show the corresponding successful \(t\)-rows. Receipts are under `.../initial_population/task_{023,024,027}/adapter_failure.json` and `heartbeat.json`; the aggregate ledger is [rejects_ledger.csv](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/parallel_search_o/output/model/joint_nested_overnight/search/rejects_ledger.csv).

The solver starts at the prior-price guess, expands a bracket by factor 1.18 for at most 18 expansions, then bisects a signed demand-minus-supply bracket. A candidate returns only when its relative residual is within tolerance; otherwise, after `max_iter`, it evaluates the final midpoint and raises the recorded error. See [run_dynamic_population_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_dynamic_population_transition.py:530) and especially lines 564–615. Thus, each receipt’s final error establishes a successful sign bracket followed by exhausted bisection at the relevant tolerance—not a failed bracket.

The open-population path first applies a \(5\times10^{-5}\) target with `min(market_max_iter,18)`, and catches only this exact non-clearing error to retry the same date at the declared \(2\times10^{-4}\) tolerance and 60 iterations. The code explicitly notes that deterministic tenure choices on a finite wealth grid can make off-stationary demand locally discontinuous. See [run_e5f_open_population_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_open_population_transition.py:1152) and lines 1166–1188.

Established versus conjecture:

- Established: none of these are failed-bracket events; no receipt shows nonfinite evaluations; all three exhausted the fallback bisection without satisfying the unchanged gate.
- Established: the saved artifacts do not retain bracket endpoints, midpoint prices, signed excess demand, or per-iteration residuals. Therefore they cannot distinguish a jump/discontinuity from slow convergence within a smooth schedule, a narrow oscillatory grid effect, or another numerical mechanism.
- Conjecture only: the documented finite-grid deterministic-choice discontinuity is a plausible explanation, but these receipts do not prove it is operative in any particular case. They also provide no evidence of economic market nonexistence or a coding error.

Smallest resolving future diagnostic: save a compact trace for the failed date and fallback only—initial bracket endpoints plus each bisection midpoint’s price, signed excess demand, and relative residual. That single trace would distinguish a discontinuous jump across zero from smooth residual decay/max-iteration exhaustion, without changing gates or targets.

For context only, the one-hour-cap rejects stopped after actual completed progress of \(t=01\) (case 005) and \(t=02\) (case 010); neither is evidence about market clearing.

Lead reconciliation: the cited bisection and retry code and saved failure records were inspected. The report correctly distinguishes a successful sign bracket from failure to meet the residual gate. Sixty bisections at ordinary housing-price magnitudes normally reach floating-point price resolution, so merely adding iterations is not a supported remedy. A saved signed-demand trace is needed to establish whether a discontinuity or finite-precision sensitivity is responsible. No code or numerical gate was changed.
