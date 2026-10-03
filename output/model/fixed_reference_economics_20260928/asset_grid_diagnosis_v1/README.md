# Asset-grid diagnosis and model control

## Verified facts

- Reference: authenticated selected soft-financing checkpoint **chain16/case0046**, original transaction timing, at fixed price **0.7266387868818555**. The canonical saved-array SHA-256 is `e154b3430b970390bb4db43eddd26919a11e1f9e15bcba16a37394f31700eab0`. The baseline replay matches all 11 checked native arrays exactly.
- Financial assets \(b\) are net of financial debt and exclude the separately tracked house. The model already expresses these assets in mean annual gross-earnings units. The 120-node grid spans \([-12,3000]\); about 99.8665% of pre-tenure household mass lies in \([-6.4,33.6605536290589]\).
- Both the compiled Bellman interpolation (`engine/kernels.py:interp_scalar`) and the forward maps (`engine/utils.py:interp_indices`) clip queries to grid endpoints. Extending the maximum to 6000 changes high-wealth policies but changes **none** of the tested tenure-averaged policies at originally occupied states or population aggregates.
- Spacing matters. Refinement strictly above the exact core endpoint changes mean inherited financial assets by 0.3155%. Refining the core to 214 and then 402 total nodes changes this mean by 0.7627% and 0.9619% relative to the baseline. Ownership is 66.6495%, 66.4100%, and 66.3334%, respectively.
- All economic inputs, income states, entry wealth atoms and probabilities, purchase rule, financed share, prices and supply coefficient are held fixed. No credit, transfer, demographic or preference fallback was introduced. These are numerical fixed-price tests; no grid or calibration has been adopted.

## Complete population comparison

Consumption is a four-year flow. Financial assets and changes are in mean annual gross-earnings units; ownership is percent. Asset change is mean \(b'\) minus mean inherited \(b\), including housing transaction cash flows; it is not national-account saving.

| Numerical case | Nodes | Mean consumption | Mean inherited assets | Mean next assets | Mean asset change | Mean rooms | Ownership (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Baseline | 120 | 2.403257705 | 1.770789328 | 1.681880949 | -0.088908379 | 5.957167757 | 66.649532706 |
| Maximum extended to 6000 | 124 | 2.403257705 | 1.770789328 | 1.681880949 | -0.088908379 | 5.957167757 | 66.649532706 |
| Tail plus boundary interval | 138 | 2.403851446 | 1.781483275 | 1.692842841 | -0.088640433 | 5.955843780 | 66.626637346 |
| Strict tail refinement | 137 | 2.403574277 | 1.776376819 | 1.687612788 | -0.088764031 | 5.956722648 | 66.643188640 |
| Core intervals halved | 214 | 2.403338508 | 1.784296006 | 1.695380898 | -0.088915108 | 5.953624628 | 66.410040447 |
| Core intervals quartered | 402 | 2.403320328 | 1.787823273 | 1.698840272 | -0.088983001 | 5.952842393 | 66.333384322 |

Changes shrink between the two core refinements: inherited assets change another 0.1977%, and ownership another -0.07666 percentage points. Individual policies are less stable than these means: the largest consumption difference between the 214- and 402-node grids is 0.56166 at a state with at least \(10^{-12}\) of original population mass; the population-weighted mean absolute difference is 0.006748. The tests establish measurable resolution error, **not** a full convergence certificate for individual policies, calibration targets or equilibrium prices.

## Evidence and reproduction

- [Initial four-case results](summary.json), [strict-tail and second-core results](followup_run2/summary.json). Each of the six solutions has the unchanged standard 17-figure diagnostic packet in its `diagnostics/` subfolder. Policy comparisons average tenure choices at inherited asset nodes, using the native transaction map and separate owner-stayer arrays, and weight states by the baseline beginning distribution.
- [Supplemental raw consumption comparison](supplemental_grid_policy_comparison.png) shows all six cases at age 30, income state 5, for the childless renter branch, with the occupied range and upper tail shown separately. Raw initialized entries at infeasible nodes are not occupied households.
- The initial tail test used the rounded cutoff 33.66055 and included the interval ending at 33.6605536290589. The follow-up isolates intervals whose **lower** endpoint is at or above that exact node. The original result is retained and explicitly labeled above.
- The first follow-up failed before solving because the minimal saved archive omitted an equivalent wealth-distribution array. Recovery reconstructs that alias from the pinned native copy semantics and authenticates the cached baseline against all 11 canonical arrays, prior replay receipt, reference price and supply coefficient. The reproducible input fingerprint excludes only the temporary evidence-directory path; that path is not an economic input. Failure evidence remains in `followup_run.log`; successful recovery is `followup_run2.log`.
- The full initial packet took 40.7 seconds; the follow-up took 54.9 seconds. Six solves were sequential on one numerical thread with 10-minute packet and 180-second solve limits. Solve times were approximately 6.3, 6.5, 7.5, 14.7, 7.3 and 43.1 seconds. macOS did not expose current virtual-memory size, so the requested 24-GiB address-space cap was recorded but not enforced.

From the repository root:

```bash
MPLBACKEND=Agg code/model/.venv/bin/python code/model/tools/diagnose_playground_asset_grid.py
MPLBACKEND=Agg code/model/.venv/bin/python code/model/tools/diagnose_playground_asset_grid.py --followup
```

The driver is [diagnose_playground_asset_grid.py](../../../../code/model/tools/diagnose_playground_asset_grid.py). Its commands regenerate this diagnostic packet. Input identity and estimates are pinned by [soft_selected.json](../soft_timing_review_v1/soft_selected.json); all economic estimates remain unchanged in these grid tests.

## Interactive model control

The [Python playground guide](../../../../code/model/tools/MODEL_PLAYGROUND.md) documents editable parameters, fixed-price solves, raw policy plots, population aggregates and baseline comparisons. The explorer at <http://127.0.0.1:8765/> now separates conditional policy tenure from inherited tenure and provides actual asset values, an equivalent earnings-unit label, and node indices. Its whole-population section is independent of household slice selectors and uses the native realized distribution with separate owner-stayer policies.

The [beta verification](../model_control_beta_check_v1/verification.json) raises annual beta from 0.9671931106 to 0.98 with the same fixed price. One solve plus plots took 7.43 seconds; the four-year discount factor is 0.92236816. This is a parameter experiment, not a recalibrated equilibrium.
