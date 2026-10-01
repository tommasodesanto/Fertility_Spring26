# Three-hour normalized calibration restarts from verified v1 candidates

This is a new set of 24 Nelder–Mead restarts from six postcheck-verified v1
candidates, with one exact center and three small deterministic joint neighbors
per center. It is not a saved optimizer continuation: v1 did not serialize its
final simplex, and an interrupted SciPy optimizer cannot be reconstructed from
selected-point receipts. V2 saves final simplex coordinates and objectives only
when SciPy returns normally; budget interruptions make no reconstruction claim.

The shared native gate uses verified v1 chain 20, case `0028_nm`: original-base
loss 30.371887956158005, price 0.7167873404451099, derived $H_0$
6.851575289344519 and $\psi=0.17198899419542374$. Its complete ROOT/REPEAT
verification is in `deployment/best_v1_candidate/`.

All economics and numerical equilibrium gates are unchanged from
`normalized_calibration_v1`: actual national benchmark population $N_0=1$,
derived internally calibrated $H_0=h(p)/(up/\bar r)^\xi$ with original accepted-
root bounds $[0.2,80]$, free searched $\psi$, physical first-child housing floor,
$D=0$, nonnegative mean-preserving entry, 120×9 grid and owned rooms 2/4/6/8/10.
The same ten scored targets, four other reported rows, all original weights,
coordinate bounds, completed-fertility normalization and topcode correction are
retained. There is no price/rent target or external $H_0$ anchor. The housing
supply coefficient is derived at actual native demand before reports and
checkpoints; the mean occupied rooms target remains scored rather than forced.

The numerical changes are starts, simplex size, and search budget. Each chain
has 10,800 seconds from actual launcher start, at most 150 objective calls,
the unchanged 900-second final native reserve, one CPU/24 GiB/one thread and
unchanged native deadline/32-lifecycle-call guards. Maximum production work is
3600 exploratory GEs plus 24 selected native postchecks. V1 completed about
35–47 objective calls in 92–93 minutes. With unchanged guards, the new window
allows roughly 153 minutes of exploration, suggesting about 58–78 calls per
chain (roughly 1400–1900 across all chains). The 150-call limit is a hard cap,
not an expected completed count. No retries, budget extensions or tolerance
relaxation are introduced.

`plan.json` records the six authenticated centers, every starting vector,
selected starting prices, RNG seed, physical perturbation scales and clipping.
The physical initial simplex steps are 25% of v1's steps. Original bounds and
weights are pinned and checked. `center/` contains the exact global-best v1
ROOT tables and verification for the shared native gate.

Verification commands:

```sh
PYTHONDONTWRITEBYTECODE=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/normalized_calibration_v2/test_configuration.py
python FULL_DEPLOYED_PACKET/run_psi.py --chain 0 --out NEW_NATIVE_GATE_DIRECTORY --deadline-epoch ACTUAL_START_PLUS_10800 --smoke-only
```

The first command exercises all 24 actual mocked 150-call loops with zero model
solves, plus the existing normalization, source-variant, mass, bounds, reserve
and fatal checks. The shared Torch gate uses the exact best verified v1 point
and its selected price, requires all 14 target values and all 31 parameter
estimates including $H_0$ to match its authenticated ROOT tables within
$10^{-10}$, and retains the full native ROOT/REPEAT arrays, tables, closure and
17 plot-hash checks. Searches use `--fast-objective`; selected postchecks use
`--verify-only SEARCH_DIRECTORY/search_completed.json`. Source and gate results
must pass before dependent production starts. No native solve or submission is
performed by preparation; `deployment/` owns collection, staging and launch.
