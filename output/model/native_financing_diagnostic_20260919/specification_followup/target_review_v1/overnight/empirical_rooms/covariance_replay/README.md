# First-birth rooms covariance replay staging

This folder contains a review-only clone of the frozen first-birth rooms
regression. The real PSID replay has **not** been launched. The original source
hash is `f588d6946d772155b75cfdf272b8ae46e5e362bc3b6f977fccd16aa57d18f36c`;
`replay_v1.diff` shows every change from that source.

The clone preserves the source data, sample construction, room alignment and
coding, weights, controls, fixed effects, estimator, and target calculation. It
changes only the output directory and Stata's processor cap (8 to 4). Diagnostic
flags count aligned room values 0, 9, 98, and 99 before room/education deletion,
after room/education deletion but before later cohort restrictions, and in `e(sample)`; they do not filter or refit.
The minimum nonmissing room count is recorded at each stage, including zero if
present.

`eventstudyinteract.ado` constructs `e(b_interact)` as cohort by event time. Its
`e(V_interact)` is only a cohort-by-event marginal-variance array: the ado takes
`diagonal(e(V))` before reshaping. The full clustered interaction covariance is
the leading interaction block of `e(V)`. The ado creates that block in event-time
outer, cohort inner order (`rel_time_list`, then `cohort_list`). Each exported numeric scalar uses 17 significant digits. Before export, the clone asserts that every mapped cohort/event coefficient and `e(V_interact)` marginal variance agrees with the corresponding entry in `e(b)` and the diagonal of `e(V)` to `1e-14`. The covariance CSV
uses compact numeric row/column indices; the coefficient CSV maps each index to
cohort, event time, and Stata interaction-variable name. `e(b_iw)` and its
full `e(V_iw)` are separately exported for the existing interaction-weighted
aggregation. `e(V_iw)` includes the package's cohort-share variance addition.
Fixed-common-cohort contrast uncertainty computed from the interaction block
will be conditional on the chosen cohort weights; this packet does not claim
sampling uncertainty for externally fixed weights.

The staged outputs contain aggregate cohort/event support and coefficient
matrices only; never-treated controls are aggregated with event time missing
because event time is undefined for that group. No person-level PSID records are exported. `eventstudyinteract_replay.ster`
will be written only by the future reviewed run. The gate CSV reports (without
stopping or changing the regression) whether target `0.7202462623815278`, SE
`0.0852600513385958`, 49,457 observations, and 4,112 people reproduce. Target and standard-error reproduction gates use absolute tolerance `1e-9`, matching the nine-decimal result display in the original fit log; the full-precision expected and observed values are retained, so absolute gaps can be calculated in `reproduction_gates.csv`. Observation and person counts require exact equality. These are empirical reproduction checks and do not change any model gate.
A failed gate is reported as a discrepancy, with no automatic retry or relaxed
tolerance.

Preflight found Stata 17.0 and the installed `eventstudyinteract` and `svmat2`
ado files. Their paths and SHA-256 hashes, plus the frozen-source and staged-file
hashes, are in `source_package_manifest.json`. A synthetic Stata-only smoke passed
for four named cohort/event coefficients and all 16 cells of their full covariance
export. It used no PSID records. The smoke script is
`output/covariance_export_smoke.do`; its synthetic CSVs are retained as the small
smoke receipt.

The original local runtime was 511.832 seconds with eight Stata processors. A
single reviewed replay at the staged four-processor cap has a 20-minute wall-clock
budget, with no retries or specification changes. Expected command after lead
line-by-line approval:

```bash
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -q -b /absolute/path/to/replay_v1.do
```
