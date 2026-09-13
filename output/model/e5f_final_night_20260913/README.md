# Authorized final-night work

Latest evidence: September13, approximately10:00UTC. Deadline18:00UTC.
The cluster runs and collector are independent of the laptop. Idle system sleep
is prevented locally until the deadline. Last account check:70% weekly remaining.

## Numerical repair under verification

Trace17642567 reproduced the native2011 mass failure. The loss comes from
float32 tenure-probability row sums. `probability_mass_native_verification.json`
records an exact old-output replay and two identical corrected outputs; relative
mass error falls from1.067e-8 to5.69e-13, below the unchanged1e-8 gate. Both the
stationary Markov KFE and the transition Markov KFE now normalize these rows in
float64. The first isolated corrected initial job17651740 reached its root but
failed stationary nesting because it only corrected the transition copy. A new
two-repetition full initial check is being prepared. No corrected calibration
is accepted yet; all original source folders and rejected outputs are retained.

`fiscal_polish_Aplus_6_summary.json` records probe17650020: the stalled A+ forecast
passes every market/fiscal/replay check in three evaluations (245seconds), with
all seven asset prices held at the prior market-admissible values. Its2007 TFR
is1.9788233781 against1.974875, within0.005. This is one forecast without historical
carry. An opt-in automatic early polish passes13 wiring tests and preserves the
full residual vector,24-evaluation total cap, deadline and exact replay. A full
native combined-controller run is still required before adoption.

Remote batch: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/`.
The requested plan is `docs/model/e5f_two_closure_overnight_plan.md`.

## Current work

| Jobs | Work | Evidence/status |
|---|---|---|
|17613033|A0/A+ six-date histories and rebated policy|A0 first window accepted and carried into2011; A+ fitting|
|17613034|A0/A+ 24-date histories and rebated policy|Equilibrium iterations in progress|
|17632922|A0/A+ 100-date histories and rebated policy|Submitted; queue counts against deadline|
|17635971/17635972|One-forecast Jacobian warm-start comparison|A0 passed; A+ rejected by mass gate|
|17641919|Native saved-state restart validation|Every inherited state field exactly verified, no solve|
|17642567|Instrumented reproduction of2011mass failure|One forecast,30-minute cap; main jobs unchanged|
|17607147|Five-minute cluster receipt collector|Running; discovers refit arrays automatically|
|17603133|Initial calibration recovery|Two exact numerical repetitions verified|
|17597260|Extra initial fertility-age pilot|Numerical outputs verified by separate report recovery|

No complete fitted history, horizon certificate or completed rebated policy is
claimed. A0_6 first-window target1.974875 is matched by1.9731122320, a gap of
−0.001762768 within the0.005 fit requirement. Preference0.1289153142. The market,
PAYGO pension, property-tax rebate and exact replay checks pass. The next fitted
surprise is2011; households carry forward between accepted windows.

## Calibration and scientific contract

The recovered initial loss is179.2984252281 versus182.6491468669 for the rebated
seed. This is a verified candidate, not a converged optimizer. The original search
stopped after repeated20-evaluation numerical root limits. Complete13-row target
and17-row parameter/restriction tables: `initial_search_recovered/README.md`.
The displayed beta bound is the enforced annual0.99; raw scorer metadata retains
its old0.9995 separately. Housing moments remain weak.

The extra age pilot barely improves its augmented loss306.5873 to306.1929;
its age component worsens123.9381 to126.0059. It does not establish better age
fit and is not promoted. All original13, extra6 and augmented18-row tables are
in `age_pilot_recovered/README.md`. The extra synthetic5% scales are not empirical
standard errors. Recovery checked source, checkpoints, fiscal gaps and numerical
repetitions without another solve; serialized checkpoint hashes are provenance,
not numerical outputs that must reproduce byte for byte.

Twelve scored initial moments, nine free coordinates, unchanged weights, annual
beta cap0.99 and separate fertility normalization2.1 remain. Four unexpected
permanent preference changes are fitted at2007/2011/2015/2019; each vintage
expects current preferences to persist, with preferences fixed after2023.
Every property-tax comparison returns revenue equally per current household
head. PAYGO pensions balance separately at every accepted date. A0 removes all
post2023 migration; A+ retains the supplied migration sensitivity. Historical
head-age conditioning through2023 remains an imposed and disclosed bridge.

B0/B+ histories have not been launched: the child/person/household formation
conversion remains unresolved, and B+ also needs signed migrant-state allocation.
The demographic operator preflight17592167 passes both accounting identities;
its mass-preserving coefficient0.5647956 is diagnostic and not adopted.

The finite boundary evaluates remaining-lifetime household values at constant
boundary conditions. Boundary prices, pensions and rebates are solved on the
actual carried households. No stationary population reset is imported. Fiscal
feasibility beyond the boundary and insensitivity to horizon remain unverified.

## Numerical checks and budgets

Native exact-policy cache probe17598785 matched all dated values, policies,
household distributions, economic rows and accounting residuals exactly:
183.9926seconds uncached versus43.5013seconds cached, with11/12 repeated calls
reused. This is a constant-input mapping comparison, not a universal speedup.
All solver arguments enter the key, unsupported arguments bypass the cache,
and cache hits return fresh arrays. Cache bounds are6GiB for6/24dates and24GiB
for100dates. Main jobs request16/32/48GiB respectively and use one thread each.

The root cap is24 evaluations, with unchanged market/fiscal/exact-replay gates.
Numerical guesses shift one date between surprise vintages. Each history allows
at most24 primary forecast trials and one bounded alternative start; the shared
deadline and two-hour policy reserve still govern. The total worker ceiling is18.
Progress, latest completed cases and best-so-far receipts remain on the cluster.

Jacobian probes use the root solver's saved derivative matrix for the immediately
following preference trial in the same vintage. The opt-in default remains off;
failed roots and vintage changes reset reuse, and policies never inherit this
matrix. Ten focused root/resume tests pass. A0 native comparison passes:7 versus11 evaluations, maximum forecast fertility difference4.36e-7, exact replays zero. A+ hits the unchanged mass gate; cause remains under investigation. Sources:
`code/model/tools/run_e5f_final_rebated_history.py` and
`code/model/tools/run_e5f_forecast_jacobian_probe.py`. Their new snapshot is
`jacobian_source`; the running histories remain frozen in `history_source_cached`.

## Readout and traceability

`collect_e5f_final_history_readout.py` prepares the complete2023 observations
from accepted native2019/2023 snapshots without solving again. Compilation and
schema checks pass; native extraction awaits an accepted final historical window.
The stable17-graph packets are generated alongside successful forecasts. The
saved A+ first-trial contact sheet was visually inspected. Its legacy filename
`lifecycle_2023.csv` does not make that2007 snapshot a2023 result.

Source and job receipts: `jobs.json`, `resource_resubmission.json`,
`initial_search_recovered/`, `age_pilot_recovered/`; remote submissions remain
under `histories_refit`, `histories_refit_100_v2` and `jacobian_probes`.

Superseded jobs17592542/17592728/17593512/17593865 failed or were replaced before
any promotion. Initial17594990 passed the complete scoring loop in about six
minutes. Search17596347 stopped at numerical limits; recovery17603133 passed.
Handoff17605279 submitted refits; entirely pending17608564 was split into memory
classes, and entirely pending17613035 was replaced by17632922. All uncached
17595967 tasks were retired after the cached replacements passed. All old outputs
are retained. No result is promoted from failed harnesses17597051/17597204/17598460.

Verified restart support loads a pinned contiguous accepted history prefix and compares every inherited state field to the selected saved forecast. Native validation17641919 passed in1.89seconds without solving. Trace17642567 replays only the first2011forecast from that exact state with the existing stage profiler; it stops before any fallback, next trial or historical carry. The earlier mass trace17554347 had no reproduced failure. Float32 row sums are a hypothesis until the captured-cohort comparison establishes the source.
