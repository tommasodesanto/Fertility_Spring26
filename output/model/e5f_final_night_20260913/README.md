# Authorized final-night work

Latest evidence: September 13, 10:40 UTC. Fixed deadline: 18:00 UTC.
The cluster jobs and collector continue independently of the laptop. Local idle
sleep prevention lasts until the deadline. Latest account check: 58% weekly
remaining; retain the author's 20% floor.

## Verified progress

The corrected initial equilibrium passed two exact numerical reproductions in
job 17655042. Full target and parameter tables: `corrected_initial/README.md`.
The probability repair is numerical: it does not change the structural search
coordinates, empirical targets, weights or acceptance tolerances. It changes the
loss by less than 0.000001 and does not improve the economic fit.

Both six-period cases passed their first historical window and carried households
to 2011. With ordinary starts, the A0 fertility gap is −0.001750752 and the A+
gap is +0.003827671 against target 1.974875, within the 0.005 requirement.
All market, PAYGO, rebate, household and exact-replay gates pass. Price-seeded
counterparts also pass in two mappings, about 174 seconds; ordinary starts took
14 mappings, about 1,223 seconds. Different nodes prevent a controlled timing claim.

## Live jobs

| Job | Work | State at update |
|---|---|---|
| 17658836 | A0/A+, 6 periods, ordinary starts | Fitting subsequent surprises |
| 17661737 | A0/A+, 6 periods, pinned numerical price starts | Fitting subsequent surprises |
| 17663940 | A0/A+, 24 periods, ordinary starts | Submitted by verified handoff |
| 17663986 | A0/A+, 100 periods | Submitted by verified handoff |
| 17664449 | A0/A+, 24 periods, pinned numerical price starts | Submitted |
| 17607147 | Five-minute receipt collector | Running |

Handoff 17658318 completed after verifying both six-period first-window gates.
The superseded old-source arrays 17613033, 17613034 and 17632922 were cancelled;
all their outputs remain. No seeded 100-period duplicate was submitted.

No complete four-shock history, horizon certificate or completed policy is claimed.
The fiscal and mass problems have verified numerical repairs; the economic fit
and full historical/policy exercise still require assessment.

Remote batch: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/`.
Plan: `docs/model/e5f_two_closure_overnight_plan.md`.

## Calibration and scientific contract

The corrected initial loss is 179.2984242480 versus 182.6491468669 for the rebated
seed. This is a verified candidate, not a converged optimizer. The original search
stopped after repeated20-evaluation numerical root limits. Complete13-row target
and17-row parameter/restriction tables: `corrected_initial/README.md`.
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
The target contract maps decision vintages2007/2011/2015/2019 to published TFR
means for2008–2011/2012–2015/2016–2019/2020–2023 respectively. Fit plots must show
these four-year observation windows explicitly; a2023 decision-vintage fertility
flow belongs to the forecast period, rather than another fitted historical row.
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

The native corrected first forecasts now verify the combined cache, row
normalization, optional Jacobian reuse and fiscal-polish controller. Numerical
price seeding reuses only pinned root coordinates; it never reuses old-source
household distributions or Jacobians. Fourteen focused seed/driver tests pass.

## Readout and traceability

`verified_initial_readout.pdf` contains the three-page assessment, complete
initial tables and unchanged 17-figure appendix. It is explicitly an initial
readout while history/policy work runs. All 12 rendered pages were visually
reviewed. Regenerate without a model solve using the bundled Python runtime:

```sh
python code/model/tools/build_e5f_final_night_report.py --packet output/model/e5f_final_night_20260913 --output output/model/e5f_final_night_20260913/verified_initial_readout.pdf --as-of "13 September 2026, 10:45 UTC"
```


`collect_e5f_final_history_readout.py` prepares the complete2023 observations
from accepted native2019/2023 snapshots without solving again. Compilation and
schema checks pass; native extraction awaits an accepted final historical window.
The stable17-graph packets are generated alongside successful forecasts. The
saved A+ first-trial contact sheet was visually inspected. Its legacy filename
`lifecycle_2023.csv` does not make that2007 snapshot a2023 result.

Source and job receipts are in `jobs.json`, the corrected-source receipts,
`corrected_initial/`, and the native remote batch. Historical failed and superseded
outputs remain preserved. No result is promoted from a failed numerical harness.
