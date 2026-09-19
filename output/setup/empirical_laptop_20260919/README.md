# Empirical laptop replay — September 19, 2026

Both selected real-data pipelines passed on the Apple M5 Pro. This verifies
these workflows, not every empirical script in the repository. No production
source, estimator, sample, weight, event-time definition or saved output changed.

| Pipeline | Elapsed | Peak RSS | Reproduction |
| --- | ---: | ---: | --- |
| R PSID tenure/liquid-wealth gradient | 8 s (whole-second wrapper) | 1.17 GiB | Both CSVs match within 2e-15 |
| Stata first-birth housing event study | 363.45 s | 2.90 GiB | All 19 event-study rows match exactly; target receipt numeric gap below 7e-15 |

The source data are the existing 5.9 GiB
`/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta`.
SHA-256: `f1c5d48d5ef5357c40e743895dd25c73cd45789f779d1d575408f671fe637029`.
The input was hashed once and the fingerprint reused for the two pipelines.

## R

Source: `code/data/psid_followup_mar2026/build_tenure_liquid_wealth_gradient.R`.
Native R 4.6.1; original 500 bootstrap draws and seed retained. Primary and
sensitivity samples contain 2,933 and 2,931 individuals. The copied script
changes only repository/output path assignments; thread limits were set in
the execution environment. The full PSID input was used.

The two output CSVs are in `r/run_output/`. `r/comparison_report.json` records
equal dimensions, headers and missingness and no numeric differences above
1e-9. Independent lead comparison found maximum absolute gap 2e-15. Runtime
and peak RSS come from `r/run.status` and `r/run.stderr`. An initial wrapper
attempt could not invoke unavailable `timeout`; no R pipeline ran in that
attempt. The direct run then completed in eight seconds.

## Stata

Source: `code/data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do`.
StataMP 17, two cores instead of the source's eight. The only other change to
the isolated script is the output-root path. `stata/relocation.diff` records
both substitutions. The full sample and estimator were retained: 49,457
estimation observations, 4,112 individuals, weighted Sun–Abraham event study,
person/year fixed effects, age/education controls, individual-ID clustering.

The replay's full output tables are under
`stata/sa_rooms_first_birth_household_aligned_v1/`. `stata/comparison.json`
checks all 120 numeric cells across the 19-row event-study table and one-row
target receipt, plus all text fields. Only runtime is excluded. The largest
numeric gap is 6.88e-15, below the 1e-9 absolute tolerance. The baseline table
hashes remain unchanged. `stata/execution.json` records successful process
exit, elapsed time and peak child RSS; successful completion is also confirmed
by the final target receipt and closed Stata log. The run was capped at 12 minutes.

The saved older Stata receipt reports 511.832 seconds with eight cores. This
is contextual timing evidence, not a controlled same-core hardware benchmark.

## Reproduction and safety

The relocated source copies and source hashes are preserved with the receipts.
Use a new output directory for any rerun. The Stata script intentionally refuses
to overwrite a completed target receipt. Do not run either original builder
merely to repeat this check, because its normal destination is production output.
No local empirical process remained after validation.
