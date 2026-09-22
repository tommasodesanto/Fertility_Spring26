# ACS Twin1 / SameSex2 national housing IV — final results (2026-09-22)

**Job 18274624, COMPLETED 0:0, elapsed 00:42:50, batch MaxRSS 77809928K (~74.21GiB). 0 warnings in this national run's log** (the 168 fixest VCOV-not-PSD warnings seen in the prior VT recovery smoke, job 18271160, were a separate, smaller-sample observation and are not claimed resolved, reproduced, or harmless here).

Code: `twins_samesex_iv_recovery_0f239725` (checksummed, read-only). All 18 design×outcome fits reached `full_fit` and pass internal consistency checks (nobs rf==fs==iv, AR error-free, full V symmetric and dimension-matched to named coefficients, V diagonal matches reported SE, primary-receipt coefficients equal final-receipt coefficients). Full per-fit check log: [`national_128g_results/review_summary.json`](national_128g_results/review_summary.json).

## Sample construction (reproduces prior construction, not new biology/social-link proof)

51 state+DC source partitions (2005–2019 ACS 1-year product, `SAMPLE==YEAR*100+1`): 59,046,776 raw rows read, 12,672,840 excluded as out-of-year, 0 excluded as non-1yr product, 46,373,936 kept. Mothers age 21–35 with oldest linked child <18: 4,024,238 unique mothers, 3,783,862 households. `MOMLOC` links are coresident/social, not confirmed biological histories. These are **contemporaneous cross-sections**, not longitudinal or pre-birth histories.

- **Twin1**: age-only twin-like proxy (extract27 has no `BIRTHQTR` — not a confirmed-twin test), event age = oldest linked child's age (pooled 0:5), treatment = ≥2 linked children. 857,607 eligible, 16,496 instrument-positive.
- **SameSex2**: oldest-two linked children same sex, event age = second-oldest child's age (pooled 0:5), treatment = ≥3 linked children — a **third-child margin, a different population from Twin1's second-child margin**. 656,986 primary-eligible (of 1,132,295 pool; excludes 32,186 primary-age ties), 329,974 same-sex positive (173,813 both-boys, 156,161 both-girls).

## Methods

Mother rows are unique by `(YEAR, SAMPLE, SERIAL, PERNUM)`, weighted by mother `PERWT`. Standard errors are clustered by household, `(YEAR, SAMPLE, SERIAL)`. Control formula: `i(mat_age_at_event) + i(RACE) + i(survey_year) + i(event_age)` — indicator (dummy) sets for mother's inferred age at the event (`mat_age_at_event = mother's ACS-interview age − event age`), race, ACS survey year, and event age itself (0:5). Confidence intervals are the normal-quantile approximation (±1.96×clustered SE) from the reported clustered covariance. The reported F is a single-instrument cluster-robust Wald statistic, `F = (FS coefficient / clustered SE)^2` — this is **not** a Kleibergen–Paap statistic. Outcomes: `ROOMS` valid codes 1:27,30, capped at 9 for the primary outcome; `OWNERSHP` 1=owner/2=renter, recoded to 0/1; `BEDROOMS` (diagnostic only) valid codes 1:22, recoded to (code−1) capped at 5. Sample effective size (ESS) was **not saved** by the collected receipts and is not reported below — not fabricated, not recomputed from a heavy read.

## Table 1 — Reduced form + first stage (PRIMARY objects, pooled event age 0:5)

RF = effect of the instrument (Twin1 age-only proxy / SameSex2 oldest-two-same-sex) on the housing outcome; interpret as **instrument–outcome association**, not an established causal additional-child effect. FS = effect of the instrument on the additional-child treatment probability, in **percentage points**.

| Design | Outcome | N | HH | Z-pos | RF coef | RF SE | RF 95% CI | Unit | FS coef (pp) | FS SE (pp) | First-stage F |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Twin1 | Rooms | 857607 | 855366 | 16496 | 0.3286 | 0.0172 | [0.2950, 0.3622] | rooms | 66.2141 | 0.2380 | 77396.4 |
| Twin1 | Ownership | 857607 | 855366 | 16496 | 3.1060 | 0.4608 | [2.2029, 4.0092] | pp owner | 66.2141 | 0.2380 | 77396.4 |
| SameSex2 | Rooms | 656985 | 656159 | 329974 | -0.0376 | 0.0051 | [-0.0476, -0.0276] | rooms | 3.2745 | 0.1224 | 715.2 |
| SameSex2 | Ownership | 656986 | 656160 | 329974 | -0.2457 | 0.1475 | [-0.5348, 0.0434] | pp owner | 3.2745 | 0.1224 | 715.2 |

## Table 2 — 2SLS and Anderson–Rubin (ASSUMPTION-DEPENDENT diagnostics, not causally certified)

2SLS is the assumed instrumented **effect of the additional-child treatment** on the outcome. Reported only because the first stage is non-degenerate; the exclusion restriction is **unresolved** for both designs (see caveats below) so this is a diagnostic, not a causal estimate.

| Design | Outcome | 2SLS coef | 2SLS SE | 2SLS 95% CI | Unit | AR accepted set (tested grid) | AR honesty flag |
|---|---|---|---|---|---|---|---|
| Twin1 | Rooms | 0.4962 | 0.0259 | [0.4454, 0.5470] | rooms | [0.5000, 0.5000] | SINGLE TESTED GRID POINT (0.5000) accepted -- not a zero-width CI, interior points untested between grid steps |
| Twin1 | Ownership | 4.6909 | 0.6955 | [3.3278, 6.0540] | pp owner | [4.0000, 6.0000] | interior-bounded, 1 component(s), grid-truncated approximation |
| SameSex2 | Rooms | -1.1491 | 0.1638 | [-1.4702, -0.8280] | rooms | [-1.4000, -0.9000] | interior-bounded, 1 component(s), grid-truncated approximation |
| SameSex2 | Ownership | -7.5029 | 4.5130 | [-16.3481, 1.3423] | pp owner | [-16.0000, 0.0000] | interior-bounded, 1 component(s), grid-truncated approximation |

## Descriptive interpretation

Twin1 and SameSex2 show **opposite-signed** room-count RF (Twin1 positive, SameSex2 negative). This is a descriptive contrast across **different samples and margins** (second-child vs. third-child; different instrument construction and different populations), not evidence that an additional child causally shrinks the home — the SameSex2 result is also consistent with a direct room-sharing/allocation response to same-sex composition unrelated to any third-child effect (the exclusion restriction concern above). Neither RF is a like-for-like counterfactual of the other.

## Identification caveats

- **RF and FS are the primary, most defensible objects.** 2SLS is an assumption-dependent diagnostic, not a causally certified estimate.
- **Exclusion restriction is unresolved for both designs.** Same-sex composition of the oldest two children may directly change room-sharing/housing demand independent of any third-child effect. Twin-like status associates with maternal health conditions and birth spacing that independently affect housing. Neither is addressed by this run.
- **Twin1 is an age-only proxy**, not a confirmed-twin indicator — extract27 has no `BIRTHQTR` to bridge to child birth quarter.
- **Twin1 and SameSex2 identify different populations/margins** (second-child vs. third-child) and must not be pooled or compared as a single estimate.
- **AR intervals are finite-grid approximations.** Some cells' accepted set is a single tested grid point (see the `ar_honesty_note` column in the enhanced CSV) — this is not a zero-width continuous confidence interval; untested points between grid steps are not certified rejected or accepted.
- No causal certification and no calibration adoption is implied by this run.

## Full results table (all 18: 2 designs × 3 outcomes × 3 windows, no cherry-picking)

[`national_128g_results/national_18row_table.csv`](national_128g_results/national_18row_table.csv) -- original file, **native probability units**, unmodified. [`national_128g_results/national_18row_table_with_pp.csv`](national_128g_results/national_18row_table_with_pp.csv) -- same data plus added percentage-point and AR-honesty columns (original columns untouched).

Brief event-age-3/5 read (from the full CSV, no selection by significance): Twin1 rooms RF is positive at both event age 3 and 5, consistent in sign with the pooled estimate; SameSex2 rooms RF is negative at both, also consistent with the pooled estimate. Magnitudes and precision vary across the smaller event-age-specific subsamples; see the CSV for exact values rather than a restated table here.

Figure: [`national_128g_results/national_main_results.png`](national_128g_results/national_main_results.png) / [`.pdf`](national_128g_results/national_main_results.pdf)

## Prior history (preserved, not overwritten)

- VT recovery-smoke review (job 18271160, computational gate only): [`recovery_smoke_review/`](recovery_smoke_review/)
- Original national OOM failure and diagnosis (job 18247804, superseded by this run): preserved in git history.
