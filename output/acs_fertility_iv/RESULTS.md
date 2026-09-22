# ACS Twin1 / SameSex2 national housing IV — final results (2026-09-22)
**Job 18274624, COMPLETED 0:0, 42m50s, batch MaxRSS 77,809,928K (~74.2GiB) at 128GB. Driver-reported timer: 2566.8s. 0 warnings in this national run's log** (the 168 fixest VCOV-not-PSD warnings seen in the prior VT recovery smoke, job 18271160, were a separate, smaller-sample observation and are not claimed resolved, reproduced, or harmless here).
Code: `code_snapshots/twins_samesex_iv_recovery_0f239725` (checksummed, read-only). All 18 design×outcome fits reached `full_fit`; all pass internal consistency checks (nobs rf==fs==iv, AR error-free, full V symmetric and dimension-matched to named coefficients, V diagonal matches reported SE, primary-receipt coefficients equal final-receipt coefficients). See `review_summary.json` for the complete per-fit check log.
## Sample construction (reproduces prior construction, not new biology/social-link proof)
51 state+DC source partitions (2005–2019 ACS 1-year product, `SAMPLE==YEAR*100+1`): 59,046,776 raw rows read, 12,672,840 excluded as out-of-year, 0 excluded as non-1yr product, 46,373,936 kept. Mothers age 21–35 with oldest linked child <18: 4,024,238 unique mothers, 3,783,862 households. `MOMLOC` links are coresident/social, not confirmed biological histories.
- **Twin1**: age-only twin-like proxy (extract27 has no `BIRTHQTR` — not a confirmed-twin test), event age = oldest linked child's age (pooled 0:5), treatment = ≥2 linked children. 857,607 eligible, 16,496 instrument-positive.
- **SameSex2**: oldest-two linked children same sex, event age = second-oldest child's age (pooled 0:5), treatment = ≥3 linked children — a **third-child margin, a different population from Twin1's second-child margin**. 656,986 primary-eligible (of 1,132,295 pool; excludes 32,186 primary-age ties), 329,974 same-sex positive (173,813 both-boys, 156,161 both-girls).
## Main results — pooled event age 0:5 (RF and FS primary; 2SLS assumption-dependent)
| Design | Outcome | N | HH | Z-positive | FS coef (extra-child pp) | FS clustered SE | First-stage F | RF coef | RF 95% CI | 2SLS coef | 2SLS 95% CI | AR (finite-grid) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Twin1 | Rooms (rooms (capped 9)) | 857607 | 855366 | 16496 | 0.6621 | 0.0024 | 77396.4 | 0.3286 | [0.2950, 0.3622] | 0.4962 | [0.4454, 0.5470] | [0.5, 0.5] bounded, 1 component(s), 0 grid errors |
| Twin1 | Ownership (pp owner) | 857607 | 855366 | 16496 | 0.6621 | 0.0024 | 77396.4 | 0.0311 | [0.0220, 0.0401] | 0.0469 | [0.0333, 0.0605] | [0.04, 0.06] bounded, 1 component(s), 0 grid errors |
| SameSex2 | Rooms (rooms (capped 9)) | 656985 | 656159 | 329974 | 0.0327 | 0.0012 | 715.2 | -0.0376 | [-0.0476, -0.0276] | -1.1491 | [-1.4702, -0.8280] | [-1.4, -0.9] bounded, 1 component(s), 0 grid errors |
| SameSex2 | Ownership (pp owner) | 656986 | 656160 | 329974 | 0.0327 | 0.0012 | 715.2 | -0.0025 | [-0.0053, 0.0004] | -0.0750 | [-0.1635, 0.0134] | [-0.16, 0] bounded, 1 component(s), 0 grid errors |

Units: ROOMS outcome is raw room count capped at 9; OWNERSHP outcome is owner=1/renter=0 (coefficients here are in probability units, i.e. percentage-point effects ×100 for pp). FS coefficient is the effect of the instrument on the probability of the additional-child treatment. AR intervals are the finite-grid Anderson–Rubin accepted set (Twin1: grid [-3,3] by 0.1 for rooms, [-0.5,0.5] by 0.02 for ownership; SameSex2 same grids) — a grid-truncated approximation, not an exact analytic boundary.
## Interpretation and identification caveats
- **RF and FS are the primary, most defensible objects.** 2SLS is reported as an assumption-dependent diagnostic, not a causally certified estimate.
- **Exclusion restriction is unresolved for both designs.** Same-sex composition of the oldest two children may directly change room-sharing/housing demand independent of any third-child effect. Twin-like status associates with maternal health conditions and birth spacing that independently affect housing. Neither is addressed by this run.
- **Twin1 is an age-only proxy**, not a confirmed-twin indicator — extract27 has no `BIRTHQTR` to bridge to child birth quarter.
- **Twin1 and SameSex2 identify different populations/margins** (second-child vs. third-child) and must not be pooled or compared as a single estimate.
- No causal certification and no calibration adoption is implied by this run.
## Full 18-row table (all outcomes × pooled + event-age 3/5, incl. Bedrooms diagnostic)
See `national_18row_table.csv` for the complete, non-cherry-picked table (18 rows: 2 designs × 3 outcomes [Rooms, Ownership, Bedrooms] × 3 windows [pooled 0:5, event age 3, event age 5]).
## Prior history (preserved, not overwritten)
- VT recovery-smoke review (job 18271160, computational gate only): `output/acs_fertility_iv/recovery_smoke_review/`
- Original national OOM failure and diagnosis (job 18247804, superseded by this run): prior content of this file is preserved in git history at commit prior to this rewrite.
