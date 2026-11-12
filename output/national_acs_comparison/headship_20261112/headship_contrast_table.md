| Design | Outcome | Raw level at -1 | +3 minus -1 (SE) | Regression N | Persons | Households | Anchors |
|---|---|---:|---:|---:|---:|---:|---:|
| A baseline (all positions) | rooms9 | 5.41 | 0.488 (0.008) | 14,155,689 | 3,283,266 | 3,155,246 | 130,601 |
| A baseline (all positions) | bedrooms5 | 2.60 | 0.337 (0.004) | 14,155,711 | 3,283,266 | 3,155,246 | 130,601 |
| A baseline (all positions) | ownership_lw | 57.78 | 3.291 (0.215) | 14,155,711 | 3,283,266 | 3,155,246 | 130,601 |
| B householder/spouse/partner | rooms9 | 5.21 | 0.593 (0.008) | 11,254,007 | 2,477,463 | 2,459,662 | 121,309 |
| B householder/spouse/partner | bedrooms5 | 2.46 | 0.407 (0.005) | 11,254,029 | 2,477,463 | 2,459,662 | 121,309 |
| B householder/spouse/partner | ownership_lw | 54.59 | 4.605 (0.225) | 11,254,029 | 2,477,463 | 2,459,662 | 121,309 |
| C householder/spouse | rooms9 | 5.25 | 0.588 (0.009) | 10,797,894 | 2,305,783 | 2,295,167 | 116,278 |
| C householder/spouse | bedrooms5 | 2.48 | 0.404 (0.005) | 10,797,916 | 2,305,783 | 2,295,167 | 116,278 |
| C householder/spouse | ownership_lw | 55.81 | 4.483 (0.228) | 10,797,916 | 2,305,783 | 2,295,167 | 116,278 |
| D households only (no group quarters) | rooms9 | 5.41 | 0.489 (0.008) | 14,155,689 | 3,221,864 | 3,093,844 | 130,601 |
| D households only (no group quarters) | bedrooms5 | 2.60 | 0.337 (0.004) | 14,155,711 | 3,221,864 | 3,093,844 | 130,601 |
| D households only (no group quarters) | ownership_lw | 57.68 | 3.378 (0.215) | 14,155,711 | 3,221,864 | 3,093,844 | 130,601 |
| E donor-only household-position restriction | rooms9 | 5.15 | 0.697 (0.008) | 11,815,841 | 2,611,969 | 2,589,196 | 130,601 |
| E donor-only household-position restriction | bedrooms5 | 2.43 | 0.475 (0.004) | 11,815,863 | 2,611,969 | 2,589,196 | 130,601 |
| E donor-only household-position restriction | ownership_lw | 52.61 | 7.227 (0.219) | 11,815,863 | 2,611,969 | 2,589,196 | 130,601 |
| A, first birth at 30-45 | rooms9 | 5.50 | 0.691 (0.011) | 8,287,809 | 1,758,340 | 1,728,579 | 73,832 |
| A, first birth at 30-45 | bedrooms5 | 2.63 | 0.436 (0.006) | 8,287,823 | 1,758,340 | 1,728,579 | 73,832 |
| A, first birth at 30-45 | ownership_lw | 62.21 | 7.764 (0.286) | 8,287,823 | 1,758,340 | 1,728,579 | 73,832 |
| B, first birth at 30-45 | rooms9 | 5.38 | 0.789 (0.012) | 7,497,375 | 1,503,967 | 1,496,538 | 70,206 |
| B, first birth at 30-45 | bedrooms5 | 2.54 | 0.500 (0.006) | 7,497,389 | 1,503,967 | 1,496,538 | 70,206 |
| B, first birth at 30-45 | ownership_lw | 60.49 | 9.012 (0.297) | 7,497,389 | 1,503,967 | 1,496,538 | 70,206 |
| B, anchors with donors at every pre-birth lag | rooms9 | 5.25 | 0.546 (0.009) | 9,331,335 | 2,154,853 | 2,146,891 | 57,061 |
| B, anchors with donors at every pre-birth lag | bedrooms5 | 2.47 | 0.395 (0.005) | 9,331,357 | 2,154,853 | 2,146,891 | 57,061 |
| B, anchors with donors at every pre-birth lag | ownership_lw | 54.21 | 3.246 (0.240) | 9,331,357 | 2,154,853 | 2,146,891 | 57,061 |

PSID references (not like for like; sample and windows stated): women-only household sample (women reference person or spouse/partner each interview) rooms 1.025 (0.073), +3/+4 vs -3/-2; main sample (men and women, reference person or spouse 2-3 years before the birth) rooms 1.465 (0.050), +3/+4 vs -3/-2. Source: code/data/psid_followup_mar2026/output/sa_rooms_first_birth_v2/summary.csv (H, A2h).

A, anchors with donors at every pre-birth lag (A_fullpre; Torch jobs 22861217 and rerun 22869863 both aborted with R heap corruption on entering the bedrooms fits, after rooms and ownership completed; values from the run log, full spec): rooms 0.4438 (0.0086); owner-occupied dwelling 2.09 pp (0.23); regression N 11,918,138; 2,588,298 source households. Bedrooms not available.
