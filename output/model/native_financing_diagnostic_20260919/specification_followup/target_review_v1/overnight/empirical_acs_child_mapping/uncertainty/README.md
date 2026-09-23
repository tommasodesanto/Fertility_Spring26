# Paired metro-bootstrap uncertainty for family rooms

This receipt quantifies uncertainty for the existing ACS `family_rooms` contrast and the parent-linked under-18 candidate from the adjacent child-mapping diagnostic. It keeps the same ACS 2005/2006 head sample, 42 MET2013 metros, `HHWT` weighting, room cap at 9, and fixed group memberships within each metro. It adds no geography or target choice, and the candidate remains diagnostic only.

## Authoritative bootstrap method

The source is `output/model/e5f_matched_pf_20260909a/design_research/housing/summarize_early_housing.py` and its saved `early_housing_candidate_receipt.json`. The builder's exact method is:

- 1,000 replications, NumPy `default_rng(20260910)`.
- For each replication, draw 42 metro indices with replacement from the 42 active metro IDs and convert the draw to metro-frequency counts.
- Apply the same frequency vector to all within-metro household and group totals. For the two dates pooled into the 2005–06 target, sum metro sufficient statistics across both years first.
- Recompute each group mean as the frequency-weighted capped-room total divided by the frequency-weighted `HHWT` total, then subtract the 1–2-child mean from the 3+-child mean.
- Report the sample standard deviation/covariance across replications (`ddof=1`). This is empirical metro-cluster resampling uncertainty, not an ACS official design or replicate-weight SE.

Group membership is held fixed within each metro for all replications. The bootstrap varies the metro multiplicities; it does not resample households or re-estimate the age/child grouping. The candidate retains the same head sample and uses children linked to the head by `MOMLOC` or `POPLOC`, counting those under age 18.

The prior child-mapping output had only overall candidate totals, so one sequential memory-mapped pass over ACS 2005–06 was needed to form candidate sufficient statistics by metro. It wrote no person- or household-level records. The all-age parent-linked count matches `NCHILD` bins capped at 3 in all 275,009 households; because `NCHILD` is top-coded, this does not validate exact counts above 3.

## Results

All contrasts are in rooms after the `min(ROOMS,9)` cap.

| Contrast | Point | Metro-bootstrap SE | Inverse-variance weight |
|---|---:|---:|---:|
| Existing `NCHILD` 3+ minus 1–2 | 0.347067 | 0.059705 | 280.528084 |
| Parent-linked own children under 18, 3+ minus 1–2 | 0.336220 | 0.058386 | 293.350607 |
| Candidate minus existing | −0.010847 | 0.008077 | — |

The paired bootstrap covariance matrix, ordered as existing then candidate, is

\[
\begin{pmatrix}
0.0035647055 & 0.0034541806\\
0.0034541806 & 0.0034088902
\end{pmatrix},
\]

with correlation 0.990893. The percentile 95% interval for candidate minus existing is [−0.030346, 0.001362] rooms. This resampling result describes sensitivity across the fixed 42-metro sample; it is not a test of whether either measurement is adequate or an ACS design-based confidence interval.

The old point reproduces exactly, the old SE differs from the saved builder SE by less than \(10^{-12}\), and the inverse-variance weight differs from the saved housing receipt by less than \(10^{-9}\). The candidate point reproduces the prior aggregate-only value within \(10^{-9}\).

## Artifacts

- `metro_bootstrap_family_rooms.py` — exact-method bootstrap plus a single-pass candidate metro-stat builder; asserts source and prior aggregate reproduction.
- `metro_sufficient_stats.csv` — per-metro old/new group counts, `HHWT` totals, capped-room weighted totals, family counts, and linked-bin checks. No microdata rows.
- `paired_bootstrap_draws.csv` — 1,000 paired old, candidate, and difference estimates.
- `bootstrap_uncertainty.json` — sample, method, results, covariance, hashes, and saved-value reproduction.

Input SHA-256 values and source identity are in the JSON. The raw extract hash is reused from the canonical housing receipt rather than recomputed; file size and modification time match that receipt.
