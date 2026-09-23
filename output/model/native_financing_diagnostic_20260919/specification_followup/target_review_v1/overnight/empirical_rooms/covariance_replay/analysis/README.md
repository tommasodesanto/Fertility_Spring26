# Post-replay covariance analysis

`analyze_replay.py` reads only the single reviewed Stata run under the sibling
`output/` folder and the prior common-cohort point-estimate review. It refuses
analytical calculations until `execution.json` is terminal with exit code zero,
the executed `.do` hash matches the lead-approved hash, the target and SE are
within the approved (10^{-9}) reproduction tolerance, and observation/person
counts are exactly 49,457/4,112. It also checks the event-study baseline and full
interaction covariance dimensions, symmetry, and agreement between covariance
diagonals and exported interaction variances.

After those gates pass, the script computes pre-weighted and post-weighted
common-cohort (-1\to+3) contrasts using the exact cohort shares from the same
regression's `e(sample)` support. It also computes separate pre-weighted and
post-weighted common-cohort (-2\to+2) candidates. For the latter, event time
(-2) is the omitted `F2event` reference, set to zero only for cohorts with
positive `e(sample)` support at both (-2) and (+2).

Each variance is (a'Va), with (V) the full interaction block from `e(V)` and
(a) the fixed-weight cohort/event contrast vector. An independent cohort
contrast covariance calculation checks the same linear form. These standard
errors are conditional on the selected cohort weights; they omit sampling
uncertainty in estimated cohort shares. The outputs are candidate aggregations,
not an adopted target, calibration change, or causal claim.

When the replay is complete, run:

```bash
python3 output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/covariance_replay/analysis/analyze_replay.py
```

The script writes `replay_analysis.json` and `cohort_weights.csv` beside itself.
It has no regression or person-level data access. If invoked while the replay is
still running, it prints the pending status and exits without producing an
analysis.

## Completed replay receipt

The single reviewed replay finished in 330.1 seconds at four processors with
Stata exit code 0. It reproduced the frozen target as 0.720246262381528
(absolute gap (2.22\times10^{-16})) and its SE as 0.0852600513385958 (gap 0).
The target and SE checks use the approved (10^{-9}) tolerance. The original
sample counts matched exactly: 49,457 woman-household-years and 4,112 women.
The replay preserved the event-study curve, including the explicit zero at the
omitted (-2) period.

| Candidate | Estimate | Conditional SE | Cohorts | Endpoint cohorts | Retained weight at start/end |
|---|---:|---:|---:|---:|---:|
| Common (-1\to+3), weights fixed at (-1) | 0.81139749 | 0.07660629 | 27 | 34 / 32 | 80.89% / 88.31% |
| Common (-1\to+3), weights fixed at (+3) | 0.79775926 | 0.07545040 | 27 | 34 / 32 | 80.89% / 88.31% |
| Common (-2\to+2), weights fixed at (-2) | 0.73459943 | 0.07562861 | 27 | 34 / 33 | 83.60% / 87.51% |
| Common (-2\to+2), weights fixed at (+2) | 0.70389142 | 0.07542286 | 27 | 34 / 33 | 83.60% / 87.51% |

The two (-1\to+3) point estimates reproduce the earlier saved-coefficient
candidates within (6.7\times10^{-16}). The full interaction covariance has
900 coefficients and 810,000 cells; its maximum symmetry discrepancy and its
maximum diagonal-versus-marginal-variance discrepancy are both zero. The
independent cohort-contrast and (a'Va) calculations also agree to below
(9\times10^{-19}) for the reported candidates. Each SE is conditional on its
fixed cohort weights and excludes sampling variation in estimated cohort
shares. These aggregations remain candidates, without adoption or causal claims.

The original fit left room value 0 unchanged. There are 100 such rows in
`e(sample)`, and the minimum nonmissing `rooms` value is 0. The `e(sample)` raw
value counts are 2,266 for code 9, 1 for code 98, and 0 for code 99. The source
uses year-specific rules for codes 9, 98, and 99, so their raw counts do not
classify every observation as missing. The room-code receipt preserves the
requested counts without filtering or refitting. This flags the 100 zero-room
observations for author review; it does not revise the frozen estimate.
