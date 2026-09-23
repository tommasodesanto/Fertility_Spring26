# Bounded CPS and NCHS fertility check

This packet independently reproduces June 2004/2006 CPS fertility-stock
summaries and the 2003–2006 NCHS first-birth age convention. It reads local
inputs and the saved model-B diagnostic receipt only. It does not edit the
current builder, target values or weights, and it launches no model solve.

Run from the repository root with:

```bash
python3 output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_fertility/reproduce_cps_nchs.py
```

The script streams the local IPUMS CPS fixed-width gzip extract because the
uncompressed `.dat` path named in the September 10 receipt is absent. The gzip
extract reproduces both prior June-partition hashes and the exact valid sample
counts: 5,600 women in 2004 and 5,272 in 2006. Its sample is June interview,
female, ages 40–44, valid `FREVER` 0–20, and positive `FRSUPPWT`; the loader's
weight scaling is retained. The JSON records source hashes, the full June
partition hashes, invalid/missing fertility-code counts, weighted denominators,
and every requested moment. CSVs provide the CPS summary and annual NCHS
comparison.

For pooled CPS women, the weighted shares of 0, 1, 2 and 3+ children are
19.83%, 17.13%, 34.42% and 28.63%. The uncapped mean is 1.8784; the means
capped at three and five are 1.7184 and 1.8566. Among women with 3+ children,
the uncapped mean is 3.5589 and the cap-five mean is 3.4828. The cap-five
objects reproduce the earlier fertility-contract receipt; the cap-three
statistic is included because it was requested for this check. Replacing the
3+ group by the model's inherited 3.602359422009 while retaining this sample's
observed shares gives a mean of 1.8908. The separate identity diagnostic applies a hypothetical
2.1 mean to this **same CPS population**: with its observed childlessness and
one-child share, the implied mean among 2+ women is 3.0594, which corresponds
to a 66.11% 3+ share among 2+ women if the 3+ group is assigned 3.602359422009. This is
only within-sample arithmetic. It does not establish an inconsistency or a
loss lower bound between CPS ages 40–44 and the model's distinct age-cell
normalization.

The saved B receipt makes that distinction concrete. `solve_old_steady_state`
calls the active chain's `extract_moments`; the chain is loaded from
`run_e1_chain.py`, whose `extract_moments` uses the frozen solver's
`parity_dist`. That distribution sums ages from `P.A_f_end` onward. The active
four-year profile has age start 18, fertility through the period cell starting
at 42 (`A_f_end=7`), so normalization applies to the post-fertility ages 46+
population. Its saved shares are 14.83%, 17.92%, 31.30% and 35.94% for 0, 1,
2 and 3+. Applying 3.602359422009 to that same 46+ distribution gives
2.1000029; applying literal 3 gives 1.8835100, the saved legacy
`mean_completed_fertility`. This verifies the 2.1 identity on its own
population. The separate `target_age_fertility_moments` routine reads
`evaluation.g_current` at the cell indexed by age 42, but it is called by the
transition measurement path, not by `solve_old_steady_state`. The initial-
fertility observer also saves `g_post_fertility` age-specific shares; its
age-42 cell implies 2.0998737 with the same top-bin value. These distinct
measurement phases are retained separately.

The separate B [40,45) projection assigns shares 16.89%, 21.21%, 34.86% and
27.03% to 0, 1, 2 and 3+ under uniform birth-time interpolation. Its
constant-post-cell alternative gives 15.70%, 19.26%, 33.09% and 31.95%. These
are descriptive model outputs with the observer's household-proxy and age
projection caveats. No arithmetic comparison between that projected stock
and the post-fertility 46+ normalization identifies an incompatibility or a loss bound.
No arbitrary percentage-point cutoff is applied.

For NCHS pooled first births in 2003–2006, the raw recorded single-age mean is
25.1608 years. Adding 0.5 treats each age category as a single-year interval
and uses its midpoint; it is an assumption, not an exact continuous-age mean.
That convention gives 25.6608. The current four-year model-cell mapping gives
25.9763, a further 0.3155 years. The under-22 first-bin mapping contributes
0.2415 years in total, but that includes ages 18–21 inside model support. The
JSON and CSV separate the contribution and birth count/share for ages under 18,
18–21, 22–41, 42–45 and over 45. The 2006 raw recorded-age mean is 25.0476,
consistent with the cited NCHS report's rounded 25.0; adding the half-year
midpoint yields 25.5476. These are period count-weighted conventions, not
cohort fertility or individual-level maternal exposure rates.

The source table in the prior review cites [NCHS Data Brief No. 21](https://www.cdc.gov/nchs/products/databriefs/db21.htm)
for the published 2006 mean. The NCHS cache builder README documents the
birth-certificate age-code handling and links the [2003 NCHS detail-file
documentation](https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/Nat2003doc.pdf).
The CPS fixed-width fields and supplement weight scaling follow the local
IPUMS `loader.do`; the earlier CPS extraction receipt and current fertility
contract are retained as local provenance references. The active contract
specifies a cap of five for the CPS mean and conditional 3+ diagnostic; the
model top-bin representative 3.602359422009 is that cap-five conditional mean
from the June 2024 CPS vintage, not the 2004/2006 cap-five conditional mean.

Files:

- `reproduce_cps_nchs.py`: deterministic read-only extraction and calculations.
- `empirical_reproduction.json`: exact sample counts, source hashes, CPS
  moments, NCHS decomposition, and the saved B observer readout.
- `cps_weighted_moments.csv`: per-year and pooled weighted CPS table.
- `nchs_midpoint_comparison.csv`: annual NCHS single-age and cell-midpoint
  calculations.

No CPS standard errors are recomputed. The old/current empirical builders and
targets remain unchanged. Detailed data and source limitations are recorded in
the JSON receipt.
