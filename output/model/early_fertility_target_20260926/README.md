# Early fertility target from the historical CPS supplement

This packet estimates cumulative live births among women recorded at exact age
25 in pooled June 2004 and June 2006 CPS fertility supplements. It does not
change calibration targets or weights. Ages 24, 26, and pooled 24--26 are
diagnostics; they do not replace the exact-age object.

Verified result (Torch job 18619369, completed in 11 seconds):

| Age | Records | Raw mean births | Mean capped at three | Capped-three bootstrap SE |
|---|---:|---:|---:|---:|
| 24 | 1,662 | 0.684 | 0.670 | 0.025 |
| **25** | **1,774** | **0.857** | **0.810** | **0.028** |
| 26 | 1,691 | 0.959 | 0.923 | 0.030 |
| 24--26, diagnostic | 5,127 | 0.832 | 0.799 | 0.016 |

The exact-age-25 model-matched value is `0.8095276384290021`; the raw value is
`0.8569601514801901`. Both original June-partition hashes and the schema hash
match. The control replays 10,872 women ages 40--44, childlessness
`0.19827875100684098`, and exactly one child conditional on motherhood
`0.2136553252201454`. A separate AWK implementation reproduces age-25 count,
weight, raw mean, and capped-three mean to floating-point precision; its
receipt is `independent_awk_check.json`.

The model-matched candidate is the weighted mean of `min(FREVER, 3)`, because
the model's highest child-count state is 3+. The uncapped mean and a cap-at-five
sensitivity are reported separately. Age 25 describes the completed age at the
interview, hence interval [25, 26); it does not measure births exactly on the
25th birthday. This is a stock of births per woman, not an annual birth rate.

The builder is `code/data/cps_fertility/build_early_fertility_target.py`. It
reuses the schema and exact June-partition hashes in
`output/model/e5f_matched_pf_20260909a/parameter_target_audit/fertility/fertility_availability.json`.
Women have `SEX=2`, valid `FREVER` 0--20, and positive `FRSUPPWT`; the latter is
divided by 10,000 per the source loader. Pooling sums weighted numerators and
denominators across records, rather than averaging annual means. A replay of
the historical ages 40--44 childlessness and exactly-one-among-mothers targets
must pass before any estimate is reported.

`early_fertility_target.json` contains point estimates, annual diagnostics,
sample counts, partition validation, control replay, and uncertainty.
`early_fertility_target.csv` is its compact age comparison.

Uncertainty is a 2,000-draw person bootstrap stratified by survey year with
fixed annual sample counts, retaining original weights (seed 20260926, plus
age for individual-age estimates). It is **not CPS design-consistent**: no
PSU/stratum or replicate-weight correction or cross-year covariance is applied.
These standard errors are diagnostic and do not automatically establish an
SMM weight. Annual samples and pooled-age sensitivity are shown separately.

All parsing, decompression, hashing, estimation, and resampling run on Torch.
The existing compressed source was transferred without local decompression.
`run_torch.sh` runs the builder on one CPU with a ten-minute limit at
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/early_fertility_target_20260926`.
The first pass writes the two authenticated June partitions. Subsequent runs
can omit `--compressed-source` to read just those retained partitions.
