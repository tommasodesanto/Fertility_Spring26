# June 2024 CPS Fertility Stock Targets

This folder is the durable source for the calibration's completed-fertility
stock targets. It constructs the targets from the Census June 2024 Current
Population Survey Fertility Supplement public-use file.

The model key `tfr` is retained for compatibility, but the measured object is
not a period total fertility rate. It is mean children ever born among women
ages 40--44, a completed-fertility stock.

## Source and sample

- Source: June 2024 CPS Fertility Supplement `jun24pub.csv`.
- Primary sample: women (`PESEX=2`) ages 40--44 with a valid fertility
  supplement response (`PTSF1` in 0--5) and positive `PWSSWGT`.
- The age band corresponds approximately to birth years 1979--1984 depending
  on interview date and birthday. The public file does not provide an exact
  respondent birth year, so this is an age-window definition rather than an
  exact birth-cohort filter.
- Weight: `PWSSWGT`, labeled by Census as the second-stage, rake-6 final-step
  weight and suggested by Census for `PTSF1`. Its scale is immaterial for means
  and shares; dividing its raw integer value by 10,000 yields population units.
- `PTSF1=5` means five or more children. The build makes no uncapping
  adjustment.

The raw 139 MB CSV is not committed. `source_manifest.json` records its Census
URL, SHA-256 checksum, byte size, and acquisition date. The builder downloads
the file to the ignored `cache/` directory and verifies both size and checksum
before using it.

The committed `requirements.txt` pins the NumPy version used for the bootstrap.
The checked build used Python 3.9.6 and NumPy 2.0.2.

## Estimates

The committed target table contains:

- completed fertility, stored under the legacy model key `tfr`;
- childlessness;
- literal parity shares for 0, 1, 2, and 3+ children;
- the capped top-bin mean
  \(E[\min(N,5)\mid N\geq3]\).

The `model_key` column maps the capped top-bin mean to
`tfr_top_bin_weight`; the descriptive `moment_key` remains
`capped_top_bin_mean` so the censoring convention cannot be missed.

The primary sample reproduces completed fertility `1.918424850350` and
childlessness `0.188179593021`. These round to the Census-published values
`1.918` and `0.188`.

The builder uses a seeded person-level bootstrap with replacement. Each
resampled person retains her original `PWSSWGT`. The default is 1,000 draws
with seed `20260724`; the script refuses fewer than 1,000 draws. The covariance
file covers all seven reported estimates. Because childlessness and the
zero-parity share are the same statistic, their covariance rows are
intentionally identical and the full matrix is singular.

The sensitivity table reports the 40--44 primary window and the 41--45 and
43--47 shifts, each with weighted and unweighted estimates.

## First-birth timing variables

The 2024 supplement contains both `PTSF2` (year first child born) and `PTSAYFC`
(age of mother at first birth). The builder checks and reports their presence
but does not use either variable. Mean age at first birth and the share of first
births at age 30 or older remain owned by the E-strand NCHS cohort-table build.

## Reproduction

From the repository root:

```bash
python3 -m pip install -r code/data/cps_fertility/requirements.txt
python3 code/data/cps_fertility/build_cps_fertility_targets.py
```

This downloads and verifies the source if the cache is empty, then runs the
1,000-draw bootstrap and rewrites the committed outputs.

To force a fresh verified download:

```bash
python3 code/data/cps_fertility/build_cps_fertility_targets.py \
  --refresh-download
```

To use a previously downloaded file outside the managed cache:

```bash
python3 code/data/cps_fertility/build_cps_fertility_targets.py \
  --raw-file=/absolute/path/to/jun24pub.csv
```

The default outputs are:

- `output/cps_fertility_targets.csv`
- `output/cps_fertility_covariance.csv`
- `output/cps_fertility_sensitivities.csv`
- `output/build_diagnostics.json`

The E strand should consume `output/cps_fertility_targets.csv`. It should not
source timing targets from this folder.
