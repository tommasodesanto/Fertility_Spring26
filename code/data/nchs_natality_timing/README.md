# NCHS natality first-birth timing targets

This directory builds first-birth timing targets from the local NCHS natality microdata archive. The archive is not part of this repository and is read-only:
`/Users/tommasodesanto/Desktop/Projects/Datasets/natality_data/YYYY/`.
Only national `us` files are used; `ps` territory and `_cwed` files are ignored. The target contract is deliberately frozen to data years 1987--2023. This start date observes every assigned-cohort age cell from age 12 through the age attainable in 2023 for both the 1975--1980 comparison cohorts and the 1979--1984 target cohorts; it does not remove right censoring for the younger cohorts. Later files cannot silently change the contract.

## Reproduce

From the repository root, run:

```sh
Rscript code/data/nchs_natality_timing/build_first_birth_timing.R
```

Set `NCHS_NATALITY_ROOT` to use another local archive. Add `--force` to rebuild the cache. Otherwise the builder skips the heavy microdata pass only when `first_birth_counts_manifest.csv` exactly matches the current US-file list, byte sizes, and raw-file SHA-256 hashes.

## Construction

A first birth is a record with live-birth order equal to 1. The builder retains mothers aged 12--49, excludes missing/unknown live-birth order (code 99 in legacy `dlivord`; code 9 in `lbo_rec`), and counts first births by calendar year and single-year age. It reads only the two required columns from each `.dta` file through `haven::read_dta(..., col_select = ...)`; it never loads a full natality record file.

| Data years | Mother's single-year age | Live-birth order |
| --- | --- | --- |
| 1987--2002 | `dmage` | `dlivord` |
| 2003 | `mager41`, decoded from the official 41-category recode | `lbo_rec` |
| 2004--2023 | `mager` | `lbo_rec` |

The 1989 certificate revision changes the source of maternal age toward completed age derived from date of birth when available; it does not change the completed-age concept or the definition of detail live-birth order. The annual first-birth distributions and unknown-order shares are smooth across this source-item transition, so no adjustment is applied. The relevant official documentation receipts are `Nat1987doc.pdf` (SHA-256 `145d8496d764ee66583ef1520ae56804e1c3373bea35d3a570ecbc1426bd5c79`) and `Nat1989doc.pdf` (SHA-256 `92dab8115baec71eec3633239cbd042b2079ad6b80bd1b3a3a43c3276ac3a7cb`).

The cohort convention is `cohort = data_year - age`. Because both date of birth and age are measured within calendar years, this is an integer approximation with a mechanical +/-1-year cohort-assignment blur.

The inherited cache previously treated 2003 `MAGER41` codes as literal ages. That was wrong: the official NCHS 2003 detail-file documentation defines code 01 as under age 15, code 02 as age 15, and subsequent codes in single years. The corrected cache decodes codes 02--41 as `code + 13`; code 01 is placed at 14 only to preserve its below-first-cutoff status. The target cohorts are ages 19--28 in 2003, where the recode is single-year exact. The official file is [Nat2003doc.pdf](https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/Nat2003doc.pdf), SHA-256 `82246197e30d54c56a69314bfdcb8f553e6ca0a0d509f8ce8112c2e996b5b2f5`.

`timing_targets.csv` preserves the original exact-age output for backward compatibility, and `timing_targets_exact_age.csv` gives the same reference moments plus their uncertainty receipt. Exact-age moments are not directly comparable with a model that records births only in four-year cells.

## Model-comparable timing contract

The recommended minimal-change operator uses every observed first birth. Exact ages below 22 map to the model's first cell (start age 18); ages 22--25 map to 22; ages 26--29 to 26; ages 30--33 to 30; ages 34--37 to 34; ages 38--41 to 38; and exact ages 42 or older to the terminal cell (start age 42). This is called `boundary_collapsed` in the machine-readable files.

Both data and the live model measurement are labeled by the continuous-age interval midpoints 20, 24, ..., 44. Relative to an archived period-start report, midpoint coding adds exactly 2 to both sides and therefore changes neither the model--data gap nor the loss. The age-30 threshold is aligned with a bin boundary: `bin_start >= 30` is identically the same classification as exact age at least 30 under both reported operators.

| Contract | Cohorts | Mean first-birth age | Share age 30+ | First births used |
| --- | --- | ---: | ---: | ---: |
| Recommended boundary collapse | 1979--1984 | 26.044627257483 | 0.260327401667 | 10578871 |
| Boundary-collapse comparison window | 1975--1980 | 25.919040268685 | 0.250397460163 | 10029307 |
| Common support, ages 18--45 | 1979--1984 | 26.665069712650 | 0.287048349933 | 9594098 |
| Common-support comparison window | 1975--1980 | 26.605067479802 | 0.279424130173 | 8976612 |
| Exact-age reference | 1979--1984 | 25.162728423477 | 0.260327401667 | 10578871 |
| Exact-age comparison window | 1975--1980 | 24.993733963872 | 0.250397460163 | 10029307 |

Conditioning on the literal common support 18--45 is retained as a sensitivity, not the recommendation. It excludes 984773 observed first births (9.3089 percent) from the primary window, almost entirely teen births, and therefore changes the age-30 share's denominator from 0.260327 to 0.287048. Boundary collapse instead preserves the unconditional observed-mothers denominator used alongside completed fertility and childlessness. Its cost is transparent coarsening: teen births are dated in the first model cell and births at 42 or older in the last. A future model with parity at age-18 entry could treat the teen-birth stock directly; the current target does not claim to do so.

The primary-versus-comparison-window spread is recomputed after applying each operator. Both cohort windows are observed from age 12 onward in the 1987--2023 cache, so this comparison is not contaminated by differential left censoring. Under the recommended boundary collapse the spread is 0.125587 year for the mean and 0.009930 for the age-30 share. Following the existing conservative convention, the mean spread is rounded up to the nearest 0.05 year and the share spread to the nearest 0.001, giving standard errors 0.150 and 0.010 and inverse-variance weights 44.444444 and 10000.000000. Under common support the corresponding standard errors are 0.100 and 0.008, with weights 100.000000 and 15625.000000. These are window-stability measures, not sampling standard errors.

Machine-readable artifacts:

- `timing_target_contract.csv`: long-form activation contract; the recommended rows are labeled `recommended_model_comparable`.
- `timing_targets_model_comparable.csv`: wide target, support, tail, uncertainty, and weight receipt for both model operators.
- `timing_age_bin_counts.csv`: every window-by-operator-by-bin count and share.
- `timing_targets_exact_age.csv`: exact-age reference and uncertainty receipt.
- `timing_targets.csv`: backward-compatible six-column exact-age reference.
- `timing_target_metadata.json`: fail-closed operator, source-document, file-hash, source-bundle, and contract-bundle provenance.

Activation requires changing the target contract and its fingerprint: use the recommended primary target 26.044627257483 for the midpoint-labeled mean and 0.260327401667 for the age-30 share, with standard errors 0.150 and 0.010 and weights 44.4444444444 and 10000. The live model measurement is already midpoint-coded and uses the aligned age-30 boundary. No structural model or fertility profile is changed by this builder.

## Order exclusions

The table below reports unknown-order exclusions among all age-12--49 birth records. The same values are available in machine-readable form in `order_exclusion_shares_by_year.csv`.

| Year | Age variable | Order variable | Unknown order records | Share of age-12--49 records |
| --- | --- | --- | ---: | ---: |
| 1987 | dmage | dlivord | 16012 | 0.4199% |
| 1988 | dmage | dlivord | 22050 | 0.5634% |
| 1989 | dmage | dlivord | 22384 | 0.5533% |
| 1990 | dmage | dlivord | 26577 | 0.6384% |
| 1991 | dmage | dlivord | 20722 | 0.5035% |
| 1992 | dmage | dlivord | 18489 | 0.4543% |
| 1993 | dmage | dlivord | 18110 | 0.4522% |
| 1994 | dmage | dlivord | 22917 | 0.5792% |
| 1995 | dmage | dlivord | 27993 | 0.7172% |
| 1996 | dmage | dlivord | 23253 | 0.5970% |
| 1997 | dmage | dlivord | 20680 | 0.5324% |
| 1998 | dmage | dlivord | 26308 | 0.6669% |
| 1999 | dmage | dlivord | 18265 | 0.4609% |
| 2000 | dmage | dlivord | 17725 | 0.4362% |
| 2001 | dmage | dlivord | 13940 | 0.3458% |
| 2002 | dmage | dlivord | 10256 | 0.2547% |
| 2003 | mager41 | lbo_rec | 11986 | 0.2926% |
| 2004 | mager | lbo_rec | 18692 | 0.4539% |
| 2005 | mager | lbo_rec | 17530 | 0.4229% |
| 2006 | mager | lbo_rec | 23748 | 0.5558% |
| 2007 | mager | lbo_rec | 21309 | 0.4929% |
| 2008 | mager | lbo_rec | 26599 | 0.6252% |
| 2009 | mager | lbo_rec | 26591 | 0.6427% |
| 2010 | mager | lbo_rec | 29129 | 0.7270% |
| 2011 | mager | lbo_rec | 28444 | 0.7182% |
| 2012 | mager | lbo_rec | 19933 | 0.5033% |
| 2013 | mager | lbo_rec | 18523 | 0.4701% |
| 2014 | mager | lbo_rec | 20012 | 0.5006% |
| 2015 | mager | lbo_rec | 20513 | 0.5144% |
| 2016 | mager | lbo_rec | 15956 | 0.4034% |
| 2017 | mager | lbo_rec | 9593 | 0.2483% |
| 2018 | mager | lbo_rec | 10724 | 0.2822% |
| 2019 | mager | lbo_rec | 9022 | 0.2402% |
| 2020 | mager | lbo_rec | 8311 | 0.2297% |
| 2021 | mager | lbo_rec | 12339 | 0.3363% |
| 2022 | mager | lbo_rec | 12016 | 0.3270% |
| 2023 | mager | lbo_rec | 9985 | 0.2771% |

## Truncation

For each target cohort, `cohort_truncation.csv` reports the maximum observed first-birth age. Its estimated missed tail is the share of first births above that age in the pooled 1970--1974 reference cohorts. Under the minimal 1987 start, that reference pool lacks some early-age cells, so this is an uncorrected right-tail diagnostic rather than an input to either target or weight. The primary boundary-collapse estimate is 0.6891 percent; the corresponding common-support diagnostic is 0.6651 percent. Neither target is corrected or imputed. An observed count of zero above age 45 in the primary window does not mean the true tail is zero: later ages for these cohorts occur after the fixed 2023 data endpoint.
