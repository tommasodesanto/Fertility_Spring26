# NCHS natality first-birth timing targets

This directory builds first-birth timing targets from the local NCHS natality microdata archive. The archive is not part of this repository and is read-only:
`/Users/tommasodesanto/Desktop/Projects/Datasets/natality_data/YYYY/natalityYYYYus.dta`.
Only `us` files are used; `ps` territory files are ignored. The currently detected archive spans 1994 through 2023.

## Reproduce

From the repository root, run:

```sh
Rscript code/data/nchs_natality_timing/build_first_birth_timing.R
```

Set `NCHS_NATALITY_ROOT` to use another local archive. Add `--force` to rebuild the cache. Otherwise the builder skips the heavy microdata pass only when `first_birth_counts_manifest.csv` exactly matches the current US-file list and byte sizes.

## Construction

A first birth is a record with live-birth order equal to 1. The builder retains mothers aged 12--49, excludes missing/unknown live-birth order (legacy code 99 or a missing recode), and counts first births by calendar year and single-year age. It reads only the two required columns from each `.dta` file through `haven::read_dta(..., col_select = ...)`; it never loads a full natality record file.

| Data years | Mother's single-year age | Live-birth order |
| --- | --- | --- |
| 1994--2002 | `dmage` | `dlivord` |
| 2003 | `mager41` | `lbo_rec` |
| 2004--2023 | `mager` | `lbo_rec` |

The cohort convention is `cohort = data_year - age`. Because both date of birth and age are measured within calendar years, this is an integer approximation with a mechanical +/-1-year cohort-assignment blur.

`timing_targets.csv` reports pooled mothers-only first-birth moments for the primary 1979--1984 cohorts and the 1975--1980 sensitivity window. `mean_age_first_birth` and `share_30plus` are weighted by the first-birth counts in `first_birth_counts_year_age.csv`.

## Order exclusions

The table below reports unknown-order exclusions among all age-12--49 birth records. The same values are available in machine-readable form in `order_exclusion_shares_by_year.csv`.

| Year | Age variable | Order variable | Unknown order records | Share of age-12--49 records |
| --- | --- | --- | ---: | ---: |
| 1994 | dmage | dlivord | 22917 | 0.5792% |
| 1995 | dmage | dlivord | 27993 | 0.7172% |
| 1996 | dmage | dlivord | 23253 | 0.5970% |
| 1997 | dmage | dlivord | 20680 | 0.5324% |
| 1998 | dmage | dlivord | 26308 | 0.6669% |
| 1999 | dmage | dlivord | 18265 | 0.4609% |
| 2000 | dmage | dlivord | 17725 | 0.4362% |
| 2001 | dmage | dlivord | 13940 | 0.3458% |
| 2002 | dmage | dlivord | 10256 | 0.2547% |
| 2003 | mager41 | lbo_rec | 0 | 0.0000% |
| 2004 | mager | lbo_rec | 0 | 0.0000% |
| 2005 | mager | lbo_rec | 0 | 0.0000% |
| 2006 | mager | lbo_rec | 0 | 0.0000% |
| 2007 | mager | lbo_rec | 0 | 0.0000% |
| 2008 | mager | lbo_rec | 0 | 0.0000% |
| 2009 | mager | lbo_rec | 0 | 0.0000% |
| 2010 | mager | lbo_rec | 0 | 0.0000% |
| 2011 | mager | lbo_rec | 0 | 0.0000% |
| 2012 | mager | lbo_rec | 0 | 0.0000% |
| 2013 | mager | lbo_rec | 0 | 0.0000% |
| 2014 | mager | lbo_rec | 0 | 0.0000% |
| 2015 | mager | lbo_rec | 0 | 0.0000% |
| 2016 | mager | lbo_rec | 0 | 0.0000% |
| 2017 | mager | lbo_rec | 0 | 0.0000% |
| 2018 | mager | lbo_rec | 0 | 0.0000% |
| 2019 | mager | lbo_rec | 0 | 0.0000% |
| 2020 | mager | lbo_rec | 0 | 0.0000% |
| 2021 | mager | lbo_rec | 0 | 0.0000% |
| 2022 | mager | lbo_rec | 0 | 0.0000% |
| 2023 | mager | lbo_rec | 0 | 0.0000% |

## Truncation

For each target cohort, `cohort_truncation.csv` reports the maximum observed first-birth age. Its estimated missed tail is the share of first births above that age in the pooled 1970--1974 reference cohorts. The window-level `truncation_estimated_share` in `timing_targets.csv` is the observed-first-birth-count-weighted average of those cohort estimates. This is a transparent reference-tail approximation, not a correction applied to the reported timing moments.
