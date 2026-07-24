#!/usr/bin/env Rscript

# Build cohort first-birth timing targets from the local NCHS natality archive.
# The only microdata read is the age/order pair needed for this construction.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(haven))

script_file <- function() {
  arg <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", arg, value = TRUE)
  if (length(hit) == 1L) return(normalizePath(sub("^--file=", "", hit)))
  normalizePath("build_first_birth_timing.R")
}

out_dir <- dirname(script_file())
archive_root <- Sys.getenv(
  "NCHS_NATALITY_ROOT",
  "/Users/tommasodesanto/Desktop/Projects/Datasets/natality_data"
)
force_rebuild <- "--force" %in% commandArgs(trailingOnly = TRUE)
years_required <- 1994L
ages <- 12:49

cache_path <- file.path(out_dir, "first_birth_counts_year_age.csv")
manifest_path <- file.path(out_dir, "first_birth_counts_manifest.csv")
exclusions_path <- file.path(out_dir, "order_exclusion_shares_by_year.csv")
targets_path <- file.path(out_dir, "timing_targets.csv")
truncation_path <- file.path(out_dir, "cohort_truncation.csv")
readme_path <- file.path(out_dir, "README.md")

find_us_files <- function(root) {
  candidates <- list.files(
    root, recursive = TRUE, full.names = TRUE,
    pattern = "^natality[0-9]{4}us\\.dta$"
  )
  year <- as.integer(sub("^natality([0-9]{4})us\\.dta$", "\\1", basename(candidates)))
  keep <- !is.na(year) & year >= years_required
  ans <- data.frame(year = year[keep], path = normalizePath(candidates[keep]),
                    stringsAsFactors = FALSE)
  ans <- ans[order(ans$year, ans$path), , drop = FALSE]
  if (anyDuplicated(ans$year)) stop("More than one US natality file found for a year.")
  if (!nrow(ans)) stop("No natalityYYYYus.dta files found from 1994 onward under: ", root)
  ans
}

make_manifest <- function(files) {
  data.frame(
    year = files$year,
    file = files$path,
    # Store byte sizes as text so a CSV round trip cannot change their type or
    # precision (and thereby spuriously invalidate the cache).
    size_bytes = format(file.info(files$path)$size, scientific = FALSE, trim = TRUE),
    stringsAsFactors = FALSE
  )
}

same_manifest <- function(current) {
  if (!file.exists(manifest_path) || !file.exists(cache_path) || !file.exists(exclusions_path)) return(FALSE)
  old <- read.csv(
    manifest_path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = c("integer", "character", "character")
  )
  identical(old, current)
}

variables_for_file <- function(path, year) {
  header <- read_dta(path, n_max = 0)
  available <- names(header)
  age <- if ("dmage" %in% available) "dmage" else if ("mager" %in% available) "mager" else if ("mager41" %in% available) "mager41" else NA_character_
  order <- if ("dlivord" %in% available) "dlivord" else if ("lbo_rec" %in% available) "lbo_rec" else NA_character_
  if (is.na(age) || is.na(order)) stop("Missing required age/order variable in ", path, " (year ", year, ").")
  c(age = age, order = order)
}

build_cache <- function(files, manifest) {
  counts <- list()
  exclusions <- list()
  for (i in seq_len(nrow(files))) {
    year <- files$year[i]
    path <- files$path[i]
    vars <- variables_for_file(path, year)
    message(sprintf("Reading %d (%s, %s) ...", year, vars["age"], vars["order"]))
    # `col_select` means haven reads only these two columns, not the full record.
    d <- read_dta(path, col_select = tidyselect::all_of(unname(vars)))
    age <- as.numeric(d[[vars["age"]]])
    order <- as.numeric(d[[vars["order"]]])
    valid_age <- !is.na(age) & age %in% ages
    unknown_order <- is.na(order) | order == 99
    first <- valid_age & !unknown_order & order == 1
    by_age <- tabulate(match(age[first], ages), nbins = length(ages))
    counts[[i]] <- data.frame(year = year, age = ages, n_first_births = by_age)
    exclusions[[i]] <- data.frame(
      year = year,
      n_birth_records_age_12_49 = sum(valid_age),
      n_unknown_live_birth_order = sum(valid_age & unknown_order),
      unknown_live_birth_order_share = mean(unknown_order[valid_age]),
      age_variable = vars["age"],
      order_variable = vars["order"],
      stringsAsFactors = FALSE
    )
    rm(d, age, order)
    gc(verbose = FALSE)
  }
  write.csv(do.call(rbind, counts), cache_path, row.names = FALSE)
  write.csv(do.call(rbind, exclusions), exclusions_path, row.names = FALSE)
  write.csv(manifest, manifest_path, row.names = FALSE)
}

weighted_mean_or_na <- function(x, w) if (sum(w) == 0) NA_real_ else sum(x * w) / sum(w)

make_targets <- function(counts) {
  cohort_counts <- transform(counts, cohort = year - age)
  # The integer assignment is conventional: births in year t to an age-a mother
  # are assigned to t-a. Dates within a calendar year create an unavoidable +/-1 blur.
  reference <- subset(cohort_counts, cohort >= 1970 & cohort <= 1974)
  windows <- data.frame(
    window = c("primary_1979_1984", "sensitivity_1975_1980"),
    cohort_min = c(1979L, 1975L), cohort_max = c(1984L, 1980L)
  )
  targets <- list()
  truncation <- list()
  for (i in seq_len(nrow(windows))) {
    w <- windows[i, ]
    x <- subset(cohort_counts, cohort >= w$cohort_min & cohort <= w$cohort_max)
    cohort_rows <- lapply(w$cohort_min:w$cohort_max, function(c) {
      z <- subset(x, cohort == c)
      max_age <- if (any(z$n_first_births > 0)) max(z$age[z$n_first_births > 0]) else NA_integer_
      above <- if (is.na(max_age)) NA_real_ else sum(reference$n_first_births[reference$age > max_age])
      ref_total <- sum(reference$n_first_births)
      data.frame(
        window = w$window, cohort = c, max_observed_age = max_age,
        n_observed_first_births = sum(z$n_first_births),
        reference_n_first_births = ref_total,
        reference_n_first_births_above_max_age = above,
        estimated_share_plausibly_missed = above / ref_total
      )
    })
    tr <- do.call(rbind, cohort_rows)
    truncation[[i]] <- tr
    targets[[i]] <- data.frame(
      window = w$window,
      cohort_range = sprintf("%d-%d", w$cohort_min, w$cohort_max),
      mean_age_first_birth = weighted_mean_or_na(x$age, x$n_first_births),
      share_30plus = if (sum(x$n_first_births) == 0) NA_real_ else sum(x$n_first_births[x$age >= 30]) / sum(x$n_first_births),
      n_first_births = sum(x$n_first_births),
      truncation_estimated_share = weighted_mean_or_na(
        tr$estimated_share_plausibly_missed, tr$n_observed_first_births
      )
    )
  }
  write.csv(do.call(rbind, targets), targets_path, row.names = FALSE)
  write.csv(do.call(rbind, truncation), truncation_path, row.names = FALSE)
}

write_readme <- function(files) {
  exclusions <- read.csv(exclusions_path, stringsAsFactors = FALSE)
  fmt <- function(x) format(x, scientific = FALSE, trim = TRUE)
  rows <- vapply(seq_len(nrow(exclusions)), function(i) sprintf(
    "| %d | %s | %s | %s | %.4f%% |",
    exclusions$year[i], exclusions$age_variable[i], exclusions$order_variable[i],
    fmt(exclusions$n_unknown_live_birth_order[i]),
    100 * exclusions$unknown_live_birth_order_share[i]
  ), character(1))
  lines <- c(
    "# NCHS natality first-birth timing targets",
    "",
    "This directory builds first-birth timing targets from the local NCHS natality microdata archive. The archive is not part of this repository and is read-only:",
    "`/Users/tommasodesanto/Desktop/Projects/Datasets/natality_data/YYYY/natalityYYYYus.dta`.",
    "Only `us` files are used; `ps` territory files are ignored. The currently detected archive spans 1994 through 2023.",
    "",
    "## Reproduce",
    "",
    "From the repository root, run:",
    "",
    "```sh",
    "Rscript code/data/nchs_natality_timing/build_first_birth_timing.R",
    "```",
    "",
    "Set `NCHS_NATALITY_ROOT` to use another local archive. Add `--force` to rebuild the cache. Otherwise the builder skips the heavy microdata pass only when `first_birth_counts_manifest.csv` exactly matches the current US-file list and byte sizes.",
    "",
    "## Construction",
    "",
    "A first birth is a record with live-birth order equal to 1. The builder retains mothers aged 12--49, excludes missing/unknown live-birth order (legacy code 99 or a missing recode), and counts first births by calendar year and single-year age. It reads only the two required columns from each `.dta` file through `haven::read_dta(..., col_select = ...)`; it never loads a full natality record file.",
    "",
    "| Data years | Mother's single-year age | Live-birth order |",
    "| --- | --- | --- |",
    "| 1994--2002 | `dmage` | `dlivord` |",
    "| 2003 | `mager41` | `lbo_rec` |",
    "| 2004--2023 | `mager` | `lbo_rec` |",
    "",
    "The cohort convention is `cohort = data_year - age`. Because both date of birth and age are measured within calendar years, this is an integer approximation with a mechanical +/-1-year cohort-assignment blur.",
    "",
    "`timing_targets.csv` reports pooled mothers-only first-birth moments for the primary 1979--1984 cohorts and the 1975--1980 sensitivity window. `mean_age_first_birth` and `share_30plus` are weighted by the first-birth counts in `first_birth_counts_year_age.csv`.",
    "",
    "## Order exclusions",
    "",
    "The table below reports unknown-order exclusions among all age-12--49 birth records. The same values are available in machine-readable form in `order_exclusion_shares_by_year.csv`.",
    "",
    "| Year | Age variable | Order variable | Unknown order records | Share of age-12--49 records |",
    "| --- | --- | --- | ---: | ---: |",
    rows,
    "",
    "## Truncation",
    "",
    "For each target cohort, `cohort_truncation.csv` reports the maximum observed first-birth age. Its estimated missed tail is the share of first births above that age in the pooled 1970--1974 reference cohorts. The window-level `truncation_estimated_share` in `timing_targets.csv` is the observed-first-birth-count-weighted average of those cohort estimates. This is a transparent reference-tail approximation, not a correction applied to the reported timing moments."
  )
  writeLines(lines, readme_path)
}

files <- find_us_files(archive_root)
manifest <- make_manifest(files)
if (force_rebuild || !same_manifest(manifest)) {
  message("Cache missing or archive manifest changed: building year-by-age counts.")
  build_cache(files, manifest)
} else {
  message("Archive manifest matches cache: skipping microdata pass.")
}
counts <- read.csv(cache_path, stringsAsFactors = FALSE)
make_targets(counts)
write_readme(files)

targets <- read.csv(targets_path, stringsAsFactors = FALSE)
primary <- subset(targets, window == "primary_1979_1984")
cat("\nPRIMARY 1979-1984 COHORT TARGETS\n")
cat(sprintf("Mean age at first birth: %.6f\n", primary$mean_age_first_birth))
cat(sprintf("Share of first births at age >= 30: %.6f\n", primary$share_30plus))
