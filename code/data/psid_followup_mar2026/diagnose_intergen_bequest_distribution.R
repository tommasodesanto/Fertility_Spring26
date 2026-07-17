#!/usr/bin/env Rscript

# Diagnostic decomposition of the three internally calibrated estate moments.
# These outputs are not calibration targets.  They separate the wealth
# numerator, income denominator, housing component, fertility composition, and
# marital composition on the exact PSID samples used by the target audit.

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

weighted_quantile_safe <- function(x, w, prob) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= prob)[1L]]
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

weighted_cov_safe <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  y <- y[ok]
  w <- w[ok]
  mx <- weighted.mean(x, w)
  my <- weighted.mean(y, w)
  sum(w * (x - mx) * (y - my)) / sum(w)
}

weighted_cor_safe <- function(x, y, w) {
  cxy <- weighted_cov_safe(x, y, w)
  vx <- weighted_cov_safe(x, x, w)
  vy <- weighted_cov_safe(y, y, w)
  if (!is.finite(cxy) || vx <= 0 || vy <= 0) return(NA_real_)
  cxy / sqrt(vx * vy)
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[startsWith(args, "--file=")]
  if (length(file_arg) == 0L) stop("Run this file with Rscript.")
  dirname(normalizePath(sub("--file=", "", file_arg[1L], fixed = TRUE)))
}

root <- normalizePath(file.path(script_dir(), "..", "..", ".."))
psid_path <- file.path(dirname(root), "PSID", "PSIDSHELF_MOBILITY.dta")
outdir <- file.path(
  root, "code", "data", "psid_followup_mar2026", "output",
  "intergen_bequest_distribution_diagnostic"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(psid_path)) stop("Missing PSID shelf: ", psid_path)

message("Reading the minimal PSID diagnostic columns...")
raw <- as.data.table(read_dta(
  psid_path,
  col_select = c(
    "ID", "year", "AGEREP", "DEATHYEAR", "RELTOHEAD_", "RELCHINUM",
    "INCFAMR", "NETWORTHR", "NETWORTH2R", "HOMEEQUITYR", "IW",
    "FAMMARRIED"
  )
))

dt <- raw[, .(
  id = as.numeric(ID),
  year = as.integer(year),
  age = as.numeric(AGEREP),
  death_year = as.numeric(DEATHYEAR),
  reference_person = as.numeric(RELTOHEAD_) == 10,
  n_children = as.numeric(RELCHINUM),
  income = as.numeric(INCFAMR),
  total_estate = as.numeric(NETWORTHR),
  nonhousing = as.numeric(NETWORTH2R),
  housing = as.numeric(HOMEEQUITYR),
  weight = as.numeric(IW),
  married = as.numeric(FAMMARRIED)
)]
rm(raw)

# Use one common complete-case sample for every diagnostic so numerator and
# denominator comparisons cannot change merely because the sample changed.
dt <- dt[
  year >= 1984 & year <= 2019 &
    (is.na(death_year) | year <= death_year) &
    reference_person & age >= 65 & age <= 84 &
    is.finite(n_children) & is.finite(income) & income > 1000 &
    is.finite(total_estate) & is.finite(nonhousing) & is.finite(housing) &
    is.finite(weight) & weight > 0
]
dt[, estate_income_ratio := total_estate / income]

windows <- list(ages_65_75 = c(65, 75), ages_76_84 = c(76, 84))
groups <- list(
  all = function(x) rep(TRUE, nrow(x)),
  children_0 = function(x) x$n_children == 0,
  children_1 = function(x) x$n_children == 1,
  children_2plus = function(x) x$n_children >= 2,
  married = function(x) x$married == 1,
  not_married = function(x) x$married == 0,
  children_1_married = function(x) x$n_children == 1 & x$married == 1,
  children_2plus_married = function(x) x$n_children >= 2 & x$married == 1,
  children_1_not_married = function(x) x$n_children == 1 & x$married == 0,
  children_2plus_not_married = function(x) x$n_children >= 2 & x$married == 0
)

summarize_group <- function(block, window_name, group_name, benchmark_income_median, window_weight) {
  w <- block$weight
  q <- function(x, p) weighted_quantile_safe(x, w, p)
  estate_p90 <- q(block$total_estate, 0.90)
  top <- is.finite(block$total_estate) & block$total_estate >= estate_p90
  top_housing_share <- if (any(top)) {
    sum(w[top] * block$housing[top]) / sum(w[top] * block$total_estate[top])
  } else {
    NA_real_
  }
  data.table(
    dataset = "PSID",
    window = window_name,
    group = group_name,
    n_rows = nrow(block),
    n_people = uniqueN(block$id),
    mass_or_weight = sum(w),
    share_of_window_weight = sum(w) / window_weight,
    benchmark_income_median = benchmark_income_median,
    income_p10 = q(block$income, 0.10),
    income_p50 = q(block$income, 0.50),
    income_p90 = q(block$income, 0.90),
    total_estate_p50 = q(block$total_estate, 0.50),
    total_estate_p75 = q(block$total_estate, 0.75),
    total_estate_p90 = estate_p90,
    nonhousing_p50 = q(block$nonhousing, 0.50),
    nonhousing_p75 = q(block$nonhousing, 0.75),
    nonhousing_p90 = q(block$nonhousing, 0.90),
    housing_p50 = q(block$housing, 0.50),
    housing_p75 = q(block$housing, 0.75),
    housing_p90 = q(block$housing, 0.90),
    estate_income_ratio_p50 = q(block$estate_income_ratio, 0.50),
    estate_income_ratio_p75 = q(block$estate_income_ratio, 0.75),
    estate_income_ratio_p90 = q(block$estate_income_ratio, 0.90),
    estate_income_ratio_p90_p50 =
      q(block$estate_income_ratio, 0.90) / q(block$estate_income_ratio, 0.50),
    total_estate_income_covariance = weighted_cov_safe(block$total_estate, block$income, w),
    normalized_estate_income_covariance =
      weighted_cov_safe(block$total_estate, block$income, w) / benchmark_income_median^2,
    estate_income_correlation = weighted_cor_safe(block$total_estate, block$income, w),
    asinh_estate_log_income_correlation = weighted_cor_safe(
      asinh(block$total_estate / benchmark_income_median),
      log(block$income / benchmark_income_median),
      w
    ),
    mean_housing_share =
      weighted_mean_safe(block$housing, w) / weighted_mean_safe(block$total_estate, w),
    top_estate_decile_housing_share = top_housing_share
  )
}

rows <- list()
sample_rows <- list()
for (window_name in names(windows)) {
  bounds <- windows[[window_name]]
  window <- dt[age >= bounds[1] & age <= bounds[2]]
  benchmark_income_median <- weighted_quantile_safe(window$income, window$weight, 0.50)
  window_weight <- sum(window$weight)
  sample_rows[[window_name]] <- data.table(
    window = window_name,
    age_min = bounds[1],
    age_max = bounds[2],
    n_rows = nrow(window),
    n_people = uniqueN(window$id),
    weighted_mean_age = weighted_mean_safe(window$age, window$weight),
    benchmark_income_median = benchmark_income_median
  )
  for (group_name in names(groups)) {
    keep <- groups[[group_name]](window)
    keep[is.na(keep)] <- FALSE
    block <- window[keep]
    if (nrow(block) == 0L) next
    rows[[length(rows) + 1L]] <- summarize_group(
      block, window_name, group_name, benchmark_income_median, window_weight
    )
  }
}

result <- rbindlist(rows, fill = TRUE)
sample_counts <- rbindlist(sample_rows, fill = TRUE)

# Add comparable level normalizations.  Raw PSID columns remain in 2022 USD.
level_columns <- c(
  "income_p10", "income_p50", "income_p90",
  "total_estate_p50", "total_estate_p75", "total_estate_p90",
  "nonhousing_p50", "nonhousing_p75", "nonhousing_p90",
  "housing_p50", "housing_p75", "housing_p90"
)
for (column in level_columns) {
  result[, paste0(column, "_over_benchmark_income") :=
    get(column) / benchmark_income_median]
}

fwrite(result, file.path(outdir, "psid_distribution_decomposition.csv"))
fwrite(sample_counts, file.path(outdir, "sample_counts.csv"))
writeLines(c(
  "# PSID estate-distribution diagnostic",
  "",
  "These outputs are diagnostics, not proposed calibration moments.",
  "",
  "Sample: PSID reference persons, 1984--2019, alive in the observation year,",
  "real family income above $1,000, positive longitudinal weight, and complete",
  "total net worth, nonhousing net worth, home equity, income, and completed-child",
  "records. Dollar values are real 2022 USD. The two age windows exactly match the",
  "internally calibrated estate targets. Pooled family-years and IW weights match",
  "the existing target construction.",
  "",
  "Quantiles of total, nonhousing, and housing wealth are reported separately but",
  "are not added: component quantiles do not mechanically sum to total-wealth",
  "quantiles. `top_estate_decile_housing_share` is the aggregate housing share",
  "within the top total-estate decile."
), file.path(outdir, "README.md"))

message("Wrote PSID decomposition to ", outdir)
