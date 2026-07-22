#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
})

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

weighted_quantile_safe <- function(x, w, prob) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= prob)[1]]
}

script_dir <- function() {
  flag <- commandArgs(trailingOnly = FALSE)
  flag <- flag[startsWith(flag, "--file=")]
  if (length(flag) == 0L) stop("Run with Rscript.")
  dirname(normalizePath(sub("--file=", "", flag[1L], fixed = TRUE)))
}

root <- normalizePath(file.path(script_dir(), "..", "..", ".."))
extract_path <- file.path(root, "code", "data", "Spatial_aggregate_withmicrodata", "raw_data", "extract27.dta")
lookup_dir <- file.path(root, "code", "data", "mms_center_periphery", "data")
out_dir <- file.path(root, "code", "data", "mms_center_periphery", "output_intergen_calibration_mapping")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lookup <- bind_rows(
  read_csv(file.path(lookup_dir, "puma_mms_lookup_2010.csv"), show_col_types = FALSE) %>%
    transmute(lookup_period = "pre2022", statefip, puma, met2013 = cbsacode, mms_location),
  read_csv(file.path(lookup_dir, "puma_mms_lookup_2020.csv"), show_col_types = FALSE) %>%
    transmute(lookup_period = "post2021", statefip, puma, met2013 = cbsacode, mms_location)
)

message("Reading ACS/IPUMS columns...")
df <- read_dta(
  extract_path,
  col_select = c(year, statefip, puma, met2013, gq, pernum, hhwt,
                 ownershp, rent, hhincome, nchild, age)
) %>%
  transmute(
    year = as.integer(year), statefip = as.integer(statefip), puma = as.integer(puma),
    met2013 = as.integer(met2013), gq = as.integer(gq), pernum = as.integer(pernum),
    hhwt = as.numeric(hhwt), ownershp = as.integer(ownershp), rent = as.numeric(rent),
    income = as.numeric(hhincome), nchild = as.numeric(nchild), age = as.integer(age),
    lookup_period = if_else(year <= 2021, "pre2022", "post2021")
  ) %>%
  filter(year >= 2012, met2013 > 0, gq %in% c(1, 2), pernum == 1,
         age >= 30, age <= 55, ownershp == 2, nchild == 0,
         is.finite(hhwt), hhwt > 0, rent > 0, income > 1000) %>%
  left_join(lookup, by = c("lookup_period", "statefip", "puma", "met2013")) %>%
  mutate(
    mms_location = if_else(mms_location == "middle", "center", mms_location),
    annual_contract_rent = 12 * rent,
    contract_rent_share = annual_contract_rent / income
  ) %>%
  filter(mms_location %in% c("center", "periphery"))

income_lo <- weighted_quantile_safe(df$income, df$hhwt, 0.01)
income_hi <- weighted_quantile_safe(df$income, df$hhwt, 0.99)
share_lo <- weighted_quantile_safe(df$contract_rent_share, df$hhwt, 0.01)
share_hi <- weighted_quantile_safe(df$contract_rent_share, df$hhwt, 0.99)

est <- df %>%
  filter(between(income, income_lo, income_hi),
         between(contract_rent_share, share_lo, share_hi))

fit <- lm(contract_rent_share ~ log(income), data = est, weights = hhwt)
fit_year <- lm(contract_rent_share ~ log(income) + factor(year), data = est, weights = hhwt)

cuts <- vapply(1:4, function(q) weighted_quantile_safe(df$income, df$hhwt, q / 5), numeric(1))
df$income_quintile <- findInterval(df$income, cuts) + 1L
cells <- df %>%
  group_by(income_quintile) %>%
  summarise(
    n = n(), weight = sum(hhwt),
    mean_income = weighted_mean_safe(income, hhwt),
    median_income = weighted_quantile_safe(income, hhwt, 0.5),
    mean_contract_rent = weighted_mean_safe(annual_contract_rent, hhwt),
    mean_individual_share = weighted_mean_safe(contract_rent_share, hhwt),
    median_individual_share = weighted_quantile_safe(contract_rent_share, hhwt, 0.5),
    aggregate_share = sum(hhwt * annual_contract_rent) / sum(hhwt * income),
    .groups = "drop"
  )

moments <- tibble(
  moment = c(
    "childless_renter_contract_rent_share_aggregate",
    "childless_renter_contract_rent_share_mean",
    "childless_renter_contract_rent_share_median",
    "childless_renter_contract_rent_share_top_income_quintile",
    "childless_renter_contract_rent_share_log_income_slope",
    "childless_renter_contract_rent_share_log_income_slope_year_fe"
  ),
  value = c(
    sum(df$hhwt * df$annual_contract_rent) / sum(df$hhwt * df$income),
    weighted_mean_safe(df$contract_rent_share, df$hhwt),
    weighted_quantile_safe(df$contract_rent_share, df$hhwt, 0.5),
    cells$aggregate_share[cells$income_quintile == 5],
    unname(coef(fit)["log(income)"]),
    unname(coef(fit_year)["log(income)"])
  ),
  n = c(rep(nrow(df), 4), rep(nrow(est), 2)),
  weight = c(rep(sum(df$hhwt), 4), rep(sum(est$hhwt), 2)),
  sample = "ACS/IPUMS 2012-2023; MMS metros; middle->center; household heads 30-55; childless renters; positive contract rent; household income >1000",
  definition = c(
    "sum(weight*12*monthly contract rent)/sum(weight*household income)",
    "weighted mean of household 12*contract rent/income",
    "weighted median of household 12*contract rent/income",
    "aggregate contract-rent share in weighted top household-income quintile",
    "WLS coefficient in contract-rent share on log household income; income and share trimmed p1-p99",
    "same WLS coefficient with survey-year fixed effects"
  ),
  status = "candidate; contract rent excludes utilities and is not total housing expenditure"
)

write_csv(moments, file.path(out_dir, "acs_childless_renter_engel_moments.csv"))
write_csv(cells, file.path(out_dir, "acs_childless_renter_engel_income_cells.csv"))

writeLines(c(
  "# Childless-renter contract-rent Engel targets",
  "",
  "These are candidate calibration moments for the Stone--Geary housing block.",
  "The ACS variable is monthly contract rent, not gross rent: utilities are excluded.",
  "The aggregate or high-income share can discipline alpha_cons; the negative income",
  "gradient can discipline c_bar_0 jointly with the childless-renter rooms target.",
  "Do not label these objects total housing-expenditure shares.",
  "",
  "Outputs:",
  "- acs_childless_renter_engel_moments.csv",
  "- acs_childless_renter_engel_income_cells.csv"
), file.path(out_dir, "README.md"))

message("Wrote childless-renter Engel targets to ", out_dir)
