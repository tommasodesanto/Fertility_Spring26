#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", args[startsWith(args, file_flag)])
  if (length(script_path) == 0) {
    stop("Run this script with Rscript so --file= is available.")
  }
  dirname(normalizePath(script_path[1]))
}

safe_wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

fmt_num <- function(x) sprintf("%.3f", x)
fmt_pct <- function(x) sprintf("%.2f%%", 100 * x)

load_hpi_controls <- function(path) {
  if (!file.exists(path)) {
    warning(sprintf("Missing HPI file: %s", path))
    return(data.table(geo_code = character(), year = integer(), hpi = numeric()))
  }
  x <- as.data.table(read_dta(path))
  if (!all(c("GEOID", "year", "HPI") %in% names(x))) {
    stop(sprintf("HPI file has unexpected schema: %s", path))
  }
  x[, geo_code := sprintf("%05d", as.integer(GEOID))]
  x[, year := as.integer(year)]
  x <- x[year >= 2005 & year <= 2012 & is.finite(HPI)]
  x[, .(hpi = mean(HPI, na.rm = TRUE)), by = .(geo_code, year)]
}

load_unemp_controls <- function(path) {
  if (!file.exists(path)) {
    warning(sprintf("Missing unemployment file: %s", path))
    return(data.table(geo_code = character(), year = integer(), unemp_rate = numeric()))
  }
  x <- as.data.table(read_dta(path))
  if (!all(c("GEOID", "year", "annual_value") %in% names(x))) {
    stop(sprintf("Unemployment file has unexpected schema: %s", path))
  }
  x[, geo_code := sprintf("%05d", as.integer(GEOID))]
  x[, year := as.integer(year)]
  x <- x[year >= 2005 & year <= 2012 & is.finite(annual_value)]
  x[, .(unemp_rate = mean(annual_value, na.rm = TRUE)), by = .(geo_code, year)]
}

aggregate_acs_year <- function(path, year_value) {
  if (!file.exists(path)) {
    stop(sprintf("Missing ACS yearly file: %s", path))
  }
  message(sprintf("Reading ACS %d: %s", year_value, path))
  dt <- as.data.table(readRDS(path))
  required <- c(
    "year", "met2013", "gq", "ownershp", "valueh", "rooms", "pernum",
    "hhwt", "nchild", "relate", "age", "sex", "fertyr", "perwt"
  )
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(sprintf("ACS file missing columns: %s", paste(missing, collapse = ", ")))
  }

  hh <- dt[
    year == year_value &
      pernum == 1L &
      relate == 1L &
      gq == 1L &
      ownershp %in% c(1L, 2L) &
      age >= 20 & age <= 54 &
      !is.na(met2013) &
      met2013 > 0
  ]
  hh[, geo_code := sprintf("%05d", as.integer(met2013))]
  hh[, age_group := fcase(
    age >= 20 & age <= 24, "20_24",
    age >= 25 & age <= 34, "25_34",
    age >= 35 & age <= 44, "35_44",
    age >= 45 & age <= 54, "45_54"
  )]
  hh[, owner := as.integer(ownershp == 1L)]
  hh[, rooms_clean := fifelse(is.finite(rooms) & rooms > 0 & rooms < 90, as.numeric(rooms), NA_real_)]
  hh[, valid_home_value := fifelse(
    owner == 1L & is.finite(valueh) & valueh > 0 & valueh < 9999998,
    as.numeric(valueh),
    NA_real_
  )]

  housing <- hh[, .(
    panel = "housing_heads",
    observations = .N,
    weight = sum(hhwt, na.rm = TRUE),
    owner_rate = safe_wmean(owner, hhwt),
    mean_rooms = safe_wmean(rooms_clean, hhwt),
    share_rooms_ge_6 = safe_wmean(as.numeric(rooms_clean >= 6), hhwt),
    median_home_value = weighted_median(valid_home_value, hhwt)
  ), by = .(geo_code, year, age_group)]

  women <- dt[
    year == year_value &
      gq == 1L &
      sex == 2L &
      age >= 20 & age <= 44 &
      fertyr %in% c(1L, 2L) &
      !is.na(met2013) &
      met2013 > 0
  ]
  women[, geo_code := sprintf("%05d", as.integer(met2013))]
  women[, age_group := fcase(
    age >= 20 & age <= 24, "20_24",
    age >= 25 & age <= 34, "25_34",
    age >= 35 & age <= 44, "35_44"
  )]
  women[, birth_past_year := as.integer(fertyr == 2L)]
  fertility <- women[, .(
    panel = "fertility_women",
    observations = .N,
    weight = sum(perwt, na.rm = TRUE),
    birth_past_year_rate = safe_wmean(birth_past_year, perwt)
  ), by = .(geo_code, year, age_group)]

  list(housing = housing, fertility = fertility)
}

event_reg <- function(dt, panel_value, age_group_value, outcome, spec_name, treatment_var) {
  work <- dt[
    panel == panel_value &
      age_group == age_group_value &
      year >= 2005 & year <= 2012 &
      is.finite(get(outcome)) &
      is.finite(get(treatment_var)) &
      is.finite(weight) &
      weight > 0 &
      is.finite(log_hpi_2007_norm) &
      is.finite(unemp_rate)
  ]
  if (nrow(work) < 50) {
    return(data.table())
  }
  event_years <- setdiff(sort(unique(work$year)), 2007L)
  for (yy in event_years) {
    work[, paste0("treat_x_", yy) := get(treatment_var) * as.integer(year == yy)]
  }
  rhs <- paste(c(paste0("treat_x_", event_years), "log_hpi_2007_norm", "unemp_rate"), collapse = " + ")
  if (requireNamespace("fixest", quietly = TRUE)) {
    fit <- fixest::feols(
      as.formula(sprintf("%s ~ %s | geo_code + year", outcome, rhs)),
      data = work,
      weights = ~weight,
      cluster = ~geo_code,
      warn = FALSE,
      notes = FALSE
    )
    coefs <- as.data.table(fixest::coeftable(fit), keep.rownames = "term")
    note_text <- "Geography and year FE; weighted; controls for HPI relative to 2007 and unemployment; SE clustered by geography."
  } else {
    fit <- lm(
      as.formula(sprintf("%s ~ %s + factor(geo_code) + factor(year)", outcome, rhs)),
      data = work,
      weights = weight
    )
    coefs <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
    note_text <- "Geography and year FE; weighted; controls included; SE not clustered because fixest is unavailable."
  }
  coefs <- coefs[startsWith(term, "treat_x_")]
  if (!nrow(coefs)) {
    return(data.table())
  }
  setnames(coefs, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    c("estimate", "std_error", "t_value", "p_value")
  )
  coefs[, `:=`(
    panel = panel_value,
    age_group = age_group_value,
    spec = spec_name,
    treatment = treatment_var,
    outcome = outcome,
    baseline_year = 2007L,
    event_year = as.integer(sub("treat_x_", "", term)),
    observations = nrow(work),
    geographies = uniqueN(work$geo_code),
    note = note_text
  )]
  coefs[, .(
    panel, age_group, spec, treatment, outcome, baseline_year, event_year,
    estimate, std_error, t_value, p_value, observations, geographies, note
  )]
}

script_dir <- get_script_dir()
out_dir <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
acs_dir <- Sys.getenv(
  "ACS_YEARLY_RDS_DIR",
  unset = file.path(repo_root, "code", "data", "Spatial_aggregate_withmicrodata", "processed_data", "yearly_rds_v7")
)
year_text <- Sys.getenv("YEARS", unset = "2005,2006,2007,2008,2009,2010,2011,2012")
years <- as.integer(trimws(strsplit(year_text, ",")[[1]]))
if (any(is.na(years))) {
  stop(sprintf("Invalid YEARS value: %s", year_text))
}

hpi_path <- Sys.getenv(
  "HPI_PANEL",
  unset = "/Users/tommasodesanto/Desktop/Projects/Datasets/BirthsHybridMSAHpiPop_collapsed.dta"
)
unemp_path <- Sys.getenv(
  "UNEMP_PANEL",
  unset = "/Users/tommasodesanto/Desktop/Projects/Datasets/unemployment_rates/msa_unemployment_clean_nodups.dta"
)

year_data <- lapply(years, function(yy) {
  aggregate_acs_year(file.path(acs_dir, sprintf("fertility_microdata_%d.rds", yy)), yy)
})
housing <- rbindlist(lapply(year_data, `[[`, "housing"), fill = TRUE)
fertility <- rbindlist(lapply(year_data, `[[`, "fertility"), fill = TRUE)

exposure <- housing[
  year == 2007L &
    age_group == "25_34" &
    is.finite(median_home_value) &
    median_home_value > 0,
  .(
    median_home_value_2007 = weighted_median(median_home_value, weight),
    owner_observations_2007 = sum(observations, na.rm = TRUE),
    owner_weight_2007 = sum(weight, na.rm = TRUE)
  ),
  by = geo_code
]
exposure <- exposure[is.finite(median_home_value_2007) & median_home_value_2007 > 0]
exposure[, fthc_share_2007 := 8000 / median_home_value_2007]
exposure[, fthc_share_10pp := fthc_share_2007 / 0.10]
cutoff <- as.numeric(quantile(exposure$fthc_share_2007, 0.75, na.rm = TRUE))
exposure[, high_fthc_exposure := as.integer(fthc_share_2007 >= cutoff)]
fwrite(exposure, file.path(out_dir, "fthc_exposure_2007.csv"))

panel <- rbindlist(list(housing, fertility), fill = TRUE)
panel <- merge(panel, exposure, by = "geo_code", all.x = FALSE)
hpi <- load_hpi_controls(hpi_path)
unemp <- load_unemp_controls(unemp_path)
panel <- merge(panel, hpi, by = c("geo_code", "year"), all.x = TRUE)
panel <- merge(panel, unemp, by = c("geo_code", "year"), all.x = TRUE)
panel[, hpi_2007 := hpi[year == 2007][1], by = geo_code]
panel[, log_hpi_2007_norm := fifelse(
  is.finite(hpi) & is.finite(hpi_2007) & hpi > 0 & hpi_2007 > 0,
  log(hpi / hpi_2007),
  NA_real_
)]
setorder(panel, panel, age_group, geo_code, year)
fwrite(panel, file.path(out_dir, "fthc_acs_screen_panel_2005_2012.csv"))

results <- rbindlist(list(
  rbindlist(lapply(c("owner_rate", "mean_rooms", "share_rooms_ge_6"), event_reg,
    dt = panel,
    panel_value = "housing_heads",
    age_group_value = "25_34",
    spec_name = "continuous_credit_share_per_10pp",
    treatment_var = "fthc_share_10pp"
  ), fill = TRUE),
  rbindlist(lapply(c("owner_rate", "mean_rooms", "share_rooms_ge_6"), event_reg,
    dt = panel,
    panel_value = "housing_heads",
    age_group_value = "25_34",
    spec_name = "high_credit_share_top_quartile",
    treatment_var = "high_fthc_exposure"
  ), fill = TRUE),
  event_reg(
    dt = panel,
    panel_value = "fertility_women",
    age_group_value = "25_34",
    outcome = "birth_past_year_rate",
    spec_name = "continuous_credit_share_per_10pp",
    treatment_var = "fthc_share_10pp"
  ),
  event_reg(
    dt = panel,
    panel_value = "fertility_women",
    age_group_value = "25_34",
    outcome = "birth_past_year_rate",
    spec_name = "high_credit_share_top_quartile",
    treatment_var = "high_fthc_exposure"
  ),
  event_reg(
    dt = panel,
    panel_value = "fertility_women",
    age_group_value = "35_44",
    outcome = "birth_past_year_rate",
    spec_name = "placebo_35_44_high_credit_share_top_quartile",
    treatment_var = "high_fthc_exposure"
  )
), fill = TRUE)
fwrite(results, file.path(out_dir, "fthc_acs_event_study.csv"))

summary_bins <- panel[
  age_group == "25_34",
  .(
    geographies = uniqueN(geo_code),
    observations = sum(observations, na.rm = TRUE),
    weight = sum(weight, na.rm = TRUE),
    owner_rate = safe_wmean(owner_rate, weight),
    mean_rooms = safe_wmean(mean_rooms, weight),
    share_rooms_ge_6 = safe_wmean(share_rooms_ge_6, weight),
    birth_past_year_rate = safe_wmean(birth_past_year_rate, weight),
    median_home_value_2007 = weighted_median(median_home_value_2007, weight),
    fthc_share_2007 = weighted_median(fthc_share_2007, weight),
    log_hpi_2007_norm = safe_wmean(log_hpi_2007_norm, weight),
    unemp_rate = safe_wmean(unemp_rate, weight)
  ),
  by = .(panel, year, high_fthc_exposure)
]
fwrite(summary_bins, file.path(out_dir, "fthc_acs_summary_bins.csv"))

display_bins <- copy(summary_bins[panel %in% c("housing_heads", "fertility_women")])
display_bins[, high_fthc_exposure := fifelse(
  high_fthc_exposure == 1L,
  "high_credit_share",
  "other"
)]
display_bins[, `:=`(
  owner_rate = fmt_pct(owner_rate),
  mean_rooms = fmt_num(mean_rooms),
  share_rooms_ge_6 = fmt_pct(share_rooms_ge_6),
  birth_past_year_rate = fmt_pct(birth_past_year_rate),
  median_home_value_2007 = sprintf("%.0f", median_home_value_2007),
  fthc_share_2007 = fmt_pct(fthc_share_2007),
  log_hpi_2007_norm = fmt_num(log_hpi_2007_norm),
  unemp_rate = fmt_num(unemp_rate)
)]
setorder(display_bins, panel, high_fthc_exposure, year)

display_results <- copy(results[
  spec %in% c("continuous_credit_share_per_10pp", "high_credit_share_top_quartile", "placebo_35_44_high_credit_share_top_quartile")
])
display_results[, `:=`(
  estimate = fmt_num(estimate),
  std_error = fmt_num(std_error),
  p_value = fmt_num(p_value)
)]
setorder(display_results, panel, age_group, spec, outcome, event_year)

summary_lines <- c(
  "# First-Time Homebuyer Tax Credit ACS Screen",
  "",
  sprintf("Built at: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("ACS source directory: `%s`", acs_dir),
  sprintf("Years: `%s`", paste(years, collapse = ", ")),
  "",
  "## Design",
  "",
  "- Exposure is `8000 / median owner-reported home value in 2007`, computed from ACS owner household heads ages `25_34` by CBSA.",
  "- Continuous treatment is scaled per `10` percentage points of home value.",
  sprintf("- High exposure is the top quartile of the credit-share measure; cutoff = `%.4f`.", cutoff),
  "- Outcomes use ACS household heads ages `25_34` for housing and ACS women ages `25_34` for birth in past year.",
  "- Regressions include geography and year fixed effects, ACS weights, HPI relative to `2007`, unemployment, and geography-clustered SEs when `fixest` is available.",
  "",
  "## Read",
  "",
  "- Treat this as a screening design. The policy is national and coincides with the Great Recession.",
  "- A credible mechanism would need positive ownership or housing-space effects before using fertility as evidence.",
  "- ACS birth-in-past-year is a timing outcome, not completed fertility.",
  "",
  "## Ages 25-34 Summary",
  "",
  paste(capture.output(print(display_bins)), collapse = "\n"),
  "",
  "## Event-Study Diagnostics",
  "",
  paste(capture.output(print(display_results)), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "fthc_acs_summary.md"))

message("Wrote:")
message("  ", file.path(out_dir, "fthc_exposure_2007.csv"))
message("  ", file.path(out_dir, "fthc_acs_screen_panel_2005_2012.csv"))
message("  ", file.path(out_dir, "fthc_acs_summary_bins.csv"))
message("  ", file.path(out_dir, "fthc_acs_event_study.csv"))
message("  ", file.path(out_dir, "fthc_acs_summary.md"))
