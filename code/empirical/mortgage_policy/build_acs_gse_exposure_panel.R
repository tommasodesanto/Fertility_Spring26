#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
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

fmt_pct <- function(x) {
  sprintf("%.2f%%", 100 * x)
}

fmt_num <- function(x) {
  sprintf("%.3f", x)
}

aggregate_acs_year <- function(path, year_value) {
  if (!file.exists(path)) {
    stop(sprintf("Missing ACS yearly file: %s", path))
  }
  message(sprintf("Reading ACS %d: %s", year_value, path))
  dt <- as.data.table(readRDS(path))
  required <- c(
    "year", "met2013", "gq", "ownershp", "rent", "hhincome", "rooms",
    "pernum", "hhwt", "nchild", "relate", "age"
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
  if (!nrow(hh)) {
    stop(sprintf("No eligible household heads in ACS %d.", year_value))
  }

  hh[, geo_code := sprintf("%05d", as.integer(met2013))]
  hh[, age_group := fcase(
    age >= 20 & age <= 24, "20_24",
    age >= 25 & age <= 34, "25_34",
    age >= 35 & age <= 44, "35_44",
    age >= 45 & age <= 54, "45_54"
  )]
  hh[, owner := as.integer(ownershp == 1L)]
  hh[, renter := as.integer(ownershp == 2L)]
  hh[, has_children := as.integer(!is.na(nchild) & nchild > 0)]
  hh[, rooms_clean := fifelse(is.finite(rooms) & rooms > 0 & rooms < 90, as.numeric(rooms), NA_real_)]
  hh[, rent_per_room := fifelse(
    renter == 1L & is.finite(rent) & rent > 0 & is.finite(rooms_clean) & rooms_clean > 0,
    rent / rooms_clean,
    NA_real_
  )]

  out <- hh[, .(
    households = .N,
    hh_weight = sum(hhwt, na.rm = TRUE),
    owner_rate = safe_wmean(owner, hhwt),
    renter_rate = safe_wmean(renter, hhwt),
    parent_share = safe_wmean(has_children, hhwt),
    mean_rooms = safe_wmean(rooms_clean, hhwt),
    mean_rooms_owner = safe_wmean(rooms_clean[owner == 1L], hhwt[owner == 1L]),
    mean_rooms_renter = safe_wmean(rooms_clean[renter == 1L], hhwt[renter == 1L]),
    share_rooms_ge_6 = safe_wmean(as.numeric(rooms_clean >= 6), hhwt),
    share_rooms_ge_7 = safe_wmean(as.numeric(rooms_clean >= 7), hhwt),
    share_rooms_ge_8 = safe_wmean(as.numeric(rooms_clean >= 8), hhwt),
    mean_rent_per_room = safe_wmean(rent_per_room, hhwt)
  ), by = .(geo_code, year, age_group)]
  out[]
}

estimate_event_study <- function(panel, outcome) {
  work <- panel[
    age_group == "25_34" &
      year >= 2005 & year <= 2012 &
      is.finite(get(outcome)) &
      is.finite(gse_log_change_2009_vs_2007_proxy_max) &
      is.finite(hh_weight) &
      hh_weight > 0
  ]
  if (nrow(work) < 50) {
    return(data.table())
  }
  event_years <- setdiff(sort(unique(work$year)), 2007L)
  for (yy in event_years) {
    work[, paste0("exposure_x_", yy) := gse_log_change_2009_vs_2007_proxy_max * as.integer(year == yy)]
  }
  rhs <- paste(c(
    paste0("exposure_x_", event_years)
  ), collapse = " + ")
  if (requireNamespace("fixest", quietly = TRUE)) {
    fit <- fixest::feols(
      as.formula(sprintf("%s ~ %s | geo_code + year", outcome, rhs)),
      data = work,
      weights = ~hh_weight,
      cluster = ~geo_code,
      warn = FALSE,
      notes = FALSE
    )
    coefs <- as.data.table(fixest::coeftable(fit), keep.rownames = "term")
    note_text <- "Weighted metro and year FE diagnostic; standard errors clustered by metro."
  } else {
    lm_rhs <- paste(c("factor(geo_code)", "factor(year)", paste0("exposure_x_", event_years)), collapse = " + ")
    fit <- lm(
      as.formula(sprintf("%s ~ %s", outcome, lm_rhs)),
      data = work,
      weights = hh_weight
    )
    coefs <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
    note_text <- "Weighted metro and year FE diagnostic; standard errors are not clustered."
  }
  coefs <- coefs[startsWith(term, "exposure_x_")]
  if (!nrow(coefs)) {
    return(data.table())
  }
  setnames(coefs, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    c("estimate", "std_error", "t_value", "p_value")
  )
  coefs[, `:=`(
    outcome = outcome,
    event_year = as.integer(sub("exposure_x_", "", term)),
    baseline_year = 2007L,
    note = note_text
  )]
  coefs[, .(outcome, baseline_year, event_year, estimate, std_error, t_value, p_value, note)]
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

exposure_path <- file.path(out_dir, "loan_limit_exposure_cbsa_2008_2009.csv")
if (!file.exists(exposure_path)) {
  stop("Missing CBSA exposure panel. Run build_loan_limit_exposure_panel.R first.")
}
exposure <- fread(exposure_path, colClasses = c(msa_code = "character"))
positive_exposure <- exposure[gse_log_change_2009_vs_2007_proxy_max > 0, gse_log_change_2009_vs_2007_proxy_max]
if (!length(positive_exposure)) {
  stop("No positive GSE CBSA exposure values found.")
}
treated_median <- median(positive_exposure, na.rm = TRUE)
exposure[, exposure_bin := fcase(
  gse_log_change_2009_vs_2007_proxy_max <= 0, "untreated",
  gse_log_change_2009_vs_2007_proxy_max < treated_median, "treated_below_median",
  default = "treated_above_median"
)]

year_panels <- lapply(years, function(yy) {
  aggregate_acs_year(file.path(acs_dir, sprintf("fertility_microdata_%d.rds", yy)), yy)
})
acs_panel <- rbindlist(year_panels, fill = TRUE)

merged <- exposure[acs_panel, on = c(msa_code = "geo_code"), nomatch = 0]
setnames(merged, "msa_code", "geo_code")
setorder(merged, geo_code, age_group, year)
fwrite(merged, file.path(out_dir, "acs_gse_exposure_panel_2005_2012.csv"))

summary <- merged[, .(
  geographies = uniqueN(geo_code),
  households = sum(households, na.rm = TRUE),
  hh_weight = sum(hh_weight, na.rm = TRUE),
  owner_rate = safe_wmean(owner_rate, hh_weight),
  renter_rate = safe_wmean(renter_rate, hh_weight),
  parent_share = safe_wmean(parent_share, hh_weight),
  mean_rooms = safe_wmean(mean_rooms, hh_weight),
  mean_rooms_owner = safe_wmean(mean_rooms_owner, hh_weight),
  mean_rooms_renter = safe_wmean(mean_rooms_renter, hh_weight),
  share_rooms_ge_6 = safe_wmean(share_rooms_ge_6, hh_weight),
  share_rooms_ge_7 = safe_wmean(share_rooms_ge_7, hh_weight),
  share_rooms_ge_8 = safe_wmean(share_rooms_ge_8, hh_weight),
  mean_rent_per_room = safe_wmean(mean_rent_per_room, hh_weight),
  mean_gse_log_exposure = mean(gse_log_change_2009_vs_2007_proxy_max, na.rm = TRUE)
), by = .(year, age_group, exposure_bin)]
setorder(summary, age_group, exposure_bin, year)
fwrite(summary, file.path(out_dir, "acs_gse_exposure_summary_2005_2012.csv"))

event_outcomes <- c("owner_rate", "mean_rooms", "mean_rooms_owner", "mean_rooms_renter", "share_rooms_ge_6")
event_results <- rbindlist(lapply(event_outcomes, function(v) estimate_event_study(merged, v)), fill = TRUE)
fwrite(event_results, file.path(out_dir, "acs_gse_event_study_diagnostics.csv"))

display <- copy(summary[age_group == "25_34", .(
  year, exposure_bin, geographies, households,
  owner_rate = fmt_pct(owner_rate),
  mean_rooms = fmt_num(mean_rooms),
  mean_rooms_owner = fmt_num(mean_rooms_owner),
  mean_rooms_renter = fmt_num(mean_rooms_renter),
  share_rooms_ge_6 = fmt_pct(share_rooms_ge_6)
)])

reg_display <- copy(event_results[outcome %in% c("owner_rate", "mean_rooms", "share_rooms_ge_6")])
reg_display[, `:=`(
  estimate = fmt_num(estimate),
  std_error = fmt_num(std_error),
  p_value = fmt_num(p_value)
)]

summary_lines <- c(
  "# ACS GSE Exposure Outcomes",
  "",
  sprintf("Built at: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("ACS source directory: `%s`", acs_dir),
  sprintf("Years: `%s`", paste(years, collapse = ", ")),
  sprintf("Matched ACS geographies: `%d`", uniqueN(merged$geo_code)),
  sprintf("Treated median GSE log exposure cutoff: `%.4f`", treated_median),
  "",
  "## Main Use",
  "",
  "- This is the ACS housing-outcome panel for the GSE loan-limit shock.",
  "- The primary age group is household heads ages `25_34`.",
  "- The event-study diagnostics include metro and year fixed effects, weighted by ACS household weight.",
  "- Standard errors are clustered by metro when `fixest` is available; use the coefficients as screening diagnostics only.",
  "",
  "## Ages 25-34 Exposure-Bin Summary",
  "",
  paste(capture.output(print(display)), collapse = "\n"),
  "",
  "## Diagnostic FE Event-Study Coefficients",
  "",
  paste(capture.output(print(reg_display)), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "acs_gse_exposure_summary.md"))

message("Wrote:")
message("  ", file.path(out_dir, "acs_gse_exposure_panel_2005_2012.csv"))
message("  ", file.path(out_dir, "acs_gse_exposure_summary_2005_2012.csv"))
message("  ", file.path(out_dir, "acs_gse_event_study_diagnostics.csv"))
message("  ", file.path(out_dir, "acs_gse_exposure_summary.md"))
