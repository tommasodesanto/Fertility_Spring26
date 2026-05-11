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

fmt_pct <- function(x) sprintf("%.2f%%", 100 * x)
fmt_num <- function(x) sprintf("%.3f", x)

safe_wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

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
  out <- x[, .(hpi = mean(HPI, na.rm = TRUE)), by = .(geo_code, year)]
  setorder(out, geo_code, year)
  out
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
  out <- x[, .(unemp_rate = mean(annual_value, na.rm = TRUE)), by = .(geo_code, year)]
  setorder(out, geo_code, year)
  out
}

event_reg <- function(dt, outcome, spec_name, treatment_var, use_controls = TRUE) {
  work <- dt[
    age_group == "25_34" &
      year >= 2005 & year <= 2012 &
      is.finite(get(outcome)) &
      is.finite(get(treatment_var)) &
      is.finite(hh_weight) &
      hh_weight > 0
  ]
  if (use_controls) {
    work <- work[is.finite(log_hpi_2007_norm) & is.finite(unemp_rate)]
  }
  if (nrow(work) < 30) {
    return(data.table())
  }
  event_years <- setdiff(sort(unique(work$year)), 2007L)
  for (yy in event_years) {
    work[, paste0("treat_x_", yy) := get(treatment_var) * as.integer(year == yy)]
  }
  rhs <- paste(paste0("treat_x_", event_years), collapse = " + ")
  if (use_controls) {
    rhs <- paste(rhs, "+ log_hpi_2007_norm + unemp_rate")
  }

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
    note_text <- "Metro and year FE; ACS weights; SE clustered by metro; controls included where specified."
  } else {
    lm_rhs <- paste(rhs, "+ factor(geo_code) + factor(year)")
    fit <- lm(
      as.formula(sprintf("%s ~ %s", outcome, lm_rhs)),
      data = work,
      weights = hh_weight
    )
    coefs <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
    note_text <- "Metro and year FE; ACS weights; SE not clustered because fixest is unavailable."
  }
  coefs <- coefs[startsWith(term, "treat_x_")]
  if (!nrow(coefs)) {
    return(data.table())
  }
  setnames(coefs, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    c("estimate", "std_error", "t_value", "p_value")
  )
  coefs[, `:=`(
    spec = spec_name,
    treatment = treatment_var,
    outcome = outcome,
    event_year = as.integer(sub("treat_x_", "", term)),
    baseline_year = 2007L,
    observations = nrow(work),
    metros = uniqueN(work$geo_code),
    controls = use_controls,
    note = note_text
  )]
  coefs[, .(
    spec, treatment, outcome, baseline_year, event_year, estimate, std_error,
    t_value, p_value, observations, metros, controls, note
  )]
}

script_dir <- get_script_dir()
out_dir <- file.path(script_dir, "output")

acs_path <- file.path(out_dir, "acs_gse_exposure_panel_2005_2012.csv")
hmda_path <- file.path(out_dir, "hmda_cbsa_first_stage.csv")
if (!file.exists(acs_path)) {
  stop("Missing ACS exposure panel. Run build_acs_gse_exposure_panel.R first.")
}
if (!file.exists(hmda_path)) {
  stop("Missing HMDA CBSA first-stage file. Run build_hmda_first_stage_county.R first.")
}

hpi_path <- Sys.getenv(
  "HPI_PANEL",
  unset = "/Users/tommasodesanto/Desktop/Projects/Datasets/BirthsHybridMSAHpiPop_collapsed.dta"
)
unemp_path <- Sys.getenv(
  "UNEMP_PANEL",
  unset = "/Users/tommasodesanto/Desktop/Projects/Datasets/unemployment_rates/msa_unemployment_clean_nodups.dta"
)

acs <- fread(acs_path, colClasses = c(geo_code = "character"))
hmda <- fread(hmda_path, colClasses = c(msa_code = "character"))

hmda[, gse_new_band_gse_sold_share := gse_new_band_gse_sold_originations / originations]
hmda[, gse_new_band_share := gse_new_band_originations / originations]
hmda <- hmda[, .(
  msa_name = paste(sort(unique(msa_name)), collapse = " / "),
  counties = sum(counties, na.rm = TRUE),
  originations = sum(originations, na.rm = TRUE),
  gse_new_band_gse_sold_originations = sum(gse_new_band_gse_sold_originations, na.rm = TRUE),
  gse_new_band_originations = sum(gse_new_band_originations, na.rm = TRUE),
  max_gse_log_change_2009_vs_2007_proxy = max(max_gse_log_change_2009_vs_2007_proxy, na.rm = TRUE)
), by = .(msa_code, year)]
hmda[, gse_new_band_gse_sold_share := gse_new_band_gse_sold_originations / originations]
hmda[, gse_new_band_share := gse_new_band_originations / originations]
hmda_wide <- dcast(
  hmda,
  msa_code + msa_name + max_gse_log_change_2009_vs_2007_proxy ~ year,
  value.var = c("gse_new_band_gse_sold_share", "gse_new_band_share", "originations")
)
hmda_wide[, `:=`(
  takeup_delta_2008_vs_2007 =
    gse_new_band_gse_sold_share_2008 - gse_new_band_gse_sold_share_2007,
  takeup_delta_2009_vs_2007 =
    gse_new_band_gse_sold_share_2009 - gse_new_band_gse_sold_share_2007,
  takeup_delta_post_vs_2007 =
    rowMeans(cbind(gse_new_band_gse_sold_share_2008, gse_new_band_gse_sold_share_2009), na.rm = TRUE) -
      gse_new_band_gse_sold_share_2007,
  positive_gse_exposure = max_gse_log_change_2009_vs_2007_proxy > 0
)]

treated_deltas <- hmda_wide[
  positive_gse_exposure == TRUE & is.finite(takeup_delta_post_vs_2007),
  takeup_delta_post_vs_2007
]
if (!length(treated_deltas)) {
  stop("No treated HMDA takeup deltas found.")
}
takeup_cutoff <- as.numeric(quantile(treated_deltas, 0.75, na.rm = TRUE))
hmda_wide[, high_takeup := as.integer(
  positive_gse_exposure == TRUE &
    is.finite(takeup_delta_post_vs_2007) &
    takeup_delta_post_vs_2007 >= takeup_cutoff
)]

design <- hmda_wide[, .(
  geo_code = msa_code,
  msa_name,
  positive_gse_exposure,
  max_gse_log_change_2009_vs_2007_proxy,
  gse_new_band_gse_sold_share_2007,
  gse_new_band_gse_sold_share_2008,
  gse_new_band_gse_sold_share_2009,
  takeup_delta_2008_vs_2007,
  takeup_delta_2009_vs_2007,
  takeup_delta_post_vs_2007,
  high_takeup
)]
fwrite(design, file.path(out_dir, "hmda_refined_takeup_design.csv"))

hpi <- load_hpi_controls(hpi_path)
unemp <- load_unemp_controls(unemp_path)

panel <- merge(acs, design, by = "geo_code", all.x = TRUE)
panel <- merge(panel, hpi, by = c("geo_code", "year"), all.x = TRUE)
panel <- merge(panel, unemp, by = c("geo_code", "year"), all.x = TRUE)

panel[, hpi_2007 := hpi[year == 2007][1], by = geo_code]
panel[, log_hpi_2007_norm := fifelse(
  is.finite(hpi) & is.finite(hpi_2007) & hpi > 0 & hpi_2007 > 0,
  log(hpi / hpi_2007),
  NA_real_
)]
panel[, high_takeup := fifelse(is.na(high_takeup), 0L, high_takeup)]
panel[, takeup_delta_post_vs_2007 := fifelse(
  is.finite(takeup_delta_post_vs_2007),
  takeup_delta_post_vs_2007,
  0
)]

fwrite(panel, file.path(out_dir, "acs_refined_takeup_panel_2005_2012.csv"))

sample_positive <- panel[positive_gse_exposure == TRUE]
sample_high_low <- sample_positive[
  is.finite(takeup_delta_post_vs_2007) &
    (high_takeup == 1L | takeup_delta_post_vs_2007 < takeup_cutoff)
]

outcomes <- c("owner_rate", "mean_rooms", "mean_rooms_owner", "mean_rooms_renter", "share_rooms_ge_6")
results <- rbindlist(list(
  rbindlist(lapply(outcomes, event_reg,
    dt = sample_high_low,
    spec_name = "positive_exposure_high_takeup_binary",
    treatment_var = "high_takeup",
    use_controls = TRUE
  ), fill = TRUE),
  rbindlist(lapply(outcomes, event_reg,
    dt = sample_positive,
    spec_name = "positive_exposure_continuous_takeup_delta",
    treatment_var = "takeup_delta_post_vs_2007",
    use_controls = TRUE
  ), fill = TRUE),
  rbindlist(lapply(outcomes, event_reg,
    dt = panel,
    spec_name = "all_metros_continuous_takeup_delta",
    treatment_var = "takeup_delta_post_vs_2007",
    use_controls = TRUE
  ), fill = TRUE)
), fill = TRUE)
fwrite(results, file.path(out_dir, "acs_refined_takeup_event_study.csv"))

summary_bins <- sample_high_low[age_group == "25_34", .(
  geographies = uniqueN(geo_code),
  households = sum(households, na.rm = TRUE),
  hh_weight = sum(hh_weight, na.rm = TRUE),
  owner_rate = safe_wmean(owner_rate, hh_weight),
  mean_rooms = safe_wmean(mean_rooms, hh_weight),
  mean_rooms_owner = safe_wmean(mean_rooms_owner, hh_weight),
  mean_rooms_renter = safe_wmean(mean_rooms_renter, hh_weight),
  share_rooms_ge_6 = safe_wmean(share_rooms_ge_6, hh_weight),
  log_hpi_2007_norm = safe_wmean(log_hpi_2007_norm, hh_weight),
  unemp_rate = safe_wmean(unemp_rate, hh_weight)
), by = .(year, high_takeup)]
setorder(summary_bins, high_takeup, year)
fwrite(summary_bins, file.path(out_dir, "acs_refined_takeup_summary_bins.csv"))

display_bins <- copy(summary_bins)
display_bins[, `:=`(
  high_takeup = fifelse(high_takeup == 1L, "high_takeup", "other_positive_exposure"),
  owner_rate = fmt_pct(owner_rate),
  mean_rooms = fmt_num(mean_rooms),
  mean_rooms_owner = fmt_num(mean_rooms_owner),
  mean_rooms_renter = fmt_num(mean_rooms_renter),
  share_rooms_ge_6 = fmt_pct(share_rooms_ge_6),
  log_hpi_2007_norm = fmt_num(log_hpi_2007_norm),
  unemp_rate = fmt_num(unemp_rate)
)]

display_results <- copy(results[
  spec %in% c("positive_exposure_high_takeup_binary", "positive_exposure_continuous_takeup_delta") &
    outcome %in% c("owner_rate", "mean_rooms", "share_rooms_ge_6")
])
display_results[, `:=`(
  estimate = fmt_num(estimate),
  std_error = fmt_num(std_error),
  p_value = fmt_num(p_value)
)]
setorder(display_results, spec, outcome, event_year)

high_n <- design[high_takeup == 1L, .N]
positive_n <- design[positive_gse_exposure == TRUE, .N]
summary_lines <- c(
  "# Refined ACS GSE Takeup Design",
  "",
  sprintf("Built at: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  sprintf("- positive-exposure metros in HMDA file: `%d`", positive_n),
  sprintf("- high-takeup cutoff: post-2008/2009 average GSE-sold new-band share minus 2007 share >= `%.4f`", takeup_cutoff),
  sprintf("- high-takeup metros: `%d`", high_n),
  "- controls: log HPI relative to 2007 and MSA unemployment rate",
  "- estimators: metro and year FE, ACS household weights, metro-clustered SEs via `fixest`",
  "",
  "## Read",
  "",
  "- This is still a screening design, but it is narrower than the broad exposure screen.",
  "- If the mortgage-access story is in the ACS housing margins, it should appear in the high-takeup binary or continuous-takeup specifications.",
  "- Do not treat the high-takeup definition as exogenous; it is an empirical first-stage filter.",
  "",
  "## Ages 25-34 Summary Within Positive-Exposure Metros",
  "",
  paste(capture.output(print(display_bins)), collapse = "\n"),
  "",
  "## Event-Study Diagnostics",
  "",
  paste(capture.output(print(display_results)), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "acs_refined_takeup_summary.md"))

message("Wrote:")
message("  ", file.path(out_dir, "hmda_refined_takeup_design.csv"))
message("  ", file.path(out_dir, "acs_refined_takeup_panel_2005_2012.csv"))
message("  ", file.path(out_dir, "acs_refined_takeup_summary_bins.csv"))
message("  ", file.path(out_dir, "acs_refined_takeup_event_study.csv"))
message("  ", file.path(out_dir, "acs_refined_takeup_summary.md"))
