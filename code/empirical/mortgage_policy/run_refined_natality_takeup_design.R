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

fmt_num <- function(x) sprintf("%.3f", x)

load_unemp_controls <- function(path) {
  if (!file.exists(path)) {
    warning(sprintf("Missing unemployment file: %s", path))
    return(data.table(geo_code = character(), year = integer(), unemp_rate = numeric()))
  }
  x <- as.data.table(read_dta(path))
  x[, geo_code := sprintf("%05d", as.integer(GEOID))]
  x[, year := as.integer(year)]
  x <- x[year >= 2005 & year <= 2012 & is.finite(annual_value)]
  x[, .(unemp_rate = mean(annual_value, na.rm = TRUE)), by = .(geo_code, year)]
}

event_reg <- function(dt, age_group_value, spec_name, treatment_var) {
  work <- dt[
    analysis_age_group == age_group_value &
      year >= 2005 & year <= 2012 &
      is.finite(birth_rate_per_1000) &
      is.finite(get(treatment_var)) &
      is.finite(population) &
      population > 0 &
      is.finite(log_hpi_2007_norm) &
      is.finite(unemp_rate)
  ]
  if (nrow(work) < 30) {
    return(data.table())
  }
  event_years <- setdiff(sort(unique(work$year)), 2007L)
  for (yy in event_years) {
    work[, paste0("treat_x_", yy) := get(treatment_var) * as.integer(year == yy)]
  }
  rhs <- paste(c(paste0("treat_x_", event_years), "log_hpi_2007_norm", "unemp_rate"), collapse = " + ")
  if (requireNamespace("fixest", quietly = TRUE)) {
    fit <- fixest::feols(
      as.formula(sprintf("birth_rate_per_1000 ~ %s | geo_code + year", rhs)),
      data = work,
      weights = ~population,
      cluster = ~geo_code,
      warn = FALSE,
      notes = FALSE
    )
    coefs <- as.data.table(fixest::coeftable(fit), keep.rownames = "term")
    note_text <- "Metro and year FE; population weights; SE clustered by metro; controls included."
  } else {
    fit <- lm(
      as.formula(sprintf("birth_rate_per_1000 ~ %s + factor(geo_code) + factor(year)", rhs)),
      data = work,
      weights = population
    )
    coefs <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
    note_text <- "Metro and year FE; population weights; SE not clustered because fixest is unavailable."
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
    analysis_age_group = age_group_value,
    event_year = as.integer(sub("treat_x_", "", term)),
    baseline_year = 2007L,
    observations = nrow(work),
    metros = uniqueN(work$geo_code),
    note = note_text
  )]
  coefs[, .(
    spec, treatment, analysis_age_group, baseline_year, event_year, estimate,
    std_error, t_value, p_value, observations, metros, note
  )]
}

script_dir <- get_script_dir()
out_dir <- file.path(script_dir, "output")

natality_path <- file.path(out_dir, "natality_gse_exposure_panel_2005_2012.csv")
design_path <- file.path(out_dir, "hmda_refined_takeup_design.csv")
if (!file.exists(natality_path)) {
  stop("Missing natality exposure panel. Run build_natality_exposure_panel.R first.")
}
if (!file.exists(design_path)) {
  stop("Missing refined HMDA takeup design. Run run_refined_acs_takeup_design.R first.")
}

unemp_path <- Sys.getenv(
  "UNEMP_PANEL",
  unset = "/Users/tommasodesanto/Desktop/Projects/Datasets/unemployment_rates/msa_unemployment_clean_nodups.dta"
)

nat <- fread(natality_path, colClasses = c(geo_code = "character"))
design <- fread(design_path, colClasses = c(geo_code = "character"))
unemp <- load_unemp_controls(unemp_path)

nat <- merge(nat, design[, .(
  geo_code,
  positive_gse_exposure,
  high_takeup,
  takeup_delta_post_vs_2007,
  max_gse_log_change_2009_vs_2007_proxy
)], by = "geo_code", all.x = TRUE)
nat <- merge(nat, unemp, by = c("geo_code", "year"), all.x = TRUE)
nat[, hpi_2007 := hpi[year == 2007][1], by = geo_code]
nat[, log_hpi_2007_norm := fifelse(
  is.finite(hpi) & is.finite(hpi_2007) & hpi > 0 & hpi_2007 > 0,
  log(hpi / hpi_2007),
  NA_real_
)]
nat[, high_takeup := fifelse(is.na(high_takeup), 0L, high_takeup)]
nat[, takeup_delta_post_vs_2007 := fifelse(is.finite(takeup_delta_post_vs_2007), takeup_delta_post_vs_2007, 0)]
fwrite(nat, file.path(out_dir, "natality_refined_takeup_panel_2005_2012.csv"))

positive <- nat[positive_gse_exposure == TRUE]
high_low <- positive[high_takeup == 1L | takeup_delta_post_vs_2007 < unique(design[high_takeup == 1L, min(takeup_delta_post_vs_2007)])]
age_groups <- c("20_24", "25_34", "35_44")

results <- rbindlist(list(
  rbindlist(lapply(age_groups, event_reg,
    dt = high_low,
    spec_name = "positive_exposure_high_takeup_binary",
    treatment_var = "high_takeup"
  ), fill = TRUE),
  rbindlist(lapply(age_groups, event_reg,
    dt = positive,
    spec_name = "positive_exposure_continuous_takeup_delta",
    treatment_var = "takeup_delta_post_vs_2007"
  ), fill = TRUE),
  rbindlist(lapply(age_groups, event_reg,
    dt = nat,
    spec_name = "all_metros_continuous_takeup_delta",
    treatment_var = "takeup_delta_post_vs_2007"
  ), fill = TRUE)
), fill = TRUE)
fwrite(results, file.path(out_dir, "natality_refined_takeup_event_study.csv"))

summary_bins <- high_low[analysis_age_group == "25_34", .(
  geographies = uniqueN(geo_code),
  births = sum(births, na.rm = TRUE),
  population = sum(population, na.rm = TRUE),
  birth_rate_per_1000 = sum(births, na.rm = TRUE) / sum(population, na.rm = TRUE) * 1000,
  log_hpi_2007_norm = weighted.mean(log_hpi_2007_norm, population, na.rm = TRUE),
  unemp_rate = weighted.mean(unemp_rate, population, na.rm = TRUE)
), by = .(year, high_takeup)]
setorder(summary_bins, high_takeup, year)
fwrite(summary_bins, file.path(out_dir, "natality_refined_takeup_summary_bins.csv"))

display_bins <- copy(summary_bins)
display_bins[, `:=`(
  high_takeup = fifelse(high_takeup == 1L, "high_takeup", "other_positive_exposure"),
  birth_rate_per_1000 = fmt_num(birth_rate_per_1000),
  log_hpi_2007_norm = fmt_num(log_hpi_2007_norm),
  unemp_rate = fmt_num(unemp_rate)
)]

display_results <- copy(results[
  spec %in% c("positive_exposure_high_takeup_binary", "positive_exposure_continuous_takeup_delta") &
    analysis_age_group == "25_34"
])
display_results[, `:=`(
  estimate = fmt_num(estimate),
  std_error = fmt_num(std_error),
  p_value = fmt_num(p_value)
)]
setorder(display_results, spec, event_year)

summary_lines <- c(
  "# Refined Natality GSE Takeup Design",
  "",
  sprintf("Built at: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Read",
  "",
  "- This applies the same high-takeup GSE screen used for ACS outcomes to the existing natality scaffold.",
  "- Treat this as a diagnostic only: the natality panel itself still needs audit before any fertility claim.",
  "- The main group shown below is ages `25_34`.",
  "",
  "## Ages 25-34 Summary Within Positive-Exposure Metros",
  "",
  paste(capture.output(print(display_bins)), collapse = "\n"),
  "",
  "## Event-Study Diagnostics, Ages 25-34",
  "",
  paste(capture.output(print(display_results)), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "natality_refined_takeup_summary.md"))

message("Wrote:")
message("  ", file.path(out_dir, "natality_refined_takeup_panel_2005_2012.csv"))
message("  ", file.path(out_dir, "natality_refined_takeup_summary_bins.csv"))
message("  ", file.path(out_dir, "natality_refined_takeup_event_study.csv"))
message("  ", file.path(out_dir, "natality_refined_takeup_summary.md"))
