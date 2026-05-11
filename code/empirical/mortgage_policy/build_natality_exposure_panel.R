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

fmt_rate <- function(x) {
  sprintf("%.2f", x)
}

script_dir <- get_script_dir()
out_dir <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

natality_path <- Sys.getenv(
  "NATALITY_PANEL",
  unset = "/Users/tommasodesanto/Desktop/Projects/Datasets/BirthsHybridMSAHpiPop_collapsed.dta"
)
if (!file.exists(natality_path)) {
  stop(sprintf("Missing natality panel: %s", natality_path))
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

message(sprintf("Reading natality panel: %s", natality_path))
nat <- as.data.table(read_dta(natality_path))

if ("divisioncode" %in% names(nat)) {
  nat[, geo_code := sprintf("%05d", as.integer(divisioncode))]
} else if ("cbsacode" %in% names(nat)) {
  nat[, geo_code := sprintf("%05d", as.integer(cbsacode))]
} else if ("GEOID" %in% names(nat)) {
  nat[, geo_code := sprintf("%05d", as.integer(GEOID))]
} else {
  stop("Natality panel has no recognized MSA/CBSA code column.")
}

needed <- c("geo_code", "year", "age_group", "recwt", "pop")
missing <- setdiff(needed, names(nat))
if (length(missing)) {
  stop(sprintf("Natality panel missing columns: %s", paste(missing, collapse = ", ")))
}

nat[, year := as.integer(year)]
nat[, age_group := as.integer(age_group)]
nat <- nat[year >= 2005 & year <= 2012 & !is.na(geo_code) & !is.na(age_group)]
nat <- nat[is.finite(recwt) & is.finite(pop) & pop > 0]

nat[, analysis_age_group := fcase(
  age_group == 2L, "20_24",
  age_group %in% c(3L, 4L), "25_34",
  age_group %in% c(5L, 6L), "35_44",
  default = NA_character_
)]
nat <- nat[!is.na(analysis_age_group)]

nat_agg <- nat[, .(
  births = sum(recwt, na.rm = TRUE),
  population = sum(pop, na.rm = TRUE),
  hpi = mean(HPI, na.rm = TRUE),
  fmr2br = if ("fmr2br_" %in% names(.SD)) mean(fmr2br_, na.rm = TRUE) else NA_real_
), by = .(geo_code, year, analysis_age_group)]
nat_agg[, birth_rate_per_1000 := births / population * 1000]

merged <- exposure[nat_agg, on = c(msa_code = "geo_code"), nomatch = 0]
setnames(merged, "msa_code", "geo_code")
setorder(merged, geo_code, analysis_age_group, year)
fwrite(merged, file.path(out_dir, "natality_gse_exposure_panel_2005_2012.csv"))

summary <- merged[, .(
  geographies = uniqueN(geo_code),
  births = sum(births, na.rm = TRUE),
  population = sum(population, na.rm = TRUE),
  birth_rate_per_1000 = sum(births, na.rm = TRUE) / sum(population, na.rm = TRUE) * 1000,
  mean_gse_log_exposure = mean(gse_log_change_2009_vs_2007_proxy_max, na.rm = TRUE)
), by = .(year, analysis_age_group, exposure_bin)]
setorder(summary, analysis_age_group, exposure_bin, year)
fwrite(summary, file.path(out_dir, "natality_gse_exposure_summary_2005_2012.csv"))

target_display <- summary[analysis_age_group == "25_34" & year %in% 2005:2012]
target_display[, birth_rate_per_1000_fmt := fmt_rate(birth_rate_per_1000)]
target_display <- target_display[, .(
  year, exposure_bin, geographies, births, population,
  birth_rate_per_1000 = birth_rate_per_1000_fmt
)]

matched_geos <- uniqueN(merged$geo_code)
source_geos <- uniqueN(nat_agg$geo_code)
summary_lines <- c(
  "# Natality Exposure Merge",
  "",
  sprintf("Built at: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  sprintf("Natality source: `%s`", natality_path),
  sprintf("Matched geographies: `%d` of `%d` natality geographies in 2005-2012 target age groups.", matched_geos, source_geos),
  sprintf("Treated median GSE log exposure cutoff: `%.4f`", treated_median),
  "",
  "## Main Use",
  "",
  "- This is a merged outcome panel, not a causal estimate.",
  "- Start with ages `25_34`, where first-time-buyer / family-formation exposure is most plausible.",
  "- Use ages `35_44` as a weaker placebo / comparison group, not as a definitive falsification test.",
  "- Add MSA/division and year fixed effects before making any fertility claim.",
  "",
  "## Birth Rates Per 1,000, Ages 25-34",
  "",
  paste(capture.output(print(target_display)), collapse = "\n")
)

writeLines(summary_lines, file.path(out_dir, "natality_gse_exposure_summary.md"))

message("Wrote:")
message("  ", file.path(out_dir, "natality_gse_exposure_panel_2005_2012.csv"))
message("  ", file.path(out_dir, "natality_gse_exposure_summary_2005_2012.csv"))
message("  ", file.path(out_dir, "natality_gse_exposure_summary.md"))
