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

fmt_pct <- function(x) {
  sprintf("%.2f%%", 100 * x)
}

script_dir <- get_script_dir()
out_dir <- file.path(script_dir, "output")
county_path <- file.path(out_dir, "hmda_county_first_stage.csv")
if (!file.exists(county_path)) {
  stop("Missing HMDA county first-stage file. Run build_hmda_first_stage_county.R first.")
}

county <- fread(
  county_path,
  colClasses = c(county_fips = "character", msa_code = "character")
)

positive_exposure <- unique(county[gse_log_change_2009_vs_2007_proxy > 0, gse_log_change_2009_vs_2007_proxy])
if (!length(positive_exposure)) {
  stop("No positive GSE exposure values found.")
}
treated_median <- median(positive_exposure, na.rm = TRUE)

county[, exposure_bin := fcase(
  gse_log_change_2009_vs_2007_proxy <= 0, "untreated",
  gse_log_change_2009_vs_2007_proxy < treated_median, "treated_below_median",
  default = "treated_above_median"
)]

bin_summary <- county[, .(
  counties = uniqueN(county_fips),
  originations = sum(originations, na.rm = TRUE),
  fha_originations = sum(fha_originations, na.rm = TRUE),
  conventional_gse_sold_originations = sum(conventional_gse_sold_originations, na.rm = TRUE),
  gse_new_band_originations = sum(gse_new_band_originations, na.rm = TRUE),
  gse_new_band_gse_sold_originations = sum(gse_new_band_gse_sold_originations, na.rm = TRUE),
  fha_new_band_originations = sum(fha_new_band_originations, na.rm = TRUE),
  fha_new_band_fha_originations = sum(fha_new_band_fha_originations, na.rm = TRUE),
  mean_gse_log_change_2009_vs_2007_proxy = mean(gse_log_change_2009_vs_2007_proxy, na.rm = TRUE),
  max_gse_log_change_2009_vs_2007_proxy = max(gse_log_change_2009_vs_2007_proxy, na.rm = TRUE)
), by = .(year, exposure_bin)]

bin_summary[, `:=`(
  fha_share = fha_originations / originations,
  conventional_gse_sold_share = conventional_gse_sold_originations / originations,
  gse_new_band_share = gse_new_band_originations / originations,
  gse_new_band_gse_sold_share = gse_new_band_gse_sold_originations / originations,
  fha_new_band_share = fha_new_band_originations / originations,
  fha_new_band_fha_share = fha_new_band_fha_originations / originations
)]
setorder(bin_summary, exposure_bin, year)
fwrite(bin_summary, file.path(out_dir, "hmda_first_stage_exposure_bins.csv"))

baseline <- bin_summary[year == 2007, .(
  exposure_bin,
  base_fha_share = fha_share,
  base_conventional_gse_sold_share = conventional_gse_sold_share,
  base_gse_new_band_gse_sold_share = gse_new_band_gse_sold_share,
  base_fha_new_band_fha_share = fha_new_band_fha_share
)]

bin_delta <- merge(bin_summary, baseline, by = "exposure_bin", all.x = TRUE)
bin_delta[, `:=`(
  delta_fha_share_vs_2007 = fha_share - base_fha_share,
  delta_conventional_gse_sold_share_vs_2007 = conventional_gse_sold_share - base_conventional_gse_sold_share,
  delta_gse_new_band_gse_sold_share_vs_2007 =
    gse_new_band_gse_sold_share - base_gse_new_band_gse_sold_share,
  delta_fha_new_band_fha_share_vs_2007 =
    fha_new_band_fha_share - base_fha_new_band_fha_share
)]
setorder(bin_delta, exposure_bin, year)
fwrite(bin_delta, file.path(out_dir, "hmda_first_stage_exposure_bin_deltas.csv"))

display <- copy(bin_summary)
display[, `:=`(
  fha_share = fmt_pct(fha_share),
  conventional_gse_sold_share = fmt_pct(conventional_gse_sold_share),
  gse_new_band_share = fmt_pct(gse_new_band_share),
  gse_new_band_gse_sold_share = fmt_pct(gse_new_band_gse_sold_share),
  fha_new_band_fha_share = fmt_pct(fha_new_band_fha_share)
)]
display <- display[, .(
  year, exposure_bin, counties, originations,
  fha_share,
  conventional_gse_sold_share,
  gse_new_band_share,
  gse_new_band_gse_sold_share,
  fha_new_band_fha_share
)]

high <- bin_summary[exposure_bin == "treated_above_median"]
high_2007 <- high[year == 2007, gse_new_band_gse_sold_share]
high_2008 <- high[year == 2008, gse_new_band_gse_sold_share]
high_2009 <- high[year == 2009, gse_new_band_gse_sold_share]

summary_lines <- c(
  "# HMDA First-Stage Diagnostic",
  "",
  sprintf("Built at: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Main Read",
  "",
  sprintf(
    "- In above-median treated counties, conventional loans sold to Fannie/Freddie in the newly eligible GSE band rose from `%s` of originations in 2007 to `%s` in 2008 and `%s` in 2009.",
    fmt_pct(high_2007), fmt_pct(high_2008), fmt_pct(high_2009)
  ),
  "- This is the first-stage object to track. Total originations fell sharply during the crisis, so levels alone are misleading.",
  "- The broad `gse_new_band_share` includes jumbo / portfolio loans that existed before eligibility; use the product-specific `gse_new_band_gse_sold_share` for the policy channel.",
  "- FHA 2009-vs-2008 band counts are small because the 2008 stimulus had already expanded FHA limits.",
  "",
  sprintf("Treated median GSE log exposure cutoff: `%.4f`", treated_median),
  "",
  "## Exposure-Bin Summary",
  "",
  paste(capture.output(print(display)), collapse = "\n")
)

writeLines(summary_lines, file.path(out_dir, "hmda_first_stage_summary.md"))

message("Wrote:")
message("  ", file.path(out_dir, "hmda_first_stage_exposure_bins.csv"))
message("  ", file.path(out_dir, "hmda_first_stage_exposure_bin_deltas.csv"))
message("  ", file.path(out_dir, "hmda_first_stage_summary.md"))
