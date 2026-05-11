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

to_num <- function(x) {
  suppressWarnings(as.numeric(trimws(x)))
}

trim_substr <- function(x, first, last) {
  trimws(substr(x, first, last))
}

state_fips_lookup <- function() {
  data.table(
    state = c(
      "AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "DC", "FL",
      "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME",
      "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH",
      "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI",
      "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI",
      "WY", "AS", "GU", "MP", "PR", "VI"
    ),
    state_fips = c(
      "01", "02", "04", "05", "06", "08", "09", "10", "11", "12",
      "13", "15", "16", "17", "18", "19", "20", "21", "22", "23",
      "24", "25", "26", "27", "28", "29", "30", "31", "32", "33",
      "34", "35", "36", "37", "38", "39", "40", "41", "42", "44",
      "45", "46", "47", "48", "49", "50", "51", "53", "54", "55",
      "56", "60", "66", "69", "72", "78"
    )
  )
}

parse_hud_limit_file <- function(path, year, product) {
  if (!file.exists(path)) {
    stop(sprintf("Missing HUD loan-limit file: %s", path))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nchar(lines, type = "chars", allowNA = FALSE, keepNA = FALSE) >= 163]
  dt <- data.table(
    msa_code = trim_substr(lines, 1, 5),
    md_code = trim_substr(lines, 6, 10),
    msa_name = trim_substr(lines, 11, 60),
    soa_code = trim_substr(lines, 61, 65),
    limit_type = trim_substr(lines, 66, 66),
    median_price = to_num(substr(lines, 67, 73)),
    limit_1unit = to_num(substr(lines, 74, 80)),
    limit_2unit = to_num(substr(lines, 81, 87)),
    limit_3unit = to_num(substr(lines, 88, 94)),
    limit_4unit = to_num(substr(lines, 95, 101)),
    state = trim_substr(lines, 102, 103),
    county_code = trim_substr(lines, 104, 106),
    state_name = trim_substr(lines, 107, 132),
    county_name = trim_substr(lines, 133, 147),
    county_tx_date = trim_substr(lines, 148, 155),
    limit_tx_date = trim_substr(lines, 156, 163),
    median_price_determining_limit = to_num(substr(lines, 164, 170)),
    median_price_year = to_num(substr(lines, 171, 175))
  )

  state_map <- state_fips_lookup()
  dt <- state_map[dt, on = "state"]
  dt <- dt[!is.na(state_fips) & county_code != "" & is.finite(limit_1unit)]
  dt[, county_fips := paste0(state_fips, county_code)]
  dt[, `:=`(year = as.integer(year), product = product)]
  setcolorder(dt, c(
    "product", "year", "county_fips", "state", "county_code", "county_name",
    "state_name", "msa_code", "md_code", "msa_name", "soa_code", "limit_type",
    "median_price", "limit_1unit", "limit_2unit", "limit_3unit", "limit_4unit",
    "county_tx_date", "limit_tx_date", "median_price_determining_limit",
    "median_price_year"
  ))
  unique(dt, by = c("product", "year", "county_fips"))
}

pick_limit <- function(dt, product_name, year_value, value_name) {
  out <- dt[product == product_name & year == year_value, .(
    county_fips,
    value = limit_1unit,
    limit_type = limit_type,
    median_price = median_price
  )]
  setnames(out, c("value", "limit_type", "median_price"), c(
    value_name,
    paste0(value_name, "_type"),
    paste0(value_name, "_median_price")
  ))
  out
}

safe_log_change <- function(new, old) {
  fifelse(is.finite(new) & is.finite(old) & new > 0 & old > 0, log(new) - log(old), NA_real_)
}

script_dir <- get_script_dir()
raw_dir <- file.path(script_dir, "raw")
out_dir <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sources <- rbindlist(list(
  parse_hud_limit_file(file.path(raw_dir, "cy2008-forward-limits.txt"), 2008, "fha_forward"),
  parse_hud_limit_file(file.path(raw_dir, "cy2009-forward-limits.txt"), 2009, "fha_forward"),
  parse_hud_limit_file(file.path(raw_dir, "cy2008-gse-limits.txt"), 2008, "gse"),
  parse_hud_limit_file(file.path(raw_dir, "cy2009-gse-limits.txt"), 2009, "gse")
), fill = TRUE)

fwrite(sources, file.path(out_dir, "hud_loan_limits_county_long_2008_2009.csv"))

base <- unique(sources[year == 2009, .(
  county_fips, state, county_code, county_name, state_name, msa_code, md_code, msa_name
)], by = "county_fips")

pieces <- list(
  base,
  pick_limit(sources, "fha_forward", 2008, "fha_2008_limit"),
  pick_limit(sources, "fha_forward", 2009, "fha_2009_limit"),
  pick_limit(sources, "gse", 2008, "gse_2008_limit"),
  pick_limit(sources, "gse", 2009, "gse_2009_limit")
)

panel <- Reduce(function(x, y) merge(x, y, by = "county_fips", all = TRUE), pieces)

statutory_high_cost_states <- c("AK", "HI", "GU", "VI")
panel[, gse_2007_proxy_limit := fifelse(state %in% statutory_high_cost_states, 625500, 417000)]
panel[, fha_2009_minus_2008 := fha_2009_limit - fha_2008_limit]
panel[, gse_2009_minus_2008 := gse_2009_limit - gse_2008_limit]
panel[, gse_2009_minus_2007_proxy := gse_2009_limit - gse_2007_proxy_limit]
panel[, fha_log_change_2009_vs_2008 := safe_log_change(fha_2009_limit, fha_2008_limit)]
panel[, gse_log_change_2009_vs_2008 := safe_log_change(gse_2009_limit, gse_2008_limit)]
panel[, gse_log_change_2009_vs_2007_proxy := safe_log_change(gse_2009_limit, gse_2007_proxy_limit)]
panel[, gse_newly_eligible_band_low := pmin(gse_2007_proxy_limit, gse_2009_limit, na.rm = TRUE)]
panel[, gse_newly_eligible_band_high := pmax(gse_2007_proxy_limit, gse_2009_limit, na.rm = TRUE)]
panel[, fha_change_sign_2009_vs_2008 := fifelse(fha_2009_minus_2008 > 0, "increase",
  fifelse(fha_2009_minus_2008 < 0, "decrease", "flat")
)]
panel[, gse_change_sign_2009_vs_2007_proxy := fifelse(gse_2009_minus_2007_proxy > 0, "increase",
  fifelse(gse_2009_minus_2007_proxy < 0, "decrease", "flat")
)]

setorder(panel, state, county_name)
fwrite(panel, file.path(out_dir, "loan_limit_exposure_county_2008_2009.csv"))

cbsa <- panel[, .(
  counties = .N,
  states = paste(sort(unique(state)), collapse = " "),
  fha_2008_limit_max = max(fha_2008_limit, na.rm = TRUE),
  fha_2009_limit_max = max(fha_2009_limit, na.rm = TRUE),
  gse_2007_proxy_limit_max = max(gse_2007_proxy_limit, na.rm = TRUE),
  gse_2008_limit_max = max(gse_2008_limit, na.rm = TRUE),
  gse_2009_limit_max = max(gse_2009_limit, na.rm = TRUE),
  fha_log_change_2009_vs_2008_max = max(fha_log_change_2009_vs_2008, na.rm = TRUE),
  gse_log_change_2009_vs_2007_proxy_max = max(gse_log_change_2009_vs_2007_proxy, na.rm = TRUE),
  share_counties_fha_increase_2009_vs_2008 = mean(fha_2009_minus_2008 > 0, na.rm = TRUE),
  share_counties_gse_increase_2009_vs_2007_proxy = mean(gse_2009_minus_2007_proxy > 0, na.rm = TRUE)
), by = .(msa_code, msa_name)]
setorder(cbsa, -gse_log_change_2009_vs_2007_proxy_max, msa_name)
fwrite(cbsa, file.path(out_dir, "loan_limit_exposure_cbsa_2008_2009.csv"))

top_counties <- panel[order(-gse_log_change_2009_vs_2007_proxy)][
  1:min(.N, 25),
  .(
    county_fips, state, county_name, msa_name,
    gse_2007_proxy_limit, gse_2009_limit,
    gse_log_change_2009_vs_2007_proxy,
    fha_2008_limit, fha_2009_limit,
    fha_log_change_2009_vs_2008
  )
]

old_width <- getOption("width")
options(width = 200)
summary_lines <- c(
  "# Loan-Limit Exposure Build",
  "",
  sprintf("Built at: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## County Counts",
  "",
  sprintf("- counties in panel: `%s`", format(nrow(panel), big.mark = ",")),
  sprintf("- counties with GSE 2009 limit above 2007 proxy: `%s`",
    format(panel[, sum(gse_2009_minus_2007_proxy > 0, na.rm = TRUE)], big.mark = ",")
  ),
  sprintf("- counties with FHA 2009 limit above FHA 2008 limit: `%s`",
    format(panel[, sum(fha_2009_minus_2008 > 0, na.rm = TRUE)], big.mark = ",")
  ),
  "",
  "## Interpretation",
  "",
  "- `gse_2007_proxy_limit` is set to `$417,000` outside statutory high-cost areas and `$625,500` in AK/HI/GU/VI.",
  "- The GSE 2009-vs-2007 measure is the cleanest available public proxy for exposure to the crisis-era conforming-limit expansion.",
  "- FHA 2009-vs-2008 is a narrower comparison because the 2008 Economic Stimulus Act had already expanded FHA limits.",
  "- Use HMDA to validate that these exposure measures actually predict loan originations in newly eligible bands before running fertility outcomes.",
  "",
  "## Top Counties By GSE 2009-vs-2007 Proxy Log Change",
  "",
  paste(capture.output(print(top_counties)), collapse = "\n")
)
options(width = old_width)

writeLines(summary_lines, file.path(out_dir, "loan_limit_exposure_summary.md"))

message("Wrote:")
message("  ", file.path(out_dir, "hud_loan_limits_county_long_2008_2009.csv"))
message("  ", file.path(out_dir, "loan_limit_exposure_county_2008_2009.csv"))
message("  ", file.path(out_dir, "loan_limit_exposure_cbsa_2008_2009.csv"))
message("  ", file.path(out_dir, "loan_limit_exposure_summary.md"))
