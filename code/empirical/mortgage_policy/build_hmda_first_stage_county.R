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

clean_names <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

as_int <- function(x) {
  suppressWarnings(as.integer(x))
}

zip_csv_cmd <- function(zip_path) {
  listing <- unzip(zip_path, list = TRUE)
  csv_name <- listing$Name[grepl("\\.csv$", listing$Name, ignore.case = TRUE)][1]
  if (is.na(csv_name)) {
    stop(sprintf("No CSV found inside %s", zip_path))
  }
  sprintf("unzip -p %s %s", shQuote(zip_path), shQuote(csv_name))
}

read_hmda_year <- function(zip_path, exposure, year_value) {
  if (!file.exists(zip_path)) {
    stop(sprintf(
      "Missing HMDA zip for %s: %s\nRun DOWNLOAD_HMDA=1 Rscript download_mortgage_policy_sources.R first.",
      year_value, zip_path
    ))
  }
  cmd <- zip_csv_cmd(zip_path)
  header <- names(fread(cmd = cmd, nrows = 0, showProgress = FALSE))
  clean_header <- clean_names(header)
  select_clean <- c(
    "as_of_year", "loan_type", "property_type", "loan_purpose", "occupancy",
    "owner_occupancy", "loan_amount_000s", "action_type", "action_taken",
    "msa_md", "msamd", "state_code", "county_code",
    "applicant_income_000s", "purchaser_type", "lien_status"
  )
  select_cols <- header[clean_header %in% select_clean]
  if (!length(select_cols)) {
    stop(sprintf("Could not match HMDA columns in %s", zip_path))
  }
  dt <- fread(cmd = cmd, select = select_cols, showProgress = TRUE)
  setnames(dt, clean_names(names(dt)))
  if ("owner_occupancy" %in% names(dt) && !"occupancy" %in% names(dt)) {
    setnames(dt, "owner_occupancy", "occupancy")
  }
  if ("action_taken" %in% names(dt) && !"action_type" %in% names(dt)) {
    setnames(dt, "action_taken", "action_type")
  }
  if ("msamd" %in% names(dt) && !"msa_md" %in% names(dt)) {
    setnames(dt, "msamd", "msa_md")
  }

  required <- c("loan_type", "loan_purpose", "action_type", "state_code", "county_code", "loan_amount_000s")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(sprintf("Missing required HMDA columns: %s", paste(missing, collapse = ", ")))
  }

  for (v in intersect(names(dt), select_clean)) {
    dt[, (v) := as_int(get(v))]
  }
  if (!"as_of_year" %in% names(dt)) {
    dt[, as_of_year := as.integer(year_value)]
  }
  if (!"property_type" %in% names(dt)) {
    dt[, property_type := 1L]
  }
  if (!"occupancy" %in% names(dt)) {
    dt[, occupancy := 1L]
  }
  if (!"lien_status" %in% names(dt)) {
    dt[, lien_status := 1L]
  }
  if (!"purchaser_type" %in% names(dt)) {
    dt[, purchaser_type := NA_integer_]
  }
  if (!"applicant_income_000s" %in% names(dt)) {
    dt[, applicant_income_000s := NA_integer_]
  }

  dt <- dt[
    action_type == 1L &
      loan_purpose == 1L &
      property_type == 1L &
      occupancy == 1L &
      lien_status == 1L &
      !is.na(state_code) &
      !is.na(county_code) &
      !is.na(loan_amount_000s)
  ]
  dt[, county_fips := sprintf("%02d%03d", state_code, county_code)]
  dt[, loan_amount := loan_amount_000s * 1000]
  dt <- exposure[dt, on = "county_fips", nomatch = 0]

  dt[, is_fha := loan_type == 2L]
  dt[, is_conventional_gse_sold := loan_type == 1L & purchaser_type %in% c(1L, 3L)]
  dt[, is_gse_newly_eligible_2009_vs_2007 := loan_amount > gse_2007_proxy_limit & loan_amount <= gse_2009_limit]
  dt[, is_fha_2009_vs_2008_band := loan_amount > fha_2008_limit & loan_amount <= fha_2009_limit]
  dt[, is_near_gse_old_limit := loan_amount > 0.90 * gse_2007_proxy_limit & loan_amount <= gse_2007_proxy_limit]

  dt[, .(
    originations = .N,
    fha_originations = sum(is_fha, na.rm = TRUE),
    conventional_gse_sold_originations = sum(is_conventional_gse_sold, na.rm = TRUE),
    gse_new_band_originations = sum(is_gse_newly_eligible_2009_vs_2007, na.rm = TRUE),
    gse_new_band_gse_sold_originations = sum(
      is_conventional_gse_sold & is_gse_newly_eligible_2009_vs_2007,
      na.rm = TRUE
    ),
    fha_new_band_originations = sum(is_fha_2009_vs_2008_band, na.rm = TRUE),
    fha_new_band_fha_originations = sum(is_fha & is_fha_2009_vs_2008_band, na.rm = TRUE),
    near_gse_old_limit_originations = sum(is_near_gse_old_limit, na.rm = TRUE),
    mean_loan_amount = mean(loan_amount, na.rm = TRUE),
    median_loan_amount = median(loan_amount, na.rm = TRUE),
    mean_applicant_income_000s = mean(applicant_income_000s, na.rm = TRUE)
  ), by = .(
    year = as_of_year,
    county_fips,
    state,
    county_name,
    msa_code,
    msa_name,
    gse_2007_proxy_limit,
    gse_2009_limit,
    fha_2008_limit,
    fha_2009_limit,
    gse_log_change_2009_vs_2007_proxy,
    fha_log_change_2009_vs_2008
  )]
}

script_dir <- get_script_dir()
raw_dir <- file.path(script_dir, "raw")
out_dir <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

exposure_path <- file.path(out_dir, "loan_limit_exposure_county_2008_2009.csv")
if (!file.exists(exposure_path)) {
  stop("Missing loan-limit exposure panel. Run build_loan_limit_exposure_panel.R first.")
}
exposure <- fread(
  exposure_path,
  colClasses = c(
    county_fips = "character",
    county_code = "character",
    msa_code = "character",
    md_code = "character"
  )
)

year_text <- Sys.getenv("YEARS", unset = "2007,2008,2009")
years <- as.integer(strsplit(year_text, ",")[[1]])
if (any(is.na(years))) {
  stop(sprintf("Invalid YEARS value: %s", year_text))
}

county_list <- lapply(years, function(yy) {
  zip_path <- file.path(
    raw_dir,
    sprintf("hmda_%d_nationwide_first-lien-owner-occupied-1-4-family-records_labels.zip", yy)
  )
  message(sprintf("Aggregating HMDA %d", yy))
  read_hmda_year(zip_path, exposure, yy)
})

county <- rbindlist(county_list, fill = TRUE)
setorder(county, county_fips, year)
fwrite(county, file.path(out_dir, "hmda_county_first_stage.csv"))

cbsa <- county[, .(
  counties = uniqueN(county_fips),
  originations = sum(originations, na.rm = TRUE),
  fha_originations = sum(fha_originations, na.rm = TRUE),
  conventional_gse_sold_originations = sum(conventional_gse_sold_originations, na.rm = TRUE),
  gse_new_band_originations = sum(gse_new_band_originations, na.rm = TRUE),
  gse_new_band_gse_sold_originations = sum(gse_new_band_gse_sold_originations, na.rm = TRUE),
  fha_new_band_originations = sum(fha_new_band_originations, na.rm = TRUE),
  fha_new_band_fha_originations = sum(fha_new_band_fha_originations, na.rm = TRUE),
  near_gse_old_limit_originations = sum(near_gse_old_limit_originations, na.rm = TRUE),
  mean_gse_log_change_2009_vs_2007_proxy = mean(gse_log_change_2009_vs_2007_proxy, na.rm = TRUE),
  max_gse_log_change_2009_vs_2007_proxy = max(gse_log_change_2009_vs_2007_proxy, na.rm = TRUE),
  mean_fha_log_change_2009_vs_2008 = mean(fha_log_change_2009_vs_2008, na.rm = TRUE)
), by = .(year, msa_code, msa_name)]
setorder(cbsa, msa_code, year)
fwrite(cbsa, file.path(out_dir, "hmda_cbsa_first_stage.csv"))

message("Wrote:")
message("  ", file.path(out_dir, "hmda_county_first_stage.csv"))
message("  ", file.path(out_dir, "hmda_cbsa_first_stage.csv"))
