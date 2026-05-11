#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) == 0) {
  stop("Could not determine script path from commandArgs().")
}
script_dir <- dirname(normalizePath(script_path))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(repo_root, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sample_mode <- Sys.getenv("ROOM_AUDIT_SAMPLE", unset = "full_extract")
bin_mode <- Sys.getenv("ROOM_BIN_MODE", unset = "coarse")

acs_path <- file.path(
  repo_root,
  "code",
  "data",
  "Spatial_aggregate_withmicrodata",
  "processed_data",
  "yearly_rds_v7",
  "fertility_microdata_2023.rds"
)
lookup_path <- file.path(
  repo_root,
  "code",
  "data",
  "mms_center_periphery",
  "data",
  "puma_mms_lookup_2020.csv"
)

if (!file.exists(acs_path)) {
  stop(sprintf("Missing ACS extract: %s", acs_path))
}

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[keep]
  w <- w[keep]
  if (!length(x)) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(
    probs,
    function(p) x[which(cw >= p)[1]],
    numeric(1)
  )
}

summarize_var <- function(dt, value_var, weight_var) {
  x <- dt[[value_var]]
  w <- dt[[weight_var]]
  q <- weighted_quantile(x, w, c(0.25, 0.50, 0.75))
  data.table(
    households = nrow(dt),
    hh_weight = sum(w, na.rm = TRUE),
    mean = weighted.mean(x, w, na.rm = TRUE),
    p25 = q[1],
    p50 = q[2],
    p75 = q[3],
    share_ge_7 = weighted.mean(x >= 7, w, na.rm = TRUE),
    share_ge_8 = weighted.mean(x >= 8, w, na.rm = TRUE),
    share_ge_11 = weighted.mean(x >= 11, w, na.rm = TRUE)
  )
}

bin_rooms <- function(x) {
  if (identical(bin_mode, "split_5_6")) {
    return(cut(
      x,
      breaks = c(-Inf, 4, 5, 6, 8, 10, Inf),
      labels = c("<=4", "5", "6", "7-8", "9-10", "11+"),
      right = TRUE
    ))
  }

  cut(
    x,
    breaks = c(-Inf, 4, 6, 8, 10, Inf),
    labels = c("<=4", "5-6", "7-8", "9-10", "11+"),
    right = TRUE
  )
}

bin_bedrooms <- function(x) {
  cut(
    x,
    breaks = c(-Inf, 1, 2, 3, Inf),
    labels = c("<=1", "2", "3", "4+"),
    right = TRUE
  )
}

summarize_bins <- function(dt, value_var, weight_var, bin_fun) {
  work <- copy(dt)
  work[, bin := bin_fun(get(value_var))]
  work <- work[!is.na(bin)]
  out <- work[, .(hh_weight = sum(get(weight_var), na.rm = TRUE)), by = bin]
  total_w <- sum(out$hh_weight)
  out[, share := hh_weight / total_w]
  setorder(out, bin)
  out[]
}

build_household_file <- function(dt) {
  dt <- dt[
    gq == 1 &
      ownershp %in% c(1, 2) &
      pernum == 1 &
      relate == 1 &
      !is.na(hhwt) &
      hhwt > 0
  ]
  dt[, tenure := fifelse(ownershp == 1, "Owner", "Renter")]
  dt[, child_bin := fifelse(
    is.na(nchild), NA_character_,
    fifelse(nchild <= 0, "0", fifelse(nchild == 1, "1", "2+"))
  )]
  dt[, age_window := "all"]
  dt_age <- dt[age >= 25 & age <= 45]
  dt_age[, age_window := "25_45"]
  base <- rbindlist(list(dt, dt_age), use.names = TRUE)
  base_all <- copy(base)
  base_all[, child_bin := "all"]
  rbindlist(list(base_all, base), use.names = TRUE)
}

apply_sample_filter <- function(dt) {
  if (identical(sample_mode, "full_extract")) {
    return(dt)
  }
  if (!identical(sample_mode, "mms_big_metros")) {
    stop(sprintf("Unknown ROOM_AUDIT_SAMPLE: %s", sample_mode))
  }
  if (!file.exists(lookup_path)) {
    stop(sprintf("Missing MMS lookup: %s", lookup_path))
  }

  lookup <- fread(lookup_path)[
    ,
    .(
      statefip = as.integer(statefip),
      puma = as.integer(puma),
      met2013 = as.integer(cbsacode),
      mms_location
    )
  ]
  lookup <- unique(lookup)
  dt[, `:=`(
    statefip = as.integer(statefip),
    puma = as.integer(puma),
    met2013 = as.integer(met2013)
  )]
  merge(dt, lookup, by = c("statefip", "puma", "met2013"), all = FALSE)
}

write_text_report <- function(path, rooms_all, rooms_2545, bedrooms_all, rooms_bins_2545) {
  fid <- file(path, open = "wt")
  on.exit(close(fid), add = TRUE)

  cat("ACS room-distribution audit\n", file = fid)
  cat(sprintf("date = %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = fid)
  cat(sprintf("source = %s\n\n", acs_path), file = fid)

  cat("Sample definition\n", file = fid)
  if (identical(sample_mode, "mms_big_metros")) {
    cat("- 2023 ACS household-head records restricted to matched MMS large-metro PUMAs\n", file = fid)
    cat(sprintf("- lookup = %s\n", lookup_path), file = fid)
  } else {
    cat("- 2023 ACS household-head records from the full local extract\n", file = fid)
  }
  cat("- occupied housing units only: gq == 1, ownershp in {1,2}\n", file = fid)
  cat("- household weights: hhwt\n", file = fid)
  cat("- child count uses household-head nchild from the extract\n\n", file = fid)
  cat("Binning\n", file = fid)
  cat(sprintf("- ROOM_BIN_MODE = %s\n\n", bin_mode), file = fid)

  cat("All-household room summaries by tenure\n", file = fid)
  cat(paste(capture.output(print(rooms_all)), collapse = "\n"), file = fid)
  cat("\n", file = fid)

  cat("All-household bedroom summaries by tenure\n", file = fid)
  cat(paste(capture.output(print(bedrooms_all)), collapse = "\n"), file = fid)
  cat("\n", file = fid)

  cat("Age 25-45 room summaries by tenure x child bin\n", file = fid)
  cat(paste(capture.output(print(rooms_2545)), collapse = "\n"), file = fid)
  cat("\n", file = fid)

  cat("Age 25-45 room-bin shares by tenure x child bin\n", file = fid)
  cat(paste(capture.output(print(rooms_bins_2545)), collapse = "\n"), file = fid)
}

acs_raw <- as.data.table(readRDS(acs_path))
acs_raw <- apply_sample_filter(acs_raw)
acs_hh <- build_household_file(acs_raw)
rm(acs_raw)
gc()

rooms_dt <- acs_hh[is.finite(rooms) & rooms > 0 & rooms < 90]
bedrooms_dt <- acs_hh[is.finite(bedrooms) & bedrooms >= 0 & bedrooms < 90]

rooms_summary <- rooms_dt[, summarize_var(.SD, "rooms", "hhwt"), by = .(age_window, tenure, child_bin)]
bedrooms_summary <- bedrooms_dt[, summarize_var(.SD, "bedrooms", "hhwt"), by = .(age_window, tenure, child_bin)]

rooms_bins <- rooms_dt[, summarize_bins(.SD, "rooms", "hhwt", bin_rooms), by = .(age_window, tenure, child_bin)]
bedrooms_bins <- bedrooms_dt[, summarize_bins(.SD, "bedrooms", "hhwt", bin_bedrooms), by = .(age_window, tenure, child_bin)]

rooms_all <- rooms_summary[age_window == "all" & child_bin %in% c("all", "0", "1", "2+")]
rooms_2545 <- rooms_summary[age_window == "25_45" & child_bin %in% c("all", "0", "1", "2+")]
bedrooms_all <- bedrooms_summary[age_window == "all" & child_bin %in% c("all", "0", "1", "2+")]
rooms_bins_2545 <- rooms_bins[age_window == "25_45" & child_bin %in% c("all", "0", "1", "2+")]

setorder(rooms_summary, age_window, tenure, child_bin)
setorder(bedrooms_summary, age_window, tenure, child_bin)
setorder(rooms_bins, age_window, tenure, child_bin, bin)
setorder(bedrooms_bins, age_window, tenure, child_bin, bin)

fwrite(rooms_summary, file.path(out_dir, "acs_2023_rooms_summary.csv"))
fwrite(bedrooms_summary, file.path(out_dir, "acs_2023_bedrooms_summary.csv"))
fwrite(rooms_bins, file.path(out_dir, "acs_2023_rooms_bins.csv"))
fwrite(bedrooms_bins, file.path(out_dir, "acs_2023_bedrooms_bins.csv"))

write_text_report(
  file.path(out_dir, "acs_2023_room_distribution_audit.txt"),
  rooms_all,
  rooms_2545,
  bedrooms_all,
  rooms_bins_2545
)

cat(file.path(out_dir, "acs_2023_room_distribution_audit.txt"), "\n")
