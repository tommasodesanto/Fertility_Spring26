#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", file_arg[startsWith(file_arg, file_flag)])
  if (length(script_path) == 0) {
    stop("Run this script with Rscript so --file= is available.")
  }
  dirname(normalizePath(script_path[1]))
}

weighted_mean_safe <- function(x, w) {
  ok <- is.na(x) == FALSE & is.na(w) == FALSE
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

weighted_median_safe <- function(x, w) {
  ok <- is.na(x) == FALSE & is.na(w) == FALSE
  if (!any(ok)) {
    return(NA_real_)
  }
  x_ok <- x[ok]
  w_ok <- w[ok]
  ord <- order(x_ok)
  x_ok <- x_ok[ord]
  w_ok <- w_ok[ord]
  cw <- cumsum(w_ok) / sum(w_ok)
  x_ok[which(cw >= 0.5)[1]]
}

ownership_table <- function(dt, label, age_min, age_max) {
  sub <- dt[age >= age_min & age <= age_max]
  no_children_rate <- sub[nchild == 0 & is.na(nchild) == FALSE, weighted.mean(owner, perwt)]
  parent_rate <- sub[nchild > 0 & is.na(nchild) == FALSE, weighted.mean(owner, perwt)]
  newparent_rate <- sub[eldch < 4 & eldch != 99 & nchild > 0 & is.na(eldch) == FALSE & is.na(nchild) == FALSE, weighted.mean(owner, perwt)]
  by_loc <- sub[, .(owner_rate = weighted.mean(owner, perwt)), by = mms_location]
  center_rate <- by_loc[mms_location == "center", owner_rate]
  periphery_rate <- by_loc[mms_location == "periphery", owner_rate]
  data.table(
    moment_group = "ownership",
    window = label,
    overall = sub[, weighted.mean(owner, perwt)],
    center = center_rate,
    periphery = periphery_rate,
    p_minus_c = periphery_rate - center_rate,
    parent_minus_nochildren = parent_rate - no_children_rate,
    newparent_minus_nochildren = newparent_rate - no_children_rate
  )
}

rent_table <- function(dt, label, age_min, age_max) {
  sub <- dt[age >= age_min & age <= age_max & renter == TRUE & is.na(rent) == FALSE & rent > 0]
  sub <- sub[is.na(rooms) == FALSE & rooms > 0]
  sub[, rent_per_room := rent / rooms]
  by_loc <- sub[, .(
    mean_rent = weighted.mean(rent, perwt),
    median_rent = weighted_median_safe(rent, perwt),
    mean_rent_per_room = weighted.mean(rent_per_room, perwt),
    median_rent_per_room = weighted_median_safe(rent_per_room, perwt)
  ), by = mms_location]
  center_mean <- by_loc[mms_location == "center", mean_rent]
  periphery_mean <- by_loc[mms_location == "periphery", mean_rent]
  center_median <- by_loc[mms_location == "center", median_rent]
  periphery_median <- by_loc[mms_location == "periphery", median_rent]
  center_mean_rr <- by_loc[mms_location == "center", mean_rent_per_room]
  periphery_mean_rr <- by_loc[mms_location == "periphery", mean_rent_per_room]
  center_median_rr <- by_loc[mms_location == "center", median_rent_per_room]
  periphery_median_rr <- by_loc[mms_location == "periphery", median_rent_per_room]
  data.table(
    moment_group = "rent",
    window = label,
    mean_center = center_mean,
    mean_periphery = periphery_mean,
    mean_ratio_c_to_p = center_mean / periphery_mean,
    median_center = center_median,
    median_periphery = periphery_median,
    median_ratio_c_to_p = center_median / periphery_median,
    mean_rent_per_room_center = center_mean_rr,
    mean_rent_per_room_periphery = periphery_mean_rr,
    mean_rent_per_room_ratio_c_to_p = center_mean_rr / periphery_mean_rr,
    median_rent_per_room_center = center_median_rr,
    median_rent_per_room_periphery = periphery_median_rr,
    median_rent_per_room_ratio_c_to_p = center_median_rr / periphery_median_rr
  )
}

tfr_table <- function(dt) {
  fert <- dt[sex == 2 & age >= 15 & age <= 44 & fertyr %in% c(1, 2)]
  fert[, had_birth := as.integer(fertyr == 2)]
  fert[, age_group := cut(
    age,
    breaks = c(14, 19, 24, 29, 34, 39, 44),
    labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44")
  )]
  asfr <- fert[, .(asfr = weighted.mean(had_birth, perwt) * 1000), by = .(mms_location, age_group)]
  tfr <- asfr[, .(tfr = sum(asfr) * 5 / 1000), by = mms_location]
  birth_rate <- fert[, .(birth_rate = weighted.mean(had_birth, perwt) * 1000), by = mms_location]
  center_tfr <- tfr[mms_location == "center", tfr]
  periphery_tfr <- tfr[mms_location == "periphery", tfr]
  center_br <- birth_rate[mms_location == "center", birth_rate]
  periphery_br <- birth_rate[mms_location == "periphery", birth_rate]
  data.table(
    moment_group = "fertility",
    window = "women_15_44",
    tfr_center = center_tfr,
    tfr_periphery = periphery_tfr,
    tfr_p_minus_c = periphery_tfr - center_tfr,
    birth_rate_center = center_br,
    birth_rate_periphery = periphery_br
  )
}

move_table <- function(dt, label, age_min, age_max) {
  sub <- copy(dt[age >= age_min & age <= age_max])
  sub[, moved1y := as.integer(migrate1 >= 2 & is.na(migrate1) == FALSE)]
  sub[, same_msa := as.integer(moved1y == 1 & met2013 > 0 & migmet131 > 0 & met2013 == migmet131)]
  sub[, diff_msa := as.integer(moved1y == 1 & met2013 > 0 & migmet131 > 0 & met2013 != migmet131)]
  movers <- sub[moved1y == 1]
  data.table(
    moment_group = "moves",
    window = label,
    any_move_rate = weighted.mean(sub$moved1y, sub$perwt),
    within_cbsa_move_rate = weighted.mean(sub$same_msa, sub$perwt),
    across_cbsa_move_rate = weighted.mean(sub$diff_msa, sub$perwt),
    within_share_among_movers = weighted.mean(movers$same_msa, movers$perwt)
  )
}

script_dir <- get_script_dir()
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_suffix <- Sys.getenv("MMS_OUTPUT_SUFFIX", "")
out_dir <- file.path(script_dir, paste0("output", out_suffix))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

middle_target <- tolower(Sys.getenv("MMS_MIDDLE_TARGET", ""))
if (middle_target == "") {
  include_middle_in_center <- tolower(Sys.getenv("MMS_INCLUDE_MIDDLE_IN_CENTER", "true")) %in% c("1", "true", "yes", "y")
  middle_target <- if (include_middle_in_center) "center" else "drop"
}
if (!(middle_target %in% c("drop", "center", "periphery"))) {
  stop("MMS_MIDDLE_TARGET must be one of: drop, center, periphery")
}

lookup_2010 <- as.data.table(fread(file.path(data_dir, "puma_mms_lookup_2010.csv")))
lookup_2020 <- as.data.table(fread(file.path(data_dir, "puma_mms_lookup_2020.csv")))
coverage <- as.data.table(fread(file.path(out_dir, "mms_move_coverage_summary.csv")))
origin_transitions <- as.data.table(fread(file.path(out_dir, "mms_origin_transition_summary.csv")))

year_dir <- "code/data/Spatial_aggregate_withmicrodata/processed_data/yearly_rds_v7"
years <- 2012:2023
yearly_list <- lapply(years, function(yy) {
  path <- file.path(year_dir, sprintf("fertility_microdata_%d.rds", yy))
  x <- as.data.table(readRDS(path))
  x <- x[
    age >= 15 & age <= 55 &
      gq %in% c(1, 2),
    .(
      year = as.integer(year),
      statefip = as.integer(statefip),
      puma = as.integer(puma),
      metro = as.integer(metro),
      met2013 = as.integer(met2013),
      age = as.integer(age),
      sex = as.integer(sex),
      perwt = as.numeric(perwt),
      gq = as.integer(gq),
      ownershp = as.integer(ownershp),
      rent = as.numeric(rent),
      rooms = as.numeric(rooms),
      nchild = as.numeric(nchild),
      eldch = as.numeric(eldch),
      fertyr = as.integer(fertyr),
      migrate1 = as.integer(migrate1),
      migmet131 = as.integer(migmet131)
    )
  ]
  x
})

dt <- rbindlist(yearly_list, use.names = TRUE)
rm(yearly_list)
gc()
dt[, owner := ownershp == 1L]
dt[, renter := ownershp == 2L]

pre <- merge(
  dt[year <= 2021],
  lookup_2010,
  by.x = c("statefip", "puma", "met2013"),
  by.y = c("statefip", "puma", "cbsacode"),
  all.x = TRUE,
  sort = FALSE
)
post <- merge(
  dt[year >= 2022],
  lookup_2020,
  by.x = c("statefip", "puma", "met2013"),
  by.y = c("statefip", "puma", "cbsacode"),
  all.x = TRUE,
  sort = FALSE
)
rm(dt)
gc()

dt_mms <- rbindlist(list(pre, post), fill = TRUE, use.names = TRUE)
rm(pre, post)
gc()

if (middle_target == "center") {
  dt_mms[mms_location == "middle", mms_location := "center"]
} else if (middle_target == "periphery") {
  dt_mms[mms_location == "middle", mms_location := "periphery"]
}
dt_mms <- dt_mms[mms_location %in% c("center", "periphery")]

ownership_results <- rbindlist(list(
  ownership_table(dt_mms, "22_45", 22, 45),
  ownership_table(dt_mms, "30_45", 30, 45),
  ownership_table(dt_mms, "30_55", 30, 55)
), fill = TRUE)

rent_results <- rbindlist(list(
  rent_table(dt_mms, "22_45", 22, 45),
  rent_table(dt_mms, "22_55", 22, 55),
  rent_table(dt_mms, "30_55", 30, 55)
), fill = TRUE)

fertility_results <- tfr_table(dt_mms)

move_results <- rbindlist(list(
  move_table(dt_mms, "22_45", 22, 45),
  move_table(dt_mms, "25_45", 25, 45),
  move_table(dt_mms, "30_55", 30, 55)
), fill = TRUE)

all_matched <- coverage[metric == "All households with MMS-matched destination", weighted_n]
within_origin_matched <- coverage[metric == "Within-CBSA movers with bridge-matched origin mass", weighted_n]
switch_mass <- origin_transitions[origin_label != dest_label, sum(pop_weight)]
move_switch_results <- data.table(
  moment_group = "moves_cp_switch",
  window = "22_45_existing_outputs",
  cp_switch_rate_all_matched = switch_mass / all_matched,
  cp_switch_share_among_within = switch_mass / within_origin_matched
)

fwrite(ownership_results, file.path(out_dir, "reconstructed_ownership_targets.csv"))
fwrite(rent_results, file.path(out_dir, "reconstructed_rent_targets.csv"))
fwrite(fertility_results, file.path(out_dir, "reconstructed_fertility_targets.csv"))
fwrite(move_results, file.path(out_dir, "reconstructed_move_targets.csv"))
fwrite(move_switch_results, file.path(out_dir, "reconstructed_cp_switch_targets.csv"))

cat("=== OWNERSHIP ===\n")
print(ownership_results)
cat("\n=== RENT ===\n")
print(rent_results)
cat("\n=== FERTILITY ===\n")
print(fertility_results)
cat("\n=== MOVES ===\n")
print(move_results)
cat("\n=== CP SWITCH ===\n")
print(move_switch_results)
