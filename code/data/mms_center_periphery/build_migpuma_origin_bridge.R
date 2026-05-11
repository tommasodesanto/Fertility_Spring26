#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(tidyr)
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

read_rel_2010 <- function(path) {
  read_excel(
    path,
    skip = 2,
    col_names = c("statefip", "puma", "migplac1", "migpuma1")
  ) %>%
    mutate(across(everything(), as.integer))
}

read_rel_2020 <- function(path) {
  x <- read_excel(path)
  transmute(
    x,
    statefip = as.integer(.data[["State of Residence (ST)"]]),
    puma = as.integer(.data[["PUMA"]]),
    migplac1 = as.integer(.data[["Place of Work State (PWSTATE2) or Migration State (MIGPLAC1)"]]),
    migpuma1 = as.integer(.data[["PWPUMA00 or MIGPUMA1"]])
  )
}

safe_share <- function(x, denom) {
  ifelse(denom > 0, x / denom, NA_real_)
}

collapse_mms_location <- function(x, include_middle_in_center) {
  if (!include_middle_in_center) {
    return(x)
  }
  case_when(
    x == "middle" ~ "center",
    TRUE ~ x
  )
}

script_dir <- get_script_dir()
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
crosswalk_dir <- file.path(script_dir, "data", "ipums_crosswalks")
include_middle_in_center <- tolower(Sys.getenv("MMS_INCLUDE_MIDDLE_IN_CENTER", "false")) %in% c("1", "true", "yes", "y")

lookup_2010 <- read_csv(file.path(data_dir, "puma_mms_lookup_2010.csv"), show_col_types = FALSE)
lookup_2020 <- read_csv(file.path(data_dir, "puma_mms_lookup_2020.csv"), show_col_types = FALSE)

rel_2010 <- read_rel_2010(file.path(crosswalk_dir, "puma_migpuma1_pwpuma00_2010.xls")) %>%
  mutate(period = "pre2022")
rel_2020 <- read_rel_2020(file.path(crosswalk_dir, "puma_migpuma1_pwpuma00_2020.xls")) %>%
  mutate(period = "post2021")

lookup_stack <- bind_rows(
  lookup_2010 %>%
    transmute(
      period = "pre2022",
      statefip,
      puma,
      cbsacode,
      tract_pop2000,
      mms_location
    ),
  lookup_2020 %>%
    transmute(
      period = "post2021",
      statefip,
      puma,
      cbsacode,
      tract_pop2000,
      mms_location
    )
  ) %>%
  mutate(mms_location = collapse_mms_location(mms_location, include_middle_in_center))

rel_stack <- bind_rows(rel_2010, rel_2020)

bridge_base <- rel_stack %>%
  left_join(lookup_stack, by = c("period", "statefip", "puma")) %>%
  filter(!is.na(cbsacode), !is.na(mms_location))

bridge <- bridge_base %>%
  group_by(period, migplac1, migpuma1, cbsacode) %>%
  summarise(
    matched_pop_weight = sum(tract_pop2000, na.rm = TRUE),
    n_puma = n(),
    n_location = n_distinct(mms_location),
    center_share = safe_share(sum(tract_pop2000[mms_location == "center"], na.rm = TRUE), matched_pop_weight),
    periphery_share = safe_share(sum(tract_pop2000[mms_location == "periphery"], na.rm = TRUE), matched_pop_weight),
    middle_share = safe_share(sum(tract_pop2000[mms_location == "middle"], na.rm = TRUE), matched_pop_weight),
    .groups = "drop"
  ) %>%
  mutate(
    center_share = coalesce(center_share, 0),
    periphery_share = coalesce(periphery_share, 0),
    middle_share = coalesce(middle_share, 0),
    center_periphery_share = center_share + periphery_share,
    dominant_location = case_when(
      center_share >= periphery_share & center_share >= middle_share ~ "center",
      periphery_share >= center_share & periphery_share >= middle_share ~ "periphery",
      TRUE ~ "middle"
    ),
    dominant_share = pmax(center_share, periphery_share, middle_share)
  )

diag_summary <- bridge %>%
  summarise(
    n_bridge_groups = n(),
    pure_groups = sum(dominant_share == 1, na.rm = TRUE),
    dominant_ge_90 = mean(dominant_share >= 0.9, na.rm = TRUE),
    dominant_ge_75 = mean(dominant_share >= 0.75, na.rm = TRUE),
    cp_only_groups = mean(center_periphery_share == 1, na.rm = TRUE),
    cp_mixed_groups = mean(center_share > 0 & periphery_share > 0, na.rm = TRUE)
  )

state_consistency <- rel_stack %>%
  summarise(equal_share = mean(statefip == migplac1, na.rm = TRUE))

write_csv(bridge, file.path(data_dir, "migpuma_mms_origin_bridge.csv"))
write_csv(diag_summary, file.path(data_dir, "migpuma_mms_origin_bridge_diagnostics.csv"))

summary_lines <- c(
  "MIGPUMA-to-MMS origin bridge",
  sprintf("Middle absorbed into center: %s", ifelse(include_middle_in_center, "yes", "no")),
  sprintf("Bridge groups: %s", format(diag_summary$n_bridge_groups, big.mark = ",")),
  sprintf("Pure single-location groups: %s", format(diag_summary$pure_groups, big.mark = ",")),
  sprintf("Dominant share >= 0.90: %.1f%%", 100 * diag_summary$dominant_ge_90),
  sprintf("Dominant share >= 0.75: %.1f%%", 100 * diag_summary$dominant_ge_75),
  sprintf("Center/periphery only groups: %.1f%%", 100 * diag_summary$cp_only_groups),
  sprintf("Mixed center-periphery groups: %.1f%%", 100 * diag_summary$cp_mixed_groups),
  sprintf("Relationship-file state == MIGPLAC1 share: %.1f%%", 100 * state_consistency$equal_share)
)

writeLines(summary_lines, file.path(data_dir, "migpuma_mms_origin_bridge_summary.txt"))
message("Saved bridge: ", file.path(data_dir, "migpuma_mms_origin_bridge.csv"))
