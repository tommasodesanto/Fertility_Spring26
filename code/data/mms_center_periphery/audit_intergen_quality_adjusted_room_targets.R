#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
})

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", file_arg[startsWith(file_arg, file_flag)])
  if (length(script_path) == 0) {
    stop("Run this script with Rscript so --file= is available.")
  }
  dirname(normalizePath(script_path[1]))
}

collapse_mms_location <- function(x, middle_target) {
  if (middle_target == "drop") {
    return(x)
  }
  if (middle_target == "center") {
    return(case_when(x == "middle" ~ "center", TRUE ~ x))
  }
  if (middle_target == "periphery") {
    return(case_when(x == "middle" ~ "periphery", TRUE ~ x))
  }
  stop("Invalid MMS_MIDDLE_TARGET: ", middle_target)
}

summarise_tenure <- function(data, group_vars) {
  data %>%
    group_by(across(all_of(group_vars)), tenure) %>%
    summarise(
      n = n(),
      weight = sum(hhwt, na.rm = TRUE),
      mean_rooms = weighted_mean_safe(rooms, hhwt),
      mean_bedrooms = weighted_mean_safe(bedrooms_count, hhwt),
      mean_rent_to_income = weighted_mean_safe(rent_to_income, hhwt),
      .groups = "drop"
    ) %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(weight_share_in_group = weight / sum(weight, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(tenure) %>%
    mutate(weight_share_within_tenure = weight / sum(weight, na.rm = TRUE)) %>%
    ungroup()
}

make_raw_gap_row <- function(data, raw_owner_mean, raw_renter_mean, model_gap,
                             model_owner_mean, model_renter_mean) {
  owner_weight <- sum(data$hhwt[data$owner], na.rm = TRUE)
  renter_weight <- sum(data$hhwt[data$renter], na.rm = TRUE)
  tibble(
    spec = "raw_current_target",
    description = "Raw weighted owner and renter means; reproduces the current ACS room target.",
    target_object = "owner-renter mean room gap among age 30-55 childless household heads",
    target_owner_mean_rooms = raw_owner_mean,
    model_owner_mean_rooms = model_owner_mean,
    model_minus_target_owner_mean_rooms = model_owner_mean - raw_owner_mean,
    target_renter_mean_rooms = raw_renter_mean,
    model_renter_mean_rooms = model_renter_mean,
    model_minus_target_renter_mean_rooms = model_renter_mean - raw_renter_mean,
    target_owner_minus_renter_rooms = raw_owner_mean - raw_renter_mean,
    model_owner_minus_renter_rooms = model_gap,
    model_minus_target_gap = model_gap - (raw_owner_mean - raw_renter_mean),
    raw_target_owner_mean_rooms = raw_owner_mean,
    raw_target_renter_mean_rooms = raw_renter_mean,
    raw_target_owner_minus_renter_rooms = raw_owner_mean - raw_renter_mean,
    change_from_raw_gap = 0,
    n = nrow(data),
    weight = sum(data$hhwt, na.rm = TRUE),
    matched_cell_weight = sum(data$hhwt, na.rm = TRUE),
    matched_cell_weight_share = 1,
    matched_owner_weight_share = 1,
    matched_renter_weight_share = 1,
    contrast_weighted_gap = raw_owner_mean - raw_renter_mean,
    controls = "none",
    status = "raw current target; not quality adjusted"
  )
}

compute_cell_gap <- function(data, spec, description, group_vars, baseline_status,
                             raw_owner_mean, raw_renter_mean, model_gap,
                             model_owner_mean, model_renter_mean) {
  fit_data <- data
  for (g in group_vars) {
    fit_data <- fit_data %>% filter(!is.na(.data[[g]]))
  }
  cell_by_tenure <- fit_data %>%
    group_by(across(all_of(group_vars)), tenure) %>%
    summarise(
      n = n(),
      weight = sum(hhwt, na.rm = TRUE),
      mean_rooms = weighted_mean_safe(rooms, hhwt),
      .groups = "drop"
    )
  owner_cells <- cell_by_tenure %>%
    filter(tenure == "owner") %>%
    select(all_of(group_vars), owner_n = n, owner_weight = weight, owner_mean_rooms = mean_rooms)
  renter_cells <- cell_by_tenure %>%
    filter(tenure == "renter") %>%
    select(all_of(group_vars), renter_n = n, renter_weight = weight, renter_mean_rooms = mean_rooms)
  matched <- inner_join(owner_cells, renter_cells, by = group_vars) %>%
    mutate(
      cell_gap = owner_mean_rooms - renter_mean_rooms,
      total_cell_weight = owner_weight + renter_weight,
      contrast_weight = owner_weight * renter_weight / (owner_weight + renter_weight)
    )

  target_owner_mean <- weighted_mean_safe(matched$owner_mean_rooms, matched$total_cell_weight)
  target_renter_mean <- weighted_mean_safe(matched$renter_mean_rooms, matched$total_cell_weight)
  target_gap <- target_owner_mean - target_renter_mean
  contrast_gap <- weighted_mean_safe(matched$cell_gap, matched$contrast_weight)

  tibble(
    spec = spec,
    description = description,
    target_object = "common-composition owner-renter room gap among age 30-55 childless household heads",
    target_owner_mean_rooms = target_owner_mean,
    model_owner_mean_rooms = model_owner_mean,
    model_minus_target_owner_mean_rooms = model_owner_mean - target_owner_mean,
    target_renter_mean_rooms = target_renter_mean,
    model_renter_mean_rooms = model_renter_mean,
    model_minus_target_renter_mean_rooms = model_renter_mean - target_renter_mean,
    target_owner_minus_renter_rooms = target_gap,
    model_owner_minus_renter_rooms = model_gap,
    model_minus_target_gap = model_gap - target_gap,
    raw_target_owner_mean_rooms = raw_owner_mean,
    raw_target_renter_mean_rooms = raw_renter_mean,
    raw_target_owner_minus_renter_rooms = raw_owner_mean - raw_renter_mean,
    change_from_raw_gap = target_gap - (raw_owner_mean - raw_renter_mean),
    n = sum(matched$owner_n + matched$renter_n, na.rm = TRUE),
    weight = sum(fit_data$hhwt, na.rm = TRUE),
    matched_cell_weight = sum(matched$total_cell_weight, na.rm = TRUE),
    matched_cell_weight_share = sum(matched$total_cell_weight, na.rm = TRUE) / sum(fit_data$hhwt, na.rm = TRUE),
    matched_owner_weight_share = sum(matched$owner_weight, na.rm = TRUE) / sum(fit_data$hhwt[fit_data$owner], na.rm = TRUE),
    matched_renter_weight_share = sum(matched$renter_weight, na.rm = TRUE) / sum(fit_data$hhwt[fit_data$renter], na.rm = TRUE),
    contrast_weighted_gap = contrast_gap,
    controls = paste(group_vars, collapse = " x "),
    status = baseline_status
  )
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "../../.."))
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_dir <- file.path(script_dir, "output_intergen_quality_adjusted_room_targets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

middle_target <- tolower(Sys.getenv("MMS_MIDDLE_TARGET", ""))
if (middle_target == "") {
  include_middle_in_center <- tolower(Sys.getenv("MMS_INCLUDE_MIDDLE_IN_CENTER", "true")) %in% c("1", "true", "yes", "y")
  middle_target <- if (include_middle_in_center) "center" else "drop"
}

extract_path <- file.path(repo_root, "code", "data", "Spatial_aggregate_withmicrodata", "raw_data", "extract27.dta")
lookup_2010_path <- file.path(data_dir, "puma_mms_lookup_2010.csv")
lookup_2020_path <- file.path(data_dir, "puma_mms_lookup_2020.csv")

if (!file.exists(extract_path)) {
  stop("Missing microdata extract: ", extract_path)
}
if (!file.exists(lookup_2010_path) || !file.exists(lookup_2020_path)) {
  stop("Missing MMS lookup files in ", data_dir)
}

lookup_2010 <- read_csv(lookup_2010_path, show_col_types = FALSE)
lookup_2020 <- read_csv(lookup_2020_path, show_col_types = FALSE)
dest_lookup <- bind_rows(
  lookup_2010 %>%
    transmute(
      lookup_period = "pre2022",
      statefip,
      puma,
      met2013 = cbsacode,
      cbsatitle,
      mms_location
    ),
  lookup_2020 %>%
    transmute(
      lookup_period = "post2021",
      statefip,
      puma,
      met2013 = cbsacode,
      cbsatitle,
      mms_location
    )
)

message("Reading ACS/IPUMS household-head extract for quality-adjusted room audit...")
df <- read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "met2013", "gq", "pernum", "hhwt",
    "ownershp", "rooms", "rent", "hhincome", "nchild", "yngch", "age",
    "builtyr2", "unitsstr", "bedrooms"
  )
) %>%
  mutate(
    unitsstr_label = as.character(as_factor(unitsstr)),
    builtyr2_label = as.character(as_factor(builtyr2))
  ) %>%
  transmute(
    year = as.integer(year),
    statefip = as.integer(statefip),
    puma = as.integer(puma),
    met2013 = as.integer(met2013),
    gq = as.integer(gq),
    pernum = as.integer(pernum),
    hhwt = as.numeric(hhwt),
    ownershp = as.integer(ownershp),
    rooms = as.numeric(rooms),
    rent = as.numeric(rent),
    hhincome = as.numeric(hhincome),
    nchild = as.numeric(nchild),
    yngch = as.numeric(yngch),
    age = as.integer(age),
    builtyr2 = as.integer(builtyr2),
    builtyr2_label = builtyr2_label,
    unitsstr = as.integer(unitsstr),
    unitsstr_label = unitsstr_label,
    bedrooms = as.numeric(bedrooms),
    lookup_period = if_else(year <= 2021, "pre2022", "post2021")
  ) %>%
  filter(
    year >= 2012,
    met2013 > 0,
    gq %in% c(1, 2),
    pernum == 1,
    age >= 18,
    !is.na(hhwt),
    hhwt > 0,
    !is.na(rooms),
    rooms > 0
  ) %>%
  left_join(dest_lookup, by = c("lookup_period", "statefip", "puma", "met2013")) %>%
  mutate(
    mms_location = collapse_mms_location(mms_location, middle_target),
    in_mms_sample = mms_location %in% c("center", "periphery"),
    owner = ownershp == 1,
    renter = ownershp == 2,
    owner_int = as.integer(owner),
    tenure = case_when(owner ~ "owner", renter ~ "renter", TRUE ~ "other"),
    childless = nchild == 0 & !is.na(nchild),
    prime_30_55 = age >= 30 & age <= 55,
    age_bin_5 = cut(age, breaks = c(29, 34, 39, 44, 49, 55), include.lowest = TRUE, right = TRUE),
    rent_to_income = if_else(renter & rent > 0 & hhincome > 1000, 12 * rent / hhincome, NA_real_),
    bedrooms_count = case_when(
      bedrooms == 1 ~ 0,
      bedrooms >= 2 & bedrooms <= 13 ~ bedrooms - 2,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(in_mms_sample, tenure %in% c("owner", "renter"))

audit_sample <- df %>%
  filter(prime_30_55, childless)

owner_childless <- audit_sample %>% filter(owner)
renter_childless <- audit_sample %>% filter(renter)

raw_owner_mean <- weighted_mean_safe(owner_childless$rooms, owner_childless$hhwt)
raw_renter_mean <- weighted_mean_safe(renter_childless$rooms, renter_childless$hhwt)

# Current diagnostic model values from the June 23 global room-gap readout,
# best 14-moment/scalar diagnostic point. These are comparison values only.
model_owner_mean <- 6.052569140341261
model_renter_mean <- 4.559938806863727
model_gap <- 1.4926303334775346

specs <- list(
  list(
    spec = "same_year_metro_mms_age",
    description = "Common-composition gap within year, CBSA, MMS center/periphery, and five-year age cells.",
    controls = c("year", "met2013", "mms_location", "age_bin_5"),
    status = "candidate location/year/age-adjusted diagnostic"
  ),
  list(
    spec = "same_year_metro_mms_age_structure",
    description = "Adds structure type to the common-composition cells.",
    controls = c("year", "met2013", "mms_location", "age_bin_5", "unitsstr_label"),
    status = "candidate location-and-structure-adjusted diagnostic"
  ),
  list(
    spec = "same_year_metro_mms_age_structure_vintage",
    description = "Adds structure type and building vintage to the common-composition cells.",
    controls = c("year", "met2013", "mms_location", "age_bin_5", "unitsstr_label", "builtyr2_label"),
    status = "aggressive quality/location adjustment diagnostic; do not use as target without review"
  ),
  list(
    spec = "diagnostic_same_year_metro_mms_age_structure_vintage_bedrooms",
    description = "Adds bedroom count to the common-composition cells; diagnostic only because bedrooms are a physical size margin.",
    controls = c("year", "met2013", "mms_location", "age_bin_5", "unitsstr_label", "builtyr2_label", "bedrooms_count"),
    status = "diagnostic only; bedrooms over-control the housing-size object"
  )
)

gap_results <- bind_rows(
  make_raw_gap_row(
    audit_sample,
    raw_owner_mean = raw_owner_mean,
    raw_renter_mean = raw_renter_mean,
    model_gap = model_gap,
    model_owner_mean = model_owner_mean,
    model_renter_mean = model_renter_mean
  ),
  lapply(specs, function(s) {
    compute_cell_gap(
      audit_sample,
      spec = s$spec,
      description = s$description,
      group_vars = s$controls,
      baseline_status = s$status,
      raw_owner_mean = raw_owner_mean,
      raw_renter_mean = raw_renter_mean,
      model_gap = model_gap,
      model_owner_mean = model_owner_mean,
      model_renter_mean = model_renter_mean
    )
  })
)

location_summary <- summarise_tenure(audit_sample, c("mms_location"))
structure_summary <- summarise_tenure(audit_sample, c("unitsstr_label"))
vintage_summary <- summarise_tenure(audit_sample, c("builtyr2_label"))

sample_summary <- tibble(
  statistic = c(
    "sample_n",
    "sample_weight",
    "owner_n",
    "renter_n",
    "owner_weight",
    "renter_weight",
    "owner_weight_share",
    "renter_weight_share",
    "raw_owner_mean_rooms",
    "raw_renter_mean_rooms",
    "raw_owner_minus_renter_mean_rooms",
    "model_owner_mean_rooms",
    "model_renter_mean_rooms",
    "model_owner_minus_renter_mean_rooms"
  ),
  value = c(
    nrow(audit_sample),
    sum(audit_sample$hhwt, na.rm = TRUE),
    nrow(owner_childless),
    nrow(renter_childless),
    sum(owner_childless$hhwt, na.rm = TRUE),
    sum(renter_childless$hhwt, na.rm = TRUE),
    sum(owner_childless$hhwt, na.rm = TRUE) / sum(audit_sample$hhwt, na.rm = TRUE),
    sum(renter_childless$hhwt, na.rm = TRUE) / sum(audit_sample$hhwt, na.rm = TRUE),
    raw_owner_mean,
    raw_renter_mean,
    raw_owner_mean - raw_renter_mean,
    model_owner_mean,
    model_renter_mean,
    model_gap
  )
)

write_csv(gap_results, file.path(out_dir, "intergen_quality_adjusted_room_gap_targets.csv"))
write_csv(location_summary, file.path(out_dir, "intergen_quality_adjusted_room_location_summary.csv"))
write_csv(structure_summary, file.path(out_dir, "intergen_quality_adjusted_room_structure_summary.csv"))
write_csv(vintage_summary, file.path(out_dir, "intergen_quality_adjusted_room_vintage_summary.csv"))
write_csv(sample_summary, file.path(out_dir, "intergen_quality_adjusted_room_sample_summary.csv"))

readme <- c(
  "# Intergen Quality-Adjusted Room Target Audit",
  "",
  "These outputs are diagnostic candidate target objects, not approved SMM moments.",
  "",
  "Main file:",
  "- `intergen_quality_adjusted_room_gap_targets.csv`: raw and common-composition owner-renter room-gap targets for prime-age childless household heads, with current diagnostic model values next to targets.",
  "",
  "Supporting files:",
  "- `intergen_quality_adjusted_room_location_summary.csv`: tenure room means by MMS center/periphery.",
  "- `intergen_quality_adjusted_room_structure_summary.csv`: tenure room means by structure type.",
  "- `intergen_quality_adjusted_room_vintage_summary.csv`: tenure room means by building vintage.",
  "- `intergen_quality_adjusted_room_sample_summary.csv`: raw sample and model comparison values.",
  "",
  "Baseline sample matches the existing intergen ACS housing target builder: ACS/IPUMS extract27 household heads, 2012-2023, matched MMS metro sample, middle collapsed to center, age 30-55, childless, owner or renter.",
  "",
  "The most aggressive quality/location adjustment is `same_year_metro_mms_age_structure_vintage`: owner and renter room means compared within exact year, CBSA, MMS center/periphery, five-year age, structure type, and building-vintage cells, then averaged over common-support cell population weights.",
  "",
  "Bedroom-count specifications are diagnostic only: bedroom count is itself a housing-size margin.",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote quality-adjusted room audit to ", out_dir)
