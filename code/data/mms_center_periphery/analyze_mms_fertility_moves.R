#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(fixest)
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

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

safe_percent <- function(x, digits = 0.1) {
  ifelse(is.na(x), "NA", percent(x, accuracy = digits))
}

extract_newparent <- function(model, label) {
  tibble(
    outcome = label,
    estimate = unname(coef(model)["newparent"]),
    se = unname(se(model)["newparent"]),
    p_value = unname(pvalue(model)["newparent"])
  )
}

extract_parent_coef <- function(model, term, label) {
  tibble(
    outcome = label,
    estimate = unname(coef(model)[term]),
    se = unname(se(model)[term]),
    p_value = unname(pvalue(model)[term])
  )
}

run_newparent_reg <- function(data, outcome, label, weight_var = "perwt") {
  y <- data[[outcome]]
  if (nrow(data) == 0 || sum(!is.na(y)) == 0 || dplyr::n_distinct(y[!is.na(y)]) < 2) {
    return(tibble(outcome = label, estimate = NA_real_, se = NA_real_, p_value = NA_real_))
  }
  model <- feols(
    as.formula(paste(outcome, "~ newparent + i(sex) + i(race) | age + year")),
    data = data,
    weights = as.formula(paste0("~", weight_var))
  )
  extract_newparent(model, label)
}

run_anyparent_reg <- function(data, outcome, label, weight_var = "perwt") {
  y <- data[[outcome]]
  if (nrow(data) == 0 || sum(!is.na(y)) == 0 || dplyr::n_distinct(y[!is.na(y)]) < 2) {
    return(tibble(outcome = label, estimate = NA_real_, se = NA_real_, p_value = NA_real_))
  }
  model <- feols(
    as.formula(paste(outcome, "~ any_parent + i(sex) + i(race) | age + year")),
    data = data,
    weights = as.formula(paste0("~", weight_var))
  )
  extract_parent_coef(model, "any_parent", label)
}

script_dir <- get_script_dir()
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_suffix <- Sys.getenv("MMS_OUTPUT_SUFFIX", "")
out_dir <- file.path(script_dir, paste0("output", out_suffix))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
middle_target <- tolower(Sys.getenv("MMS_MIDDLE_TARGET", ""))
if (middle_target == "") {
  include_middle_in_center <- tolower(Sys.getenv("MMS_INCLUDE_MIDDLE_IN_CENTER", "false")) %in% c("1", "true", "yes", "y")
  middle_target <- if (include_middle_in_center) "center" else "drop"
}
if (!(middle_target %in% c("drop", "center", "periphery"))) {
  stop("MMS_MIDDLE_TARGET must be one of: drop, center, periphery")
}
core_pop_share <- as.numeric(Sys.getenv("MMS_CORE_POP_SHARE", "0.10"))
newparent_eldch_cutoff <- as.numeric(Sys.getenv("MMS_NEWPARENT_ELDCH_CUTOFF", "4"))

lookup_2010_path <- file.path(data_dir, "puma_mms_lookup_2010.csv")
lookup_2020_path <- file.path(data_dir, "puma_mms_lookup_2020.csv")
origin_bridge_path <- file.path(data_dir, "migpuma_mms_origin_bridge.csv")
origin_geo_path <- "code/data/Spatial_aggregate_withmicrodata/raw_data/extract28_origin_geo_2012_2023_age22_45_households.csv"

if (!file.exists(lookup_2010_path) || !file.exists(lookup_2020_path) || !file.exists(origin_bridge_path)) {
  stop("Run build_mms_geography.R and build_migpuma_origin_bridge.R before analyze_mms_fertility_moves.R.")
}

extract_path <- "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
if (!file.exists(extract_path)) {
  stop("Missing microdata extract: ", extract_path)
}
if (!file.exists(origin_geo_path)) {
  stop("Missing extract28 origin supplement: ", origin_geo_path)
}

lookup_2010 <- read_csv(lookup_2010_path, show_col_types = FALSE)
lookup_2020 <- read_csv(lookup_2020_path, show_col_types = FALSE)
origin_bridge <- read_csv(origin_bridge_path, show_col_types = FALSE)
origin_geo <- read_csv(origin_geo_path, show_col_types = FALSE) %>%
  transmute(
    year = as.integer(year),
    sample = as.integer(sample),
    serial = as.numeric(serial),
    pernum = as.integer(pernum),
    migplac1 = as.integer(migplac1),
    migmetro1 = as.integer(migmetro1)
  )

dup_origin_keys <- origin_geo %>%
  count(year, sample, serial, pernum, name = "n") %>%
  filter(n > 1)
if (nrow(dup_origin_keys) > 0) {
  stop("Origin supplement has duplicate person keys.")
}

dest_lookup <- bind_rows(
  lookup_2010 %>%
    transmute(
      lookup_period = "pre2022",
      statefip,
      puma,
      cbsacode,
      dest_mms_location = mms_location
    ),
  lookup_2020 %>%
    transmute(
      lookup_period = "post2021",
      statefip,
      puma,
      cbsacode,
      dest_mms_location = mms_location
    )
)

message("Loading extract27 columns needed for MMS fertility/move analysis...")
df <- read_dta(
  extract_path,
  col_select = c(
    "year", "sample", "serial", "pernum", "statefip", "puma", "migpuma1", "metro", "met2013", "migmet131",
    "age", "sex", "race", "perwt", "gq", "ownershp", "rooms", "rent", "hhincome",
    "nchild", "nchlt5", "eldch", "migrate1", "migrate1d"
  )
) %>%
  transmute(
    year = as.integer(year),
    sample = as.integer(sample),
    serial = as.numeric(serial),
    pernum = as.integer(pernum),
    statefip = as.integer(statefip),
    puma = as.integer(puma),
    migpuma1 = as.integer(migpuma1),
    metro = as.integer(metro),
    met2013 = as.integer(met2013),
    migmet131 = as.integer(migmet131),
    age = as.integer(age),
    sex = as.integer(sex),
    race = as.integer(race),
    perwt = as.numeric(perwt),
    gq = as.integer(gq),
    ownershp = as.integer(ownershp),
    rooms = as.numeric(rooms),
    rent = as.numeric(rent),
    hhincome = as.numeric(hhincome),
    nchild = as.numeric(nchild),
    nchlt5 = as.numeric(nchlt5),
    eldch = as.numeric(eldch),
    migrate1 = as.integer(migrate1),
    migrate1d = as.integer(migrate1d)
  ) %>%
  filter(
    year >= 2012,
    age >= 22,
    age <= 45,
    gq %in% c(1, 2)
  ) %>%
  mutate(
    owner = ownershp == 1,
    renter = ownershp == 2,
    has_youngchild = nchlt5 > 0 & !is.na(nchlt5),
    has_children = nchild > 0 & !is.na(nchild),
    childless = nchild == 0 & !is.na(nchild),
    newparent = as.integer(eldch < newparent_eldch_cutoff & eldch != 99 & nchild > 0),
    any_parent = as.integer(has_children),
    parent_status = case_when(
      childless ~ "Non-Parents",
      newparent == 1 ~ "New Parents",
      has_children ~ "Older Parents",
      TRUE ~ NA_character_
    ),
    parent_compare = case_when(
      childless ~ "Non-Parents",
      newparent == 1 ~ "New Parents",
      TRUE ~ NA_character_
    ),
    parent_compare_all = case_when(
      childless ~ "Non-Parents",
      has_children ~ "All Parents",
      TRUE ~ NA_character_
    ),
    moved1y = migrate1 >= 2 & !is.na(migrate1),
    same_msa = moved1y & met2013 > 0 & migmet131 > 0 & met2013 == migmet131,
    diff_msa = moved1y & met2013 > 0 & migmet131 > 0 & met2013 != migmet131,
    within_state_move = migrate1d %in% c(23, 24, 25),
    rti = ifelse(rent > 0 & hhincome > 0, 12 * rent / hhincome, NA_real_),
    lookup_period = if_else(year <= 2021, "pre2022", "post2021")
  ) %>%
  left_join(origin_geo, by = c("year", "sample", "serial", "pernum"))

df_mms <- df %>%
  left_join(dest_lookup, by = c("lookup_period", "statefip", "puma", "met2013" = "cbsacode")) %>%
  left_join(
    origin_bridge,
    by = c(
      "lookup_period" = "period",
      "migplac1" = "migplac1",
      "migpuma1" = "migpuma1",
      "migmet131" = "cbsacode"
    )
  ) %>%
  mutate(
    dest_mms_location = collapse_mms_location(dest_mms_location, middle_target),
    dest_mms_location = if_else(dest_mms_location %in% c("center", "periphery"), dest_mms_location, NA_character_),
    dest_label = case_when(
      dest_mms_location == "center" ~ "Center",
      dest_mms_location == "periphery" ~ "Periphery",
      TRUE ~ NA_character_
    ),
    center_share = coalesce(center_share, 0),
    periphery_share = coalesce(periphery_share, 0),
    middle_share = coalesce(middle_share, 0),
    center_share = if (middle_target == "center") center_share + middle_share else center_share,
    periphery_share = if (middle_target == "periphery") periphery_share + middle_share else periphery_share,
    middle_share = if (middle_target == "drop") middle_share else 0,
    origin_cp_share = center_share + periphery_share,
    origin_match = !is.na(matched_pop_weight),
    origin_center_weight = perwt * center_share,
    origin_periphery_weight = perwt * periphery_share
  ) %>%
  select(-lookup_period)

origin_match_rate <- mean(!is.na(df_mms$migplac1))
if (origin_match_rate < 0.99) {
  stop(sprintf("Origin supplement merge rate too low: %.3f", origin_match_rate))
}

rm(df)
invisible(gc())

message("Computing MMS fertility/location summaries...")

location_parent_shares <- df_mms %>%
  filter(!is.na(parent_status), !is.na(dest_label)) %>%
  group_by(parent_status, dest_label) %>%
  summarise(pop_weight = sum(perwt, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(share = pop_weight / sum(pop_weight)) %>%
  ungroup()

location_parent_summary <- df_mms %>%
  filter(!is.na(parent_status), !is.na(dest_label)) %>%
  group_by(parent_status) %>%
  summarise(
    weighted_n = sum(perwt, na.rm = TRUE),
    center_share = weighted_mean_safe(dest_mms_location == "center", perwt),
    owner_rate = weighted_mean_safe(owner, perwt),
    has_youngchild_rate = weighted_mean_safe(has_youngchild, perwt),
    mean_rooms = weighted_mean_safe(rooms, perwt),
    mean_rti = weighted_mean_safe(rti, perwt),
    .groups = "drop"
  )

age_fertility_profiles <- df_mms %>%
  filter(!is.na(dest_label)) %>%
  group_by(dest_label, age) %>%
  summarise(
    weighted_n = sum(perwt, na.rm = TRUE),
    has_youngchild_rate = weighted_mean_safe(has_youngchild, perwt),
    has_children_rate = weighted_mean_safe(has_children, perwt),
    newparent_rate = weighted_mean_safe(newparent == 1, perwt),
    owner_rate = weighted_mean_safe(owner, perwt),
    .groups = "drop"
  )

write_csv(location_parent_shares, file.path(out_dir, "mms_location_by_parent_shares.csv"))
write_csv(location_parent_summary, file.path(out_dir, "mms_location_by_parent_summary.csv"))
write_csv(age_fertility_profiles, file.path(out_dir, "mms_age_fertility_profiles.csv"))

message("Computing MMS mover diagnostics...")

movers <- df_mms %>%
  filter(
    moved1y,
    met2013 > 0,
    migmet131 > 0,
    !is.na(parent_compare),
    !is.na(dest_label)
  ) %>%
  mutate(
    center_dest = as.integer(dest_mms_location == "center"),
    move_type = case_when(
      same_msa & dest_mms_location == "center" ~ "Within CBSA -> Center",
      same_msa & dest_mms_location == "periphery" ~ "Within CBSA -> Periphery",
      diff_msa & dest_mms_location == "center" ~ "Across CBSA -> Center",
      diff_msa & dest_mms_location == "periphery" ~ "Across CBSA -> Periphery",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(move_type))

move_type_summary <- movers %>%
  group_by(parent_compare) %>%
  summarise(
    weighted_n = sum(perwt, na.rm = TRUE),
    within_cbsa_share = weighted_mean_safe(same_msa, perwt),
    across_cbsa_share = weighted_mean_safe(diff_msa, perwt),
    center_dest_share_all = weighted_mean_safe(dest_mms_location == "center", perwt),
    center_dest_share_within = weighted_mean_safe(dest_mms_location[same_msa] == "center", perwt[same_msa]),
    owner_rate = weighted_mean_safe(owner, perwt),
    .groups = "drop"
  )

four_way_shares <- movers %>%
  group_by(parent_compare, move_type) %>%
  summarise(pop_weight = sum(perwt, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(share = pop_weight / sum(pop_weight)) %>%
  ungroup()

ownership_by_move_type <- movers %>%
  group_by(parent_compare, move_type) %>%
  summarise(
    weighted_n = sum(perwt, na.rm = TRUE),
    owner_rate = weighted_mean_safe(owner, perwt),
    mean_rooms = weighted_mean_safe(rooms, perwt),
    mean_rti = weighted_mean_safe(rti, perwt),
    .groups = "drop"
  )

within_movers <- movers %>%
  filter(same_msa)

within_movers_origin <- within_movers %>%
  filter(origin_cp_share > 0)

move_regressions <- bind_rows(
  run_newparent_reg(movers, "diff_msa", "Across-CBSA mover"),
  run_newparent_reg(within_movers, "center_dest", "Center destination among within-CBSA movers")
) %>%
  mutate(
    estimate_pp = 100 * estimate,
    se_pp = 100 * se
  )

four_way_reg_data <- movers %>%
  mutate(
    within_center = as.integer(move_type == "Within CBSA -> Center"),
    within_periphery = as.integer(move_type == "Within CBSA -> Periphery"),
    across_center = as.integer(move_type == "Across CBSA -> Center"),
    across_periphery = as.integer(move_type == "Across CBSA -> Periphery")
  )

four_way_regressions <- bind_rows(
  run_newparent_reg(four_way_reg_data, "within_center", "Within CBSA -> Center"),
  run_newparent_reg(four_way_reg_data, "within_periphery", "Within CBSA -> Periphery"),
  run_newparent_reg(four_way_reg_data, "across_center", "Across CBSA -> Center"),
  run_newparent_reg(four_way_reg_data, "across_periphery", "Across CBSA -> Periphery")
) %>%
  mutate(
    estimate_pp = 100 * estimate,
    se_pp = 100 * se
  )

origin_transition_cells <- bind_rows(
  within_movers_origin %>%
    transmute(
      parent_compare,
      origin_label = "Center",
      dest_label,
      flow_weight = origin_center_weight
    ),
  within_movers_origin %>%
    transmute(
      parent_compare,
      origin_label = "Periphery",
      dest_label,
      flow_weight = origin_periphery_weight
    )
) %>%
  filter(flow_weight > 0, !is.na(dest_label))

origin_transition_summary <- origin_transition_cells %>%
  group_by(parent_compare, origin_label, dest_label) %>%
  summarise(pop_weight = sum(flow_weight, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(share = pop_weight / sum(pop_weight)) %>%
  ungroup()

origin_weight_summary <- within_movers %>%
  group_by(parent_compare) %>%
  summarise(
    weighted_n = sum(perwt, na.rm = TRUE),
    origin_bridge_share = weighted_mean_safe(origin_match, perwt),
    origin_cp_share = weighted_mean_safe(origin_cp_share, perwt),
    center_origin_mass = sum(origin_center_weight, na.rm = TRUE),
    periphery_origin_mass = sum(origin_periphery_weight, na.rm = TRUE),
    .groups = "drop"
  )

origin_regressions <- bind_rows(
  run_newparent_reg(
    within_movers_origin %>% mutate(reg_weight = origin_center_weight),
    "center_dest",
    "Center destination among within-CBSA movers, center-origin mass",
    weight_var = "reg_weight"
  ),
  run_newparent_reg(
    within_movers_origin %>% mutate(reg_weight = origin_periphery_weight),
    "center_dest",
    "Center destination among within-CBSA movers, periphery-origin mass",
    weight_var = "reg_weight"
  )
) %>%
  mutate(
    estimate_pp = 100 * estimate,
    se_pp = 100 * se
  )

# ---------- Parallel "Non-Parents vs All Parents" block ----------
movers_all <- df_mms %>%
  filter(
    moved1y,
    met2013 > 0,
    migmet131 > 0,
    !is.na(parent_compare_all),
    !is.na(dest_label)
  ) %>%
  mutate(
    center_dest = as.integer(dest_mms_location == "center"),
    move_type = case_when(
      same_msa & dest_mms_location == "center" ~ "Within CBSA -> Center",
      same_msa & dest_mms_location == "periphery" ~ "Within CBSA -> Periphery",
      diff_msa & dest_mms_location == "center" ~ "Across CBSA -> Center",
      diff_msa & dest_mms_location == "periphery" ~ "Across CBSA -> Periphery",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(move_type))

move_type_summary_all <- movers_all %>%
  group_by(parent_compare_all) %>%
  summarise(
    weighted_n = sum(perwt, na.rm = TRUE),
    within_cbsa_share = weighted_mean_safe(same_msa, perwt),
    across_cbsa_share = weighted_mean_safe(diff_msa, perwt),
    center_dest_share_all = weighted_mean_safe(dest_mms_location == "center", perwt),
    center_dest_share_within = weighted_mean_safe(dest_mms_location[same_msa] == "center", perwt[same_msa]),
    owner_rate = weighted_mean_safe(owner, perwt),
    .groups = "drop"
  )

four_way_shares_all <- movers_all %>%
  group_by(parent_compare_all, move_type) %>%
  summarise(pop_weight = sum(perwt, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(share = pop_weight / sum(pop_weight)) %>%
  ungroup()

within_movers_all <- movers_all %>% filter(same_msa)
within_movers_origin_all <- within_movers_all %>% filter(origin_cp_share > 0)

four_way_reg_data_all <- movers_all %>%
  mutate(
    within_center = as.integer(move_type == "Within CBSA -> Center"),
    within_periphery = as.integer(move_type == "Within CBSA -> Periphery"),
    across_center = as.integer(move_type == "Across CBSA -> Center"),
    across_periphery = as.integer(move_type == "Across CBSA -> Periphery")
  )

move_regressions_all <- bind_rows(
  run_anyparent_reg(movers_all, "diff_msa", "Across-CBSA mover"),
  run_anyparent_reg(within_movers_all, "center_dest", "Center destination among within-CBSA movers")
) %>% mutate(estimate_pp = 100 * estimate, se_pp = 100 * se)

four_way_regressions_all <- bind_rows(
  run_anyparent_reg(four_way_reg_data_all, "within_center", "Within CBSA -> Center"),
  run_anyparent_reg(four_way_reg_data_all, "within_periphery", "Within CBSA -> Periphery"),
  run_anyparent_reg(four_way_reg_data_all, "across_center", "Across CBSA -> Center"),
  run_anyparent_reg(four_way_reg_data_all, "across_periphery", "Across CBSA -> Periphery")
) %>% mutate(estimate_pp = 100 * estimate, se_pp = 100 * se)

origin_transition_cells_all <- bind_rows(
  within_movers_origin_all %>% transmute(parent_compare_all, origin_label = "Center", dest_label, flow_weight = origin_center_weight),
  within_movers_origin_all %>% transmute(parent_compare_all, origin_label = "Periphery", dest_label, flow_weight = origin_periphery_weight)
) %>% filter(flow_weight > 0, !is.na(dest_label))

origin_transition_summary_all <- origin_transition_cells_all %>%
  group_by(parent_compare_all, origin_label, dest_label) %>%
  summarise(pop_weight = sum(flow_weight, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(share = pop_weight / sum(pop_weight)) %>%
  ungroup()

origin_regressions_all <- bind_rows(
  run_anyparent_reg(
    within_movers_origin_all %>% mutate(reg_weight = origin_center_weight),
    "center_dest",
    "Center destination among within-CBSA movers, center-origin mass",
    weight_var = "reg_weight"
  ),
  run_anyparent_reg(
    within_movers_origin_all %>% mutate(reg_weight = origin_periphery_weight),
    "center_dest",
    "Center destination among within-CBSA movers, periphery-origin mass",
    weight_var = "reg_weight"
  )
) %>% mutate(estimate_pp = 100 * estimate, se_pp = 100 * se)

coverage_summary <- tibble(
  metric = c(
    "All households with MMS-matched destination",
    "All households with extract28 origin supplement",
    "Movers with MMS-matched destination and previous CBSA",
    "Within-CBSA movers with bridge-matched origin mass"
  ),
  weighted_n = c(
    sum(df_mms$perwt[!is.na(df_mms$dest_label)], na.rm = TRUE),
    sum(df_mms$perwt[!is.na(df_mms$migplac1)], na.rm = TRUE),
    sum(movers$perwt, na.rm = TRUE),
    sum(within_movers_origin$perwt, na.rm = TRUE)
  ),
  raw_n = c(
    sum(!is.na(df_mms$dest_label)),
    sum(!is.na(df_mms$migplac1)),
    nrow(movers),
    nrow(within_movers_origin)
  )
)

write_csv(coverage_summary, file.path(out_dir, "mms_move_coverage_summary.csv"))
write_csv(move_type_summary, file.path(out_dir, "mms_move_type_summary.csv"))
write_csv(four_way_shares, file.path(out_dir, "mms_four_way_move_shares.csv"))
write_csv(ownership_by_move_type, file.path(out_dir, "mms_ownership_by_move_type.csv"))
write_csv(move_regressions, file.path(out_dir, "mms_move_regressions.csv"))
write_csv(four_way_regressions, file.path(out_dir, "mms_four_way_regressions.csv"))
write_csv(origin_transition_summary, file.path(out_dir, "mms_origin_transition_summary.csv"))
write_csv(origin_weight_summary, file.path(out_dir, "mms_origin_weight_summary.csv"))
write_csv(origin_regressions, file.path(out_dir, "mms_origin_regressions.csv"))

# Parallel all-parent outputs
write_csv(move_type_summary_all, file.path(out_dir, "mms_move_type_summary_allparent.csv"))
write_csv(four_way_shares_all, file.path(out_dir, "mms_four_way_move_shares_allparent.csv"))
write_csv(move_regressions_all, file.path(out_dir, "mms_move_regressions_allparent.csv"))
write_csv(four_way_regressions_all, file.path(out_dir, "mms_four_way_regressions_allparent.csv"))
write_csv(origin_transition_summary_all, file.path(out_dir, "mms_origin_transition_summary_allparent.csv"))
write_csv(origin_regressions_all, file.path(out_dir, "mms_origin_regressions_allparent.csv"))

location_plot <- location_parent_shares %>%
  mutate(
    parent_status = factor(parent_status, levels = c("Non-Parents", "New Parents", "Older Parents")),
    dest_label = factor(dest_label, levels = c("Center", "Periphery"))
  ) %>%
  ggplot(aes(x = parent_status, y = share, fill = dest_label)) +
  geom_col(position = "fill", width = 0.68) +
  geom_text(
    aes(label = safe_percent(share)),
    position = position_fill(vjust = 0.5),
    color = "white",
    size = 4
  ) +
  scale_fill_manual(values = c("Center" = "#1b4965", "Periphery" = "#ca6702")) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "",
    y = "Share",
    fill = "",
    title = "Current MMS Location by Parenthood Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "mms_location_by_parent.png"), location_plot, width = 8, height = 5, dpi = 150)

move_bar_data <- move_type_summary %>%
  select(parent_compare, within_cbsa_share, across_cbsa_share) %>%
  pivot_longer(cols = c(within_cbsa_share, across_cbsa_share), names_to = "move_margin", values_to = "share") %>%
  mutate(
    move_margin = recode(move_margin, within_cbsa_share = "Within CBSA", across_cbsa_share = "Across CBSA"),
    parent_compare = factor(parent_compare, levels = c("Non-Parents", "New Parents"))
  )

move_bar_plot <- ggplot(move_bar_data, aes(x = parent_compare, y = share, fill = move_margin)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(
    aes(label = safe_percent(share)),
    position = position_dodge(width = 0.6),
    vjust = -0.35,
    size = 4
  ) +
  scale_fill_manual(values = c("Within CBSA" = "#1b4965", "Across CBSA" = "#9c6644")) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "",
    y = "Share of movers",
    fill = "",
    title = "Move Type by Parenthood Status Under MMS Geography"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "mms_bar_move_type_by_parent.png"), move_bar_plot, width = 7.5, height = 5, dpi = 150)

four_way_plot_data <- four_way_shares %>%
  mutate(
    parent_compare = factor(parent_compare, levels = c("Non-Parents", "New Parents")),
    move_type = factor(
      move_type,
      levels = c(
        "Within CBSA -> Center",
        "Within CBSA -> Periphery",
        "Across CBSA -> Center",
        "Across CBSA -> Periphery"
      )
    )
  )

four_way_plot <- ggplot(four_way_plot_data, aes(x = move_type, y = share, fill = parent_compare)) +
  geom_col(position = "dodge", width = 0.68) +
  geom_text(
    aes(label = safe_percent(share)),
    position = position_dodge(width = 0.68),
    vjust = -0.35,
    size = 3.5
  ) +
  scale_fill_manual(values = c("Non-Parents" = "#1b4965", "New Parents" = "#ae2012")) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "",
    y = "Share of movers",
    fill = "",
    title = "MMS Destination Decomposition for Movers"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )
ggsave(file.path(out_dir, "mms_grouped_4way_decomposition.png"), four_way_plot, width = 9, height = 5.5, dpi = 150)

fertility_profile_plot <- age_fertility_profiles %>%
  mutate(dest_label = factor(dest_label, levels = c("Center", "Periphery"))) %>%
  ggplot(aes(x = age, y = has_youngchild_rate, color = dest_label)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Center" = "#1b4965", "Periphery" = "#ca6702")) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Age",
    y = "Share with child < 5",
    color = "",
    title = "Fertility Gradient Under MMS Geography"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "mms_age_fertility_gradient.png"), fertility_profile_plot, width = 8, height = 5, dpi = 150)

origin_flow_plot_data <- origin_transition_summary %>%
  mutate(
    parent_compare = factor(parent_compare, levels = c("Non-Parents", "New Parents")),
    origin_label = factor(origin_label, levels = c("Center", "Periphery")),
    dest_label = factor(dest_label, levels = c("Center", "Periphery"))
  )

origin_flow_plot <- ggplot(origin_flow_plot_data, aes(x = dest_label, y = share, fill = parent_compare)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(
    aes(label = safe_percent(share)),
    position = position_dodge(width = 0.65),
    vjust = -0.35,
    size = 3.4
  ) +
  facet_wrap(~origin_label) +
  scale_fill_manual(values = c("Non-Parents" = "#1b4965", "New Parents" = "#ae2012")) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "",
    y = "Destination share within origin group",
    fill = "",
    title = "Within-CBSA MMS Destination by Origin Location"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "mms_within_cbsa_origin_destination.png"), origin_flow_plot, width = 8.5, height = 5.5, dpi = 150)

# ---------- Parallel all-parent plots ----------
four_way_plot_data_all <- four_way_shares_all %>%
  mutate(
    parent_compare_all = factor(parent_compare_all, levels = c("Non-Parents", "All Parents")),
    move_type = factor(
      move_type,
      levels = c(
        "Within CBSA -> Center",
        "Within CBSA -> Periphery",
        "Across CBSA -> Center",
        "Across CBSA -> Periphery"
      )
    )
  )

four_way_plot_all <- ggplot(four_way_plot_data_all, aes(x = move_type, y = share, fill = parent_compare_all)) +
  geom_col(position = "dodge", width = 0.68) +
  geom_text(
    aes(label = safe_percent(share)),
    position = position_dodge(width = 0.68),
    vjust = -0.35,
    size = 3.5
  ) +
  scale_fill_manual(values = c("Non-Parents" = "#1b4965", "All Parents" = "#ae2012")) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "",
    y = "Share of movers",
    fill = "",
    title = "MMS Destination Decomposition for Movers"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )
ggsave(file.path(out_dir, "mms_grouped_4way_decomposition_allparent.png"), four_way_plot_all, width = 9, height = 5.5, dpi = 150)

origin_flow_plot_data_all <- origin_transition_summary_all %>%
  mutate(
    parent_compare_all = factor(parent_compare_all, levels = c("Non-Parents", "All Parents")),
    origin_label = factor(origin_label, levels = c("Center", "Periphery")),
    dest_label = factor(dest_label, levels = c("Center", "Periphery"))
  )

origin_flow_plot_all <- ggplot(origin_flow_plot_data_all, aes(x = dest_label, y = share, fill = parent_compare_all)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(
    aes(label = safe_percent(share)),
    position = position_dodge(width = 0.65),
    vjust = -0.35,
    size = 3.4
  ) +
  facet_wrap(~origin_label) +
  scale_fill_manual(values = c("Non-Parents" = "#1b4965", "All Parents" = "#ae2012")) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "",
    y = "Destination share within origin group",
    fill = "",
    title = "Within-CBSA MMS Destination by Origin Location"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "mms_within_cbsa_origin_destination_allparent.png"), origin_flow_plot_all, width = 8.5, height = 5.5, dpi = 150)

np_move <- move_type_summary %>% filter(parent_compare == "Non-Parents")
par_move <- move_type_summary %>% filter(parent_compare == "New Parents")
np_origin <- origin_transition_summary %>% filter(parent_compare == "Non-Parents")
par_origin <- origin_transition_summary %>% filter(parent_compare == "New Parents")

summary_lines <- c(
  "# MMS Fertility and Move Summary",
  "",
  "## Coverage",
  sprintf("- Core tract population target share: %.2f.", core_pop_share),
  sprintf("- Middle target: %s.", middle_target),
  sprintf("- New-parent cutoff: eldest child < %.0f.", newparent_eldch_cutoff),
  sprintf(
    "- Current-location MMS match: %.0f weighted observations across %d raw records.",
    coverage_summary$weighted_n[coverage_summary$metric == "All households with MMS-matched destination"],
    coverage_summary$raw_n[coverage_summary$metric == "All households with MMS-matched destination"]
  ),
  sprintf(
    "- Movers with current and previous CBSA plus MMS destination: %.0f weighted observations across %d raw records.",
    coverage_summary$weighted_n[coverage_summary$metric == "Movers with MMS-matched destination and previous CBSA"],
    coverage_summary$raw_n[coverage_summary$metric == "Movers with MMS-matched destination and previous CBSA"]
  ),
  sprintf(
    "- Within-CBSA movers with bridge-matched center/periphery origin mass: %.0f weighted observations across %d raw records.",
    coverage_summary$weighted_n[coverage_summary$metric == "Within-CBSA movers with bridge-matched origin mass"],
    coverage_summary$raw_n[coverage_summary$metric == "Within-CBSA movers with bridge-matched origin mass"]
  ),
  "",
  "## Current location by parent status",
  sprintf(
    "- Center share: %s for non-parents, %s for new parents, %s for older parents.",
    safe_percent(location_parent_summary$center_share[location_parent_summary$parent_status == "Non-Parents"]),
    safe_percent(location_parent_summary$center_share[location_parent_summary$parent_status == "New Parents"]),
    safe_percent(location_parent_summary$center_share[location_parent_summary$parent_status == "Older Parents"])
  ),
  sprintf(
    "- Ownership rate: %s for non-parents, %s for new parents, %s for older parents.",
    safe_percent(location_parent_summary$owner_rate[location_parent_summary$parent_status == "Non-Parents"]),
    safe_percent(location_parent_summary$owner_rate[location_parent_summary$parent_status == "New Parents"]),
    safe_percent(location_parent_summary$owner_rate[location_parent_summary$parent_status == "Older Parents"])
  ),
  "",
  "## Movers: non-parents vs new parents",
  sprintf(
    "- Within-CBSA share: %s for non-parents vs %s for new parents.",
    safe_percent(np_move$within_cbsa_share),
    safe_percent(par_move$within_cbsa_share)
  ),
  sprintf(
    "- Among within-CBSA movers, destination center share: %s for non-parents vs %s for new parents.",
    safe_percent(np_move$center_dest_share_within),
    safe_percent(par_move$center_dest_share_within)
  ),
  sprintf(
    "- New-parent coefficient on across-CBSA moving: %.2f pp (se %.2f pp, p = %.4f).",
    move_regressions$estimate_pp[move_regressions$outcome == "Across-CBSA mover"],
    move_regressions$se_pp[move_regressions$outcome == "Across-CBSA mover"],
    move_regressions$p_value[move_regressions$outcome == "Across-CBSA mover"]
  ),
  sprintf(
    "- New-parent coefficient on destination center among within-CBSA movers: %.2f pp (se %.2f pp, p = %.4f).",
    move_regressions$estimate_pp[move_regressions$outcome == "Center destination among within-CBSA movers"],
    move_regressions$se_pp[move_regressions$outcome == "Center destination among within-CBSA movers"],
    move_regressions$p_value[move_regressions$outcome == "Center destination among within-CBSA movers"]
  ),
  "",
  "## Local origin-destination refinement",
  sprintf(
    "- Average bridge-classified center/periphery origin mass among within-CBSA movers: %s for non-parents vs %s for new parents.",
    safe_percent(origin_weight_summary$origin_cp_share[origin_weight_summary$parent_compare == "Non-Parents"]),
    safe_percent(origin_weight_summary$origin_cp_share[origin_weight_summary$parent_compare == "New Parents"])
  ),
  sprintf(
    "- From center-origin mass, destination center share is %s for non-parents vs %s for new parents.",
    safe_percent(np_origin$share[np_origin$origin_label == "Center" & np_origin$dest_label == "Center"]),
    safe_percent(par_origin$share[par_origin$origin_label == "Center" & par_origin$dest_label == "Center"])
  ),
  sprintf(
    "- From periphery-origin mass, destination center share is %s for non-parents vs %s for new parents.",
    safe_percent(np_origin$share[np_origin$origin_label == "Periphery" & np_origin$dest_label == "Center"]),
    safe_percent(par_origin$share[par_origin$origin_label == "Periphery" & par_origin$dest_label == "Center"])
  ),
  sprintf(
    "- New-parent coefficient on destination center within center-origin mass: %.2f pp (se %.2f pp, p = %.4f).",
    origin_regressions$estimate_pp[origin_regressions$outcome == "Center destination among within-CBSA movers, center-origin mass"],
    origin_regressions$se_pp[origin_regressions$outcome == "Center destination among within-CBSA movers, center-origin mass"],
    origin_regressions$p_value[origin_regressions$outcome == "Center destination among within-CBSA movers, center-origin mass"]
  ),
  sprintf(
    "- New-parent coefficient on destination center within periphery-origin mass: %.2f pp (se %.2f pp, p = %.4f).",
    origin_regressions$estimate_pp[origin_regressions$outcome == "Center destination among within-CBSA movers, periphery-origin mass"],
    origin_regressions$se_pp[origin_regressions$outcome == "Center destination among within-CBSA movers, periphery-origin mass"],
    origin_regressions$p_value[origin_regressions$outcome == "Center destination among within-CBSA movers, periphery-origin mass"]
  ),
  "",
  "## Files",
  "- `mms_location_by_parent.png`",
  "- `mms_bar_move_type_by_parent.png`",
  "- `mms_grouped_4way_decomposition.png`",
  "- `mms_age_fertility_gradient.png`",
  "- `mms_within_cbsa_origin_destination.png`"
)

writeLines(summary_lines, file.path(out_dir, "mms_fertility_moves_summary.md"))

message("MMS fertility/move analysis complete.")
