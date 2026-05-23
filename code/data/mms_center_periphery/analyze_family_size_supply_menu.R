#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
  library(tidyr)
  library(ggplot2)
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
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

weighted_quantile_safe <- function(x, w, prob = 0.5) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= prob)[1]]
}

coef_table <- function(model, model_name) {
  ct <- as.data.frame(coeftable(model))
  ct$term <- rownames(ct)
  rownames(ct) <- NULL
  names(ct) <- c("estimate", "se", "t_stat", "p_value", "term")
  ct %>%
    transmute(
      model = model_name,
      term,
      estimate,
      se,
      t_stat,
      p_value
    )
}

script_dir <- get_script_dir()
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_dir <- file.path(script_dir, Sys.getenv("FAMILY_SIZE_OUTPUT_DIR", "output_family_size_supply"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

middle_target <- tolower(Sys.getenv("MMS_MIDDLE_TARGET", ""))
if (middle_target == "") {
  include_middle_in_center <- tolower(Sys.getenv("MMS_INCLUDE_MIDDLE_IN_CENTER", "true")) %in% c("1", "true", "yes", "y")
  middle_target <- if (include_middle_in_center) "center" else "drop"
}
if (!(middle_target %in% c("drop", "center", "periphery"))) {
  stop("MMS_MIDDLE_TARGET must be one of: drop, center, periphery")
}

lookup_2010_path <- file.path(data_dir, "puma_mms_lookup_2010.csv")
lookup_2020_path <- file.path(data_dir, "puma_mms_lookup_2020.csv")
extract_path <- "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"

if (!file.exists(lookup_2010_path) || !file.exists(lookup_2020_path)) {
  stop("Missing MMS lookup files. Run build_mms_geography.R first.")
}
if (!file.exists(extract_path)) {
  stop("Missing microdata extract: ", extract_path)
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

message("Loading ACS/IPUMS household-head sample for family-size supply packet...")
df <- read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "met2013", "gq", "pernum", "hhwt", "perwt",
    "ownershp", "rooms", "bedrooms", "rent", "hhincome", "nchild", "nchlt5",
    "eldch", "yngch", "age", "sex", "race", "educ"
  )
) %>%
  transmute(
    year = as.integer(year),
    statefip = as.integer(statefip),
    puma = as.integer(puma),
    met2013 = as.integer(met2013),
    gq = as.integer(gq),
    pernum = as.integer(pernum),
    hhwt = as.numeric(hhwt),
    perwt = as.numeric(perwt),
    ownershp = as.integer(ownershp),
    rooms = as.numeric(rooms),
    bedrooms = as.numeric(bedrooms),
    rent = as.numeric(rent),
    hhincome = as.numeric(hhincome),
    nchild = as.numeric(nchild),
    nchlt5 = as.numeric(nchlt5),
    eldch = as.numeric(eldch),
    yngch = as.numeric(yngch),
    age = as.integer(age),
    sex = as.integer(sex),
    race = as.integer(race),
    educ = as.integer(educ),
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
    mms_location = if_else(mms_location %in% c("center", "periphery"), mms_location, NA_character_),
    center = as.integer(mms_location == "center"),
    owner = ownershp == 1,
    renter = ownershp == 2,
    tenure = case_when(owner ~ "owner", renter ~ "renter", TRUE ~ "other"),
    rent_monthly = if_else(renter & rent > 0, rent, NA_real_),
    rent_per_room = rent_monthly / rooms,
    size_bin = case_when(
      rooms <= 4 ~ "S_1_4",
      rooms <= 6 ~ "M_5_6",
      rooms >= 7 ~ "L_7plus",
      TRUE ~ NA_character_
    ),
    family_capable = size_bin %in% c("M_5_6", "L_7plus"),
    moderate_family = size_bin == "M_5_6",
    has_children = nchild > 0 & !is.na(nchild),
    childless = nchild == 0 & !is.na(nchild),
    has_child_u18 = nchild > 0 & !is.na(yngch) & yngch != 99 & yngch < 18,
    newparent = eldch < 4 & eldch != 99 & nchild > 0,
    parent_status = case_when(
      has_child_u18 ~ "parent_u18",
      childless ~ "childless",
      TRUE ~ NA_character_
    ),
    young_sample = age >= 22 & age <= 45,
    log_hhincome = if_else(hhincome > 0, log(pmax(hhincome, 1)), NA_real_),
    college = as.integer(educ >= 10),
    sex_factor = factor(sex),
    race_factor = factor(race),
    tenure_factor = factor(tenure),
    income_tercile = NA_integer_
  ) %>%
  filter(!is.na(mms_location), !is.na(size_bin))

if (nrow(df) == 0) {
  stop("No matched MMS household-head records.")
}

df <- df %>%
  group_by(year) %>%
  mutate(
    income_tercile = if_else(
      !is.na(log_hhincome),
      as.integer(ntile(log_hhincome, 3)),
      NA_integer_
    ),
    log_hhincome_std = if_else(
      !is.na(log_hhincome),
      as.numeric(scale(log_hhincome)),
      NA_real_
    )
  ) %>%
  ungroup()

message("Computing stock, rent-premium, sorting, and regression outputs...")

cell_stock <- df %>%
  group_by(met2013, cbsatitle, mms_location, size_bin, tenure) %>%
  summarise(
    n_hh = n(),
    hh_weight = sum(hhwt, na.rm = TRUE),
    mean_rooms = weighted_mean_safe(rooms, hhwt),
    mean_rent = weighted_mean_safe(rent_monthly, hhwt),
    median_rent = weighted_quantile_safe(rent_monthly, hhwt, 0.5),
    mean_rent_per_room = weighted_mean_safe(rent_per_room, hhwt),
    median_rent_per_room = weighted_quantile_safe(rent_per_room, hhwt, 0.5),
    .groups = "drop"
  ) %>%
  group_by(met2013, mms_location) %>%
  mutate(share_within_location = hh_weight / sum(hh_weight, na.rm = TRUE)) %>%
  ungroup()

location_size_summary <- df %>%
  group_by(mms_location, size_bin, tenure) %>%
  summarise(
    n_hh = n(),
    hh_weight = sum(hhwt, na.rm = TRUE),
    mean_rooms = weighted_mean_safe(rooms, hhwt),
    mean_rent = weighted_mean_safe(rent_monthly, hhwt),
    median_rent = weighted_quantile_safe(rent_monthly, hhwt, 0.5),
    mean_rent_per_room = weighted_mean_safe(rent_per_room, hhwt),
    median_rent_per_room = weighted_quantile_safe(rent_per_room, hhwt, 0.5),
    .groups = "drop"
  ) %>%
  group_by(mms_location) %>%
  mutate(share_within_location = hh_weight / sum(hh_weight, na.rm = TRUE)) %>%
  ungroup()

renter_price_cells <- df %>%
  filter(renter, !is.na(rent_monthly), rent_monthly > 0) %>%
  group_by(met2013, cbsatitle, mms_location, size_bin) %>%
  summarise(
    renter_weight = sum(hhwt, na.rm = TRUE),
    mean_rent = weighted_mean_safe(rent_monthly, hhwt),
    median_rent = weighted_quantile_safe(rent_monthly, hhwt, 0.5),
    mean_rent_per_room = weighted_mean_safe(rent_per_room, hhwt),
    median_rent_per_room = weighted_quantile_safe(rent_per_room, hhwt, 0.5),
    .groups = "drop"
  )

metro_base <- df %>%
  distinct(met2013, cbsatitle)

price_cell <- function(loc, bin, prefix) {
  renter_price_cells %>%
    filter(mms_location == loc, size_bin == bin) %>%
    transmute(
      met2013,
      cbsatitle,
      !!paste0(prefix, "_renter_weight") := renter_weight,
      !!paste0(prefix, "_median_rent") := median_rent,
      !!paste0(prefix, "_mean_rent") := mean_rent,
      !!paste0(prefix, "_median_rent_per_room") := median_rent_per_room,
      !!paste0(prefix, "_mean_rent_per_room") := mean_rent_per_room
    )
}

central_supply <- cell_stock %>%
  filter(mms_location == "center", tenure %in% c("owner", "renter")) %>%
  group_by(met2013, cbsatitle) %>%
  summarise(
    center_stock_s = sum(if_else(size_bin == "S_1_4", hh_weight, 0), na.rm = TRUE),
    center_stock_m = sum(if_else(size_bin == "M_5_6", hh_weight, 0), na.rm = TRUE),
    center_stock_l = sum(if_else(size_bin == "L_7plus", hh_weight, 0), na.rm = TRUE),
    center_renter_s = sum(if_else(size_bin == "S_1_4" & tenure == "renter", hh_weight, 0), na.rm = TRUE),
    center_renter_m = sum(if_else(size_bin == "M_5_6" & tenure == "renter", hh_weight, 0), na.rm = TRUE),
    center_renter_l = sum(if_else(size_bin == "L_7plus" & tenure == "renter", hh_weight, 0), na.rm = TRUE),
    center_m_total = sum(if_else(size_bin == "M_5_6", hh_weight, 0), na.rm = TRUE),
    .groups = "drop"
  )

location_supply <- cell_stock %>%
  filter(tenure %in% c("owner", "renter")) %>%
  group_by(met2013, cbsatitle, mms_location) %>%
  summarise(
    stock_s = sum(if_else(size_bin == "S_1_4", hh_weight, 0), na.rm = TRUE),
    stock_m = sum(if_else(size_bin == "M_5_6", hh_weight, 0), na.rm = TRUE),
    stock_l = sum(if_else(size_bin == "L_7plus", hh_weight, 0), na.rm = TRUE),
    stock_family = stock_m + stock_l,
    stock_total = stock_s + stock_family,
    family_stock_share = stock_family / stock_total,
    small_to_family_ratio = stock_s / stock_family,
    .groups = "drop"
  ) %>%
  select(met2013, cbsatitle, mms_location, stock_s, stock_m, stock_l, stock_family, stock_total, family_stock_share, small_to_family_ratio) %>%
  pivot_wider(
    names_from = mms_location,
    values_from = c(stock_s, stock_m, stock_l, stock_family, stock_total, family_stock_share, small_to_family_ratio)
  )

metro_premia <- metro_base %>%
  left_join(price_cell("center", "S_1_4", "center_s"), by = c("met2013", "cbsatitle")) %>%
  left_join(price_cell("periphery", "S_1_4", "periphery_s"), by = c("met2013", "cbsatitle")) %>%
  left_join(price_cell("center", "M_5_6", "center_m"), by = c("met2013", "cbsatitle")) %>%
  left_join(price_cell("periphery", "M_5_6", "periphery_m"), by = c("met2013", "cbsatitle")) %>%
  left_join(price_cell("center", "L_7plus", "center_l"), by = c("met2013", "cbsatitle")) %>%
  left_join(price_cell("periphery", "L_7plus", "periphery_l"), by = c("met2013", "cbsatitle")) %>%
  left_join(central_supply, by = c("met2013", "cbsatitle")) %>%
  left_join(location_supply, by = c("met2013", "cbsatitle")) %>%
  mutate(
    log_center_small = log(center_s_median_rent),
    log_periphery_small = log(periphery_s_median_rent),
    log_center_m = log(center_m_median_rent),
    log_periphery_m = log(periphery_m_median_rent),
    log_center_l = log(center_l_median_rent),
    log_periphery_l = log(periphery_l_median_rent),
    family_size_central_premium_m = (log_center_m - log_periphery_m) - (log_center_small - log_periphery_small),
    family_size_central_premium_l = (log_center_l - log_periphery_l) - (log_center_small - log_periphery_small),
    unit_log_center_small = log(center_s_median_rent_per_room),
    unit_log_periphery_small = log(periphery_s_median_rent_per_room),
    unit_log_center_m = log(center_m_median_rent_per_room),
    unit_log_periphery_m = log(periphery_m_median_rent_per_room),
    unit_family_size_central_premium_m = (unit_log_center_m - unit_log_periphery_m) - (unit_log_center_small - unit_log_periphery_small),
    central_family_rental_weight = center_renter_m + center_renter_l,
    central_total_rental_weight = center_renter_s + center_renter_m + center_renter_l,
    central_m_rental_share_of_rentals = center_renter_m / central_total_rental_weight,
    central_family_rental_share_of_rentals = central_family_rental_weight / central_total_rental_weight,
    central_m_units_rental_share = center_renter_m / center_m_total,
    central_family_stock_share = (center_stock_m + center_stock_l) /
      (center_stock_s + center_stock_m + center_stock_l),
    family_space_central_scarcity = log(small_to_family_ratio_center) - log(small_to_family_ratio_periphery),
    family_space_central_access = family_stock_share_center - family_stock_share_periphery,
    center_stock_weight = center_stock_s + center_stock_m + center_stock_l
  )

parent_sorting <- df %>%
  filter(young_sample, parent_status %in% c("childless", "parent_u18")) %>%
  group_by(parent_status, income_tercile, mms_location, tenure) %>%
  summarise(
    hh_weight = sum(hhwt, na.rm = TRUE),
    mean_rooms = weighted_mean_safe(rooms, hhwt),
    mean_rent = weighted_mean_safe(rent_monthly, hhwt),
    .groups = "drop"
  ) %>%
  group_by(parent_status, income_tercile) %>%
  mutate(
    share = hh_weight / sum(hh_weight, na.rm = TRUE)
  ) %>%
  ungroup()

parent_centrality_metro <- df %>%
  filter(young_sample, parent_status %in% c("childless", "parent_u18")) %>%
  group_by(met2013, cbsatitle, parent_status) %>%
  summarise(
    hh_weight = sum(hhwt, na.rm = TRUE),
    center_share = weighted_mean_safe(center, hhwt),
    owner_share = weighted_mean_safe(owner, hhwt),
    mean_rooms = weighted_mean_safe(rooms, hhwt),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = parent_status,
    values_from = c(hh_weight, center_share, owner_share, mean_rooms)
  ) %>%
  mutate(
    parent_center_gap = center_share_parent_u18 - center_share_childless,
    childless_parent_center_gap = center_share_childless - center_share_parent_u18,
    parent_owner_gap = owner_share_parent_u18 - owner_share_childless,
    parent_rooms_gap = mean_rooms_parent_u18 - mean_rooms_childless,
    analysis_weight = coalesce(hh_weight_parent_u18, 0) + coalesce(hh_weight_childless, 0)
  ) %>%
  left_join(
    metro_premia %>%
      select(
        met2013,
        family_size_central_premium_m,
        family_size_central_premium_l,
        unit_family_size_central_premium_m,
        central_m_rental_share_of_rentals,
        central_family_rental_share_of_rentals,
        central_m_units_rental_share,
        central_family_stock_share,
        family_space_central_scarcity,
        family_space_central_access,
        center_stock_weight
      ),
    by = "met2013"
  )

young_reg_data <- df %>%
  filter(
    young_sample,
    parent_status %in% c("childless", "parent_u18"),
    !is.na(log_hhincome_std),
    !is.na(college),
    tenure %in% c("owner", "renter")
  ) %>%
  mutate(
    any_parent = as.integer(parent_status == "parent_u18")
  ) %>%
  left_join(
    metro_premia %>%
      select(met2013, family_size_central_premium_m, central_m_rental_share_of_rentals),
    by = "met2013"
  ) %>%
  filter(!is.na(family_size_central_premium_m))

model_within_metro <- feols(
  center ~ any_parent * log_hhincome_std + any_parent * college + tenure_factor + sex_factor + race_factor | met2013 + year + age,
  data = young_reg_data,
  weights = ~hhwt,
  cluster = ~met2013
)

model_cross_metro <- feols(
  center ~ any_parent * log_hhincome_std +
    any_parent:family_size_central_premium_m +
    any_parent:log_hhincome_std:family_size_central_premium_m +
    any_parent:central_m_rental_share_of_rentals +
    college + tenure_factor + sex_factor + race_factor | year + age,
  data = young_reg_data,
  weights = ~hhwt,
  cluster = ~met2013
)

metro_reg_data <- parent_centrality_metro %>%
  filter(
    !is.na(parent_center_gap),
    !is.na(family_size_central_premium_m),
    !is.na(family_space_central_scarcity),
    !is.na(central_m_rental_share_of_rentals),
    analysis_weight > 0
  )

metro_gap_model <- lm(
  parent_center_gap ~ family_size_central_premium_m + central_m_rental_share_of_rentals + central_family_stock_share,
  data = metro_reg_data,
  weights = analysis_weight
)

metro_access_model <- lm(
  childless_parent_center_gap ~ family_space_central_scarcity + family_size_central_premium_m + central_family_stock_share,
  data = metro_reg_data,
  weights = analysis_weight
)

regression_table <- bind_rows(
  coef_table(model_within_metro, "household_center_with_metro_fe"),
  coef_table(model_cross_metro, "household_center_family_premium")
)

metro_gap_table <- as.data.frame(summary(metro_gap_model)$coefficients)
metro_gap_table$term <- rownames(metro_gap_table)
rownames(metro_gap_table) <- NULL
names(metro_gap_table) <- c("estimate", "se", "t_stat", "p_value", "term")
metro_gap_table <- metro_gap_table %>%
  transmute(model = "metro_parent_center_gap", term, estimate, se, t_stat, p_value)

metro_access_table <- as.data.frame(summary(metro_access_model)$coefficients)
metro_access_table$term <- rownames(metro_access_table)
rownames(metro_access_table) <- NULL
names(metro_access_table) <- c("estimate", "se", "t_stat", "p_value", "term")
metro_access_table <- metro_access_table %>%
  transmute(model = "metro_childless_parent_gap_access", term, estimate, se, t_stat, p_value)

write_csv(cell_stock, file.path(out_dir, "acs_family_size_supply_cells.csv"))
write_csv(location_size_summary, file.path(out_dir, "acs_location_size_tenure_summary.csv"))
write_csv(renter_price_cells, file.path(out_dir, "acs_renter_price_cells.csv"))
write_csv(metro_premia, file.path(out_dir, "acs_family_size_premia_by_metro.csv"))
write_csv(parent_sorting, file.path(out_dir, "acs_parent_income_sorting_summary.csv"))
write_csv(parent_centrality_metro, file.path(out_dir, "acs_parent_centrality_by_metro.csv"))
write_csv(bind_rows(regression_table, metro_gap_table, metro_access_table), file.path(out_dir, "acs_family_size_regressions.csv"))

location_plot_data <- location_size_summary %>%
  group_by(mms_location, size_bin) %>%
  summarise(hh_weight = sum(hh_weight, na.rm = TRUE), .groups = "drop") %>%
  group_by(mms_location) %>%
  mutate(share = hh_weight / sum(hh_weight)) %>%
  ungroup()

p_stock <- ggplot(location_plot_data, aes(x = mms_location, y = share, fill = size_bin)) +
  geom_col(width = 0.65) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "", y = "Share of housing units", fill = "Rooms", title = "Housing Menu by Center/Periphery") +
  theme_minimal(base_size = 12)

p_gap <- ggplot(metro_reg_data, aes(x = family_size_central_premium_m, y = parent_center_gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_point(aes(size = analysis_weight), alpha = 0.55, color = "#2C5F8A") +
  geom_smooth(aes(weight = analysis_weight), method = "lm", se = FALSE, color = "#9A3B3B", linewidth = 0.8) +
  scale_size_continuous(guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Central family-size rent premium, M vs S",
    y = "Parent minus childless center share",
    title = "Parent Centrality Gap and Family-Size Premium"
  ) +
  theme_minimal(base_size = 12)

p_sort <- parent_sorting %>%
  filter(mms_location == "center", tenure %in% c("owner", "renter")) %>%
  group_by(parent_status, income_tercile) %>%
  summarise(center_share = sum(share, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(income_tercile), y = center_share, fill = parent_status)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Household income tercile",
    y = "Center share",
    fill = "",
    title = "Center Sorting by Parenthood and Income"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "acs_family_size_stock_menu.png"), p_stock, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_parent_gap_vs_family_premium.png"), p_gap, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_center_sorting_by_parent_income.png"), p_sort, width = 7.5, height = 4.5, dpi = 180)

p_access_stock <- ggplot(metro_reg_data, aes(x = family_space_central_scarcity, y = childless_parent_center_gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_point(aes(size = analysis_weight), alpha = 0.55, color = "#2C5F8A") +
  geom_smooth(aes(weight = analysis_weight), method = "lm", se = FALSE, color = "#9A3B3B", linewidth = 0.8) +
  scale_size_continuous(guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Central small/family stock ratio relative to periphery",
    y = "Childless minus parent center share",
    title = "Parent Centrality Gap and Family-Space Stock Scarcity"
  ) +
  theme_minimal(base_size = 12)

p_access_price <- ggplot(metro_reg_data, aes(x = family_size_central_premium_m, y = childless_parent_center_gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_point(aes(size = analysis_weight), alpha = 0.55, color = "#2C5F8A") +
  geom_smooth(aes(weight = analysis_weight), method = "lm", se = FALSE, color = "#9A3B3B", linewidth = 0.8) +
  scale_size_continuous(guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Central family-size rent premium, M vs S",
    y = "Childless minus parent center share",
    title = "Parent Centrality Gap and Family-Space Price Premium"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "acs_childless_parent_gap_vs_stock_scarcity.png"), p_access_stock, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_childless_parent_gap_vs_price_premium.png"), p_access_price, width = 7.5, height = 4.5, dpi = 180)

aggregate_size <- location_plot_data %>%
  mutate(label = paste0(mms_location, "_", size_bin)) %>%
  select(label, share) %>%
  tidyr::pivot_wider(names_from = label, values_from = share)

premium_summary <- metro_premia %>%
  filter(!is.na(family_size_central_premium_m), center_stock_weight > 0) %>%
  summarise(
    n_metros = n(),
    weighted_mean_family_premium_m = weighted_mean_safe(family_size_central_premium_m, center_stock_weight),
    median_family_premium_m = median(family_size_central_premium_m, na.rm = TRUE),
    weighted_mean_unit_family_premium_m = weighted_mean_safe(unit_family_size_central_premium_m, center_stock_weight),
    weighted_mean_central_m_rental_share = weighted_mean_safe(central_m_rental_share_of_rentals, center_stock_weight),
    weighted_mean_central_family_rental_share = weighted_mean_safe(central_family_rental_share_of_rentals, center_stock_weight)
  )

gap_summary <- parent_centrality_metro %>%
  filter(!is.na(parent_center_gap), analysis_weight > 0) %>%
  summarise(
    n_metros = n(),
    weighted_parent_center_gap = weighted_mean_safe(parent_center_gap, analysis_weight),
    weighted_parent_owner_gap = weighted_mean_safe(parent_owner_gap, analysis_weight),
    weighted_parent_rooms_gap = weighted_mean_safe(parent_rooms_gap, analysis_weight)
  )

coef_lookup <- bind_rows(regression_table, metro_gap_table, metro_access_table)
get_coef_line <- function(model_name, term_name) {
  row <- coef_lookup %>% filter(model == model_name, term == term_name)
  if (nrow(row) == 0) {
    return("not estimated")
  }
  sprintf("%.4f (SE %.4f, p=%.3g)", row$estimate[1], row$se[1], row$p_value[1])
}

md <- c(
  "# ACS/MMS Family-Size Supply and Sorting Packet",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  "- Source: ACS/IPUMS `extract27.dta` joined to the MMS center/periphery PUMA lookup.",
  "- Household unit sample: household heads (`pernum == 1`), non-group-quarter households, years 2012-2023.",
  "- Sorting sample: household heads ages 22-45, childless or with children under 18.",
  sprintf("- Middle PUMAs are assigned to: `%s`.", middle_target),
  "- Size bins: `S_1_4` = 1-4 rooms, `M_5_6` = 5-6 rooms, `L_7plus` = 7+ rooms.",
  "",
  "## Main Descriptives",
  "",
  sprintf("- Matched household-head records: `%s`.", comma(nrow(df))),
  sprintf("- Metros with identified family-size rent premium: `%s`.", comma(premium_summary$n_metros[1])),
  sprintf("- Weighted mean central M-size rent premium relative to small units: `%.3f` log points.", premium_summary$weighted_mean_family_premium_m[1]),
  sprintf("- Median central M-size rent premium relative to small units: `%.3f` log points.", premium_summary$median_family_premium_m[1]),
  sprintf("- Weighted mean central M rental share among central rentals: `%.3f`.", premium_summary$weighted_mean_central_m_rental_share[1]),
  sprintf("- Weighted mean central family-capable rental share among central rentals: `%.3f`.", premium_summary$weighted_mean_central_family_rental_share[1]),
  sprintf("- Weighted parent-minus-childless center-share gap across metros: `%.3f`.", gap_summary$weighted_parent_center_gap[1]),
  sprintf("- Weighted childless-minus-parent center-share gap across metros: `%.3f`.", -gap_summary$weighted_parent_center_gap[1]),
  sprintf("- Weighted parent-minus-childless rooms gap across metros: `%.3f` rooms.", gap_summary$weighted_parent_rooms_gap[1]),
  "",
  "## Regression Reads",
  "",
  sprintf("- Within-metro parent effect on central residence: `%s`.", get_coef_line("household_center_with_metro_fe", "any_parent")),
  sprintf("- Within-metro parent x income interaction: `%s`.", get_coef_line("household_center_with_metro_fe", "any_parent:log_hhincome_std")),
  sprintf("- Household-level parent x family-size premium interaction: `%s`.", get_coef_line("household_center_family_premium", "any_parent:family_size_central_premium_m")),
  sprintf("- Household-level parent x income x family-size premium interaction: `%s`.", get_coef_line("household_center_family_premium", "any_parent:log_hhincome_std:family_size_central_premium_m")),
  sprintf("- Metro parent centrality gap slope on family-size premium: `%s`.", get_coef_line("metro_parent_center_gap", "family_size_central_premium_m")),
  sprintf("- Metro parent centrality gap slope on central M rental share: `%s`.", get_coef_line("metro_parent_center_gap", "central_m_rental_share_of_rentals")),
  sprintf("- Metro childless-minus-parent gap slope on central family-space stock scarcity: `%s`.", get_coef_line("metro_childless_parent_gap_access", "family_space_central_scarcity")),
  sprintf("- Metro childless-minus-parent gap slope on family-size price premium: `%s`.", get_coef_line("metro_childless_parent_gap_access", "family_size_central_premium_m")),
  "",
  "## Outputs",
  "",
  "- `acs_family_size_supply_cells.csv`: metro x location x size x tenure stock and rent cells.",
  "- `acs_location_size_tenure_summary.csv`: aggregate center/periphery housing menu.",
  "- `acs_family_size_premia_by_metro.csv`: family-size central premia, central/periphery family-stock shares, and stock-scarcity measures.",
  "- `acs_parent_income_sorting_summary.csv`: parent/childless center sorting by income tercile.",
  "- `acs_parent_centrality_by_metro.csv`: metro parent centrality gaps and joined supply measures.",
  "- `acs_family_size_regressions.csv`: household and metro regression coefficients.",
  "- PNG figures for the housing menu, parent gap vs premium, center sorting by parent/income, and the childless-parent gap versus stock scarcity and price premium.",
  "",
  "## Interpretation Guardrails",
  "",
  "These are descriptive ACS cross-sections. They can document the housing menu, type composition, and parent centrality gaps. They do not by themselves identify within-household birth-induced moves or causal fertility timing responses."
)

writeLines(md, file.path(out_dir, "ACS_MMS_FAMILY_SIZE_SUPPLY_PACKET.md"))
message("Wrote family-size supply packet to: ", out_dir)
