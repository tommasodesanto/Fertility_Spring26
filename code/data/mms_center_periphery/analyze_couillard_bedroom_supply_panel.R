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
  ct %>% transmute(model = model_name, term, estimate, se, t_stat, p_value)
}

script_dir <- get_script_dir()
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_dir <- file.path(script_dir, Sys.getenv("COUILLARD_BEDROOM_OUTPUT_DIR", "output_couillard_bedroom_supply"))
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

message("Loading ACS/IPUMS household-head sample for Couillard bedroom supply panel...")
df <- read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "met2013", "gq", "pernum", "hhwt",
    "ownershp", "rooms", "bedrooms", "rent", "nchild", "yngch", "age"
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
    ownershp = as.integer(ownershp),
    rooms = as.numeric(rooms),
    bedrooms = as.numeric(bedrooms),
    rent = as.numeric(rent),
    nchild = as.numeric(nchild),
    yngch = as.numeric(yngch),
    age = as.integer(age),
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
    !is.na(bedrooms),
    bedrooms >= 0,
    bedrooms < 99
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
    rent_per_bedroom = if_else(renter & rent > 0 & bedrooms > 0, rent / bedrooms, NA_real_),
    bedroom_bin = case_when(
      bedrooms <= 1 ~ "B_0_1",
      bedrooms == 2 ~ "B_2",
      bedrooms >= 3 ~ "B_3plus",
      TRUE ~ NA_character_
    ),
    bedroom_family = bedroom_bin == "B_3plus",
    parent_status = case_when(
      nchild > 0 & !is.na(yngch) & yngch != 99 & yngch < 18 ~ "parent_u18",
      nchild == 0 & !is.na(nchild) ~ "childless",
      TRUE ~ NA_character_
    ),
    young_sample = age >= 22 & age <= 45
  ) %>%
  filter(!is.na(mms_location), !is.na(bedroom_bin))

if (nrow(df) == 0) {
  stop("No matched MMS household-head records.")
}

message("Computing Couillard-style bedroom stock, price, sorting, and panel outputs...")

stock_cells <- df %>%
  group_by(met2013, cbsatitle, mms_location, bedroom_bin, tenure) %>%
  summarise(
    n_hh = n(),
    hh_weight = sum(hhwt, na.rm = TRUE),
    mean_bedrooms = weighted_mean_safe(bedrooms, hhwt),
    mean_rooms = weighted_mean_safe(rooms, hhwt),
    median_rent = weighted_quantile_safe(rent_monthly, hhwt, 0.5),
    median_rent_per_bedroom = weighted_quantile_safe(rent_per_bedroom, hhwt, 0.5),
    .groups = "drop"
  ) %>%
  group_by(met2013, mms_location) %>%
  mutate(share_within_location = hh_weight / sum(hh_weight, na.rm = TRUE)) %>%
  ungroup()

stock_panel <- df %>%
  group_by(year, met2013, cbsatitle, mms_location, bedroom_bin, tenure) %>%
  summarise(
    hh_weight = sum(hhwt, na.rm = TRUE),
    median_rent = weighted_quantile_safe(rent_monthly, hhwt, 0.5),
    .groups = "drop"
  )

location_summary <- df %>%
  group_by(mms_location, bedroom_bin, tenure) %>%
  summarise(
    n_hh = n(),
    hh_weight = sum(hhwt, na.rm = TRUE),
    median_rent = weighted_quantile_safe(rent_monthly, hhwt, 0.5),
    .groups = "drop"
  ) %>%
  group_by(mms_location) %>%
  mutate(share_within_location = hh_weight / sum(hh_weight, na.rm = TRUE)) %>%
  ungroup()

renter_price <- df %>%
  filter(renter, !is.na(rent_monthly), rent_monthly > 0) %>%
  group_by(met2013, cbsatitle, mms_location, bedroom_bin) %>%
  summarise(
    renter_weight = sum(hhwt, na.rm = TRUE),
    median_rent = weighted_quantile_safe(rent_monthly, hhwt, 0.5),
    median_rent_per_bedroom = weighted_quantile_safe(rent_per_bedroom, hhwt, 0.5),
    .groups = "drop"
  )

price_cell <- function(loc, bin, prefix) {
  renter_price %>%
    filter(mms_location == loc, bedroom_bin == bin) %>%
    transmute(
      met2013,
      cbsatitle,
      !!paste0(prefix, "_renter_weight") := renter_weight,
      !!paste0(prefix, "_median_rent") := median_rent,
      !!paste0(prefix, "_median_rent_per_bedroom") := median_rent_per_bedroom
    )
}

location_supply <- stock_cells %>%
  filter(tenure %in% c("owner", "renter")) %>%
  group_by(met2013, cbsatitle, mms_location) %>%
  summarise(
    stock_small = sum(if_else(bedroom_bin == "B_0_1", hh_weight, 0), na.rm = TRUE),
    stock_two = sum(if_else(bedroom_bin == "B_2", hh_weight, 0), na.rm = TRUE),
    stock_family = sum(if_else(bedroom_bin == "B_3plus", hh_weight, 0), na.rm = TRUE),
    stock_total = stock_small + stock_two + stock_family,
    family_stock_share = stock_family / stock_total,
    small_to_family_ratio = stock_small / stock_family,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = mms_location,
    values_from = c(stock_small, stock_two, stock_family, stock_total, family_stock_share, small_to_family_ratio)
  )

premia <- df %>%
  distinct(met2013, cbsatitle) %>%
  left_join(price_cell("center", "B_0_1", "center_small"), by = c("met2013", "cbsatitle")) %>%
  left_join(price_cell("periphery", "B_0_1", "periphery_small"), by = c("met2013", "cbsatitle")) %>%
  left_join(price_cell("center", "B_3plus", "center_family"), by = c("met2013", "cbsatitle")) %>%
  left_join(price_cell("periphery", "B_3plus", "periphery_family"), by = c("met2013", "cbsatitle")) %>%
  left_join(location_supply, by = c("met2013", "cbsatitle")) %>%
  mutate(
    bedroom_family_price_premium = (log(center_family_median_rent) - log(periphery_family_median_rent)) -
      (log(center_small_median_rent) - log(periphery_small_median_rent)),
    bedroom_family_price_premium_per_bedroom = (log(center_family_median_rent_per_bedroom) - log(periphery_family_median_rent_per_bedroom)) -
      (log(center_small_median_rent_per_bedroom) - log(periphery_small_median_rent_per_bedroom)),
    bedroom_family_stock_scarcity = log(small_to_family_ratio_center) - log(small_to_family_ratio_periphery),
    bedroom_family_access = family_stock_share_center - family_stock_share_periphery,
    center_stock_weight = stock_total_center
  )

parent_centrality <- df %>%
  filter(young_sample, parent_status %in% c("childless", "parent_u18")) %>%
  group_by(met2013, cbsatitle, parent_status) %>%
  summarise(
    hh_weight = sum(hhwt, na.rm = TRUE),
    center_share = weighted_mean_safe(center, hhwt),
    owner_share = weighted_mean_safe(owner, hhwt),
    mean_bedrooms = weighted_mean_safe(bedrooms, hhwt),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = parent_status, values_from = c(hh_weight, center_share, owner_share, mean_bedrooms)) %>%
  mutate(
    childless_parent_center_gap = center_share_childless - center_share_parent_u18,
    parent_owner_gap = owner_share_parent_u18 - owner_share_childless,
    parent_bedroom_gap = mean_bedrooms_parent_u18 - mean_bedrooms_childless,
    analysis_weight = coalesce(hh_weight_parent_u18, 0) + coalesce(hh_weight_childless, 0)
  ) %>%
  left_join(
    premia %>%
      select(
        met2013,
        bedroom_family_price_premium,
        bedroom_family_price_premium_per_bedroom,
        bedroom_family_stock_scarcity,
        bedroom_family_access,
        family_stock_share_center,
        center_stock_weight
      ),
    by = "met2013"
  )

reg_data <- parent_centrality %>%
  filter(
    analysis_weight > 0,
    !is.na(childless_parent_center_gap),
    !is.na(bedroom_family_stock_scarcity),
    !is.na(bedroom_family_price_premium)
  )

gap_model <- lm(
  childless_parent_center_gap ~ bedroom_family_stock_scarcity + bedroom_family_price_premium + family_stock_share_center,
  data = reg_data,
  weights = analysis_weight
)

early_late <- df %>%
  filter(year %in% c(2012, 2013, 2014, 2021, 2022, 2023), tenure %in% c("owner", "renter")) %>%
  mutate(period = case_when(
    year %in% c(2012, 2013, 2014) ~ "early_2012_2014",
    year %in% c(2021, 2022, 2023) ~ "late_2021_2023",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  group_by(met2013, cbsatitle, mms_location, period) %>%
  summarise(
    stock_total = sum(hhwt, na.rm = TRUE),
    stock_family = sum(hhwt[bedroom_bin == "B_3plus"], na.rm = TRUE),
    stock_small = sum(hhwt[bedroom_bin == "B_0_1"], na.rm = TRUE),
    family_stock_share = stock_family / stock_total,
    small_to_family_ratio = stock_small / stock_family,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = c(mms_location, period),
    values_from = c(stock_total, stock_family, stock_small, family_stock_share, small_to_family_ratio)
  ) %>%
  mutate(
    center_family_stock_growth_log = log(stock_family_center_late_2021_2023) - log(stock_family_center_early_2012_2014),
    center_family_share_change = family_stock_share_center_late_2021_2023 - family_stock_share_center_early_2012_2014,
    center_small_family_ratio_change = log(small_to_family_ratio_center_late_2021_2023) - log(small_to_family_ratio_center_early_2012_2014)
  ) %>%
  left_join(
    premia %>%
      select(met2013, bedroom_family_price_premium, bedroom_family_stock_scarcity),
    by = "met2013"
  )

supply_response_data <- early_late %>%
  filter(
    !is.na(center_family_share_change),
    !is.na(bedroom_family_price_premium),
    is.finite(center_family_share_change),
    is.finite(bedroom_family_price_premium)
  )

supply_response_model <- lm(
  center_family_share_change ~ bedroom_family_price_premium + bedroom_family_stock_scarcity,
  data = supply_response_data
)

regression_table <- bind_rows(
  coef_table(gap_model, "bedroom_parent_gap"),
  coef_table(supply_response_model, "bedroom_supply_response")
)

write_csv(stock_cells, file.path(out_dir, "acs_bedroom_supply_cells.csv"))
write_csv(stock_panel, file.path(out_dir, "acs_bedroom_supply_panel.csv"))
write_csv(location_summary, file.path(out_dir, "acs_bedroom_location_summary.csv"))
write_csv(renter_price, file.path(out_dir, "acs_bedroom_renter_price_cells.csv"))
write_csv(premia, file.path(out_dir, "acs_bedroom_premia_by_metro.csv"))
write_csv(parent_centrality, file.path(out_dir, "acs_bedroom_parent_centrality_by_metro.csv"))
write_csv(early_late, file.path(out_dir, "acs_bedroom_early_late_supply_response.csv"))
write_csv(regression_table, file.path(out_dir, "acs_bedroom_supply_regressions.csv"))

bedroom_levels <- c("B_0_1", "B_2", "B_3plus")

p_menu <- location_summary %>%
  group_by(mms_location, bedroom_bin) %>%
  summarise(hh_weight = sum(hh_weight, na.rm = TRUE), .groups = "drop") %>%
  group_by(mms_location) %>%
  mutate(
    share = hh_weight / sum(hh_weight, na.rm = TRUE),
    bedroom_bin = factor(bedroom_bin, levels = bedroom_levels)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = mms_location, y = share, fill = bedroom_bin)) +
  geom_col(width = 0.65) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Share of housing units", fill = "Bedrooms", title = "Bedroom Menu by Center/Periphery") +
  theme_minimal(base_size = 12)

p_gap_stock <- ggplot(reg_data, aes(x = bedroom_family_stock_scarcity, y = childless_parent_center_gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_point(aes(size = analysis_weight), alpha = 0.55, color = "#2C5F8A") +
  geom_smooth(aes(weight = analysis_weight), method = "lm", se = FALSE, color = "#9A3B3B", linewidth = 0.8) +
  scale_size_continuous(guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Central 0-1/3+ bedroom stock ratio relative to periphery",
    y = "Childless minus parent center share",
    title = "Parent Centrality Gap and 3+ Bedroom Stock Scarcity"
  ) +
  theme_minimal(base_size = 12)

p_gap_price <- ggplot(reg_data, aes(x = bedroom_family_price_premium, y = childless_parent_center_gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_point(aes(size = analysis_weight), alpha = 0.55, color = "#2C5F8A") +
  geom_smooth(aes(weight = analysis_weight), method = "lm", se = FALSE, color = "#9A3B3B", linewidth = 0.8) +
  scale_size_continuous(guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Central 3+ bedroom rent premium relative to 0-1 bedrooms",
    y = "Childless minus parent center share",
    title = "Parent Centrality Gap and 3+ Bedroom Price Premium"
  ) +
  theme_minimal(base_size = 12)

p_supply <- ggplot(supply_response_data, aes(x = bedroom_family_price_premium, y = center_family_share_change)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_point(alpha = 0.65, color = "#2C5F8A") +
  geom_smooth(method = "lm", se = FALSE, color = "#9A3B3B", linewidth = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    x = "Current central 3+ bedroom rent premium",
    y = "Change in central 3+ bedroom stock share, early to late",
    title = "Descriptive Supply Response: Family Stock Share"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "acs_bedroom_menu_center_periphery.png"), p_menu, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_bedroom_parent_gap_vs_stock_scarcity.png"), p_gap_stock, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_bedroom_parent_gap_vs_price_premium.png"), p_gap_price, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_bedroom_supply_response.png"), p_supply, width = 7.5, height = 4.5, dpi = 180)

summary_stats <- tibble(
  statistic = c(
    "matched_records",
    "metros_with_bedroom_premium",
    "weighted_childless_parent_center_gap",
    "weighted_parent_owner_gap",
    "weighted_parent_bedroom_gap",
    "weighted_mean_bedroom_stock_scarcity",
    "weighted_mean_bedroom_price_premium",
    "weighted_mean_center_family_share_change"
  ),
  value = c(
    nrow(df),
    sum(!is.na(premia$bedroom_family_price_premium)),
    weighted_mean_safe(reg_data$childless_parent_center_gap, reg_data$analysis_weight),
    weighted_mean_safe(reg_data$parent_owner_gap, reg_data$analysis_weight),
    weighted_mean_safe(reg_data$parent_bedroom_gap, reg_data$analysis_weight),
    weighted_mean_safe(premia$bedroom_family_stock_scarcity, premia$center_stock_weight),
    weighted_mean_safe(premia$bedroom_family_price_premium, premia$center_stock_weight),
    weighted_mean_safe(supply_response_data$center_family_share_change, supply_response_data$stock_total_center_late_2021_2023)
  )
)
write_csv(summary_stats, file.path(out_dir, "acs_bedroom_summary_stats.csv"))

coef_line <- function(model_name, term_name) {
  row <- regression_table %>% filter(model == model_name, term == term_name)
  if (nrow(row) == 0) {
    return("not estimated")
  }
  sprintf("%.4f (SE %.4f, p=%.3g)", row$estimate[1], row$se[1], row$p_value[1])
}

md <- c(
  "# ACS/MMS Couillard-Style Bedroom Supply Packet",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  "- Source: ACS/IPUMS `extract27.dta` joined to the MMS center/periphery PUMA lookup.",
  "- Household-head sample, years 2012-2023, non-group-quarter households.",
  "- Bedroom bins follow Couillard's size contrast: `B_0_1`, `B_2`, and `B_3plus`.",
  sprintf("- Middle PUMAs are assigned to: `%s`.", middle_target),
  "- Supply-response diagnostic compares early 2012-2014 to late 2021-2023. It is descriptive, not causal.",
  "",
  "## Main Facts",
  "",
  sprintf("- Matched household-head records: `%s`.", comma(summary_stats$value[summary_stats$statistic == "matched_records"])),
  sprintf("- Metros with identified 3+ bedroom price premium: `%s`.", comma(summary_stats$value[summary_stats$statistic == "metros_with_bedroom_premium"])),
  sprintf("- Weighted childless-minus-parent center-share gap: `%.3f`.", summary_stats$value[summary_stats$statistic == "weighted_childless_parent_center_gap"]),
  sprintf("- Weighted parent-minus-childless owner-share gap: `%.3f`.", summary_stats$value[summary_stats$statistic == "weighted_parent_owner_gap"]),
  sprintf("- Weighted parent-minus-childless bedroom gap: `%.3f` bedrooms.", summary_stats$value[summary_stats$statistic == "weighted_parent_bedroom_gap"]),
  sprintf("- Weighted mean central 0-1/3+ bedroom stock scarcity: `%.3f`.", summary_stats$value[summary_stats$statistic == "weighted_mean_bedroom_stock_scarcity"]),
  sprintf("- Weighted mean central 3+ bedroom price premium: `%.3f` log points.", summary_stats$value[summary_stats$statistic == "weighted_mean_bedroom_price_premium"]),
  sprintf("- Weighted mean change in central 3+ bedroom stock share, early to late: `%.3f`.", summary_stats$value[summary_stats$statistic == "weighted_mean_center_family_share_change"]),
  "",
  "## Regression Reads",
  "",
  sprintf("- Parent centrality gap slope on 3+ bedroom stock scarcity: `%s`.", coef_line("bedroom_parent_gap", "bedroom_family_stock_scarcity")),
  sprintf("- Parent centrality gap slope on 3+ bedroom price premium: `%s`.", coef_line("bedroom_parent_gap", "bedroom_family_price_premium")),
  sprintf("- Central 3+ bedroom stock-share change slope on price premium: `%s`.", coef_line("bedroom_supply_response", "bedroom_family_price_premium")),
  sprintf("- Central 3+ bedroom stock-share change slope on stock scarcity: `%s`.", coef_line("bedroom_supply_response", "bedroom_family_stock_scarcity")),
  "",
  "## Outputs",
  "",
  "- `acs_bedroom_supply_cells.csv`",
  "- `acs_bedroom_supply_panel.csv`",
  "- `acs_bedroom_location_summary.csv`",
  "- `acs_bedroom_renter_price_cells.csv`",
  "- `acs_bedroom_premia_by_metro.csv`",
  "- `acs_bedroom_parent_centrality_by_metro.csv`",
  "- `acs_bedroom_early_late_supply_response.csv`",
  "- `acs_bedroom_supply_regressions.csv`",
  "- `acs_bedroom_summary_stats.csv`",
  "- PNG figures for the bedroom menu, parent centrality gap, and descriptive stock-share response.",
  "",
  "## Interpretation Guardrails",
  "",
  "This packet uses bedrooms, which is closer to Couillard than the room-bin packet. It is still descriptive. The supply-response exercise uses observed ACS stock changes, not an exogenous supply shifter, and should not be read as a causal elasticity."
)

writeLines(md, file.path(out_dir, "ACS_MMS_COUILLARD_BEDROOM_SUPPLY_PACKET.md"))
message("Wrote Couillard-style bedroom supply packet to: ", out_dir)
