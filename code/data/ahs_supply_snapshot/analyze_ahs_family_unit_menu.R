#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
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

download_if_missing <- function(url, path) {
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    message("Downloading: ", url)
    utils::download.file(url, path, mode = "wb", quiet = FALSE)
  }
}

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

weighted_median_safe <- function(x, w) {
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
  x[which(cw >= 0.5)[1]]
}

weighted_quantile_safe <- function(x, w, prob) {
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

summarise_weighted <- function(data, group_vars) {
  data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(units = sum(weight, na.rm = TRUE), n = n(), .groups = "drop") %>%
    group_by(across(all_of(setdiff(group_vars, tail(group_vars, 1))))) %>%
    mutate(share_within_parent = units / sum(units, na.rm = TRUE)) %>%
    ungroup()
}

script_dir <- get_script_dir()
raw_dir <- file.path(script_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

sample_name <- tolower(Sys.getenv("AHS_SAMPLE", "metro"))
if (!(sample_name %in% c("metro", "national"))) {
  stop("AHS_SAMPLE must be either 'metro' or 'national'.")
}

out_dir <- file.path(
  script_dir,
  Sys.getenv("AHS_OUTPUT_DIR", paste0("output_ahs_family_unit_menu_", sample_name))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://www2.census.gov/programs-surveys/ahs/2023"
labels_url <- file.path(base_url, "AHS%202023%20Value%20Labels%20Package.zip")
labels_zip <- file.path(raw_dir, "ahs_2023_value_labels.zip")
download_if_missing(labels_url, labels_zip)

if (sample_name == "metro") {
  puf_url <- file.path(base_url, "AHS%202023%20Metropolitan%20PUF%20v1.1%20CSV.zip")
  puf_zip <- file.path(raw_dir, "ahs_2023_metro_puf_v1_1_csv.zip")
} else {
  puf_url <- file.path(base_url, "AHS%202023%20National%20PUF%20v1.1%20CSV.zip")
  puf_zip <- file.path(raw_dir, "ahs_2023_national_puf_v1_1_csv.zip")
}
download_if_missing(puf_url, puf_zip)

tmp <- tempfile("ahs_2023_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

utils::unzip(puf_zip, files = "household.csv", exdir = tmp, overwrite = TRUE)
utils::unzip(labels_zip, files = "AHS 2023 Value Labels.csv", exdir = tmp, overwrite = TRUE)

household_path <- file.path(tmp, "household.csv")
labels_path <- file.path(tmp, "AHS 2023 Value Labels.csv")

vars <- c(
  "CONTROL", "INTSTATUS", "VACANCY", "TENURE", "WEIGHT", "OMB13CBSA",
  "BLD", "BEDROOMS", "TOTROOMS", "FINROOMS", "UNITSIZE", "YRBUILT",
  "RENT", "MARKETVAL", "ADEQUACY", "NUMPEOPLE", "NUMYNGKIDS", "NUMOLDKIDS",
  "HHYNGKIDS", "HHOLDKIDS", "HSHLDTYPE"
)

message("Reading AHS 2023 ", sample_name, " household PUF...")
hh <- read_csv(
  household_path,
  col_select = any_of(vars),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  mutate(across(where(is.character), ~ str_remove_all(.x, "^'|'$")))

missing_vars <- setdiff(vars, names(hh))
if (length(missing_vars) > 0) {
  for (var in missing_vars) {
    hh[[var]] <- NA_character_
  }
}

labels <- read_csv(labels_path, col_types = cols(.default = col_character()), show_col_types = FALSE)

get_labels <- function(var) {
  labels %>%
    filter(Table %in% c("HOUSEHOLD", "AHS2023"), NAME == var) %>%
    distinct(Value, Label) %>%
    rename(value = Value, label = Label)
}

cbsa_labels <- get_labels("OMB13CBSA")
bld_labels <- get_labels("BLD")

numeric_clean <- function(x) {
  parse_number(if_else(x %in% c(".N", "N", ".M", "M", ""), NA_character_, x))
}

unit_size_midpoint <- function(x) {
  case_when(
    x == "1" ~ 400,
    x == "2" ~ 625,
    x == "3" ~ 875,
    x == "4" ~ 1250,
    x == "5" ~ 1750,
    x == "6" ~ 2250,
    x == "7" ~ 2750,
    x == "8" ~ 3500,
    x == "9" ~ 4500,
    TRUE ~ NA_real_
  )
}

hh <- hh %>%
  mutate(
    weight = numeric_clean(WEIGHT),
    int_status = case_when(
      INTSTATUS == "1" ~ "occupied",
      INTSTATUS == "2" ~ "usual_residence_elsewhere",
      INTSTATUS == "3" ~ "vacant",
      TRUE ~ "unknown"
    ),
    tenure = case_when(
      TENURE == "1" ~ "owner",
      TENURE == "2" ~ "renter",
      TENURE == "3" ~ "occupied_no_rent",
      int_status == "vacant" ~ "vacant",
      TRUE ~ "unknown"
    ),
    cbsa_code = OMB13CBSA,
    bedrooms = numeric_clean(BEDROOMS),
    rooms = numeric_clean(TOTROOMS),
    finrooms = numeric_clean(FINROOMS),
    rent = numeric_clean(RENT),
    market_value = numeric_clean(MARKETVAL),
    unit_size_mid = unit_size_midpoint(UNITSIZE),
    year_built = numeric_clean(YRBUILT),
    num_people = numeric_clean(NUMPEOPLE),
    num_young_kids = numeric_clean(NUMYNGKIDS),
    num_old_kids = numeric_clean(NUMOLDKIDS),
    any_kids = coalesce(num_young_kids, 0) + coalesce(num_old_kids, 0) > 0,
    bedroom_bin = case_when(
      bedrooms <= 1 ~ "0-1 bedrooms",
      bedrooms == 2 ~ "2 bedrooms",
      bedrooms == 3 ~ "3 bedrooms",
      bedrooms >= 4 ~ "4+ bedrooms",
      TRUE ~ "missing"
    ),
    couillard_size = case_when(
      bedrooms <= 1 ~ "small_0_1_bed",
      bedrooms == 2 ~ "mid_2_bed",
      bedrooms >= 3 ~ "large_3plus_bed",
      TRUE ~ "missing"
    ),
    room_bin = case_when(
      rooms <= 4 ~ "S: <=4 rooms",
      rooms %in% 5:6 ~ "M: 5-6 rooms",
      rooms >= 7 ~ "L: 7+ rooms",
      TRUE ~ "missing"
    ),
    structure = case_when(
      BLD == "01" ~ "mobile_home",
      BLD == "02" ~ "detached_single_family",
      BLD == "03" ~ "attached_single_family",
      BLD %in% c("04", "05") ~ "small_multifamily_2_4",
      BLD %in% c("06", "07") ~ "mid_multifamily_5_19",
      BLD %in% c("08", "09") ~ "large_multifamily_20plus",
      TRUE ~ "other_or_missing"
    ),
    missing_middle_unit = bedroom_bin %in% c("2 bedrooms", "3 bedrooms") &
      structure %in% c("attached_single_family", "small_multifamily_2_4", "mid_multifamily_5_19"),
    family_sized = bedrooms >= 3,
    family_rental = tenure == "renter" & family_sized,
    occupied = int_status == "occupied",
    rental_rent_observed = tenure == "renter" & !is.na(rent) & rent > 0 & rent < 29998,
    rent_per_bedroom = if_else(rental_rent_observed & bedrooms > 0, rent / bedrooms, NA_real_),
    rent_per_sqft = if_else(rental_rent_observed & !is.na(unit_size_mid) & unit_size_mid > 0, rent / unit_size_mid, NA_real_)
  ) %>%
  filter(!is.na(weight), weight > 0, bedroom_bin != "missing")

hh <- hh %>%
  left_join(cbsa_labels, by = c("cbsa_code" = "value")) %>%
  rename(cbsa_name = label)

if (!("cbsa_name" %in% names(hh))) {
  hh$cbsa_name <- NA_character_
}

hh <- hh %>%
  mutate(cbsa_name = coalesce(cbsa_name, cbsa_code))

overall_units <- sum(hh$weight, na.rm = TRUE)
occupied_units <- sum(hh$weight[hh$occupied], na.rm = TRUE)
renter_units <- sum(hh$weight[hh$tenure == "renter"], na.rm = TRUE)
family_sized_units <- sum(hh$weight[hh$family_sized], na.rm = TRUE)
family_sized_rentals <- sum(hh$weight[hh$family_rental], na.rm = TRUE)
missing_middle_units <- sum(hh$weight[hh$missing_middle_unit], na.rm = TRUE)

stock_by_bedroom <- hh %>%
  group_by(bedroom_bin) %>%
  summarise(units = sum(weight, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(share = units / sum(units, na.rm = TRUE))

stock_by_bedroom_tenure <- hh %>%
  group_by(tenure, bedroom_bin) %>%
  summarise(units = sum(weight, na.rm = TRUE), n = n(), .groups = "drop") %>%
  group_by(tenure) %>%
  mutate(share_within_tenure = units / sum(units, na.rm = TRUE)) %>%
  ungroup()

stock_by_bedroom_structure <- hh %>%
  group_by(structure, bedroom_bin) %>%
  summarise(units = sum(weight, na.rm = TRUE), n = n(), .groups = "drop") %>%
  group_by(structure) %>%
  mutate(share_within_structure = units / sum(units, na.rm = TRUE)) %>%
  ungroup()

stock_by_room_tenure <- hh %>%
  group_by(tenure, room_bin) %>%
  summarise(units = sum(weight, na.rm = TRUE), n = n(), .groups = "drop") %>%
  group_by(tenure) %>%
  mutate(share_within_tenure = units / sum(units, na.rm = TRUE)) %>%
  ungroup()

missing_middle_by_tenure <- hh %>%
  group_by(tenure) %>%
  summarise(
    units = sum(weight, na.rm = TRUE),
    missing_middle_units = sum(weight[missing_middle_unit], na.rm = TRUE),
    missing_middle_share = missing_middle_units / units,
    family_sized_units = sum(weight[family_sized], na.rm = TRUE),
    family_sized_share = family_sized_units / units,
    .groups = "drop"
  )

family_sized_composition <- hh %>%
  filter(family_sized) %>%
  group_by(tenure, structure) %>%
  summarise(units = sum(weight, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(share_of_family_sized = units / sum(units, na.rm = TRUE))

renter_price_by_bedroom <- hh %>%
  filter(rental_rent_observed) %>%
  group_by(bedroom_bin) %>%
  summarise(
    units = sum(weight, na.rm = TRUE),
    n = n(),
    mean_rent = weighted_mean_safe(rent, weight),
    median_rent = weighted_median_safe(rent, weight),
    p25_rent = weighted_quantile_safe(rent, weight, 0.25),
    p75_rent = weighted_quantile_safe(rent, weight, 0.75),
    mean_rent_per_bedroom = weighted_mean_safe(rent_per_bedroom, weight),
    mean_rent_per_sqft = weighted_mean_safe(rent_per_sqft, weight),
    .groups = "drop"
  )

renter_price_by_bedroom_structure <- hh %>%
  filter(rental_rent_observed) %>%
  group_by(structure, bedroom_bin) %>%
  summarise(
    units = sum(weight, na.rm = TRUE),
    n = n(),
    median_rent = weighted_median_safe(rent, weight),
    mean_rent_per_sqft = weighted_mean_safe(rent_per_sqft, weight),
    .groups = "drop"
  )

occupied_by_kids_and_bedrooms <- hh %>%
  filter(occupied) %>%
  group_by(any_kids, bedroom_bin) %>%
  summarise(units = sum(weight, na.rm = TRUE), n = n(), .groups = "drop") %>%
  group_by(any_kids) %>%
  mutate(share_within_kid_status = units / sum(units, na.rm = TRUE)) %>%
  ungroup()

children_space_mismatch <- hh %>%
  filter(occupied, any_kids) %>%
  group_by(tenure) %>%
  summarise(
    units_with_children = sum(weight, na.rm = TRUE),
    share_0_1_bed = sum(weight[bedrooms <= 1], na.rm = TRUE) / units_with_children,
    share_0_2_bed = sum(weight[bedrooms <= 2], na.rm = TRUE) / units_with_children,
    mean_bedrooms = weighted_mean_safe(bedrooms, weight),
    mean_rooms = weighted_mean_safe(rooms, weight),
    .groups = "drop"
  )

metro_menu <- tibble()
if (any(!is.na(hh$cbsa_code) & !(hh$cbsa_code %in% c("99998", "99999")))) {
  metro_menu <- hh %>%
    filter(!is.na(cbsa_code), !(cbsa_code %in% c("99998", "99999"))) %>%
    group_by(cbsa_code, cbsa_name) %>%
    summarise(
      units = sum(weight, na.rm = TRUE),
      renter_units = sum(weight[tenure == "renter"], na.rm = TRUE),
      family_sized_share = sum(weight[family_sized], na.rm = TRUE) / units,
      family_sized_rental_share_of_all = sum(weight[family_rental], na.rm = TRUE) / units,
      family_sized_share_within_rentals = sum(weight[family_rental], na.rm = TRUE) / renter_units,
      missing_middle_share = sum(weight[missing_middle_unit], na.rm = TRUE) / units,
      rental_missing_middle_share = sum(weight[missing_middle_unit & tenure == "renter"], na.rm = TRUE) / renter_units,
      renter_child_share_0_2_bed = sum(weight[occupied & tenure == "renter" & any_kids & bedrooms <= 2], na.rm = TRUE) /
        sum(weight[occupied & tenure == "renter" & any_kids], na.rm = TRUE),
      median_rent_1bed = weighted_median_safe(rent[tenure == "renter" & bedroom_bin == "0-1 bedrooms" & rental_rent_observed], weight[tenure == "renter" & bedroom_bin == "0-1 bedrooms" & rental_rent_observed]),
      median_rent_2bed = weighted_median_safe(rent[tenure == "renter" & bedroom_bin == "2 bedrooms" & rental_rent_observed], weight[tenure == "renter" & bedroom_bin == "2 bedrooms" & rental_rent_observed]),
      median_rent_3bed = weighted_median_safe(rent[tenure == "renter" & bedroom_bin == "3 bedrooms" & rental_rent_observed], weight[tenure == "renter" & bedroom_bin == "3 bedrooms" & rental_rent_observed]),
      .groups = "drop"
    ) %>%
    mutate(
      rent_3bed_vs_1bed_log = log(median_rent_3bed) - log(median_rent_1bed),
      rent_3bed_vs_2bed_log = log(median_rent_3bed) - log(median_rent_2bed),
      family_rental_scarcity_index = 1 - family_sized_share_within_rentals
    ) %>%
    arrange(desc(family_rental_scarcity_index))
}

write_csv(stock_by_bedroom, file.path(out_dir, "ahs_stock_by_bedroom.csv"))
write_csv(stock_by_bedroom_tenure, file.path(out_dir, "ahs_stock_by_bedroom_tenure.csv"))
write_csv(stock_by_bedroom_structure, file.path(out_dir, "ahs_stock_by_bedroom_structure.csv"))
write_csv(stock_by_room_tenure, file.path(out_dir, "ahs_stock_by_room_tenure.csv"))
write_csv(missing_middle_by_tenure, file.path(out_dir, "ahs_missing_middle_by_tenure.csv"))
write_csv(family_sized_composition, file.path(out_dir, "ahs_family_sized_composition.csv"))
write_csv(renter_price_by_bedroom, file.path(out_dir, "ahs_renter_price_by_bedroom.csv"))
write_csv(renter_price_by_bedroom_structure, file.path(out_dir, "ahs_renter_price_by_bedroom_structure.csv"))
write_csv(occupied_by_kids_and_bedrooms, file.path(out_dir, "ahs_occupied_by_kids_and_bedrooms.csv"))
write_csv(children_space_mismatch, file.path(out_dir, "ahs_children_space_mismatch.csv"))
write_csv(metro_menu, file.path(out_dir, "ahs_metro_family_unit_menu.csv"))

bedroom_levels <- c("0-1 bedrooms", "2 bedrooms", "3 bedrooms", "4+ bedrooms")
room_levels <- c("S: <=4 rooms", "M: 5-6 rooms", "L: 7+ rooms")

plot_stock_bedroom <- stock_by_bedroom %>%
  mutate(bedroom_bin = factor(bedroom_bin, levels = bedroom_levels))

p0 <- ggplot(plot_stock_bedroom, aes(x = bedroom_bin, y = share)) +
  geom_col(fill = "#4c78a8") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Share of housing units", title = "AHS 2023 physical bedroom menu") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "ahs_bedroom_menu_overall.png"), p0, width = 7.2, height = 4.2, dpi = 180)

plot_stock_rooms <- hh %>%
  group_by(room_bin) %>%
  summarise(units = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    share = units / sum(units, na.rm = TRUE),
    room_bin = factor(room_bin, levels = room_levels)
  )

p_rooms <- ggplot(plot_stock_rooms, aes(x = room_bin, y = share)) +
  geom_col(fill = "#54a24b") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Share of housing units", title = "AHS 2023 room-size menu") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "ahs_room_menu_overall.png"), p_rooms, width = 7.2, height = 4.2, dpi = 180)

plot_structure_bedroom <- hh %>%
  group_by(structure, bedroom_bin) %>%
  summarise(units = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    share_total = units / sum(units, na.rm = TRUE),
    bedroom_bin = factor(bedroom_bin, levels = bedroom_levels)
  ) %>%
  group_by(structure) %>%
  mutate(structure_units = sum(units, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(structure = reorder(structure, structure_units))

p_structure_bedroom <- ggplot(plot_structure_bedroom, aes(x = share_total, y = structure, fill = bedroom_bin)) +
  geom_col() +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "0-1 bedrooms" = "#4c78a8",
    "2 bedrooms" = "#72b7b2",
    "3 bedrooms" = "#f58518",
    "4+ bedrooms" = "#e45756"
  )) +
  labs(x = "Share of all housing units", y = NULL, fill = NULL, title = "Physical stock by structure and bedrooms") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(out_dir, "ahs_structure_by_bedroom_share.png"), p_structure_bedroom, width = 7.8, height = 4.8, dpi = 180)

plot_bedroom <- stock_by_bedroom_tenure %>%
  filter(tenure %in% c("owner", "renter", "vacant")) %>%
  mutate(
    bedroom_bin = factor(bedroom_bin, levels = bedroom_levels),
    tenure = factor(tenure, levels = c("renter", "owner", "vacant"))
  )

p1 <- ggplot(plot_bedroom, aes(x = bedroom_bin, y = share_within_tenure, fill = tenure)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(renter = "#4c78a8", owner = "#f58518", vacant = "#54a24b")) +
  labs(x = NULL, y = "Share within tenure", fill = NULL, title = "AHS 2023 bedroom menu by tenure") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(out_dir, "ahs_bedroom_menu_by_tenure.png"), p1, width = 7.5, height = 4.5, dpi = 180)

plot_kids_bedroom <- occupied_by_kids_and_bedrooms %>%
  mutate(
    bedroom_bin = factor(bedroom_bin, levels = bedroom_levels),
    child_status = if_else(any_kids, "Household has children", "No children in household")
  )

p_kids <- ggplot(plot_kids_bedroom, aes(x = bedroom_bin, y = share_within_kid_status, fill = child_status)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Household has children" = "#e45756", "No children in household" = "#4c78a8")) +
  labs(x = NULL, y = "Share within household type", fill = NULL, title = "Bedroom consumption by child status") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(out_dir, "ahs_bedroom_menu_by_children.png"), p_kids, width = 7.5, height = 4.5, dpi = 180)

plot_rent_gradient <- renter_price_by_bedroom %>%
  mutate(bedroom_bin = factor(bedroom_bin, levels = bedroom_levels))

p_rent <- ggplot(plot_rent_gradient, aes(x = bedroom_bin, y = median_rent, group = 1)) +
  geom_line(color = "#b279a2", linewidth = 0.8) +
  geom_point(color = "#b279a2", size = 2.5) +
  scale_y_continuous(labels = dollar_format(accuracy = 1)) +
  labs(x = NULL, y = "Weighted median monthly rent", title = "Rental price gradient by bedroom count") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "ahs_rent_gradient_by_bedroom.png"), p_rent, width = 7.2, height = 4.2, dpi = 180)

plot_structure <- family_sized_composition %>%
  filter(tenure %in% c("owner", "renter")) %>%
  mutate(structure = reorder(structure, units))

p2 <- ggplot(plot_structure, aes(x = units, y = structure, fill = tenure)) +
  geom_col() +
  scale_x_continuous(labels = comma) +
  scale_fill_manual(values = c(renter = "#4c78a8", owner = "#f58518")) +
  labs(x = "Weighted units", y = NULL, fill = NULL, title = "Where 3+ bedroom units sit in the stock") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(out_dir, "ahs_family_sized_units_by_structure_tenure.png"), p2, width = 7.5, height = 4.8, dpi = 180)

top_scarce <- metro_menu %>%
  filter(!is.na(family_sized_share_within_rentals), renter_units > 0) %>%
  slice_max(order_by = family_rental_scarcity_index, n = 12) %>%
  mutate(cbsa_name = str_replace(cbsa_name, ",.*$", ""))

if (nrow(top_scarce) > 0) {
  p3 <- ggplot(top_scarce, aes(x = family_rental_scarcity_index, y = reorder(cbsa_name, family_rental_scarcity_index))) +
    geom_col(fill = "#b279a2") +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "1 - share of rentals with 3+ bedrooms", y = NULL, title = "Selected metros with scarce family-sized rentals") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "ahs_family_rental_scarcity_metros.png"), p3, width = 7.5, height = 4.8, dpi = 180)
}

pdf(file.path(out_dir, "AHS_2023_FAMILY_UNIT_FIGURE_PACKET.pdf"), width = 8, height = 5)
print(p0)
print(p_rooms)
print(p_structure_bedroom)
print(p1)
print(p_kids)
print(p_rent)
print(p2)
if (exists("p3")) {
  print(p3)
}
dev.off()

fmt_pct <- function(x) percent(x, accuracy = 0.1)
fmt_num <- function(x) comma(x, accuracy = 1)
fmt_dol <- function(x) dollar(x, accuracy = 1)

renter_bed <- stock_by_bedroom_tenure %>% filter(tenure == "renter")
owner_bed <- stock_by_bedroom_tenure %>% filter(tenure == "owner")
renter_small_share <- renter_bed %>% filter(bedroom_bin == "0-1 bedrooms") %>% pull(share_within_tenure)
renter_family_share <- renter_bed %>% filter(bedroom_bin %in% c("3 bedrooms", "4+ bedrooms")) %>% summarise(x = sum(share_within_tenure, na.rm = TRUE)) %>% pull(x)
owner_family_share <- owner_bed %>% filter(bedroom_bin %in% c("3 bedrooms", "4+ bedrooms")) %>% summarise(x = sum(share_within_tenure, na.rm = TRUE)) %>% pull(x)

kids_renter_mismatch <- children_space_mismatch %>% filter(tenure == "renter") %>% pull(share_0_2_bed)
kids_owner_mismatch <- children_space_mismatch %>% filter(tenure == "owner") %>% pull(share_0_2_bed)

rent_1 <- renter_price_by_bedroom %>% filter(bedroom_bin == "0-1 bedrooms") %>% pull(median_rent)
rent_2 <- renter_price_by_bedroom %>% filter(bedroom_bin == "2 bedrooms") %>% pull(median_rent)
rent_3 <- renter_price_by_bedroom %>% filter(bedroom_bin == "3 bedrooms") %>% pull(median_rent)

mm_renter_share <- missing_middle_by_tenure %>% filter(tenure == "renter") %>% pull(missing_middle_share)
mm_owner_share <- missing_middle_by_tenure %>% filter(tenure == "owner") %>% pull(missing_middle_share)

scarce_lines <- "- Metro table not available in this AHS sample."
if (nrow(metro_menu) > 0) {
  scarce_lines <- metro_menu %>%
    filter(!is.na(family_sized_share_within_rentals), renter_units > 0) %>%
    slice_max(order_by = family_rental_scarcity_index, n = 8) %>%
    transmute(line = sprintf(
      "- %s: 3+ bedroom share of rentals `%s`; rental missing-middle share `%s`; median 3-bed rent `%s`.",
      cbsa_name,
      fmt_pct(family_sized_share_within_rentals),
      fmt_pct(rental_missing_middle_share),
      fmt_dol(median_rent_3bed)
    )) %>%
    pull(line)
}

md <- c(
  "# AHS 2023 Family-Unit Housing Menu",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  sprintf("- Source: AHS 2023 `%s` public-use household CSV, version 1.1.", sample_name),
  "- Weight: `WEIGHT`.",
  "- Unit-size bins follow the Couillard-style bedroom contrast: small `0-1`, middle `2`, and large/family-sized `3+` bedrooms.",
  "- Missing-middle proxy: 2-3 bedroom units in attached single-family, 2-4 unit multifamily, or 5-19 unit multifamily structures.",
  "- Family-sized rental scarcity proxy: `1 - Pr(3+ bedrooms | rental stock)` within metro.",
  "- This is a stock/menu snapshot. It is not causal evidence that family-size scarcity changes fertility.",
  "",
  "## National/Metro Stock Snapshot",
  "",
  sprintf("- Weighted housing units in this PUF read: `%s`; occupied units: `%s`.", fmt_num(overall_units), fmt_num(occupied_units)),
  sprintf("- Family-sized units, 3+ bedrooms: `%s` of all units.", fmt_pct(family_sized_units / overall_units)),
  sprintf("- Family-sized rentals, 3+ bedrooms and rented: `%s` of all units; `%s` of 3+ bedroom units.", fmt_pct(family_sized_rentals / overall_units), fmt_pct(family_sized_rentals / family_sized_units)),
  sprintf("- Missing-middle proxy units: `%s` of all units.", fmt_pct(missing_middle_units / overall_units)),
  sprintf("- Missing-middle share within rentals: `%s`; within owners: `%s`.", fmt_pct(mm_renter_share), fmt_pct(mm_owner_share)),
  "",
  "## Tenure Contrast",
  "",
  sprintf("- Rental stock share with 0-1 bedrooms: `%s`.", fmt_pct(renter_small_share)),
  sprintf("- Rental stock share with 3+ bedrooms: `%s`.", fmt_pct(renter_family_share)),
  sprintf("- Owner stock share with 3+ bedrooms: `%s`.", fmt_pct(owner_family_share)),
  sprintf("- Renter households with children in 0-2 bedroom units: `%s`; owner households with children in 0-2 bedroom units: `%s`.", fmt_pct(kids_renter_mismatch), fmt_pct(kids_owner_mismatch)),
  "",
  "## Rental Price Gradient",
  "",
  sprintf("- Median rent, 0-1 bedrooms: `%s`; 2 bedrooms: `%s`; 3 bedrooms: `%s`.", fmt_dol(rent_1), fmt_dol(rent_2), fmt_dol(rent_3)),
  sprintf("- Log 3-bedroom versus 1-bedroom median-rent gradient: `%.3f`.", log(rent_3) - log(rent_1)),
  sprintf("- Log 3-bedroom versus 2-bedroom median-rent gradient: `%.3f`.", log(rent_3) - log(rent_2)),
  "",
  "## Selected Metro Scarcity Reads",
  "",
  scarce_lines,
  "",
  "## Outputs",
  "",
  "- `ahs_stock_by_bedroom.csv`",
  "- `ahs_stock_by_bedroom_tenure.csv`",
  "- `ahs_stock_by_bedroom_structure.csv`",
  "- `ahs_stock_by_room_tenure.csv`",
  "- `ahs_missing_middle_by_tenure.csv`",
  "- `ahs_family_sized_composition.csv`",
  "- `ahs_renter_price_by_bedroom.csv`",
  "- `ahs_renter_price_by_bedroom_structure.csv`",
  "- `ahs_occupied_by_kids_and_bedrooms.csv`",
  "- `ahs_children_space_mismatch.csv`",
  "- `ahs_metro_family_unit_menu.csv`",
  "- PNG figures and `AHS_2023_FAMILY_UNIT_FIGURE_PACKET.pdf` for the physical bedroom/room menu, structure, tenure, child-status consumption, rent gradient, and selected metro family-rental scarcity.",
  "",
  "## Interpretation",
  "",
  "The missing-middle question is not whether large units exist. It is whether family-capable units are available outside the detached-owner bundle. A Couillard-style 3+ bedroom policy margin maps naturally to the family-sized rental share and the 2-3 bedroom attached/multifamily share. If these objects are low, the fertility-relevant housing constraint is a tenure/structure menu constraint, not just an average rent level."
)

writeLines(md, file.path(out_dir, "AHS_2023_FAMILY_UNIT_MENU.md"))

message("Wrote AHS family-unit menu packet to: ", out_dir)
