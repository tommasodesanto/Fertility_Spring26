#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
  library(ggplot2)
  library(tidyr)
  library(scales)
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

lookup_2010_path <- file.path(data_dir, "puma_mms_lookup_2010.csv")
lookup_2020_path <- file.path(data_dir, "puma_mms_lookup_2020.csv")

if (!file.exists(lookup_2010_path) || !file.exists(lookup_2020_path)) {
  stop("Run build_mms_geography.R before compute_mms_descriptives.R.")
}

extract_path <- "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
if (!file.exists(extract_path)) {
  stop("Missing microdata extract: ", extract_path)
}

lookup_2010 <- read_csv(lookup_2010_path, show_col_types = FALSE)
lookup_2020 <- read_csv(lookup_2020_path, show_col_types = FALSE)

message("Loading extract27 columns needed for MMS descriptives...")
df <- read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "metro", "met2013", "age", "sex", "perwt", "gq",
    "ownershp", "rooms", "rent", "hhincome", "nchild", "nchlt5", "eldch", "yngch"
  )
) %>%
  transmute(
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
    rooms = as.numeric(rooms),
    rent = as.numeric(rent),
    hhincome = as.numeric(hhincome),
    nchild = as.numeric(nchild),
    nchlt5 = as.numeric(nchlt5),
    eldch = as.numeric(eldch),
    yngch = as.numeric(yngch)
  ) %>%
  filter(
    year >= 2012,
    age >= 22,
    gq %in% c(1, 2)
  ) %>%
  mutate(
    owner = ownershp == 1,
    renter = ownershp == 2,
    has_youngchild = nchlt5 > 0 & !is.na(nchlt5),
    has_child_u18 = nchild > 0 & !is.na(yngch) & yngch != 99 & yngch < 18,
    has_children = nchild > 0 & !is.na(nchild),
    childless = nchild == 0 & !is.na(nchild),
    newparent = eldch < 4 & eldch != 99 & nchild > 0,
    rti = ifelse(rent > 0 & hhincome > 0, 12 * rent / hhincome, NA_real_)
  )

df_2010 <- df %>%
  filter(year <= 2021) %>%
  left_join(
    lookup_2010,
    by = c("statefip", "puma", "met2013" = "cbsacode")
  )

df_2020 <- df %>%
  filter(year >= 2022) %>%
  left_join(
    lookup_2020,
    by = c("statefip", "puma", "met2013" = "cbsacode")
  )

df_mms <- bind_rows(df_2010, df_2020) %>%
  mutate(mms_location = collapse_mms_location(mms_location, middle_target)) %>%
  filter(mms_location %in% c("center", "periphery"))

df_mms_window <- df_mms %>%
  filter(age <= 45)

if (nrow(df_mms) == 0) {
  stop("No matched households after applying MMS lookup.")
}

message("Computing descriptive outputs...")
location_summary <- df_mms_window %>%
  group_by(mms_location) %>%
  summarise(
    obs = n(),
    pop_weight = sum(perwt, na.rm = TRUE),
    owner_rate = weighted_mean_safe(owner, perwt),
    renter_rate = weighted_mean_safe(renter, perwt),
    mean_rooms = weighted_mean_safe(rooms, perwt),
    mean_rent = weighted_mean_safe(rent, perwt),
    mean_rti = weighted_mean_safe(rti, perwt),
    has_child_u18_rate = weighted_mean_safe(has_child_u18, perwt),
    has_youngchild_rate = weighted_mean_safe(has_youngchild, perwt),
    has_children_rate = weighted_mean_safe(has_children, perwt),
    childless_rate = weighted_mean_safe(childless, perwt),
    newparent_rate = weighted_mean_safe(newparent, perwt),
    .groups = "drop"
  ) %>%
  mutate(pop_share = pop_weight / sum(pop_weight))

tenure_parent_summary <- df_mms_window %>%
  mutate(
    tenure = case_when(
      owner ~ "Owner",
      renter ~ "Renter",
      TRUE ~ "Other"
    ),
    parent_status = case_when(
      has_children ~ "With children",
      childless ~ "No children",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(tenure %in% c("Owner", "Renter"), parent_status != "Unknown") %>%
  group_by(mms_location, tenure, parent_status) %>%
  summarise(
    pop_weight = sum(perwt, na.rm = TRUE),
    mean_rooms = weighted_mean_safe(rooms, perwt),
    mean_rent = weighted_mean_safe(rent, perwt),
    mean_rti = weighted_mean_safe(rti, perwt),
    has_youngchild_rate = weighted_mean_safe(has_youngchild, perwt),
    .groups = "drop"
  )

age_profiles <- df_mms_window %>%
  group_by(mms_location, age) %>%
  summarise(
    pop_weight = sum(perwt, na.rm = TRUE),
    owner_rate = weighted_mean_safe(owner, perwt),
    has_child_u18_rate = weighted_mean_safe(has_child_u18, perwt),
    has_children_rate = weighted_mean_safe(has_children, perwt),
    has_youngchild_rate = weighted_mean_safe(has_youngchild, perwt),
    mean_rooms = weighted_mean_safe(rooms, perwt),
    .groups = "drop"
  )

age_profiles_full <- df_mms %>%
  group_by(mms_location, age) %>%
  summarise(
    pop_weight = sum(perwt, na.rm = TRUE),
    owner_rate = weighted_mean_safe(owner, perwt),
    has_child_u18_rate = weighted_mean_safe(has_child_u18, perwt),
    has_children_rate = weighted_mean_safe(has_children, perwt),
    has_youngchild_rate = weighted_mean_safe(has_youngchild, perwt),
    mean_rooms = weighted_mean_safe(rooms, perwt),
    .groups = "drop"
  )

city_summary <- df_mms_window %>%
  group_by(cbsatitle, met2013, mms_location) %>%
  summarise(
    pop_weight = sum(perwt, na.rm = TRUE),
    owner_rate = weighted_mean_safe(owner, perwt),
    has_youngchild_rate = weighted_mean_safe(has_youngchild, perwt),
    mean_rooms = weighted_mean_safe(rooms, perwt),
    .groups = "drop"
  )

write_csv(location_summary, file.path(out_dir, "mms_location_summary.csv"))
write_csv(tenure_parent_summary, file.path(out_dir, "mms_tenure_parent_summary.csv"))
write_csv(age_profiles, file.path(out_dir, "mms_age_profiles.csv"))
write_csv(age_profiles_full, file.path(out_dir, "mms_age_profiles_full.csv"))
write_csv(city_summary, file.path(out_dir, "mms_city_summary.csv"))

p1 <- ggplot(location_summary, aes(x = mms_location, y = owner_rate, fill = mms_location)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = percent(owner_rate, accuracy = 0.1)), vjust = -0.4, size = 4) +
  scale_y_continuous(labels = percent_format(), limits = c(0, max(location_summary$owner_rate) * 1.15)) +
  labs(x = "", y = "Ownership rate", title = "MMS Geography: Ownership by Location") +
  theme_minimal(base_size = 13)

p2 <- ggplot(location_summary, aes(x = mms_location, y = has_youngchild_rate, fill = mms_location)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = percent(has_youngchild_rate, accuracy = 0.1)), vjust = -0.4, size = 4) +
  scale_y_continuous(labels = percent_format(), limits = c(0, max(location_summary$has_youngchild_rate) * 1.15)) +
  labs(x = "", y = "Share with child < 5", title = "MMS Geography: Fertility Proxy by Location") +
  theme_minimal(base_size = 13)

p3 <- ggplot(tenure_parent_summary, aes(x = tenure, y = mean_rooms, fill = parent_status)) +
  geom_col(position = "dodge", width = 0.6) +
  facet_wrap(~mms_location) +
  labs(x = "", y = "Mean rooms", fill = "", title = "Rooms by Tenure, Parenthood, and MMS Location") +
  theme_minimal(base_size = 13)

png(file.path(out_dir, "mms_location_profiles.png"), width = 1800, height = 600, res = 150)
layout(matrix(1:3, nrow = 1))
print(p1)
print(p2)
print(p3)
dev.off()

summary_lines <- c(
  "# MMS Geography Summary",
  "",
  sprintf("Core tract population target share: %.2f", core_pop_share),
  sprintf("Middle target: %s", middle_target),
  sprintf("Matched sample years: %d-%d", min(df_mms$year), max(df_mms$year)),
  sprintf("Matched weighted sample: %.0f", sum(df_mms$perwt, na.rm = TRUE)),
  "",
  "## By location",
  ""
)

for (i in seq_len(nrow(location_summary))) {
  row <- location_summary[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "- %s: pop share %.1f%%, owner %.1f%%, child<5 %.1f%%, mean rooms %.2f, mean RTI %.3f",
      row$mms_location,
      100 * row$pop_share,
      100 * row$owner_rate,
      100 * row$has_youngchild_rate,
      row$mean_rooms,
      row$mean_rti
    )
  )
}

writeLines(summary_lines, file.path(out_dir, "mms_summary.md"))

message("MMS descriptives complete.")
