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

post_moment <- function(rows, moment, value, sample, source, formula, status, n = NA_real_, weight = NA_real_) {
  rows[[length(rows) + 1]] <- tibble(
    moment = moment,
    value = as.numeric(value),
    sample = sample,
    n = as.numeric(n),
    weight = as.numeric(weight),
    source = source,
    formula = formula,
    status = status
  )
  rows
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

script_dir <- get_script_dir()
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_dir <- file.path(script_dir, "output_intergen_one_market_targets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

middle_target <- tolower(Sys.getenv("MMS_MIDDLE_TARGET", ""))
if (middle_target == "") {
  include_middle_in_center <- tolower(Sys.getenv("MMS_INCLUDE_MIDDLE_IN_CENTER", "true")) %in% c("1", "true", "yes", "y")
  middle_target <- if (include_middle_in_center) "center" else "drop"
}

extract_path <- file.path("code", "data", "Spatial_aggregate_withmicrodata", "raw_data", "extract27.dta")
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

message("Reading ACS/IPUMS household-head extract for one-market intergen targets...")
df <- read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "met2013", "gq", "pernum", "hhwt",
    "ownershp", "rooms", "rent", "hhincome", "nchild", "yngch", "age"
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
    rent = as.numeric(rent),
    hhincome = as.numeric(hhincome),
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
    age <= 85,
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
    tenure = case_when(owner ~ "owner", renter ~ "renter", TRUE ~ "other"),
    childless = nchild == 0 & !is.na(nchild),
    parent_u18 = nchild > 0 & !is.na(yngch) & yngch != 99 & yngch < 18,
    prime_30_55 = age >= 30 & age <= 55,
    rent_monthly = if_else(renter & rent > 0, rent, NA_real_),
    rent_to_income = if_else(renter & rent > 0 & hhincome > 1000, 12 * rent / hhincome, NA_real_),
    room_bin = case_when(
      rooms <= 4 ~ "S_1_4",
      rooms <= 6 ~ "M_5_6",
      rooms >= 7 ~ "L_7plus",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(in_mms_sample, tenure %in% c("owner", "renter"))

prime_childless <- df %>% filter(prime_30_55, childless)
prime_parent <- df %>% filter(prime_30_55, parent_u18)

rows <- list()
source <- "ACS/IPUMS extract27 household heads, 2012-2023, matched MMS metro sample, middle collapsed to center"
status <- "candidate; re-audit source/sample/formula before SMM use"

rows <- post_moment(rows, "aggregate_mean_occupied_rooms_18_85",
  weighted_mean_safe(df$rooms, df$hhwt),
  "age 18-85, household head, owner or renter, positive literal ROOMS", source,
  "sum(HHWT * ROOMS) / sum(HHWT)", status,
  nrow(df), sum(df$hhwt, na.rm = TRUE))

owner_childless <- prime_childless %>% filter(owner)
renter_childless <- prime_childless %>% filter(renter)
owner_parent <- prime_parent %>% filter(owner)
renter_parent <- prime_parent %>% filter(renter)

own_mean <- weighted_mean_safe(owner_childless$rooms, owner_childless$hhwt)
rent_mean <- weighted_mean_safe(renter_childless$rooms, renter_childless$hhwt)

rows <- post_moment(rows, "prime30_55_childless_owner_mean_rooms", own_mean,
  "age 30-55, household head, childless, owner", source,
  "weighted mean rooms", status, nrow(owner_childless), sum(owner_childless$hhwt, na.rm = TRUE))
rows <- post_moment(rows, "prime30_55_childless_renter_mean_rooms", rent_mean,
  "age 30-55, household head, childless, renter", source,
  "weighted mean rooms", status, nrow(renter_childless), sum(renter_childless$hhwt, na.rm = TRUE))
rows <- post_moment(rows, "prime30_55_childless_owner_minus_renter_mean_rooms", own_mean - rent_mean,
  "age 30-55, household head, childless", source,
  "owner weighted mean rooms minus renter weighted mean rooms", status, nrow(prime_childless), sum(prime_childless$hhwt, na.rm = TRUE))

for (threshold in c(5, 6, 7)) {
  rows <- post_moment(rows, paste0("prime30_55_childless_owner_share_rooms_ge_", threshold),
    weighted_mean_safe(as.numeric(owner_childless$rooms >= threshold), owner_childless$hhwt),
    "age 30-55, household head, childless, owner", source,
    paste0("weighted share with rooms >= ", threshold), status,
    nrow(owner_childless), sum(owner_childless$hhwt, na.rm = TRUE))
  rows <- post_moment(rows, paste0("prime30_55_childless_renter_share_rooms_ge_", threshold),
    weighted_mean_safe(as.numeric(renter_childless$rooms >= threshold), renter_childless$hhwt),
    "age 30-55, household head, childless, renter", source,
    paste0("weighted share with rooms >= ", threshold), status,
    nrow(renter_childless), sum(renter_childless$hhwt, na.rm = TRUE))
}

rows <- post_moment(rows, "prime30_55_childless_renter_rent_to_income_mean",
  weighted_mean_safe(renter_childless$rent_to_income, renter_childless$hhwt),
  "age 30-55, household head, childless, renter, rent>0, income>1000", source,
  "weighted mean 12*rent/hhincome", status,
  sum(!is.na(renter_childless$rent_to_income)), sum(renter_childless$hhwt[!is.na(renter_childless$rent_to_income)], na.rm = TRUE))
rows <- post_moment(rows, "prime30_55_childless_renter_rent_to_income_median",
  weighted_quantile_safe(renter_childless$rent_to_income, renter_childless$hhwt, 0.5),
  "age 30-55, household head, childless, renter, rent>0, income>1000", source,
  "weighted median 12*rent/hhincome", status,
  sum(!is.na(renter_childless$rent_to_income)), sum(renter_childless$hhwt[!is.na(renter_childless$rent_to_income)], na.rm = TRUE))

rows <- post_moment(rows, "prime30_55_parent_owner_minus_renter_mean_rooms",
  weighted_mean_safe(owner_parent$rooms, owner_parent$hhwt) - weighted_mean_safe(renter_parent$rooms, renter_parent$hhwt),
  "age 30-55, household head, parent with child under 18", source,
  "owner weighted mean rooms minus renter weighted mean rooms", status,
  nrow(prime_parent), sum(prime_parent$hhwt, na.rm = TRUE))

targets <- bind_rows(rows)
write_csv(targets, file.path(out_dir, "intergen_one_market_acs_housing_targets.csv"))

size_cells <- df %>%
  filter(prime_30_55, childless) %>%
  group_by(tenure, room_bin) %>%
  summarise(
    n = n(),
    weight = sum(hhwt, na.rm = TRUE),
    mean_rooms = weighted_mean_safe(rooms, hhwt),
    mean_rent_to_income = weighted_mean_safe(rent_to_income, hhwt),
    .groups = "drop"
  ) %>%
  group_by(tenure) %>%
  mutate(share_within_tenure = weight / sum(weight, na.rm = TRUE)) %>%
  ungroup()
write_csv(size_cells, file.path(out_dir, "intergen_one_market_acs_childless_size_cells.csv"))

readme <- c(
  "# Intergen One-Market ACS Housing Candidate Targets",
  "",
  "These outputs are candidate replacement moments, not approved SMM targets.",
  "They must be reaudited for sample, weights, formula, room convention, and model-object match before use.",
  "",
  "- `intergen_one_market_acs_housing_targets.csv`: collapsed candidate moments.",
  "- `intergen_one_market_acs_childless_size_cells.csv`: room-bin shares by tenure for prime-age childless household heads.",
  "",
  paste0("Source: ", source),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote ACS candidate targets to ", out_dir)
