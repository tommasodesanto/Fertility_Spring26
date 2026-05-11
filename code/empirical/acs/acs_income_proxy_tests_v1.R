#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
  library(fixest)
  library(ggplot2)
  library(scales)
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

interaction_term_name <- function(model, v1, v2) {
  nms <- names(coef(model))
  hit <- nms[grepl(v1, nms, fixed = TRUE) & grepl(v2, nms, fixed = TRUE)]
  if (length(hit) == 0) {
    stop("Missing interaction term for ", v1, " and ", v2)
  }
  hit[1]
}

lincom_affine <- function(model, base_term, inter_term, q) {
  b <- coef(model)
  v <- vcov(model)
  est <- b[[base_term]] + q * b[[inter_term]]
  vv <- v[base_term, base_term] + q^2 * v[inter_term, inter_term] + 2 * q * v[base_term, inter_term]
  se <- sqrt(vv)
  list(estimate = unname(est), se = unname(se))
}

extract_parent_penalty_by_income <- function(model, parent_term, income_term, label, q_low, q_high) {
  inter_name <- interaction_term_name(model, parent_term, income_term)
  low <- lincom_affine(model, parent_term, inter_name, q_low)
  high <- lincom_affine(model, parent_term, inter_name, q_high)
  tibble(
    outcome = label,
    parent_main = unname(coef(model)[[parent_term]]),
    parent_main_se = sqrt(vcov(model)[parent_term, parent_term]),
    parent_income_interaction = unname(coef(model)[[inter_name]]),
    parent_income_interaction_se = sqrt(vcov(model)[inter_name, inter_name]),
    income_p25 = q_low,
    income_p75 = q_high,
    parent_effect_p25 = low$estimate,
    parent_effect_p25_se = low$se,
    parent_effect_p75 = high$estimate,
    parent_effect_p75_se = high$se,
    n = nobs(model)
  )
}

extract_income_slope_by_parent <- function(model, income_term, parent_term, label) {
  inter_name <- interaction_term_name(model, income_term, parent_term)
  childless_slope <- unname(coef(model)[[income_term]])
  childless_se <- sqrt(vcov(model)[income_term, income_term])
  parent <- lincom_affine(model, income_term, inter_name, 1)
  tibble(
    outcome = label,
    childless_slope = childless_slope,
    childless_slope_se = childless_se,
    parent_interaction = unname(coef(model)[[inter_name]]),
    parent_interaction_se = sqrt(vcov(model)[inter_name, inter_name]),
    parent_slope = parent$estimate,
    parent_slope_se = parent$se,
    n = nobs(model)
  )
}

script_dir <- get_script_dir()
spring26_root <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(spring26_root, "code", "data", "mms_center_periphery", "data")
out_dir <- file.path(script_dir, "output", "acs_income_proxy_tests_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

middle_target <- "center"
lookup_2010_path <- file.path(data_dir, "puma_mms_lookup_2010.csv")
lookup_2020_path <- file.path(data_dir, "puma_mms_lookup_2020.csv")
extract_path <- file.path(spring26_root, "code", "data", "Spatial_aggregate_withmicrodata", "raw_data", "extract27.dta")

lookup_2010 <- read_csv(lookup_2010_path, show_col_types = FALSE)
lookup_2020 <- read_csv(lookup_2020_path, show_col_types = FALSE)

dest_lookup <- bind_rows(
  lookup_2010 %>%
    transmute(
      lookup_period = "pre2022",
      statefip,
      puma,
      cbsacode,
      cbsatitle,
      dest_mms_location = mms_location
    ),
  lookup_2020 %>%
    transmute(
      lookup_period = "post2021",
      statefip,
      puma,
      cbsacode,
      cbsatitle,
      dest_mms_location = mms_location
    )
)

df <- read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "metro", "met2013", "migmet131",
    "age", "sex", "race", "perwt", "gq", "ownershp", "rooms", "rent", "hhincome",
    "nchild", "nchlt5", "eldch", "migrate1"
  )
) %>%
  transmute(
    year = as.integer(year),
    statefip = as.integer(statefip),
    puma = as.integer(puma),
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
    migrate1 = as.integer(migrate1)
  ) %>%
  filter(
    year >= 2012,
    age >= 22,
    age <= 45,
    gq %in% c(1, 2)
  ) %>%
  mutate(
    owner = as.integer(ownershp == 1),
    renter = as.integer(ownershp == 2),
    has_children = nchild > 0 & !is.na(nchild),
    childless = nchild == 0 & !is.na(nchild),
    newparent = as.integer(eldch < 4 & eldch != 99 & nchild > 0),
    any_parent = as.integer(has_children),
    parent_compare = case_when(
      childless ~ "Non-Parents",
      newparent == 1 ~ "New Parents",
      TRUE ~ NA_character_
    ),
    moved1y = migrate1 >= 2 & !is.na(migrate1),
    same_msa = moved1y & met2013 > 0 & migmet131 > 0 & met2013 == migmet131,
    lookup_period = if_else(year <= 2021, "pre2022", "post2021"),
    ln_income = ifelse(hhincome > 1000, log(hhincome), NA_real_),
    ln_rooms = ifelse(rooms > 0, log(rooms), NA_real_),
    rent_per_room = ifelse(renter == 1 & rent > 0 & rooms > 0, rent / rooms, NA_real_)
  )

q_lo <- quantile(df$hhincome[df$hhincome > 1000], probs = 0.01, na.rm = TRUE)
q_hi <- quantile(df$hhincome[df$hhincome > 1000], probs = 0.99, na.rm = TRUE)

df_mms <- df %>%
  filter(is.na(hhincome) | (hhincome >= q_lo & hhincome <= q_hi)) %>%
  left_join(dest_lookup, by = c("lookup_period", "statefip", "puma", "met2013" = "cbsacode")) %>%
  mutate(
    dest_mms_location = collapse_mms_location(dest_mms_location, middle_target),
    dest_mms_location = if_else(dest_mms_location %in% c("center", "periphery"), dest_mms_location, NA_character_),
    dest_label = case_when(
      dest_mms_location == "center" ~ "Center",
      dest_mms_location == "periphery" ~ "Periphery",
      TRUE ~ NA_character_
    ),
    center_live = as.integer(dest_label == "Center"),
    center_dest = as.integer(dest_label == "Center")
  ) %>%
  filter(!is.na(dest_label), !is.na(cbsatitle), !is.na(ln_income))

restricted_np <- df_mms %>%
  filter(!is.na(parent_compare))

within_cbsa <- restricted_np %>%
  filter(same_msa)

income_q25 <- quantile(restricted_np$ln_income, 0.25, na.rm = TRUE)
income_q75 <- quantile(restricted_np$ln_income, 0.75, na.rm = TRUE)

location_model <- feols(
  center_live ~ newparent * ln_income + i(sex) + i(race) | age + year,
  data = restricted_np,
  weights = ~perwt
)

dest_model <- feols(
  center_dest ~ newparent * ln_income + i(sex) + i(race) | age + year,
  data = within_cbsa,
  weights = ~perwt
)

rooms_model <- feols(
  ln_rooms ~ ln_income * any_parent + i(sex) + i(race) | age + year,
  data = df_mms %>% filter(!is.na(ln_rooms)),
  weights = ~perwt
)

own_model <- feols(
  owner ~ ln_income * any_parent + i(sex) + i(race) | age + year,
  data = df_mms,
  weights = ~perwt
)

parent_penalty_summary <- bind_rows(
  extract_parent_penalty_by_income(location_model, "newparent", "ln_income", "center_live_newparent_income", income_q25, income_q75),
  extract_parent_penalty_by_income(dest_model, "newparent", "ln_income", "center_dest_newparent_income", income_q25, income_q75)
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

housing_slope_summary <- bind_rows(
  extract_income_slope_by_parent(rooms_model, "ln_income", "any_parent", "ln_rooms_anyparent_income"),
  extract_income_slope_by_parent(own_model, "ln_income", "any_parent", "own_anyparent_income")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

write_csv(parent_penalty_summary, file.path(out_dir, "acs_income_parent_penalty_v1.csv"))
write_csv(housing_slope_summary, file.path(out_dir, "acs_income_housing_slopes_v1.csv"))

income_bin_summary <- bind_rows(
  restricted_np %>%
    mutate(
      income_tercile = ntile(ln_income, 3),
      sample = "all_mms",
      outcome = center_live
    ) %>%
    group_by(sample, parent_compare, income_tercile) %>%
    summarise(
      share = weighted_mean_safe(outcome, perwt),
      n = n(),
      .groups = "drop"
    ),
  within_cbsa %>%
    mutate(
      income_tercile = ntile(ln_income, 3),
      sample = "within_cbsa_movers",
      outcome = center_dest
    ) %>%
    group_by(sample, parent_compare, income_tercile) %>%
    summarise(
      share = weighted_mean_safe(outcome, perwt),
      n = n(),
      .groups = "drop"
    )
) %>%
  mutate(share = round(share, 6))

write_csv(income_bin_summary, file.path(out_dir, "acs_income_proxy_terciles_v1.csv"))

p_income <- income_bin_summary %>%
  mutate(
    income_tercile = factor(income_tercile, levels = 1:3, labels = c("Low", "Mid", "High")),
    sample = recode(sample, all_mms = "Center residence", within_cbsa_movers = "Center destination among within-CBSA movers")
  ) %>%
  ggplot(aes(x = income_tercile, y = share, color = parent_compare, group = parent_compare)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~sample) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Household income tercile", y = "Center share", color = "", title = "ACS income proxy for resource-moderated family flight") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "acs_income_proxy_terciles_v1.png"), p_income, width = 10, height = 5, dpi = 180)

metro_rent <- df_mms %>%
  filter(renter == 1, !is.na(rent_per_room), rent_per_room > 0) %>%
  group_by(cbsatitle, met2013, dest_label) %>%
  summarise(
    rent_per_room = weighted_mean_safe(rent_per_room, perwt),
    rent_obs = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = dest_label,
    values_from = c(rent_per_room, rent_obs),
    names_glue = "{.value}_{dest_label}"
  ) %>%
  mutate(
    rent_gap_ratio = rent_per_room_Center / rent_per_room_Periphery,
    log_rent_gap = log(rent_gap_ratio)
  )

metro_center_gap <- restricted_np %>%
  mutate(parent_key = if_else(parent_compare == "Non-Parents", "nonparent", "newparent")) %>%
  group_by(cbsatitle, met2013, parent_key) %>%
  summarise(
    center_share = weighted_mean_safe(center_live, perwt),
    group_obs = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = parent_key,
    values_from = c(center_share, group_obs),
    names_glue = "{.value}_{parent_key}"
  ) %>%
  transmute(
    cbsatitle,
    met2013,
    center_share_nonparent = center_share_nonparent,
    center_share_newparent = center_share_newparent,
    gap_pp = 100 * (center_share_nonparent - center_share_newparent),
    metro_obs = group_obs_nonparent + group_obs_newparent,
    np_obs = group_obs_nonparent,
    newparent_obs = group_obs_newparent
  )

metro_dest_gap <- within_cbsa %>%
  mutate(parent_key = if_else(parent_compare == "Non-Parents", "nonparent", "newparent")) %>%
  group_by(cbsatitle, met2013, parent_key) %>%
  summarise(
    center_dest_share = weighted_mean_safe(center_dest, perwt),
    group_obs = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = parent_key,
    values_from = c(center_dest_share, group_obs),
    names_glue = "{.value}_{parent_key}"
  ) %>%
  transmute(
    cbsatitle,
    met2013,
    center_dest_nonparent = center_dest_share_nonparent,
    center_dest_newparent = center_dest_share_newparent,
    gap_pp = 100 * (center_dest_nonparent - center_dest_newparent),
    mover_obs = group_obs_nonparent + group_obs_newparent,
    np_obs = group_obs_nonparent,
    newparent_obs = group_obs_newparent
  )

metro_sorting <- metro_center_gap %>%
  inner_join(metro_rent, by = c("cbsatitle", "met2013")) %>%
  filter(
    np_obs >= 200,
    newparent_obs >= 50,
    rent_obs_Center >= 100,
    rent_obs_Periphery >= 100,
    is.finite(log_rent_gap)
  )

metro_dest_sorting <- metro_dest_gap %>%
  inner_join(metro_rent, by = c("cbsatitle", "met2013")) %>%
  filter(
    np_obs >= 50,
    newparent_obs >= 25,
    rent_obs_Center >= 100,
    rent_obs_Periphery >= 100,
    is.finite(log_rent_gap)
  )

metro_sorting_model <- feols(gap_pp ~ log_rent_gap, data = metro_sorting, weights = ~metro_obs)
metro_dest_model <- feols(gap_pp ~ log_rent_gap, data = metro_dest_sorting, weights = ~mover_obs)

metro_model_summary <- tibble(
  outcome = c("center_residence_gap_pp", "within_cbsa_center_dest_gap_pp"),
  slope = c(unname(coef(metro_sorting_model)[["log_rent_gap"]]), unname(coef(metro_dest_model)[["log_rent_gap"]])),
  slope_se = c(unname(se(metro_sorting_model)[["log_rent_gap"]]), unname(se(metro_dest_model)[["log_rent_gap"]])),
  p_value = c(unname(pvalue(metro_sorting_model)[["log_rent_gap"]]), unname(pvalue(metro_dest_model)[["log_rent_gap"]])),
  n_metros = c(nobs(metro_sorting_model), nobs(metro_dest_model))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

write_csv(metro_sorting, file.path(out_dir, "acs_cross_metro_center_gap_v1.csv"))
write_csv(metro_dest_sorting, file.path(out_dir, "acs_cross_metro_dest_gap_v1.csv"))
write_csv(metro_model_summary, file.path(out_dir, "acs_cross_metro_regressions_v1.csv"))

p_metro1 <- ggplot(metro_sorting, aes(x = log_rent_gap, y = gap_pp, size = metro_obs)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "firebrick") +
  labs(
    x = "Log center/periphery rent-per-room gap",
    y = "Non-parent minus new-parent center share (pp)",
    size = "Obs.",
    title = "Cross-metro sorting gap vs rent gap"
  ) +
  theme_minimal(base_size = 12)

p_metro2 <- ggplot(metro_dest_sorting, aes(x = log_rent_gap, y = gap_pp, size = mover_obs)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "navy") +
  labs(
    x = "Log center/periphery rent-per-room gap",
    y = "Non-parent minus new-parent center-destination gap (pp)",
    size = "Obs.",
    title = "Cross-metro within-CBSA mover gap vs rent gap"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "acs_cross_metro_center_gap_v1.png"), p_metro1, width = 8, height = 5, dpi = 180)
ggsave(file.path(out_dir, "acs_cross_metro_dest_gap_v1.png"), p_metro2, width = 8, height = 5, dpi = 180)

cat("ACS income proxy and cross-metro tests complete.\n")
cat(file.path(out_dir, "acs_income_parent_penalty_v1.csv"), "\n")
cat(file.path(out_dir, "acs_income_housing_slopes_v1.csv"), "\n")
cat(file.path(out_dir, "acs_cross_metro_regressions_v1.csv"), "\n")
