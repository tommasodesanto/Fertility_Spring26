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

coef_table <- function(model, model_name) {
  ct <- as.data.frame(coeftable(model))
  ct$term <- rownames(ct)
  rownames(ct) <- NULL
  names(ct) <- c("estimate", "se", "t_stat", "p_value", "term")
  ct %>%
    transmute(model = model_name, term, estimate, se, t_stat, p_value)
}

decompose_gap <- function(data, outcome) {
  cells <- data %>%
    filter(!is.na(.data[[outcome]]), !is.na(income_quintile), mms_location %in% c("center", "periphery")) %>%
    group_by(mms_location, income_quintile) %>%
    summarise(
      weight = sum(perwt, na.rm = TRUE),
      mean_outcome = weighted_mean_safe(.data[[outcome]], perwt),
      .groups = "drop"
    ) %>%
    group_by(mms_location) %>%
    mutate(share = weight / sum(weight, na.rm = TRUE)) %>%
    ungroup()

  wide <- cells %>%
    select(mms_location, income_quintile, mean_outcome, share) %>%
    pivot_wider(
      names_from = mms_location,
      values_from = c(mean_outcome, share)
    ) %>%
    filter(
      !is.na(mean_outcome_center),
      !is.na(mean_outcome_periphery),
      !is.na(share_center),
      !is.na(share_periphery)
    )

  if (nrow(wide) == 0) {
    return(tibble())
  }

  center_mean <- sum(wide$share_center * wide$mean_outcome_center, na.rm = TRUE)
  periphery_mean <- sum(wide$share_periphery * wide$mean_outcome_periphery, na.rm = TRUE)
  within_component <- sum(wide$share_center * (wide$mean_outcome_center - wide$mean_outcome_periphery), na.rm = TRUE)
  composition_component <- sum((wide$share_center - wide$share_periphery) * wide$mean_outcome_periphery, na.rm = TRUE)

  tibble(
    outcome = outcome,
    center_mean = center_mean,
    periphery_mean = periphery_mean,
    center_minus_periphery = center_mean - periphery_mean,
    within_income_component = within_component,
    income_composition_component = composition_component,
    residual = (center_mean - periphery_mean) - within_component - composition_component
  )
}

script_dir <- get_script_dir()
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_dir <- file.path(script_dir, Sys.getenv("INCOME_FERTILITY_OUTPUT_DIR", "output_income_fertility_cross_section"))
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

message("Loading ACS/IPUMS women sample for income-fertility cross-section...")
df <- read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "met2013", "gq", "perwt", "ownershp", "rooms",
    "rent", "hhincome", "nchild", "nchlt5", "eldch", "yngch", "fertyr",
    "sex", "age", "race", "educ", "marst", "empstat", "incwage", "inctot"
  )
) %>%
  transmute(
    year = as.integer(year),
    statefip = as.integer(statefip),
    puma = as.integer(puma),
    met2013 = as.integer(met2013),
    gq = as.integer(gq),
    perwt = as.numeric(perwt),
    ownershp = as.integer(ownershp),
    rooms = as.numeric(rooms),
    rent = as.numeric(rent),
    hhincome = as.numeric(hhincome),
    nchild = as.numeric(nchild),
    nchlt5 = as.numeric(nchlt5),
    eldch = as.numeric(eldch),
    yngch = as.numeric(yngch),
    fertyr = as.integer(fertyr),
    sex = as.integer(sex),
    age = as.integer(age),
    race = as.integer(race),
    educ = as.integer(educ),
    marst = as.integer(marst),
    empstat = as.integer(empstat),
    incwage = as.numeric(incwage),
    inctot = as.numeric(inctot),
    lookup_period = if_else(year <= 2021, "pre2022", "post2021")
  ) %>%
  filter(
    year >= 2012,
    met2013 > 0,
    gq %in% c(1, 2),
    sex == 2,
    age >= 22,
    age <= 45,
    !is.na(perwt),
    perwt > 0
  ) %>%
  left_join(dest_lookup, by = c("lookup_period", "statefip", "puma", "met2013")) %>%
  mutate(
    mms_location = collapse_mms_location(mms_location, middle_target),
    mms_location = if_else(mms_location %in% c("center", "periphery"), mms_location, NA_character_),
    center = as.integer(mms_location == "center"),
    owner = ownershp == 1,
    renter = ownershp == 2,
    recent_birth = as.integer(fertyr == 2),
    fertyr_observed = fertyr %in% c(1, 2),
    has_child_u18 = nchild > 0 & !is.na(yngch) & yngch != 99 & yngch < 18,
    childless_in_household = nchild == 0 & !is.na(nchild),
    mean_children_household = nchild,
    parity_bin = case_when(
      nchild == 0 ~ "0",
      nchild == 1 ~ "1",
      nchild == 2 ~ "2",
      nchild >= 3 ~ "3plus",
      TRUE ~ NA_character_
    ),
    age_group = case_when(
      age <= 29 ~ "22_29",
      age <= 34 ~ "30_34",
      age <= 39 ~ "35_39",
      TRUE ~ "40_45"
    ),
    owner_label = case_when(owner ~ "owner", renter ~ "renter", TRUE ~ "other"),
    college = as.integer(educ >= 10),
    married = as.integer(marst %in% c(1, 2)),
    employed = as.integer(empstat == 1),
    hhincome_clean = if_else(hhincome > 1000 & hhincome < 2000000, hhincome, NA_real_),
    log_hhincome = if_else(!is.na(hhincome_clean), log(hhincome_clean), NA_real_),
    incwage_clean = if_else(incwage > 0 & incwage < 1000000, incwage, NA_real_),
    log_incwage = if_else(!is.na(incwage_clean), log(incwage_clean), NA_real_),
    rooms_clean = if_else(rooms >= 1 & rooms <= 12, rooms, NA_real_),
    family_capable = rooms_clean >= 5,
    m_size = rooms_clean >= 5 & rooms_clean <= 6,
    rent_burden = if_else(renter & rent > 0 & hhincome_clean > 0, 12 * rent / hhincome_clean, NA_real_)
  ) %>%
  filter(!is.na(mms_location), !is.na(log_hhincome))

if (nrow(df) == 0) {
  stop("No matched MMS women records with valid household income.")
}

df <- df %>%
  group_by(year) %>%
  mutate(
    income_quintile = as.integer(ntile(log_hhincome, 5)),
    income_tercile = as.integer(ntile(log_hhincome, 3))
  ) %>%
  ungroup()

income_resid_model <- feols(
  log_hhincome ~ i(age) + i(year) + i(educ) + i(marst),
  data = df,
  weights = ~perwt
)
df$income_residual <- resid(income_resid_model)
df <- df %>%
  mutate(
    income_resid_quintile = as.integer(ntile(income_residual, 5))
  )

message("Computing income-fertility summaries, decompositions, and regressions...")

fertility_by_income <- df %>%
  group_by(income_quintile) %>%
  summarise(
    n = n(),
    weight = sum(perwt, na.rm = TRUE),
    median_hhincome = weighted_median_safe(hhincome_clean, perwt),
    recent_birth_rate = weighted_mean_safe(if_else(fertyr_observed, recent_birth, NA_integer_), perwt),
    parent_u18_rate = weighted_mean_safe(has_child_u18, perwt),
    childless_rate = weighted_mean_safe(childless_in_household, perwt),
    mean_nchild = weighted_mean_safe(mean_children_household, perwt),
    center_share = weighted_mean_safe(center, perwt),
    owner_rate = weighted_mean_safe(owner, perwt),
    renter_rate = weighted_mean_safe(renter, perwt),
    mean_rooms = weighted_mean_safe(rooms_clean, perwt),
    family_capable_share = weighted_mean_safe(family_capable, perwt),
    m_size_share = weighted_mean_safe(m_size, perwt),
    rent_burden = weighted_mean_safe(rent_burden, perwt),
    college_rate = weighted_mean_safe(college, perwt),
    married_rate = weighted_mean_safe(married, perwt),
    employed_rate = weighted_mean_safe(employed, perwt),
    .groups = "drop"
  )

fertility_by_income_location <- df %>%
  group_by(mms_location, income_quintile) %>%
  summarise(
    n = n(),
    weight = sum(perwt, na.rm = TRUE),
    recent_birth_rate = weighted_mean_safe(if_else(fertyr_observed, recent_birth, NA_integer_), perwt),
    parent_u18_rate = weighted_mean_safe(has_child_u18, perwt),
    childless_rate = weighted_mean_safe(childless_in_household, perwt),
    mean_nchild = weighted_mean_safe(mean_children_household, perwt),
    owner_rate = weighted_mean_safe(owner, perwt),
    mean_rooms = weighted_mean_safe(rooms_clean, perwt),
    family_capable_share = weighted_mean_safe(family_capable, perwt),
    .groups = "drop"
  )

fertility_by_income_age <- df %>%
  group_by(age_group, income_quintile) %>%
  summarise(
    n = n(),
    weight = sum(perwt, na.rm = TRUE),
    recent_birth_rate = weighted_mean_safe(if_else(fertyr_observed, recent_birth, NA_integer_), perwt),
    parent_u18_rate = weighted_mean_safe(has_child_u18, perwt),
    childless_rate = weighted_mean_safe(childless_in_household, perwt),
    mean_nchild = weighted_mean_safe(mean_children_household, perwt),
    .groups = "drop"
  )

fertility_by_income_tenure <- df %>%
  filter(owner_label %in% c("owner", "renter")) %>%
  group_by(owner_label, income_quintile) %>%
  summarise(
    n = n(),
    weight = sum(perwt, na.rm = TRUE),
    recent_birth_rate = weighted_mean_safe(if_else(fertyr_observed, recent_birth, NA_integer_), perwt),
    parent_u18_rate = weighted_mean_safe(has_child_u18, perwt),
    childless_rate = weighted_mean_safe(childless_in_household, perwt),
    mean_nchild = weighted_mean_safe(mean_children_household, perwt),
    center_share = weighted_mean_safe(center, perwt),
    family_capable_share = weighted_mean_safe(family_capable, perwt),
    .groups = "drop"
  )

fertility_by_residual_income <- df %>%
  group_by(income_resid_quintile) %>%
  summarise(
    n = n(),
    weight = sum(perwt, na.rm = TRUE),
    recent_birth_rate = weighted_mean_safe(if_else(fertyr_observed, recent_birth, NA_integer_), perwt),
    parent_u18_rate = weighted_mean_safe(has_child_u18, perwt),
    childless_rate = weighted_mean_safe(childless_in_household, perwt),
    mean_nchild = weighted_mean_safe(mean_children_household, perwt),
    center_share = weighted_mean_safe(center, perwt),
    owner_rate = weighted_mean_safe(owner, perwt),
    mean_rooms = weighted_mean_safe(rooms_clean, perwt),
    college_rate = weighted_mean_safe(college, perwt),
    married_rate = weighted_mean_safe(married, perwt),
    .groups = "drop"
  )

income_distribution_by_fertility <- df %>%
  mutate(
    fertility_group = case_when(
      fertyr_observed & recent_birth == 1 ~ "recent_birth",
      has_child_u18 ~ "parent_u18",
      childless_in_household ~ "childless",
      TRUE ~ "other"
    )
  ) %>%
  filter(fertility_group %in% c("recent_birth", "parent_u18", "childless")) %>%
  group_by(fertility_group, income_quintile) %>%
  summarise(
    weight = sum(perwt, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(fertility_group) %>%
  mutate(share_within_group = weight / sum(weight, na.rm = TRUE)) %>%
  ungroup()

parity_income_distribution <- df %>%
  filter(!is.na(parity_bin)) %>%
  group_by(parity_bin, income_quintile) %>%
  summarise(
    weight = sum(perwt, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(parity_bin) %>%
  mutate(share_within_parity = weight / sum(weight, na.rm = TRUE)) %>%
  ungroup()

income_composition_by_location_parent <- df %>%
  mutate(parent_status = if_else(has_child_u18, "parent_u18", "not_parent_u18")) %>%
  group_by(mms_location, parent_status, income_quintile) %>%
  summarise(weight = sum(perwt, na.rm = TRUE), n = n(), .groups = "drop") %>%
  group_by(mms_location, parent_status) %>%
  mutate(share_within_location_parent = weight / sum(weight, na.rm = TRUE)) %>%
  ungroup()

gap_decomposition <- bind_rows(
  decompose_gap(df, "recent_birth"),
  decompose_gap(df, "has_child_u18"),
  decompose_gap(df, "childless_in_household"),
  decompose_gap(df, "mean_children_household")
)

reg_data <- df %>%
  filter(fertyr_observed, owner_label %in% c("owner", "renter")) %>%
  mutate(
    log_hhincome_std = as.numeric(scale(log_hhincome)),
    age_factor = factor(age),
    year_factor = factor(year),
    race_factor = factor(race),
    educ_factor = factor(educ),
    tenure_factor = factor(owner_label)
  )

recent_birth_reg <- feols(
  recent_birth ~ log_hhincome_std + I(log_hhincome_std^2) + college + married + tenure_factor + center | age + year + met2013,
  data = reg_data,
  weights = ~perwt,
  cluster = ~met2013
)

parent_reg <- feols(
  has_child_u18 ~ log_hhincome_std + I(log_hhincome_std^2) + college + married + tenure_factor + center | age + year + met2013,
  data = reg_data,
  weights = ~perwt,
  cluster = ~met2013
)

center_parent_reg <- feols(
  center ~ has_child_u18 * log_hhincome_std + has_child_u18 * college + tenure_factor | age + year + met2013,
  data = reg_data,
  weights = ~perwt,
  cluster = ~met2013
)

regression_table <- bind_rows(
  coef_table(recent_birth_reg, "recent_birth_income"),
  coef_table(parent_reg, "parent_u18_income"),
  coef_table(center_parent_reg, "center_parent_income")
)

write_csv(fertility_by_income, file.path(out_dir, "acs_fertility_by_income_quintile.csv"))
write_csv(fertility_by_income_location, file.path(out_dir, "acs_fertility_by_income_location.csv"))
write_csv(fertility_by_income_age, file.path(out_dir, "acs_fertility_by_income_age.csv"))
write_csv(fertility_by_income_tenure, file.path(out_dir, "acs_fertility_by_income_tenure.csv"))
write_csv(fertility_by_residual_income, file.path(out_dir, "acs_fertility_by_residual_income_quintile.csv"))
write_csv(income_distribution_by_fertility, file.path(out_dir, "acs_income_distribution_by_fertility_group.csv"))
write_csv(parity_income_distribution, file.path(out_dir, "acs_income_distribution_by_parity.csv"))
write_csv(income_composition_by_location_parent, file.path(out_dir, "acs_income_composition_by_location_parent.csv"))
write_csv(gap_decomposition, file.path(out_dir, "acs_center_periphery_fertility_income_decomposition.csv"))
write_csv(regression_table, file.path(out_dir, "acs_income_fertility_regressions.csv"))

p_recent <- ggplot(fertility_by_income, aes(x = income_quintile, y = recent_birth_rate)) +
  geom_line(linewidth = 0.8, color = "#2C5F8A") +
  geom_point(size = 2.5, color = "#2C5F8A") +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(x = "Household income quintile", y = "Recent birth rate", title = "Recent Births Across the Income Distribution") +
  theme_minimal(base_size = 12)

p_parent <- ggplot(fertility_by_income, aes(x = income_quintile)) +
  geom_line(aes(y = parent_u18_rate, color = "Parent with child <18"), linewidth = 0.8) +
  geom_point(aes(y = parent_u18_rate, color = "Parent with child <18"), size = 2.3) +
  geom_line(aes(y = childless_rate, color = "No own child in household"), linewidth = 0.8) +
  geom_point(aes(y = childless_rate, color = "No own child in household"), size = 2.3) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Household income quintile", y = "Share", color = "", title = "Parenthood and Childlessness by Income") +
  theme_minimal(base_size = 12)

p_loc <- ggplot(fertility_by_income_location, aes(x = income_quintile, y = parent_u18_rate, color = mms_location)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Household income quintile", y = "Parent with child <18", color = "", title = "Income-Fertility Gradient by Center/Periphery") +
  theme_minimal(base_size = 12)

p_dist <- income_distribution_by_fertility %>%
  filter(fertility_group %in% c("recent_birth", "childless")) %>%
  ggplot(aes(x = income_quintile, y = share_within_group, fill = fertility_group)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Household income quintile", y = "Share within group", fill = "", title = "Where Recent Births and Childlessness Sit in Income") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "acs_recent_birth_by_income.png"), p_recent, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_parent_childless_by_income.png"), p_parent, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_parent_by_income_location.png"), p_loc, width = 7.5, height = 4.5, dpi = 180)
ggsave(file.path(out_dir, "acs_income_distribution_birth_childless.png"), p_dist, width = 7.5, height = 4.5, dpi = 180)

get_coef_line <- function(model_name, term_name) {
  row <- regression_table %>% filter(model == model_name, term == term_name)
  if (nrow(row) == 0) {
    return("not estimated")
  }
  sprintf("%.4f (SE %.4f, p=%.3g)", row$estimate[1], row$se[1], row$p_value[1])
}

q1 <- fertility_by_income %>% filter(income_quintile == 1)
q5 <- fertility_by_income %>% filter(income_quintile == 5)
recent_dist <- income_distribution_by_fertility %>%
  filter(fertility_group == "recent_birth") %>%
  arrange(income_quintile)
childless_dist <- income_distribution_by_fertility %>%
  filter(fertility_group == "childless") %>%
  arrange(income_quintile)
decomp_recent <- gap_decomposition %>% filter(outcome == "recent_birth")
decomp_parent <- gap_decomposition %>% filter(outcome == "has_child_u18")

md <- c(
  "# ACS/MMS Income-Fertility Cross-Section",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  "- Source: ACS/IPUMS `extract27.dta` joined to the MMS center/periphery PUMA lookup.",
  "- Sample: women ages 22-45, non-group-quarter households, years 2012-2023, positive household income.",
  "- Fertility objects are ACS cross-sectional measures: recent birth from `fertyr == 2`, own children in household from `nchild`, and childlessness as no own child in household. They are not completed fertility.",
  sprintf("- Middle PUMAs are assigned to: `%s`.", middle_target),
  "- Income bins are within-year household-income quintiles. A residual-income quintile table additionally residualizes log household income on age, year, education, and marital status.",
  "",
  "## Main Cross-Section Facts",
  "",
  sprintf("- Matched women records: `%s`.", comma(nrow(df))),
  sprintf("- Recent birth rate, bottom income quintile: `%.3f`; top income quintile: `%.3f`.", q1$recent_birth_rate[1], q5$recent_birth_rate[1]),
  sprintf("- Parent-with-child-under-18 rate, bottom income quintile: `%.3f`; top income quintile: `%.3f`.", q1$parent_u18_rate[1], q5$parent_u18_rate[1]),
  sprintf("- Childless-in-household rate, bottom income quintile: `%.3f`; top income quintile: `%.3f`.", q1$childless_rate[1], q5$childless_rate[1]),
  sprintf("- Owner rate, bottom income quintile: `%.3f`; top income quintile: `%.3f`.", q1$owner_rate[1], q5$owner_rate[1]),
  sprintf("- Mean rooms, bottom income quintile: `%.3f`; top income quintile: `%.3f`.", q1$mean_rooms[1], q5$mean_rooms[1]),
  sprintf("- Share of recent births in income quintiles 1-5: %s.", paste(sprintf("%.3f", recent_dist$share_within_group), collapse = ", ")),
  sprintf("- Share of childless women in income quintiles 1-5: %s.", paste(sprintf("%.3f", childless_dist$share_within_group), collapse = ", ")),
  "",
  "## Center-Periphery Income Decomposition",
  "",
  sprintf("- Center-minus-periphery recent-birth gap: `%.4f`; within-income component `%.4f`; income-composition component `%.4f`.",
          decomp_recent$center_minus_periphery[1], decomp_recent$within_income_component[1], decomp_recent$income_composition_component[1]),
  sprintf("- Center-minus-periphery parent-with-child-under-18 gap: `%.4f`; within-income component `%.4f`; income-composition component `%.4f`.",
          decomp_parent$center_minus_periphery[1], decomp_parent$within_income_component[1], decomp_parent$income_composition_component[1]),
  "",
  "## Regression Reads",
  "",
  sprintf("- Recent birth slope on standardized log household income: `%s`.", get_coef_line("recent_birth_income", "log_hhincome_std")),
  sprintf("- Recent birth quadratic in standardized log household income: `%s`.", get_coef_line("recent_birth_income", "I(log_hhincome_std^2)")),
  sprintf("- Parent-with-child-under-18 slope on standardized log household income: `%s`.", get_coef_line("parent_u18_income", "log_hhincome_std")),
  sprintf("- Parent centrality interaction with standardized log household income: `%s`.", get_coef_line("center_parent_income", "has_child_u18TRUE:log_hhincome_std")),
  "",
  "## Outputs",
  "",
  "- `acs_fertility_by_income_quintile.csv`",
  "- `acs_fertility_by_income_location.csv`",
  "- `acs_fertility_by_income_age.csv`",
  "- `acs_fertility_by_income_tenure.csv`",
  "- `acs_fertility_by_residual_income_quintile.csv`",
  "- `acs_income_distribution_by_fertility_group.csv`",
  "- `acs_income_distribution_by_parity.csv`",
  "- `acs_income_composition_by_location_parent.csv`",
  "- `acs_center_periphery_fertility_income_decomposition.csv`",
  "- `acs_income_fertility_regressions.csv`",
  "- PNG figures for recent births, parenthood/childlessness, location gradients, and income distributions.",
  "",
  "## Interpretation Guardrails",
  "",
  "Household income is endogenous to marriage, labor supply, and fertility. The residual-income table is a type-composition diagnostic, not an instrument. ACS `nchild` measures own children in the household, not completed fertility."
)

writeLines(md, file.path(out_dir, "ACS_MMS_INCOME_FERTILITY_CROSS_SECTION.md"))
message("Wrote income-fertility cross-section packet to: ", out_dir)
