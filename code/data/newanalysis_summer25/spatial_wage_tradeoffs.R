# =============================================================================
# SPATIAL WAGE TRADEOFFS: DO PARENTS ACCEPT WAGE CUTS FOR CHEAPER HOUSING?
# =============================================================================
# This script explores whether new parents systematically move to lower-wage
# cities, accepting wage cuts in exchange for lower housing costs.
#
# Analyses:
#   1. City tier transition matrices by parenthood status
#   2. Direct wage change regressions for MSA-to-MSA movers
#   3. Superstar city outflow analysis
#   4. Wage-housing tradeoff visualization
# =============================================================================

library(tidyverse)
library(tidycensus)
library(haven)
library(fixest)
library(knitr)
library(scales)

# Set paths
OUTPUT_DIR <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/spatial_wage_analysis"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Set theme
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, color = "gray40")))

# =============================================================================
# PART 0: LOAD AND EXPLORE DATA
# =============================================================================

cat("Loading ACS microdata...\n")
ACS_PATH <- "/Users/tommasodesanto/Desktop/Projects/Datasets/extract_15/acs_longer_larger.dta"
acs_raw <- read_dta(ACS_PATH)

cat("Dataset dimensions:", nrow(acs_raw), "rows,", ncol(acs_raw), "cols\n")

# Check year range and migration variable availability
cat("\nYear distribution:\n")
print(table(acs_raw$year))

cat("\nMigration variable availability by year:\n")
acs_raw %>%
  group_by(year) %>%
  summarise(
    n = n(),
    pct_migmet1_available = mean(!is.na(migmet1) & migmet1 > 0) * 100,
    pct_migrate1d_available = mean(!is.na(migrate1d)) * 100
  ) %>%
  print(n = 30)

# =============================================================================
# PART 1: GET MSA CHARACTERISTICS FROM CENSUS API
# =============================================================================

cat("\n\nFetching MSA characteristics from Census API...\n")

get_msa_characteristics <- function(year = 2019) {

  econ_vars <- c(
    "B19013_001",  # Median household income
    "B25077_001",  # Median home value
    "B25064_001",  # Median gross rent
    "B01003_001",  # Total population
    "B20002_001"   # Median earnings
  )

  msa_data <- get_acs(
    geography = "metropolitan statistical area/micropolitan statistical area",
    variables = econ_vars,
    year = year,
    survey = "acs5",
    output = "wide",
    cache_table = TRUE
  ) %>%
    transmute(
      GEOID = GEOID,
      msa_name = str_extract(NAME, "^[^,]+"),
      median_income = B19013_001E,
      median_home_value = B25077_001E,
      median_rent = B25064_001E,
      population = B01003_001E,
      median_earnings = B20002_001E,
      price_to_income = median_home_value / median_income,
      rent_to_income = (median_rent * 12) / median_income
    ) %>%
    filter(
      !is.na(median_income),
      !is.na(median_rent),
      population > 100000
    )

  return(msa_data)
}

msa_chars <- get_msa_characteristics(2019)
cat("Retrieved characteristics for", nrow(msa_chars), "MSAs\n")

# Create city tiers based on wages
msa_chars <- msa_chars %>%
  mutate(
    wage_tercile = ntile(median_earnings, 3),
    city_tier = case_when(
      wage_tercile == 1 ~ "Peripheral",
      wage_tercile == 2 ~ "Secondary",
      wage_tercile == 3 ~ "Superstar"
    ),
    city_tier = factor(city_tier, levels = c("Peripheral", "Secondary", "Superstar")),
    log_earnings = log(median_earnings),
    log_rent = log(median_rent)
  )

cat("\nCity tier distribution:\n")
print(table(msa_chars$city_tier))

cat("\nMSA characteristics by tier:\n")
msa_chars %>%
  group_by(city_tier) %>%
  summarise(
    n = n(),
    mean_earnings = mean(median_earnings, na.rm = TRUE),
    mean_rent = mean(median_rent, na.rm = TRUE),
    mean_price_to_income = mean(price_to_income, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# PART 2: CREATE CROSSWALK FOR IPUMS METRO CODES TO CBSA
# =============================================================================

# IPUMS metarea codes need to be mapped to CBSA codes
# migmet1 = metarea * 10 based on our exploration

cat("\n\nCreating metro code crosswalk...\n")

# We'll work with metarea codes directly and create our own wage/rent measures
# by aggregating the ACS microdata itself

# First, let's compute MSA-level characteristics from the microdata
msa_from_microdata <- acs_raw %>%
  filter(
    year >= 2010,  # Recent years with good coverage
    metarea > 0,
    !is.na(inctot),
    inctot > 0,
    inctot < 9999999,
    age >= 25, age <= 64  # Prime working age
  ) %>%
  group_by(metarea) %>%
  summarise(
    n_obs = n(),
    mean_income = weighted.mean(inctot, perwt, na.rm = TRUE),
    median_income_approx = median(inctot, na.rm = TRUE),
    mean_rent = weighted.mean(rent[rent > 0], perwt[rent > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_obs >= 500)  # Only MSAs with enough observations

cat("Computed characteristics for", nrow(msa_from_microdata), "metro areas from microdata\n")

# Create tiers based on microdata income
msa_from_microdata <- msa_from_microdata %>%
  mutate(
    wage_tercile = ntile(mean_income, 3),
    city_tier = case_when(
      wage_tercile == 1 ~ "Peripheral",
      wage_tercile == 2 ~ "Secondary",
      wage_tercile == 3 ~ "Superstar"
    ),
    city_tier = factor(city_tier, levels = c("Peripheral", "Secondary", "Superstar")),
    log_income = log(mean_income),
    log_rent = log(mean_rent)
  )

cat("\nCity tiers from microdata:\n")
print(table(msa_from_microdata$city_tier))

# =============================================================================
# PART 3: PREPARE INDIVIDUAL MOVER DATA
# =============================================================================

cat("\n\nPreparing individual mover data...\n")

# Focus on recent years with migration data
movers <- acs_raw %>%
  filter(
    year >= 2010,
    age >= 22 & age <= 50,
    # Moved across MSAs (not same house, not within same metro)
    migrate1d %in% c(24, 31, 32),  # Between PUMAs within state, or between states
    metarea > 0,
    migmet1 > 0
  ) %>%
  mutate(
    # Convert migmet1 to metarea scale (migmet1 = metarea * 10)
    origin_metarea = migmet1 / 10,
    dest_metarea = metarea,

    # Parent status
    new_parent = case_when(
      nchild == 0 ~ "Non-Parent",
      eldch <= 3 ~ "New Parent",
      TRUE ~ "Older Parent"
    )
  ) %>%
  filter(new_parent %in% c("Non-Parent", "New Parent"))

cat("Initial movers sample:", nrow(movers), "\n")

# Merge origin MSA characteristics
movers <- movers %>%
  left_join(
    msa_from_microdata %>%
      select(metarea,
             origin_income = mean_income,
             origin_rent = mean_rent,
             origin_tier = city_tier,
             origin_log_income = log_income,
             origin_log_rent = log_rent),
    by = c("origin_metarea" = "metarea")
  )

# Merge destination MSA characteristics
movers <- movers %>%
  left_join(
    msa_from_microdata %>%
      select(metarea,
             dest_income = mean_income,
             dest_rent = mean_rent,
             dest_tier = city_tier,
             dest_log_income = log_income,
             dest_log_rent = log_rent),
    by = c("dest_metarea" = "metarea")
  )

# Keep only movers where we have both origin and destination data
movers <- movers %>%
  filter(!is.na(origin_tier), !is.na(dest_tier))

cat("Movers with origin and destination data:", nrow(movers), "\n")
cat("By parent status:\n")
print(table(movers$new_parent))

# Compute changes
movers <- movers %>%
  mutate(
    delta_log_income = dest_log_income - origin_log_income,
    delta_log_rent = dest_log_rent - origin_log_rent,
    pct_income_change = (exp(delta_log_income) - 1) * 100,
    pct_rent_change = (exp(delta_log_rent) - 1) * 100,

    # Direction of move
    moved_down = (origin_tier == "Superstar" & dest_tier %in% c("Secondary", "Peripheral")) |
                 (origin_tier == "Secondary" & dest_tier == "Peripheral"),
    moved_up = (origin_tier == "Peripheral" & dest_tier %in% c("Secondary", "Superstar")) |
               (origin_tier == "Secondary" & dest_tier == "Superstar")
  )

# =============================================================================
# PART 4: CITY TIER TRANSITION ANALYSIS
# =============================================================================

cat("\n")
cat("========================================\n")
cat("ANALYSIS 1: CITY TIER TRANSITIONS\n")
cat("========================================\n")

# Transition matrix by parent status
tier_transitions <- movers %>%
  group_by(new_parent, origin_tier, dest_tier) %>%
  summarise(
    n = n(),
    weighted_n = sum(perwt),
    .groups = "drop"
  ) %>%
  group_by(new_parent, origin_tier) %>%
  mutate(pct = weighted_n / sum(weighted_n) * 100) %>%
  ungroup()

# Display transition matrices
cat("\n--- Transition Matrix: NON-PARENTS ---\n")
cat("(Rows = Origin, Columns = Destination)\n\n")
tier_transitions %>%
  filter(new_parent == "Non-Parent") %>%
  select(origin_tier, dest_tier, pct) %>%
  pivot_wider(names_from = dest_tier, values_from = pct, values_fill = 0) %>%
  kable(digits = 1, format = "simple") %>%
  print()

cat("\n--- Transition Matrix: NEW PARENTS ---\n")
cat("(Rows = Origin, Columns = Destination)\n\n")
tier_transitions %>%
  filter(new_parent == "New Parent") %>%
  select(origin_tier, dest_tier, pct) %>%
  pivot_wider(names_from = dest_tier, values_from = pct, values_fill = 0) %>%
  kable(digits = 1, format = "simple") %>%
  print()

# Summary: probability of moving down
downward_summary <- movers %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    n_weighted = sum(perwt),
    pct_moved_down = weighted.mean(moved_down, perwt) * 100,
    pct_moved_up = weighted.mean(moved_up, perwt) * 100,
    pct_same_tier = 100 - weighted.mean(moved_down, perwt) * 100 - weighted.mean(moved_up, perwt) * 100,
    .groups = "drop"
  )

cat("\n--- Summary: Direction of Move ---\n")
print(kable(downward_summary, digits = 1, format = "simple"))

# Statistical test
cat("\n--- Chi-squared Test: Moved Down by Parent Status ---\n")
chisq_result <- chisq.test(table(movers$new_parent, movers$moved_down))
print(chisq_result)

# Save
write_csv(tier_transitions, file.path(OUTPUT_DIR, "tier_transitions.csv"))

# Plot
p_transitions <- tier_transitions %>%
  ggplot(aes(x = dest_tier, y = pct, fill = new_parent)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~origin_tier, labeller = labeller(origin_tier = function(x) paste("From:", x))) +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  labs(
    title = "City Tier Transitions by Parenthood Status",
    subtitle = "Where do movers from each tier end up?",
    x = "Destination Tier",
    y = "Percent of Movers",
    fill = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "tier_transitions.png"), p_transitions,
       width = 10, height = 6, dpi = 300, bg = "white")

# =============================================================================
# PART 5: WAGE CHANGE ANALYSIS
# =============================================================================

cat("\n")
cat("========================================\n")
cat("ANALYSIS 2: WAGE/INCOME CHANGES\n")
cat("========================================\n")

# Summary statistics
income_change_summary <- movers %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    mean_pct_income_change = weighted.mean(pct_income_change, perwt, na.rm = TRUE),
    median_pct_income_change = median(pct_income_change, na.rm = TRUE),
    mean_pct_rent_change = weighted.mean(pct_rent_change, perwt, na.rm = TRUE),
    median_pct_rent_change = median(pct_rent_change, na.rm = TRUE),
    share_income_cut = weighted.mean(pct_income_change < 0, perwt, na.rm = TRUE) * 100,
    share_rent_cut = weighted.mean(pct_rent_change < 0, perwt, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n--- Income and Rent Changes by Parent Status ---\n")
print(kable(income_change_summary, digits = 2, format = "simple"))

# Regressions
cat("\n--- Regression: Delta Log Income ~ Parent Status ---\n")

reg1 <- feols(delta_log_income ~ i(new_parent, ref = "Non-Parent"),
              data = movers, weights = ~perwt)

reg2 <- feols(delta_log_income ~ i(new_parent, ref = "Non-Parent") +
                age + I(age^2) + factor(sex) + factor(educd),
              data = movers, weights = ~perwt)

reg3 <- feols(delta_log_income ~ i(new_parent, ref = "Non-Parent") +
                age + I(age^2) + factor(sex) + factor(educd) | year,
              data = movers, weights = ~perwt)

reg4 <- feols(delta_log_income ~ i(new_parent, ref = "Non-Parent") +
                age + I(age^2) + factor(sex) + factor(educd) | year + origin_tier,
              data = movers, weights = ~perwt)

cat("\nModel 1: Unconditional\n")
print(summary(reg1, se = "hetero"))

cat("\nModel 2: With demographics\n")
print(summary(reg2, se = "hetero"))

cat("\nModel 3: + Year FE\n")
print(summary(reg3, se = "cluster"))

cat("\nModel 4: + Year FE + Origin Tier FE\n")
print(summary(reg4, se = "cluster"))

# Rent change regression
cat("\n--- Regression: Delta Log Rent ~ Parent Status ---\n")
reg_rent <- feols(delta_log_rent ~ i(new_parent, ref = "Non-Parent") +
                    age + I(age^2) + factor(sex) + factor(educd) | year + origin_tier,
                  data = movers, weights = ~perwt)
print(summary(reg_rent, se = "cluster"))

# =============================================================================
# PART 6: SUPERSTAR CITY OUTFLOWS
# =============================================================================

cat("\n")
cat("========================================\n")
cat("ANALYSIS 3: SUPERSTAR OUTFLOWS\n")
cat("========================================\n")

superstar_leavers <- movers %>%
  filter(origin_tier == "Superstar")

cat("People leaving Superstar cities:", nrow(superstar_leavers), "\n")

if (nrow(superstar_leavers) > 100) {

  # Where do they go?
  superstar_dest <- superstar_leavers %>%
    group_by(new_parent, dest_tier) %>%
    summarise(
      n = n(),
      weighted_n = sum(perwt),
      .groups = "drop"
    ) %>%
    group_by(new_parent) %>%
    mutate(pct = weighted_n / sum(weighted_n) * 100)

  cat("\n--- Destinations of Superstar Leavers ---\n")
  superstar_dest %>%
    select(new_parent, dest_tier, pct) %>%
    pivot_wider(names_from = dest_tier, values_from = pct, values_fill = 0) %>%
    kable(digits = 1, format = "simple") %>%
    print()

  # Income/rent changes
  superstar_changes <- superstar_leavers %>%
    group_by(new_parent) %>%
    summarise(
      n = n(),
      mean_income_change = weighted.mean(pct_income_change, perwt, na.rm = TRUE),
      mean_rent_change = weighted.mean(pct_rent_change, perwt, na.rm = TRUE),
      share_to_peripheral = weighted.mean(dest_tier == "Peripheral", perwt) * 100,
      .groups = "drop"
    )

  cat("\n--- Income/Rent Changes for Superstar Leavers ---\n")
  print(kable(superstar_changes, digits = 2, format = "simple"))

  # Plot
  p_superstar <- superstar_dest %>%
    ggplot(aes(x = dest_tier, y = pct, fill = new_parent)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
    labs(
      title = "Where Do Superstar City Leavers Go?",
      subtitle = "Destination tier by parenthood status",
      x = "Destination Tier",
      y = "Percent of Leavers",
      fill = NULL
    ) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "superstar_outflows.png"), p_superstar,
         width = 8, height = 5, dpi = 300, bg = "white")
}

# =============================================================================
# PART 7: WAGE-HOUSING TRADEOFF VISUALIZATION
# =============================================================================

cat("\n")
cat("========================================\n")
cat("ANALYSIS 4: WAGE-HOUSING TRADEOFF\n")
cat("========================================\n")

# Quadrant analysis
quadrant_summary <- movers %>%
  mutate(
    quadrant = case_when(
      pct_income_change >= 0 & pct_rent_change >= 0 ~ "Higher income, Higher rent",
      pct_income_change >= 0 & pct_rent_change < 0 ~ "Higher income, Lower rent",
      pct_income_change < 0 & pct_rent_change >= 0 ~ "Lower income, Higher rent",
      pct_income_change < 0 & pct_rent_change < 0 ~ "Lower income, Lower rent"
    )
  ) %>%
  group_by(new_parent, quadrant) %>%
  summarise(
    n = n(),
    weighted_n = sum(perwt),
    .groups = "drop"
  ) %>%
  group_by(new_parent) %>%
  mutate(pct = weighted_n / sum(weighted_n) * 100)

cat("\n--- Wage-Housing Tradeoff Quadrants ---\n")
quadrant_summary %>%
  select(new_parent, quadrant, pct) %>%
  pivot_wider(names_from = new_parent, values_from = pct, values_fill = 0) %>%
  kable(digits = 1, format = "simple") %>%
  print()

# Key finding
cat("\n--- KEY FINDING: Share accepting 'Lower Income, Lower Rent' ---\n")
quadrant_summary %>%
  filter(quadrant == "Lower income, Lower rent") %>%
  select(new_parent, pct) %>%
  print()

# Scatter plot
p_tradeoff <- movers %>%
  filter(abs(pct_income_change) < 50, abs(pct_rent_change) < 50) %>%
  ggplot(aes(x = pct_income_change, y = pct_rent_change, color = new_parent)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  labs(
    title = "The Wage-Housing Tradeoff",
    subtitle = "Change in local income vs. change in local rents for movers",
    x = "Change in MSA Mean Income (%)",
    y = "Change in MSA Mean Rent (%)",
    color = NULL
  ) +
  annotate("text", x = -25, y = -35, label = "Lower income,\nLower rent",
           color = "gray40", size = 3) +
  annotate("text", x = 25, y = 35, label = "Higher income,\nHigher rent",
           color = "gray40", size = 3) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "wage_housing_tradeoff.png"), p_tradeoff,
       width = 8, height = 7, dpi = 300, bg = "white")

# Density plot
p_density <- movers %>%
  filter(abs(pct_income_change) < 50) %>%
  ggplot(aes(x = pct_income_change, fill = new_parent, color = new_parent)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  scale_color_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  labs(
    title = "Distribution of MSA Income Changes for Movers",
    subtitle = "Do new parents show a leftward shift (more 'wage cuts')?",
    x = "Change in MSA Mean Income (%)",
    y = "Density",
    fill = NULL, color = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "income_change_density.png"), p_density,
       width = 8, height = 5, dpi = 300, bg = "white")

# =============================================================================
# PART 8: SUMMARY REPORT
# =============================================================================

cat("\n")
cat("========================================\n")
cat("SUMMARY REPORT\n")
cat("========================================\n")

cat("\nSAMPLE:\n")
cat("- Total movers (cross-MSA, age 22-50, 2010+):", nrow(movers), "\n")
cat("- Non-Parents:", sum(movers$new_parent == "Non-Parent"), "\n")
cat("- New Parents:", sum(movers$new_parent == "New Parent"), "\n")

cat("\nKEY FINDINGS:\n")

# Finding 1: Downward mobility
down_nonparent <- downward_summary$pct_moved_down[downward_summary$new_parent == "Non-Parent"]
down_newparent <- downward_summary$pct_moved_down[downward_summary$new_parent == "New Parent"]
cat(sprintf("1. Probability of moving to lower-wage tier:\n"))
cat(sprintf("   - Non-Parents: %.1f%%\n", down_nonparent))
cat(sprintf("   - New Parents: %.1f%%\n", down_newparent))
cat(sprintf("   - Difference: %.1f pp\n", down_newparent - down_nonparent))

# Finding 2: Income change
inc_nonparent <- income_change_summary$mean_pct_income_change[income_change_summary$new_parent == "Non-Parent"]
inc_newparent <- income_change_summary$mean_pct_income_change[income_change_summary$new_parent == "New Parent"]
cat(sprintf("\n2. Mean MSA income change when moving:\n"))
cat(sprintf("   - Non-Parents: %.2f%%\n", inc_nonparent))
cat(sprintf("   - New Parents: %.2f%%\n", inc_newparent))
cat(sprintf("   - Difference: %.2f pp\n", inc_newparent - inc_nonparent))

# Finding 3: Regression coefficient
coef_newparent <- coef(reg4)["new_parent::New Parent"]
cat(sprintf("\n3. Regression (with Year + Origin Tier FE):\n"))
cat(sprintf("   - New Parent coefficient on delta log income: %.4f\n", coef_newparent))
cat(sprintf("   - Interpretation: New parents move to MSAs with %.1f%% lower income\n",
            (exp(coef_newparent) - 1) * 100))

cat("\nOUTPUT FILES:\n")
cat("- ", file.path(OUTPUT_DIR, "tier_transitions.csv"), "\n")
cat("- ", file.path(OUTPUT_DIR, "tier_transitions.png"), "\n")
cat("- ", file.path(OUTPUT_DIR, "superstar_outflows.png"), "\n")
cat("- ", file.path(OUTPUT_DIR, "wage_housing_tradeoff.png"), "\n")
cat("- ", file.path(OUTPUT_DIR, "income_change_density.png"), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
