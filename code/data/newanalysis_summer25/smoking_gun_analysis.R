# =============================================================================
# SMOKING GUN ANALYSIS: Parents Trade Wages for Cheaper Housing
# =============================================================================
# This analysis directly shows that parents who move systematically
# end up in places with:
#   1. Lower average wages
#   2. Lower housing costs (rents)
#   3. More housing space
# =============================================================================

library(tidyverse)
library(haven)
library(fixest)
library(knitr)
library(scales)

OUTPUT_DIR <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/spatial_wage_analysis"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 14)))

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading ACS data...\n")
acs <- read_dta("/Users/tommasodesanto/Desktop/Projects/Datasets/extract_15/acs_longer_larger.dta")

# Focus on prime working age, recent years
acs <- acs %>%
  filter(
    year >= 2010,
    age >= 25 & age <= 45,
    inctot > 5000 & inctot < 500000  # Working people with reasonable income
  ) %>%
  mutate(
    # Parent status - clean definition
    new_parent = case_when(
      nchild == 0 ~ "Non-Parent",
      eldch <= 3 ~ "New Parent",
      TRUE ~ "Older Parent"
    ),

    # Mover status (1-year migration)
    moved = migrate1d %in% c(23, 24, 25, 31, 32, 40),

    # Metro status
    in_metro = metro %in% c(2, 3, 4),
    in_principal_city = metro == 2,

    # Housing outcomes
    owner = ownershp == 1,

    # Origin metro type
    origin_principal_city = migtype1 == 3,
    origin_suburb = migtype1 == 2,
    origin_in_metro = migtype1 %in% c(2, 3)
  )

cat("Sample size:", nrow(acs), "\n")

# =============================================================================
# COMPUTE MSA-LEVEL CHARACTERISTICS
# =============================================================================

cat("\n=== Computing MSA characteristics ===\n")

# Use met2013 if available, else metarea
acs <- acs %>%
  mutate(msa_id = ifelse(!is.na(met2013) & met2013 > 0, met2013, metarea))

msa_chars <- acs %>%
  filter(msa_id > 0) %>%
  group_by(msa_id, year) %>%
  summarise(
    # Labor market
    msa_mean_wage = weighted.mean(incwage[incwage > 0 & incwage < 500000],
                                   perwt[incwage > 0 & incwage < 500000], na.rm = TRUE),
    msa_median_wage = median(incwage[incwage > 0 & incwage < 500000], na.rm = TRUE),
    msa_p75_wage = quantile(incwage[incwage > 0 & incwage < 500000], 0.75, na.rm = TRUE),

    # Housing costs
    msa_mean_rent = weighted.mean(rent[rent > 0], perwt[rent > 0], na.rm = TRUE),
    msa_median_rent = median(rent[rent > 0], na.rm = TRUE),
    msa_mean_value = weighted.mean(valueh[valueh > 0 & valueh < 9999998],
                                    perwt[valueh > 0 & valueh < 9999998], na.rm = TRUE),

    # Housing size
    msa_mean_rooms = weighted.mean(rooms, perwt, na.rm = TRUE),

    # Population
    msa_pop = sum(perwt),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(msa_pop > 100000, n_obs > 500)  # Only decent-sized MSAs

cat("MSAs with characteristics:", n_distinct(msa_chars$msa_id), "\n")

# Create terciles based on pooled characteristics
msa_terciles <- msa_chars %>%
  group_by(year) %>%
  mutate(
    wage_tercile = ntile(msa_mean_wage, 3),
    rent_tercile = ntile(msa_mean_rent, 3),
    wage_label = case_when(
      wage_tercile == 1 ~ "Low-Wage MSA",
      wage_tercile == 2 ~ "Mid-Wage MSA",
      wage_tercile == 3 ~ "High-Wage MSA"
    ),
    rent_label = case_when(
      rent_tercile == 1 ~ "Low-Cost MSA",
      rent_tercile == 2 ~ "Mid-Cost MSA",
      rent_tercile == 3 ~ "High-Cost MSA"
    )
  ) %>%
  ungroup()

# =============================================================================
# ANALYSIS: WHERE DO MOVERS END UP?
# =============================================================================

cat("\n=== SMOKING GUN: Destination Characteristics ===\n")

# Merge MSA characteristics to individuals
movers <- acs %>%
  filter(
    moved,
    new_parent %in% c("Non-Parent", "New Parent"),
    msa_id > 0
  ) %>%
  left_join(msa_terciles, by = c("msa_id", "year"))

cat("Movers with MSA data:", nrow(movers %>% filter(!is.na(msa_mean_wage))), "\n")

# Key comparison: Destination MSA characteristics
dest_summary <- movers %>%
  filter(!is.na(msa_mean_wage)) %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    # Weighted means of destination MSA characteristics
    dest_mean_wage = weighted.mean(msa_mean_wage, perwt, na.rm = TRUE),
    dest_mean_rent = weighted.mean(msa_mean_rent, perwt, na.rm = TRUE),
    dest_mean_value = weighted.mean(msa_mean_value, perwt, na.rm = TRUE),
    dest_mean_rooms = weighted.mean(msa_mean_rooms, perwt, na.rm = TRUE),

    # Distribution across terciles
    pct_low_wage = weighted.mean(wage_tercile == 1, perwt, na.rm = TRUE) * 100,
    pct_high_wage = weighted.mean(wage_tercile == 3, perwt, na.rm = TRUE) * 100,
    pct_low_rent = weighted.mean(rent_tercile == 1, perwt, na.rm = TRUE) * 100,
    pct_high_rent = weighted.mean(rent_tercile == 3, perwt, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n--- DESTINATION MSA CHARACTERISTICS BY PARENT STATUS ---\n")
print(kable(dest_summary, digits = 2, format = "simple"))

# Compute differences
diff_wage <- dest_summary$dest_mean_wage[dest_summary$new_parent == "New Parent"] -
             dest_summary$dest_mean_wage[dest_summary$new_parent == "Non-Parent"]
diff_rent <- dest_summary$dest_mean_rent[dest_summary$new_parent == "New Parent"] -
             dest_summary$dest_mean_rent[dest_summary$new_parent == "Non-Parent"]

cat("\n=== KEY FINDING ===\n")
cat(sprintf("New parents move to MSAs with:\n"))
cat(sprintf("  - $%.0f LOWER average wages\n", -diff_wage))
cat(sprintf("  - $%.0f LOWER average rents\n", -diff_rent))

# =============================================================================
# REGRESSION: DESTINATION MSA WAGE CONTROLLING FOR CHARACTERISTICS
# =============================================================================

cat("\n=== REGRESSION: Destination MSA Wage ===\n")

reg_data <- movers %>%
  filter(!is.na(msa_mean_wage), !is.na(educd))

reg_dest_wage <- feols(
  log(msa_mean_wage) ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) + log(inctot) | year,
  data = reg_data,
  weights = ~perwt
)
print(summary(reg_dest_wage))

cat("\n=== REGRESSION: Destination MSA Rent ===\n")
reg_dest_rent <- feols(
  log(msa_mean_rent) ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) + log(inctot) | year,
  data = reg_data,
  weights = ~perwt
)
print(summary(reg_dest_rent))

# =============================================================================
# VISUALIZATION 1: BAR CHART OF DESTINATION CHARACTERISTICS
# =============================================================================

# Create comparison data for plotting
plot_data <- dest_summary %>%
  select(new_parent, dest_mean_wage, dest_mean_rent) %>%
  pivot_longer(cols = starts_with("dest_"), names_to = "metric", values_to = "value") %>%
  mutate(
    metric_label = case_when(
      metric == "dest_mean_wage" ~ "Avg Wage in\nDestination MSA",
      metric == "dest_mean_rent" ~ "Avg Rent in\nDestination MSA"
    )
  )

# Separate plots for wage and rent
p_wage <- dest_summary %>%
  ggplot(aes(x = new_parent, y = dest_mean_wage / 1000, fill = new_parent)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("$%.1fK", dest_mean_wage/1000)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  scale_y_continuous(labels = dollar_format(prefix = "$", suffix = "K"),
                     limits = c(0, max(dest_summary$dest_mean_wage/1000) * 1.15)) +
  labs(
    title = "Destination MSA Average Wage",
    subtitle = "Among movers aged 25-45",
    x = NULL,
    y = "Average Wage ($K)",
    fill = NULL
  ) +
  theme(legend.position = "none")

p_rent <- dest_summary %>%
  ggplot(aes(x = new_parent, y = dest_mean_rent, fill = new_parent)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("$%.0f", dest_mean_rent)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  scale_y_continuous(labels = dollar_format(),
                     limits = c(0, max(dest_summary$dest_mean_rent) * 1.15)) +
  labs(
    title = "Destination MSA Average Rent",
    subtitle = "Among movers aged 25-45",
    x = NULL,
    y = "Average Monthly Rent ($)",
    fill = NULL
  ) +
  theme(legend.position = "none")

# Combined figure
library(patchwork)
p_combined <- p_wage + p_rent +
  plot_annotation(
    title = "The Wage-Housing Tradeoff: Where Do New Parents Move?",
    subtitle = "New parents who move choose lower-wage, lower-cost metro areas",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

ggsave(file.path(OUTPUT_DIR, "smoking_gun_destination_comparison.png"), p_combined,
       width = 10, height = 5, dpi = 300, bg = "white")

# =============================================================================
# VISUALIZATION 2: TERCILE DISTRIBUTION
# =============================================================================

tercile_data <- movers %>%
  filter(!is.na(wage_tercile)) %>%
  group_by(new_parent, wage_label) %>%
  summarise(n = sum(perwt), .groups = "drop") %>%
  group_by(new_parent) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_tercile <- tercile_data %>%
  mutate(wage_label = factor(wage_label, levels = c("Low-Wage MSA", "Mid-Wage MSA", "High-Wage MSA"))) %>%
  ggplot(aes(x = wage_label, y = pct, fill = new_parent)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  labs(
    title = "Distribution of Movers Across MSA Wage Terciles",
    subtitle = "New parents disproportionately move to lower-wage metros",
    x = "Destination MSA Wage Level",
    y = "Percent of Movers",
    fill = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "smoking_gun_wage_terciles.png"), p_tercile,
       width = 8, height = 5, dpi = 300, bg = "white")

# =============================================================================
# ANALYSIS: CITY VS SUBURB - DETAILED
# =============================================================================

cat("\n=== CITY VS SUBURB ANALYSIS ===\n")

# Among movers currently in metro areas
metro_movers <- acs %>%
  filter(
    moved,
    in_metro,
    new_parent %in% c("Non-Parent", "New Parent")
  )

city_suburb_summary <- metro_movers %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    pct_principal_city = weighted.mean(in_principal_city, perwt) * 100,
    pct_suburb = weighted.mean(metro == 3, perwt) * 100,
    .groups = "drop"
  )

cat("\n--- Where Movers End Up (City vs Suburb) ---\n")
print(kable(city_suburb_summary, digits = 2, format = "simple"))

# Regression
cat("\n--- Regression: P(Principal City) ~ New Parent ---\n")
reg_city <- feols(
  in_principal_city ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) + log(inctot) | year,
  data = metro_movers %>% filter(inctot > 0),
  weights = ~perwt
)
print(summary(reg_city))

# Bar chart
p_city <- city_suburb_summary %>%
  select(new_parent, pct_principal_city, pct_suburb) %>%
  pivot_longer(cols = starts_with("pct_"), names_to = "location", values_to = "pct") %>%
  mutate(location_label = ifelse(location == "pct_principal_city", "Principal City", "Suburb")) %>%
  ggplot(aes(x = location_label, y = pct, fill = new_parent)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  labs(
    title = "Suburbanization: Where Do Movers End Up?",
    subtitle = "New parents much less likely to be in principal city",
    x = "Metro Location",
    y = "Percent of Movers",
    fill = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "smoking_gun_city_suburb.png"), p_city,
       width = 7, height = 5, dpi = 300, bg = "white")

# =============================================================================
# SUMMARY TABLE FOR PAPER
# =============================================================================

cat("\n" %+% strrep("=", 60) %+% "\n")
cat("SMOKING GUN SUMMARY\n")
cat(strrep("=", 60) %+% "\n")

# Key statistics
np_wage <- dest_summary$dest_mean_wage[dest_summary$new_parent == "Non-Parent"]
p_wage <- dest_summary$dest_mean_wage[dest_summary$new_parent == "New Parent"]
np_rent <- dest_summary$dest_mean_rent[dest_summary$new_parent == "Non-Parent"]
p_rent <- dest_summary$dest_mean_rent[dest_summary$new_parent == "New Parent"]

np_city <- city_suburb_summary$pct_principal_city[city_suburb_summary$new_parent == "Non-Parent"]
p_city <- city_suburb_summary$pct_principal_city[city_suburb_summary$new_parent == "New Parent"]

cat("\n1. DESTINATION MSA WAGES:\n")
cat(sprintf("   Non-Parents move to MSAs with avg wage: $%.0f\n", np_wage))
cat(sprintf("   New Parents move to MSAs with avg wage: $%.0f\n", p_wage))
cat(sprintf("   Difference: -$%.0f (%.1f%% lower)\n", np_wage - p_wage, (np_wage - p_wage)/np_wage * 100))

cat("\n2. DESTINATION MSA RENTS:\n")
cat(sprintf("   Non-Parents move to MSAs with avg rent: $%.0f\n", np_rent))
cat(sprintf("   New Parents move to MSAs with avg rent:  $%.0f\n", p_rent))
cat(sprintf("   Difference: -$%.0f (%.1f%% lower)\n", np_rent - p_rent, (np_rent - p_rent)/np_rent * 100))

cat("\n3. CITY VS SUBURB:\n")
cat(sprintf("   Non-Parents in principal city: %.1f%%\n", np_city))
cat(sprintf("   New Parents in principal city: %.1f%%\n", p_city))
cat(sprintf("   Difference: %.1f pp\n", p_city - np_city))

cat("\n4. REGRESSION COEFFICIENTS (log destination MSA wage):\n")
cat(sprintf("   New Parent effect: %.4f (%.1f%% lower wages)\n",
            coef(reg_dest_wage)["new_parent::New Parent"],
            (1 - exp(coef(reg_dest_wage)["new_parent::New Parent"])) * 100))

cat("\n=== Files saved to:", OUTPUT_DIR, "===\n")
