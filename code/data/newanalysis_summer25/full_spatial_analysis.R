# =============================================================================
# FULL SPATIAL ANALYSIS: MSA Changes + Suburbanization
# =============================================================================
# 1. MSA wage changes (origin vs destination)
# 2. Suburbanization decomposed by same-MSA vs different-MSA
# =============================================================================

library(tidyverse)
library(haven)
library(fixest)
library(knitr)
library(scales)
library(patchwork)

OUTPUT_DIR <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/spatial_wage_analysis"

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 14)))

cat("Loading ACS data...\n")
acs <- read_dta("/Users/tommasodesanto/Desktop/Projects/Datasets/extract_15/acs_longer_larger.dta")

# =============================================================================
# STEP 1: Create county-to-MSA crosswalk
# =============================================================================

cat("\n=== Creating county-to-MSA crosswalk ===\n")

county_to_msa <- acs %>%
  filter(
    year >= 2010, year <= 2021,
    !is.na(statefip), statefip > 0,
    !is.na(countyfip), countyfip > 0,
    !is.na(met2013), met2013 > 0
  ) %>%
  group_by(statefip, countyfip) %>%
  summarise(
    msa = as.numeric(names(which.max(table(met2013)))),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n > 100)

cat("Counties mapped to MSAs:", nrow(county_to_msa), "\n")

# =============================================================================
# STEP 2: Compute MSA characteristics
# =============================================================================

cat("\n=== Computing MSA characteristics ===\n")

msa_chars <- acs %>%
  filter(
    year >= 2010, year <= 2021,
    !is.na(met2013), met2013 > 0,
    age >= 25, age <= 55,
    incwage > 0, incwage < 500000
  ) %>%
  group_by(met2013) %>%
  summarise(
    msa_mean_wage = weighted.mean(incwage, perwt, na.rm = TRUE),
    msa_mean_rent = weighted.mean(rent[rent > 0], perwt[rent > 0], na.rm = TRUE),
    msa_pop = sum(perwt),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(n_obs > 1000)

cat("MSAs with characteristics:", nrow(msa_chars), "\n")

# =============================================================================
# STEP 3: Process movers
# =============================================================================

cat("\n=== Processing movers ===\n")

movers <- acs %>%
  filter(
    year >= 2012, year <= 2021,
    migrate1d %in% c(23, 24, 25, 31, 32),  # All movers (including same county)
    age >= 25, age <= 45
  ) %>%
  mutate(
    new_parent = case_when(
      nchild == 0 ~ "Non-Parent",
      eldch <= 3 ~ "New Parent",
      TRUE ~ "Older Parent"
    ),

    # Current location
    in_metro = metro %in% c(2, 3, 4),
    in_principal_city = metro == 2,
    in_suburb = metro == 3,

    # Origin location (from migtype1)
    # 1 = not in metro, 2 = metro not principal, 3 = principal city, 4 = metro not principal, 5 = abroad
    origin_in_metro = migtype1 %in% c(2, 3, 4),
    origin_principal_city = migtype1 == 3,
    origin_suburb = migtype1 %in% c(2, 4),

    # Origin county for MSA mapping
    origin_statefip = as.numeric(migplac1),
    origin_countyfip = as.numeric(migcounty1),

    # County changer?
    county_changer = migrate1d %in% c(24, 25, 31, 32)
  )

cat("Total movers:", nrow(movers), "\n")

# =============================================================================
# STEP 4: Map to MSAs and identify MSA changes
# =============================================================================

cat("\n=== Mapping origin county to MSA ===\n")

# Map origin county to origin MSA
movers <- movers %>%
  left_join(
    county_to_msa %>%
      rename(origin_msa = msa, origin_statefip = statefip, origin_countyfip = countyfip),
    by = c("origin_statefip", "origin_countyfip")
  ) %>%
  mutate(
    # Determine if changed MSA
    changed_msa = case_when(
      is.na(origin_msa) | is.na(met2013) ~ NA,
      origin_msa != met2013 ~ TRUE,
      TRUE ~ FALSE
    ),
    same_msa = !is.na(changed_msa) & !changed_msa
  )

cat("Movers with MSA mapping:", sum(!is.na(movers$origin_msa)), "\n")

# =============================================================================
# PART A: MSA WAGE CHANGES
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PART A: MSA WAGE CHANGES (Origin vs Destination)\n")
cat(strrep("=", 60), "\n")

# Merge MSA characteristics
msa_changers <- movers %>%
  filter(
    county_changer,
    !is.na(origin_msa),
    !is.na(met2013), met2013 > 0,
    changed_msa == TRUE,
    new_parent %in% c("Non-Parent", "New Parent")
  ) %>%
  left_join(msa_chars %>% select(met2013, dest_wage = msa_mean_wage, dest_rent = msa_mean_rent),
            by = "met2013") %>%
  left_join(msa_chars %>% select(met2013, origin_wage = msa_mean_wage, origin_rent = msa_mean_rent) %>%
              rename(origin_msa = met2013),
            by = "origin_msa") %>%
  filter(!is.na(dest_wage), !is.na(origin_wage)) %>%
  mutate(
    wage_change = dest_wage - origin_wage,
    wage_pct_change = (dest_wage - origin_wage) / origin_wage * 100,
    moved_to_lower_wage = dest_wage < origin_wage
  )

cat("MSA changers with wage data:", nrow(msa_changers), "\n")

# Summary
msa_summary <- msa_changers %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    weighted_n = sum(perwt),
    mean_origin_wage = weighted.mean(origin_wage, perwt),
    mean_dest_wage = weighted.mean(dest_wage, perwt),
    mean_wage_change = weighted.mean(wage_change, perwt),
    pct_lower_wage = weighted.mean(moved_to_lower_wage, perwt) * 100,
    mean_origin_rent = weighted.mean(origin_rent, perwt),
    mean_dest_rent = weighted.mean(dest_rent, perwt),
    .groups = "drop"
  )

cat("\n--- MSA WAGE CHANGES ---\n")
print(kable(msa_summary, digits = 1, format = "simple"))

# Regression
cat("\n--- Regression: P(Moved to Lower-Wage MSA) ---\n")
reg_wage <- feols(
  moved_to_lower_wage ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = msa_changers,
  weights = ~perwt
)
cat("New Parent coefficient:", round(coef(reg_wage)["new_parent::New Parent"], 4),
    "p-value:", round(fixest::pvalue(reg_wage)["new_parent::New Parent"], 4), "\n")

# =============================================================================
# PART B: SUBURBANIZATION - DECOMPOSED
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PART B: SUBURBANIZATION (Same-MSA vs Different-MSA)\n")
cat(strrep("=", 60), "\n")

# Focus on people who started in principal city
city_starters <- movers %>%
  filter(
    origin_principal_city,  # Started in principal city
    in_metro,               # Currently in metro
    new_parent %in% c("Non-Parent", "New Parent")
  )

cat("\nPeople who started in principal city:", nrow(city_starters), "\n")
cat("With MSA change info:", sum(!is.na(city_starters$changed_msa)), "\n")

# Decomposition
suburb_decomp <- city_starters %>%
  filter(!is.na(changed_msa)) %>%
  mutate(
    move_category = case_when(
      same_msa & in_principal_city ~ "Same MSA: Stay in City",
      same_msa & in_suburb ~ "Same MSA: City → Suburb",
      changed_msa & in_principal_city ~ "Diff MSA: → City",
      changed_msa & in_suburb ~ "Diff MSA: → Suburb",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(new_parent, move_category) %>%
  summarise(n = sum(perwt), .groups = "drop") %>%
  group_by(new_parent) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

cat("\n--- SUBURBANIZATION DECOMPOSITION (from Principal City) ---\n")
suburb_wide <- suburb_decomp %>%
  select(-n) %>%
  pivot_wider(names_from = new_parent, values_from = pct, values_fill = 0)
print(kable(suburb_wide, digits = 1, format = "simple"))

# Key rates
cat("\n--- KEY SUBURBANIZATION RATES ---\n")
key_rates <- suburb_decomp %>%
  filter(move_category %in% c("Same MSA: City → Suburb", "Diff MSA: → Suburb")) %>%
  group_by(new_parent) %>%
  summarise(
    same_msa_suburb = sum(pct[move_category == "Same MSA: City → Suburb"]),
    diff_msa_suburb = sum(pct[move_category == "Diff MSA: → Suburb"]),
    total_suburb = sum(pct),
    .groups = "drop"
  )
print(kable(key_rates, digits = 1, format = "simple"))

# Regression: P(End up in suburb | started in city)
cat("\n--- Regression: P(Suburb) ~ New Parent (among city starters) ---\n")
reg_suburb <- feols(
  in_suburb ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = city_starters,
  weights = ~perwt
)
cat("New Parent coefficient:", round(coef(reg_suburb)["new_parent::New Parent"], 4),
    "p-value:", round(fixest::pvalue(reg_suburb)["new_parent::New Parent"], 4), "\n")

# Regression with MSA change interaction
cat("\n--- Regression: P(Suburb) ~ New Parent * Changed MSA ---\n")
reg_suburb_int <- feols(
  in_suburb ~ i(new_parent, ref = "Non-Parent") * changed_msa +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = city_starters %>% filter(!is.na(changed_msa)),
  weights = ~perwt
)
print(summary(reg_suburb_int, keep = c("new_parent", "changed_msa")))

# =============================================================================
# PART C: REVERSE DIRECTION - Suburb starters
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("PART C: SUBURB STARTERS - Where do they go?\n")
cat(strrep("=", 60), "\n")

suburb_starters <- movers %>%
  filter(
    origin_suburb,
    in_metro,
    new_parent %in% c("Non-Parent", "New Parent")
  )

cat("People who started in suburb:", nrow(suburb_starters), "\n")

suburb_to_city <- suburb_starters %>%
  filter(!is.na(changed_msa)) %>%
  mutate(
    move_category = case_when(
      same_msa & in_suburb ~ "Same MSA: Stay in Suburb",
      same_msa & in_principal_city ~ "Same MSA: Suburb → City",
      changed_msa & in_suburb ~ "Diff MSA: → Suburb",
      changed_msa & in_principal_city ~ "Diff MSA: → City",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(new_parent, move_category) %>%
  summarise(n = sum(perwt), .groups = "drop") %>%
  group_by(new_parent) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

cat("\n--- WHERE DO SUBURB STARTERS GO? ---\n")
suburb_start_wide <- suburb_to_city %>%
  select(-n) %>%
  pivot_wider(names_from = new_parent, values_from = pct, values_fill = 0)
print(kable(suburb_start_wide, digits = 1, format = "simple"))

# =============================================================================
# VISUALIZATIONS
# =============================================================================

cat("\n=== Creating visualizations ===\n")

# Plot 1: MSA Wage Change
p1 <- msa_summary %>%
  select(new_parent, mean_origin_wage, mean_dest_wage) %>%
  pivot_longer(-new_parent, names_to = "type", values_to = "wage") %>%
  mutate(type = ifelse(type == "mean_origin_wage", "Origin MSA", "Destination MSA")) %>%
  ggplot(aes(x = type, y = wage/1000, fill = new_parent)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("$%.1fK", wage/1000)),
            position = position_dodge(0.6), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  scale_y_continuous(labels = dollar_format(suffix = "K"),
                     limits = c(0, max(msa_summary$mean_dest_wage)/1000 * 1.15)) +
  labs(title = "MSA Wage: Origin vs Destination",
       subtitle = "Among people who changed MSA",
       x = NULL, y = "Average MSA Wage ($K)", fill = NULL) +
  theme(legend.position = "bottom")

# Plot 2: Suburbanization decomposition
p2 <- suburb_decomp %>%
  filter(pct > 1) %>%
  ggplot(aes(x = reorder(move_category, pct), y = pct, fill = new_parent)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(0.7), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  labs(title = "Where Do City Movers Go?",
       subtitle = "Decomposed by same vs different MSA",
       x = NULL, y = "Percent", fill = NULL) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, max(suburb_decomp$pct) * 1.15))

# Combined
p_combined <- p1 + p2 +
  plot_annotation(
    title = "Spatial Sorting: New Parents vs Non-Parents",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave(file.path(OUTPUT_DIR, "spatial_sorting_combined.png"), p_combined,
       width = 14, height = 6, dpi = 300, bg = "white")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\n1. MSA WAGE CHANGES (among MSA changers):\n")
np <- msa_summary %>% filter(new_parent == "Non-Parent")
p <- msa_summary %>% filter(new_parent == "New Parent")
cat(sprintf("   Non-Parents: $%.0f → $%.0f (change: $%.0f)\n",
            np$mean_origin_wage, np$mean_dest_wage, np$mean_wage_change))
cat(sprintf("   New Parents: $%.0f → $%.0f (change: $%.0f)\n",
            p$mean_origin_wage, p$mean_dest_wage, p$mean_wage_change))
cat(sprintf("   DIFFERENCE: $%.0f\n", p$mean_wage_change - np$mean_wage_change))
cat(sprintf("   P(lower wage): Non-Parents %.1f%%, New Parents %.1f%% (diff: %.1f pp)\n",
            np$pct_lower_wage, p$pct_lower_wage, p$pct_lower_wage - np$pct_lower_wage))

cat("\n2. SUBURBANIZATION (among city starters):\n")
np_sub <- key_rates %>% filter(new_parent == "Non-Parent")
p_sub <- key_rates %>% filter(new_parent == "New Parent")
cat(sprintf("   Non-Parents: Same-MSA suburb %.1f%%, Diff-MSA suburb %.1f%%, Total %.1f%%\n",
            np_sub$same_msa_suburb, np_sub$diff_msa_suburb, np_sub$total_suburb))
cat(sprintf("   New Parents: Same-MSA suburb %.1f%%, Diff-MSA suburb %.1f%%, Total %.1f%%\n",
            p_sub$same_msa_suburb, p_sub$diff_msa_suburb, p_sub$total_suburb))

cat("\n=== Output saved to:", OUTPUT_DIR, "===\n")
