# =============================================================================
# SPATIAL WAGE ANALYSIS V2: ALTERNATIVE APPROACHES
# =============================================================================
# The MSA-to-MSA analysis was limited by migmet1 availability (2005-2011 only).
# This script tries different approaches with better data coverage.
#
# Analyses:
#   1. Within-metro moves: Principal city → suburbs
#   2. Individual outcomes for movers: rooms, rent burden, income
#   3. Exit from metro areas entirely
#   4. Destination characteristics conditional on moving
# =============================================================================

library(tidyverse)
library(haven)
library(fixest)
library(knitr)
library(scales)

# Set paths
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

# Focus on recent years with good migration data
acs <- acs %>%
  filter(
    year >= 2010,
    age >= 22 & age <= 50
  ) %>%
  mutate(
    # Parent status
    new_parent = case_when(
      nchild == 0 ~ "Non-Parent",
      eldch <= 3 ~ "New Parent",
      TRUE ~ "Older Parent"
    ),

    # Mover status (1-year)
    moved = migrate1d %in% c(23, 24, 25, 31, 32, 40),
    moved_across_state = migrate1d %in% c(31, 32),
    moved_within_state = migrate1d %in% c(23, 24, 25),

    # Metro status
    in_metro = metro %in% c(2, 3, 4),
    in_principal_city = metro == 2,
    in_metro_not_principal = metro == 3,
    not_in_metro = metro == 1,

    # Origin metro status (1 year ago)
    origin_in_metro = migtype1 %in% c(2, 3),
    origin_principal_city = migtype1 == 3,
    origin_not_principal = migtype1 == 2,
    origin_not_metro = migtype1 %in% c(1, 4, 5),

    # Rent burden
    rent_burden = ifelse(rent > 0 & hhincome > 0 & hhincome < 9999999,
                         (rent * 12) / hhincome, NA),

    # Owner
    owner = ownershp == 1
  )

cat("Sample size after filters:", nrow(acs), "\n")
cat("By parent status:\n")
print(table(acs$new_parent))

# =============================================================================
# ANALYSIS 1: WITHIN-METRO MOVES (Principal City <-> Suburbs)
# =============================================================================

cat("\n========================================\n")
cat("ANALYSIS 1: WITHIN-METRO TRANSITIONS\n")
cat("========================================\n")

# Filter to people in metro areas with valid migration type
metro_movers <- acs %>%
  filter(
    in_metro,
    !is.na(migtype1),
    migtype1 > 0,
    new_parent %in% c("Non-Parent", "New Parent")
  )

cat("\nSample with metro migration info:", nrow(metro_movers), "\n")

# Create transition categories
metro_movers <- metro_movers %>%
  mutate(
    transition = case_when(
      origin_principal_city & in_principal_city ~ "City → City",
      origin_principal_city & in_metro_not_principal ~ "City → Suburb",
      origin_not_principal & in_principal_city ~ "Suburb → City",
      origin_not_principal & in_metro_not_principal ~ "Suburb → Suburb",
      origin_not_metro & in_metro ~ "Non-metro → Metro",
      TRUE ~ "Other"
    )
  ) %>%
  filter(transition != "Other")

# Transition probabilities by parent status
transition_probs <- metro_movers %>%
  group_by(new_parent, transition) %>%
  summarise(n = n(), weighted_n = sum(perwt), .groups = "drop") %>%
  group_by(new_parent) %>%
  mutate(pct = weighted_n / sum(weighted_n) * 100) %>%
  ungroup()

cat("\n--- Within-Metro Transitions by Parent Status ---\n")
transition_probs %>%
  select(new_parent, transition, pct) %>%
  pivot_wider(names_from = new_parent, values_from = pct, values_fill = 0) %>%
  kable(digits = 1, format = "simple") %>%
  print()

# Key comparison: City → Suburb rate
city_to_suburb <- transition_probs %>%
  filter(transition == "City → Suburb") %>%
  select(new_parent, pct)

cat("\n--- KEY: City → Suburb Rate ---\n")
print(city_to_suburb)

# Statistical test
suburbanization_test <- metro_movers %>%
  filter(origin_principal_city) %>%
  mutate(went_suburb = in_metro_not_principal)

cat("\n--- Chi-squared: Suburbanization by Parent Status ---\n")
print(chisq.test(table(suburbanization_test$new_parent, suburbanization_test$went_suburb)))

# Regression
cat("\n--- Regression: P(City → Suburb) ~ New Parent ---\n")
reg_suburb <- feols(
  in_metro_not_principal ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = metro_movers %>% filter(origin_principal_city),
  weights = ~perwt
)
print(summary(reg_suburb))

# Plot
p_transitions <- transition_probs %>%
  ggplot(aes(x = transition, y = pct, fill = new_parent)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
  coord_flip() +
  labs(
    title = "Within-Metro Transitions by Parenthood Status",
    subtitle = "Among movers currently in metro areas",
    x = NULL,
    y = "Percent",
    fill = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "within_metro_transitions.png"), p_transitions,
       width = 9, height = 5, dpi = 300, bg = "white")

# =============================================================================
# ANALYSIS 2: EXIT FROM METRO AREAS
# =============================================================================

cat("\n========================================\n")
cat("ANALYSIS 2: EXIT FROM METRO AREAS\n")
cat("========================================\n")

# People who were in metro 1 year ago
metro_origin <- acs %>%
  filter(
    origin_in_metro,
    new_parent %in% c("Non-Parent", "New Parent")
  )

cat("\nPeople originating from metro:", nrow(metro_origin), "\n")

# Did they leave?
exit_rates <- metro_origin %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    weighted_n = sum(perwt),
    pct_left_metro = weighted.mean(not_in_metro, perwt) * 100,
    pct_stayed_metro = weighted.mean(in_metro, perwt) * 100,
    .groups = "drop"
  )

cat("\n--- Metro Exit Rates by Parent Status ---\n")
print(kable(exit_rates, digits = 2, format = "simple"))

# Regression
cat("\n--- Regression: P(Left Metro) ~ New Parent ---\n")
reg_exit <- feols(
  not_in_metro ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = metro_origin,
  weights = ~perwt
)
print(summary(reg_exit))

# =============================================================================
# ANALYSIS 3: INDIVIDUAL OUTCOMES FOR MOVERS
# =============================================================================

cat("\n========================================\n")
cat("ANALYSIS 3: INDIVIDUAL OUTCOMES FOR MOVERS\n")
cat("========================================\n")

# Focus on people who moved
movers_only <- acs %>%
  filter(
    moved,
    new_parent %in% c("Non-Parent", "New Parent")
  )

cat("\nMovers sample:", nrow(movers_only), "\n")

# Outcomes by parent status
mover_outcomes <- movers_only %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    # Housing
    mean_rooms = weighted.mean(rooms, perwt, na.rm = TRUE),
    mean_bedrooms = weighted.mean(bedrooms, perwt, na.rm = TRUE),
    pct_owner = weighted.mean(owner, perwt, na.rm = TRUE) * 100,
    # Rent burden (for renters)
    mean_rent_burden = weighted.mean(rent_burden[!owner], perwt[!owner], na.rm = TRUE),
    median_rent_burden = median(rent_burden[!owner], na.rm = TRUE),
    # Location
    pct_principal_city = weighted.mean(in_principal_city, perwt, na.rm = TRUE) * 100,
    pct_not_metro = weighted.mean(not_in_metro, perwt, na.rm = TRUE) * 100,
    # Income
    mean_income = weighted.mean(inctot[inctot > 0 & inctot < 9999999],
                                perwt[inctot > 0 & inctot < 9999999], na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Outcomes for Movers by Parent Status ---\n")
print(kable(mover_outcomes, digits = 2, format = "simple"))

# Regression: Rooms conditional on moving
cat("\n--- Regression: Rooms ~ New Parent (among movers) ---\n")
reg_rooms <- feols(
  rooms ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) + log(inctot + 1) | year,
  data = movers_only %>% filter(inctot > 0, inctot < 9999999),
  weights = ~perwt
)
print(summary(reg_rooms))

# Regression: Rent burden conditional on moving (renters only)
cat("\n--- Regression: Rent Burden ~ New Parent (among renter movers) ---\n")
reg_burden <- feols(
  rent_burden ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = movers_only %>% filter(!owner, !is.na(rent_burden), rent_burden < 1),
  weights = ~perwt
)
print(summary(reg_burden))

# Regression: Principal city conditional on moving
cat("\n--- Regression: P(Principal City) ~ New Parent (among movers) ---\n")
reg_principal <- feols(
  in_principal_city ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = movers_only %>% filter(in_metro),
  weights = ~perwt
)
print(summary(reg_principal))

# =============================================================================
# ANALYSIS 4: DESTINATION CHARACTERISTICS BY PARENT STATUS
# =============================================================================

cat("\n========================================\n")
cat("ANALYSIS 4: WHERE DO MOVERS END UP?\n")
cat("========================================\n")

# Compute MSA-level characteristics from ALL data (not just movers)
msa_chars <- acs %>%
  filter(metarea > 0, inctot > 0, inctot < 9999999) %>%
  group_by(metarea) %>%
  summarise(
    msa_mean_income = weighted.mean(inctot, perwt, na.rm = TRUE),
    msa_mean_rent = weighted.mean(rent[rent > 0], perwt[rent > 0], na.rm = TRUE),
    msa_median_rooms = median(rooms, na.rm = TRUE),
    msa_pct_owner = weighted.mean(ownershp == 1, perwt, na.rm = TRUE) * 100,
    msa_pop = sum(perwt),
    .groups = "drop"
  ) %>%
  filter(msa_pop > 50000) %>%  # Only reasonably sized MSAs
  mutate(
    msa_income_tercile = ntile(msa_mean_income, 3),
    msa_rent_tercile = ntile(msa_mean_rent, 3)
  )

# Merge to movers
movers_with_msa <- movers_only %>%
  filter(metarea > 0) %>%
  left_join(msa_chars, by = "metarea")

cat("\nMovers with MSA characteristics:", nrow(movers_with_msa %>% filter(!is.na(msa_mean_income))), "\n")

# Destination characteristics
dest_chars <- movers_with_msa %>%
  filter(!is.na(msa_mean_income)) %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    mean_dest_income = weighted.mean(msa_mean_income, perwt, na.rm = TRUE),
    mean_dest_rent = weighted.mean(msa_mean_rent, perwt, na.rm = TRUE),
    pct_to_low_income_msa = weighted.mean(msa_income_tercile == 1, perwt, na.rm = TRUE) * 100,
    pct_to_high_income_msa = weighted.mean(msa_income_tercile == 3, perwt, na.rm = TRUE) * 100,
    pct_to_low_rent_msa = weighted.mean(msa_rent_tercile == 1, perwt, na.rm = TRUE) * 100,
    pct_to_high_rent_msa = weighted.mean(msa_rent_tercile == 3, perwt, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n--- Destination MSA Characteristics by Parent Status ---\n")
print(kable(dest_chars, digits = 2, format = "simple"))

# Regression: destination MSA income
cat("\n--- Regression: Dest MSA Income ~ New Parent ---\n")
reg_dest_income <- feols(
  log(msa_mean_income) ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = movers_with_msa %>% filter(!is.na(msa_mean_income)),
  weights = ~perwt
)
print(summary(reg_dest_income))

# Regression: destination MSA rent
cat("\n--- Regression: Dest MSA Rent ~ New Parent ---\n")
reg_dest_rent <- feols(
  log(msa_mean_rent) ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = movers_with_msa %>% filter(!is.na(msa_mean_rent)),
  weights = ~perwt
)
print(summary(reg_dest_rent))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n========================================\n")
cat("SUMMARY OF FINDINGS\n")
cat("========================================\n")

cat("\n1. WITHIN-METRO TRANSITIONS:\n")
city_suburb_np <- transition_probs %>% filter(new_parent == "Non-Parent", transition == "City → Suburb") %>% pull(pct)
city_suburb_p <- transition_probs %>% filter(new_parent == "New Parent", transition == "City → Suburb") %>% pull(pct)
cat(sprintf("   City → Suburb rate: Non-Parents %.1f%%, New Parents %.1f%% (diff: %.1f pp)\n",
            city_suburb_np, city_suburb_p, city_suburb_p - city_suburb_np))

cat("\n2. METRO EXIT:\n")
exit_np <- exit_rates$pct_left_metro[exit_rates$new_parent == "Non-Parent"]
exit_p <- exit_rates$pct_left_metro[exit_rates$new_parent == "New Parent"]
cat(sprintf("   Left metro rate: Non-Parents %.2f%%, New Parents %.2f%% (diff: %.2f pp)\n",
            exit_np, exit_p, exit_p - exit_np))

cat("\n3. MOVERS' INDIVIDUAL OUTCOMES:\n")
rooms_np <- mover_outcomes$mean_rooms[mover_outcomes$new_parent == "Non-Parent"]
rooms_p <- mover_outcomes$mean_rooms[mover_outcomes$new_parent == "New Parent"]
cat(sprintf("   Mean rooms: Non-Parents %.2f, New Parents %.2f (diff: %.2f)\n",
            rooms_np, rooms_p, rooms_p - rooms_np))

principal_np <- mover_outcomes$pct_principal_city[mover_outcomes$new_parent == "Non-Parent"]
principal_p <- mover_outcomes$pct_principal_city[mover_outcomes$new_parent == "New Parent"]
cat(sprintf("   In principal city: Non-Parents %.1f%%, New Parents %.1f%% (diff: %.1f pp)\n",
            principal_np, principal_p, principal_p - principal_np))

cat("\n4. DESTINATION MSA CHARACTERISTICS:\n")
dest_inc_np <- dest_chars$mean_dest_income[dest_chars$new_parent == "Non-Parent"]
dest_inc_p <- dest_chars$mean_dest_income[dest_chars$new_parent == "New Parent"]
cat(sprintf("   Mean dest MSA income: Non-Parents $%.0f, New Parents $%.0f (diff: $%.0f)\n",
            dest_inc_np, dest_inc_p, dest_inc_p - dest_inc_np))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output saved to:", OUTPUT_DIR, "\n")
