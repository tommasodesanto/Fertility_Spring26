# =============================================================================
# DETAILED MOVES ANALYSIS
# =============================================================================
# Breaking down:
# 1. Suburbanization: same-MSA vs different-MSA
# 2. Job changes for movers
# 3. Origin vs destination comparison (not just destination levels)
# =============================================================================

library(tidyverse)
library(haven)
library(fixest)
library(knitr)
library(scales)

OUTPUT_DIR <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/spatial_wage_analysis"

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 14)))

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading ACS data...\n")
acs <- read_dta("/Users/tommasodesanto/Desktop/Projects/Datasets/extract_15/acs_longer_larger.dta")

acs <- acs %>%
  filter(
    year >= 2010,
    age >= 25 & age <= 45
  ) %>%
  mutate(
    # Parent status
    new_parent = case_when(
      nchild == 0 ~ "Non-Parent",
      eldch <= 3 ~ "New Parent",
      TRUE ~ "Older Parent"
    ),

    # Mover status
    moved = migrate1d %in% c(23, 24, 25, 31, 32, 40),
    moved_diff_county_same_state = migrate1d %in% c(24, 25),
    moved_diff_state = migrate1d %in% c(31, 32),

    # Current location
    in_metro = metro %in% c(2, 3, 4),
    in_principal_city = metro == 2,
    in_suburb = metro == 3,  # metro, not in principal city

    # Origin metro type (1 year ago)
    origin_principal_city = migtype1 == 3,
    origin_suburb = migtype1 == 2,
    origin_in_metro = migtype1 %in% c(2, 3),
    origin_not_metro = migtype1 %in% c(1, 4, 5),

    # MSA identifier
    msa_id = ifelse(!is.na(met2013) & met2013 > 0, met2013, metarea),

    # Origin MSA (need to convert - migmet1 = metarea * 10)
    origin_msa = ifelse(!is.na(migmet1) & migmet1 > 0, migmet1 / 10, NA),

    # Did they change MSA?
    changed_msa = !is.na(origin_msa) & !is.na(msa_id) & (origin_msa != msa_id),
    same_msa = !is.na(origin_msa) & !is.na(msa_id) & (origin_msa == msa_id),

    # Work location
    work_in_diff_metro = !is.na(pwmetro) & pwmetro > 0 & pwmetro != metro,

    # Employment
    employed = empstat == 1,

    # Owner
    owner = ownershp == 1
  )

cat("Sample size:", nrow(acs), "\n")

# =============================================================================
# ANALYSIS 1: DECOMPOSE SUBURBANIZATION
# =============================================================================

cat("\n" , strrep("=", 60), "\n")
cat("ANALYSIS 1: DECOMPOSING SUBURBANIZATION\n")
cat(strrep("=", 60), "\n")

# Focus on movers who started in a metro area (principal city specifically)
city_origin <- acs %>%
  filter(
    moved,
    origin_principal_city,  # Started in principal city
    new_parent %in% c("Non-Parent", "New Parent")
  )

cat("\nMovers who started in principal city:", nrow(city_origin), "\n")

# Check how many have MSA change info
cat("With MSA change info:", sum(!is.na(city_origin$changed_msa)), "\n")

# Decomposition: Where do they end up?
decomp_all <- city_origin %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    # Current location
    pct_still_principal = weighted.mean(in_principal_city, perwt, na.rm = TRUE) * 100,
    pct_now_suburb = weighted.mean(in_suburb, perwt, na.rm = TRUE) * 100,
    pct_left_metro = weighted.mean(!in_metro, perwt, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n--- Where do Principal City Movers End Up? ---\n")
print(kable(decomp_all, digits = 1, format = "simple"))

# Now decompose by same vs different MSA (where we have data)
decomp_msa <- city_origin %>%
  filter(!is.na(changed_msa)) %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    # Same MSA moves
    pct_same_msa = weighted.mean(same_msa, perwt, na.rm = TRUE) * 100,
    pct_diff_msa = weighted.mean(changed_msa, perwt, na.rm = TRUE) * 100,

    # Among same-MSA: where?
    pct_same_msa_to_suburb = weighted.mean(same_msa & in_suburb, perwt, na.rm = TRUE) * 100,
    pct_same_msa_stay_city = weighted.mean(same_msa & in_principal_city, perwt, na.rm = TRUE) * 100,

    # Among diff-MSA: where?
    pct_diff_msa_to_suburb = weighted.mean(changed_msa & in_suburb, perwt, na.rm = TRUE) * 100,
    pct_diff_msa_to_city = weighted.mean(changed_msa & in_principal_city, perwt, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n--- MSA Change Decomposition (where data available) ---\n")
print(kable(decomp_msa, digits = 1, format = "simple"))

# Better breakdown
cat("\n--- DETAILED BREAKDOWN ---\n")
detailed_decomp <- city_origin %>%
  filter(!is.na(changed_msa)) %>%
  mutate(
    move_type = case_when(
      same_msa & in_principal_city ~ "Same MSA: City → City",
      same_msa & in_suburb ~ "Same MSA: City → Suburb",
      same_msa & !in_metro ~ "Same MSA: City → Non-metro",
      changed_msa & in_principal_city ~ "Diff MSA: → City",
      changed_msa & in_suburb ~ "Diff MSA: → Suburb",
      changed_msa & !in_metro ~ "Diff MSA: → Non-metro",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(new_parent, move_type) %>%
  summarise(n = sum(perwt), .groups = "drop") %>%
  group_by(new_parent) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = new_parent, values_from = pct, values_fill = 0)

print(kable(detailed_decomp, digits = 1, format = "simple"))

# Key comparison
cat("\n--- KEY: Suburbanization Rate by Move Type ---\n")
suburb_rates <- city_origin %>%
  filter(!is.na(changed_msa), in_metro) %>%  # Stay in metro
  group_by(new_parent, same_msa) %>%
  summarise(
    n = n(),
    pct_suburb = weighted.mean(in_suburb, perwt) * 100,
    .groups = "drop"
  ) %>%
  mutate(msa_type = ifelse(same_msa, "Same MSA", "Different MSA"))

print(kable(suburb_rates %>% select(new_parent, msa_type, pct_suburb, n), digits = 1, format = "simple"))

# =============================================================================
# ANALYSIS 2: JOB CHANGES
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("ANALYSIS 2: DID MOVERS CHANGE JOBS?\n")
cat(strrep("=", 60), "\n")

# Look at employment and work location for movers
movers <- acs %>%
  filter(
    moved,
    new_parent %in% c("Non-Parent", "New Parent"),
    age >= 25 & age <= 45
  )

job_outcomes <- movers %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    # Employment
    pct_employed = weighted.mean(employed, perwt, na.rm = TRUE) * 100,

    # Among employed: work in different metro than residence?
    pct_work_diff_metro = weighted.mean(work_in_diff_metro[employed], perwt[employed], na.rm = TRUE) * 100,

    # Income (proxy for job quality)
    mean_income = weighted.mean(inctot[inctot > 0 & inctot < 500000],
                                perwt[inctot > 0 & inctot < 500000], na.rm = TRUE),
    mean_wage = weighted.mean(incwage[incwage > 0 & incwage < 500000],
                              perwt[incwage > 0 & incwage < 500000], na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Employment Outcomes for Movers ---\n")
print(kable(job_outcomes, digits = 1, format = "simple"))

# Regression: income among movers
cat("\n--- Regression: Income ~ New Parent (among movers) ---\n")
reg_income <- feols(
  log(incwage) ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = movers %>% filter(incwage > 0, incwage < 500000),
  weights = ~perwt
)
print(summary(reg_income))

# =============================================================================
# ANALYSIS 3: ORIGIN VS DESTINATION (Better wage comparison)
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("ANALYSIS 3: ORIGIN VS DESTINATION MSA WAGES\n")
cat(strrep("=", 60), "\n")

# Compute MSA wage characteristics
msa_wages <- acs %>%
  filter(msa_id > 0, incwage > 0, incwage < 500000, employed) %>%
  group_by(msa_id, year) %>%
  summarise(
    msa_wage = weighted.mean(incwage, perwt, na.rm = TRUE),
    msa_rent = weighted.mean(rent[rent > 0], perwt[rent > 0], na.rm = TRUE),
    msa_n = n(),
    .groups = "drop"
  ) %>%
  filter(msa_n > 200)

# Get destination characteristics
movers_dest <- movers %>%
  filter(msa_id > 0) %>%
  left_join(msa_wages %>% rename(dest_wage = msa_wage, dest_rent = msa_rent),
            by = c("msa_id", "year"))

# Get origin characteristics (using origin_msa)
movers_both <- movers_dest %>%
  filter(!is.na(origin_msa)) %>%
  left_join(msa_wages %>%
              rename(origin_wage = msa_wage, origin_rent = msa_rent, origin_msa = msa_id),
            by = c("origin_msa", "year"))

cat("Movers with both origin and destination MSA data:",
    nrow(movers_both %>% filter(!is.na(dest_wage), !is.na(origin_wage))), "\n")

# Compute wage change
wage_change <- movers_both %>%
  filter(!is.na(dest_wage), !is.na(origin_wage)) %>%
  mutate(
    wage_diff = dest_wage - origin_wage,
    wage_pct_change = (dest_wage - origin_wage) / origin_wage * 100,
    moved_to_lower_wage = dest_wage < origin_wage,
    rent_diff = dest_rent - origin_rent,
    rent_pct_change = (dest_rent - origin_rent) / origin_rent * 100
  )

cat("\nSample with wage change:", nrow(wage_change), "\n")

# Summary by parent status
wage_change_summary <- wage_change %>%
  group_by(new_parent) %>%
  summarise(
    n = n(),
    # Origin
    mean_origin_wage = weighted.mean(origin_wage, perwt, na.rm = TRUE),
    mean_origin_rent = weighted.mean(origin_rent, perwt, na.rm = TRUE),
    # Destination
    mean_dest_wage = weighted.mean(dest_wage, perwt, na.rm = TRUE),
    mean_dest_rent = weighted.mean(dest_rent, perwt, na.rm = TRUE),
    # Change
    mean_wage_change = weighted.mean(wage_diff, perwt, na.rm = TRUE),
    mean_wage_pct_change = weighted.mean(wage_pct_change, perwt, na.rm = TRUE),
    pct_moved_lower_wage = weighted.mean(moved_to_lower_wage, perwt, na.rm = TRUE) * 100,
    mean_rent_change = weighted.mean(rent_diff, perwt, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- ORIGIN VS DESTINATION: Wage and Rent Changes ---\n")
print(kable(wage_change_summary, digits = 1, format = "simple"))

# Regression: Did they move to lower-wage MSA?
cat("\n--- Regression: P(Moved to Lower-Wage MSA) ~ New Parent ---\n")
reg_lower_wage <- feols(
  moved_to_lower_wage ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) + log(inctot + 1) | year,
  data = wage_change %>% filter(inctot > 0),
  weights = ~perwt
)
print(summary(reg_lower_wage))

# Regression: MSA wage percent change
cat("\n--- Regression: MSA Wage % Change ~ New Parent ---\n")
reg_wage_change <- feols(
  wage_pct_change ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) + log(inctot + 1) | year,
  data = wage_change %>% filter(inctot > 0, abs(wage_pct_change) < 100),  # Remove outliers
  weights = ~perwt
)
print(summary(reg_wage_change))

# =============================================================================
# VISUALIZATION
# =============================================================================

# Plot 1: Suburbanization decomposition
if (nrow(detailed_decomp) > 0) {
  plot_decomp <- detailed_decomp %>%
    pivot_longer(cols = -move_type, names_to = "parent_status", values_to = "pct") %>%
    filter(pct > 0) %>%
    ggplot(aes(x = reorder(move_type, pct), y = pct, fill = parent_status)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              position = position_dodge(width = 0.7), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
    labs(
      title = "Where Do Principal City Movers Go?",
      subtitle = "Breakdown by same vs different MSA",
      x = NULL,
      y = "Percent of Movers",
      fill = NULL
    ) +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, max(detailed_decomp[,-1], na.rm = TRUE) * 1.2))

  ggsave(file.path(OUTPUT_DIR, "suburbanization_decomposed.png"), plot_decomp,
         width = 10, height = 6, dpi = 300, bg = "white")
}

# Plot 2: Wage change comparison
if (nrow(wage_change_summary) > 0) {
  plot_wage <- wage_change_summary %>%
    select(new_parent, mean_origin_wage, mean_dest_wage) %>%
    pivot_longer(cols = -new_parent, names_to = "type", values_to = "wage") %>%
    mutate(type_label = ifelse(type == "mean_origin_wage", "Origin MSA", "Destination MSA")) %>%
    ggplot(aes(x = type_label, y = wage/1000, fill = new_parent, group = new_parent)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = sprintf("$%.1fK", wage/1000)),
              position = position_dodge(width = 0.6), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("Non-Parent" = "#2E86AB", "New Parent" = "#A23B72")) +
    scale_y_continuous(labels = dollar_format(prefix = "$", suffix = "K"),
                       limits = c(0, max(wage_change_summary$mean_dest_wage,
                                         wage_change_summary$mean_origin_wage)/1000 * 1.15)) +
    labs(
      title = "Origin vs Destination MSA Wages",
      subtitle = "Comparing where movers came from vs where they went",
      x = NULL,
      y = "Average MSA Wage ($K)",
      fill = NULL
    ) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "origin_vs_destination_wages.png"), plot_wage,
         width = 8, height = 5, dpi = 300, bg = "white")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\n1. SUBURBANIZATION DECOMPOSITION:\n")
if (exists("detailed_decomp")) {
  cat("   See table above for same-MSA vs different-MSA breakdown\n")
}

cat("\n2. WAGE CHANGES (Origin → Destination):\n")
if (nrow(wage_change_summary) > 0) {
  np_change <- wage_change_summary$mean_wage_change[wage_change_summary$new_parent == "Non-Parent"]
  p_change <- wage_change_summary$mean_wage_change[wage_change_summary$new_parent == "New Parent"]
  cat(sprintf("   Non-Parents: MSA wage change = $%.0f\n", np_change))
  cat(sprintf("   New Parents: MSA wage change = $%.0f\n", p_change))
  cat(sprintf("   Difference: $%.0f\n", p_change - np_change))

  np_pct <- wage_change_summary$pct_moved_lower_wage[wage_change_summary$new_parent == "Non-Parent"]
  p_pct <- wage_change_summary$pct_moved_lower_wage[wage_change_summary$new_parent == "New Parent"]
  cat(sprintf("\n   Percent moved to lower-wage MSA:\n"))
  cat(sprintf("   Non-Parents: %.1f%%\n", np_pct))
  cat(sprintf("   New Parents: %.1f%%\n", p_pct))
}

cat("\n=== Files saved to:", OUTPUT_DIR, "===\n")
