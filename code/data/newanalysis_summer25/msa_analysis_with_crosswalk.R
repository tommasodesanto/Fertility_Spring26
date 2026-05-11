# =============================================================================
# MSA ANALYSIS WITH COUNTY-TO-MSA CROSSWALK
# =============================================================================
# Map origin county to origin MSA, then compare origin vs destination
# =============================================================================

library(tidyverse)
library(haven)
library(fixest)
library(knitr)

OUTPUT_DIR <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/spatial_wage_analysis"

cat("Loading ACS data...\n")
acs <- read_dta("/Users/tommasodesanto/Desktop/Projects/Datasets/extract_15/acs_longer_larger.dta")

# =============================================================================
# STEP 1: Create county-to-MSA crosswalk from the data itself
# =============================================================================
# We can infer county-to-MSA mapping from the data: for each county,
# what MSA are people in?

cat("\n=== Creating county-to-MSA crosswalk ===\n")

# Build crosswalk: for each state-county combo, find the most common MSA
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
  filter(n > 100)  # Only keep counties with enough obs

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
    msa_median_wage = median(incwage, na.rm = TRUE),
    msa_mean_rent = weighted.mean(rent[rent > 0], perwt[rent > 0], na.rm = TRUE),
    msa_pop = sum(perwt),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(n_obs > 1000)

cat("MSAs with characteristics:", nrow(msa_chars), "\n")

# =============================================================================
# STEP 3: Identify movers with origin county info
# =============================================================================

cat("\n=== Processing movers ===\n")

movers <- acs %>%
  filter(
    year >= 2012, year <= 2021,  # After migmet1 stops
    migrate1d %in% c(24, 25, 31, 32),  # County changers
    age >= 25, age <= 45,
    !is.na(migcounty1), migcounty1 > 0,
    !is.na(migplac1), migplac1 > 0, migplac1 < 100,  # US states only
    !is.na(met2013), met2013 > 0  # Currently in an MSA
  ) %>%
  mutate(
    new_parent = case_when(
      nchild == 0 ~ "Non-Parent",
      eldch <= 3 ~ "New Parent",
      TRUE ~ "Older Parent"
    ),
    # Create origin state-county FIPS
    origin_statefip = as.numeric(migplac1),
    origin_countyfip = as.numeric(migcounty1)
  )

cat("Movers with origin county:", nrow(movers), "\n")

# =============================================================================
# STEP 4: Map origin county to origin MSA
# =============================================================================

cat("\n=== Mapping origin county to origin MSA ===\n")

movers_mapped <- movers %>%
  left_join(
    county_to_msa %>%
      rename(origin_msa = msa, origin_statefip = statefip, origin_countyfip = countyfip),
    by = c("origin_statefip", "origin_countyfip")
  )

cat("Movers with origin MSA mapped:", sum(!is.na(movers_mapped$origin_msa)), "\n")

# =============================================================================
# STEP 5: Merge MSA characteristics for both origin and destination
# =============================================================================

movers_full <- movers_mapped %>%
  filter(!is.na(origin_msa)) %>%
  # Destination MSA characteristics
  left_join(
    msa_chars %>% rename(dest_msa = met2013, dest_wage = msa_mean_wage,
                         dest_rent = msa_mean_rent),
    by = c("met2013" = "dest_msa")
  ) %>%
  # Origin MSA characteristics
  left_join(
    msa_chars %>% rename(origin_wage = msa_mean_wage, origin_rent = msa_mean_rent) %>%
      select(met2013, origin_wage, origin_rent),
    by = c("origin_msa" = "met2013")
  ) %>%
  filter(!is.na(dest_wage), !is.na(origin_wage)) %>%
  mutate(
    changed_msa = origin_msa != met2013,
    wage_change = dest_wage - origin_wage,
    wage_pct_change = (dest_wage - origin_wage) / origin_wage * 100,
    moved_to_lower_wage = dest_wage < origin_wage,
    rent_change = dest_rent - origin_rent,
    rent_pct_change = (dest_rent - origin_rent) / origin_rent * 100
  )

cat("\nMovers with full origin+destination MSA data:", nrow(movers_full), "\n")
cat("Of which changed MSA:", sum(movers_full$changed_msa), "\n")

# =============================================================================
# STEP 6: ANALYSIS - Compare parents vs non-parents
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("ANALYSIS: ORIGIN VS DESTINATION MSA WAGES\n")
cat(strrep("=", 60), "\n")

# Focus on MSA changers
msa_changers <- movers_full %>%
  filter(changed_msa, new_parent %in% c("Non-Parent", "New Parent"))

cat("\nMSA changers for analysis:", nrow(msa_changers), "\n")

# Summary by parent status
summary_table <- msa_changers %>%
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
    mean_wage_change = weighted.mean(wage_change, perwt, na.rm = TRUE),
    mean_wage_pct = weighted.mean(wage_pct_change, perwt, na.rm = TRUE),
    pct_moved_lower_wage = weighted.mean(moved_to_lower_wage, perwt, na.rm = TRUE) * 100,
    mean_rent_change = weighted.mean(rent_change, perwt, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- ORIGIN VS DESTINATION: MSA Changers ---\n")
print(kable(summary_table, digits = 1, format = "simple"))

# Key statistics
cat("\n=== KEY FINDINGS ===\n")
np <- summary_table %>% filter(new_parent == "Non-Parent")
p <- summary_table %>% filter(new_parent == "New Parent")

cat(sprintf("\nNon-Parents who change MSA:\n"))
cat(sprintf("  Origin MSA wage: $%.0f → Dest MSA wage: $%.0f (change: $%.0f)\n",
            np$mean_origin_wage, np$mean_dest_wage, np$mean_wage_change))
cat(sprintf("  %.1f%% moved to lower-wage MSA\n", np$pct_moved_lower_wage))

cat(sprintf("\nNew Parents who change MSA:\n"))
cat(sprintf("  Origin MSA wage: $%.0f → Dest MSA wage: $%.0f (change: $%.0f)\n",
            p$mean_origin_wage, p$mean_dest_wage, p$mean_wage_change))
cat(sprintf("  %.1f%% moved to lower-wage MSA\n", p$pct_moved_lower_wage))

cat(sprintf("\nDifference in wage change: $%.0f\n", p$mean_wage_change - np$mean_wage_change))
cat(sprintf("Difference in pct moved to lower-wage MSA: %.1f pp\n",
            p$pct_moved_lower_wage - np$pct_moved_lower_wage))

# =============================================================================
# STEP 7: Regressions
# =============================================================================

cat("\n--- Regression: P(Moved to Lower-Wage MSA) ~ New Parent ---\n")
reg_lower <- feols(
  moved_to_lower_wage ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = msa_changers,
  weights = ~perwt
)
print(summary(reg_lower))

cat("\n--- Regression: MSA Wage % Change ~ New Parent ---\n")
reg_change <- feols(
  wage_pct_change ~ i(new_parent, ref = "Non-Parent") +
    age + I(age^2) + factor(sex) + factor(educd) | year,
  data = msa_changers %>% filter(abs(wage_pct_change) < 100),
  weights = ~perwt
)
print(summary(reg_change))

# =============================================================================
# STEP 8: Save results
# =============================================================================

write_csv(summary_table, file.path(OUTPUT_DIR, "msa_change_summary.csv"))

cat("\n=== Analysis complete ===\n")
cat("Output saved to:", OUTPUT_DIR, "\n")
