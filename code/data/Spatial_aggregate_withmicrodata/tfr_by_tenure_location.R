# ======================================================================================
#    TFR / Recent Birth Rate by Tenure x Location
#
#    Goal: Measure ACTUAL fertility gradient for owners vs renters
#    Using recent births (fertyr), not children in household (nchild)
# ======================================================================================

cat("=== TFR BY TENURE x LOCATION ===\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(fixest)
})

# Load data
data_path <- "processed_data/fertility_microdata_clean_10pct_stratified_sample_v7.rds"
if (!file.exists(data_path)) {
  data_path <- "processed_data/fertility_microdata_clean_5pct_stratified_sample.rds"
}
cat(sprintf("Loading: %s\n", data_path))
acs_data <- readRDS(data_path)
setDT(acs_data)
setnames(acs_data, old = names(acs_data), new = tolower(names(acs_data)))

# Create variables
acs_data[, met2013 := as.character(met2013)]
if (!"metro" %in% names(acs_data)) {
  acs_data[, metro := ifelse(met2013 == "0", 0, 1)]
}

# MSA population
msa_pop <- acs_data[met2013 != "0", .(msa_population = sum(perwt, na.rm = TRUE)), by = .(met2013, year)]
acs_data <- merge(acs_data, msa_pop, by = c("met2013", "year"), all.x = TRUE)
acs_data[met2013 == "0", msa_population := 0]

# Location type
acs_data[, location_type := fcase(
  metro == 0 | msa_population < 250000, "1_Peripheral",
  msa_population >= 250000 & msa_population < 2000000, "2_Secondary",
  msa_population >= 2000000, "3_Superstar",
  default = NA_character_
)]

# Tenure
if ("ownershp" %in% names(acs_data)) {
  acs_data[, is_owner := (ownershp == 1)]
  acs_data[, tenure := ifelse(is_owner, "Owner", "Renter")]
}

# Recent birth indicator (fertyr == 2 means birth in last 12 months)
if ("fertyr" %in% names(acs_data)) {
  acs_data[, had_birth := (fertyr == 2)]
} else {
  cat("WARNING: fertyr not in data, looking for is_recent_mother...\n")
  if ("is_recent_mother" %in% names(acs_data)) {
    acs_data[, had_birth := (is_recent_mother == 1)]
  }
}

# Check what we have
cat("\nVariables available:\n")
cat(sprintf("  fertyr: %s\n", "fertyr" %in% names(acs_data)))
cat(sprintf("  had_birth: %s\n", "had_birth" %in% names(acs_data)))
cat(sprintf("  tenure: %s\n", "tenure" %in% names(acs_data)))

# Filter to women of childbearing age (15-44)
women_fertile <- acs_data[sex == 2 & age >= 15 & age <= 44 & !is.na(location_type) & !is.na(tenure) & !is.na(had_birth)]

cat(sprintf("\nSample size: %d women aged 15-44\n", nrow(women_fertile)))
cat(sprintf("  With recent birth: %d (%.1f%%)\n",
            sum(women_fertile$had_birth),
            100 * mean(women_fertile$had_birth)))

# ======================================================================================
# KEY TABLE: Birth rate by tenure x location
# ======================================================================================
cat("\n===========================================\n")
cat("  BIRTH RATE (per 1000 women 15-44) BY TENURE x LOCATION\n")
cat("===========================================\n\n")

birth_rate_by_tenure_loc <- women_fertile[, .(
  birth_rate = weighted.mean(had_birth, perwt, na.rm = TRUE) * 1000,
  n_births = sum(had_birth * perwt, na.rm = TRUE),
  n_women = sum(perwt, na.rm = TRUE),
  n_obs = .N
), by = .(tenure, location_type)][order(tenure, location_type)]

# Print nicely
cat("Tenure     | Location    | Birth Rate | N births  | N women\n")
cat("-----------|-------------|------------|-----------|----------\n")
for (i in 1:nrow(birth_rate_by_tenure_loc)) {
  row <- birth_rate_by_tenure_loc[i]
  cat(sprintf("%-10s | %-11s |   %5.1f    | %9.0f | %9.0f\n",
              row$tenure, row$location_type, row$birth_rate, row$n_births, row$n_women))
}

# ======================================================================================
# Compute TFR approximation: sum of age-specific fertility rates * 5 / 1000
# ======================================================================================
cat("\n\n===========================================\n")
cat("  AGE-SPECIFIC FERTILITY RATES BY TENURE x LOCATION\n")
cat("===========================================\n\n")

# Create age groups
women_fertile[, age_group := cut(age, breaks = c(14, 19, 24, 29, 34, 39, 44),
                                  labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44"))]

asfr_by_tenure_loc <- women_fertile[, .(
  asfr = weighted.mean(had_birth, perwt, na.rm = TRUE) * 1000,
  n_obs = .N
), by = .(tenure, location_type, age_group)][order(tenure, location_type, age_group)]

# Compute TFR = sum of ASFR * 5 / 1000 (5-year age groups)
tfr_by_tenure_loc <- asfr_by_tenure_loc[, .(
  tfr = sum(asfr) * 5 / 1000
), by = .(tenure, location_type)]

cat("Tenure     | Location    | TFR\n")
cat("-----------|-------------|------\n")
for (i in 1:nrow(tfr_by_tenure_loc)) {
  row <- tfr_by_tenure_loc[i]
  cat(sprintf("%-10s | %-11s | %.2f\n",
              row$tenure, row$location_type, row$tfr))
}

# ======================================================================================
# GRADIENT COMPUTATION
# ======================================================================================
cat("\n\n===========================================\n")
cat("  TFR GRADIENT (Peripheral - Superstar)\n")
cat("===========================================\n\n")

gradient_by_tenure <- tfr_by_tenure_loc %>%
  pivot_wider(id_cols = tenure, names_from = location_type, values_from = tfr) %>%
  mutate(gradient = `1_Peripheral` - `3_Superstar`)

print(gradient_by_tenure)

cat("\n\nInterpretation:\n")
cat("  - Positive gradient = higher fertility in peripheral (cheaper) areas\n")
cat("  - Model predicts: Renters should show POSITIVE gradient\n")
cat("  - Model predicts: Owners may show ZERO or NEGATIVE gradient (selection)\n")

# ======================================================================================
# REGRESSION
# ======================================================================================
cat("\n\n===========================================\n")
cat("  REGRESSION: Birth Rate ~ Tenure x Location\n")
cat("===========================================\n\n")

women_fertile[, loc_superstar := (location_type == "3_Superstar")]
women_fertile[, loc_peripheral := (location_type == "1_Peripheral")]

# Model: birth ~ tenure * location
m1 <- feglm(had_birth ~ is_owner * loc_superstar + factor(age),
            data = women_fertile, family = "binomial", weights = ~perwt)

cat("Logit: had_birth ~ is_owner * loc_superstar + age FE\n\n")

# Get marginal effects
cat("Coefficients (log-odds):\n")
print(coef(summary(m1))[1:4, ])

# Alternative: compare peripheral vs superstar
m2 <- feglm(had_birth ~ is_owner * loc_peripheral + factor(age),
            data = women_fertile, family = "binomial", weights = ~perwt)

cat("\n\nLogit: had_birth ~ is_owner * loc_peripheral + age FE\n\n")
cat("Coefficients (log-odds):\n")
print(coef(summary(m2))[1:4, ])

cat("\n\nKey interpretation:\n")
cat("  loc_peripheral coef > 0 means: Renters have HIGHER fertility in peripheral\n")
cat("  is_owner:loc_peripheral coef tells us if owners differ\n")

# ======================================================================================
# SAVE RESULTS
# ======================================================================================
write.csv(tfr_by_tenure_loc, "calibration_targets_output/tfr_by_tenure_location.csv", row.names = FALSE)
write.csv(birth_rate_by_tenure_loc, "calibration_targets_output/birth_rate_by_tenure_location.csv", row.names = FALSE)

cat("\n\nResults saved to calibration_targets_output/\n")
cat("Done!\n")
