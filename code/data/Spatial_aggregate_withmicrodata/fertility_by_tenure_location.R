# ======================================================================================
#    Quick Analysis: Fertility by Tenure x Location
#
#    Goal: Test the model prediction that renters should show fertility gradient
#          but owners should not (or reverse)
# ======================================================================================

cat("=== FERTILITY BY TENURE x LOCATION ===\n\n")

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

# Filter to women 45+ (completed fertility)
women_completed <- acs_data[sex == 2 & age >= 45 & age <= 65 & !is.na(nchild) & !is.na(location_type) & !is.na(tenure)]

cat(sprintf("\nSample size: %d women aged 45-65\n", nrow(women_completed)))

# ======================================================================================
# KEY TABLE: Mean fertility by tenure x location
# ======================================================================================
cat("\n===========================================\n")
cat("  MEAN COMPLETED FERTILITY BY TENURE x LOCATION\n")
cat("===========================================\n\n")

fert_by_tenure_loc <- women_completed[, .(
  mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
  pct_childless = weighted.mean(nchild == 0, perwt, na.rm = TRUE) * 100,
  n_obs = .N,
  pop_weight = sum(perwt, na.rm = TRUE)
), by = .(tenure, location_type)][order(tenure, location_type)]

# Print nicely
cat("Tenure     | Location    | Mean N | Childless | N obs\n")
cat("-----------|-------------|--------|-----------|-------\n")
for (i in 1:nrow(fert_by_tenure_loc)) {
  row <- fert_by_tenure_loc[i]
  cat(sprintf("%-10s | %-11s | %5.2f  |   %5.1f%%  | %6d\n",
              row$tenure, row$location_type, row$mean_nchild, row$pct_childless, row$n_obs))
}

# Compute gradients
cat("\n\n===========================================\n")
cat("  FERTILITY GRADIENT (Peripheral - Superstar)\n")
cat("===========================================\n\n")

gradient_by_tenure <- fert_by_tenure_loc %>%
  pivot_wider(id_cols = tenure, names_from = location_type, values_from = mean_nchild) %>%
  mutate(gradient = `1_Peripheral` - `3_Superstar`)

print(gradient_by_tenure)

# Save results
output <- list(
  by_tenure_loc = fert_by_tenure_loc,
  gradients = gradient_by_tenure
)
write.csv(fert_by_tenure_loc, "calibration_targets_output/fertility_by_tenure_location.csv", row.names = FALSE)

cat("\n\n===========================================\n")
cat("  REGRESSION: Tenure x Location Interaction\n")
cat("===========================================\n\n")

# Regression with tenure x location interaction
women_completed[, loc_superstar := (location_type == "3_Superstar")]

# Simple regression
m1 <- feols(nchild ~ is_owner * loc_superstar | age, data = women_completed, weights = ~perwt)
cat("Model: nchild ~ is_owner * loc_superstar | age FE\n\n")
print(summary(m1))

cat("\n\nInterpretation:\n")
cat("  - Coefficient on loc_superstar = effect for RENTERS\n")
cat("  - Coefficient on is_owner:loc_superstar = DIFFERENCE in effect for owners\n")
cat("  - If model prediction holds:\n")
cat("      * loc_superstar should be NEGATIVE (renters have fewer kids in superstars)\n")
cat("      * is_owner:loc_superstar should be POSITIVE (owners do not show this pattern)\n")

cat("\n\nDone!\n")
