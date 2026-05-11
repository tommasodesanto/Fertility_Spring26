# ======================================================================================
#    CALIBRATION TARGETS: RTI-BASED LOCATION CLASSIFICATION
#
#    This script computes calibration targets using rent-to-income (RTI) terciles
#    to classify locations, rather than arbitrary population thresholds.
#
#    Rationale: The model's mechanism is about housing costs affecting fertility.
#    Classifying by housing affordability (RTI) directly captures this mechanism.
#
#    Author: Tommaso De Santo
#    Date: January 2026
# ======================================================================================

suppressPackageStartupMessages({
  library(data.table)
})

cat("======================================================================\n")
cat("  CALIBRATION TARGETS: RTI-BASED LOCATION CLASSIFICATION\n")
cat("======================================================================\n\n")

# Load data
data_path <- "processed_data/fertility_microdata_clean_10pct_stratified_sample_v7.rds"
if (!file.exists(data_path)) {
  data_path <- "processed_data/fertility_microdata_clean_5pct_stratified_sample.rds"
}
cat("Loading:", data_path, "\n\n")
acs_data <- readRDS(data_path)
setDT(acs_data)
setnames(acs_data, old = names(acs_data), new = tolower(names(acs_data)))

acs_data[, met2013 := as.character(met2013)]

# =============================================================================
# STEP 1: Compute MSA-level median INCOME and classify into terciles
# =============================================================================
cat("--- 1. LOCATION CLASSIFICATION (WAGE/INCOME Terciles) ---\n\n")

# Clean income variable
acs_data[hhincome <= 0 | hhincome >= 9999999, hhincome_clean := NA_real_]
acs_data[hhincome > 0 & hhincome < 9999999, hhincome_clean := hhincome]

# MSA-level median income
msa_stats <- acs_data[met2013 != "0" & !is.na(hhincome_clean), .(
  msa_income = median(hhincome_clean, na.rm = TRUE),
  msa_pop = sum(perwt, na.rm = TRUE)
), by = met2013]

# Custom split: 40% Peripheral, 40% Secondary, 20% Superstar (by MSA count)
# Top 20% of MSAs by wage = Superstar (smaller, truly high-wage)
msa_stats <- msa_stats[order(msa_income)]
n_msa <- nrow(msa_stats)
msa_stats[, msa_rank := 1:.N]

# Bottom 40%, middle 40% (40-80), top 20% (80-100)
q40_rank <- floor(n_msa * 0.40)
q80_rank <- floor(n_msa * 0.80)

q33 <- msa_stats[msa_rank == q40_rank, msa_income]  # 40th percentile
q67 <- msa_stats[msa_rank == q80_rank, msa_income]  # 80th percentile

cat("INCOME tercile cutoffs (population-weighted):\n")
cat(sprintf("  Peripheral (low wage):   median_income < $%.0f\n", q33))
cat(sprintf("  Secondary (medium):      median_income $%.0f - $%.0f\n", q33, q67))
cat(sprintf("  Superstar (high wage):   median_income >= $%.0f\n\n", q67))

# Classify MSAs by INCOME (Superstar = HIGH wage)
msa_stats[, location_type := fcase(
  msa_income < q33, "1_Peripheral",
  msa_income >= q33 & msa_income < q67, "2_Secondary",
  msa_income >= q67, "3_Superstar"
)]

# Merge to individual data
acs_data <- merge(acs_data, msa_stats[, .(met2013, msa_income, location_type)],
                  by = "met2013", all.x = TRUE)
# NOTE: Non-metro (met2013 == "0") is EXCLUDED - MSA-only analysis
# acs_data[met2013 == "0", location_type := "1_Peripheral"]  # Commented out

# Prepare variables
acs_data[, had_birth := (fertyr == 2)]
acs_data[, is_owner := (ownershp == 1)]
acs_data[, is_female := (sex == 2)]
acs_data[, has_children := (nchild > 0)]

# =============================================================================
# STEP 2: Population Shares
# =============================================================================
cat("--- 2. POPULATION SHARES ---\n\n")
pop_shares <- acs_data[!is.na(location_type), .(pop = sum(perwt)), by = location_type]
pop_shares <- pop_shares[order(location_type)]
pop_shares[, pct := round(100 * pop / sum(pop), 1)]
print(pop_shares[, .(location_type, pct)])

# =============================================================================
# STEP 3: TFR by Location
# =============================================================================
cat("\n--- 3. TOTAL FERTILITY RATE (TFR) BY LOCATION ---\n\n")
women_fertile <- acs_data[is_female == TRUE & age >= 15 & age <= 44 &
                          !is.na(location_type) & !is.na(had_birth)]
women_fertile[, age_group := cut(age, breaks = c(14, 19, 24, 29, 34, 39, 44),
                                  labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44"))]

asfr <- women_fertile[, .(asfr = weighted.mean(had_birth, perwt, na.rm = TRUE) * 1000),
                      by = .(location_type, age_group)]
tfr <- asfr[, .(tfr = round(sum(asfr) * 5 / 1000, 3)), by = location_type]
tfr <- tfr[order(location_type)]
print(tfr)

tfr_gradient <- tfr[location_type == "1_Peripheral", tfr] -
                tfr[location_type == "3_Superstar", tfr]
cat(sprintf("\nTFR GRADIENT (Peripheral - Superstar): +%.3f\n", tfr_gradient))

# =============================================================================
# STEP 4: TFR by Tenure x Location
# =============================================================================
cat("\n--- 4. TFR BY TENURE x LOCATION ---\n\n")
women_fertile[, tenure := ifelse(ownershp == 1, "Owner", "Renter")]
asfr_tenure <- women_fertile[!is.na(tenure),
                             .(asfr = weighted.mean(had_birth, perwt, na.rm = TRUE) * 1000),
                             by = .(tenure, location_type, age_group)]
tfr_tenure <- asfr_tenure[, .(tfr = round(sum(asfr) * 5 / 1000, 3)),
                          by = .(tenure, location_type)]
tfr_tenure <- tfr_tenure[order(tenure, location_type)]
print(tfr_tenure)

tfr_grad_owner <- tfr_tenure[tenure == "Owner" & location_type == "1_Peripheral", tfr] -
                  tfr_tenure[tenure == "Owner" & location_type == "3_Superstar", tfr]
tfr_grad_renter <- tfr_tenure[tenure == "Renter" & location_type == "1_Peripheral", tfr] -
                   tfr_tenure[tenure == "Renter" & location_type == "3_Superstar", tfr]
cat(sprintf("\nTFR Gradient - Owners:  +%.3f\n", tfr_grad_owner))
cat(sprintf("TFR Gradient - Renters: +%.3f\n", tfr_grad_renter))

# =============================================================================
# STEP 5: NCHILD (women 40-50) - Stock measure
# =============================================================================
cat("\n--- 5. MEAN NCHILD (Women 40-50) - Stock Measure ---\n\n")
women_40_50 <- acs_data[is_female == TRUE & age >= 40 & age <= 50 & !is.na(location_type)]
nchild_stats <- women_40_50[, .(
  mean_nchild = round(weighted.mean(nchild, perwt, na.rm = TRUE), 3),
  pct_childless = round(100 * weighted.mean(nchild == 0, perwt, na.rm = TRUE), 1)
), by = location_type]
nchild_stats <- nchild_stats[order(location_type)]
print(nchild_stats)

nchild_gradient <- nchild_stats[location_type == "1_Peripheral", mean_nchild] -
                   nchild_stats[location_type == "3_Superstar", mean_nchild]
cat(sprintf("\nNCHILD Gradient (P - S): %.3f\n", nchild_gradient))
cat("NOTE: OPPOSITE sign to TFR! This is selection + children leaving home.\n")

# =============================================================================
# STEP 6: Homeownership
# =============================================================================
cat("\n--- 6. HOMEOWNERSHIP RATE ---\n\n")
own_by_loc <- acs_data[!is.na(location_type) & !is.na(is_owner), .(
  ownership_rate = round(100 * weighted.mean(is_owner, perwt, na.rm = TRUE), 1)
), by = location_type]
own_by_loc <- own_by_loc[order(location_type)]
print(own_by_loc)

own_gradient <- own_by_loc[location_type == "1_Peripheral", ownership_rate] -
                own_by_loc[location_type == "3_Superstar", ownership_rate]
cat(sprintf("\nOwnership Gradient (P - S): +%.1f pp\n", own_gradient))

# =============================================================================
# STEP 7: Ownership by Family Status
# =============================================================================
cat("\n--- 7. OWNERSHIP BY FAMILY STATUS x LOCATION ---\n\n")
own_family <- acs_data[!is.na(location_type) & !is.na(is_owner) & age >= 25 & age <= 55, .(
  ownership_rate = round(100 * weighted.mean(is_owner, perwt, na.rm = TRUE), 1)
), by = .(location_type, has_children)]
own_family <- own_family[order(location_type, has_children)]
print(dcast(own_family, location_type ~ has_children, value.var = "ownership_rate"))

# Family-childless gap
family_gap <- acs_data[!is.na(location_type) & !is.na(is_owner) & age >= 25 & age <= 55, .(
  gap = round(100 * (weighted.mean(is_owner[has_children == TRUE], perwt[has_children == TRUE], na.rm = TRUE) -
               weighted.mean(is_owner[has_children == FALSE], perwt[has_children == FALSE], na.rm = TRUE)), 1)
), by = location_type]
family_gap <- family_gap[order(location_type)]
cat("\nFamily-Childless Ownership Gap by Location:\n")
print(family_gap)

# =============================================================================
# STEP 8: Relative Income (wage proxy)
# =============================================================================
cat("\n--- 8. MEDIAN HOUSEHOLD INCOME BY LOCATION ---\n\n")
acs_data[hhincome <= 0 | hhincome >= 9999999, hhincome := NA]
income_by_loc <- acs_data[!is.na(location_type) & !is.na(hhincome), .(
  median_income = median(hhincome, na.rm = TRUE)
), by = location_type]
income_by_loc <- income_by_loc[order(location_type)]
income_by_loc[, relative := round(median_income /
                                  income_by_loc[location_type == "2_Secondary", median_income], 3)]
print(income_by_loc)

# =============================================================================
# STEP 9: Median Rent by Location (for rent ratios)
# =============================================================================
cat("\n--- 9. MEDIAN RENT BY LOCATION ---\n\n")
# Use monthly rent for renters only
rent_by_loc <- acs_data[!is.na(location_type) & rent > 0 & ownershp != 1, .(
  median_rent = median(rent, na.rm = TRUE)
), by = location_type]
rent_by_loc <- rent_by_loc[order(location_type)]
rent_by_loc[, relative := round(median_rent /
                                rent_by_loc[location_type == "2_Secondary", median_rent], 3)]
print(rent_by_loc)

cat("\nRENT RATIOS (vs Secondary = 1.00):\n")
cat(sprintf("  rent_ratio_P = %.3f  (Peripheral / Secondary)\n",
            rent_by_loc[location_type == "1_Peripheral", relative]))
cat(sprintf("  rent_ratio_X = %.3f  (Superstar / Secondary)\n",
            rent_by_loc[location_type == "3_Superstar", relative]))

# =============================================================================
# SUMMARY TABLE
# =============================================================================
cat("\n\n======================================================================\n")
cat("                    SUMMARY: CALIBRATION TARGETS                      \n")
cat("              (WAGE/INCOME-Based Location Classification)             \n")
cat("======================================================================\n\n")

cat("POPULATION SHARES:\n")
cat(sprintf("  Peripheral:  %.0f%%\n", pop_shares[location_type == "1_Peripheral", pct]))
cat(sprintf("  Secondary:   %.0f%%\n", pop_shares[location_type == "2_Secondary", pct]))
cat(sprintf("  Superstar:   %.0f%%\n\n", pop_shares[location_type == "3_Superstar", pct]))

cat("TFR (Total Fertility Rate):\n")
cat(sprintf("  Peripheral:  %.2f\n", tfr[location_type == "1_Peripheral", tfr]))
cat(sprintf("  Secondary:   %.2f\n", tfr[location_type == "2_Secondary", tfr]))
cat(sprintf("  Superstar:   %.2f\n", tfr[location_type == "3_Superstar", tfr]))
cat(sprintf("  GRADIENT:    +%.2f\n\n", tfr_gradient))

cat("TFR GRADIENT BY TENURE:\n")
cat(sprintf("  Owners:      +%.2f\n", tfr_grad_owner))
cat(sprintf("  Renters:     +%.2f\n\n", tfr_grad_renter))

cat("HOMEOWNERSHIP:\n")
cat(sprintf("  Peripheral:  %.1f%%\n", own_by_loc[location_type == "1_Peripheral", ownership_rate]))
cat(sprintf("  Secondary:   %.1f%%\n", own_by_loc[location_type == "2_Secondary", ownership_rate]))
cat(sprintf("  Superstar:   %.1f%%\n", own_by_loc[location_type == "3_Superstar", ownership_rate]))
cat(sprintf("  GRADIENT:    +%.1f pp\n\n", own_gradient))

cat("FAMILY-CHILDLESS OWNERSHIP GAP:\n")
cat(sprintf("  Average:     +%.0f pp\n\n", mean(family_gap$gap)))

cat("RELATIVE INCOME (wage ratios, vs Secondary = 1.00):\n")
cat(sprintf("  wage_ratio_P:  %.3f  (Peripheral / Secondary)\n", income_by_loc[location_type == "1_Peripheral", relative]))
cat(sprintf("  wage_ratio_X:  %.3f  (Superstar / Secondary)\n\n", income_by_loc[location_type == "3_Superstar", relative]))

cat("RELATIVE RENT (rent ratios, vs Secondary = 1.00):\n")
cat(sprintf("  rent_ratio_P:  %.3f  (Peripheral / Secondary)\n", rent_by_loc[location_type == "1_Peripheral", relative]))
cat(sprintf("  rent_ratio_X:  %.3f  (Superstar / Secondary)\n\n", rent_by_loc[location_type == "3_Superstar", relative]))

cat("CHILDLESSNESS (external): 15%% (Census Bureau vital statistics)\n\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
output_dir <- "calibration_targets_output"
dir.create(output_dir, showWarnings = FALSE)

write.csv(tfr, file.path(output_dir, "tfr_by_location_rti.csv"), row.names = FALSE)
write.csv(tfr_tenure, file.path(output_dir, "tfr_by_tenure_location_rti.csv"), row.names = FALSE)
write.csv(own_by_loc, file.path(output_dir, "ownership_by_location_rti.csv"), row.names = FALSE)
write.csv(pop_shares, file.path(output_dir, "population_shares_rti.csv"), row.names = FALSE)
write.csv(income_by_loc, file.path(output_dir, "income_by_location_rti.csv"), row.names = FALSE)
write.csv(rent_by_loc, file.path(output_dir, "rent_by_location_rti.csv"), row.names = FALSE)
write.csv(nchild_stats, file.path(output_dir, "nchild_by_location_rti.csv"), row.names = FALSE)
write.csv(family_gap, file.path(output_dir, "family_ownership_gap_rti.csv"), row.names = FALSE)

cat("Outputs saved to:", output_dir, "\n")
cat("\nDone!\n")
