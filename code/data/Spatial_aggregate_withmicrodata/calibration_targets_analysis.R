# ======================================================================================
#
#    Calibration Targets Analysis for Spatial Fertility Model
#
#    Purpose: Generate empirically-grounded calibration targets for the structural model
#             with three location types (Peripheral / Secondary / Superstar)
#
#    Key Outputs:
#    1. Completed fertility by location type (women 40-50)
#    2. Childlessness rate by location type (women 45-50)
#    3. Fertility gradient (unconditional and conditional on education/income)
#    4. Homeownership rate by location type
#
#    Author: Tommaso De Santo
#    Date: January 2026
#
# ======================================================================================

# ======================================================================================
# 0. SETUP
# ======================================================================================
cat("=== CALIBRATION TARGETS ANALYSIS ===\n")
cat("Setting up environment...\n\n")

# Load packages
required_packages <- c("tidyverse", "data.table", "fixest", "modelsummary",
                       "matrixStats", "kableExtra", "haven")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Output directory
OUTPUT_DIR <- "calibration_targets_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================================================================================
# 1. DATA LOADING
# ======================================================================================
cat("--- 1. LOADING DATA ---\n")

# Try different possible data paths
data_paths <- c(
  "processed_data/fertility_microdata_clean_10pct_stratified_sample_v7.rds",
  "processed_data/fertility_microdata_clean_5pct_stratified_sample.rds",
  "output/validated_microdata_sample.rds"
)

data_loaded <- FALSE
for (path in data_paths) {
  if (file.exists(path)) {
    cat(sprintf("Loading data from: %s\n", path))
    acs_data <- readRDS(path)
    data_loaded <- TRUE
    break
  }
}

if (!data_loaded) {
  stop("No data file found. Please check data paths.")
}

setDT(acs_data)
setnames(acs_data, old = names(acs_data), new = tolower(names(acs_data)))
cat(sprintf("Data loaded: %s observations, %s variables\n\n",
            format(nrow(acs_data), big.mark=","), ncol(acs_data)))

# ======================================================================================
# 2. VARIABLE PREPARATION
# ======================================================================================
cat("--- 2. PREPARING VARIABLES ---\n")

# Ensure key variables exist and are properly coded
acs_data[, met2013 := as.character(met2013)]

# Create metro status if not exists
if (!"metro" %in% names(acs_data)) {
  # metro = 0 means not in metro area, met2013 = "0" also indicates non-metro
  acs_data[, metro := ifelse(met2013 == "0", 0, 1)]
}

# Calculate MSA population by year
cat("  Calculating MSA populations...\n")
msa_pop <- acs_data[met2013 != "0",
                    .(msa_population = sum(perwt, na.rm = TRUE)),
                    by = .(met2013, year)]

# Merge back to main data
acs_data <- merge(acs_data, msa_pop, by = c("met2013", "year"), all.x = TRUE)
acs_data[met2013 == "0", msa_population := 0]

# ======================================================================================
# 3. DEFINE LOCATION TYPES
# ======================================================================================
cat("  Defining location types (Peripheral / Secondary / Superstar)...\n")

# Location type based on metro status and population
acs_data[, location_type := fcase(
  metro == 0 | msa_population < 250000, "1_Peripheral",
  msa_population >= 250000 & msa_population < 2000000, "2_Secondary",
  msa_population >= 2000000, "3_Superstar",
  default = NA_character_
)]

# Check distribution
cat("\n  Location type distribution:\n")
loc_dist <- acs_data[!is.na(location_type),
                     .(n = .N,
                       pop = sum(perwt, na.rm=TRUE)),
                     by = location_type][order(location_type)]
loc_dist[, pct := pop / sum(pop) * 100]
print(loc_dist)

# ======================================================================================
# 4. PREPARE ANALYSIS SAMPLES
# ======================================================================================
cat("\n--- 3. PREPARING ANALYSIS SAMPLES ---\n")

# Ensure sex variable is coded correctly (2 = female in IPUMS)
if (max(acs_data$sex, na.rm=TRUE) == 2) {
  acs_data[, is_female := (sex == 2)]
} else {
  acs_data[, is_female := (sex == 1)]  # Adjust if coding is different
}

# Education factor (create if not exists)
if (!"educ_factor" %in% names(acs_data)) {
  if ("educd" %in% names(acs_data)) {
    acs_data[, educ_factor := fcase(
      educd < 62, "1_Less than HS",
      educd >= 62 & educd < 65, "2_HS Grad",
      educd >= 65 & educd < 101, "3_Some College",
      educd >= 101 & educd < 114, "4_Bachelors",
      educd >= 114, "5_Graduate",
      default = NA_character_
    )]
  } else if ("educ" %in% names(acs_data)) {
    acs_data[, educ_factor := fcase(
      educ < 6, "1_Less than HS",
      educ == 6, "2_HS Grad",
      educ >= 7 & educ < 10, "3_Some College",
      educ == 10, "4_Bachelors",
      educ >= 11, "5_Graduate",
      default = NA_character_
    )]
  }
}

# Clean income
acs_data[hhincome <= 0 | hhincome == 9999999, hhincome := NA]
acs_data[!is.na(hhincome) & hhincome > 0, log_hhincome := log(hhincome)]

# Ownership (ownershp = 1 means owned in IPUMS)
if ("ownershp" %in% names(acs_data)) {
  acs_data[, is_owner := (ownershp == 1)]
} else if ("ownershpd" %in% names(acs_data)) {
  acs_data[, is_owner := (ownershpd %in% c(12, 13))]  # Owned free and clear, or with mortgage
}

# Create analysis samples
women_all <- acs_data[is_female == TRUE & !is.na(location_type)]
women_40_50 <- acs_data[is_female == TRUE & age %between% c(40, 50) & !is.na(location_type)]
women_45_50 <- acs_data[is_female == TRUE & age %between% c(45, 50) & !is.na(location_type)]
women_15_50 <- acs_data[is_female == TRUE & age %between% c(15, 50) & !is.na(location_type)]

cat(sprintf("  Women 40-50 sample: %s observations\n", format(nrow(women_40_50), big.mark=",")))
cat(sprintf("  Women 45-50 sample: %s observations\n", format(nrow(women_45_50), big.mark=",")))

# ======================================================================================
# 5. COMPLETED FERTILITY BY LOCATION TYPE
# ======================================================================================
cat("\n--- 4. COMPLETED FERTILITY BY LOCATION TYPE ---\n")
cat("  (Using women aged 40-50 as proxy for completed fertility)\n\n")

# Summary statistics
fertility_by_location <- women_40_50[, .(
  mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
  median_nchild = tryCatch(weightedMedian(nchild, perwt, na.rm = TRUE), error = function(e) NA_real_),
  sd_nchild = sqrt(weighted.mean((nchild - weighted.mean(nchild, perwt, na.rm=TRUE))^2, perwt, na.rm = TRUE)),
  pct_0_children = weighted.mean(nchild == 0, perwt, na.rm = TRUE) * 100,
  pct_1_child = weighted.mean(nchild == 1, perwt, na.rm = TRUE) * 100,
  pct_2_children = weighted.mean(nchild == 2, perwt, na.rm = TRUE) * 100,
  pct_3plus = weighted.mean(nchild >= 3, perwt, na.rm = TRUE) * 100,
  n_obs = .N,
  pop_weight = sum(perwt, na.rm = TRUE)
), by = location_type][order(location_type)]

cat("=== MEAN CHILDREN IN HOUSEHOLD (Women 40-50) ===\n")
print(fertility_by_location[, .(location_type, mean_nchild, sd_nchild, n_obs)])

# Calculate gradients
peripheral_fertility <- fertility_by_location[location_type == "1_Peripheral", mean_nchild]
secondary_fertility <- fertility_by_location[location_type == "2_Secondary", mean_nchild]
superstar_fertility <- fertility_by_location[location_type == "3_Superstar", mean_nchild]

cat(sprintf("\n  Fertility gradient (Peripheral - Superstar): %.3f children\n",
            peripheral_fertility - superstar_fertility))
cat(sprintf("  Fertility gradient (Peripheral - Secondary): %.3f children\n",
            peripheral_fertility - secondary_fertility))

# ======================================================================================
# 6. CHILDLESSNESS RATE BY LOCATION TYPE
# ======================================================================================
cat("\n--- 5. CHILDLESSNESS RATE BY LOCATION TYPE ---\n")
cat("  (Using women aged 45-50, NCHILD == 0 as proxy)\n")
cat("  NOTE: This is a LOWER BOUND - children may have left home\n\n")

childlessness_by_location <- women_45_50[, .(
  pct_childless = weighted.mean(nchild == 0, perwt, na.rm = TRUE) * 100,
  n_obs = .N,
  pop_weight = sum(perwt, na.rm = TRUE)
), by = location_type][order(location_type)]

cat("=== CHILDLESSNESS RATE (Women 45-50, NCHILD == 0) ===\n")
print(childlessness_by_location)

# Calculate gradient
peripheral_childless <- childlessness_by_location[location_type == "1_Peripheral", pct_childless]
superstar_childless <- childlessness_by_location[location_type == "3_Superstar", pct_childless]

cat(sprintf("\n  Childlessness gradient (Superstar - Peripheral): %.1f percentage points\n",
            superstar_childless - peripheral_childless))

# ======================================================================================
# 7. REGRESSION ANALYSIS: FERTILITY GRADIENT
# ======================================================================================
cat("\n--- 6. REGRESSION ANALYSIS: FERTILITY GRADIENT ---\n")

# Set reference category
women_40_50[, location_type := relevel(factor(location_type), ref = "2_Secondary")]

# Model 1: Unconditional (only age FE)
cat("\n  Model 1: Unconditional (age FE only)\n")
m1_unconditional <- feols(nchild ~ location_type | age,
                          data = women_40_50,
                          weights = ~perwt,
                          vcov = "hetero")

# Model 2: Conditional on education
cat("  Model 2: Conditional on education\n")
m2_educ <- feols(nchild ~ location_type + educ_factor | age,
                 data = women_40_50[!is.na(educ_factor)],
                 weights = ~perwt,
                 vcov = "hetero")

# Model 3: Conditional on education + income
cat("  Model 3: Conditional on education + income\n")
m3_full <- feols(nchild ~ location_type + educ_factor + log_hhincome | age,
                 data = women_40_50[!is.na(educ_factor) & !is.na(log_hhincome)],
                 weights = ~perwt,
                 vcov = "hetero")

# Model 4: With year FE
cat("  Model 4: With year fixed effects\n")
m4_year_fe <- feols(nchild ~ location_type + educ_factor + log_hhincome | age + year,
                    data = women_40_50[!is.na(educ_factor) & !is.na(log_hhincome)],
                    weights = ~perwt,
                    vcov = "hetero")

# Display results
cat("\n=== REGRESSION RESULTS: NUMBER OF CHILDREN (Women 40-50) ===\n")
cat("Reference category: Secondary (medium metros 250K-2M)\n\n")

etable(m1_unconditional, m2_educ, m3_full, m4_year_fe,
       headers = c("Unconditional", "+ Education", "+ Income", "+ Year FE"),
       signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
       fitstat = ~ n + r2)

# Save regression table
etable(m1_unconditional, m2_educ, m3_full, m4_year_fe,
       headers = c("Unconditional", "+ Education", "+ Income", "+ Year FE"),
       signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
       fitstat = ~ n + r2,
       file = file.path(OUTPUT_DIR, "fertility_gradient_regressions.tex"),
       tex = TRUE)

# ======================================================================================
# 8. CHILDLESSNESS REGRESSION (LOGIT)
# ======================================================================================
cat("\n--- 7. CHILDLESSNESS REGRESSION (LOGIT) ---\n")

women_45_50[, location_type := relevel(factor(location_type), ref = "2_Secondary")]
women_45_50[, is_childless := as.integer(nchild == 0)]

# Logit models
m_childless_1 <- feglm(is_childless ~ location_type | age,
                       data = women_45_50,
                       family = "binomial",
                       weights = ~perwt)

m_childless_2 <- feglm(is_childless ~ location_type + educ_factor | age,
                       data = women_45_50[!is.na(educ_factor)],
                       family = "binomial",
                       weights = ~perwt)

m_childless_3 <- feglm(is_childless ~ location_type + educ_factor + log_hhincome | age,
                       data = women_45_50[!is.na(educ_factor) & !is.na(log_hhincome)],
                       family = "binomial",
                       weights = ~perwt)

cat("\n=== LOGIT REGRESSION: CHILDLESSNESS (Women 45-50) ===\n")
cat("Reference category: Secondary. Coefficients are log-odds.\n\n")

etable(m_childless_1, m_childless_2, m_childless_3,
       headers = c("Unconditional", "+ Education", "+ Income"),
       signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
       fitstat = ~ n + pr2)

# ======================================================================================
# 9. HOMEOWNERSHIP BY LOCATION TYPE
# ======================================================================================
cat("\n--- 8. HOMEOWNERSHIP RATE BY LOCATION TYPE ---\n")

if ("is_owner" %in% names(acs_data)) {
  # By location type (all households)
  ownership_by_location <- acs_data[!is.na(location_type) & !is.na(is_owner), .(
    ownership_rate = weighted.mean(is_owner, perwt, na.rm = TRUE) * 100,
    n_obs = .N
  ), by = location_type][order(location_type)]

  cat("\n=== HOMEOWNERSHIP RATE BY LOCATION TYPE (All) ===\n")
  print(ownership_by_location)

  # By location type and age group
  acs_data[, age_group := fcase(
    age %between% c(25, 34), "25-34",
    age %between% c(35, 44), "35-44",
    age %between% c(45, 54), "45-54",
    age %between% c(55, 64), "55-64",
    default = NA_character_
  )]

  ownership_by_loc_age <- acs_data[!is.na(location_type) & !is.na(is_owner) & !is.na(age_group), .(
    ownership_rate = weighted.mean(is_owner, perwt, na.rm = TRUE) * 100,
    n_obs = .N
  ), by = .(location_type, age_group)][order(location_type, age_group)]

  cat("\n=== HOMEOWNERSHIP RATE BY LOCATION TYPE AND AGE ===\n")
  print(dcast(ownership_by_loc_age, age_group ~ location_type, value.var = "ownership_rate"))
}

# ======================================================================================
# 10. FERTILITY BY AGE AND LOCATION (FOR TFR-STYLE ANALYSIS)
# ======================================================================================
cat("\n--- 9. FERTILITY RATE BY AGE AND LOCATION ---\n")

# Recent birth rate by age group and location
if ("fertyr" %in% names(acs_data) | "is_recent_mother" %in% names(acs_data)) {

  if (!"is_recent_mother" %in% names(acs_data)) {
    acs_data[, is_recent_mother := as.numeric(fertyr == 2 & sex == 2)]
  }

  women_15_50[, age_group_5yr := cut(age,
                                      breaks = c(14, 19, 24, 29, 34, 39, 44, 50),
                                      labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-50"))]

  asfr_by_location <- women_15_50[!is.na(age_group_5yr), .(
    birth_rate = weighted.mean(is_recent_mother, perwt, na.rm = TRUE) * 1000,
    n_obs = .N
  ), by = .(location_type, age_group_5yr)][order(location_type, age_group_5yr)]

  cat("\n=== AGE-SPECIFIC FERTILITY RATE (births per 1000 women) ===\n")
  print(dcast(asfr_by_location, age_group_5yr ~ location_type, value.var = "birth_rate"))

  # Calculate TFR by location
  tfr_by_location <- asfr_by_location[, .(
    tfr = 5 * sum(birth_rate / 1000, na.rm = TRUE)
  ), by = location_type]

  cat("\n=== TOTAL FERTILITY RATE BY LOCATION ===\n")
  print(tfr_by_location)
}

# ======================================================================================
# 11. SUMMARY TABLE FOR CALIBRATION
# ======================================================================================
cat("\n\n")
cat("================================================================================\n")
cat("                    SUMMARY: CALIBRATION TARGETS                               \n")
cat("================================================================================\n\n")

# Create summary table
summary_table <- data.table(
  Moment = c(
    "Mean children (women 40-50)",
    "Mean children (women 40-50)",
    "Mean children (women 40-50)",
    "Fertility gradient (Periph - Super)",
    "",
    "Childlessness rate (women 45-50)",
    "Childlessness rate (women 45-50)",
    "Childlessness rate (women 45-50)",
    "Childlessness gradient (Super - Periph)",
    "",
    "Homeownership rate",
    "Homeownership rate",
    "Homeownership rate"
  ),
  Location = c(
    "Peripheral", "Secondary", "Superstar", "Gradient", "",
    "Peripheral", "Secondary", "Superstar", "Gradient", "",
    "Peripheral", "Secondary", "Superstar"
  ),
  Value = c(
    sprintf("%.3f", fertility_by_location[location_type == "1_Peripheral", mean_nchild]),
    sprintf("%.3f", fertility_by_location[location_type == "2_Secondary", mean_nchild]),
    sprintf("%.3f", fertility_by_location[location_type == "3_Superstar", mean_nchild]),
    sprintf("%.3f", peripheral_fertility - superstar_fertility),
    "",
    sprintf("%.1f%%", childlessness_by_location[location_type == "1_Peripheral", pct_childless]),
    sprintf("%.1f%%", childlessness_by_location[location_type == "2_Secondary", pct_childless]),
    sprintf("%.1f%%", childlessness_by_location[location_type == "3_Superstar", pct_childless]),
    sprintf("%.1f pp", superstar_childless - peripheral_childless),
    "",
    ifelse("is_owner" %in% names(acs_data),
           sprintf("%.1f%%", ownership_by_location[location_type == "1_Peripheral", ownership_rate]), "N/A"),
    ifelse("is_owner" %in% names(acs_data),
           sprintf("%.1f%%", ownership_by_location[location_type == "2_Secondary", ownership_rate]), "N/A"),
    ifelse("is_owner" %in% names(acs_data),
           sprintf("%.1f%%", ownership_by_location[location_type == "3_Superstar", ownership_rate]), "N/A")
  ),
  Notes = c(
    "NCHILD in household (lower bound)", "", "", "Unconditional", "",
    "NCHILD == 0 (lower bound)", "", "", "", "",
    "", "", ""
  )
)

print(summary_table)

# Save summary
write.csv(summary_table, file.path(OUTPUT_DIR, "calibration_targets_summary.csv"), row.names = FALSE)

# ======================================================================================
# 12. REGRESSION COEFFICIENTS FOR MODEL
# ======================================================================================
cat("\n\n=== KEY REGRESSION COEFFICIENTS FOR MODEL CALIBRATION ===\n")

cat("\nUnconditional fertility gradient (vs Secondary):\n")
cat(sprintf("  Peripheral: %+.4f (SE: %.4f)\n",
            coef(m1_unconditional)["location_type1_Peripheral"],
            sqrt(vcov(m1_unconditional)["location_type1_Peripheral", "location_type1_Peripheral"])))
cat(sprintf("  Superstar:  %+.4f (SE: %.4f)\n",
            coef(m1_unconditional)["location_type3_Superstar"],
            sqrt(vcov(m1_unconditional)["location_type3_Superstar", "location_type3_Superstar"])))

cat("\nConditional fertility gradient (controlling for educ + income):\n")
cat(sprintf("  Peripheral: %+.4f (SE: %.4f)\n",
            coef(m3_full)["location_type1_Peripheral"],
            sqrt(vcov(m3_full)["location_type1_Peripheral", "location_type1_Peripheral"])))
cat(sprintf("  Superstar:  %+.4f (SE: %.4f)\n",
            coef(m3_full)["location_type3_Superstar"],
            sqrt(vcov(m3_full)["location_type3_Superstar", "location_type3_Superstar"])))

# ======================================================================================
# 13. SAVE ALL OUTPUTS
# ======================================================================================
cat("\n--- 10. SAVING OUTPUTS ---\n")

# Save detailed tables
write.csv(fertility_by_location, file.path(OUTPUT_DIR, "fertility_by_location.csv"), row.names = FALSE)
write.csv(childlessness_by_location, file.path(OUTPUT_DIR, "childlessness_by_location.csv"), row.names = FALSE)
if ("is_owner" %in% names(acs_data)) {
  write.csv(ownership_by_location, file.path(OUTPUT_DIR, "ownership_by_location.csv"), row.names = FALSE)
}
if (exists("asfr_by_location")) {
  write.csv(asfr_by_location, file.path(OUTPUT_DIR, "asfr_by_location.csv"), row.names = FALSE)
  write.csv(tfr_by_location, file.path(OUTPUT_DIR, "tfr_by_location.csv"), row.names = FALSE)
}

# Save regression models
saveRDS(list(
  m1_unconditional = m1_unconditional,
  m2_educ = m2_educ,
  m3_full = m3_full,
  m4_year_fe = m4_year_fe,
  m_childless_1 = m_childless_1,
  m_childless_2 = m_childless_2,
  m_childless_3 = m_childless_3
), file.path(OUTPUT_DIR, "regression_models.rds"))

cat(sprintf("\nAll outputs saved to: %s\n", OUTPUT_DIR))
cat("\n=== ANALYSIS COMPLETE ===\n")
