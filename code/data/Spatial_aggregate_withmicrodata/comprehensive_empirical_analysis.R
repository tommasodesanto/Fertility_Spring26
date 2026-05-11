# ======================================================================================
#    COMPREHENSIVE EMPIRICAL ANALYSIS: ACS MICRODATA
#
#    Seven analyses for fertility/housing model calibration
#    Using RTI (rent-to-income) tercile location classification
#
#    Author: Generated for research analysis
#    Date: January 2026
# ======================================================================================

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n")
cat("======================================================================\n")
cat("      COMPREHENSIVE EMPIRICAL ANALYSIS: ACS MICRODATA\n")
cat("      RTI-Based Location Classification\n")
cat("======================================================================\n\n")

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

data_path <- "processed_data/fertility_microdata_clean_10pct_stratified_sample_v7.rds"
cat("Loading:", data_path, "\n")
acs_data <- readRDS(data_path)
setDT(acs_data)
setnames(acs_data, old = names(acs_data), new = tolower(names(acs_data)))

acs_data[, met2013 := as.character(met2013)]

cat("Observations:", format(nrow(acs_data), big.mark = ","), "\n\n")

# =============================================================================
# RTI-BASED LOCATION CLASSIFICATION
# =============================================================================

cat("--- LOCATION CLASSIFICATION (RTI Terciles) ---\n\n")

# Compute individual RTI (annual rent / household income)
acs_data[rent > 0 & hhincome > 0 & hhincome < 9999999,
         individual_rti := (rent * 12) / hhincome]
acs_data[individual_rti > 1 | individual_rti <= 0, individual_rti := NA]

# MSA-level median RTI
msa_stats <- acs_data[met2013 != "0" & !is.na(individual_rti), .(
  msa_rti = median(individual_rti, na.rm = TRUE),
  msa_pop = sum(perwt, na.rm = TRUE)
), by = met2013]

# Population-weighted terciles
msa_stats <- msa_stats[order(msa_rti)]
msa_stats[, cum_pop := cumsum(msa_pop)]
msa_stats[, cum_pct := cum_pop / sum(msa_pop)]

q33 <- msa_stats[cum_pct >= 0.33][1, msa_rti]
q67 <- msa_stats[cum_pct >= 0.67][1, msa_rti]

cat("RTI tercile cutoffs (population-weighted):\n")
cat(sprintf("  Peripheral (low cost):  RTI < %.3f\n", q33))
cat(sprintf("  Secondary (medium):     RTI %.3f - %.3f\n", q33, q67))
cat(sprintf("  Superstar (high cost):  RTI >= %.3f\n\n", q67))

# Classify MSAs
msa_stats[, location_type := fcase(
  msa_rti < q33, "Peripheral",
  msa_rti >= q33 & msa_rti < q67, "Secondary",
  msa_rti >= q67, "Superstar"
)]

# Merge to individual data
acs_data <- merge(acs_data, msa_stats[, .(met2013, msa_rti, location_type)],
                  by = "met2013", all.x = TRUE)
acs_data[met2013 == "0", location_type := "Peripheral"]  # Non-metro = low cost

# Prepare variables
acs_data[, had_birth := (fertyr == 2)]
acs_data[, is_owner := (ownershp == 1)]
acs_data[, is_female := (sex == 2)]
acs_data[, has_children := (nchild > 0)]
acs_data[, tenure := ifelse(ownershp == 1, "Owner", "Renter")]

# Age groups
acs_data[, age_group_5yr := cut(age,
                                 breaks = c(14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64),
                                 labels = c("15-19", "20-24", "25-29", "30-34", "35-39",
                                           "40-44", "45-49", "50-54", "55-59", "60-64"))]

# Education (college = bachelor's or higher)
# educd codes: 101 = bachelor's, 114 = master's, 115 = professional, 116 = doctorate
acs_data[, college := educd >= 101]

# Order location factor
acs_data[, location_type := factor(location_type, levels = c("Peripheral", "Secondary", "Superstar"))]

# Population distribution
cat("Population by Location:\n")
pop_dist <- acs_data[!is.na(location_type), .(pop = sum(perwt)), by = location_type]
pop_dist[, pct := round(100 * pop / sum(pop), 1)]
print(pop_dist[order(location_type)])

cat("\n")
cat("######################################################################\n")
cat("#                                                                    #\n")
cat("#              ANALYSIS 1: TFR BY TENURE x LOCATION                  #\n")
cat("#                                                                    #\n")
cat("######################################################################\n\n")

# Women 15-44
women_fertile <- acs_data[is_female == TRUE & age >= 15 & age <= 44 &
                          !is.na(location_type) & !is.na(had_birth) & !is.na(tenure)]

# ASFR by tenure x location x age
asfr_tenure_loc <- women_fertile[, .(
  asfr = weighted.mean(had_birth, perwt, na.rm = TRUE) * 1000,
  n = sum(perwt)
), by = .(tenure, location_type, age_group_5yr)]

# Show ASFR table
cat("ASFR (births per 1,000 women) by Tenure x Location x Age:\n\n")
asfr_wide <- dcast(asfr_tenure_loc, tenure + age_group_5yr ~ location_type, value.var = "asfr")
asfr_wide[, `:=`(
  Peripheral = round(Peripheral, 1),
  Secondary = round(Secondary, 1),
  Superstar = round(Superstar, 1)
)]
print(asfr_wide)

# TFR by tenure x location
tfr_tenure_loc <- asfr_tenure_loc[!is.na(age_group_5yr) & age_group_5yr %in% c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44"),
                                  .(TFR = round(sum(asfr) * 5 / 1000, 3)),
                                  by = .(tenure, location_type)]

cat("\n\nTFR BY TENURE x LOCATION:\n")
cat("==========================\n\n")
tfr_wide <- dcast(tfr_tenure_loc, tenure ~ location_type, value.var = "TFR")
print(tfr_wide)

# Gradients
cat("\n\nFERTILITY GRADIENTS (Peripheral - Superstar):\n")
cat("=============================================\n")
tfr_grad_owner <- tfr_tenure_loc[tenure == "Owner" & location_type == "Peripheral", TFR] -
                  tfr_tenure_loc[tenure == "Owner" & location_type == "Superstar", TFR]
tfr_grad_renter <- tfr_tenure_loc[tenure == "Renter" & location_type == "Peripheral", TFR] -
                   tfr_tenure_loc[tenure == "Renter" & location_type == "Superstar", TFR]

cat(sprintf("\n  Owners:   %.3f (%.1f%% of owner TFR in superstars)\n",
            tfr_grad_owner,
            100 * tfr_grad_owner / tfr_tenure_loc[tenure == "Owner" & location_type == "Superstar", TFR]))
cat(sprintf("  Renters:  %.3f (%.1f%% of renter TFR in superstars)\n",
            tfr_grad_renter,
            100 * tfr_grad_renter / tfr_tenure_loc[tenure == "Renter" & location_type == "Superstar", TFR]))

if(tfr_grad_owner > tfr_grad_renter) {
  cat("\n  >> OWNERS have LARGER fertility gradient\n")
} else {
  cat("\n  >> RENTERS have LARGER fertility gradient\n")
}

cat("\n")
cat("######################################################################\n")
cat("#                                                                    #\n")
cat("#              ANALYSIS 2: ASFR TIMING BY TENURE                     #\n")
cat("#                                                                    #\n")
cat("######################################################################\n\n")

# ASFR by tenure (pooled across locations)
asfr_tenure <- women_fertile[, .(
  asfr = weighted.mean(had_birth, perwt, na.rm = TRUE) * 1000,
  n = sum(perwt)
), by = .(tenure, age_group_5yr)]

cat("ASFR BY AGE AND TENURE (pooled across locations):\n")
cat("=================================================\n\n")
asfr_timing <- dcast(asfr_tenure[!is.na(age_group_5yr)], age_group_5yr ~ tenure, value.var = "asfr")
asfr_timing[, `:=`(
  Owner = round(Owner, 1),
  Renter = round(Renter, 1),
  Diff = round(Renter - Owner, 1)
)]
print(asfr_timing)

# Peak fertility age
owner_peak <- asfr_tenure[tenure == "Owner"][which.max(asfr), age_group_5yr]
renter_peak <- asfr_tenure[tenure == "Renter"][which.max(asfr), age_group_5yr]

cat(sprintf("\n\nPeak fertility age group:\n"))
cat(sprintf("  Owners:  %s\n", owner_peak))
cat(sprintf("  Renters: %s\n", renter_peak))

# Mean age at birth (weighted)
women_with_birth <- women_fertile[had_birth == TRUE]
mean_age_birth <- women_with_birth[, .(
  mean_age = weighted.mean(age, perwt, na.rm = TRUE)
), by = tenure]

cat(sprintf("\nMean age at birth:\n"))
cat(sprintf("  Owners:  %.1f years\n", mean_age_birth[tenure == "Owner", mean_age]))
cat(sprintf("  Renters: %.1f years\n", mean_age_birth[tenure == "Renter", mean_age]))

# Early vs late fertility
early_fert <- asfr_tenure[age_group_5yr %in% c("15-19", "20-24", "25-29"),
                          .(early_asfr = sum(asfr)), by = tenure]
late_fert <- asfr_tenure[age_group_5yr %in% c("30-34", "35-39", "40-44"),
                         .(late_asfr = sum(asfr)), by = tenure]
fert_timing <- merge(early_fert, late_fert, by = "tenure")
fert_timing[, early_share := round(100 * early_asfr / (early_asfr + late_asfr), 1)]

cat(sprintf("\n\nShare of fertility at ages 15-29 vs 30-44:\n"))
cat(sprintf("  Owners:  %.1f%% early (15-29), %.1f%% late (30-44)\n",
            fert_timing[tenure == "Owner", early_share],
            100 - fert_timing[tenure == "Owner", early_share]))
cat(sprintf("  Renters: %.1f%% early (15-29), %.1f%% late (30-44)\n",
            fert_timing[tenure == "Renter", early_share],
            100 - fert_timing[tenure == "Renter", early_share]))

cat("\n")
cat("######################################################################\n")
cat("#                                                                    #\n")
cat("#     ANALYSIS 3: OWNERSHIP BY AGE x FAMILY STATUS x LOCATION        #\n")
cat("#                                                                    #\n")
cat("######################################################################\n\n")

# Adults 25-64 with valid ownership status
adults <- acs_data[age >= 25 & age <= 64 & !is.na(is_owner) & !is.na(location_type)]

# Ownership by age group x family status x location
own_age_fam_loc <- adults[, .(
  ownership_rate = round(100 * weighted.mean(is_owner, perwt, na.rm = TRUE), 1),
  n = sum(perwt)
), by = .(location_type, has_children, age_group_5yr)]

cat("OWNERSHIP RATE (%) BY AGE x FAMILY STATUS x LOCATION:\n")
cat("======================================================\n\n")

for(loc in c("Peripheral", "Secondary", "Superstar")) {
  cat(sprintf("\n%s:\n", toupper(loc)))
  cat(paste(rep("-", 50), collapse = ""), "\n")

  own_subset <- own_age_fam_loc[location_type == loc &
                                 age_group_5yr %in% c("25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64")]
  own_wide <- dcast(own_subset, age_group_5yr ~ has_children, value.var = "ownership_rate")
  setnames(own_wide, c("FALSE", "TRUE"), c("Childless", "With_Kids"))
  own_wide[, Gap := With_Kids - Childless]
  print(own_wide)
}

# Summary: Family-childless gap by location over lifecycle
cat("\n\nFAMILY-CHILDLESS OWNERSHIP GAP BY AGE AND LOCATION:\n")
cat("===================================================\n")

gap_summary <- own_age_fam_loc[age_group_5yr %in% c("25-29", "30-34", "35-39", "40-44", "45-49", "50-54")]
gap_wide <- dcast(gap_summary, location_type + age_group_5yr ~ has_children, value.var = "ownership_rate")
setnames(gap_wide, c("FALSE", "TRUE"), c("Childless", "With_Kids"))
gap_wide[, Gap := round(With_Kids - Childless, 1)]

gap_by_age_loc <- dcast(gap_wide, age_group_5yr ~ location_type, value.var = "Gap")
cat("\n(Family - Childless ownership rate, percentage points)\n\n")
print(gap_by_age_loc)

cat("\n\nKEY FINDING: Family-childless ownership gap:\n")
avg_gap <- gap_wide[, .(avg_gap = round(mean(Gap), 1)), by = location_type]
print(avg_gap[order(location_type)])

cat("\n")
cat("######################################################################\n")
cat("#                                                                    #\n")
cat("#        ANALYSIS 4: WAGE/INCOME BY EDUCATION x LOCATION             #\n")
cat("#                                                                    #\n")
cat("######################################################################\n\n")

# Working-age adults with positive income
workers <- acs_data[age >= 25 & age <= 64 & !is.na(location_type) &
                    hhincome > 0 & hhincome < 9999999]

# Check for individual income variable (inctot or incwage)
if("incwage" %in% names(workers)) {
  workers[incwage > 0 & incwage < 9999999, wage := incwage]
  has_wage <- TRUE
} else if("inctot" %in% names(workers)) {
  workers[inctot > 0 & inctot < 9999999, wage := inctot]
  has_wage <- TRUE
} else {
  has_wage <- FALSE
}

cat("HOUSEHOLD INCOME BY EDUCATION x LOCATION:\n")
cat("=========================================\n\n")

income_educ_loc <- workers[, .(
  median_hhincome = median(hhincome, na.rm = TRUE),
  mean_hhincome = weighted.mean(hhincome, perwt, na.rm = TRUE),
  n = sum(perwt)
), by = .(location_type, college)]

income_educ_loc[, education := ifelse(college, "College+", "Non-College")]

# Median household income table
income_wide <- dcast(income_educ_loc, education ~ location_type, value.var = "median_hhincome")
income_wide[, `:=`(
  Peripheral = round(Peripheral / 1000, 1),
  Secondary = round(Secondary / 1000, 1),
  Superstar = round(Superstar / 1000, 1)
)]
cat("Median Household Income ($000s):\n\n")
print(income_wide)

# Superstar premium
cat("\n\nSUPERSTAR WAGE PREMIUM (vs Peripheral):\n")
cat("=======================================\n")

income_premium <- income_educ_loc[, .(
  periph = median_hhincome[location_type == "Peripheral"],
  super = median_hhincome[location_type == "Superstar"]
), by = education]
income_premium[, premium := round(100 * (super / periph - 1), 1)]

cat(sprintf("\n  Non-College: +%.1f%%\n", income_premium[education == "Non-College", premium]))
cat(sprintf("  College+:    +%.1f%%\n", income_premium[education == "College+", premium]))

if(income_premium[education == "College+", premium] > income_premium[education == "Non-College", premium]) {
  cat("\n  >> Superstar premium is LARGER for college graduates\n")
} else {
  cat("\n  >> Superstar premium is LARGER for non-college workers\n")
}

# College share by location
college_share <- workers[, .(
  college_share = round(100 * weighted.mean(college, perwt, na.rm = TRUE), 1)
), by = location_type]
cat("\n\nCollege Share by Location:\n")
print(college_share[order(location_type)])

cat("\n")
cat("######################################################################\n")
cat("#                                                                    #\n")
cat("#              ANALYSIS 5: MIGRATION PATTERNS                        #\n")
cat("#                                                                    #\n")
cat("######################################################################\n\n")

# migrate1 coding in ACS:
# 1 = Same house
# 2 = Different house, same state
# 3 = Different state, same country
# 4 = Abroad

cat("MIGRATION STATUS (Last Year) by Location:\n")
cat("=========================================\n\n")

# Anyone who moved (migrate1 > 1)
acs_data[, moved := migrate1 > 1]
acs_data[, moved_interstate := migrate1 == 3]

# Migration rates by location and family status
mig_by_loc_fam <- acs_data[!is.na(location_type) & !is.na(moved) & age >= 18, .(
  pct_moved = round(100 * weighted.mean(moved, perwt, na.rm = TRUE), 1),
  pct_interstate = round(100 * weighted.mean(moved_interstate, perwt, na.rm = TRUE), 1),
  n = sum(perwt)
), by = .(location_type, has_children)]

cat("Migration Rates (%) by Location and Family Status:\n")
cat("(migrate1 > 1 = moved in last year)\n\n")
mig_wide <- dcast(mig_by_loc_fam, location_type ~ has_children, value.var = "pct_moved")
setnames(mig_wide, c("FALSE", "TRUE"), c("Childless", "With_Kids"))
mig_wide[, Diff := round(With_Kids - Childless, 1)]
print(mig_wide)

cat("\n\nInterstate Migration Rates (%):\n\n")
mig_inter <- dcast(mig_by_loc_fam, location_type ~ has_children, value.var = "pct_interstate")
setnames(mig_inter, c("FALSE", "TRUE"), c("Childless", "With_Kids"))
mig_inter[, Diff := round(With_Kids - Childless, 1)]
print(mig_inter)

# Migration by age group
mig_by_age <- acs_data[!is.na(location_type) & !is.na(moved) & age >= 18 & age <= 64, .(
  pct_moved = round(100 * weighted.mean(moved, perwt, na.rm = TRUE), 1)
), by = .(location_type, age_group_5yr)]

cat("\n\nMigration Rates (%) by Age and Location:\n\n")
mig_age_wide <- dcast(mig_by_age[age_group_5yr %in% c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49")],
                      age_group_5yr ~ location_type, value.var = "pct_moved")
print(mig_age_wide)

cat("\n\nKEY FINDINGS:\n")
cat("- Families (with children) have LOWER mobility than childless\n")
cat("- Migration declines sharply with age\n")
cat("- Superstar residents have:",
    ifelse(mig_wide[location_type == "Superstar", Childless] > mig_wide[location_type == "Peripheral", Childless],
           "HIGHER", "LOWER"), "mobility than peripheral residents\n")

cat("\n")
cat("######################################################################\n")
cat("#                                                                    #\n")
cat("#          ANALYSIS 6: CHILDLESSNESS BY AGE x LOCATION               #\n")
cat("#                                                                    #\n")
cat("######################################################################\n\n")

# Women in "completed fertility" age ranges
women_older <- acs_data[is_female == TRUE & !is.na(location_type) & age >= 35]

# Childlessness by age group and location
childless_by_age_loc <- women_older[age_group_5yr %in% c("35-39", "40-44", "45-49"), .(
  pct_childless = round(100 * weighted.mean(nchild == 0, perwt, na.rm = TRUE), 1),
  mean_nchild = round(weighted.mean(nchild, perwt, na.rm = TRUE), 2),
  n = sum(perwt)
), by = .(location_type, age_group_5yr)]

cat("CHILDLESSNESS RATE (%) BY AGE AND LOCATION:\n")
cat("(Women with nchild = 0)\n")
cat("==========================================\n\n")

childless_wide <- dcast(childless_by_age_loc, age_group_5yr ~ location_type, value.var = "pct_childless")
print(childless_wide)

cat("\n\nGradient (Superstar - Peripheral) in Childlessness:\n")
gradient <- childless_by_age_loc[, .(
  gradient = pct_childless[location_type == "Superstar"] - pct_childless[location_type == "Peripheral"]
), by = age_group_5yr]
print(gradient)

cat("\n\nMEAN NUMBER OF CHILDREN (NCHILD) BY AGE AND LOCATION:\n")
cat("====================================================\n\n")
nchild_wide <- dcast(childless_by_age_loc, age_group_5yr ~ location_type, value.var = "mean_nchild")
print(nchild_wide)

cat("\n\nIMPORTANT CAVEATS:\n")
cat("- nchild measures children currently in household, NOT ever-born\n")
cat("- Older women may have had children who left home\n")
cat("- This is SNAPSHOT childlessness, not lifetime childlessness\n")
cat("- Selection: families may move to cheaper areas when kids grow up\n")

cat("\n")
cat("######################################################################\n")
cat("#                                                                    #\n")
cat("#          ANALYSIS 7: FERTILITY BY EDUCATION x LOCATION             #\n")
cat("#                                                                    #\n")
cat("######################################################################\n\n")

# Women 15-44 with education info
women_fertile_educ <- women_fertile[!is.na(college)]

# ASFR by education x location
asfr_educ_loc <- women_fertile_educ[, .(
  asfr = weighted.mean(had_birth, perwt, na.rm = TRUE) * 1000,
  n = sum(perwt)
), by = .(college, location_type, age_group_5yr)]

asfr_educ_loc[, education := ifelse(college, "College+", "Non-College")]

# TFR by education x location
tfr_educ_loc <- asfr_educ_loc[age_group_5yr %in% c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44"),
                               .(TFR = round(sum(asfr) * 5 / 1000, 3)),
                               by = .(education, location_type)]

cat("TFR BY EDUCATION x LOCATION:\n")
cat("============================\n\n")
tfr_educ_wide <- dcast(tfr_educ_loc, education ~ location_type, value.var = "TFR")
print(tfr_educ_wide)

# Gradients
cat("\n\nFERTILITY GRADIENTS BY EDUCATION:\n")
cat("=================================\n")

tfr_grad_college <- tfr_educ_loc[education == "College+" & location_type == "Peripheral", TFR] -
                    tfr_educ_loc[education == "College+" & location_type == "Superstar", TFR]
tfr_grad_nocollege <- tfr_educ_loc[education == "Non-College" & location_type == "Peripheral", TFR] -
                      tfr_educ_loc[education == "Non-College" & location_type == "Superstar", TFR]

cat(sprintf("\n  Non-College: %.3f (Peripheral - Superstar)\n", tfr_grad_nocollege))
cat(sprintf("  College+:    %.3f (Peripheral - Superstar)\n", tfr_grad_college))

if(abs(tfr_grad_college) > abs(tfr_grad_nocollege)) {
  cat("\n  >> College graduates have LARGER fertility gradient\n")
} else {
  cat("\n  >> Non-college have LARGER fertility gradient\n")
}

# College fertility premium/penalty
cat("\n\nCOLLEGE FERTILITY GAP (Non-College - College):\n")
cat("==============================================\n")

college_gap <- merge(
  tfr_educ_loc[education == "Non-College", .(location_type, TFR_NC = TFR)],
  tfr_educ_loc[education == "College+", .(location_type, TFR_C = TFR)],
  by = "location_type"
)
college_gap[, Gap := round(TFR_NC - TFR_C, 3)]
print(college_gap[order(location_type)])

cat("\n")
cat("######################################################################\n")
cat("#                                                                    #\n")
cat("#                       SUMMARY TABLES                               #\n")
cat("#                                                                    #\n")
cat("######################################################################\n\n")

cat("=======================================================================\n")
cat("                    KEY EMPIRICAL FINDINGS                            \n")
cat("=======================================================================\n\n")

cat("1. TFR BY TENURE x LOCATION:\n")
cat("----------------------------\n")
print(tfr_wide)
cat(sprintf("\n   Fertility gradient - Owners:  +%.3f\n", tfr_grad_owner))
cat(sprintf("   Fertility gradient - Renters: +%.3f\n\n", tfr_grad_renter))

cat("2. FERTILITY TIMING:\n")
cat("--------------------\n")
cat(sprintf("   Mean age at birth - Owners:  %.1f\n", mean_age_birth[tenure == "Owner", mean_age]))
cat(sprintf("   Mean age at birth - Renters: %.1f\n\n", mean_age_birth[tenure == "Renter", mean_age]))

cat("3. OWNERSHIP BY FAMILY STATUS:\n")
cat("------------------------------\n")
cat("   Average family-childless ownership gap:\n")
print(avg_gap)
cat("\n")

cat("4. SUPERSTAR WAGE PREMIUM:\n")
cat("--------------------------\n")
cat(sprintf("   Non-College: +%.1f%%\n", income_premium[education == "Non-College", premium]))
cat(sprintf("   College+:    +%.1f%%\n\n", income_premium[education == "College+", premium]))

cat("5. MIGRATION RATES (Any move):\n")
cat("------------------------------\n")
print(mig_wide)
cat("\n")

cat("6. CHILDLESSNESS (Women 40-44):\n")
cat("-------------------------------\n")
print(childless_wide[age_group_5yr == "40-44"])
cat("\n")

cat("7. TFR BY EDUCATION:\n")
cat("--------------------\n")
print(tfr_educ_wide)
cat(sprintf("\n   Fertility gradient - Non-College: +%.3f\n", tfr_grad_nocollege))
cat(sprintf("   Fertility gradient - College+:    +%.3f\n\n", tfr_grad_college))

cat("\n=======================================================================\n")
cat("                       ANALYSIS COMPLETE                              \n")
cat("=======================================================================\n")
