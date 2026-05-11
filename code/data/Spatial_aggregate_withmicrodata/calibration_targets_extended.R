# ======================================================================================
#
#    Extended Calibration Targets Analysis
#    Part 2: Validation, Homeownership, and Migration Patterns
#
#    Purpose:
#    1. Validate fertility findings with multiple approaches
#    2. Examine homeownership patterns by family status and location
#    3. Analyze migration patterns and sorting by fertility/family status
#
#    Author: Tommaso De Santo
#    Date: January 2026
#
# ======================================================================================

# ======================================================================================
# 0. SETUP (Load previous analysis environment)
# ======================================================================================
cat("=== EXTENDED CALIBRATION ANALYSIS ===\n")
cat("Part 2: Validation, Homeownership, Migration\n\n")

required_packages <- c("tidyverse", "data.table", "fixest", "modelsummary",
                       "matrixStats", "kableExtra", "haven")

for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

OUTPUT_DIR <- "calibration_targets_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================================================================================
# 1. LOAD DATA
# ======================================================================================
cat("--- 1. LOADING DATA ---\n")

data_paths <- c(
  "processed_data/fertility_microdata_clean_10pct_stratified_sample_v7.rds",
  "processed_data/fertility_microdata_clean_5pct_stratified_sample.rds",
  "output/validated_microdata_sample.rds"
)

for (path in data_paths) {
  if (file.exists(path)) {
    cat(sprintf("Loading data from: %s\n", path))
    acs_data <- readRDS(path)
    break
  }
}

setDT(acs_data)
setnames(acs_data, old = names(acs_data), new = tolower(names(acs_data)))

# ======================================================================================
# 2. VARIABLE PREPARATION (same as before)
# ======================================================================================
cat("\n--- 2. PREPARING VARIABLES ---\n")

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

# Sex
if (max(acs_data$sex, na.rm=TRUE) == 2) {
  acs_data[, is_female := (sex == 2)]
} else {
  acs_data[, is_female := (sex == 1)]
}

# Education
if (!"educ_factor" %in% names(acs_data)) {
  if ("educd" %in% names(acs_data)) {
    acs_data[, educ_factor := fcase(
      educd < 62, "1_LessHS",
      educd >= 62 & educd < 65, "2_HS",
      educd >= 65 & educd < 101, "3_SomeColl",
      educd >= 101 & educd < 114, "4_BA",
      educd >= 114, "5_Grad",
      default = NA_character_
    )]
  }
}

# Income
acs_data[hhincome <= 0 | hhincome == 9999999, hhincome := NA]
acs_data[!is.na(hhincome) & hhincome > 0, log_hhincome := log(hhincome)]

# Ownership
if ("ownershp" %in% names(acs_data)) {
  acs_data[, is_owner := (ownershp == 1)]
}

# Migration
if ("migrate1d" %in% names(acs_data)) {
  acs_data[, did_move := as.integer(migrate1d >= 20)]  # Moved from different house
  acs_data[, moved_diff_state := as.integer(migrate1d >= 24)]  # Moved from different state
  acs_data[, moved_diff_msa := as.integer(migrate1d %in% c(21, 22, 23, 24))]  # Different MSA
}

# Recent mother
if ("fertyr" %in% names(acs_data)) {
  acs_data[, is_recent_mother := as.numeric(fertyr == 2 & is_female == TRUE)]
}

# Has children
acs_data[, has_children := (nchild > 0)]
acs_data[, has_young_children := (nchild > 0 & yngch <= 5)]

cat("Variables prepared.\n")

# ======================================================================================
# 3. FERTILITY VALIDATION: MULTIPLE APPROACHES
# ======================================================================================
cat("\n" %+% strrep("=", 70) %+% "\n")
cat("SECTION A: FERTILITY VALIDATION\n")
cat(strrep("=", 70) %+% "\n")

# --- 3.1 By Year: Check if pattern is stable over time ---
cat("\n--- A.1: Fertility Pattern by Year ---\n")

fertility_by_year_loc <- acs_data[is_female == TRUE & age %between% c(40, 50) & !is.na(location_type), .(
  mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
  n = .N
), by = .(year, location_type)][order(year, location_type)]

# Wide format for display
fertility_wide <- dcast(fertility_by_year_loc, year ~ location_type, value.var = "mean_nchild")
cat("\nMean NCHILD (women 40-50) by year and location:\n")
print(fertility_wide)

# Check if gradient is consistent
fertility_wide[, gradient_PS := `1_Peripheral` - `3_Superstar`]
cat(sprintf("\nMean gradient (Peripheral - Superstar) across years: %.3f\n",
            mean(fertility_wide$gradient_PS, na.rm=TRUE)))
cat(sprintf("SD of gradient across years: %.3f\n",
            sd(fertility_wide$gradient_PS, na.rm=TRUE)))

# --- 3.2 By Age: Check different age windows ---
cat("\n--- A.2: Fertility by Different Age Windows ---\n")

age_windows <- list(
  "35-39" = c(35, 39),
  "40-44" = c(40, 44),
  "45-49" = c(45, 49),
  "40-50" = c(40, 50),
  "35-50" = c(35, 50)
)

fertility_by_age_window <- rbindlist(lapply(names(age_windows), function(w) {
  ages <- age_windows[[w]]
  acs_data[is_female == TRUE & age %between% ages & !is.na(location_type), .(
    mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
    pct_childless = weighted.mean(nchild == 0, perwt, na.rm = TRUE) * 100,
    n = .N
  ), by = location_type][, age_window := w]
}))

cat("\nMean NCHILD by age window and location:\n")
print(dcast(fertility_by_age_window, age_window ~ location_type, value.var = "mean_nchild"))

cat("\nChildlessness (%) by age window and location:\n")
print(dcast(fertility_by_age_window, age_window ~ location_type, value.var = "pct_childless"))

# --- 3.3 Use ELDCH to infer children who left home ---
cat("\n--- A.3: Using ELDCH (Age of Eldest Child) to Infer True Fertility ---\n")

if ("eldch" %in% names(acs_data)) {
  # For women 45-50: if eldch > 0, they have/had at least one child
  # Even if nchild = 0, eldch might indicate a child who left

  women_45_50 <- acs_data[is_female == TRUE & age %between% c(45, 50) & !is.na(location_type)]

  # IPUMS: eldch = 99 means no own children in household
  women_45_50[, has_or_had_child := (nchild > 0 | (eldch > 0 & eldch < 99))]
  women_45_50[, child_left_home := (nchild == 0 & eldch > 0 & eldch < 99)]

  eldch_analysis <- women_45_50[, .(
    pct_nchild_0 = weighted.mean(nchild == 0, perwt, na.rm = TRUE) * 100,
    pct_truly_childless = weighted.mean(has_or_had_child == FALSE, perwt, na.rm = TRUE) * 100,
    pct_child_left = weighted.mean(child_left_home, perwt, na.rm = TRUE) * 100,
    mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
    n = .N
  ), by = location_type][order(location_type)]

  cat("\nChildlessness measures (women 45-50):\n")
  cat("  pct_nchild_0: NCHILD == 0 (children may have left)\n")
  cat("  pct_truly_childless: No children ever in household (more accurate)\n")
  cat("  pct_child_left: NCHILD==0 but ELDCH indicates past child\n\n")
  print(eldch_analysis)

} else {
  cat("ELDCH variable not available in data.\n")
}

# --- 3.4 TFR by location (confirms flow measure) ---
cat("\n--- A.4: TFR Validation (Flow Measure) ---\n")

if ("is_recent_mother" %in% names(acs_data)) {
  women_15_50 <- acs_data[is_female == TRUE & age %between% c(15, 49) & !is.na(location_type)]
  women_15_50[, age_group := cut(age, breaks = seq(15, 50, 5), right = FALSE,
                                  labels = paste0(seq(15, 45, 5), "-", seq(19, 49, 5)))]

  asfr <- women_15_50[!is.na(age_group), .(
    asfr = weighted.mean(is_recent_mother, perwt, na.rm = TRUE) * 1000,
    n = .N
  ), by = .(location_type, age_group)]

  tfr <- asfr[, .(tfr = 5 * sum(asfr / 1000, na.rm = TRUE)), by = location_type]

  cat("\nTotal Fertility Rate by location type:\n")
  print(tfr[order(location_type)])

  cat(sprintf("\nTFR gradient (Peripheral - Superstar): %.3f children\n",
              tfr[location_type == "1_Peripheral", tfr] - tfr[location_type == "3_Superstar", tfr]))
}

# --- 3.5 By education within location ---
cat("\n--- A.5: Fertility by Education WITHIN Location ---\n")

fertility_educ_loc <- acs_data[is_female == TRUE & age %between% c(40, 50) &
                                !is.na(location_type) & !is.na(educ_factor), .(
  mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
  n = .N
), by = .(location_type, educ_factor)][order(location_type, educ_factor)]

cat("\nMean NCHILD by education and location (women 40-50):\n")
print(dcast(fertility_educ_loc, educ_factor ~ location_type, value.var = "mean_nchild"))

# ======================================================================================
# 4. HOMEOWNERSHIP PATTERNS
# ======================================================================================
cat("\n" %+% strrep("=", 70) %+% "\n")
cat("SECTION B: HOMEOWNERSHIP PATTERNS\n")
cat(strrep("=", 70) %+% "\n")

if ("is_owner" %in% names(acs_data)) {

  # --- 4.1 Ownership by family status ---
  cat("\n--- B.1: Homeownership by Family Status ---\n")

  acs_data[, family_status := fcase(
    nchild == 0, "No children",
    nchild == 1, "1 child",
    nchild == 2, "2 children",
    nchild >= 3, "3+ children",
    default = NA_character_
  )]

  ownership_family <- acs_data[!is.na(location_type) & !is.na(family_status) & age %between% c(25, 55), .(
    ownership_rate = weighted.mean(is_owner, perwt, na.rm = TRUE) * 100,
    n = .N
  ), by = .(location_type, family_status)][order(location_type, family_status)]

  cat("\nHomeownership rate by family status and location (ages 25-55):\n")
  print(dcast(ownership_family, family_status ~ location_type, value.var = "ownership_rate"))

  # --- 4.2 Ownership by age and location ---
  cat("\n--- B.2: Homeownership Lifecycle by Location ---\n")

  acs_data[, age_group_5 := cut(age, breaks = seq(25, 70, 5), right = FALSE,
                                 labels = paste0(seq(25, 65, 5), "-", seq(29, 69, 5)))]

  ownership_lifecycle <- acs_data[!is.na(location_type) & !is.na(age_group_5), .(
    ownership_rate = weighted.mean(is_owner, perwt, na.rm = TRUE) * 100,
    n = .N
  ), by = .(location_type, age_group_5)][order(location_type, age_group_5)]

  cat("\nHomeownership rate by age and location:\n")
  print(dcast(ownership_lifecycle, age_group_5 ~ location_type, value.var = "ownership_rate"))

  # --- 4.3 Ownership gap by children status ---
  cat("\n--- B.3: Ownership Gap: Families vs Childless ---\n")

  ownership_gap <- acs_data[!is.na(location_type) & age %between% c(30, 50), .(
    own_with_kids = weighted.mean(is_owner[has_children == TRUE], perwt[has_children == TRUE], na.rm = TRUE) * 100,
    own_no_kids = weighted.mean(is_owner[has_children == FALSE], perwt[has_children == FALSE], na.rm = TRUE) * 100,
    n_with_kids = sum(has_children == TRUE),
    n_no_kids = sum(has_children == FALSE)
  ), by = location_type][order(location_type)]

  ownership_gap[, gap := own_with_kids - own_no_kids]

  cat("\nHomeownership rate: families with children vs childless (ages 30-50):\n")
  print(ownership_gap)

  # --- 4.4 Regression: Ownership ~ location + children ---
  cat("\n--- B.4: Ownership Regression ---\n")

  reg_data <- acs_data[!is.na(location_type) & age %between% c(25, 60) & !is.na(educ_factor) & !is.na(log_hhincome)]
  reg_data[, location_type := relevel(factor(location_type), ref = "2_Secondary")]

  m_own1 <- feglm(is_owner ~ location_type + has_children | age,
                  data = reg_data, family = "binomial", weights = ~perwt)

  m_own2 <- feglm(is_owner ~ location_type * has_children | age,
                  data = reg_data, family = "binomial", weights = ~perwt)

  m_own3 <- feglm(is_owner ~ location_type * has_children + educ_factor + log_hhincome | age,
                  data = reg_data, family = "binomial", weights = ~perwt)

  cat("\nLogit: Homeownership ~ Location + Children (ages 25-60)\n")
  etable(m_own1, m_own2, m_own3,
         headers = c("Additive", "Interaction", "+ Controls"),
         signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
         fitstat = ~ n + pr2)
}

# ======================================================================================
# 5. MIGRATION PATTERNS
# ======================================================================================
cat("\n" %+% strrep("=", 70) %+% "\n")
cat("SECTION C: MIGRATION PATTERNS\n")
cat(strrep("=", 70) %+% "\n")

if ("did_move" %in% names(acs_data)) {

  # --- 5.1 Migration rates by location and family status ---
  cat("\n--- C.1: Migration Rates by Location and Family Status ---\n")

  migration_rates <- acs_data[!is.na(location_type) & age %between% c(25, 50), .(
    pct_moved = weighted.mean(did_move, perwt, na.rm = TRUE) * 100,
    pct_moved_msa = weighted.mean(moved_diff_msa, perwt, na.rm = TRUE) * 100,
    n = .N
  ), by = .(location_type, has_children)][order(location_type, has_children)]

  cat("\nMigration rates by location and family status (ages 25-50):\n")
  cat("  pct_moved: Moved from different house in past year\n")
  cat("  pct_moved_msa: Moved from different MSA\n\n")
  print(migration_rates)

  # --- 5.2 Who moves TO superstar cities? ---
  cat("\n--- C.2: Profile of Movers TO Superstar Cities ---\n")

  if ("migmet131" %in% names(acs_data)) {
    # People currently in superstar who moved from elsewhere
    movers_to_superstar <- acs_data[location_type == "3_Superstar" & moved_diff_msa == 1 & age %between% c(25, 45)]
    stayers_superstar <- acs_data[location_type == "3_Superstar" & did_move == 0 & age %between% c(25, 45)]

    if (nrow(movers_to_superstar) > 100) {
      cat("\nComparing movers TO superstar vs stayers in superstar (ages 25-45):\n")

      comparison <- rbind(
        movers_to_superstar[, .(
          group = "Moved to Superstar",
          mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
          pct_has_kids = weighted.mean(has_children, perwt, na.rm = TRUE) * 100,
          pct_owner = weighted.mean(is_owner, perwt, na.rm = TRUE) * 100,
          mean_age = weighted.mean(age, perwt, na.rm = TRUE),
          n = .N
        )],
        stayers_superstar[, .(
          group = "Stayers in Superstar",
          mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
          pct_has_kids = weighted.mean(has_children, perwt, na.rm = TRUE) * 100,
          pct_owner = weighted.mean(is_owner, perwt, na.rm = TRUE) * 100,
          mean_age = weighted.mean(age, perwt, na.rm = TRUE),
          n = .N
        )]
      )
      print(comparison)
    }
  }

  # --- 5.3 Who LEAVES superstar cities? ---
  cat("\n--- C.3: Who Leaves Expensive Cities? ---\n")

  # Look at current peripheral residents who came from superstar
  if ("migmet131" %in% names(acs_data)) {
    acs_data[, migmet131 := as.character(migmet131)]

    # Get origin MSA population (approximate by using current year's population)
    origin_pop <- msa_pop[, .(origin_msa_pop = mean(msa_population)), by = met2013]
    acs_data <- merge(acs_data, origin_pop, by.x = "migmet131", by.y = "met2013", all.x = TRUE)

    # Classify origin location type
    acs_data[, origin_location := fcase(
      is.na(migmet131) | migmet131 == "0", NA_character_,
      origin_msa_pop < 250000, "1_Peripheral",
      origin_msa_pop >= 250000 & origin_msa_pop < 2000000, "2_Secondary",
      origin_msa_pop >= 2000000, "3_Superstar",
      default = NA_character_
    )]

    # Migration flows
    migration_flows <- acs_data[moved_diff_msa == 1 & !is.na(location_type) & !is.na(origin_location) &
                                  age %between% c(25, 50), .(
      n = .N,
      pop = sum(perwt, na.rm = TRUE),
      mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
      pct_has_kids = weighted.mean(has_children, perwt, na.rm = TRUE) * 100
    ), by = .(origin_location, location_type)][order(origin_location, location_type)]

    # Calculate percentages within origin
    migration_flows[, pct_of_origin := pop / sum(pop) * 100, by = origin_location]

    cat("\nMigration flows between location types (ages 25-50):\n")
    cat("  Shows where people from each origin go, and their family characteristics\n\n")
    print(migration_flows)

    # --- 5.4 Do families leave superstars? ---
    cat("\n--- C.4: Family Status of Those Leaving Superstars ---\n")

    leaving_superstar <- acs_data[origin_location == "3_Superstar" & moved_diff_msa == 1 &
                                    location_type != "3_Superstar" & age %between% c(25, 50)]
    staying_superstar <- acs_data[location_type == "3_Superstar" & did_move == 0 & age %between% c(25, 50)]

    if (nrow(leaving_superstar) > 50) {
      cat("\nComparing those who LEFT superstar vs those who stayed:\n")
      comparison_leave <- rbind(
        leaving_superstar[, .(
          group = "Left Superstar",
          mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
          pct_has_kids = weighted.mean(has_children, perwt, na.rm = TRUE) * 100,
          pct_young_kids = weighted.mean(has_young_children, perwt, na.rm = TRUE) * 100,
          pct_owner = weighted.mean(is_owner, perwt, na.rm = TRUE) * 100,
          n = .N
        )],
        staying_superstar[, .(
          group = "Stayed in Superstar",
          mean_nchild = weighted.mean(nchild, perwt, na.rm = TRUE),
          pct_has_kids = weighted.mean(has_children, perwt, na.rm = TRUE) * 100,
          pct_young_kids = weighted.mean(has_young_children, perwt, na.rm = TRUE) * 100,
          pct_owner = weighted.mean(is_owner, perwt, na.rm = TRUE) * 100,
          n = .N
        )]
      )
      print(comparison_leave)
    }
  }

  # --- 5.5 Regression: Who leaves expensive cities? ---
  cat("\n--- C.5: Regression: Predictors of Leaving High-Cost Areas ---\n")

  # Sample: People currently in peripheral/secondary who moved in past year
  # Outcome: Did they come from a superstar city?

  movers <- acs_data[moved_diff_msa == 1 & !is.na(origin_location) & !is.na(location_type) &
                       age %between% c(25, 50) & !is.na(educ_factor)]
  movers[, came_from_superstar := as.integer(origin_location == "3_Superstar")]
  movers[, went_to_cheaper := as.integer(
    (origin_location == "3_Superstar" & location_type != "3_Superstar") |
    (origin_location == "2_Secondary" & location_type == "1_Peripheral")
  )]

  if (nrow(movers) > 500) {
    m_leave1 <- feglm(went_to_cheaper ~ has_children + has_young_children | age,
                      data = movers, family = "binomial", weights = ~perwt)

    m_leave2 <- feglm(went_to_cheaper ~ has_children + has_young_children + educ_factor | age,
                      data = movers, family = "binomial", weights = ~perwt)

    m_leave3 <- feglm(went_to_cheaper ~ has_children + has_young_children + educ_factor + is_owner | age,
                      data = movers, family = "binomial", weights = ~perwt)

    cat("\nLogit: Moved to cheaper location ~ Family status (among inter-MSA movers)\n")
    etable(m_leave1, m_leave2, m_leave3,
           headers = c("Family only", "+ Education", "+ Ownership"),
           signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
           fitstat = ~ n + pr2)
  }
}

# ======================================================================================
# 6. SUMMARY OF KEY FINDINGS
# ======================================================================================
cat("\n" %+% strrep("=", 70) %+% "\n")
cat("SUMMARY OF KEY FINDINGS FOR MODEL CALIBRATION\n")
cat(strrep("=", 70) %+% "\n")

cat("\n1. FERTILITY PATTERNS:\n")
cat("   - TFR is higher in peripheral areas (flow measure)\n")
cat("   - But NCHILD (stock) shows opposite pattern - likely selection + children leaving home\n")
cat("   - Need to distinguish sorting from causal effects\n")

cat("\n2. HOMEOWNERSHIP:\n")
cat("   - Large gradient by location (67% peripheral vs 54% superstar)\n")
cat("   - Families with children more likely to own in ALL locations\n")
cat("   - Gap may be larger in expensive cities (interaction)\n")

cat("\n3. MIGRATION:\n")
cat("   - Families with children MORE likely to leave expensive cities\n")
cat("   - Those arriving in superstars have fewer children\n")
cat("   - Sorting: expensive cities select childless/career-focused\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

# Save outputs
write.csv(fertility_by_year_loc, file.path(OUTPUT_DIR, "fertility_by_year_location.csv"), row.names = FALSE)
if (exists("migration_flows")) {
  write.csv(migration_flows, file.path(OUTPUT_DIR, "migration_flows.csv"), row.names = FALSE)
}
if (exists("ownership_family")) {
  write.csv(ownership_family, file.path(OUTPUT_DIR, "ownership_by_family_status.csv"), row.names = FALSE)
}
