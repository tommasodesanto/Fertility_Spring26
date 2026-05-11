# =============================================================================
# MSA PANEL SCRIPT (V6.2 - Robust Column Handling in Summaries)
# - ACS 2009: Calculated Fertility (B01001/B13002), B-tables for Pop, CWR, Income, Housing. Education/Family likely NA.
# - ACS 2010+: S1301 Fertility, DP02 for Education/Family (if available), B-tables for Pop, CWR, Income, Housing.
# =============================================================================

# PART 0: SETUP & LIBRARIES
# =============================================================================
cat("PART 0: Setting up libraries and environment...\n")
suppressPackageStartupMessages({
  if (!require(tidycensus)) install.packages("tidycensus"); library(tidycensus)
  if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
  if (!require(scales)) install.packages("scales"); library(scales) # Keep for potential future plots
  if (!require(knitr)) install.packages("knitr"); library(knitr)   # Keep for potential future tables
  if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2) # Keep for potential future plots
})

if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("Census API Key not found. Please set it using census_api_key('YOUR_KEY', install=TRUE).")
}
options(tigris_use_cache = TRUE)

# PART 1: CONFIGURATION
# =============================================================================
cat("PART 1: Configuration...\n")
YEAR_2009_CALC <- 2009
VARS_B01001_DENOM_2009_CWR <- c( # Used for CWR denom and 2009 fertility denom
  f_15_17 = "B01001_030E", f_18_19 = "B01001_031E", f_20 = "B01001_032E",
  f_21 = "B01001_033E", f_22_24 = "B01001_034E", f_25_29 = "B01001_035E",
  f_30_34 = "B01001_036E", f_35_39 = "B01001_037E", f_40_44 = "B01001_038E",
  f_45_49 = "B01001_039E" 
)
VAR_B13002_NUM_2009 <- c(women_with_birth_15_50 = "B13002_002E") # Women 15-50 with a birth

START_YEAR_ACS_S1301_DP02 <- 2010 
LATEST_YEAR_ACS <- 2022 
ALL_TARGET_YEARS <- c(YEAR_2009_CALC, START_YEAR_ACS_S1301_DP02:LATEST_YEAR_ACS)

S1301_FERTILITY_RATE_VAR <- "S1301_C04_001" 

# PART 2: CORRECTED DATA FETCHING WRAPPER
# =============================================================================
cat("PART 2: Defining corrected data fetching wrapper safe_fetch_acs()...\n")
safe_fetch_acs <- function(year_to_fetch, geography_type, variables_list = NULL, table_name = NULL, 
                           survey_endpoint, 
                           cache_table = TRUE, output_format = "wide") {
  tryCatch({
    if (year_to_fetch == 2009 && (survey_endpoint == "acs5/subject" || survey_endpoint == "acs5/profile")) {
      cat(paste0("  Skipping Table '", table_name, "' (", survey_endpoint, ") for 2009 ACS product due to API endpoint issues.\n"))
      return(NULL)
    }
    if (!is.null(table_name)) {
      get_acs(geography = geography_type, table = table_name, year = year_to_fetch, 
              survey = survey_endpoint, cache_table = cache_table, show_call = FALSE)
    } else {
      get_acs(geography = geography_type, variables = variables_list, year = year_to_fetch, 
              survey = survey_endpoint, cache_table = cache_table, output = output_format, show_call = FALSE)
    }
  }, error = function(e) {
    cat(paste0("  Error fetching ACS data (Year: ", year_to_fetch, 
               ", Survey Endpoint: ", survey_endpoint, 
               ", Table/Vars: ", ifelse(!is.null(table_name), table_name, paste(names(variables_list)[1], "...", sep="")),
               "): ", e$message, "\n"))
    return(NULL)
  })
}

# PART 3: MAIN DATA COLLECTION FUNCTION FOR PANEL
# =============================================================================
cat("PART 3: Defining main data collection function get_msa_data_panel()...\n")
get_msa_data_panel <- function(year_to_fetch) {
  cat(paste0("\n--- Fetching data for MSA panel, ACS Product Year: ", year_to_fetch, " ---\n"))
  data_list_for_year <- list()
  base_geographies_df <- NULL

  cat("  Fetching Population (B01003) & Sex/Age for CWR (B01001)...\n")
  vars_pop_cwr_denom <- c(
      VARS_B01001_DENOM_2009_CWR, 
      pop_total_b01003 = "B01003_001E",
      children_0_4_male_b01001 = "B01001_003E",
      children_0_4_female_b01001 = "B01001_027E"
  )
  pop_cwr_data_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, 
                                 geography_type = "metropolitan statistical area/micropolitan statistical area",
                                 variables_list = vars_pop_cwr_denom, survey_endpoint = "acs5", output_format = "wide")

  if (!is.null(pop_cwr_data_raw) && nrow(pop_cwr_data_raw) > 0) {
    base_geographies_df <- pop_cwr_data_raw %>% select(GEOID, NAME) %>% distinct()
    
    # Ensure all expected columns from vars_pop_cwr_denom are present before calculations
    present_pop_vars <- names(vars_pop_cwr_denom)[names(vars_pop_cwr_denom) %in% names(pop_cwr_data_raw)]
    missing_pop_vars <- setdiff(names(vars_pop_cwr_denom), present_pop_vars)
    if(length(missing_pop_vars) > 0) {
        cat(paste0("    Warning: Missing population/CWR variables for year ", year_to_fetch, ": ", paste(missing_pop_vars, collapse=", "), ". CWR might be NA.\n"))
        for(m_var in missing_pop_vars) pop_cwr_data_raw[[m_var]] <- NA_real_ # Add missing as NA
    }

    pop_cwr_processed <- pop_cwr_data_raw %>%
      mutate(
        total_population = as.numeric(pop_total_b01003),
        women_15_49 = rowSums(select(., all_of(names(VARS_B01001_DENOM_2009_CWR)[names(VARS_B01001_DENOM_2009_CWR) %in% names(.)])), na.rm = FALSE),
        children_0_4 = ifelse(all(c("children_0_4_male_b01001", "children_0_4_female_b01001") %in% names(.)),
                              as.numeric(children_0_4_male_b01001) + as.numeric(children_0_4_female_b01001), NA_real_),
        cwr = ifelse(!is.na(women_15_49) & women_15_49 > 0 & !is.na(children_0_4), (children_0_4 / women_15_49) * 1000, NA_real_)
      ) %>% select(GEOID, total_population, women_15_49, children_0_4, cwr)
    data_list_for_year[["pop_cwr"]] <- pop_cwr_processed
  } else {
    cat("    CRITICAL: Failed to fetch B01001/B01003. Attempting to get base GEOID/NAME from another source...\n")
    temp_geo_vars <- c(fallback_var = "B19013_001E") 
    temp_geo_df <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = "metropolitan statistical area/micropolitan statistical area",
                                  variables_list = temp_geo_vars, survey_endpoint = "acs5", output_format="wide")
    if(!is.null(temp_geo_df) && nrow(temp_geo_df) > 0) {
        base_geographies_df <- temp_geo_df %>% select(GEOID, NAME) %>% distinct()
        cat("      Fallback GEOID/NAME fetched.\n")
    } else {
        cat("      Fallback GEOID/NAME also failed. This year's data will likely be empty.\n")
    }
  }

  cat("  Fetching Fertility Rate...\n")
  if (year_to_fetch == YEAR_2009_CALC) {
    cat("    Calculating for 2009 from B13002 (numerator) and B01001 (denominator from pop_cwr)...\n")
    fert_num_data_2009 <- safe_fetch_acs(year_to_fetch = year_to_fetch,
                                        variables_list = VAR_B13002_NUM_2009, survey_endpoint = "acs5", output_format = "wide")
    if (!is.null(fert_num_data_2009) && "pop_cwr" %in% names(data_list_for_year) && !is.null(data_list_for_year[["pop_cwr"]])) {
      # Ensure 'women_15_49' column exists in the pop_cwr component
      if ("women_15_49" %in% names(data_list_for_year[["pop_cwr"]])) {
        fert_data <- data_list_for_year[["pop_cwr"]] %>%
          select(GEOID, women_15_49) %>%
          left_join(fert_num_data_2009 %>% select(GEOID, numerator_women_with_birth = !!sym(names(VAR_B13002_NUM_2009)[1])), by = "GEOID") %>%
          mutate(
            birth_rate_15_50 = ifelse(!is.na(women_15_49) & women_15_49 > 0 & !is.na(numerator_women_with_birth), 
                                      (as.numeric(numerator_women_with_birth) / women_15_49) * 1000, NA_real_),
            birth_rate_15_50_moe = NA_real_, fertility_data_type = "Calculated (WwB/W15-49)"
          ) %>% select(GEOID, birth_rate_15_50, birth_rate_15_50_moe, fertility_data_type)
        data_list_for_year[["fertility"]] <- fert_data
      } else {cat("    'women_15_49' column missing from pop_cwr data for 2009 fertility calculation.\n")}
    } else {cat("    Failed to get components for 2009 calculated fertility rate.\n")}
  } else { 
    cat(paste0("    Using S1301_C04_001 from acs5/subject for ", year_to_fetch, "...\n"))
    fert_s1301_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, table_name = "S1301", survey_endpoint = "acs5/subject")
    if (!is.null(fert_s1301_raw)) {
      if (S1301_FERTILITY_RATE_VAR %in% fert_s1301_raw$variable) {
        fert_data <- fert_s1301_raw %>%
          filter(variable == S1301_FERTILITY_RATE_VAR) %>%
          select(GEOID, birth_rate_15_50 = estimate, birth_rate_15_50_moe = moe) %>%
          mutate(birth_rate_15_50 = as.numeric(birth_rate_15_50), birth_rate_15_50_moe = as.numeric(birth_rate_15_50_moe),
                 fertility_data_type = "S1301 (Births/W15-50)")
        data_list_for_year[["fertility"]] <- fert_data
      } else {cat(paste0("    S1301_C04_001 not found in S1301 table for ", year_to_fetch, ".\n"))}
    } else {cat(paste0("    S1301 table fetch failed for ", year_to_fetch, ".\n"))}
  }
  
  cat("  Fetching Household/Family Structure (DP02)...\n")
  dp02_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, table_name = "DP02", survey_endpoint = "acs5/profile")
  if (!is.null(dp02_raw)) {
    dp02_vars_needed_map <- c(DP02_0002E = "total_families", DP02_0003E = "families_with_children_lt18")
    dp02_data_wide <- dp02_raw %>%
      filter(variable %in% names(dp02_vars_needed_map)) %>%
      select(GEOID, variable, estimate) %>%
      pivot_wider(names_from = variable, values_from = estimate)
      
    # Check if all expected columns exist after pivot
    all_pivoted_cols_exist <- all(names(dp02_vars_needed_map) %in% names(dp02_data_wide))
    if (nrow(dp02_data_wide) > 0 && all_pivoted_cols_exist) {
        dp02_processed <- dp02_data_wide %>%
            transmute(GEOID,
                total_families = as.numeric(DP02_0002E),
                families_with_children_lt18 = as.numeric(DP02_0003E),
                share_families_with_children_lt18 = ifelse(total_families > 0, families_with_children_lt18 / total_families, NA_real_)
            )
        data_list_for_year[["family_structure"]] <- dp02_processed
    } else { cat("    Could not process DP02 family vars for ", year_to_fetch, "; expected variables: ", paste(names(dp02_vars_needed_map), collapse=", "), " not all found after pivot. Vars present: ", paste(names(dp02_data_wide), collapse=", "), "\n") }
  } else { cat("    DP02 table not available or failed to fetch for ", year_to_fetch, ".\n") }

  cat("  Fetching Income data (B19013, B20002)...\n")
  income_vars_b_map <- c(B19013_001E = "median_household_income", B20002_001E = "median_earnings_workers")
  income_data_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, variables_list = names(income_vars_b_map), survey_endpoint = "acs5")
  if (!is.null(income_data_raw)) {
    # Rename using the map values for clarity if successful
    current_names <- names(income_data_raw)
    new_names <- current_names
    for(i in seq_along(current_names)){
        if(current_names[i] %in% names(income_vars_b_map)) {
            new_names[i] <- income_vars_b_map[current_names[i]]
        }
    }
    names(income_data_raw) <- new_names
    # Ensure target columns exist even if renamed variable wasn't fetched
    for(target_col in unname(income_vars_b_map)) if(!target_col %in% names(income_data_raw)) income_data_raw[[target_col]] <- NA_real_

    income_processed <- income_data_raw %>%
      select(GEOID, any_of(unname(income_vars_b_map))) %>% # Select by new names
      mutate(across(all_of(intersect(unname(income_vars_b_map), names(.))), as.numeric)) # Convert only existing columns
    data_list_for_year[["income"]] <- income_processed
  }
  
  cat("  Fetching Educational Attainment...\n")
  if (year_to_fetch > 2009 && !is.null(dp02_raw)) { 
      cat("    Using DP02 for Education for ", year_to_fetch, "...\n")
      edu_denom_var <- "DP02_0059E"
      edu_numer_vars_map <- if (year_to_fetch < 2013) c("DP02_0065E", "DP02_0066E") else c("DP02_0067E")
      all_edu_vars_dp02 <- c(edu_denom_var, edu_numer_vars_map)
      
      dp02_edu_wide <- dp02_raw %>%
        filter(variable %in% all_edu_vars_dp02) %>%
        select(GEOID, variable, estimate) %>%
        pivot_wider(names_from = variable, values_from = estimate)

      all_pivoted_edu_cols_exist <- all(all_edu_vars_dp02 %in% names(dp02_edu_wide))
      if (nrow(dp02_edu_wide) > 0 && all_pivoted_edu_cols_exist) {
          edu_processed <- dp02_edu_wide %>%
            mutate(
              pop_25plus = as.numeric(!!sym(edu_denom_var)),
              pop_bachelors_or_higher_val = if (year_to_fetch < 2013) {
                as.numeric(!!sym(edu_numer_vars_map[1])) + as.numeric(!!sym(edu_numer_vars_map[2]))
              } else {
                as.numeric(!!sym(edu_numer_vars_map[1]))
              }
            ) %>%
            mutate(pct_bachelors_or_higher = ifelse(pop_25plus > 0, (pop_bachelors_or_higher_val / pop_25plus) * 100, NA_real_)) %>%
            select(GEOID, pct_bachelors_or_higher)
          data_list_for_year[["education"]] <- edu_processed
      } else {cat("    Could not process DP02 education vars for ", year_to_fetch, "; expected variables: ", paste(all_edu_vars_dp02, collapse=", "), " not all found after pivot. Vars present: ", paste(names(dp02_edu_wide), collapse=", "), "\n")}
  } else {
    cat(paste0("    Educational attainment from DP02 not attempted or DP02 unavailable for ", year_to_fetch, ".\n"))
  }

  cat("  Fetching Housing data (B25077, B25064)...\n")
  housing_vars_b_map <- c(B25077_001E = "median_home_value", B25064_001E = "median_gross_rent")
  housing_data_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, variables_list = names(housing_vars_b_map), survey_endpoint = "acs5")
  if (!is.null(housing_data_raw)) {
    current_names_h <- names(housing_data_raw)
    new_names_h <- current_names_h
    for(i in seq_along(current_names_h)){
        if(current_names_h[i] %in% names(housing_vars_b_map)) {
            new_names_h[i] <- housing_vars_b_map[current_names_h[i]]
        }
    }
    names(housing_data_raw) <- new_names_h
    for(target_col_h in unname(housing_vars_b_map)) if(!target_col_h %in% names(housing_data_raw)) housing_data_raw[[target_col_h]] <- NA_real_
    
    housing_processed <- housing_data_raw %>%
      select(GEOID, any_of(unname(housing_vars_b_map))) %>%
      mutate(across(all_of(intersect(unname(housing_vars_b_map), names(.))), as.numeric))
    data_list_for_year[["housing"]] <- housing_processed
  }

  if (is.null(base_geographies_df) || nrow(base_geographies_df) == 0) {
      cat(paste0("    CRITICAL: No base geographies (GEOID, NAME) captured for year ", year_to_fetch, ". Cannot create panel row.\n"))
      return(NULL)
  }
  final_data_for_year <- base_geographies_df
  
  data_list_for_year[["pop_cwr"]] <- data_list_for_year[["pop_cwr"]] %>% select(-any_of("women_15_49")) # Drop intermediate column if not needed in final panel

  for (df_name in names(data_list_for_year)) {
    current_df <- data_list_for_year[[df_name]]
    if (!is.null(current_df) && nrow(current_df) > 0 && "GEOID" %in% names(current_df)) {
      cols_in_final <- names(final_data_for_year)
      cols_in_current <- names(current_df)
      join_by_col <- "GEOID"
      cols_to_add <- setdiff(cols_in_current, cols_in_final) 
      if (length(cols_to_add) > 0) {
         final_data_for_year <- final_data_for_year %>% 
           left_join(current_df %>% select(all_of(c(join_by_col, cols_to_add))), by = join_by_col)
      }
    }
  }
  
  final_data_for_year <- final_data_for_year %>%
    mutate(year_acs_product = as.integer(year_to_fetch))
  
  # Check if any actual data columns were added beyond GEOID, NAME, year_acs_product
  if (ncol(final_data_for_year) <= 3 && nrow(final_data_for_year) > 0) { # Only GEOID, NAME, year_acs_product
      cat(paste0("--- Only base geographies and year retained for MSA panel, ACS Product Year: ", year_to_fetch, 
             ". MSAs: ", nrow(final_data_for_year), ". This means all specific data fetches might have failed or returned no new columns.\n"))
  } else if (nrow(final_data_for_year) > 0) {
    cat(paste0("--- Successfully processed data for MSA panel, ACS Product Year: ", year_to_fetch, 
             ". MSAs: ", nrow(final_data_for_year), " ---\n"))
  } else {
    cat(paste0("--- No data successfully processed for MSA panel, ACS Product Year: ", year_to_fetch, " ---\n"))
    return(NULL)
  }
  return(final_data_for_year)
}

# PART 4: DATA AGGREGATION & PANEL CREATION (Using get_msa_data_panel)
# =============================================================================
cat("\nPART 4: Creating Full MSA Panel...\n")
cat(paste0("Target ACS Product Years for panel: ", paste(ALL_TARGET_YEARS, collapse = ", "), "\n"))

msa_panel_list_full <- lapply(ALL_TARGET_YEARS, get_msa_data_panel)
msa_panel_list_full <- msa_panel_list_full[!sapply(msa_panel_list_full, is.null)] 

if (length(msa_panel_list_full) == 0) {
  stop("FATAL: No data fetched for ANY year for the full MSA panel.")
}
msa_panel_full_raw <- bind_rows(msa_panel_list_full)

cat(paste0("\nFull raw panel created with ", nrow(msa_panel_full_raw), " MSA-year observations.\n"))
if(nrow(msa_panel_full_raw) > 0) {
    cat(paste0("ACS Product Years in full panel: ", paste(sort(unique(msa_panel_full_raw$year_acs_product)), collapse = ", "), "\n"))
}

# PART 5: FINAL CLEANING & SAVING
# =============================================================================
cat("\nPART 5: Final Cleaning, Feature Engineering & Saving...\n")

msa_panel_final <- msa_panel_full_raw 
if("year_acs_product" %in% names(msa_panel_final)) {
    msa_panel_final <- msa_panel_final %>% rename(year = year_acs_product)
}

intended_numeric_cols <- c(
  "total_population", "women_15_49", "children_0_4", "cwr", 
  "birth_rate_15_50", "birth_rate_15_50_moe",
  "total_families", "families_with_children_lt18", "share_families_with_children_lt18",
  "median_household_income", "median_earnings_workers", 
  "pct_bachelors_or_higher", 
  "median_home_value", "median_gross_rent"
)
actual_cols_to_convert <- intersect(intended_numeric_cols, names(msa_panel_final))
cat(paste0("Actual columns found in panel to convert to numeric: ", paste(actual_cols_to_convert, collapse=", "),"\n"))
if (length(actual_cols_to_convert) > 0) {
  msa_panel_final <- msa_panel_final %>%
    mutate(across(all_of(actual_cols_to_convert), ~as.numeric(as.character(.))))
}

if (all(c("median_home_value", "median_household_income") %in% names(msa_panel_final))) {
  msa_panel_final <- msa_panel_final %>%
    mutate(price_to_income_ratio = ifelse(median_household_income > 0 & median_home_value > 0, 
                                          median_home_value / median_household_income, NA_real_))
} else {
  cat("Warning: median_home_value or median_household_income missing, cannot calculate price_to_income_ratio.\n")
  if (!"price_to_income_ratio" %in% names(msa_panel_final)) msa_panel_final$price_to_income_ratio <- NA_real_
}

output_main_dir_v6_2 <- "msa_panel_output_v6_2"
dir.create(output_main_dir_v6_2, showWarnings = FALSE)
output_file_final <- file.path(output_main_dir_v6_2, "msa_panel_v6_2_2009calc_onwards.rds")
write_csv(msa_panel_final, gsub(".rds$", ".csv", output_file_final))
saveRDS(msa_panel_final, output_file_final)
cat(paste0("\nFinal MSA panel V6.2 saved to: ", output_file_final, " (and .csv)\n"))

cat("\n--- Summary of Key Variables (Sample from latest year) ---\n")
if (nrow(msa_panel_final) > 0 && "year" %in% names(msa_panel_final)) {
  latest_year_in_panel <- max(msa_panel_final$year, na.rm = TRUE)
  if(is.finite(latest_year_in_panel)){
    summary_cols <- c("GEOID", "NAME", "year", "total_population", "birth_rate_15_50", 
                      "fertility_data_type", "cwr", "share_families_with_children_lt18", 
                      "median_household_income", "median_earnings_workers", 
                      "pct_bachelors_or_higher", "median_home_value", "median_gross_rent",
                      "price_to_income_ratio")
    cols_present_for_summary <- intersect(summary_cols, names(msa_panel_final))
    if(length(cols_present_for_summary) > 2) {
        summary_sample_latest <- msa_panel_final %>% 
        filter(year == latest_year_in_panel) %>%
        select(all_of(cols_present_for_summary)) %>%
        summary()
        print(summary_sample_latest)
    } else {cat("Not enough relevant columns present for a meaningful summary of the latest year.\n")}
  } else {cat("Could not determine latest year for summary.\n")}
} else {cat("Panel is empty or 'year' column missing for summary.\n")}

cat("\n--- Check 2009 calculated fertility data specifically ---\n")
if (nrow(msa_panel_final) > 0 && YEAR_2009_CALC %in% msa_panel_final$year) {
    summary_2009_cols <- c("GEOID", "NAME", "year", "birth_rate_15_50", "fertility_data_type", 
                           "total_population", "cwr", "median_household_income", 
                           "median_earnings_workers", "pct_bachelors_or_higher") 
    cols_present_for_2009_summary <- intersect(summary_2009_cols, names(msa_panel_final))
    if(length(cols_present_for_2009_summary) > 2){
        summary_2009_calc_df <- msa_panel_final %>%
            filter(year == YEAR_2009_CALC) %>%
            select(all_of(cols_present_for_2009_summary)) %>%
            arrange(desc(total_population)) 
        print(head(summary_2009_calc_df, 10))
        
        # Count MSAs with 2009 rate if birth_rate_15_50 column exists
        if("birth_rate_15_50" %in% names(msa_panel_final)){
            n_2009_rate <- msa_panel_final %>% 
                            filter(year == YEAR_2009_CALC & !is.na(birth_rate_15_50)) %>% nrow()
            cat(paste0("Number of MSAs with calculated 2009 rate: ", n_2009_rate, "\n"))
        } else {
            cat("Column 'birth_rate_15_50' not found for 2009 NA count.\n")
        }
    } else {cat("Not enough relevant columns present for a meaningful summary of 2009.\n")}
} else {cat("No data for 2009 found in the final panel for summary.\n")}

cat("\n--- V6.2 PANEL SCRIPT COMPLETE ---\n")