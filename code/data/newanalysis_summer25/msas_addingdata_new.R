# =============================================================================
# MSA PANEL SCRIPT (V6.5.1 - Corrected B11012 Logic, Table Fetches, Children in HHs, Enhanced Affordability)
# - ACS 2009: Calculated Fertility (B13002), B-tables. DP02/DP04 attempted (may be skipped by safe_fetch_acs).
# - ACS 2010+: S1301 Fertility, DP02, DP04, B-tables.
# - Corrected survey_endpoint logic for get_acs table calls.
# - Corrected B11012 variable usage for comprehensive child shares.
# =============================================================================

# PART 0: SETUP & LIBRARIES
# =============================================================================
cat("PART 0: Setting up libraries and environment...\n")
suppressPackageStartupMessages({
  if (!require(tidycensus)) install.packages("tidycensus"); library(tidycensus)
  if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
  if (!require(scales)) install.packages("scales"); library(scales)
  if (!require(knitr)) install.packages("knitr"); library(knitr)
  if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
})

if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("Census API Key not found. Please set it using census_api_key('YOUR_KEY', install=TRUE) and restarting R.")
}
options(tigris_use_cache = TRUE)

# PART 1: CONFIGURATION
# =============================================================================
cat("PART 1: Configuration...\n")
YEAR_2009_CALC <- 2009
VARS_B01001_WOMEN_15_49_FOR_CWR <- c(
  f_15_17_cwr = "B01001_030E", f_18_19_cwr = "B01001_031E", f_20_cwr = "B01001_032E",
  f_21_cwr = "B01001_033E", f_22_24_cwr = "B01001_034E", f_25_29_cwr = "B01001_035E",
  f_30_34_cwr = "B01001_036E", f_35_39_cwr = "B01001_037E", f_40_44_cwr = "B01001_038E",
  f_45_49_cwr = "B01001_039E"
)
VARS_B13002_FERTILITY_2009 <- c( # Used for 2009 calculated fertility rate
  denom_women_15_50_fert2009 = "B13002_001E", # Total Women 15-50
  num_women_w_birth_15_50_fert2009 = "B13002_002E"  # Women 15-50 with a birth
)
START_YEAR_ACS_S1301_DP02 <- 2010
LATEST_YEAR_ACS <- 2022 #Or Sys.Date() %>% format("%Y") %>% as.numeric() - 2 for more dynamic latest year
ALL_TARGET_YEARS <- c(YEAR_2009_CALC, START_YEAR_ACS_S1301_DP02:LATEST_YEAR_ACS)
# For quicker testing:
# ALL_TARGET_YEARS <- c(2009, 2022)

S1301_FERTILITY_RATE_VAR <- "S1301_C04_001"

# PART 2: DATA FETCHING WRAPPER
# =============================================================================
cat("PART 2: Defining data fetching wrapper safe_fetch_acs()...\n")
safe_fetch_acs <- function(year_to_fetch, geography_type, variables_list = NULL, table_name = NULL,
                           survey_endpoint,
                           cache_table = TRUE, output_format = "wide") {
  tryCatch({
    if (year_to_fetch == 2009 && survey_endpoint %in% c("acs5/subject", "acs5/profile")) { # Defensive
      cat(paste0("  Skipping Table '", table_name, "' (", survey_endpoint, ") for 2009 per script logic (DP/S tables not available for 2009 product).\n"))
      return(NULL)
    }
    if (!is.null(table_name)) {
      # For DP/S tables, survey_endpoint should be "acs5" and tidycensus infers /profile or /subject
      # For B/C tables by table name, survey_endpoint is also "acs5"
      get_acs(geography = geography_type, table = table_name, year = year_to_fetch,
              survey = survey_endpoint, cache_table = cache_table, show_call = FALSE)
    } else { # Fetching specific variables
      get_acs(geography = geography_type, variables = variables_list, year = year_to_fetch,
              survey = survey_endpoint, cache_table = cache_table, output = output_format, show_call = FALSE)
    }
  }, error = function(e) {
    cat(paste0("  Error fetching (Year: ", year_to_fetch, ", Geo: ", geography_type, ", Survey: ", survey_endpoint,
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
  current_geography_type <- "metropolitan statistical area/micropolitan statistical area"
  
  base_survey_for_tables <- "acs5" # Use this for S, DP tables when 'table_name' is specified
  survey_for_variables <- "acs5"   # Use this for B-tables (and others) when 'variables_list' is specified
  
  # 1. Population & CWR Components (B01001, B01003)
  cat("  1. Population & CWR Components (B01001, B01003)...\n")
  vars_pop_cwr <- c(VARS_B01001_WOMEN_15_49_FOR_CWR,
                    pop_total_alias = "B01003_001E",
                    children_0_4_m_alias = "B01001_003E", children_0_4_f_alias = "B01001_027E")
  pop_cwr_data_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                                     variables_list = vars_pop_cwr, survey_endpoint = survey_for_variables, output_format = "wide")
  if (!is.null(pop_cwr_data_raw) && nrow(pop_cwr_data_raw) > 0) {
    base_geographies_df <- pop_cwr_data_raw %>% select(GEOID, NAME) %>% distinct()
    expected_cwr_denom_cols <- names(VARS_B01001_WOMEN_15_49_FOR_CWR)
    present_cwr_denom_cols <- intersect(expected_cwr_denom_cols, names(pop_cwr_data_raw))
    for(m_var in setdiff(expected_cwr_denom_cols, present_cwr_denom_cols)) pop_cwr_data_raw[[m_var]] <- NA_real_
    pop_cwr_processed <- pop_cwr_data_raw %>%
      mutate(
        total_population = if("pop_total_alias" %in% names(.)) as.numeric(pop_total_alias) else NA_real_,
        women_15_49_cwr_denom = rowSums(select(., all_of(present_cwr_denom_cols)), na.rm = FALSE),
        children_0_4 = if(all(c("children_0_4_m_alias", "children_0_4_f_alias") %in% names(.))) as.numeric(children_0_4_m_alias) + as.numeric(children_0_4_f_alias) else NA_real_,
        cwr = ifelse(women_15_49_cwr_denom > 0 & !is.na(children_0_4), (children_0_4 / women_15_49_cwr_denom) * 1000, NA_real_)
      ) %>% select(GEOID, total_population, children_0_4, cwr)
    data_list_for_year[["pop_cwr"]] <- pop_cwr_processed
  } else {
    cat("    CRITICAL: Pop/CWR (B01001/3) failed. Fallback for GEOID/NAME...\n")
    temp_geo_vars <- c(fb="B19013_001E")
    temp_geo_df <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                                  variables_list = temp_geo_vars, survey_endpoint = survey_for_variables, output_format="wide")
    if(!is.null(temp_geo_df) && nrow(temp_geo_df) > 0) base_geographies_df <- temp_geo_df %>% select(GEOID, NAME) %>% distinct()
  }
  
  # 2. Fertility Rate
  cat("  2. Fertility Rate...\n")
  if (year_to_fetch == YEAR_2009_CALC) {
    fert_comp_09 <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                                   variables_list = VARS_B13002_FERTILITY_2009, survey_endpoint = survey_for_variables, output_format = "wide")
    if (!is.null(fert_comp_09) && all(names(VARS_B13002_FERTILITY_2009) %in% names(fert_comp_09))) {
      data_list_for_year[["fertility"]] <- fert_comp_09 %>%
        transmute(GEOID, birth_rate_15_50 = ifelse(as.numeric(denom_women_15_50_fert2009) > 0 & !is.na(num_women_w_birth_15_50_fert2009),
                                                   (as.numeric(num_women_w_birth_15_50_fert2009) / as.numeric(denom_women_15_50_fert2009)) * 1000, NA_real_),
                  birth_rate_15_50_moe = NA_real_, fertility_data_type = "Calculated (B13002)")
    } else { cat("    Failed B13002 for 2009 fertility.\n") }
  } else {
    fert_s1301 <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                                 table_name = "S1301", survey_endpoint = base_survey_for_tables)
    if (!is.null(fert_s1301) && S1301_FERTILITY_RATE_VAR %in% fert_s1301$variable) {
      data_list_for_year[["fertility"]] <- fert_s1301 %>% filter(variable == S1301_FERTILITY_RATE_VAR) %>%
        select(GEOID, birth_rate_15_50 = estimate, birth_rate_15_50_moe = moe) %>%
        mutate(across(c(birth_rate_15_50, birth_rate_15_50_moe), as.numeric), fertility_data_type = "S1301 Direct")
    } else { cat(paste0("    S1301 or var ",S1301_FERTILITY_RATE_VAR," failed for ",year_to_fetch,".\n")) }
  }
  
  # 3. DP02 Data (Family Structure, Education)
  cat("  3. DP02 (Family Structure, Education)...\n")
  # For 2009, DP tables are not reliably available via "acs5" survey endpoint in tidycensus for table calls.
  # safe_fetch_acs handles this by returning NULL if survey_endpoint is "acs5/profile" for 2009.
  # Since we now use base_survey_for_tables ("acs5"), let's add specific logic to skip DP for 2009.
  if (year_to_fetch == YEAR_2009_CALC) {
    dp02_raw <- NULL # Explicitly set to NULL for 2009
    cat("    DP02 data fetch skipped for 2009.\n")
  } else {
    dp02_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                               table_name = "DP02", survey_endpoint = base_survey_for_tables)
  }
  
  if (!is.null(dp02_raw) && nrow(dp02_raw) > 0) {
    dp02_fam_vars <- c(DP02_0002E = "total_families", DP02_0003E = "families_w_child_lt18_dp02")
    dp02_fam_wide <- dp02_raw %>% filter(variable %in% names(dp02_fam_vars)) %>% select(GEOID, variable, estimate) %>% pivot_wider(names_from=variable, values_from=estimate, values_fn=first)
    if (nrow(dp02_fam_wide) > 0 && all(names(dp02_fam_vars) %in% names(dp02_fam_wide))) {
      data_list_for_year[["family_structure"]] <- dp02_fam_wide %>%
        transmute(GEOID, total_families = as.numeric(DP02_0002E), families_w_child_lt18 = as.numeric(DP02_0003E),
                  share_families_with_children_lt18 = ifelse(total_families > 0 & !is.na(families_w_child_lt18), (families_w_child_lt18 / total_families) * 100, NA_real_))
    } else { cat("    DP02 family vars processing failed (or vars not found after pivot).\n") }
    
    edu_denom <- "DP02_0059E"; edu_num_vars <- if (year_to_fetch < 2013) c("DP02_0065E","DP02_0066E") else "DP02_0067E"
    all_edu_vars <- c(edu_denom, edu_num_vars)
    dp02_edu_wide <- dp02_raw %>% filter(variable %in% all_edu_vars) %>% select(GEOID,variable,estimate) %>% pivot_wider(names_from=variable,values_from=estimate,values_fn=first)
    if (nrow(dp02_edu_wide) > 0 && all(all_edu_vars %in% names(dp02_edu_wide))) {
      data_list_for_year[["education"]] <- dp02_edu_wide %>%
        mutate(pop25_edu_denom = as.numeric(!!sym(edu_denom)),
               bach_higher_edu_num = if(year_to_fetch<2013) {
                 val1 <- as.numeric(!!sym(edu_num_vars[1])); val2 <- as.numeric(!!sym(edu_num_vars[2]))
                 if(is.na(val1) | is.na(val2)) NA_real_ else val1 + val2
               } else { as.numeric(!!sym(edu_num_vars)) }) %>%
        transmute(GEOID, pct_bachelors_or_higher = ifelse(pop25_edu_denom > 0 & !is.na(bach_higher_edu_num), (bach_higher_edu_num/pop25_edu_denom)*100, NA_real_))
    } else { cat("    DP02 education vars processing failed (or vars not found after pivot).\n") }
  } else if (year_to_fetch != YEAR_2009_CALC) {
    cat("    DP02 table unavailable/empty for a non-2009 year.\n")
  }
  
  # 4. Household Composition (Children) from B11012 - CORRECTED LOGIC
  cat("  4. Household Composition (Children - B11012)...\n")
  vars_hh_children_corrected <- c(
    total_households_b11012 = "B11012_001E",
    mcf_with_children_lt18 = "B11012_003E", mhh_with_children_lt18 = "B11012_009E", fhh_with_children_lt18 = "B11012_015E",
    mcf_children_lt6_only = "B11012_004E", mcf_children_lt6_and_6to17 = "B11012_005E",
    mhh_children_lt6_only = "B11012_010E", mhh_children_lt6_and_6to17 = "B11012_011E",
    fhh_children_lt6_only = "B11012_016E", fhh_children_lt6_and_6to17 = "B11012_017E"
  )
  hh_children_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                                    variables_list = vars_hh_children_corrected, survey_endpoint = survey_for_variables, output_format = "wide")
  
  if (!is.null(hh_children_raw) && nrow(hh_children_raw) > 0) {
    for(v_alias in names(vars_hh_children_corrected)) if(!v_alias %in% names(hh_children_raw)) hh_children_raw[[v_alias]] <- NA_real_
    hh_children_raw <- hh_children_raw %>% mutate(across(all_of(names(vars_hh_children_corrected)), as.numeric))
    
    # Optional Diagnostic Block
    if (year_to_fetch %in% c(ALL_TARGET_YEARS[1], ALL_TARGET_YEARS[length(ALL_TARGET_YEARS)])) {
      cat(paste0("DIAGNOSTIC (B11012 - Year: ", year_to_fetch, "): Raw values for B11012 components (up to 5 rows):\n"))
      diag_cols_subset <- c("GEOID", "NAME", "total_households_b11012", "mcf_with_children_lt18", "mhh_with_children_lt18", "fhh_with_children_lt18", "mcf_children_lt6_only", "mhh_children_lt6_only", "fhh_children_lt6_only")
      diag_cols_present <- intersect(diag_cols_subset, names(hh_children_raw))
      if (length(diag_cols_present) > ("GEOID" %in% diag_cols_present) + ("NAME" %in% diag_cols_present) + 2) { # Check if enough cols exist
        print(head(hh_children_raw %>% select(all_of(diag_cols_present)) %>% filter(!is.na(total_households_b11012)), 5))
      } else { cat("  Insufficient columns for B11012 diagnostic print.\n") }
    }
    
    data_list_for_year[["household_composition_children"]] <- hh_children_raw %>%
      mutate(
        total_hh_with_children_lt18 = rowSums(select(., all_of(c("mcf_with_children_lt18", "mhh_with_children_lt18", "fhh_with_children_lt18"))), na.rm = FALSE),
        total_hh_with_children_lt6 = rowSums(select(., all_of(c("mcf_children_lt6_only", "mcf_children_lt6_and_6to17", "mhh_children_lt6_only", "mhh_children_lt6_and_6to17", "fhh_children_lt6_only", "fhh_children_lt6_and_6to17"))), na.rm = FALSE)
      ) %>%
      transmute( GEOID,
                 share_households_with_children_lt18 = ifelse(total_households_b11012 > 0 & !is.na(total_hh_with_children_lt18), (total_hh_with_children_lt18 / total_households_b11012) * 100, NA_real_),
                 share_households_with_children_lt6 = ifelse(total_households_b11012 > 0 & !is.na(total_hh_with_children_lt6), (total_hh_with_children_lt6 / total_households_b11012) * 100, NA_real_)
      )
  } else { cat("    Household composition (B11012) fetch failed.\n") }
  
  
  # 5. Income (B19013, B20002, B25119)
  cat("  5. Income (B19013, B20002, B25119)...\n")
  inc_vars <- c(median_household_income="B19013_001E", median_earnings_workers="B20002_001E", median_renter_household_income="B25119_003E")
  inc_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                            variables_list = inc_vars, survey_endpoint = survey_for_variables, output_format = "wide")
  if (!is.null(inc_raw) && nrow(inc_raw) > 0) {
    for(v_alias in names(inc_vars)) if(!v_alias %in% names(inc_raw)) inc_raw[[v_alias]] <- NA_real_
    data_list_for_year[["income"]] <- inc_raw %>% select(GEOID, any_of(names(inc_vars))) %>%
      mutate(across(all_of(intersect(names(inc_vars), names(.))), as.numeric))
  } else { cat("    Income fetch failed.\n")}
  
  # 6. Housing Values (B25077, B25064)
  cat("  6. Housing Values (B25077, B25064)...\n")
  hous_vars <- c(median_home_value="B25077_001E", median_gross_rent="B25064_001E")
  hous_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                             variables_list = hous_vars, survey_endpoint = survey_for_variables, output_format = "wide")
  if (!is.null(hous_raw) && nrow(hous_raw) > 0) {
    for(v_alias in names(hous_vars)) if(!v_alias %in% names(hous_raw)) hous_raw[[v_alias]] <- NA_real_
    data_list_for_year[["housing_values"]] <- hous_raw %>% select(GEOID, any_of(names(hous_vars))) %>%
      mutate(across(all_of(intersect(names(hous_vars), names(.))), as.numeric))
  } else { cat("    Housing values fetch failed.\n")}
  
  # 7. Housing Cost Burden (DP04)
  cat("  7. Housing Cost Burden (DP04)...\n")
  if (year_to_fetch == YEAR_2009_CALC) {
    dp04_raw <- NULL # Explicitly set to NULL for 2009
    cat("    DP04 data fetch skipped for 2009.\n")
  } else {
    dp04_raw <- safe_fetch_acs(year_to_fetch = year_to_fetch, geography_type = current_geography_type,
                               table_name="DP04", survey_endpoint = base_survey_for_tables)
  }
  
  if(!is.null(dp04_raw) && nrow(dp04_raw) > 0){
    dp04_vars_to_fetch <- c( "DP04_0133PE", "DP04_0134PE", "DP04_0100PE", "DP04_0101PE", "DP04_0106PE", "DP04_0107PE", "DP04_0045E",  "DP04_0126E", "DP04_0095E",  "DP04_0102E" )
    dp04_burden_wide <- dp04_raw %>% filter(variable %in% dp04_vars_to_fetch) %>%
      select(GEOID, variable, estimate) %>% pivot_wider(names_from = variable, values_from = estimate, values_fn = first)
    for(v_code in dp04_vars_to_fetch) if(!v_code %in% names(dp04_burden_wide)) dp04_burden_wide[[v_code]] <- NA_real_
    
    if(nrow(dp04_burden_wide) > 0 && any(sapply(dp04_burden_wide[intersect(dp04_vars_to_fetch, names(dp04_burden_wide))], function(x) !all(is.na(x)))) ){
      dp04_burden_processed <- dp04_burden_wide %>%
        mutate(across(all_of(intersect(dp04_vars_to_fetch, names(.))), as.numeric)) %>%
        mutate(
          pct_renters_cost_burdened_30plus = if(all(c("DP04_0133PE", "DP04_0134PE") %in% names(.))) ifelse(!is.na(DP04_0133PE) & !is.na(DP04_0134PE), DP04_0133PE + DP04_0134PE, NA_real_) else NA_real_,
          pct_owners_w_mortgage_cost_burdened_30plus = if(all(c("DP04_0100PE", "DP04_0101PE") %in% names(.))) ifelse(!is.na(DP04_0100PE) & !is.na(DP04_0101PE), DP04_0100PE + DP04_0101PE, NA_real_) else NA_real_,
          pct_owners_no_mortgage_cost_burdened_30plus = if(all(c("DP04_0106PE", "DP04_0107PE") %in% names(.))) ifelse(!is.na(DP04_0106PE) & !is.na(DP04_0107PE), DP04_0106PE + DP04_0107PE, NA_real_) else NA_real_
        ) %>%
        mutate(
          num_burdened_renters = if(all(c("pct_renters_cost_burdened_30plus", "DP04_0126E") %in% names(.))) ifelse(!is.na(pct_renters_cost_burdened_30plus) & !is.na(DP04_0126E), (pct_renters_cost_burdened_30plus / 100) * DP04_0126E, NA_real_) else NA_real_,
          num_burdened_owners_w_mort = if(all(c("pct_owners_w_mortgage_cost_burdened_30plus", "DP04_0095E") %in% names(.))) ifelse(!is.na(pct_owners_w_mortgage_cost_burdened_30plus) & !is.na(DP04_0095E), (pct_owners_w_mortgage_cost_burdened_30plus / 100) * DP04_0095E, NA_real_) else NA_real_,
          num_burdened_owners_no_mort = if(all(c("pct_owners_no_mortgage_cost_burdened_30plus", "DP04_0102E") %in% names(.))) ifelse(!is.na(pct_owners_no_mortgage_cost_burdened_30plus) & !is.na(DP04_0102E), (pct_owners_no_mortgage_cost_burdened_30plus / 100) * DP04_0102E, NA_real_) else NA_real_
        ) %>%
        rowwise() %>%
        mutate( total_burdened_households = sum(c( if("num_burdened_renters" %in% names(.)) num_burdened_renters else NA, if("num_burdened_owners_w_mort" %in% names(.)) num_burdened_owners_w_mort else NA, if("num_burdened_owners_no_mort" %in% names(.)) num_burdened_owners_no_mort else NA ), na.rm = FALSE) ) %>%
        ungroup() %>%
        mutate( pct_all_households_cost_burdened_30plus = if("total_burdened_households" %in% names(.) & "DP04_0045E" %in% names(.)) ifelse( !is.na(total_burdened_households) & !is.na(DP04_0045E) & DP04_0045E > 0, (total_burdened_households / DP04_0045E) * 100, NA_real_ ) else NA_real_ ) %>%
        select(GEOID, any_of(c("pct_renters_cost_burdened_30plus", "pct_owners_w_mortgage_cost_burdened_30plus", "pct_owners_no_mortgage_cost_burdened_30plus", "pct_all_households_cost_burdened_30plus")))
      data_list_for_year[["housing_burden"]] <- dp04_burden_processed
    } else {cat("    DP04 burden vars processing failed.\n")}
  } else if (year_to_fetch != YEAR_2009_CALC) {cat("    DP04 table for burden unavailable/empty for a non-2009 year.\n")}
  
  # --- Assemble ---
  if (is.null(base_geographies_df) || nrow(base_geographies_df) == 0) {
    cat(paste0("    CRITICAL FINAL: No base_geographies_df for year ", year_to_fetch, ". Returning NULL.\n"))
    return(NULL)
  }
  final_df <- base_geographies_df
  for (df_name in names(data_list_for_year)) {
    current_processed_df <- data_list_for_year[[df_name]]
    if (!is.null(current_processed_df) && nrow(current_processed_df) > 0 && "GEOID" %in% names(current_processed_df)) {
      cols_to_join <- setdiff(names(current_processed_df), c("GEOID", if ("NAME" %in% names(final_df) && "NAME" %in% names(current_processed_df)) "NAME" else NULL))
      if(length(cols_to_join) > 0) {
        final_df <- final_df %>% left_join(current_processed_df %>% select(GEOID, all_of(cols_to_join)), by = "GEOID")
      }
    }
  }
  final_df <- final_df %>% mutate(year_acs_product = as.integer(year_to_fetch))
  cat(paste0("--- Processed Year: ", year_to_fetch, ". MSAs: ", nrow(final_df), ". Final Cols: ", ncol(final_df), " ---\n"))
  return(final_df)
}

# PART 4: DATA AGGREGATION & PANEL CREATION (Identical to V6.5.0)
# =============================================================================
cat("\nPART 4: Creating Full MSA Panel...\n")
cat(paste0("Target ACS Product Years for panel: ", paste(ALL_TARGET_YEARS, collapse = ", "), "\n"))
msa_panel_list_full <- lapply(ALL_TARGET_YEARS, get_msa_data_panel)
msa_panel_list_full <- msa_panel_list_full[!sapply(msa_panel_list_full, is.null)]
if (length(msa_panel_list_full) == 0) stop("FATAL: No data fetched for ANY year.")
msa_panel_full_raw <- bind_rows(msa_panel_list_full)
cat(paste0("\nFull raw panel created with ", nrow(msa_panel_full_raw), " MSA-year observations.\n"))
if(nrow(msa_panel_full_raw) > 0) {
  cat(paste0("ACS Product Years in panel: ", paste(sort(unique(msa_panel_full_raw$year_acs_product)), collapse = ", "), "\n"))
  cat(paste0("Columns in panel (first 10): ", paste(head(names(msa_panel_full_raw),10), collapse = ", "), "...\n"))
}

# PART 5: FINAL CLEANING, FEATURE ENGINEERING & SAVING (Identical to V6.5.0)
# =============================================================================
cat("\nPART 5: Final Cleaning, Feature Engineering & Saving...\n")
if (nrow(msa_panel_full_raw) == 0) {
  cat("Panel is empty. Skipping final cleaning.\n")
} else {
  msa_panel_final <- msa_panel_full_raw
  if("year_acs_product" %in% names(msa_panel_final)) {
    msa_panel_final <- msa_panel_final %>% rename(year = year_acs_product)
  }
  intended_numeric_cols <- c(
    "total_population", "children_0_4", "cwr", "birth_rate_15_50", "birth_rate_15_50_moe",
    "total_families", "families_w_child_lt18", "share_families_with_children_lt18", # From DP02
    "share_households_with_children_lt18", "share_households_with_children_lt6", # From B11012
    "median_household_income", "median_earnings_workers", "median_renter_household_income",
    "pct_bachelors_or_higher", "median_home_value", "median_gross_rent",
    "pct_renters_cost_burdened_30plus", "pct_owners_w_mortgage_cost_burdened_30plus",
    "pct_owners_no_mortgage_cost_burdened_30plus", "pct_all_households_cost_burdened_30plus"
  )
  actual_cols_to_convert <- intersect(intended_numeric_cols, names(msa_panel_final))
  if (length(actual_cols_to_convert) > 0) {
    msa_panel_final <- msa_panel_final %>%
      mutate(across(all_of(actual_cols_to_convert), ~as.numeric(as.character(.))))
  }
  if (all(c("median_home_value", "median_household_income") %in% names(msa_panel_final))) {
    msa_panel_final <- msa_panel_final %>%
      mutate(price_to_income_ratio = ifelse(!is.na(median_household_income) & median_household_income > 0 &
                                              !is.na(median_home_value) & median_home_value > 0,
                                            median_home_value / median_household_income, NA_real_))
  } else {
    cat("Warning: median_home_value or median_household_income missing. Cannot calculate price_to_income_ratio.\n")
    if (!"price_to_income_ratio" %in% names(msa_panel_final)) msa_panel_final$price_to_income_ratio <- NA_real_
  }
  if (all(c("median_gross_rent", "median_renter_household_income") %in% names(msa_panel_final))) {
    msa_panel_final <- msa_panel_final %>%
      mutate(median_rent_to_renter_income_ratio = ifelse(!is.na(median_renter_household_income) & median_renter_household_income > 0 &
                                                           !is.na(median_gross_rent) & median_gross_rent > 0,
                                                         (median_gross_rent * 12) / median_renter_household_income, NA_real_))
  } else {
    cat("Warning: median_gross_rent or median_renter_household_income missing. Cannot calculate median_rent_to_renter_income_ratio.\n")
    if (!"median_rent_to_renter_income_ratio" %in% names(msa_panel_final)) msa_panel_final$median_rent_to_renter_income_ratio <- NA_real_
  }
  output_main_dir <- "msa_panel_output_v6_5_1"
  dir.create(output_main_dir, showWarnings = FALSE)
  saveRDS(msa_panel_final, file.path(output_main_dir, "msa_panel_v6_5_1_full.rds"))
  write_csv(msa_panel_final, file.path(output_main_dir, "msa_panel_v6_5_1_full.csv"))
  cat(paste0("\nFinal MSA panel V6.5.1 saved to '", output_main_dir, "'\n"))
  
  cat("\n--- Summary of Key Variables (Latest year) ---\n")
  if ("year" %in% names(msa_panel_final) && nrow(msa_panel_final) > 0) {
    latest_year <- max(msa_panel_final$year, na.rm = TRUE)
    if(is.finite(latest_year)){
      summary_cols <- c("GEOID", "NAME", "year", "total_population", "birth_rate_15_50", "fertility_data_type", "cwr",
                        "share_families_with_children_lt18", # From DP02
                        "share_households_with_children_lt18", "share_households_with_children_lt6", # From B11012
                        "median_household_income", "median_renter_household_income", "pct_bachelors_or_higher",
                        "median_home_value", "price_to_income_ratio", "median_gross_rent", "median_rent_to_renter_income_ratio",
                        "pct_renters_cost_burdened_30plus", "pct_owners_w_mortgage_cost_burdened_30plus",
                        "pct_owners_no_mortgage_cost_burdened_30plus", "pct_all_households_cost_burdened_30plus")
      cols_present <- intersect(summary_cols, names(msa_panel_final))
      panel_latest <- msa_panel_final %>% filter(year == latest_year)
      if(length(cols_present) > 2 && nrow(panel_latest) > 0) {
        print(paste("Summary for latest year:", latest_year))
        print(summary(panel_latest %>% select(all_of(cols_present))))
      } else {cat("Not enough data/cols for latest year summary.\n")}
    } else {cat("Could not determine latest year.\n")}
  } else {cat("Panel empty or 'year' missing.\n")}
  
  cat("\n--- Check 2009 data specifically ---\n")
  if ("year" %in% names(msa_panel_final) && nrow(msa_panel_final) > 0 && YEAR_2009_CALC %in% msa_panel_final$year) {
    summary_2009_cols <- c("GEOID", "NAME", "year", "birth_rate_15_50", "fertility_data_type", "total_population", "cwr",
                           "share_families_with_children_lt18", # From DP02, likely NA for 2009
                           "share_households_with_children_lt18", "share_households_with_children_lt6", # From B11012
                           "median_household_income", "median_renter_household_income", "median_rent_to_renter_income_ratio",
                           "pct_bachelors_or_higher", "pct_renters_cost_burdened_30plus", "pct_all_households_cost_burdened_30plus" )
    cols_present_2009 <- intersect(summary_2009_cols, names(msa_panel_final))
    panel_2009 <- msa_panel_final %>% filter(year == YEAR_2009_CALC)
    if(length(cols_present_2009) > 2 && nrow(panel_2009) > 0){
      cat("Top MSAs by population for 2009:\n")
      df_2009_head <- panel_2009 %>% select(all_of(cols_present_2009))
      if ("total_population" %in% names(df_2009_head) && sum(!is.na(df_2009_head$total_population)) > 0) {
        df_2009_head <- df_2009_head %>% arrange(desc(total_population))
      }
      print(head(df_2009_head, 5))
      check_na_counts <- function(df, var_name) {
        if (var_name %in% names(df)) { cat(paste0("MSAs w/ 2009 ", var_name, ": ", df %>% filter(!is.na(.data[[var_name]])) %>% nrow(), "\n"))
        } else { cat(paste0("Variable ", var_name, " not present in 2009 data.\n")) }
      }
      check_na_counts(panel_2009, "birth_rate_15_50")
      check_na_counts(panel_2009, "share_households_with_children_lt18") # From B11012
      check_na_counts(panel_2009, "share_families_with_children_lt18") # From DP02
      check_na_counts(panel_2009, "median_rent_to_renter_income_ratio")
      check_na_counts(panel_2009, "pct_bachelors_or_higher")
      check_na_counts(panel_2009, "pct_all_households_cost_burdened_30plus")
    } else {cat("Not enough data/cols for 2009 summary.\n")}
  } else {cat("No 2009 data in panel.\n")}
}
cat("\n--- V6.5.1 MSA PANEL SCRIPT COMPLETE ---\n")