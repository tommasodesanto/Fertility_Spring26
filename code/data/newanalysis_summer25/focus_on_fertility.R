# =============================================================================
# SCRIPT TO FETCH FERTILITY RATES (ACS S1301) - MSA & COUNTY - CORRECTED VARIABLE
# Fetches birth rates per 1,000 women aged 15-50 for MSAs and Counties.
# Uses S1301_C04_001 for all ACS 5-year data years.
# =============================================================================

# PART 0: SETUP & LIBRARIES
# =============================================================================
cat("PART 0: Setting up libraries and environment...\n")
suppressPackageStartupMessages({
  if (!require(tidycensus)) install.packages("tidycensus"); library(tidycensus)
  if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
})

if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("Census API Key not found. Please set it using census_api_key('YOUR_KEY', install=TRUE) and restarting R.")
}
options(tigris_use_cache = TRUE)

# PART 1: CONFIGURATION
# =============================================================================
cat("PART 1: Configuration...\n")
# --- General Year Configuration ---
start_data_year_acs5 <- 2009 
latest_data_year_acs5 <- 2022 

# --- MSA Specific Configuration ---
years_to_fetch_msa <- start_data_year_acs5:latest_data_year_acs5
# Or test with a smaller set first for MSAs:
# years_to_fetch_msa <- c(2018, 2019, 2021, 2022)

# --- County Specific Configuration ---
# For fetching county data, fetching all US counties for all ~14 years can be very slow.
# To fetch for all years for counties (can be slow): county_years_to_fetch_list <- start_data_year_acs5:latest_data_year_acs5
# To fetch for the last 3 available years for counties (recommended for initial run):
county_years_to_fetch_list <- tail(start_data_year_acs5:latest_data_year_acs5, 3) 
# To fetch for specific years for counties:
# county_years_to_fetch_list <- c(2010, 2015, 2020, 2022)

# To fetch for specific states for counties (e.g., California "06" and Texas "48"):
# county_state_filter <- c("06", "48") 
county_state_filter <- NULL # Set to NULL to fetch for all US counties

# --- Corrected Fertility Rate Variable Code ---
# Based on inspection, S1301_C04_001 is "Rate per 1,000 women...Women 15 to 50 years"
# for both pre-2022 and 2022+ ACS data for table S1301.
FERTILITY_RATE_VARIABLE_CODE <- "S1301_C04_001"

# PART 2: DATA FETCHING FUNCTION
# =============================================================================
cat("PART 2: Defining data fetching function...\n")
get_fertility_s1301_generic <- function(data_year, geography_level, state_filter = NULL) {
  cat(paste0("  Fetching S1301 fertility data for year: ", data_year, 
             ", Geography: ", geography_level))
  if (!is.null(state_filter) && geography_level == "county") cat(paste0(", State(s): ", paste(state_filter, collapse=", ")))
  cat("...\n")
  cat(paste0("    Using S1301 variable code: ", FERTILITY_RATE_VARIABLE_CODE, "\n"))
  
  data_raw <- NULL
  tryCatch({
    data_raw <- get_acs(
      geography = geography_level,
      table = "S1301",
      year = data_year,
      survey = "acs5",
      state = if (geography_level == "county") state_filter else NULL, # Apply state filter only for counties
      cache_table = TRUE,
      show_call = FALSE 
    )
    
    if (is.null(data_raw) || nrow(data_raw) == 0) {
      cat(paste0("    No data returned by API for S1301, year ", data_year, ", geo ", geography_level, ".\n"))
      return(NULL)
    }
    
    if (!FERTILITY_RATE_VARIABLE_CODE %in% data_raw$variable) {
      cat(paste0("    CRITICAL: Variable '", FERTILITY_RATE_VARIABLE_CODE, 
                 "' NOT FOUND in data for year ", data_year, ", geo ", geography_level, ".\n"))
      return(NULL)
    }
    
    processed_data <- data_raw %>%
      filter(variable == FERTILITY_RATE_VARIABLE_CODE) %>%
      select(GEOID, NAME, estimate, moe) %>%
      rename(
        birth_rate_15_50 = estimate,
        birth_rate_15_50_moe = moe
      ) %>%
      mutate(
        year_acs_product = as.integer(data_year),
        birth_rate_15_50 = as.numeric(birth_rate_15_50),
        birth_rate_15_50_moe = as.numeric(birth_rate_15_50_moe),
        data_period = paste0(data_year - 4, "-", data_year)
      )
    
    cat(paste0("    Successfully processed ", nrow(processed_data), " records for '", 
               FERTILITY_RATE_VARIABLE_CODE, "' (Geo: ", geography_level, ").\n"))
    return(processed_data)
    
  }, error = function(e) {
    cat(paste0("    Error fetching/processing data for year ", data_year, ", geo ", geography_level, ": ", e$message, "\n"))
    return(NULL)
  })
}

# PART 3: FETCHING DATA FOR MSAs AND COUNTIES
# =============================================================================

# --- Fetch MSA Data ---
cat("\nPART 3A: Fetching MSA data...\n")
cat("Fetching MSA data for ACS Product years:", paste(years_to_fetch_msa, collapse=", "), "\n")
msa_fertility_list <- list()
for (yr_iter in years_to_fetch_msa) {
  Sys.sleep(0.3) 
  msa_data_for_year <- get_fertility_s1301_generic(data_year = yr_iter, 
                                                   geography_level = "metropolitan statistical area/micropolitan statistical area")
  if (!is.null(msa_data_for_year) && nrow(msa_data_for_year) > 0) {
    msa_fertility_list[[as.character(yr_iter)]] <- msa_data_for_year
  }
}
msa_fertility_panel <- bind_rows(msa_fertility_list)


# --- Fetch County Data ---
cat("\nPART 3B: Fetching County data...\n")
cat("Fetching County data for ACS Product years:", paste(county_years_to_fetch_list, collapse=", "), "\n")
if (!is.null(county_state_filter)) {
  cat(paste0("For states: ", paste(county_state_filter, collapse=", "), "\n"))
} else {
  cat("For all US states.\n")
}
cat("Note: Fetching all US counties for many years can be time-consuming.\n")

county_fertility_list <- list()
for (yr_iter in county_years_to_fetch_list) {
  Sys.sleep(0.3) 
  county_data_for_year <- get_fertility_s1301_generic(data_year = yr_iter, 
                                                      geography_level = "county",
                                                      state_filter = county_state_filter)
  if (!is.null(county_data_for_year) && nrow(county_data_for_year) > 0) {
    county_fertility_list[[as.character(yr_iter)]] <- county_data_for_year
  }
}
county_fertility_panel <- bind_rows(county_fertility_list)


# PART 4: SUMMARIZE AND OUTPUT
# =============================================================================
cat("\n\nPART 4: Summarizing and Outputting Results...\n")

output_dir <- "fertility_data_s1301_c04_001_output"
dir.create(output_dir, showWarnings = FALSE)

# --- MSA Fertility Data Summary ---
cat("\n--- MSA Fertility Data ---\n")
if (!is.null(msa_fertility_panel) && nrow(msa_fertility_panel) > 0) {
  print(paste("Total MSA-year observations:", nrow(msa_fertility_panel)))
  print(paste("ACS Product Years successfully fetched for MSAs:", paste(sort(unique(msa_fertility_panel$year_acs_product)), collapse = ", ")))
  
  cat("\nSample data (MSA):\n")
  print(head(msa_fertility_panel))
  
  msa_na_counts <- msa_fertility_panel %>% group_by(year_acs_product) %>% summarise(total_rows = n(), na_in_birth_rate = sum(is.na(birth_rate_15_50)), .groups = "drop")
  cat("\nNA counts for MSA birth_rate_15_50 by year:\n"); print(knitr::kable(msa_na_counts))
  
  cat("\nSummary of MSA birth rates by ACS Product Year (excluding NAs):\n")
  msa_summary_by_year <- msa_fertility_panel %>%
    filter(!is.na(birth_rate_15_50)) %>% 
    group_by(year_acs_product, data_period) %>%
    summarise(
      n_geographies_with_data = n(),
      avg_birth_rate = round(mean(birth_rate_15_50), 1), median_birth_rate = round(median(birth_rate_15_50), 1),
      min_birth_rate = round(min(birth_rate_15_50), 1), max_birth_rate = round(max(birth_rate_15_50), 1),
      .groups = "drop"
    ) %>% arrange(year_acs_product)
  print(knitr::kable(msa_summary_by_year, caption = "MSA Fertility Rate Summary (S1301_C04_001)"))
  
  msa_file_path <- file.path(output_dir, paste0("msa_fertility_rates_S1301_C04_001_", min(years_to_fetch_msa), "_", max(years_to_fetch_msa), ".csv"))
  write_csv(msa_fertility_panel, msa_file_path)
  cat(paste0("\nMSA fertility data saved to: ", msa_file_path, "\n"))
} else {
  cat("No MSA fertility data was successfully compiled.\n")
}

# --- County Fertility Data Summary ---
cat("\n\n--- County Fertility Data ---\n")
if (!is.null(county_fertility_panel) && nrow(county_fertility_panel) > 0) {
  print(paste("Total county-year observations:", nrow(county_fertility_panel)))
  print(paste("ACS Product Years successfully fetched for Counties:", paste(sort(unique(county_fertility_panel$year_acs_product)), collapse = ", ")))
  
  cat("\nSample data (County):\n")
  print(head(county_fertility_panel))
  
  county_na_counts <- county_fertility_panel %>% group_by(year_acs_product) %>% summarise(total_rows = n(), na_in_birth_rate = sum(is.na(birth_rate_15_50)), .groups = "drop")
  cat("\nNA counts for County birth_rate_15_50 by year:\n"); print(knitr::kable(county_na_counts))
  
  cat("\nSummary of County birth rates by ACS Product Year (excluding NAs):\n")
  county_summary_by_year <- county_fertility_panel %>%
    filter(!is.na(birth_rate_15_50)) %>%
    group_by(year_acs_product, data_period) %>%
    summarise(
      n_geographies_with_data = n(),
      avg_birth_rate = round(mean(birth_rate_15_50), 1), median_birth_rate = round(median(birth_rate_15_50), 1),
      min_birth_rate = round(min(birth_rate_15_50), 1), max_birth_rate = round(max(birth_rate_15_50), 1),
      .groups = "drop"
    ) %>% arrange(year_acs_product)
  print(knitr::kable(county_summary_by_year, caption = "County Fertility Rate Summary (S1301_C04_001)"))
  
  county_file_name_suffix <- if(!is.null(county_state_filter)) paste0("_states_", paste(county_state_filter, collapse="-")) else "_allUS"
  county_file_path <- file.path(output_dir, paste0("county_fertility_rates_S1301_C04_001_", 
                                                   min(county_years_to_fetch_list), "_", 
                                                   max(county_years_to_fetch_list),
                                                   county_file_name_suffix, ".csv"))
  write_csv(county_fertility_panel, county_file_path)
  cat(paste0("\nCounty fertility data saved to: ", county_file_path, "\n"))
} else {
  cat("No county fertility data was successfully compiled.\n")
}

cat("\n\n--- SCRIPT COMPLETE (MSA & COUNTY) ---\n")
cat(paste0("This script uses '", FERTILITY_RATE_VARIABLE_CODE, "' for all years for both geographies.\n"))