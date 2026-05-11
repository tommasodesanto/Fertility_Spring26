# =============================================================================
# SCRIPT TO FETCH MSA FERTILITY RATES (ACS 5-YEAR, 2009 PRODUCT CALCULATED + LATER S1301)
# - ACS 2009 (product year): Calculated from B01001 (denominator) and B13002 (numerator)
# - ACS 2010 onwards: Uses S1301_C04_001
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

# --- ACS 2009 (Product Year) Configuration for Calculated Rate ---
year_2009_calc <- 2009
vars_b01001_denominator_2009 <- c(
  f_15_17 = "B01001_030E", f_18_19 = "B01001_031E", f_20 = "B01001_032E",
  f_21 = "B01001_033E", f_22_24 = "B01001_034E", f_25_29 = "B01001_035E",
  f_30_34 = "B01001_036E", f_35_39 = "B01001_037E", f_40_44 = "B01001_038E",
  f_45_49 = "B01001_039E" 
  # Summing these gives women 15-49 years
)
var_b13002_numerator_2009 <- c(women_with_birth_15_50 = "B13002_002E") 
# B13002_002E: Total!!Estimate!!WOMEN 15 TO 50 YEARS WHO HAD A BIRTH IN THE PAST 12 MONTHS

# --- ACS 2010 Onwards Configuration (S1301) ---
start_data_year_acs_s1301 <- 2010 
latest_data_year_acs_s1301 <- 2022 
years_to_fetch_msa_s1301 <- start_data_year_acs_s1301:latest_data_year_acs_s1301
acs_survey_s1301 <- "acs5" # For S1301, it implies "acs5/subject"
S1301_FERTILITY_RATE_VARIABLE_CODE <- "S1301_C04_001" 

# PART 2: DATA FETCHING FUNCTIONS
# =============================================================================
cat("PART 2: Defining data fetching functions...\n")

# Function for ACS 2009 Calculated Rate
fetch_msa_fertility_acs2009_calculated <- function() {
  cat(paste0("  Fetching data for calculated ACS 2009 (product year) fertility rate (MSA)...\n"))
  cat(paste0("    Denominator from B01001 (Women 15-49), Numerator from B13002_002E (Women 15-50 with a birth).\n"))
  
  all_vars_2009 <- c(vars_b01001_denominator_2009, var_b13002_numerator_2009)
  
  data_raw_2009 <- NULL
  tryCatch({
    data_raw_2009 <- get_acs(
      geography = "metropolitan statistical area/micropolitan statistical area",
      variables = all_vars_2009, 
      year = year_2009_calc, 
      survey = "acs5", # Detailed tables endpoint
      output = "wide",
      cache_table = TRUE, show_call = FALSE
    )
    
    if (is.null(data_raw_2009) || nrow(data_raw_2009) == 0) {
      cat("    No data returned for B01001/B13002 for ACS 2009 detailed tables.\n"); return(NULL)
    }
    
    # Check if all expected columns are present
    # The 'E' is already in the variable names defined in all_vars_2009
    expected_cols_2009 <- names(all_vars_2009) 
    if (!all(expected_cols_2009 %in% names(data_raw_2009))) {
      missing_cols <- setdiff(expected_cols_2009, names(data_raw_2009))
      cat(paste0("    Not all expected B01001/B13002 variables found for 2009. Missing: ", 
                 paste(missing_cols, collapse=", "), "\n"))
      return(NULL)
    }
    
    processed_data_2009 <- data_raw_2009 %>%
      mutate(
        # Calculate denominator: Sum of women aged 15-49
        denominator_women_15_49 = rowSums(select(., all_of(names(vars_b01001_denominator_2009))), na.rm = FALSE), # Keep NA if any component is NA
        # Numerator: Women 15-50 with a birth
        numerator_women_with_birth = !!sym(names(var_b13002_numerator_2009)[1])
      ) %>%
      filter(!is.na(denominator_women_15_49) & denominator_women_15_49 > 0 & 
               !is.na(numerator_women_with_birth) & numerator_women_with_birth >= 0) %>%
      mutate(
        birth_rate_15_50 = (numerator_women_with_birth / denominator_women_15_49) * 1000,
        birth_rate_15_50_moe = NA_real_, # MOE for calculated rate is complex and omitted
        year = as.integer(year_2009_calc),
        data_period_covered = paste0(year_2009_calc - 4, "-", year_2009_calc),
        data_period_type = "ACS 5-year (Calculated from B-tables)",
        source_tables = "B01001 & B13002"
      ) %>%
      select(GEOID, NAME, year, birth_rate_15_50, birth_rate_15_50_moe, data_period_covered, data_period_type, source_tables)
    
    cat(paste0("    Successfully calculated ACS 2009 rates for ", nrow(processed_data_2009), " MSAs.\n"))
    return(processed_data_2009)
    
  }, error = function(e) {
    cat(paste0("    Error in fetch_msa_fertility_acs2009_calculated: ", e$message, "\n")); return(NULL)
  })
}

# Function for ACS 2010 Onwards (S1301)
fetch_msa_fertility_s1301 <- function(data_year) {
  cat(paste0("  Fetching S1301 (", acs_survey_s1301, ") for product year: ", data_year, " (MSA)...\n"))
  # ... (rest of this function is same as fetch_msa_fertility_acs5 from previous correct script) ...
  cat(paste0("    This covers data period: ", data_year - 4, "-", data_year, "\n"))
  cat(paste0("    Using variable code: ", S1301_FERTILITY_RATE_VARIABLE_CODE, "\n"))
  data_raw <- NULL
  tryCatch({
    data_raw <- get_acs(
      geography = "metropolitan statistical area/micropolitan statistical area",
      table = "S1301", year = data_year, survey = acs_survey_s1301, 
      cache_table = TRUE, show_call = FALSE
    )
    if (is.null(data_raw) || nrow(data_raw) == 0) {cat(paste0("    No data from S1301 (", acs_survey_s1301, "), product year ", data_year, ".\n")); return(NULL)}
    if (!S1301_FERTILITY_RATE_VARIABLE_CODE %in% data_raw$variable) {
      cat(paste0("    CRITICAL: Var '", S1301_FERTILITY_RATE_VARIABLE_CODE, 
                 "' NOT FOUND in S1301 for product year ", data_year, ".\n")); return(NULL)
    }
    processed_data <- data_raw %>%
      filter(variable == S1301_FERTILITY_RATE_VARIABLE_CODE) %>%
      select(GEOID, NAME, estimate, moe) %>%
      rename(birth_rate_15_50 = estimate, birth_rate_15_50_moe = moe) %>%
      mutate(
        year = as.integer(data_year),
        birth_rate_15_50 = as.numeric(birth_rate_15_50),
        birth_rate_15_50_moe = as.numeric(birth_rate_15_50_moe),
        data_period_covered = paste0(data_year - 4, "-", data_year), 
        data_period_type = "ACS 5-year (S1301)",
        source_tables = "S1301"
      ) %>%
      select(GEOID, NAME, year, birth_rate_15_50, birth_rate_15_50_moe, data_period_covered, data_period_type, source_tables)
    cat(paste0("    Successfully processed ", nrow(processed_data), " MSA records from S1301 for product year ", data_year, ".\n"))
    return(processed_data)
  }, error = function(e) {cat(paste0("    Error fetching/processing S1301 data for year ", data_year, ": ", e$message, "\n")); return(NULL)})
}


# PART 3: FETCHING AND COMBINING DATA
# =============================================================================
cat("\nPART 3: Fetching and combining data...\n")

# --- Fetch Calculated ACS 2009 Data ---
acs2009_calculated_fertility_msa <- fetch_msa_fertility_acs2009_calculated()

# --- Fetch ACS S1301 Data for 2010 Onwards ---
cat("\nFetching MSA data from S1301 for ACS Product years:", paste(years_to_fetch_msa_s1301, collapse=", "), "\n")
s1301_fertility_list_msa <- list()
for (yr_iter in years_to_fetch_msa_s1301) {
  Sys.sleep(0.3) 
  msa_data_for_year <- fetch_msa_fertility_s1301(data_year = yr_iter)
  if (!is.null(msa_data_for_year) && nrow(msa_data_for_year) > 0) {
    s1301_fertility_list_msa[[as.character(yr_iter)]] <- msa_data_for_year
  } else {
    cat(paste0("    -> No S1301 data successfully processed for MSAs for product year ", yr_iter, ".\n"))
  }
}
s1301_fertility_panel_msa <- bind_rows(s1301_fertility_list_msa)

# --- Combine Data ---
cat("\nCombining Calculated ACS 2009 and S1301 ACS data...\n")
all_fertility_components <- list()
if (!is.null(acs2009_calculated_fertility_msa) && nrow(acs2009_calculated_fertility_msa) > 0) {
  all_fertility_components[["acs2009_calc"]] <- acs2009_calculated_fertility_msa
} else {
  warning("Calculated ACS 2009 data not fetched successfully.")
}
if (!is.null(s1301_fertility_panel_msa) && nrow(s1301_fertility_panel_msa) > 0) {
  all_fertility_components[["acs_s1301"]] <- s1301_fertility_panel_msa
} else {
  warning("ACS S1301 data (2010 onwards) not fetched successfully.")
}

if (length(all_fertility_components) > 0) {
  combined_fertility_panel_msa <- bind_rows(all_fertility_components) %>%
    arrange(GEOID, year)
  cat(paste0("Combined panel has ", nrow(combined_fertility_panel_msa), " MSA-year observations.\n"))
} else {
  combined_fertility_panel_msa <- tibble() # Empty tibble if nothing fetched
  warning("No fertility data fetched successfully from any source.")
}


# PART 4: SUMMARIZE AND OUTPUT
# =============================================================================
cat("\n\nPART 4: Summarizing and Outputting Results (Combined Panel)...\n")
output_dir_combined <- "fertility_data_msa_acs_2009calc_onwards_output"
dir.create(output_dir_combined, showWarnings = FALSE)

cat("\n--- Combined MSA Fertility Data (ACS 2009 Calc + ACS 2010 Onwards S1301) ---\n")
if (!is.null(combined_fertility_panel_msa) && nrow(combined_fertility_panel_msa) > 0) {
  print(paste("Total MSA-year observations (Combined Panel):", nrow(combined_fertility_panel_msa)))
  print(paste("Years in combined panel (product years):", paste(sort(unique(combined_fertility_panel_msa$year)), collapse = ", ")))
  
  cat("\nSample data (Combined Panel - check 2009 carefully):\n")
  print(head(combined_fertility_panel_msa %>% filter(year == 2009)))
  print(head(combined_fertility_panel_msa %>% filter(year == 2010)))
  
  combined_na_counts <- combined_fertility_panel_msa %>% 
    group_by(year, data_period_type) %>% 
    summarise(total_msas_for_year = n(), 
              na_in_birth_rate = sum(is.na(birth_rate_15_50)),
              .groups = "drop") %>%
    arrange(year)
  cat("\nNA counts for MSA birth_rate_15_50 (Combined) by product year:\n")
  print(knitr::kable(combined_na_counts))
  
  cat("\nSummary of MSA birth rates by ACS Product Year (Combined, excluding NAs):\n")
  combined_summary_by_year <- combined_fertility_panel_msa %>%
    filter(!is.na(birth_rate_15_50)) %>% 
    group_by(year, data_period_covered, data_period_type, source_tables) %>%
    summarise(
      n_msas_with_data = n(),
      avg_birth_rate = round(mean(birth_rate_15_50), 1), 
      median_birth_rate = round(median(birth_rate_15_50), 1),
      min_birth_rate = round(min(birth_rate_15_50), 1), 
      max_birth_rate = round(max(birth_rate_15_50), 1),
      .groups = "drop"
    ) %>% arrange(year)
  print(knitr::kable(combined_summary_by_year, caption = "MSA Fertility Rate Summary (Combined)"))
  
  combined_file_path <- file.path(output_dir_combined, 
                                  paste0("msa_fertility_rates_ACS_2009calc_onwards_", 
                                         min(combined_fertility_panel_msa$year), "_", 
                                         max(combined_fertility_panel_msa$year), ".csv"))
  write_csv(combined_fertility_panel_msa, combined_file_path)
  cat(paste0("\nCombined fertility data saved to: ", combined_file_path, "\n"))
  
  if (nrow(combined_summary_by_year) > 1) {
    p_trend_combined <- ggplot(combined_summary_by_year, 
                               aes(x = year, y = avg_birth_rate, 
                                   color = data_period_type, group = 1)) +
      geom_line(linewidth = 1) +
      geom_point(aes(shape = data_period_type), size = 3) +
      geom_text(aes(label = n_msas_with_data), vjust = -0.8, size = 3, show.legend = FALSE) +
      scale_x_continuous(breaks = combined_summary_by_year$year) +
      scale_color_manual(values = c("ACS 5-year (Calculated from B-tables)" = "darkorange", "ACS 5-year (S1301)" = "darkgreen")) +
      scale_shape_manual(values = c("ACS 5-year (Calculated from B-tables)" = 15, "ACS 5-year (S1301)" = 17)) +
      labs(title = "Average MSA Fertility Rate Trend (Calculated 2009 + S1301 2010-2022)",
           subtitle = "ACS 5-Year plotted at product end-year of period. Numbers on points indicate count of MSAs with data.\nNote: 2009 rate is 'women with a birth per 1k women 15-49'; 2010+ is 'births per 1k women 15-50'.",
           x = "ACS 5-Year Product End Year", 
           y = "Average Rate per 1,000 Women",
           color = "Data Type/Source", shape = "Data Type/Source") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    print(p_trend_combined)
    ggsave(file.path(output_dir_combined, "msa_avg_fertility_trend_2009calc_onwards.png"), plot = p_trend_combined, width = 12, height = 7)
    cat("Combined trend plot saved.\n")
  }
  
} else {
  cat("No combined MSA fertility data was successfully compiled.\n")
}

cat("\n\n--- SCRIPT COMPLETE (MSA ACS 2009 Calc + 2010 Onwards S1301) ---\n")
cat("NOTE: The 2009 ACS 5-year rate is calculated as (Women 15-50 with a birth / Women 15-49) * 1000.\n")
cat("This is conceptually slightly different from the S1301 rate (Births per 1000 women 15-50) used for 2010 onwards.\n")
cat("MOEs for the calculated 2009 rate are NOT provided.\n")