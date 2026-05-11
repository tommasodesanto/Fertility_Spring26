# =============================================================================
# MSA PANEL DATA COLLECTION: ACS-1 & DECENNIAL CENSUS
# Version 1.1 - Comprehensive fertility and household analysis with robust error handling
# =============================================================================

# SETUP & CONFIGURATION
# =============================================================================
cat("Setting up environment...\n")
suppressPackageStartupMessages({
  required_packages <- c("tidycensus", "tidyverse", "purrr", "furrr", "progressr")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
})

# Enable parallel processing for efficiency
plan(multisession, workers = availableCores() - 1)

# Check Census API key
if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("Census API Key not found. Please set using census_api_key('YOUR_KEY', install=TRUE)")
}

options(tigris_use_cache = TRUE)

# CONFIGURATION
# =============================================================================
ACS_START_YEAR <- 2005  # First year of ACS-1
ACS_END_YEAR <- 2023    # Latest available
DECENNIAL_YEARS <- c(2000, 2010, 2020)
GEOGRAPHY <- "metropolitan statistical area/micropolitan statistical area"

# CONFIGURATION
# =============================================================================
# ... (other configurations like ACS_START_YEAR, etc.)

# BASE PROJECT DIRECTORY (USER: Set this to your main R project folder)
BASE_PROJECT_DIR <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/R"

# OUTPUT CONFIGURATION
output_dir <- file.path(BASE_PROJECT_DIR, "msa_panel_output") # Outputs will go here
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Data generation outputs will be saved to: %s\n", output_dir))


# Create log file
log_file <- file.path(output_dir, paste0("data_collection_log_", Sys.Date(), ".txt"))
# ... (rest of the script, including log_message function, etc.)
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
}

log_message("=== MSA PANEL DATA COLLECTION STARTED ===")
log_message(sprintf("ACS Years: %d-%d", ACS_START_YEAR, ACS_END_YEAR))
log_message(sprintf("Decennial Years: %s", paste(DECENNIAL_YEARS, collapse = ", ")))

# HELPER FUNCTIONS
# =============================================================================

# Safe data fetching wrapper with retry logic
safe_get_data <- function(fetch_function, ..., max_retries = 3) {
  for (i in 1:max_retries) {
    result <- tryCatch({
      fetch_function(...)
    }, error = function(e) {
      if (i < max_retries) {
        Sys.sleep(2^i)  # Exponential backoff
        NULL
      } else {
        log_message(paste("Error after", max_retries, "attempts:", e$message))
        NULL
      }
    })
    if (!is.null(result)) return(result)
  }
  return(NULL)
}

# Function to check data quality
check_data_quality <- function(data, year, source) {
  if (is.null(data) || nrow(data) == 0) {
    log_message(sprintf("WARNING: No data for %s year %d", source, year))
    return(FALSE)
  }
  
  # Check for critical variables
  critical_vars <- c("GEOID", "NAME")
  missing_vars <- setdiff(critical_vars, names(data))
  if (length(missing_vars) > 0) {
    log_message(sprintf("WARNING: Missing critical variables in %s year %d: %s", 
                        source, year, paste(missing_vars, collapse = ", ")))
    return(FALSE)
  }
  
  # Check for duplicate GEOIDs
  if (any(duplicated(data$GEOID))) {
    log_message(sprintf("WARNING: Duplicate GEOIDs in %s year %d", source, year))
    data <- data[!duplicated(data$GEOID), ]
  }
  
  return(TRUE)
}

# ACS DATA COLLECTION FUNCTION - IMPROVED VERSION
# =============================================================================
get_acs_data <- function(year, use_acs1 = TRUE) {
  survey_type <- if (use_acs1) "acs1" else "acs5"
  log_message(sprintf("Fetching %s data for year %d...", toupper(survey_type), year))
  
  # Initialize results list
  results <- list()
  
  # 1. POPULATION AND AGE STRUCTURE - FETCH WHOLE TABLE
  # -------------------------------------------------------------------------
  log_message("  Fetching population and age data...")
  
  # Get total population first
  total_pop <- safe_get_data(get_acs,
                             geography = GEOGRAPHY,
                             variables = "B01003_001",
                             year = year,
                             survey = survey_type)
  
  if (is.null(total_pop)) {
    log_message(sprintf("  No data available for %s year %d", survey_type, year))
    return(NULL)
  }
  
  # Try to get the entire B01001 table
  age_sex_data <- safe_get_data(get_acs,
                                geography = GEOGRAPHY,
                                table = "B01001",
                                year = year,
                                survey = survey_type,
                                cache_table = TRUE)
  
  if (!is.null(age_sex_data)) {
    # Process age/sex data
    pop_data <- age_sex_data %>%
      select(GEOID, NAME, variable, estimate) %>%
      filter(variable %in% c(
        "B01001_003", "B01001_027",  # Under 5
        "B01001_004", "B01001_028",  # 5-9
        "B01001_005", "B01001_029",  # 10-14
        "B01001_006", "B01001_030",  # 15-17
        "B01001_007", "B01001_031",  # 18-19
        "B01001_032", "B01001_033", "B01001_034", # Women 20-24
        "B01001_035", "B01001_036", "B01001_037", # Women 25-39
        "B01001_038", "B01001_039"   # Women 40-49
      )) %>%
      pivot_wider(names_from = variable, values_from = estimate)
    
    # Add total population
    pop_data <- pop_data %>%
      left_join(
        total_pop %>% 
          filter(variable == "B01003_001") %>%
          select(GEOID, estimate) %>%
          rename(total_population = estimate),
        by = "GEOID"
      )
    
    # Calculate derived variables with available columns
    available_cols <- names(pop_data)
    
    # Children calculations
    if (all(c("B01001_003", "B01001_027") %in% available_cols)) {
      pop_data$children_under_5 <- pop_data$B01001_003 + pop_data$B01001_027
    }
    
    if (all(c("B01001_004", "B01001_028") %in% available_cols)) {
      pop_data$children_5_to_9 <- pop_data$B01001_004 + pop_data$B01001_028
    }
    
    if (all(c("B01001_005", "B01001_029") %in% available_cols)) {
      pop_data$children_10_to_14 <- pop_data$B01001_005 + pop_data$B01001_029
    }
    
    if (all(c("B01001_006", "B01001_030") %in% available_cols)) {
      pop_data$children_15_to_17 <- pop_data$B01001_006 + pop_data$B01001_030
    }
    
    # Total children
    child_cols <- c("children_under_5", "children_5_to_9", "children_10_to_14", "children_15_to_17")
    child_cols_exist <- child_cols[child_cols %in% names(pop_data)]
    if (length(child_cols_exist) > 0) {
      pop_data$children_under_18 <- rowSums(pop_data[child_cols_exist], na.rm = FALSE)
    }
    
    # Women 15-49 calculations
    women_vars <- c("B01001_030", "B01001_031", "B01001_032", "B01001_033", 
                    "B01001_034", "B01001_035", "B01001_036", "B01001_037", 
                    "B01001_038", "B01001_039")
    women_vars_exist <- women_vars[women_vars %in% available_cols]
    
    if (length(women_vars_exist) > 0) {
      pop_data$women_15_49 <- rowSums(pop_data[women_vars_exist], na.rm = FALSE)
    }
    
    # Calculate ratios
    if (all(c("children_under_18", "total_population") %in% names(pop_data))) {
      pop_data$pct_under_18 <- (pop_data$children_under_18 / pop_data$total_population) * 100
    }
    
    if (all(c("children_under_5", "women_15_49") %in% names(pop_data))) {
      pop_data$child_woman_ratio <- ifelse(
        pop_data$women_15_49 > 0,
        (pop_data$children_under_5 / pop_data$women_15_49) * 1000,
        NA_real_
      )
    }
    
    # Select final columns
    keep_cols <- c("GEOID", "NAME", "total_population", "children_under_18", 
                   "pct_under_18", "children_under_5", "women_15_49", "child_woman_ratio")
    keep_cols <- keep_cols[keep_cols %in% names(pop_data)]
    
    results[["population"]] <- pop_data %>% select(all_of(keep_cols))
  }
  
  # 2. HOUSEHOLDS WITH CHILDREN (NOT families!)
  # -------------------------------------------------------------------------
  log_message("  Fetching household data...")
  
  # Try B11005 first (preferred)
  household_data <- safe_get_data(get_acs,
                                  geography = GEOGRAPHY,
                                  variables = c("B11005_001", "B11005_002"),
                                  year = year,
                                  survey = survey_type)
  
  if (!is.null(household_data)) {
    household_summary <- household_data %>%
      select(GEOID, variable, estimate) %>%
      pivot_wider(names_from = variable, values_from = estimate) %>%
      rename(
        total_households = B11005_001,
        households_with_children = B11005_002
      ) %>%
      mutate(
        pct_households_with_children = (households_with_children / total_households) * 100
      )
    
    results[["households"]] <- household_summary
  }
  
  # 3. FERTILITY DATA - B13002 Method
  # -------------------------------------------------------------------------
  log_message("  Fetching fertility data...")
  
  fertility_data <- safe_get_data(get_acs,
                                  geography = GEOGRAPHY,
                                  variables = c("B13002_001", "B13002_002"),
                                  year = year,
                                  survey = survey_type)
  
  if (!is.null(fertility_data)) {
    fertility_summary <- fertility_data %>%
      select(GEOID, variable, estimate) %>%
      pivot_wider(names_from = variable, values_from = estimate) %>%
      rename(
        women_15_50 = B13002_001,
        women_gave_birth = B13002_002
      ) %>%
      mutate(
        birth_rate_per_1000 = ifelse(
          women_15_50 > 0 & !is.na(women_gave_birth),
          (women_gave_birth / women_15_50) * 1000,
          NA_real_
        ),
        fertility_data_source = "B13002"
      )
    
    results[["fertility"]] <- fertility_summary
  }
  
  # 4. S1301 FERTILITY (if available and year >= 2010)
  # -------------------------------------------------------------------------
  if (year >= 2010 && survey_type == "acs5") {
    log_message("  Checking for S1301 fertility data...")
    
    s1301_data <- safe_get_data(get_acs,
                                geography = GEOGRAPHY,
                                table = "S1301",
                                year = year,
                                survey = survey_type,
                                cache_table = TRUE)
    
    if (!is.null(s1301_data)) {
      s1301_rate <- s1301_data %>%
        filter(variable == "S1301_C03_001") %>%
        select(GEOID, estimate) %>%
        rename(birth_rate_s1301 = estimate) %>%
        mutate(s1301_available = TRUE)
      
      results[["fertility_s1301"]] <- s1301_rate
    }
  }
  
  # 5. ECONOMIC INDICATORS
  # -------------------------------------------------------------------------
  log_message("  Fetching economic data...")
  
  econ_data <- safe_get_data(get_acs,
                             geography = GEOGRAPHY,
                             variables = c("B19013_001", "B25077_001", "B25064_001"),
                             year = year,
                             survey = survey_type)
  
  if (!is.null(econ_data)) {
    econ_summary <- econ_data %>%
      select(GEOID, variable, estimate) %>%
      pivot_wider(names_from = variable, values_from = estimate) %>%
      rename(
        median_hh_income = B19013_001,
        median_home_value = B25077_001,
        median_gross_rent = B25064_001
      ) %>%
      mutate(
        price_to_income = median_home_value / median_hh_income,
        annual_rent_to_income = (median_gross_rent * 12) / median_hh_income
      )
    
    results[["economics"]] <- econ_summary
  }
  
  # 6. EDUCATION
  # -------------------------------------------------------------------------
  log_message("  Fetching education data...")
  
  # Try to get the full B15003 table
  edu_data <- safe_get_data(get_acs,
                            geography = GEOGRAPHY,
                            table = "B15003",
                            year = year,
                            survey = survey_type,
                            cache_table = TRUE)
  
  if (!is.null(edu_data)) {
    edu_summary <- edu_data %>%
      filter(variable %in% c("B15003_001", "B15003_022", "B15003_023", 
                             "B15003_024", "B15003_025")) %>%
      select(GEOID, variable, estimate) %>%
      pivot_wider(names_from = variable, values_from = estimate)
    
    # Calculate if we have the necessary columns
    if (all(c("B15003_001", "B15003_022") %in% names(edu_summary))) {
      edu_cols <- c("B15003_022", "B15003_023", "B15003_024", "B15003_025")
      edu_cols_exist <- edu_cols[edu_cols %in% names(edu_summary)]
      
      edu_summary$college_plus <- rowSums(edu_summary[edu_cols_exist], na.rm = FALSE)
      edu_summary$pct_college_plus <- (edu_summary$college_plus / edu_summary$B15003_001) * 100
      
      results[["education"]] <- edu_summary %>%
        select(GEOID, pct_college_plus)
    }
  }
  
  # COMBINE ALL COMPONENTS
  # -------------------------------------------------------------------------
  if (length(results) == 0) {
    log_message(sprintf("WARNING: No data retrieved for %s year %d", survey_type, year))
    return(NULL)
  }
  
  # Join all data components
  combined_data <- results[[1]]
  for (i in 2:length(results)) {
    if (!is.null(results[[i]])) {
      # Remove NAME column from subsequent joins to avoid duplicates
      join_cols <- setdiff(names(results[[i]]), "NAME")
      combined_data <- combined_data %>%
        left_join(results[[i]][join_cols], by = "GEOID")
    }
  }
  
  # Add metadata
  combined_data <- combined_data %>%
    mutate(
      year = year,
      data_source = toupper(survey_type),
      collection_date = Sys.Date()
    )
  
  # Ensure NAME column exists
  if (!"NAME" %in% names(combined_data) && "NAME" %in% names(results[[1]])) {
    combined_data <- combined_data %>%
      left_join(results[[1]] %>% select(GEOID, NAME), by = "GEOID")
  }
  
  # Reorder columns
  col_order <- c("GEOID", "NAME", "year", "data_source")
  other_cols <- setdiff(names(combined_data), col_order)
  combined_data <- combined_data %>%
    select(all_of(c(col_order, other_cols)))
  
  log_message(sprintf("  Successfully retrieved %d MSAs for %s year %d", 
                      nrow(combined_data), survey_type, year))
  
  return(combined_data)
}

# DECENNIAL CENSUS DATA COLLECTION
# =============================================================================
get_decennial_data <- function(year) {
  log_message(sprintf("Fetching Decennial Census data for year %d...", year))
  
  results <- list()
  
  if (year == 2000) {
    # 2000 Census specific variables
    # Population by age
    age_vars <- c(
      "P012001",  # Total
      "P012003", "P012027",  # Under 5 (M, F)
      "P012004", "P012028",  # 5-9 (M, F)
      "P012005", "P012029",  # 10-14 (M, F)
      "P012006", "P012007", "P012030", "P012031"  # 15-17 (M, F)
    )
    
    # Households with children
    household_vars <- c(
      "P020001",  # Total households
      "P020002"   # Households with individuals under 18
    )
    
  } else if (year %in% c(2010, 2020)) {
    # 2010 and 2020 use same table structure
    age_vars <- c(
      "P12_001N",  # Total
      "P12_003N", "P12_027N",  # Under 5 (M, F)
      "P12_004N", "P12_028N",  # 5-9 (M, F)
      "P12_005N", "P12_029N",  # 10-14 (M, F)
      "P12_006N", "P12_007N", "P12_030N", "P12_031N"  # 15-17 (M, F)
    )
    
    household_vars <- c(
      "P20_001N",  # Total households
      "P20_002N"   # Households with individuals under 18
    )
  }
  
  # Get population data
  pop_data <- safe_get_data(get_decennial,
                            geography = GEOGRAPHY,
                            variables = age_vars,
                            year = year,
                            output = "wide")
  
  # Get household data
  hh_data <- safe_get_data(get_decennial,
                           geography = GEOGRAPHY,
                           variables = household_vars,
                           year = year,
                           output = "wide")
  
  # Process and combine
  combined_data <- NULL
  
  if (!is.null(pop_data) && !is.null(hh_data)) {
    # Calculate derived variables based on year
    if (year == 2000) {
      pop_processed <- pop_data %>%
        mutate(
          total_population = P012001,
          children_under_5 = P012003 + P012027,
          children_5_to_9 = P012004 + P012028,
          children_10_to_14 = P012005 + P012029,
          children_15_to_17 = P012006 + P012007 + P012030 + P012031,
          children_under_18 = children_under_5 + children_5_to_9 + 
            children_10_to_14 + children_15_to_17,
          pct_under_18 = (children_under_18 / total_population) * 100
        )
      
      hh_processed <- hh_data %>%
        mutate(
          total_households = P020001,
          households_with_children = P020002,
          pct_households_with_children = (households_with_children / total_households) * 100
        )
      
    } else {
      # 2010 and 2020
      pop_processed <- pop_data %>%
        mutate(
          total_population = P12_001N,
          children_under_5 = P12_003N + P12_027N,
          children_5_to_9 = P12_004N + P12_028N,
          children_10_to_14 = P12_005N + P12_029N,
          children_15_to_17 = P12_006N + P12_007N + P12_030N + P12_031N,
          children_under_18 = children_under_5 + children_5_to_9 + 
            children_10_to_14 + children_15_to_17,
          pct_under_18 = (children_under_18 / total_population) * 100
        )
      
      hh_processed <- hh_data %>%
        mutate(
          total_households = P20_001N,
          households_with_children = P20_002N,
          pct_households_with_children = (households_with_children / total_households) * 100
        )
    }
    
    # Combine population and household data
    combined_data <- pop_processed %>%
      select(GEOID, NAME, total_population, children_under_18, pct_under_18,
             children_under_5) %>%
      left_join(
        hh_processed %>% 
          select(GEOID, total_households, households_with_children, 
                 pct_households_with_children),
        by = "GEOID"
      ) %>%
      mutate(
        year = year,
        data_source = "DECENNIAL",
        collection_date = Sys.Date()
      )
    
    log_message(sprintf("  Successfully retrieved %d MSAs for Decennial %d", 
                        nrow(combined_data), year))
  }
  
  return(combined_data)
}

# MAIN EXECUTION
# =============================================================================
log_message("=== STARTING DATA COLLECTION ===")

# Progress tracking
with_progress({
  p <- progressor(steps = length(ACS_START_YEAR:ACS_END_YEAR) + 
                    length(DECENNIAL_YEARS))
  
  # Collect ACS data
  log_message("Collecting ACS data...")
  acs_data_list <- future_map(ACS_START_YEAR:ACS_END_YEAR, function(yr) {
    p()
    # Try ACS-1 first, fall back to ACS-5 if needed
    data <- get_acs_data(yr, use_acs1 = TRUE)
    if (is.null(data) || nrow(data) < 50) {  # If too few MSAs, try ACS-5
      log_message(sprintf("Limited ACS-1 data for %d, trying ACS-5...", yr))
      data_acs5 <- get_acs_data(yr, use_acs1 = FALSE)
      if (!is.null(data_acs5) && (is.null(data) || nrow(data_acs5) > nrow(data))) {
        data <- data_acs5
      }
    }
    return(data)
  }, .options = furrr_options(seed = TRUE))
  
  # Collect Decennial data
  log_message("Collecting Decennial Census data...")
  decennial_data_list <- map(DECENNIAL_YEARS, function(yr) {
    p()
    get_decennial_data(yr)
  })
})

# COMBINE ALL DATA
# =============================================================================
log_message("Combining all data sources...")

# Remove NULL entries
acs_data_clean <- acs_data_list[!sapply(acs_data_list, is.null)]
decennial_data_clean <- decennial_data_list[!sapply(decennial_data_list, is.null)]

# Combine
all_data <- bind_rows(acs_data_clean, decennial_data_clean) %>%
  arrange(GEOID, year)

# DATA QUALITY SUMMARY
# =============================================================================
log_message("=== DATA COLLECTION COMPLETE ===")
log_message(sprintf("Total observations: %d", nrow(all_data)))
log_message(sprintf("Unique MSAs: %d", n_distinct(all_data$GEOID)))
log_message(sprintf("Years covered: %s", paste(sort(unique(all_data$year)), collapse = ", ")))

# Create summary by year and source
summary_by_year <- all_data %>%
  group_by(year, data_source) %>%
  summarise(
    n_msas = n(),
    avg_population = mean(total_population, na.rm = TRUE),
    avg_pct_under_18 = mean(pct_under_18, na.rm = TRUE),
    avg_pct_hh_children = mean(pct_households_with_children, na.rm = TRUE),
    fertility_coverage = sum(!is.na(birth_rate_per_1000)),
    cwr_coverage = sum(!is.na(child_woman_ratio)),
    .groups = "drop"
  ) %>%
  arrange(year)

print(summary_by_year)

# SAVE OUTPUTS
# =============================================================================
log_message("Saving outputs...")

# Save main panel
saveRDS(all_data, file.path(output_dir, "msa_panel_complete.rds"))
write_csv(all_data, file.path(output_dir, "msa_panel_complete.csv"))

# Save summary statistics
write_csv(summary_by_year, file.path(output_dir, "summary_by_year.csv"))

# Create validation report
validation_report <- all_data %>%
  group_by(year, data_source) %>%
  summarise(
    across(where(is.numeric), 
           list(
             n_valid = ~sum(!is.na(.)),
             pct_missing = ~mean(is.na(.)) * 100,
             mean = ~mean(., na.rm = TRUE),
             sd = ~sd(., na.rm = TRUE),
             min = ~min(., na.rm = TRUE),
             max = ~max(., na.rm = TRUE)
           ),
           .names = "{.col}__{.fn}"),
    .groups = "drop"
  )

write_csv(validation_report, file.path(output_dir, "data_validation_report.csv"))

# Top MSAs by household children percentage (latest year)
latest_year <- max(all_data$year)
top_msas <- all_data %>%
  filter(year == latest_year) %>%
  arrange(desc(pct_households_with_children)) %>%
  select(NAME, pct_households_with_children, pct_under_18, 
         birth_rate_per_1000, child_woman_ratio, total_population) %>%
  head(20)

log_message("\nTop 20 MSAs by % Households with Children (Latest Year):")
print(top_msas)
write_csv(top_msas, file.path(output_dir, "top_msas_households_children.csv"))

# Check fertility data coverage
fertility_coverage <- all_data %>%
  group_by(year) %>%
  summarise(
    total_msas = n(),
    has_b13002 = sum(!is.na(birth_rate_per_1000)),
    has_cwr = sum(!is.na(child_woman_ratio)),
    pct_b13002 = has_b13002 / total_msas * 100,
    pct_cwr = has_cwr / total_msas * 100
  )

log_message("\nFertility Data Coverage by Year:")
print(fertility_coverage)
write_csv(fertility_coverage, file.path(output_dir, "fertility_coverage.csv"))

# Create panel balance check
panel_balance <- all_data %>%
  group_by(GEOID) %>%
  summarise(
    n_years = n(),
    years_list = paste(sort(year), collapse = ", "),
    first_year = min(year),
    last_year = max(year)
  ) %>%
  arrange(desc(n_years))

# Save balanced panel (MSAs with all years)
max_years <- max(panel_balance$n_years)
balanced_msas <- panel_balance %>%
  filter(n_years == max_years) %>%
  pull(GEOID)

balanced_panel <- all_data %>%
  filter(GEOID %in% balanced_msas)

log_message(sprintf("\nBalanced panel: %d MSAs with complete data", length(balanced_msas)))
if (nrow(balanced_panel) > 0) {
  saveRDS(balanced_panel, file.path(output_dir, "msa_panel_balanced.rds"))
}

# Final log message
log_message(sprintf("\nAll outputs saved to: %s", output_dir))
log_message("=== SCRIPT EXECUTION COMPLETE ===")

# Display final summary
cat("\n========== FINAL SUMMARY ==========\n")
cat(sprintf("Total MSA-year observations: %d\n", nrow(all_data)))
cat(sprintf("Unique MSAs: %d\n", n_distinct(all_data$GEOID)))
cat(sprintf("Time span: %d-%d\n", min(all_data$year), max(all_data$year)))
cat(sprintf("Data sources: %s\n", paste(unique(all_data$data_source), collapse = ", ")))
cat(sprintf("Average observations per MSA: %.1f\n", nrow(all_data) / n_distinct(all_data$GEOID)))
cat(sprintf("Balanced panel size: %d MSAs\n", length(balanced_msas)))
cat("====================================\n")