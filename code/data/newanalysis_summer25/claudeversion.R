# =============================================================================
# MSA PANEL SCRIPT V7.0 - COMPLETE OVERHAUL WITH FERTILITY COMPARISON
# Simplified version using stable B-tables and direct variables
# Includes multiple fertility measures for comparison:
# - S1301 direct rate (2010+)
# - B13002 calculated rate (all years)
# - Child-Woman Ratio (CWR)
# =============================================================================

# SETUP
# =============================================================================
cat("Setting up libraries and environment...\n")
suppressPackageStartupMessages({
  if (!require(tidycensus)) install.packages("tidycensus"); library(tidycensus)
  if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
  if (!require(purrr)) install.packages("purrr"); library(purrr)
})

# Check for Census API key
if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("Census API Key not found. Please set it using census_api_key('YOUR_KEY', install=TRUE)")
}
options(tigris_use_cache = TRUE)

# CONFIGURATION
# =============================================================================
START_YEAR <- 2009
END_YEAR <- 2022
GEOGRAPHY <- "metropolitan statistical area/micropolitan statistical area"

cat("\nThis script collects three fertility measures for comparison:\n")
cat("1. S1301 Direct Rate (available 2010+)\n")
cat("2. B13002 Calculated Rate (available all years)\n")
cat("3. Child-Woman Ratio (traditional demographic measure)\n\n")

# MAIN DATA COLLECTION FUNCTION
# =============================================================================
get_msa_panel_data <- function(year) {
  cat(paste0("\nFetching data for year ", year, "...\n"))
  
  # Initialize list to store all data components
  data_list <- list()
  
  # 1. BASIC DEMOGRAPHICS & POPULATION
  # =============================================================================
  cat("  1. Population and basic demographics...\n")
  pop_vars <- c(
    total_population = "B01003_001E",
    # Age breakdown for children
    pop_under_5 = "B01001_003E",  # Male under 5
    pop_under_5_f = "B01001_027E", # Female under 5
    pop_5_to_9 = "B01001_004E",    # Male 5-9
    pop_5_to_9_f = "B01001_028E",  # Female 5-9
    pop_10_to_14 = "B01001_005E",  # Male 10-14
    pop_10_to_14_f = "B01001_029E", # Female 10-14
    pop_15_to_17 = "B01001_006E",  # Male 15-17
    pop_15_to_17_f = "B01001_030E"  # Female 15-17
  )
  
  pop_data <- tryCatch({
    get_acs(geography = GEOGRAPHY, 
            variables = pop_vars, 
            year = year,
            survey = "acs5", 
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching population data: ", e$message, "\n"))
    return(NULL)
  })
  
  if (!is.null(pop_data)) {
    pop_data <- pop_data %>%
      mutate(
        children_under_5 = pop_under_5 + pop_under_5_f,
        children_5_to_9 = pop_5_to_9 + pop_5_to_9_f,
        children_10_to_14 = pop_10_to_14 + pop_10_to_14_f,
        children_15_to_17 = pop_15_to_17 + pop_15_to_17_f,
        children_under_18 = children_under_5 + children_5_to_9 + children_10_to_14 + children_15_to_17,
        pct_population_under_18 = (children_under_18 / total_population) * 100
      ) %>%
      select(GEOID, NAME, total_population, children_under_5, children_5_to_9, 
             children_10_to_14, children_15_to_17, children_under_18, pct_population_under_18)
    
    data_list[["population"]] <- pop_data
  }
  
  # 2. HOUSEHOLDS WITH CHILDREN
  # =============================================================================
  cat("  2. Households with children...\n")
  household_vars <- c(
    total_households = "B11005_001E",
    households_with_children = "B11005_002E"  # Households with one or more people under 18
  )
  
  household_data <- tryCatch({
    get_acs(geography = GEOGRAPHY,
            variables = household_vars,
            year = year,
            survey = "acs5",
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching household data: ", e$message, "\n"))
    return(NULL)
  })
  
  if (!is.null(household_data)) {
    household_data <- household_data %>%
      mutate(
        share_households_with_children = (households_with_children / total_households) * 100
      ) %>%
      select(GEOID, total_households, households_with_children, share_households_with_children)
    
    data_list[["households"]] <- household_data
  }
  
  # 3. FAMILIES BY AGE OF CHILDREN
  # =============================================================================
  cat("  3. Families by age of children...\n")
  family_vars <- c(
    total_families = "B11003_001E",
    families_with_children_under18 = "B11003_002E",
    families_with_children_under6_only = "B11003_003E",
    families_with_children_under6_and_6to17 = "B11003_004E",
    families_with_children_6to17_only = "B11003_005E"
  )
  
  family_data <- tryCatch({
    get_acs(geography = GEOGRAPHY,
            variables = family_vars,
            year = year,
            survey = "acs5",
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching family data: ", e$message, "\n"))
    return(NULL)
  })
  
  if (!is.null(family_data)) {
    family_data <- family_data %>%
      mutate(
        share_families_with_children = (families_with_children_under18 / total_families) * 100,
        share_families_children_under6_only = (families_with_children_under6_only / total_families) * 100,
        share_families_children_both_ages = (families_with_children_under6_and_6to17 / total_families) * 100,
        share_families_children_6to17_only = (families_with_children_6to17_only / total_families) * 100
      ) %>%
      select(GEOID, total_families, families_with_children_under18, share_families_with_children,
             share_families_children_under6_only, share_families_children_both_ages, 
             share_families_children_6to17_only)
    
    data_list[["families"]] <- family_data
  }
  
  # 4. DETAILED CHILDREN COUNT BY AGE
  # =============================================================================
  cat("  4. Children count by age...\n")
  children_vars <- c(
    total_own_children = "B09002_001E",
    children_under_3 = "B09002_002E",
    children_3_to_4 = "B09002_003E",
    children_5 = "B09002_004E",
    children_6_to_8 = "B09002_005E",
    children_9_to_11 = "B09002_006E",
    children_12_to_14 = "B09002_007E",
    children_15_to_17 = "B09002_008E"
  )
  
  children_data <- tryCatch({
    get_acs(geography = GEOGRAPHY,
            variables = children_vars,
            year = year,
            survey = "acs5",
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching children data: ", e$message, "\n"))
    return(NULL)
  })
  
  if (!is.null(children_data)) {
    children_data <- children_data %>%
      mutate(
        children_under_6 = children_under_3 + children_3_to_4 + children_5,
        children_6_to_17 = children_6_to_8 + children_9_to_11 + children_12_to_14 + children_15_to_17,
        pct_children_under_6 = (children_under_6 / total_own_children) * 100,
        pct_children_6_to_17 = (children_6_to_17 / total_own_children) * 100
      ) %>%
      select(GEOID, total_own_children, children_under_6, children_6_to_17,
             pct_children_under_6, pct_children_6_to_17)
    
    data_list[["children"]] <- children_data
  }
  
  # 5. FERTILITY - BOTH METHODS FOR COMPARISON
  # =============================================================================
  cat("  5. Fertility data (both methods for comparison)...\n")
  
  # Method 1: S1301 Direct Rate (2010+)
  fertility_s1301 <- NULL
  if (year >= 2010) {
    fertility_s1301 <- tryCatch({
      s1301 <- get_acs(geography = GEOGRAPHY,
                       table = "S1301",
                       year = year,
                       survey = "acs5",
                       cache_table = TRUE)
      
      # S1301_C03_001 is the fertility rate per 1,000 women
      s1301 %>%
        filter(variable == "S1301_C03_001E") %>%
        select(GEOID, estimate) %>%
        rename(birth_rate_s1301 = estimate) %>%
        mutate(
          birth_rate_s1301_moe = NA_real_, # Could get MOE if needed
          fertility_s1301_source = "S1301 Direct"
        )
    }, error = function(e) {
      cat(paste0("    Error fetching S1301 fertility data: ", e$message, "\n"))
      return(NULL)
    })
  }
  
  # Method 2: B13002 Calculated Rate (available all years)
  fertility_b13002 <- tryCatch({
    fert_vars <- c(
      women_15_50 = "B13002_001E",
      women_with_birth = "B13002_002E"
    )
    
    fert <- get_acs(geography = GEOGRAPHY,
                    variables = fert_vars,
                    year = year,
                    survey = "acs5",
                    output = "wide")
    
    fert %>%
      mutate(
        birth_rate_b13002 = ifelse(women_15_50 > 0 & !is.na(women_with_birth), 
                                   (women_with_birth / women_15_50) * 1000, 
                                   NA_real_),
        fertility_b13002_source = "B13002 Calculated"
      ) %>%
      select(GEOID, women_15_50, women_with_birth, birth_rate_b13002, fertility_b13002_source)
  }, error = function(e) {
    cat(paste0("    Error fetching B13002 fertility data: ", e$message, "\n"))
    return(NULL)
  })
  
  # Method 3: Child-Woman Ratio (CWR) - traditional demographic measure
  cwr_data <- tryCatch({
    # Women aged 15-49 for denominator
    cwr_vars <- c(
      f_15_17 = "B01001_030E", f_18_19 = "B01001_031E", f_20 = "B01001_032E",
      f_21 = "B01001_033E", f_22_24 = "B01001_034E", f_25_29 = "B01001_035E",
      f_30_34 = "B01001_036E", f_35_39 = "B01001_037E", f_40_44 = "B01001_038E",
      f_45_49 = "B01001_039E",
      # Children under 5 (already fetched in population, but need for standalone)
      children_0_4_m = "B01001_003E", 
      children_0_4_f = "B01001_027E"
    )
    
    cwr <- get_acs(geography = GEOGRAPHY,
                   variables = cwr_vars,
                   year = year,
                   survey = "acs5",
                   output = "wide")
    
    # Calculate CWR
    cwr %>%
      mutate(
        women_15_49 = rowSums(select(., f_15_17:f_45_49), na.rm = FALSE),
        children_under_5_cwr = children_0_4_m + children_0_4_f,
        child_woman_ratio = ifelse(women_15_49 > 0 & !is.na(children_under_5_cwr),
                                   (children_under_5_cwr / women_15_49) * 1000,
                                   NA_real_)
      ) %>%
      select(GEOID, women_15_49, children_under_5_cwr, child_woman_ratio)
  }, error = function(e) {
    cat(paste0("    Error calculating CWR: ", e$message, "\n"))
    return(NULL)
  })
  
  # Combine all fertility measures
  fertility_combined <- NULL
  
  # Start with B13002 (available all years)
  if (!is.null(fertility_b13002)) {
    fertility_combined <- fertility_b13002
  }
  
  # Add S1301 if available
  if (!is.null(fertility_s1301) && !is.null(fertility_combined)) {
    fertility_combined <- fertility_combined %>%
      left_join(fertility_s1301, by = "GEOID")
  } else if (!is.null(fertility_s1301)) {
    fertility_combined <- fertility_s1301
  }
  
  # Add CWR if available
  if (!is.null(cwr_data) && !is.null(fertility_combined)) {
    fertility_combined <- fertility_combined %>%
      left_join(cwr_data, by = "GEOID")
  } else if (!is.null(cwr_data)) {
    fertility_combined <- cwr_data
  }
  
  # Create a primary fertility measure (prefer S1301 when available)
  if (!is.null(fertility_combined)) {
    fertility_combined <- fertility_combined %>%
      mutate(
        # Primary measure: S1301 if available (2010+), otherwise B13002
        birth_rate_primary = if("birth_rate_s1301" %in% names(.)) {
          coalesce(birth_rate_s1301, birth_rate_b13002)
        } else {
          birth_rate_b13002
        },
        fertility_primary_source = if("birth_rate_s1301" %in% names(.)) {
          case_when(
            !is.na(birth_rate_s1301) ~ "S1301 Direct",
            !is.na(birth_rate_b13002) ~ "B13002 Calculated",
            TRUE ~ NA_character_
          )
        } else {
          if("birth_rate_b13002" %in% names(.)) {
            ifelse(!is.na(birth_rate_b13002), "B13002 Calculated", NA_character_)
          } else {
            NA_character_
          }
        }
      )
    
    data_list[["fertility"]] <- fertility_combined
  }
  
  # 6. EDUCATION
  # =============================================================================
  cat("  6. Education data...\n")
  edu_vars <- c(
    pop_25_plus = "B15003_001E",
    bachelors = "B15003_022E",
    masters = "B15003_023E",
    professional = "B15003_024E",
    doctorate = "B15003_025E"
  )
  
  edu_data <- tryCatch({
    get_acs(geography = GEOGRAPHY,
            variables = edu_vars,
            year = year,
            survey = "acs5",
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching education data: ", e$message, "\n"))
    return(NULL)
  })
  
  if (!is.null(edu_data)) {
    edu_data <- edu_data %>%
      mutate(
        bachelors_or_higher = bachelors + masters + professional + doctorate,
        pct_bachelors_or_higher = (bachelors_or_higher / pop_25_plus) * 100
      ) %>%
      select(GEOID, pop_25_plus, bachelors_or_higher, pct_bachelors_or_higher)
    
    data_list[["education"]] <- edu_data
  }
  
  # 7. INCOME
  # =============================================================================
  cat("  7. Income data...\n")
  income_vars <- c(
    median_household_income = "B19013_001E",
    median_family_income = "B19113_001E",
    median_earnings = "B20002_001E"
  )
  
  income_data <- tryCatch({
    get_acs(geography = GEOGRAPHY,
            variables = income_vars,
            year = year,
            survey = "acs5",
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching income data: ", e$message, "\n"))
    return(NULL)
  })
  
  if (!is.null(income_data)) {
    data_list[["income"]] <- income_data %>%
      select(GEOID, median_household_income, median_family_income, median_earnings)
  }
  
  # 8. HOUSING
  # =============================================================================
  cat("  8. Housing data...\n")
  housing_vars <- c(
    median_home_value = "B25077_001E",
    median_gross_rent = "B25064_001E",
    total_units = "B25001_001E",
    occupied_units = "B25002_002E",
    owner_occupied = "B25003_002E",
    renter_occupied = "B25003_003E"
  )
  
  housing_data <- tryCatch({
    get_acs(geography = GEOGRAPHY,
            variables = housing_vars,
            year = year,
            survey = "acs5",
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching housing data: ", e$message, "\n"))
    return(NULL)
  })
  
  if (!is.null(housing_data)) {
    housing_data <- housing_data %>%
      mutate(
        homeownership_rate = (owner_occupied / occupied_units) * 100,
        rental_rate = (renter_occupied / occupied_units) * 100,
        occupancy_rate = (occupied_units / total_units) * 100
      ) %>%
      select(GEOID, median_home_value, median_gross_rent, total_units,
             occupied_units, homeownership_rate, rental_rate, occupancy_rate)
    
    data_list[["housing"]] <- housing_data
  }
  
  # 9. HOUSING COST BURDEN (using B25070 for renters, B25091 for owners)
  # =============================================================================
  cat("  9. Housing cost burden...\n")
  
  # Renter cost burden
  renter_burden_vars <- c(
    renter_total = "B25070_001E",
    renter_30_34pct = "B25070_007E",
    renter_35_39pct = "B25070_008E",
    renter_40_49pct = "B25070_009E",
    renter_50plus = "B25070_010E"
  )
  
  renter_burden <- tryCatch({
    get_acs(geography = GEOGRAPHY,
            variables = renter_burden_vars,
            year = year,
            survey = "acs5",
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching renter burden data: ", e$message, "\n"))
    return(NULL)
  })
  
  # Owner cost burden
  owner_burden_vars <- c(
    owner_total = "B25091_001E",
    owner_30_34pct = "B25091_007E",
    owner_35_39pct = "B25091_008E", 
    owner_40_49pct = "B25091_009E",
    owner_50plus = "B25091_010E"
  )
  
  owner_burden <- tryCatch({
    get_acs(geography = GEOGRAPHY,
            variables = owner_burden_vars,
            year = year,
            survey = "acs5",
            output = "wide")
  }, error = function(e) {
    cat(paste0("    Error fetching owner burden data: ", e$message, "\n"))
    return(NULL)
  })
  
  # Combine burden data
  burden_data <- NULL
  if (!is.null(renter_burden) && !is.null(owner_burden)) {
    burden_data <- renter_burden %>%
      select(GEOID, starts_with("renter_")) %>%
      left_join(owner_burden %>% select(GEOID, starts_with("owner_")), by = "GEOID") %>%
      mutate(
        renter_burdened_30plus = renter_30_34pct + renter_35_39pct + renter_40_49pct + renter_50plus,
        owner_burdened_30plus = owner_30_34pct + owner_35_39pct + owner_40_49pct + owner_50plus,
        pct_renters_burdened = (renter_burdened_30plus / renter_total) * 100,
        pct_owners_burdened = (owner_burdened_30plus / owner_total) * 100
      ) %>%
      select(GEOID, pct_renters_burdened, pct_owners_burdened)
    
    data_list[["cost_burden"]] <- burden_data
  }
  
  # COMBINE ALL DATA
  # =============================================================================
  if (length(data_list) == 0) {
    cat(paste0("  No data retrieved for year ", year, "\n"))
    return(NULL)
  }
  
  # Start with the first non-null dataset
  combined_data <- NULL
  for (df in data_list) {
    if (!is.null(df) && "GEOID" %in% names(df)) {
      if (is.null(combined_data)) {
        combined_data <- df
      } else {
        # Remove NAME column from subsequent joins to avoid duplicates
        cols_to_join <- setdiff(names(df), c("NAME"))
        combined_data <- combined_data %>%
          left_join(df %>% select(all_of(cols_to_join)), by = "GEOID")
      }
    }
  }
  
  if (!is.null(combined_data)) {
    combined_data$year <- year
    
    # Add calculated fields
    if (all(c("median_home_value", "median_household_income") %in% names(combined_data))) {
      combined_data <- combined_data %>%
        mutate(price_to_income_ratio = median_home_value / median_household_income)
    }
    
    if (all(c("median_gross_rent", "median_household_income") %in% names(combined_data))) {
      combined_data <- combined_data %>%
        mutate(rent_to_income_ratio = (median_gross_rent * 12) / median_household_income)
    }
  }
  
  return(combined_data)
}

# MAIN EXECUTION
# =============================================================================
cat("\n=== STARTING MSA PANEL DATA COLLECTION ===\n")
cat(paste0("Years: ", START_YEAR, " to ", END_YEAR, "\n"))

# Collect data for all years
all_years_data <- map(START_YEAR:END_YEAR, function(yr) {
  data <- get_msa_panel_data(yr)
  if (!is.null(data)) {
    cat(paste0("  Successfully retrieved ", nrow(data), " MSAs for year ", yr, "\n"))
  }
  return(data)
})

# Remove NULL entries and combine
all_years_data <- all_years_data[!sapply(all_years_data, is.null)]

if (length(all_years_data) == 0) {
  stop("No data retrieved for any year!")
}

# Combine into panel
msa_panel <- bind_rows(all_years_data)

cat(paste0("\n=== PANEL CREATION COMPLETE ===\n"))
cat(paste0("Total observations: ", nrow(msa_panel), "\n"))
cat(paste0("Years included: ", paste(sort(unique(msa_panel$year)), collapse = ", "), "\n"))
cat(paste0("Number of MSAs: ", length(unique(msa_panel$GEOID)), "\n"))
cat(paste0("Variables: ", ncol(msa_panel), "\n"))

# SAVE OUTPUT
# =============================================================================
output_dir <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/R/msapanelv7"
dir.create(output_dir, showWarnings = FALSE)

# Save full panel
saveRDS(msa_panel, file.path(output_dir, "msa_panel_v7_full.rds"))
write_csv(msa_panel, file.path(output_dir, "msa_panel_v7_full.csv"))

# Create summary statistics for latest year
latest_year <- max(msa_panel$year)
latest_data <- msa_panel %>% filter(year == latest_year)

summary_stats <- latest_data %>%
  select(where(is.numeric), -GEOID, -year) %>%
  summarise(across(everything(), 
                   list(mean = ~mean(.x, na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE),
                        min = ~min(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE),
                        n_obs = ~sum(!is.na(.x))),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(), names_to = c("variable", "stat"), names_sep = "__") %>%
  pivot_wider(names_from = stat, values_from = value)

write_csv(summary_stats, file.path(output_dir, paste0("summary_stats_", latest_year, ".csv")))

# Top MSAs by share of households with children
top_child_msas <- latest_data %>%
  filter(!is.na(share_households_with_children)) %>%
  arrange(desc(share_households_with_children)) %>%
  select(NAME, share_households_with_children, share_families_with_children,
         pct_population_under_18, birth_rate_primary, 
         any_of(c("birth_rate_b13002", "birth_rate_s1301", "child_woman_ratio")),
         median_household_income, pct_bachelors_or_higher) %>%
  head(20)

write_csv(top_child_msas, file.path(output_dir, paste0("top_msas_children_", latest_year, ".csv")))

cat(paste0("\nOutput saved to: ", output_dir, "/\n"))

# Print sample of data
cat("\n=== SAMPLE OF FINAL DATA (Latest Year, Top 5 MSAs by Population) ===\n")
sample_data <- latest_data %>%
  arrange(desc(total_population)) %>%
  head(5) %>%
  select(NAME, year, total_population, share_households_with_children,
         share_families_with_children, pct_population_under_18,
         birth_rate_primary, fertility_primary_source,
         any_of(c("birth_rate_b13002", "birth_rate_s1301", "child_woman_ratio")),
         median_household_income, median_home_value, price_to_income_ratio)

print(sample_data)


# Save fertility comparison by year
fertility_by_year <- msa_panel %>%
  group_by(year) %>%
  summarise(
    n_total = n(),
    n_b13002 = sum(!is.na(birth_rate_b13002)),
    n_s1301 = if("birth_rate_s1301" %in% names(msa_panel)) sum(!is.na(birth_rate_s1301)) else 0,
    n_cwr = sum(!is.na(child_woman_ratio)),
    mean_b13002 = mean(birth_rate_b13002, na.rm = TRUE),
    mean_s1301 = if("birth_rate_s1301" %in% names(msa_panel)) mean(birth_rate_s1301, na.rm = TRUE) else NA_real_,
    mean_cwr = mean(child_woman_ratio, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(fertility_by_year, file.path(output_dir, "fertility_measures_by_year.csv"))

cat("\n=== SCRIPT COMPLETE ===\n")