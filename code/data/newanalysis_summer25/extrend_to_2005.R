# =============================================================================
# FINAL CONSOLIDATED SCRIPT: FERTILITY, HOUSING, AND AGGLOMERATION ANALYSIS (V5.1 - Extended Time Period)
# =============================================================================

# PART 0: SETUP & LIBRARIES
# =============================================================================
cat("PART 0: Setting up libraries and environment...\n")
suppressPackageStartupMessages({
  if (!require(tidycensus)) install.packages("tidycensus"); library(tidycensus)
  if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
  if (!require(scales)) install.packages("scales"); library(scales)
  if (!require(ggrepel)) install.packages("ggrepel"); library(ggrepel)
  if (!require(viridis)) install.packages("viridis"); library(viridis)
  if (!require(patchwork)) install.packages("patchwork"); library(patchwork)
  if (!require(knitr)) install.packages("knitr"); library(knitr)
  if (!require(modelsummary)) install.packages("modelsummary"); library(modelsummary)
  if (!require(plm)) install.packages("plm"); library(plm)
  if (!require(GGally)) install.packages("GGally"); library(GGally)
})

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0.5),
                  legend.position = "bottom",
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))

# census_api_key("YOUR_API_KEY", install = TRUE, overwrite = TRUE) # Make sure your key is set
if (Sys.getenv("CENSUS_API_KEY") == "") stop("Census API Key not found. Please set it using census_api_key('YOUR_KEY', install=TRUE).")

options(tigris_use_cache = TRUE)
output_main_dir <- "consolidated_fertility_analysis_v5_extended" # Incremented version for clarity
dir.create(output_main_dir, showWarnings = FALSE)

# =============================================================================
# PART 1: ENHANCED DATA COLLECTION FUNCTION
# =============================================================================
cat("PART 1: Defining data collection function get_msa_data_comprehensive()...\n")
get_msa_data_comprehensive <- function(year_to_fetch, survey_type = "acs5") {
  cat(paste0("\n--- Fetching data for year: ", year_to_fetch, ", survey: ", survey_type, " ---\n"))
  base_geographies_df <- NULL 
  
  safe_get_acs <- function(variables_list = NULL, table_name = NULL, processing_fun, var_label, 
                           is_base_geo_fetch = FALSE, expected_output_cols = NULL) {
    cat(paste("Attempting to fetch", var_label, "... "))
    tryCatch({
      data_raw <- if (!is.null(table_name)) {
        get_acs(geography = "metropolitan statistical area/micropolitan statistical area",
                table = table_name, year = year_to_fetch, survey = survey_type, cache_table = TRUE, show_call = FALSE)
      } else { 
        get_acs(geography = "metropolitan statistical area/micropolitan statistical area",
                variables = variables_list, year = year_to_fetch, survey = survey_type, cache_table = TRUE, output = "wide", show_call = FALSE)
      }
      if (nrow(data_raw) == 0) {cat("No data returned by API.\n"); data_to_return <- NULL}
      else {
        if (is_base_geo_fetch && ("GEOID" %in% names(data_raw)) && ("NAME" %in% names(data_raw))) {
          assign("base_geographies_df", data_raw %>% select(GEOID, NAME) %>% distinct(), envir = parent.env(environment())); cat("Base geographies captured. ")
        }
        processed_data <- processing_fun(data_raw)
        if(is.null(processed_data) || nrow(processed_data) == 0) {cat("Processing returned no data.\n"); data_to_return <- NULL}
        else {cat("Processing successful.\n"); data_to_return <- processed_data}
      }
      if (!is.null(expected_output_cols)) {
        if (is.null(data_to_return) || nrow(data_to_return) == 0) {
          if (!is.null(base_geographies_df) && nrow(base_geographies_df) > 0) {
            temp_df <- base_geographies_df %>% select(GEOID); for (col_name in expected_output_cols) { temp_df[[col_name]] <- NA_real_ }; data_to_return <- temp_df
          } else {data_to_return <- NULL}
        } else { 
          for (col_name in expected_output_cols) {if (!col_name %in% names(data_to_return)) data_to_return[[col_name]] <- NA_real_}
          if (!"GEOID" %in% names(data_to_return) && !is.null(base_geographies_df) && nrow(base_geographies_df) > 0){data_to_return <- base_geographies_df %>% select(GEOID) %>% left_join(data_to_return, by = character())}
          else if (!"GEOID" %in% names(data_to_return)) {data_to_return <- NULL}
        }
      }
      return(data_to_return)
    }, error = function(e) {
      cat(paste("Error for '", var_label, "': ", e$message, "\n"))
      if (!is.null(expected_output_cols) && !is.null(base_geographies_df) && nrow(base_geographies_df) > 0) {
        temp_df <- base_geographies_df %>% select(GEOID); for (col_name in expected_output_cols) { temp_df[[col_name]] <- NA_real_ }; return(temp_df)
      }
      return(NULL) 
    })
  }
  
  eo_population <- c("population", "NAME", "pop_lt5")
  process_population <- function(df_wide_pop) { 
    cols_needed_E <- c("B01003_001E", "B01001_003E", "B01001_027E") 
    for(col_n in cols_needed_E) { if(!col_n %in% names(df_wide_pop)) df_wide_pop[[col_n]] <- NA_real_ }
    df_wide_pop %>% transmute(GEOID, NAME, population = as.numeric(B01003_001E), pop_lt5 = as.numeric(B01001_003E) + as.numeric(B01001_027E)) 
  }
  population_data <- safe_get_acs(variables_list = c("B01003_001", "B01001_003", "B01001_027"), processing_fun = process_population, var_label = "Population", is_base_geo_fetch = TRUE, expected_output_cols = eo_population)
  if (is.null(population_data) || nrow(population_data) == 0 || !"GEOID" %in% names(population_data)) {cat("CRITICAL: Population data failed.\n"); return(NULL)}
  
  eo_households <- c("total_households", "family_households", "families_with_children_lt18", "families_with_children_lt6", "avg_household_size", "share_families_with_children_lt18", "share_families_with_children_lt6")
  process_dp02_households <- function(df_long_dp02) { 
    vars_needed_dp02_codes <- c("DP02_0001", "DP02_0002", "DP02_0003", "DP02_0009", "DP02_0016")
    df_wide <- df_long_dp02 %>% filter(variable %in% vars_needed_dp02_codes) %>% select(GEOID, variable, estimate) %>% pivot_wider(names_from = variable, values_from = estimate, values_fn = first)
    for(v_name in vars_needed_dp02_codes) { if(!v_name %in% names(df_wide)) df_wide[[v_name]] <- NA_real_ }
    df_wide %>% transmute(GEOID, total_households = as.numeric(DP02_0001), family_households = as.numeric(DP02_0002), families_with_children_lt18 = as.numeric(DP02_0003), families_with_children_lt6 = as.numeric(DP02_0009), avg_household_size = as.numeric(DP02_0016), share_families_with_children_lt18 = ifelse(!is.na(family_households) & family_households > 0, families_with_children_lt18 / family_households, NA_real_), share_families_with_children_lt6 = ifelse(!is.na(family_households) & family_households > 0, families_with_children_lt6 / family_households, NA_real_))
  }
  household_structure_data <- safe_get_acs(table_name = "DP02", processing_fun = process_dp02_households, var_label = "Household Structure (DP02)", expected_output_cols = eo_households)
  
  eo_econ <- c("median_household_income", "median_home_value", "median_gross_rent")
  econ_vars_codes <- c("B19013_001", "B25077_001", "B25064_001")
  process_econ <- function(df_wide_econ) { 
    cols_needed_econ_E <- paste0(econ_vars_codes, "E")
    for(col_n in cols_needed_econ_E) { if(!col_n %in% names(df_wide_econ)) df_wide_econ[[col_n]] <- NA_real_ }
    df_wide_econ %>% transmute(GEOID, median_household_income = as.numeric(B19013_001E), median_home_value = as.numeric(B25077_001E), median_gross_rent = as.numeric(B25064_001E))
  }
  economic_data <- safe_get_acs(variables_list = econ_vars_codes, processing_fun = process_econ, var_label = "Economic Indicators", expected_output_cols = eo_econ)
  
  eo_earnings <- c("median_earnings_workers")
  process_earnings <- function(df_wide_earn) { 
    if(!"B20002_001E" %in% names(df_wide_earn)) df_wide_earn[["B20002_001E"]] <- NA_real_
    df_wide_earn %>% transmute(GEOID, median_earnings_workers = as.numeric(B20002_001E)) 
  }
  earnings_data <- safe_get_acs(variables_list = "B20002_001", processing_fun = process_earnings, var_label = "Median Earnings", expected_output_cols = eo_earnings)
  
  eo_education <- c("college_pct")
  process_education_dp02 <- function(df_long_dp02_edu) { 
    # year_to_fetch is available from the parent function's environment
    current_year_for_codes <- year_to_fetch 
    
    if (current_year_for_codes < 2013) {
      # For ACS years before 2013 (e.g., 2009-2012 5-year data):
      # Denominator: DP02_0059 (Population 25 years and over)
      # Numerator components: DP02_0065 (Bachelor's degree) + DP02_0066 (Graduate or professional degree)
      vars_to_extract_for_pivot <- c("DP02_0059", "DP02_0065", "DP02_0066")
      cat(paste0(" [DP02 Education using codes for pre-2013 ACS (", current_year_for_codes, "): ", paste(vars_to_extract_for_pivot, collapse=", "), "] "))
    } else {
      # For ACS 2013 data and onwards:
      # Denominator: DP02_0059 (Population 25 years and over)
      # Numerator: DP02_0067 (Bachelor's degree or higher)
      vars_to_extract_for_pivot <- c("DP02_0059", "DP02_0067")
      cat(paste0(" [DP02 Education using codes for 2013+ ACS (", current_year_for_codes, "): ", paste(vars_to_extract_for_pivot, collapse=", "), "] "))
    }
    
    df_wide <- df_long_dp02_edu %>%
      filter(variable %in% vars_to_extract_for_pivot) %>%
      select(GEOID, variable, estimate) %>%
      pivot_wider(names_from = variable, values_from = estimate, values_fn = first)
    
    # Ensure all potentially needed columns exist, filling with NA_real_ if not.
    # This is important because not all columns are fetched/present for all year groups.
    all_possible_edu_vars <- c("DP02_0059", "DP02_0065", "DP02_0066", "DP02_0067")
    for(v_code in all_possible_edu_vars) {
      if(!v_code %in% names(df_wide)) {
        df_wide[[v_code]] <- NA_real_ 
      }
    }
    
    # Convert relevant columns to numeric
    df_wide <- df_wide %>%
      mutate(
        DP02_0059_num = as.numeric(DP02_0059),
        DP02_0065_num = as.numeric(DP02_0065), 
        DP02_0066_num = as.numeric(DP02_0066), 
        DP02_0067_num = as.numeric(DP02_0067)  
      )
    
    if (current_year_for_codes < 2013) {
      processed_data <- df_wide %>%
        transmute(
          GEOID,
          college_pct = ifelse(
            !is.na(DP02_0059_num) & DP02_0059_num > 0,
            ( (ifelse(is.na(DP02_0065_num), 0, DP02_0065_num) + ifelse(is.na(DP02_0066_num), 0, DP02_0066_num)) / DP02_0059_num ) * 100,
            NA_real_
          )
        )
    } else { # 2013 and later
      processed_data <- df_wide %>%
        transmute(
          GEOID,
          college_pct = ifelse(
            !is.na(DP02_0059_num) & DP02_0059_num > 0 & !is.na(DP02_0067_num), # Ensure numerator is not NA
            (DP02_0067_num / DP02_0059_num) * 100,
            NA_real_
          )
        )
    }
    return(processed_data)
  }
  education_data <- safe_get_acs(table_name = "DP02", processing_fun = process_education_dp02, var_label = "Education (DP02)", expected_output_cols = eo_education)
  
  eo_cwr <- c("cwr")
  cwr_vars_codes_for_call <- c("B01001_003", "B01001_027", "B01001_030", "B01001_031", "B01001_032", "B01001_033", "B01001_034", "B01001_035", "B01001_036", "B01001_037", "B01001_038", "B01001_039") 
  process_cwr <- function(df_wide_cwr) { 
    children_cols_E <- c("B01001_003E", "B01001_027E"); women_cols_E <- paste0(cwr_vars_codes_for_call[-(1:2)], "E"); all_needed_cols_cwr_E <- c(children_cols_E, women_cols_E)
    for(col_e in all_needed_cols_cwr_E) {if(!col_e %in% names(df_wide_cwr)) df_wide_cwr[[col_e]] <- NA_real_}
    df_wide_cwr %>% transmute(GEOID, children_0_4 = rowSums(select(., all_of(children_cols_E)), na.rm = FALSE), total_women_15_49 = rowSums(select(., all_of(women_cols_E)), na.rm = FALSE)) %>% mutate(cwr = ifelse(!is.na(total_women_15_49) & total_women_15_49 > 0 & !is.na(children_0_4) & children_0_4 >=0, (children_0_4 / total_women_15_49) * 1000, NA_real_)) %>% select(GEOID, cwr)
  }
  cwr_data <- safe_get_acs(variables_list = cwr_vars_codes_for_call, processing_fun = process_cwr, var_label = "Child-Woman Ratio", expected_output_cols = eo_cwr)
  
  eo_birth_rate <- c("birth_rate_15_50")
  process_recent_births_S1301_final <- function(df_long_s1301) { 
    current_year_for_s1301 <- year_to_fetch
    if (current_year_for_s1301 >= 2022) { # ACS 2022 5-year data (released end of 2023) and later
      var_needed_s1301_code <- "S1301_C04_001" # New code for "Rate per 1,000"
      cat(paste0(" [S1301 using new code: ", var_needed_s1301_code, " for year ", current_year_for_s1301, "] "))
    } else { # Older ACS data (e.g., up to 2021 5-year)
      var_needed_s1301_code <- "S1301_C03_001" # Old code for "Rate per 1,000"
      cat(paste0(" [S1301 using old code: ", var_needed_s1301_code, " for year ", current_year_for_s1301, "] "))
    }
    
    filtered_long <- df_long_s1301 %>%
      filter(variable == var_needed_s1301_code)
    
    if(nrow(filtered_long) == 0) { 
      unique_geoids_in_batch <- df_long_s1301 %>% distinct(GEOID)
      if (nrow(unique_geoids_in_batch) > 0) {
        return(unique_geoids_in_batch %>% mutate(birth_rate_15_50 = NA_real_))
      } else { 
        return(tibble(GEOID=character(), birth_rate_15_50=as.numeric(NA))) 
      }
    }
    
    transmuted_data <- filtered_long %>%
      transmute(GEOID, birth_rate_15_50 = as.numeric(estimate))
    
    return(transmuted_data)
  }
  recent_birth_data <- safe_get_acs(table_name = "S1301", 
                                    processing_fun = process_recent_births_S1301_final, 
                                    var_label = "Birth Rate (S1301)", 
                                    expected_output_cols = eo_birth_rate)
  
  all_data_for_year <- population_data 
  list_of_data_to_join <- list(household_structure_data, economic_data, earnings_data, education_data, cwr_data, recent_birth_data)
  for (dataset_to_join in list_of_data_to_join) {
    if (!is.null(dataset_to_join) && nrow(dataset_to_join) > 0 && "GEOID" %in% names(dataset_to_join)) {
      cols_in_current <- names(all_data_for_year); cols_in_new <- names(dataset_to_join); duplicate_cols_to_drop <- setdiff(intersect(cols_in_current, cols_in_new), "GEOID")
      if(length(duplicate_cols_to_drop) > 0) dataset_to_join <- dataset_to_join %>% select(-all_of(duplicate_cols_to_drop))
      all_data_for_year <- all_data_for_year %>% left_join(dataset_to_join, by = "GEOID")
    }
  }
  all_data_for_year <- all_data_for_year %>%
    mutate(year = as.integer(year_to_fetch), msa_name_full = NAME, msa_name = str_extract(NAME, "^[^,]+"), state_name_or_abbr = str_match(NAME, ",\\s*(.+?)(?:\\s+MSA|\\s+Micropolitan\\s+Statistical\\s+Area|\\s+Metropolitan\\s+Division|$)")[,2], cbsa_code = str_extract(GEOID, "\\d{5}$")) %>%
    filter(!str_detect(NAME, ", PR"), !is.na(population), population > 50000)
  final_schema_cols <- unique(c("GEOID", "msa_name_full", "msa_name", "state_name_or_abbr", "cbsa_code", "year", "population", "pop_lt5", "total_households", "family_households", "families_with_children_lt18", "families_with_children_lt6", "avg_household_size", "share_families_with_children_lt18", "share_families_with_children_lt6", "median_household_income", "median_home_value", "median_gross_rent", "median_earnings_workers", "college_pct", "cwr", "birth_rate_15_50"))
  for (col_s in final_schema_cols) {if (!col_s %in% names(all_data_for_year)) all_data_for_year[[col_s]] <- NA_real_}
  all_data_for_year <- all_data_for_year %>% select(all_of(final_schema_cols))
  if (nrow(all_data_for_year) == 0) {cat("WARNING: No MSAs remained for year", year_to_fetch, "\n"); return(NULL)}
  cat(paste0("--- Processed year: ", year_to_fetch, ". MSAs: ", nrow(all_data_for_year), ". Cols: ", ncol(all_data_for_year), " ---\n"))
  return(all_data_for_year)
}

# =============================================================================
# PART 2: DATA AGGREGATION & PANEL CREATION 
# =============================================================================
cat("\nPART 2: Creating MSA Panel...\n")
# MODIFIED: Extended target years to include earlier ACS 5-year data
# 2009 ACS 5-year covers 2005-2009 period
target_years_acs5 <- sort(unique(c(2009, 2011, 2013, 2014, 2016, 2018, 2022))) 
cat(paste0("Target years for ACS 5-year data: ", paste(target_years_acs5, collapse = ", "), "\n"))

msa_panel_list <- lapply(target_years_acs5, get_msa_data_comprehensive, survey_type = "acs5")
msa_panel_list <- msa_panel_list[!sapply(msa_panel_list, is.null)] 
if (length(msa_panel_list) == 0) {stop("FATAL: No data fetched for ANY year.")}
msa_panel_raw <- bind_rows(msa_panel_list)
cat(paste0("\nRaw panel created with ", nrow(msa_panel_raw), " MSA-year observations.\n"))
cat(paste0("Years in panel: ", paste(sort(unique(msa_panel_raw$year)), collapse = ", "), "\n"))

# =============================================================================
# PART 3: DATA CLEANING & FEATURE ENGINEERING 
# =============================================================================
cat("\nPART 3: Cleaning data and engineering features...\n")
fertility_vars_config <- list(
  birth_rate_15_50 = list(longname = "Birth Rate (per 1k Women 15-50)", plot_label_func = scales::comma_format(accuracy = 1)),
  share_families_with_children_lt6 = list(longname = "Share of Families w/ Children <6 yrs", plot_label_func = scales::percent_format(accuracy = 0.1)),
  cwr = list(longname = "Child-Woman Ratio (Children 0-4 per 1k W15-49)", plot_label_func = scales::comma_format(accuracy = 1))
)
fertility_vars_shortnames <- names(fertility_vars_config)
cols_to_ensure_numeric <- c("population", "total_households", "family_households", "families_with_children_lt18", "families_with_children_lt6", "avg_household_size", "share_families_with_children_lt18", "share_families_with_children_lt6", "median_household_income", "median_home_value", "median_gross_rent", "median_earnings_workers", "college_pct", "cwr", "birth_rate_15_50", "pop_lt5")
for(col_n in cols_to_ensure_numeric){if(!col_n %in% names(msa_panel_raw)) msa_panel_raw[[col_n]] <- NA_real_}

msa_panel_processed <- msa_panel_raw %>%
  mutate(across(all_of(cols_to_ensure_numeric), as.numeric)) %>%
  mutate(price_to_income_ratio = ifelse(median_household_income > 0, median_home_value / median_household_income, NA_real_), 
         annual_gross_rent = median_gross_rent * 12, 
         rent_to_income_ratio = ifelse(median_household_income > 0, annual_gross_rent / median_household_income, NA_real_), 
         log_population = ifelse(population > 0, log(population), NA_real_), 
         log_median_household_income = ifelse(median_household_income > 0, log(median_household_income), NA_real_), 
         fertility_level_category = case_when(birth_rate_15_50 >= 55 ~ "High Fertility (55+)", birth_rate_15_50 >= 50 ~ "Medium-High (50-55)", birth_rate_15_50 >= 45 ~ "Medium (45-50)", !is.na(birth_rate_15_50) ~ "Low Fertility (<45)", TRUE ~ NA_character_), 
         pct_pop_under_5 = ifelse(population > 0 & !is.na(pop_lt5), (pop_lt5 / population) * 100, NA_real_)
  )

msa_panel_final <- msa_panel_processed %>%
  filter(!is.na(population)) %>% group_by(year) %>%
  mutate(size_category = factor(case_when(population >= 5000000 ~ "1. Mega (>5M)", population >= 2500000 ~ "2. Very Large (2.5-5M)", population >= 1000000 ~ "3. Large (1-2.5M)", population >= 500000  ~ "4. Medium (0.5-1M)", population >= 250000  ~ "5. Small (0.25-0.5M)", TRUE ~ "6. Very Small (<0.25M)"), levels = c("6. Very Small (<0.25M)", "5. Small (0.25-0.5M)", "4. Medium (0.5-1M)", "3. Large (1-2.5M)", "2. Very Large (2.5-5M)", "1. Mega (>5M)"), ordered = TRUE), pop_rank = rank(-population)) %>% ungroup()

current_year_val <- NA_integer_ # Initialize
if (nrow(msa_panel_final) > 0 && any(!is.na(msa_panel_final$year))) {current_year_val <- max(msa_panel_final$year, na.rm = TRUE); msa_current_final <- msa_panel_final %>% filter(year == current_year_val); cat(paste0("\nCreated msa_current_final for year ", current_year_val, " with ", nrow(msa_current_final), " MSAs.\n"))} else {warning("msa_panel_final empty or no valid years."); msa_current_final <- tibble()}
saveRDS(msa_panel_final, file.path(output_main_dir, "msa_panel_final.rds")); saveRDS(msa_current_final, file.path(output_main_dir, "msa_current_final.rds")); cat("Saved .rds files\n")
cat("\n--- Summary of Key Variables (Current Year: ", if(!is.na(current_year_val)) current_year_val else "N/A" , ") ---\n"); current_summary_vars <- c(fertility_vars_shortnames, "population", "price_to_income_ratio", "rent_to_income_ratio", "college_pct", "pct_pop_under_5"); cols_for_summary <- intersect(current_summary_vars, names(msa_current_final)); if(length(cols_for_summary) > 0 && nrow(msa_current_final) > 0) {print(summary(msa_current_final %>% select(all_of(cols_for_summary))))} else {cat("No columns or data for summary.\n")}

analysis_plot_dir <- file.path(output_main_dir, "plots_and_tables"); dir.create(analysis_plot_dir, showWarnings = FALSE)
save_plot_obj <- function(plot_obj, filename_stem, width = 10, height = 8, dv_shortname = NULL, path = analysis_plot_dir) {filename <- if(!is.null(dv_shortname)) paste0(filename_stem, "_", dv_shortname, ".png") else paste0(filename_stem, ".png"); if (!is.null(plot_obj) && (inherits(plot_obj, "ggplot") || inherits(plot_obj, "gg"))) {tryCatch({ggsave(filename = file.path(path, filename), plot = plot_obj, width = width, height = height, dpi = 300, bg = "white"); cat(paste("Saved plot:", filename, "\n"))}, error = function(e) cat(paste("Error saving plot", filename, ":", e$message, "\n")))} else cat(paste("Skipping save for non-ggplot/NULL plot:", filename, "\n"))}
create_repel_data_smart <- function(df, x_var_sym, y_var_sym, label_var_sym = sym("msa_name"), size_var_sym = sym("population"), group_var_sym = NULL, n_overall_edges = 5, n_group_edges = 2, n_overall_size = 5) {if (nrow(df) == 0 || !as.character(x_var_sym) %in% names(df) || !as.character(y_var_sym) %in% names(df) || !as.character(size_var_sym) %in% names(df) ) return(tibble()); df_filtered <- df %>% filter(!is.na(!!x_var_sym) & !is.na(!!y_var_sym) & !is.na(!!size_var_sym)); if (nrow(df_filtered) == 0) return(tibble()); overall_largest <- df_filtered %>% arrange(desc(!!size_var_sym)) %>% slice_head(n = n_overall_size); edge_points_list <- list(df_filtered %>% arrange(!!x_var_sym) %>% slice_head(n = n_overall_edges), df_filtered %>% arrange(desc(!!x_var_sym)) %>% slice_head(n = n_overall_edges), df_filtered %>% arrange(!!y_var_sym) %>% slice_head(n = n_overall_edges), df_filtered %>% arrange(desc(!!y_var_sym)) %>% slice_head(n = n_overall_edges)); edge_points <- bind_rows(edge_points_list); group_edge_points <- tibble(); if (!is.null(group_var_sym) && as.character(group_var_sym) %in% names(df_filtered)) {group_splits <- df_filtered %>% group_by(!!group_var_sym) %>% group_split(); group_edge_points_list <- map(group_splits, function(g_df) {if (nrow(g_df) == 0) return(tibble()); bind_rows(g_df %>% arrange(!!x_var_sym) %>% slice_head(n = n_group_edges), g_df %>% arrange(desc(!!x_var_sym)) %>% slice_head(n = n_group_edges), g_df %>% arrange(!!y_var_sym) %>% slice_head(n = n_group_edges), g_df %>% arrange(desc(!!y_var_sym)) %>% slice_head(n = n_group_edges), g_df %>% arrange(desc(!!size_var_sym)) %>% slice_head(n = n_group_edges))}); group_edge_points <- bind_rows(group_edge_points_list)}; final_repel_df <- bind_rows(overall_largest, edge_points, group_edge_points) %>% distinct(!!label_var_sym, .keep_all = TRUE); return(final_repel_df)}

# =============================================================================
# PART 4: CROSS-SECTIONAL ANALYSIS (CURRENT YEAR)
# =============================================================================
if(nrow(msa_current_final) > 0 && !is.na(current_year_val)) {
  cat(paste0("\nPART 4: Cross-Sectional Analysis for Current Year (", current_year_val, ")...\n"))
  summary_by_size_current <- msa_current_final %>% filter(!is.na(size_category)) %>% group_by(size_category) %>% summarise(n_msas = n(), avg_pop_millions = weighted.mean(population / 1e6, w = population, na.rm = T), across(any_of(fertility_vars_shortnames), ~ weighted.mean(.x, w = population, na.rm = T), .names = "avg_{.col}"), avg_price_income = weighted.mean(price_to_income_ratio, w = population, na.rm = T), avg_rent_income = weighted.mean(rent_to_income_ratio, w = population, na.rm=T), avg_college_pct = weighted.mean(college_pct, w = population, na.rm = T), avg_pct_pop_under_5 = weighted.mean(pct_pop_under_5, w = population, na.rm =T), .groups = "drop") %>% arrange(size_category)
  print("--- Summary by MSA Size (Current Year) ---"); print(knitr::kable(summary_by_size_current, digits = 2)); write_csv(summary_by_size_current, file.path(analysis_plot_dir, "summary_by_size_current.csv"))
  for (dv_shortname in fertility_vars_shortnames) {
    if(!dv_shortname %in% names(msa_current_final) || sum(!is.na(msa_current_final[[dv_shortname]])) == 0) {cat(paste("DV", dv_shortname, "not found or all NA. Skipping.\n")); next}
    dv_config <- fertility_vars_config[[dv_shortname]]; dv_longname <- dv_config$longname; dv_plot_label_func <- dv_config$plot_label_func
    plot_data_current <- msa_current_final %>% filter(!is.na(.data[[dv_shortname]]) & is.finite(.data[[dv_shortname]]))
    if(nrow(plot_data_current) < 10) {cat(paste("Skipping plots for", dv_longname, "due to <10 valid obs.\n")); next}
    if("log_population" %in% names(plot_data_current) && "price_to_income_ratio" %in% names(plot_data_current)) {repel_dv_pop <- create_repel_data_smart(plot_data_current, sym("log_population"), sym(dv_shortname)); p_dv_vs_logpop <- ggplot(plot_data_current, aes(x = log_population, y = .data[[dv_shortname]])) + geom_point(aes(size = population, color = price_to_income_ratio), alpha = 0.6, na.rm = TRUE) + geom_smooth(method = "lm", aes(weight = population), color = "darkred", se = TRUE, na.rm = TRUE) + geom_text_repel(data = repel_dv_pop, aes(label = msa_name), size = 2.5, max.overlaps = Inf, na.rm = TRUE) + scale_size_continuous(name = "Population", labels = scales::comma_format(), range = c(1, 8)) + scale_color_viridis_c(name = "House Price/Income", option = "plasma", na.value = "grey80") + scale_y_continuous(labels = dv_plot_label_func) + labs(title = paste(dv_longname, "vs. City Size"), subtitle = paste("Year:", current_year_val), x = "Log Population", y = dv_longname); save_plot_obj(p_dv_vs_logpop, "fig_current_dv_vs_logpop", dv_shortname = dv_shortname)}
    if("price_to_income_ratio" %in% names(plot_data_current) && "size_category" %in% names(plot_data_current)) {repel_dv_housing <- create_repel_data_smart(plot_data_current, sym("price_to_income_ratio"), sym(dv_shortname)); p_dv_vs_housing <- ggplot(plot_data_current, aes(x = price_to_income_ratio, y = .data[[dv_shortname]])) + geom_point(aes(size = population, color = size_category), alpha = 0.6, na.rm = TRUE) + geom_smooth(method = "lm", aes(weight = population), color = "darkred", se = TRUE, na.rm = TRUE) + geom_smooth(method = "loess", aes(weight = population), color = "darkblue", se = FALSE, linetype = "dashed", span = 0.6, na.rm = TRUE) + geom_text_repel(data = repel_dv_housing, aes(label = msa_name), size = 2.5, max.overlaps = Inf, na.rm = TRUE) + scale_size_continuous(name = "Population", labels = scales::comma_format(), range = c(1, 8)) + scale_color_viridis_d(name = "City Size", direction = -1, na.translate = FALSE, drop=FALSE) + scale_y_continuous(labels = dv_plot_label_func) + labs(title = paste(dv_longname, "vs. Housing Affordability (Price)"), subtitle = paste("Year:", current_year_val), x = "House Price to Income Ratio", y = dv_longname); save_plot_obj(p_dv_vs_housing, "fig_current_dv_vs_price_housing", dv_shortname = dv_shortname)}
    if("rent_to_income_ratio" %in% names(plot_data_current) && "size_category" %in% names(plot_data_current)) {repel_dv_rent_housing <- create_repel_data_smart(plot_data_current, sym("rent_to_income_ratio"), sym(dv_shortname)); p_dv_vs_rent_housing <- ggplot(plot_data_current, aes(x = rent_to_income_ratio, y = .data[[dv_shortname]])) + geom_point(aes(size = population, color = size_category), alpha = 0.6, na.rm = TRUE) + geom_smooth(method = "lm", aes(weight = population), color = "darkred", se = TRUE, na.rm = TRUE) + geom_smooth(method = "loess", aes(weight = population), color = "darkblue", se = FALSE, linetype = "dashed", span = 0.6, na.rm = TRUE) + geom_text_repel(data = repel_dv_rent_housing, aes(label = msa_name), size = 2.5, max.overlaps = Inf, na.rm = TRUE) + scale_size_continuous(name = "Population", labels = scales::comma_format(), range = c(1, 8)) + scale_color_viridis_d(name = "City Size", direction = -1, na.translate = FALSE, drop=FALSE) + scale_y_continuous(labels = dv_plot_label_func) + labs(title = paste(dv_longname, "vs. Housing Affordability (Rent)"), subtitle = paste("Year:", current_year_val), x = "Rent to Income Ratio (Annual)", y = dv_longname); save_plot_obj(p_dv_vs_rent_housing, "fig_current_dv_vs_rent_housing", dv_shortname = dv_shortname)}
    if("college_pct" %in% names(plot_data_current) && "price_to_income_ratio" %in% names(plot_data_current)) {repel_dv_college <- create_repel_data_smart(plot_data_current, sym("college_pct"), sym(dv_shortname)); p_dv_vs_college <- ggplot(plot_data_current, aes(x = college_pct, y = .data[[dv_shortname]])) + geom_point(aes(size = population, color = price_to_income_ratio), alpha = 0.6, na.rm = TRUE) + geom_smooth(method = "lm", aes(weight = population), color = "darkred", se = TRUE, na.rm = TRUE) + geom_text_repel(data = repel_dv_college, aes(label = msa_name), size = 2.5, max.overlaps = Inf, na.rm = TRUE) + scale_size_continuous(name = "Population", labels = scales::comma_format(), range = c(1, 8)) + scale_color_viridis_c(name = "House Price/Income", option = "plasma", na.value = "grey80") + scale_x_continuous(labels = scales::percent_format(scale = 1, accuracy = 1)) + scale_y_continuous(labels = dv_plot_label_func) + labs(title = paste(dv_longname, "vs. Human Capital"), subtitle = paste("Year:", current_year_val), x = "% Adults with Bachelor's Degree or Higher", y = dv_longname); save_plot_obj(p_dv_vs_college, "fig_current_dv_vs_college", dv_shortname = dv_shortname)}
  }
  if (nrow(msa_current_final) > 10 && all(sapply(fertility_vars_shortnames, function(v) v %in% names(msa_current_final) && sum(!is.na(msa_current_final[[v]])) > 10 ))) {fert_measures_df_for_corr <- msa_current_final %>% select(all_of(fertility_vars_shortnames)) %>% filter(if_all(everything(), ~!is.na(.) & is.finite(.))); if(nrow(fert_measures_df_for_corr) > 10 && ncol(fert_measures_df_for_corr) >=2 ){p_fert_corr <- GGally::ggpairs(fert_measures_df_for_corr, title = paste("Correlation Matrix of Fertility Measures (Year:", current_year_val, ")"), upper = list(continuous = GGally::wrap("cor", size = 3, stars=FALSE, na.rm = TRUE)), lower = list(continuous = GGally::wrap("points", alpha = 0.5, size=0.8, na.rm = TRUE)), diag = list(continuous = GGally::wrap("densityDiag", alpha=0.5, na.rm = TRUE))) + theme(strip.text = element_text(size=8), axis.text = element_text(size=7, angle=0, hjust=0.5)); save_plot_obj(p_fert_corr, "fig_current_fertility_measures_corr", width = 10, height = 10)} else {cat("Skipping fertility correlation: insufficient valid data.\n")}} else {cat("Skipping fertility correlation: missing DVs or obs.\n")}
  cols_for_ols_check <- c("log_population", "price_to_income_ratio", "rent_to_income_ratio", "log_median_household_income", "college_pct", "population"); reg_data_current <- msa_current_final %>% filter(if_all(all_of(cols_for_ols_check), ~!is.na(.) & is.finite(.)))
  if (nrow(reg_data_current) >= 30) {ols_models_all_dvs <- list(); for (dv_shortname in fertility_vars_shortnames) {if(!dv_shortname %in% names(reg_data_current)) {cat(paste("DV", dv_shortname, "not for OLS.\n")); next}; dv_config <- fertility_vars_config[[dv_shortname]]; dv_longname <- dv_config$longname; current_dv_data_ols <- reg_data_current %>% filter(!is.na(.data[[dv_shortname]]) & is.finite(.data[[dv_shortname]])); if(nrow(current_dv_data_ols) < 25) {cat(paste("Skipping OLS for", dv_longname, "<25 obs.\n")); next}; models_ols <- list("LogPop Only" = lm(as.formula(paste0(dv_shortname, " ~ log_population")), data = current_dv_data_ols, weights = population), "+ Price/Inc" = lm(as.formula(paste0(dv_shortname, " ~ log_population + price_to_income_ratio")), data = current_dv_data_ols, weights = population), "+ Rent/Inc" = lm(as.formula(paste0(dv_shortname, " ~ log_population + price_to_income_ratio + rent_to_income_ratio")), data = current_dv_data_ols, weights = population), "+ College%" = lm(as.formula(paste0(dv_shortname, " ~ log_population + price_to_income_ratio + rent_to_income_ratio + college_pct")), data = current_dv_data_ols, weights = population), "Full Model" = lm(as.formula(paste0(dv_shortname, " ~ log_population + price_to_income_ratio + rent_to_income_ratio + college_pct + log_median_household_income")), data = current_dv_data_ols, weights = population)); ols_models_all_dvs[[dv_longname]] <- models_ols; coef_rename_map <- c("log_population"="Log Population", "price_to_income_ratio"="Price/Income Ratio", "rent_to_income_ratio"="Rent/Income Ratio", "college_pct"="% College Educated", "log_median_household_income"="Log Median HH Income"); msum_table_ols <- modelsummary(models_ols, stars = TRUE, gof_map = c("nobs", "r.squared"), title = paste("OLS Determinants of", dv_longname, "(Year:", current_year_val, ", Pop. Weighted)"), coef_map = coef_rename_map); print(msum_table_ols); modelsummary(models_ols, stars = TRUE, output = file.path(analysis_plot_dir, paste0("regtable_ols_", dv_shortname, ".html")), title = paste("OLS Determinants of", dv_longname, "(Year:", current_year_val, ", Pop. Weighted)"), coef_map = coef_rename_map)}} else { cat("Insufficient data for OLS (current year) after IV filtering.\n") }
} else { cat("msa_current_final empty or current_year_val not set. Skipping Part 4.\n") }

# =============================================================================
# PART 5: PANEL ANALYSIS (MULTIPLE YEARS)
# =============================================================================
if(nrow(msa_panel_final) > 0 && length(unique(msa_panel_final$year)) > 1 && !is.null(target_years_acs5) && length(target_years_acs5)>0) {
  cat("\nPART 5: Panel Analysis (Multiple Years)...\n")
  for (dv_shortname in fertility_vars_shortnames) {if(!dv_shortname %in% names(msa_panel_final) || sum(!is.na(msa_panel_final[[dv_shortname]])) == 0) {cat(paste("DV", dv_shortname, "not in panel or all NA. Skip panel plots.\n")); next}; dv_config <- fertility_vars_config[[dv_shortname]]; dv_longname <- dv_config$longname; dv_plot_label_func <- dv_config$plot_label_func; plot_data_trends <- msa_panel_final %>% filter(!is.na(.data[[dv_shortname]]) & !is.na(size_category) & !is.na(population) & is.finite(.data[[dv_shortname]])) %>% group_by(year, size_category) %>% summarise(avg_measure = weighted.mean(.data[[dv_shortname]], w = population, na.rm = TRUE), n_msas = n(), .groups = "drop") %>% filter(n_msas >= 3); if (nrow(plot_data_trends) < length(target_years_acs5) && length(unique(plot_data_trends$size_category)) < 2) {cat(paste("Skipping Time Trends for", dv_longname, "insufficient panel data.\n")); next}; p_trends <- ggplot(plot_data_trends, aes(x = year, y = avg_measure, color = size_category, group = size_category)) + geom_line(linewidth = 1.2, na.rm=TRUE) + geom_point(size = 2.5, na.rm=TRUE) + scale_y_continuous(labels = dv_plot_label_func) + scale_x_continuous(breaks = unique(msa_panel_final$year)) + scale_color_viridis_d(name = "MSA Size", direction = -1, na.translate=FALSE, drop=FALSE) + labs(title = paste("Time Trends in", dv_longname), subtitle = "Pop-weighted averages by MSA size category", x = "Year", y = paste("Avg.", dv_longname)); save_plot_obj(p_trends, "fig_panel_trends", dv_shortname = dv_shortname)}
  min_panel_year <- min(msa_panel_final$year, na.rm = TRUE); max_panel_year <- max(msa_panel_final$year, na.rm = TRUE); cols_for_changes_check <- c(fertility_vars_shortnames, "price_to_income_ratio", "rent_to_income_ratio", "population")
  msa_panel_changes <- msa_panel_final %>% filter(year %in% c(min_panel_year, max_panel_year)) %>% select(GEOID, msa_name, year, size_category, any_of(cols_for_changes_check)) %>% mutate(across(any_of(cols_for_changes_check), as.numeric)) %>% group_by(GEOID, msa_name) %>% filter(n() == 2 && all(sapply(cols_for_changes_check, function(v_col) { if(v_col %in% names(cur_data())) { sum(!is.na(cur_data()[[v_col]]) & is.finite(cur_data()[[v_col]])) == 2 } else { FALSE } } ))) %>% arrange(GEOID, year) %>% summarise(size_category_end = last(size_category), population_start = first(population), across(any_of(fertility_vars_shortnames), ~ (last(.x) - first(.x)), .names = "delta_{.col}"), delta_price_to_income = last(price_to_income_ratio) - first(price_to_income_ratio), delta_rent_to_income = last(rent_to_income_ratio) - first(rent_to_income_ratio), .groups = "drop") %>% filter(!is.na(population_start) & population_start > 0)
  if (nrow(msa_panel_changes) > 10) {for (dv_shortname in fertility_vars_shortnames) {delta_dv_colname <- paste0("delta_", dv_shortname); if(!delta_dv_colname %in% names(msa_panel_changes)) {cat(paste("Delta col for", dv_shortname, "not found. Skip change plot.\n")); next}; dv_config <- fertility_vars_config[[dv_shortname]]; dv_longname <- dv_config$longname; plot_data_delta <- msa_panel_changes %>% filter(!is.na(.data[[delta_dv_colname]]) & is.finite(.data[[delta_dv_colname]])); if(nrow(plot_data_delta) < 10 || !"delta_price_to_income" %in% names(plot_data_delta)) next; repel_delta_price <- create_repel_data_smart(plot_data_delta, sym("delta_price_to_income"), sym(delta_dv_colname), size_var_sym = sym("population_start")); p_delta_dv_vs_price <- ggplot(plot_data_delta, aes(x = delta_price_to_income, y = .data[[delta_dv_colname]])) + geom_point(aes(size = population_start, color = size_category_end), alpha = 0.7, na.rm = TRUE) + geom_smooth(method = "lm", aes(weight = population_start), color = "darkred", na.rm = TRUE) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + geom_text_repel(data = repel_delta_price, aes(label = msa_name), size = 2.5, max.overlaps = Inf, na.rm = TRUE) + scale_size_continuous(name = paste0("Pop. (", min_panel_year, ")"), labels = scales::comma, range = c(1,10)) + scale_color_viridis_d(name = paste0("Size Cat. (", max_panel_year, ")"), direction = -1, na.translate=FALSE, drop=FALSE) + labs(title = paste("Change in", dv_longname, "vs. Change in Price/Income"), subtitle = paste(min_panel_year, "to", max_panel_year), x = "Change in House Price to Income Ratio", y = paste("Change in", dv_longname)); save_plot_obj(p_delta_dv_vs_price, "fig_panel_delta_dv_vs_price", dv_shortname = dv_shortname)}} else { cat("Not enough MSAs for change analysis.\n") }
  cols_for_panel_reg_check <- c(fertility_vars_shortnames, "log_population", "price_to_income_ratio", "rent_to_income_ratio", "log_median_household_income", "college_pct", "population") # Added population for weights if needed, but plm usually doesn't use formula weights directly
  panel_reg_data_prep <- msa_panel_final %>% 
    select(GEOID, year, any_of(cols_for_panel_reg_check)) %>% 
    mutate(across(any_of(setdiff(cols_for_panel_reg_check, c("GEOID", "year"))), as.numeric)) %>% 
    group_by(GEOID) %>% 
    filter(n_distinct(year) >= 2) %>% # Keep MSAs observed in at least 2 distinct years
    # Ensure sufficient non-NA observations for key variables within each MSA
    filter(all(sapply(setdiff(cols_for_panel_reg_check, c("GEOID", "year", "population")), function(v_col) {
      if (v_col %in% names(cur_data())) { sum(!is.na(cur_data()[[v_col]]) & is.finite(cur_data()[[v_col]])) >= 2 } else { FALSE }
    }))) %>%
    ungroup()
  
  if (nrow(panel_reg_data_prep) > 0 && length(unique(panel_reg_data_prep$GEOID)) > 1 && length(unique(panel_reg_data_prep$year)) > 1) {
    panel_reg_data <- pdata.frame(panel_reg_data_prep, index = c("GEOID", "year"))
    if(nrow(panel_reg_data) > 50 && length(unique(panel_reg_data$GEOID)) > 20 && length(unique(panel_reg_data$year)) > 1) {fe_models_all_dvs <- list(); for (dv_shortname in fertility_vars_shortnames) {if(!dv_shortname %in% names(panel_reg_data)) {cat(paste("DV", dv_shortname, "not for Panel FE.\n")); next}; dv_config <- fertility_vars_config[[dv_shortname]]; dv_longname <- dv_config$longname; dv_variation_check <- panel_reg_data %>% filter(!is.na(.data[[dv_shortname]])) %>% group_by(GEOID) %>% summarise(sd_dv = sd(.data[[dv_shortname]], na.rm=TRUE), n_years_dv = n(), .groups="drop") %>% filter(!is.na(sd_dv) & sd_dv > 1e-6 & n_years_dv >=2); if(nrow(dv_variation_check) < 20){cat(paste("Skipping Panel FE for", dv_longname, "insufficient variation/obs.\n")); next}; fe_formula_str <- paste0(dv_shortname, " ~ price_to_income_ratio + rent_to_income_ratio + log_median_household_income + college_pct + log_population"); panel_data_for_model <- panel_reg_data %>% filter(GEOID %in% dv_variation_check$GEOID) %>% filter(!is.na(.data[[dv_shortname]])); if(nrow(panel_data_for_model) < 50) {cat(paste("Skipping Panel FE for", dv_longname, "insufficient obs after DV variation filter.\n")); next}; m_fe_msa <- tryCatch(plm(as.formula(fe_formula_str), data = panel_data_for_model, model = "within", effect = "individual"), error=function(e) NULL); m_fe_twoways <- tryCatch(plm(as.formula(fe_formula_str), data = panel_data_for_model, model = "within", effect = "twoways"), error=function(e) NULL); current_fe_models <- list(); if(!is.null(m_fe_msa)) current_fe_models[["MSA FE"]] <- m_fe_msa; if(!is.null(m_fe_twoways)) current_fe_models[["MSA & Year FE"]] <- m_fe_twoways; if(length(current_fe_models) > 0){fe_models_all_dvs[[dv_longname]] <- current_fe_models; coef_rename_map_fe <- c("price_to_income_ratio"="Price/Income Ratio", "rent_to_income_ratio"="Rent/Income Ratio", "log_median_household_income"="Log Median HH Income", "college_pct"="% College Educated", "log_population"="Log Population"); msum_table_fe <- modelsummary(current_fe_models, stars = TRUE, gof_map = "nobs.+", title = paste("Panel FE Determinants of", dv_longname), coef_map = coef_rename_map_fe); print(msum_table_fe); modelsummary(current_fe_models, stars = TRUE, output = file.path(analysis_plot_dir, paste0("regtable_panel_fe_", dv_shortname, ".html")), title = paste("Panel FE Determinants of", dv_longname), coef_map = coef_rename_map_fe)}}} else { cat("Insufficient data for Panel FE after initial filtering.\n") }
  } else { cat("Not enough valid MSAs or years in panel_reg_data_prep for pdata.frame conversion or modeling.\n") }
} else { cat("msa_panel_final empty or only one year. Skipping Part 5.\n") }

cat("\n--- CONSOLIDATED ANALYSIS COMPLETE ---\n")
cat(paste0("Results saved in/printed from: '", analysis_plot_dir, "'\n"))



# =============================================================================
# Quick script to re-fetch just the birth rate data and update your existing panel
# This is faster than re-running the entire analysis
# NOTE: The main script above now has robust S1301 handling. 
# This "quick script" uses S1301_C04_001 for ALL years, which is
# INCORRECT for ACS data prior to 2022. For pre-2022 ACS data,
# S1301_C03_001 was generally the correct variable for overall birth rate.
# The main analysis script's `process_recent_births_S1301_final` function handles this correctly.
# This section is retained from your original prompt but might need adjustment or
# may be superseded by running the full script above.
# =============================================================================
cat("\n--- START OF 'QUICK SCRIPT' SECTION (from original prompt) ---\n")

library(tidycensus) # Already loaded, but good for standalone context
library(tidyverse)   # Already loaded

# Load existing data (assuming it was created by a previous run of the main script)
# Make sure this path is correct if you ran the extended script above.
# It might be in "consolidated_fertility_analysis_v5_extended"
# For this example, let's assume the original path for demonstration
original_panel_path <- "consolidated_fertility_analysis_v5/msa_panel_final.rds"
if (file.exists(original_panel_path)) {
  msa_panel_final_original <- readRDS(original_panel_path)
  cat(paste("Loaded original panel data from:", original_panel_path, "\n"))
  
  # Function to fetch birth rate data for a specific year (from your quick script)
  fetch_birth_rate_fixed_quickscript <- function(year_val) {
    cat(paste0("\n[QuickScript] Fetching birth rate for year ", year_val, " using S1301_C04_001...\n"))
    cat("WARNING: [QuickScript] Using S1301_C04_001 for all years. This is incorrect for pre-2022 data.\n")
    
    tryCatch({
      s1301_data <- get_acs(
        geography = "metropolitan statistical area/micropolitan statistical area",
        table = "S1301",
        year = year_val,
        survey = "acs5"
      )
      
      # Extract birth rate from S1301_C04_001 (This is the part that's problematic for older years)
      birth_rates <- s1301_data %>%
        filter(variable == "S1301_C04_001") %>% # <<<<< Problematic for pre-2022
        transmute(
          GEOID,
          year = year_val,
          birth_rate_15_50_new_quickscript = as.numeric(estimate)
        )
      
      mean_rate <- mean(birth_rates$birth_rate_15_50_new_quickscript, na.rm = TRUE)
      cat(paste0("[QuickScript] Mean birth rate (using S1301_C04_001): ", round(mean_rate, 1), 
                 " (n=", sum(!is.na(birth_rates$birth_rate_15_50_new_quickscript)), ")\n"))
      
      return(birth_rates)
      
    }, error = function(e) {
      cat(paste0("[QuickScript] Error: ", e$message, "\n"))
      return(NULL)
    })
  }
  
  # Fetch birth rates for years to "fix" according to quick script
  years_to_fix_quickscript <- c(2014, 2016, 2018) # 2022 already "works" with S1301_C04_001
  birth_rate_updates_quickscript <- map_dfr(years_to_fix_quickscript, fetch_birth_rate_fixed_quickscript)
  
  if (!is.null(birth_rate_updates_quickscript) && nrow(birth_rate_updates_quickscript) > 0) {
    
    msa_panel_fixed_quickscript <- msa_panel_final_original %>%
      left_join(
        birth_rate_updates_quickscript %>% select(GEOID, year, birth_rate_15_50_new_quickscript),
        by = c("GEOID", "year")
      ) %>%
      mutate(
        # This logic replaces birth_rate_15_50 with the new one if available
        birth_rate_15_50_qs_updated = coalesce(birth_rate_15_50_new_quickscript, birth_rate_15_50)
      ) # Not removing the original birth_rate_15_50 for comparison here
    
    cat("\n=== [QuickScript] VERIFICATION ===\n")
    birth_rate_check_qs <- msa_panel_fixed_quickscript %>%
      group_by(year) %>%
      summarise(
        n_total = n(),
        n_with_original_br = sum(!is.na(birth_rate_15_50)),
        mean_original_br = round(mean(birth_rate_15_50, na.rm = TRUE), 1),
        n_with_qs_updated_br = sum(!is.na(birth_rate_15_50_qs_updated)),
        mean_qs_updated_br = round(mean(birth_rate_15_50_qs_updated, na.rm = TRUE), 1),
        .groups = "drop"
      )
    
    print(birth_rate_check_qs)
    
    # Save the fixed data from quick script
    quick_script_output_dir <- "consolidated_fertility_analysis_quickscript_fixed"
    dir.create(quick_script_output_dir, showWarnings = FALSE)
    saveRDS(msa_panel_fixed_quickscript, file.path(quick_script_output_dir, "msa_panel_quickscript_fixed.rds"))
    cat(paste("\n[QuickScript] Data 'fixed' by quick script saved to '", quick_script_output_dir, "' directory\n"))
    
    # Example plot from quick script logic
    cat("\n=== [QuickScript] RE-RUNNING PANEL PLOTS (using quick script logic) ===\n")
    output_dir_qs <- file.path(quick_script_output_dir, "plots_and_tables")
    dir.create(output_dir_qs, showWarnings = FALSE, recursive = TRUE)
    
    plot_data_trends_qs <- msa_panel_fixed_quickscript %>%
      filter(!is.na(birth_rate_15_50_qs_updated) & !is.na(size_category) & !is.na(population)) %>%
      group_by(year, size_category) %>%
      summarise(
        avg_measure = weighted.mean(birth_rate_15_50_qs_updated, w = population, na.rm = TRUE),
        n_msas = n(),
        .groups = "drop"
      ) %>%
      filter(n_msas >= 3)
    
    if (nrow(plot_data_trends_qs) > 0) {
      p_trends_birth_rate_qs <- ggplot(plot_data_trends_qs, 
                                       aes(x = year, y = avg_measure, 
                                           color = size_category, group = size_category)) +
        geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
        scale_y_continuous(labels = scales::comma_format(accuracy = 1)) +
        scale_x_continuous(breaks = unique(msa_panel_fixed_quickscript$year)) +
        scale_color_viridis_d(name = "MSA Size", direction = -1) +
        labs(title = "[QuickScript] Time Trends in Birth Rate (S1301_C04_001 used for some years)",
             subtitle = "Pop-weighted averages by MSA size category",
             x = "Year", y = "Avg. Birth Rate (per 1k Women 15-50)") +
        theme_minimal() + theme(legend.position = "bottom")
      
      ggsave(filename = file.path(output_dir_qs, "fig_panel_trends_birth_rate_QUICKSCRIPT_FIXED.png"),
             plot = p_trends_birth_rate_qs, width = 10, height = 8, dpi = 300, bg = "white")
      cat("[QuickScript] New birth rate trends plot (from quick script) saved!\n")
    } else {
      cat("[QuickScript] Not enough data to generate trends plot.\n")
    }
    
  } else {
    cat("[QuickScript] Failed to fetch birth rate updates using quick script logic.\n")
  }
} else {
  cat(paste("[QuickScript] Original panel data not found at:", original_panel_path, ". Skipping quick script section.\n"))
}
cat("--- END OF 'QUICK SCRIPT' SECTION ---\n")





# =============================================================================
# PART 6: TIME-VARYING COEFFICIENT ANALYSIS (FERTILITY VS. LOG POPULATION) - WEIGHTED & UNWEIGHTED
# =============================================================================
cat("\nPART 6: Time-Varying Coefficient Analysis (Fertility vs. Log Population) - Weighted & Unweighted...\n")

# Ensure msa_panel_final exists and has multiple years, and fertility_vars_shortnames is defined
if (exists("msa_panel_final") && nrow(msa_panel_final) > 0 && 
    length(unique(msa_panel_final$year)) > 1 &&
    exists("fertility_vars_shortnames") && length(fertility_vars_shortnames) > 0) {
  
  # Function to run yearly regression and extract coefficient for log_population
  run_yearly_regression_for_logpop_effect <- function(df_for_year, dv_name_str, is_weighted = TRUE) {
    
    required_cols <- c(dv_name_str, "log_population", "year")
    if (is_weighted) required_cols <- c(required_cols, "population")
    
    if (!all(required_cols %in% names(df_for_year))) {
      return(tibble(
        dv_shortname = dv_name_str,
        year = if("year" %in% names(df_for_year)) unique(df_for_year$year)[1] else NA_integer_,
        coefficient_log_pop = NA_real_,
        std_error_log_pop = NA_real_,
        p_value_log_pop = NA_real_,
        n_obs_reg = NA_integer_,
        weighted = is_weighted
      ))
    }
    
    reg_data_year_dv <- df_for_year %>%
      filter(
        !is.na(.data[[dv_name_str]]) & is.finite(.data[[dv_name_str]]),
        !is.na(log_population) & is.finite(log_population)
      )
    if (is_weighted) {
      reg_data_year_dv <- reg_data_year_dv %>% filter(!is.na(population) & population > 0)
    }
    
    current_year_val_reg <- unique(reg_data_year_dv$year)[1]
    
    if (nrow(reg_data_year_dv) < 20) { 
      return(tibble(
        dv_shortname = dv_name_str,
        year = current_year_val_reg,
        coefficient_log_pop = NA_real_,
        std_error_log_pop = NA_real_,
        p_value_log_pop = NA_real_,
        n_obs_reg = nrow(reg_data_year_dv),
        weighted = is_weighted
      ))
    }
    
    formula_obj <- as.formula(paste0(dv_name_str, " ~ log_population"))
    
    tryCatch({
      if (is_weighted) {
        model_fit <- lm(formula_obj, data = reg_data_year_dv, weights = population)
      } else {
        model_fit <- lm(formula_obj, data = reg_data_year_dv) # Unweighted
      }
      summary_fit <- summary(model_fit)
      
      coef_info <- NA
      if ("log_population" %in% rownames(summary_fit$coefficients)) {
        coef_info <- summary_fit$coefficients["log_population", ]
        results <- tibble(
          dv_shortname = dv_name_str,
          year = current_year_val_reg,
          coefficient_log_pop = coef_info["Estimate"],
          std_error_log_pop = coef_info["Std. Error"],
          p_value_log_pop = coef_info["Pr(>|t|)"],
          n_obs_reg = nobs(model_fit),
          weighted = is_weighted
        )
      } else {
        results <- tibble(
          dv_shortname = dv_name_str,
          year = current_year_val_reg,
          coefficient_log_pop = NA_real_,
          std_error_log_pop = NA_real_,
          p_value_log_pop = NA_real_,
          n_obs_reg = nobs(model_fit),
          weighted = is_weighted
        )
      }
      return(results)
    }, error = function(e) {
      cat(paste0("Error running ", if(is_weighted) "weighted" else "unweighted", " regression for DV '", dv_name_str, "' in year ", current_year_val_reg, ": ", e$message, "\n"))
      return(tibble(
        dv_shortname = dv_name_str,
        year = current_year_val_reg,
        coefficient_log_pop = NA_real_,
        std_error_log_pop = NA_real_,
        p_value_log_pop = NA_real_,
        n_obs_reg = nrow(reg_data_year_dv),
        weighted = is_weighted
      ))
    })
  }
  
  all_panel_years <- sort(unique(msa_panel_final$year))
  time_varying_coeffs_list_all <- list()
  
  for (current_dv in fertility_vars_shortnames) {
    if (!current_dv %in% names(msa_panel_final)) {
      cat(paste0("Skipping time-varying coefficient calculation for DV: '", current_dv, "' as it's not in msa_panel_final.\n"))
      next
    }
    cat(paste0("\nCalculating time-varying coefficients of log_population for DV: '", current_dv, "' (Weighted & Unweighted)...\n"))
    
    # Weighted regressions
    coeffs_for_current_dv_weighted <- map_dfr(all_panel_years, function(yr) {
      data_this_year <- msa_panel_final %>% filter(year == yr)
      run_yearly_regression_for_logpop_effect(df_for_year = data_this_year, dv_name_str = current_dv, is_weighted = TRUE)
    })
    time_varying_coeffs_list_all[[paste0(current_dv, "_weighted")]] <- coeffs_for_current_dv_weighted
    
    # Unweighted regressions
    coeffs_for_current_dv_unweighted <- map_dfr(all_panel_years, function(yr) {
      data_this_year <- msa_panel_final %>% filter(year == yr)
      run_yearly_regression_for_logpop_effect(df_for_year = data_this_year, dv_name_str = current_dv, is_weighted = FALSE)
    })
    time_varying_coeffs_list_all[[paste0(current_dv, "_unweighted")]] <- coeffs_for_current_dv_unweighted
  }
  
  all_time_varying_coeffs_df <- bind_rows(time_varying_coeffs_list_all, .id = "dv_weight_id") %>%
    separate(dv_weight_id, into = c("original_dv_shortname", "weight_status"), sep = "_", extra = "merge", fill = "right") %>%
    mutate(weight_status = ifelse(is.na(weight_status) | weight_status != "unweighted", "weighted", "unweighted"))
  
  
  if (nrow(all_time_varying_coeffs_df) > 0 && any(!is.na(all_time_varying_coeffs_df$coefficient_log_pop))) {
    write_csv(all_time_varying_coeffs_df, file.path(analysis_plot_dir, "time_varying_coefficients_log_population_weighted_unweighted.csv"))
    cat(paste0("Saved time-varying coefficients table to: ", file.path(analysis_plot_dir, "time_varying_coefficients_log_population_weighted_unweighted.csv"), "\n"))
    
    for (dv_to_plot in fertility_vars_shortnames) {
      for (is_w in c(TRUE, FALSE)) {
        
        weight_label <- if(is_w) "weighted" else "unweighted"
        data_for_plot <- all_time_varying_coeffs_df %>%
          filter(original_dv_shortname == dv_to_plot, weighted == is_w, !is.na(coefficient_log_pop))
        
        if(nrow(data_for_plot) < 2) {
          cat(paste0("Skipping plot for ", weight_label, " time-varying coefficient (DV: ", dv_to_plot, ") - insufficient data (", nrow(data_for_plot), ").\n"))
          next
        }
        
        dv_config_plot <- fertility_vars_config[[dv_to_plot]]
        dv_longname_plot <- dv_config_plot$longname
        plot_subtitle_weight_text <- if(is_w) "MSA population weighted" else "Unweighted"
        
        plot_data_with_ci <- data_for_plot %>%
          mutate(
            ci_lower = coefficient_log_pop - 1.96 * std_error_log_pop,
            ci_upper = coefficient_log_pop + 1.96 * std_error_log_pop
          )
        
        p_time_coeff <- ggplot(plot_data_with_ci, aes(x = year, y = coefficient_log_pop)) +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "dodgerblue") +
          geom_line(color = "dodgerblue", linewidth = 1.1) +
          geom_point(color = "dodgerblue", size = 2.5) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          scale_x_continuous(breaks = unique(plot_data_with_ci$year)) +
          labs(
            title = paste("Time-Varying Effect of Log Population on", dv_longname_plot),
            subtitle = paste0("Coefficient from annual OLS (DV ~ log_population, ", plot_subtitle_weight_text, ")\nError bands represent 95% CI"),
            x = "Year",
            y = "Coefficient of Log Population"
          ) +
          theme_minimal(base_size = 13) + 
          theme(
            plot.title = element_text(size = 15),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
          )
        
        save_plot_obj(p_time_coeff, 
                      filename_stem = paste0("fig_time_coeff_logpop_", weight_label), 
                      dv_shortname = dv_to_plot, 
                      width = 9, height = 7)
        print(p_time_coeff)
      }
    }
  } else {
    cat("No time-varying coefficients were successfully calculated to plot or save.\n")
  }
} else {
  cat("Skipping PART 6 (Time-Varying Coefficients): msa_panel_final is empty, has only one year, or fertility_vars_shortnames not defined.\n")
}