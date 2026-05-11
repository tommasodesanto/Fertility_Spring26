# ==============================================================================
# 05_final_validation.R
#
# Description:
# This script performs two main tasks:
# 1. It generates detailed summary statistics for the cleaned, 5% stratified
#    sample of the IPUMS microdata.
# 2. It validates the microdata by replicating key aggregate fertility metrics
#    (e.g., mean age at first birth) by city (MSA) and comparing them to
#    the results from the original aggregate analysis.
# ==============================================================================

# --- 1. SETUP & CONFIGURATION ---
cat("Setting up the environment...\n")
suppressPackageStartupMessages({
    library(data.table)
    library(here)
    library(dplyr)
    library(knitr)
})

# Define file paths
CLEAN_MICRODATA_PATH <- here("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
OUTPUT_DIR <- here("Spatial_aggregate_withmicrodata", "output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load the cleaned microdata
cat(paste("Loading cleaned microdata from:", CLEAN_MICRODATA_PATH, "\n"))
if (!file.exists(CLEAN_MICRODATA_PATH)) {
    stop("Cleaned microdata file not found. Please run script 03 first.")
}
micro_data <- readRDS(CLEAN_MICRODATA_PATH)
setDT(micro_data)
cat(paste("Successfully loaded", nrow(micro_data), "records.\n\n"))


# --- 2. GENERAL SUMMARY STATISTICS ---
cat("--- Generating General Summary Statistics ---\n")

# Select key numeric and factor variables for summary
numeric_vars <- c("age", "nchild", "age_at_first_birth", "rent", "valueh")
# 'nativity_factor' has been removed as it is not in the dataset
factor_vars <- c("educ_factor", "migration_factor", "is_recent_mother")

# Summary for numeric variables
summary_stats <- micro_data[, .(
    Variable = numeric_vars,
    Min = sapply(.SD, min, na.rm = TRUE),
    Q1 = sapply(.SD, quantile, probs = 0.25, na.rm = TRUE),
    Median = sapply(.SD, median, na.rm = TRUE),
    Mean = sapply(.SD, mean, na.rm = TRUE),
    Q3 = sapply(.SD, quantile, probs = 0.75, na.rm = TRUE),
    Max = sapply(.SD, max, na.rm = TRUE),
    NAs = sapply(.SD, function(x) sum(is.na(x)))
), .SDcols = numeric_vars]

cat("Summary for key numeric variables:\n")
print(kable(summary_stats, format = "pipe", caption = "Numeric Variable Summary"))
cat("\n")

# Summary for factor variables (as tables)
for (var_name in factor_vars) {
    cat(paste("Distribution for:", var_name, "\n"))
    freq_table <- micro_data[, .N, by = var_name]
    setnames(freq_table, "N", "Count")
    print(kable(freq_table, format = "pipe", caption = paste("Distribution of", var_name)))
    cat("\n")
}

# Save the final dataset
# ------------------------------------------------------------------------------
cat("\n\n--- Saving a copy of the validated sample ---\n")
saveRDS(micro_data, file = file.path(OUTPUT_DIR, "validated_microdata_sample.rds"))

# --- 3. AGGREGATE VALIDATION BY MSA (met2013) ---
cat("--- Aggregating Microdata by MSA for Validation ---\n")

# Ensure perwt is numeric for weighting
micro_data[, perwt := as.numeric(perwt)]
micro_data[, hhwt := as.numeric(hhwt)]

# Aggregate the microdata to the MSA level
msa_aggregated_data <- micro_data[!is.na(met2013) & met2013 != 0, .(
    
    mean_age_at_first_birth_micro = weighted.mean(age_at_first_birth, w = perwt, na.rm = TRUE),
    
    # Calculate median income. Using person-level total income as a proxy.
    median_income_micro = Hmisc::wtd.quantile(inctot, weights = perwt, probs = 0.5, na.rm = TRUE),
    
    # Calculate percentage of population that is female
    pct_female_micro = sum(perwt[sex == 2], na.rm=TRUE) / sum(perwt, na.rm=TRUE),
    
    mean_rent_micro = weighted.mean(rent, w = perwt, na.rm = TRUE),
    mean_valueh_micro = weighted.mean(valueh, w = perwt, na.rm = TRUE),
    total_women_15_44 = sum(perwt),
    total_births = sum(is_recent_mother * perwt, na.rm = TRUE),
    mean_nchild_micro = weighted.mean(nchild, w = perwt, na.rm = TRUE)
    
), by = .(met2013, year)]

# Calculate the GFR
msa_aggregated_data[, gfr_micro := (total_births / total_women_15_44) * 1000]

cat("Microdata aggregated to MSA-year level.\n")
cat("Sample of aggregated data:\n")
print(head(msa_aggregated_data))
cat("\n")

# Save the aggregated data before the comparison
AGGREGATED_OUTPUT_PATH <- file.path(OUTPUT_DIR, "microdata_msa_year_aggregates.rds")
saveRDS(msa_aggregated_data, AGGREGATED_OUTPUT_PATH)
cat(paste("Aggregated MSA-level data from microdata saved to:", AGGREGATED_OUTPUT_PATH, "\n"))


# --- 4. COMPARISON WITH ORIGINAL AGGREGATE DATA ---
cat("\n--- Loading Original Aggregate Data for Comparison ---\n")

# Define path to original aggregate data
ORIG_AGG_PATH <- here("Spatial_Aggregate_Analysis", "data_raw", "msa_panel_complete.rds")
if (!file.exists(ORIG_AGG_PATH)) {
    stop("Original aggregate data file not found at:", ORIG_AGG_PATH)
}
orig_agg_data <- readRDS(ORIG_AGG_PATH)
setDT(orig_agg_data)
cat("Successfully loaded original aggregate data.\n")
cat("\n--- Columns in Original Aggregate Data ---\n")
message(paste(names(orig_agg_data), collapse = ", "))
cat("----------------------------------------\n\n")

# Prepare for merge: standardize key column names and types
setnames(orig_agg_data, "GEOID", "met2013", skip_absent = TRUE)
setnames(orig_agg_data, "YEAR", "year", skip_absent = TRUE)
orig_agg_data[, met2013 := as.numeric(as.character(met2013))]
orig_agg_data[, year := as.numeric(as.character(year))]

# Calculate pct_female for the original aggregate data
if ("total_population" %in% names(orig_agg_data) && "B01001_002E" %in% names(orig_agg_data)) {
    orig_agg_data[, pct_female_orig := (total_population - B01001_002E) / total_population]
    cat("Calculated 'pct_female_orig' for the original aggregate data.\n")
}

# Merge the two datasets
validation_data <- merge(
    msa_aggregated_data,
    orig_agg_data,
    by = c("met2013", "year"),
    suffixes = c("_micro", "_orig")
)
cat(paste("Successfully merged", nrow(validation_data), "common MSA-year records for validation.\n\n"))

if (nrow(validation_data) == 0) {
    cat("Warning: No common MSA-year observations found between microdata and aggregate data. Cannot perform validation.\n")
} else {
    # --- 5. VISUAL & STATISTICAL VALIDATION ---
    cat("--- Performing Visual and Statistical Validation ---\n")
    
    # Define the pairs of variables to compare
    # We will look for variables that exist in both datasets after the merge
    comparison_vars <- list(
      "mean_age_at_first_birth_micro" = "mean_age_at_birth",
      "median_income_micro" = "median_hh_incomeE",
      "pct_female_micro" = "pct_female_orig" 
    )
    
    # Add ggplot2 for plotting
    if (!require("ggplot2", character.only = TRUE)) install.packages("ggplot2")
    library(ggplot2)

    # Add Hmisc for weighted median calculation
    if (!require("Hmisc", character.only = TRUE)) install.packages("Hmisc")
    library(Hmisc)

    for (micro_col in names(comparison_vars)) {
        orig_col <- comparison_vars[[micro_col]]
        
        if (!orig_col %in% names(validation_data)) {
            cat(paste("\nNOTE: Column '", orig_col, "' not found in original aggregate data. Skipping comparison for this variable.\n"))
            next
        }
        
        cat(paste("\n--- Validation for:", micro_col, "vs.", orig_col, "---\n"))
        
        temp_comp_data <- na.omit(validation_data[, .(m = .SD[[micro_col]], o = .SD[[orig_col]])])
        
        if (nrow(temp_comp_data) < 2) {
            cat("Not enough common data points to compare.\n")
            next
        }
        
        correlation <- cor(temp_comp_data$m, temp_comp_data$o)
        cat(paste("Correlation:", round(correlation, 4), "\n"))
        
        plot_title <- paste("Validation:", micro_col, "vs.", orig_col)
        plot_path <- file.path(OUTPUT_DIR, paste0("validation_plot_", micro_col, ".png"))
        
        p <- ggplot(temp_comp_data, aes(x = o, y = m)) +
            geom_point(alpha = 0.5) +
            geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
            labs(title = plot_title, x = paste("Original Aggregate:", orig_col), y = paste("Microdata Sample:", micro_col)) +
            theme_bw()
        
        ggsave(plot_path, p, width = 8, height = 6)
        cat(paste("Validation plot saved to:", plot_path, "\n"))
    }
}

cat("\n--- Script 05 Finished ---\n")
cat("Next steps would be to load the original aggregate data and merge it with our new results for a direct, row-by-row validation.\n") 