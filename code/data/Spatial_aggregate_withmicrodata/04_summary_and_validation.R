# ==============================================================================
# 04_summary_and_validation.R
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
factor_vars <- c("educ_factor", "migration_factor", "nativity_factor", "is_recent_mother")

# Summary for numeric variables - a more robust, direct approach
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
    # Create the table and pipe to kable
    freq_table <- micro_data[, .N, by = var_name]
    setnames(freq_table, "N", "Count")
    print(kable(freq_table, format = "pipe", caption = paste("Distribution of", var_name)))
    cat("\n")
}

# --- Distribution for Categorical Variables ---
cat("Distribution for: educ_factor \n\n")
print(knitr::kable(as.data.frame(table(micro_data$educ_factor)), col.names = c("educ_factor", "Count"), caption = "Distribution of educ_factor"))

cat("\nDistribution for: migration_factor \n\n")
print(knitr::kable(as.data.frame(table(micro_data$migration_factor)), col.names = c("migration_factor", "Count"), caption = "Distribution of migration_factor"))

# cat("\nDistribution for: nativity_factor \n\n")
# print(knitr::kable(as.data.frame(table(micro_data$nativity_factor)), col.names = c("nativity_factor", "Count"), caption = "Distribution of nativity_factor"))

# Save the final dataset
# ------------------------------------------------------------------------------
cat("\n\n--- Saving a copy of the validated sample ---\n")
saveRDS(micro_data, file = file.path(OUTPUT_DIR, "validated_microdata_sample.rds"))

# --- 3. AGGREGATE VALIDATION BY MSA (met2013) ---
cat("--- Aggregating Microdata by MSA for Validation ---\n")

# Ensure perwt is numeric for weighting
micro_data[, perwt := as.numeric(perwt)]

# Aggregate the microdata to the MSA level
# Here, we calculate the key metrics we want to validate
msa_aggregated_data <- micro_data[!is.na(met2013) & met2013 != 0, .(
    
    # Weighted mean age at first birth from microdata
    mean_age_at_first_birth_micro = weighted.mean(age_at_first_birth, w = perwt, na.rm = TRUE),
    
    # Weighted mean for housing variables
    mean_rent_micro = weighted.mean(rent, w = perwt, na.rm = TRUE),
    mean_valueh_micro = weighted.mean(valueh, w = perwt, na.rm = TRUE),
    
    # General Fertility Rate (GFR) from microdata
    # GFR = (Number of births / Number of women aged 15-44) * 1000
    total_women_15_44 = sum(perwt),
    total_births = sum(is_recent_mother * perwt, na.rm = TRUE),
    
    # Weighted mean number of children
    mean_nchild_micro = weighted.mean(nchild, w = perwt, na.rm = TRUE)
    
), by = .(met2013, year)]

# Calculate the GFR
msa_aggregated_data[, gfr_micro := (total_births / total_women_15_44) * 1000]

cat("Microdata aggregated to MSA-year level.\n")
cat("Sample of aggregated data:\n")
print(head(msa_aggregated_data))
cat("\n")


# --- 4. COMPARISON WITH ORIGINAL AGGREGATE DATA ---
# In a full analysis, we would load the results from the `01_get_msa_panel.R` 
# script and merge it with `msa_aggregated_data` on `met2013` and `year`
# to perform a direct comparison or plot them against each other.

# For now, we will save our aggregated results and print a summary,
# as the direct comparison requires running the other pipeline.

cat("--- Summary of Aggregated Microdata Results (Averages across all MSAs) ---\n")

# Calculate and print the overall average of our new metrics across all MSAs
overall_summary <- msa_aggregated_data[, .(
    avg_mean_age_at_first_birth = mean(mean_age_at_first_birth_micro, na.rm = TRUE),
    avg_rent = mean(mean_rent_micro, na.rm = TRUE),
    avg_valueh = mean(mean_valueh_micro, na.rm = TRUE),
    avg_gfr = mean(gfr_micro, na.rm = TRUE),
    avg_nchild = mean(mean_nchild_micro, na.rm = TRUE)
)]

cat("Overall Averages Across All MSAs from Microdata:\n")
print(overall_summary)
cat("\n")

# Save the aggregated data for future use or direct comparison
AGGREGATED_OUTPUT_PATH <- file.path(OUTPUT_DIR, "microdata_msa_year_aggregates.rds")
saveRDS(msa_aggregated_data, AGGREGATED_OUTPUT_PATH)
cat(paste("Aggregated MSA-level data from microdata saved to:", AGGREGATED_OUTPUT_PATH, "\n"))

cat("\n--- Script 04 Finished ---\n")
cat("Next steps would be to load the original aggregate data and merge it with our new results for a direct, row-by-row validation.\n") 