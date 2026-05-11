# ==============================================================================
# 06_export_to_stata.R
#
# Description:
# This script loads the final, cleaned 5% stratified sample dataset and
# exports it to a Stata (.dta) file for external inspection and analysis.
# ==============================================================================

# --- 1. SETUP & CONFIGURATION ---
cat("Setting up the environment...\n")
suppressPackageStartupMessages({
    library(here)
    library(data.table)
    library(haven)
})

# Define file paths
CLEAN_RDS_PATH <- here("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
OUTPUT_DTA_PATH <- here("Spatial_aggregate_withmicrodata", "output", "fertility_microdata_clean_5pct_sample.dta")

# --- 2. LOAD & EXPORT ---
cat("Loading the cleaned RDS dataset...\n")
if (!file.exists(CLEAN_RDS_PATH)) {
    stop("Cleaned RDS file not found. Please run scripts 03 and 04 first.")
}
final_data <- readRDS(CLEAN_RDS_PATH)
cat(paste("Loaded", nrow(final_data), "records.\n"))

cat(paste("Exporting data to Stata file at:", OUTPUT_DTA_PATH, "\n"))
haven::write_dta(final_data, OUTPUT_DTA_PATH)

cat("Export complete.\n")
cat("--- Script 06 Finished ---\n") 