# ==============================================================================
# 07_export_orig_agg_to_stata.R
#
# Description:
# This script loads the original, aggregate MSA panel data (from the
# `Spatial_aggregate_analysis` project) and exports it to a Stata (.dta)
# file for easier external inspection.
# ==============================================================================

# --- 1. SETUP & CONFIGURATION ---
cat("Setting up the environment...\n")
suppressPackageStartupMessages({
    library(here)
    library(haven)
})

# Define file paths
ORIGINAL_RDS_PATH <- here("Spatial_aggregate_analysis", "data_raw", "msa_panel_analysis_ready.rds")
OUTPUT_DTA_PATH <- here("Spatial_aggregate_withmicrodata", "output", "original_msa_panel_for_inspection.dta")

# --- 2. LOAD & EXPORT ---
cat("Loading the original aggregate RDS dataset...\n")
if (!file.exists(ORIGINAL_RDS_PATH)) {
    stop("Original aggregate RDS file not found.")
}
original_data <- readRDS(ORIGINAL_RDS_PATH)
cat(paste("Loaded", nrow(original_data), "records.\n"))

cat(paste("Exporting data to Stata file at:", OUTPUT_DTA_PATH, "\n"))
haven::write_dta(original_data, OUTPUT_DTA_PATH)

cat("Export complete.\n")
cat("--- Script 07 Finished ---\n") 