# -----------------------------------------------------------------------------
# 03_clean_and_prep_microdata.R
#
# Description:
# This script loads the yearly IPUMS microdata (.rds files), takes a 5% 
# stratified sample for each year, combines them, cleans the result,
# and saves a final, analysis-ready dataset. This approach works with the
# pre-processed yearly files to remain memory-efficient.
# -----------------------------------------------------------------------------

# --- 1. Load Packages ---
suppressPackageStartupMessages({
    library(data.table)
    library(here)
    library(dplyr)
    library(haven)
})

# --- 2. Define Paths and Parameters ---
yearly_rds_dir <- here("Spatial_aggregate_withmicrodata", "processed_data", "yearly_rds_v7")
clean_data_path <- here("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_10pct_stratified_sample_v7.rds")
YEAR_RANGE <- 2005:2023
SAMPLE_FRACTION <- 0.10
set.seed(123) # for reproducibility

# --- 3. Year-by-Year Loading and Sampling from RDS files ---
message("Starting year-by-year loading and sampling from RDS files...")

all_year_samples <- list()

for (current_year in YEAR_RANGE) {
    message(paste("Processing year:", current_year))
    
    rds_file_path <- file.path(yearly_rds_dir, paste0("fertility_microdata_", current_year, ".rds"))

    if (!file.exists(rds_file_path)) {
        message(paste("  - RDS file for year", current_year, "not found. Skipping."))
        next
    }
    
    # Load the data for the current year
    year_data <- readRDS(rds_file_path)
    setDT(year_data) # Ensure it's a data.table
    message(paste("  - Loaded", nrow(year_data), "records for", current_year))

    # --- Type Conversion Step ---
    # Convert all haven_labelled columns to numeric to ensure type consistency
    # across all yearly files before binding them together.
    labelled_cols <- names(year_data)[sapply(year_data, is.labelled)]
    if (length(labelled_cols) > 0) {
        message(paste("  - Converting", length(labelled_cols), "labelled columns to numeric."))
        year_data[, (labelled_cols) := lapply(.SD, haven::zap_labels), .SDcols = labelled_cols]
    }

    # Take a 5% subsample
    sample_size <- floor(nrow(year_data) * SAMPLE_FRACTION)
    if (sample_size > 0) {
        year_sample <- year_data[sample(.N, sample_size)]
        message(paste("  - Sampled", nrow(year_sample), "records for", current_year))
        all_year_samples[[as.character(current_year)]] <- year_sample
    } else {
        message(paste("  - Sample size is zero for year", current_year, ". Nothing to sample."))
    }
    
    # Clean up memory
    rm(year_data)
    gc()
}

# Combine all yearly samples into one data.table
message("Combining all yearly samples...")
full_sample_dt <- rbindlist(all_year_samples, use.names = TRUE, fill = TRUE)
message(paste("Total rows in final full sampled dataset:", nrow(full_sample_dt)))

# Rename columns to lowercase for consistency before saving
setnames(full_sample_dt, old = names(full_sample_dt), new = tolower(names(full_sample_dt)))

# --- 4. FEATURE ENGINEERING AND CLEANING ---
message("Starting feature engineering on the combined sample...")

# Define the condition for women of childbearing age to apply fertility-specific logic
childbearing_condition <- full_sample_dt$sex == 2 & full_sample_dt$age >= 15 & full_sample_dt$age <= 50
message(paste("Identified", sum(childbearing_condition), "women of childbearing age for feature engineering."))

# Engineer new variables ONLY for the relevant subset (women of childbearing age)
# The `:=` operator modifies the data.table in place.
full_sample_dt[childbearing_condition, `:=`(
    # is_recent_mother: Binary flag for giving birth in the last year (from `fertyr`)
    is_recent_mother = ifelse(fertyr == 2, 1, 0),
    
    # age_at_first_birth: Key variable for our analysis
    age_at_first_birth = fcase(
        fertyr == 2 & nchild == 1, as.double(age),
        default = NA_real_
    ),
    
    # educ_factor: Categorical variable for education
    educ_factor = factor(fcase(
        educd < 60, "Less than High School",
        educd %in% 60:65, "High School Grad",
        educd %in% 70:90, "Some College",
        educd == 101, "Bachelors",
        educd > 101, "Graduate",
        default = "NIU/Unknown"
    )),

    # migration_factor: Categorical variable for migration status
    migration_factor = factor(fcase(
        migrate1 == 1, "Same House",
        migrate1 %in% 2:4, "Moved within US",
        migrate1 == 5, "Moved from Abroad",
        default = "NIU/Unknown"
    ))
)]

message("Feature engineering complete. The full dataset, including men, has been retained.")


# --- 5. SAVE THE FINAL, FULL SAMPLE FOR ANALYSIS ---
# This is the final, analysis-ready dataset. It contains ALL individuals.
message(paste("Saving final analysis-ready full sample to:", clean_data_path))
saveRDS(full_sample_dt, file = clean_data_path)

message("Script 03 finished successfully, producing the final FULL sample.") 