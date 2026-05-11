# ==============================================================================
# create_datasheet_for_agent.R
#
# Description:
# This script creates a self-contained CSV file that is ideal for sharing
# with another agent or team member. It contains:
# 1. A header row with the variable names.
# 2. A second row containing the description for each variable.
# 3. The first 50 rows of data from the final dataset.
# ==============================================================================

# --- 1. SETUP & CONFIGURATION ---
cat("Setting up the environment...\n")
suppressPackageStartupMessages({
    library(haven)
    library(data.table)
    library(here)
    library(dplyr)
})

# Define file paths
RAW_DATA_PATH <- here("Spatial_aggregate_withmicrodata", "raw_data", "extract26.dta")
CLEAN_MICRODATA_PATH <- here("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
OUTPUT_CSV_PATH <- here("Spatial_aggregate_withmicrodata", "output", "data_summary_for_agent.csv")

# --- 2. PREPARE DATA ---
cat("Loading data sample and variable labels...\n")

# Load the final dataset and take the first 50 rows
full_data <- readRDS(CLEAN_MICRODATA_PATH)
data_sample <- head(full_data, 50)

# Extract variable labels from the original Stata file
stata_metadata <- read_dta(RAW_DATA_PATH, n_max = 0)
variable_labels <- sapply(stata_metadata, function(x) {
    lbl <- attr(x, "label")
    if (is.null(lbl)) "" else lbl
})

# Create a data frame for the labels/descriptions
variable_descriptions <- data.frame(
    VariableName = names(variable_labels),
    Description = variable_labels
)

# Add descriptions for the variables we created manually
new_vars_desc <- data.frame(
    VariableName = c("is_recent_mother", "age_at_first_birth", "educ_factor", "migration_factor"),
    Description = c(
        "Binary flag (1/0) for recent birth. Women 15-50 only.",
        "Mother's age at first birth. Women 15-50 only.",
        "Factor for educational attainment. Women 15-50 only.",
        "Factor for migration status. Women 15-50 only."
    )
)
full_descriptions <- bind_rows(variable_descriptions, new_vars_desc)

# Match descriptions to the order of columns in the final data
# This is crucial to ensure alignment
ordered_descriptions <- full_descriptions$Description[match(names(data_sample), full_descriptions$VariableName)]

# Create a data frame where the first row is the descriptions
# This requires converting everything to character
desc_row <- as.data.frame(t(ordered_descriptions), stringsAsFactors = FALSE)
names(desc_row) <- names(data_sample)

# Convert the data sample to character to allow binding
data_sample_char <- as.data.frame(lapply(data_sample, as.character))

# Combine the description row and the data sample
# The header will be written from data_sample_char's names
final_export_data <- bind_rows(desc_row, data_sample_char)

cat("Data prepared successfully.\n")

# --- 3. CREATE AND SAVE CSV FILE ---
cat(paste("Creating CSV file at:", OUTPUT_CSV_PATH, "\n"))

# Write to CSV. The variable names will be the header by default.
fwrite(final_export_data, file = OUTPUT_CSV_PATH)

cat("CSV file created successfully.\n")
cat("The first row of the CSV contains the variable descriptions.\n") 