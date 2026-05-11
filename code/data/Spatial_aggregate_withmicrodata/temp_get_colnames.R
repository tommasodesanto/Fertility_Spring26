suppressPackageStartupMessages({
    library(data.table)
    library(here)
})

# Load the full, correct microdata sample
CLEAN_MICRODATA_PATH <- here("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
micro_data <- readRDS(CLEAN_MICRODATA_PATH)

# Get the column names
column_names <- names(micro_data)

cat("--- Column Names in the Final 5% Stratified Sample ---\n\n")
cat(paste(column_names, collapse = "\n"))
cat("\n\n--- End of List ---\n") 