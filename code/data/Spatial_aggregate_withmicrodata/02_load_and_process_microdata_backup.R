# ==============================================================================
# LOAD, PROCESS, AND SAVE IPUMS MICRODATA (EFFICIENTLY)
#
# Description: This script reads the large Stata file ONCE, in chunks, to 
#              avoid memory issues. It processes and collects all years of data 
#              in a single pass, then writes each year to a separate RDS file.
#              This is dramatically faster than reading the source file for each year.
#
# Author: Gemini Assistant
# Date: [Current Date]
# ==============================================================================

# 1. SETUP & CONFIGURATION
# ==============================================================================
cat("Setting up the environment...\n")

# Install and load necessary packages
suppressPackageStartupMessages({
    required_packages <- c("haven", "here", "data.table", "dplyr", "purrr")
    for (pkg in required_packages) {
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg, repos = "https://cloud.r-project.org/")
            library(pkg, character.only = TRUE)
        }
    }
})

# --- Define Paths and Parameters ---
RAW_DATA_PATH <- here("Spatial_aggregate_withmicrodata", "raw_data", "extract27.dta")
PROCESSED_DATA_DIR <- here("Spatial_aggregate_withmicrodata", "processed_data", "yearly_rds_v7")
dir.create(PROCESSED_DATA_DIR, showWarnings = FALSE, recursive = TRUE)
YEAR_RANGE <- 2005:2023
CHUNK_SIZE <- 500000 # Increased chunk size for efficiency

cat(paste("Input Stata file:", RAW_DATA_PATH, "\n"))
cat(paste("Output directory for yearly RDS files:", PROCESSED_DATA_DIR, "\n\n"))

# --- Helper Function for Chunked DTA Reading ---
read_dta_chunked <- function(file, callback, chunk_size = 10000, ...) {
  skip <- 0
  while (TRUE) {
    chunk <- read_dta(file, skip = skip, n_max = chunk_size)
    if (nrow(chunk) == 0) break
    callback(chunk, skip)
    skip <- skip + chunk_size
  }
}

# --- Check which years already exist to avoid re-processing ---
existing_years <- as.numeric(gsub("fertility_microdata_|.rds", "", list.files(PROCESSED_DATA_DIR)))
years_to_process <- setdiff(YEAR_RANGE, existing_years)

if (length(years_to_process) == 0) {
    cat("All requested years have already been processed. Nothing to do.\n")
    quit(save = "no", status = 0)
}
cat(paste("Years to process:", paste(years_to_process, collapse = ", "), "\n"))


# 2. EFFICIENT DATA READING AND PROCESSING
# ==============================================================================
cat("Starting efficient single-pass processing of the Stata file...\n")

# Initialize an empty list to hold data for each year
yearly_data_agg <- list()

# --- Define the callback function for chunk processing ---
process_chunk <- function(chunk, pos) {
    chunk_dt <- as.data.table(chunk)

    # Convert all haven_labelled columns to their base types
    for (col in names(chunk_dt)) {
        if (is.labelled(chunk_dt[[col]])) {
            chunk_dt[, (col) := haven::zap_labels(get(col))]
        }
    }
    
    # Split the chunk by year
    split_by_year <- split(chunk_dt, by = "year", keep.by = TRUE)
    
    # Append the data to the correct year in our aggregation list
    for (year_str in names(split_by_year)) {
        if (is.null(yearly_data_agg[[year_str]])) {
            yearly_data_agg[[year_str]] <<- list()
        }
        yearly_data_agg[[year_str]][[length(yearly_data_agg[[year_str]]) + 1]] <<- split_by_year[[year_str]]
    }
    
    cat(paste("  - Processed chunk starting at row:", pos, "\n"))
    gc()
}

# --- Read the DTA file in chunks using the callback ---
read_dta_chunked(
    file = RAW_DATA_PATH, 
    callback = process_chunk,
    chunk_size = CHUNK_SIZE
)

cat("Finished reading and organizing data by year.\n")


# 3. SAVE EACH YEAR'S DATA TO A SEPARATE RDS FILE
# ==============================================================================
cat("Combining chunks and saving RDS file for each year...\n")

# Filter for the years the user actually wants
years_to_save <- intersect(names(yearly_data_agg), as.character(YEAR_RANGE))

# Use purrr to iterate over the years, combine chunks, and save
walk(years_to_save, ~{
    year_val <- .x
    cat(paste("  - Processing year:", year_val, "\n"))
    
    if (length(yearly_data_agg[[year_val]]) > 0) {
        full_year_data <- rbindlist(yearly_data_agg[[year_val]], use.names = TRUE, fill = TRUE)
        
        # Convert labelled columns to character for consistency
        full_year_data <- full_year_data %>%
            mutate(across(where(is.labelled), as.character))
            
        output_rds_path <- file.path(PROCESSED_DATA_DIR, paste0("fertility_microdata_", year_val, ".rds"))
        
        saveRDS(full_year_data, file = output_rds_path)
        cat(paste("    - Saved", nrow(full_year_data), "records to", basename(output_rds_path), "\n"))
    } else {
        cat(paste("    - No data found for year", year_val, ". Skipping.\n"))
    }
    
    # Free up memory
    yearly_data_agg[[year_val]] <- NULL
    gc()
})

cat("\nScript finished. All years processed and saved.\n") 