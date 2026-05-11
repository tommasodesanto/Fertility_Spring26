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


# 2. SAFER, YEAR-BY-YEAR DATA PROCESSING
# ==============================================================================
cat("\nStarting year-by-year processing of the Stata file to conserve memory...\n")

# This approach is slower as it re-reads the source DTA file for each year,
# but it is much safer for systems with memory constraints.

for (year_to_get in years_to_process) {
    cat(paste("\n--- Processing Year:", year_to_get, "---\n"))
    
    # List to hold chunks for THIS year only
    current_year_chunks <- list()

    # Define the callback function for chunk processing
    process_chunk_for_year <- function(chunk, pos) {
        chunk_dt <- as.data.table(chunk)
        
        # Filter the chunk for data belonging to the current target year
        year_subset <- chunk_dt[year == year_to_get]
        
        if (nrow(year_subset) > 0) {
            # Convert all haven_labelled columns to their base types
            for (col in names(year_subset)) {
                if (is.labelled(year_subset[[col]])) {
                    year_subset[, (col) := haven::zap_labels(get(col))]
                }
            }
            # Add the processed chunk to our list for the current year
            current_year_chunks[[length(current_year_chunks) + 1]] <<- year_subset
        }
        
        if (pos %% (CHUNK_SIZE * 10) == 0) { # Print progress less frequently
             cat(paste("  - Scanned up to row:", pos, "\n"))
        }
    }

    # Read the entire DTA file in chunks, but the callback only keeps data for the target year
    read_dta_chunked(
        file = RAW_DATA_PATH, 
        callback = process_chunk_for_year,
        chunk_size = CHUNK_SIZE
    )

    # Combine all chunks for the current year and save to a single RDS file
    if (length(current_year_chunks) > 0) {
        cat(paste("  - Combining", length(current_year_chunks), "chunks for year", year_to_get, "...\n"))
        full_year_data <- rbindlist(current_year_chunks, use.names = TRUE, fill = TRUE)
        
        # Convert any remaining labelled columns to character for consistency
        full_year_data <- full_year_data %>%
            mutate(across(where(is.labelled), as.character))
            
        output_rds_path <- file.path(PROCESSED_DATA_DIR, paste0("fertility_microdata_", year_to_get, ".rds"))
        
        saveRDS(full_year_data, file = output_rds_path)
        cat(paste("    - SUCCESS: Saved", nrow(full_year_data), "records to", basename(output_rds_path), "\n"))
    } else {
        cat(paste("    - WARNING: No data found for year", year_to_get, ". Skipping.\n"))
    }
    
    # Explicitly free up memory before starting the next year
    rm(current_year_chunks)
    if (exists("full_year_data")) {
      rm(full_year_data)
    }
    gc()
}

cat("\nScript finished. All remaining years processed and saved.\n") 