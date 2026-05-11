# ==============================================================================
# IPUMS MICRODATA EXTRACT FOR FERTILITY ANALYSIS
#
# Description: This script defines and submits an IPUMS extract request to
#              download ACS 1-year microdata for a detailed fertility analysis.
#              It includes a comprehensive set of geographic, migration, and
#              socio-economic variables for nuanced research.
#
# Author: Gemini Assistant
# Date: [Current Date]
# ==============================================================================

# 1. SETUP & CONFIGURATION
# ==============================================================================
# Ensure the ipumsr package is installed and loaded.
if (!require(ipumsr)) {
    install.packages("ipumsr")
    library(ipumsr)
}

# The `here` package helps manage file paths.
if (!require(here)) {
    install.packages("here")
    library(here)
}


# Set your IPUMS API key.
# It's recommended to set this in your .Renviron file for security.
# You can run `usethis::edit_r_environ()` to open it and add the line:
# IPUMS_API_KEY="YOUR_KEY_HERE"
# Sys.setenv(IPUMS_API_KEY = "YOUR_KEY_HERE") # Or uncomment and set here


# 2. DEFINE IPUMS EXTRACT
# ==============================================================================
# This defines the data extract we want from IPUMS USA.
# We are requesting ACS 1-year data from 2005 to 2023.

fertility_extract <- define_extract_micro(
  collection = "usa", # Specify the IPUMS collection
  description = "Fertility, Migration, and Child Age Microdata, 2005-2023",
  samples     = paste0("us", 2005:2023, "a"), # ACS 1-year samples
  variables   = c(
    # --- Core Technical Variables (Household and Person) ---
    "YEAR",      # Survey year
    "SAMPLE",    # IPUMS sample identifier
    "SERIAL",    # Household identifier
    "HHWT",      # Household weight
    "CLUSTER",   # Household cluster for variance estimation
    "STRATA",    # Household strata for variance estimation
    "PERNUM",    # Person number in household
    "PERWT",     # Person weight

    # --- Family Interrelationship Pointers (Crucial for child vars) ---
    "SPLOC",     # Location of spouse
    "MOMLOC",    # Location of mother
    "POPLOC",    # Location of father

    # --- Geographic Variables (Household Level) ---
    "REGION",    # Census region and division
    "STATEFIP",  # State FIPS code
    "COUNTYFIP", # County FIPS code (identifiable counties)
    "PUMA",      # Public Use Microdata Area
    "METRO",     # Metropolitan status
    "MET2013",   # Metropolitan area (2013 delineations)
    "DENSITY",   # Population-weighted density of PUMA

    # --- Fertility and Child Variables ---
    "FERTYR",    # Gave birth in the last year
    "NCHILD",    # Number of own children in the household
    "ELDCH",     # Age of eldest own child
    "YNGCH",     # Age of youngest own child

    # --- Individual Demographics ---
    "AGE",       # Age
    "SEX",       # Sex
    "MARST",     # Marital status
    "RACE",      # Race
    "BPL",       # Birthplace
    "EDUC",      # Educational attainment (general)
    "EDUCD",     # Educational attainment (detailed)

    # --- Individual Economic Characteristics ---
    "EMPSTAT",   # Employment status
    "INCTOT",    # Total personal income
    "INCWAGE",   # Wage and salary income
    "IND",       # Industry
    "OCC",       # Occupation

    # --- Migration and Place of Work Variables (Person Level) ---
    "MIGRATE1",  # Migration status, 1 year
    "MIGPUMA1",  # Migration PUMA, 1 year
    "MIGMET131"  # Migration metropolitan area (2013 delineations), 1 year
  ),
  data_format = "fixed_width" # Recommended format
)

# 3. SAVE EXTRACT DEFINITION
# ==============================================================================
# Save the extract definition object to a file. This is good practice for
# reproducibility.
EXTRACT_DEF_DIR <- here("Spatial_aggregate_withmicrodata", "extract_definitions")
dir.create(EXTRACT_DEF_DIR, showWarnings = FALSE, recursive = TRUE)

# Use a new name for the new extract definition
saveRDS(
  fertility_extract,
  file = here(EXTRACT_DEF_DIR, "fertility_extract_def_v3_child_age.rds")
)

cat("IPUMS extract definition has been saved to:\n")
cat(here(EXTRACT_DEF_DIR, "fertility_extract_def_v3_child_age.rds"), "\n\n")


# 4. SUBMIT EXTRACT AND DOWNLOAD DATA
# ==============================================================================
# This section submits the request and waits for the data to be prepared.
# IPUMS extracts can take time to process (from minutes to hours).
# The `wait_for_extract` function will check the status periodically.

# --- IMPORTANT ---
# Before running the code below, make sure your IPUMS API key is set.
# The following lines are commented out to prevent accidental submission.
# Uncomment them when you are ready to request the data from IPUMS.

cat("Submitting extract request to IPUMS...\n")

# Read the new definition back from the file
extract_def <- readRDS(here(EXTRACT_DEF_DIR, "fertility_extract_def_v3_child_age.rds"))

# Submit the extract and wait for it to be ready
submitted_extract <- submit_extract(extract_def)
waited_extract <- wait_for_extract(submitted_extract) # Wait for extract, no interval arg

# Define where to download the data
DOWNLOAD_DIR <- here("Spatial_aggregate_withmicrodata", "raw_data")
dir.create(DOWNLOAD_DIR, showWarnings = FALSE, recursive = TRUE)

# Download the data. This function will now wait for the extract to be ready.
download_extract(waited_extract, download_dir = DOWNLOAD_DIR, overwrite = TRUE)

cat("Download complete!\n")
cat("Your IPUMS microdata files are located in:", DOWNLOAD_DIR, "\n") 