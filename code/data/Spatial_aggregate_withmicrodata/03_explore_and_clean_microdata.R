# ==============================================================================
# EXPLORE, CLEAN, AND PREPARE IPUMS MICRODATA
#
# Description: This script loads the processed RDS file containing the IPUMS
#              microdata. It performs exploratory data analysis (EDA), cleans
#              variables, and engineers new features needed for the fertility
#              analysis, such as identifying first-time mothers.
#
# Author: Gemini Assistant
# Date: [Current Date]
# ==============================================================================

# 1. SETUP & CONFIGURATION
# ==============================================================================
cat("Setting up the environment...\n")

# Install and load necessary packages
suppressPackageStartupMessages({
    required_packages <- c("tidyverse", "here", "arrow", "haven")
    for (pkg in required_packages) {
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg, repos = "https://cloud.r-project.org/")
            library(pkg, character.only = TRUE)
        }
    }
})

# Define file paths
PROCESSED_DATA_PATH <- here("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_2005_2023.rds")
CLEANED_DATA_DIR <- here("Spatial_aggregate_withmicrodata", "processed_data")
OUTPUT_PARQUET_PATH <- here(CLEANED_DATA_DIR, "fertility_microdata_cleaned.parquet")

cat(paste("Input RDS file:", PROCESSED_DATA_PATH, "\n"))
cat(paste("Output Parquet file:", OUTPUT_PARQUET_PATH, "\n\n"))


# 2. LOAD DATA
# ==============================================================================
cat("Loading processed microdata...\n")

if (!file.exists(PROCESSED_DATA_PATH)) {
    stop("The processed RDS file does not exist. Please run the previous script.")
}
microdata <- readRDS(PROCESSED_DATA_PATH)
cat("Data loaded successfully.\n\n")

# Convert 'sex' column from haven_labelled to numeric to avoid filter errors
if (inherits(microdata$sex, "haven_labelled")) {
    microdata$sex <- as.numeric(as_factor(microdata$sex))
}


# 3. INITIAL FILTERING AND SELECTION
# ==============================================================================
cat("Filtering data to women of childbearing age and selecting columns...\n")

# We use ipums_val_labels to see what the numbers mean
# print(ipumsr::ipums_val_labels(microdata$sex)) # 1=Male, 2=Female

cleaned_data <- microdata %>%
    # Filter for women aged 15-50
    filter(sex == 2, age >= 15, age <= 50) %>%
    # Select only the columns we need for this analysis
    select(year, perwt, age, fertyr, nchild, educ, marst, race, hispan,
           migrate1, puma, metarea, statefip, rent, valueh)

cat(paste("Filtered data has", nrow(cleaned_data), "rows.\n\n"))


# 4. VARIABLE CLEANING AND FEATURE ENGINEERING
# ==============================================================================
cat("Cleaning variables and engineering new features...\n")

cleaned_data <- cleaned_data %>%
    mutate(
        # Create `gave_birth`: A binary 0/1 variable from FERTYR (2 = Yes).
        # FERTYR codes: 1=N/A, 2=Yes. We assume 1 means no.
        gave_birth = case_when(
            fertyr == 2 ~ 1,
            TRUE ~ 0
        ),

        # Create `age_at_first_birth`: Contains the age for women who had their
        # first child this year, otherwise NA. This is a powerful proxy.
        # Logic: Gave birth in the last year AND have exactly one child.
        age_at_first_birth = case_when(
            gave_birth == 1 & nchild == 1 ~ age,
            TRUE ~ NA_real_
        ),

        # Recode `EDUC` into simpler categories for analysis.
        # EDUC codes: 02-05=Some HS, 06=HS Grad, 07-09=Some College, 10=Bachelors, 11=Graduate
        educ_cat = case_when(
            educ < 6 ~ "Less than High School",
            educ == 6 ~ "High School Grad",
            educ %in% 7:9 ~ "Some College",
            educ == 10 ~ "Bachelors",
            educ >= 11 ~ "Graduate",
            TRUE ~ "NIU/Unknown"
        ) %>% factor(levels = c("Less than High School", "High School Grad", "Some College",
                                "Bachelors", "Graduate", "NIU/Unknown")),

        # Recode `MARST` into simpler categories.
        # MARST codes: 1=Married, 2=Married/sep, 3=Divorced, 4=Widowed, 5=Separated, 6=Never married
        marst_cat = case_when(
            marst %in% 1:2 ~ "Married",
            marst == 6 ~ "Never Married",
            marst %in% 3:5 ~ "Divorced/Separated/Widowed",
            TRUE ~ "NIU/Unknown"
        ) %>% factor(levels = c("Married", "Never Married", "Divorced/Separated/Widowed", "NIU/Unknown")),
        
        # Create a simple dummy for whether the person moved in the last year
        moved_last_year = case_when(
            migrate1 > 1 ~ 1,
            TRUE ~ 0
        )
    )

cat("Feature engineering complete.\n\n")


# 5. EXPLORATORY DATA ANALYSIS (EDA)
# ==============================================================================
cat("Performing brief exploratory data analysis...\n")

# Calculate the weighted mean age at first birth for each year
# This is a key validation metric to see if our proxy makes sense.
yearly_summary <- cleaned_data %>%
    filter(!is.na(age_at_first_birth)) %>%
    group_by(year) %>%
    summarise(
        mean_age_first_birth = weighted.mean(age_at_first_birth, perwt, na.rm = TRUE),
        n_first_births = n()
    )

cat("--- Yearly Mean Age at First Birth (from Microdata) ---\n")
print(yearly_summary)
cat("\nThis table is a crucial sanity check. It shows the trend in our primary new outcome variable.\n\n")

# Frequency table for education categories
cat("--- Education Level Distribution (All Years) ---\n")
print(cleaned_data %>% count(educ_cat, wt = perwt))


# 6. SAVE CLEANED DATA
# ==============================================================================
cat("\nSaving cleaned data to Parquet format...\n")

# We will save the cleaned data in Parquet format. It's efficient for large
# datasets and preserves data types well.
write_parquet(cleaned_data, sink = OUTPUT_PARQUET_PATH)

cat(paste("Cleaned data successfully saved to:", OUTPUT_PARQUET_PATH, "\n"))
cat("Script finished.\n") 