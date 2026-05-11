# ==============================================================================
# AGGREGATE MICRODATA AND VALIDATE AGAINST ORIGINAL ANALYSIS
#
# Description:
# This script loads the cleaned microdata, maps PUMAs to MSAs using
# geographic crosswalk files, and aggregates the data to the MSA-year level.
# It then merges this new data with the original aggregate analysis data
# to perform a visual validation, ensuring the trends are consistent.
# ==============================================================================

# 1. SETUP & CONFIGURATION
# ==============================================================================
cat("Setting up the environment...\n")

suppressPackageStartupMessages({
    library(data.table)
    library(here)
    library(ggplot2)
})

# Define file paths
CLEAN_MICRODATA_PATH <- here("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean.rds")
CROSSWALK_2000_PATH <- here("Spatial_aggregate_withmicrodata", "crosswalks", "geocorr_puma2000_to_cbsa2013.csv")
CROSSWALK_2010_PATH <- here("Spatial_aggregate_withmicrodata", "crosswalks", "geocorr_puma2010_to_cbsa2013.csv")
ORIGINAL_AGG_DATA_PATH <- here("Spatial_aggregate_analysis", "data_final", "msa_panel_analysis_ready.rds")
VALIDATION_OUTPUT_DIR <- here("Spatial_aggregate_withmicrodata", "validation")
dir.create(VALIDATION_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 2. LOAD DATA & CROSSWALKS
# ==============================================================================
cat("Loading cleaned microdata, crosswalks, and original aggregate data...\n")

# Load the analysis-ready microdata
if (!file.exists(CLEAN_MICRODATA_PATH)) stop("Cleaned microdata file not found.")
micro_dt <- readRDS(CLEAN_MICRODATA_PATH)

# Load crosswalks
if (!file.exists(CROSSWALK_2000_PATH)) stop("2000 PUMA crosswalk not found.")
cw_2000 <- fread(CROSSWALK_2000_PATH)

if (!file.exists(CROSSWALK_2010_PATH)) stop("2010 PUMA crosswalk not found.")
cw_2010 <- fread(CROSSWALK_2010_PATH)

# Load the final dataset from the original analysis for comparison
if (!file.exists(ORIGINAL_AGG_DATA_PATH)) stop("Original aggregate data file not found.")
agg_df <- readRDS(ORIGINAL_AGG_DATA_PATH)

cat("Data loaded successfully.\n\n")

# 3. APPLY CROSSWALK TO MAP PUMA TO CBSA
# ==============================================================================
cat("Applying geographic crosswalks...\n")

# Standardize column names in crosswalks
setnames(cw_2000, old = c("puma2k", "cbsa13", "afact"), new = c("puma", "cbsa", "afact"))
setnames(cw_2010, old = c("puma10", "cbsa13", "afact"), new = c("puma", "cbsa", "afact"))

# Split microdata by PUMA vintage
micro_2000_era <- micro_dt[year <= 2011]
micro_2010_era <- micro_dt[year >= 2012]

# Merge each part with the correct crosswalk
# We keep all records from the microdata, even if they don't match a CBSA
merged_2000 <- merge(micro_2000_era, cw_2000[, .(state, puma, cbsa, afact)], by = c("state", "puma"), all.x = TRUE)
merged_2010 <- merge(micro_2010_era, cw_2010[, .(state, puma, cbsa, afact)], by = c("state", "puma"), all.x = TRUE)

# Recombine the data
micro_dt_mapped <- rbindlist(list(merged_2000, merged_2010))

# Adjust person weight by the allocation factor. This correctly allocates
# portions of a PUMA's population to the appropriate CBSA.
micro_dt_mapped[, perwt_adj := perwt * afact]

cat("PUMA to CBSA mapping complete.\n\n")


# 4. AGGREGATE MICRODATA TO MSA-YEAR LEVEL
# ==============================================================================
cat("Aggregating microdata to MSA-year level...\n")

# Filter out non-metro areas and calculate weighted stats
msa_agg_from_micro <- micro_dt_mapped[!is.na(cbsa) & cbsa != 99999, .(
    mean_age_birth_micro = weighted.mean(age[is_recent_mother == 1], perwt_adj[is_recent_mother == 1], na.rm = TRUE),
    mean_age_first_birth_micro = weighted.mean(age_at_first_birth, perwt_adj, na.rm = TRUE),
    birth_rate_micro = (sum(is_recent_mother * perwt_adj, na.rm = TRUE) / sum(perwt_adj, na.rm = TRUE)) * 1000
), by = .(year, cbsa)]

cat("Aggregation complete.\n\n")

# 5. MERGE AND COMPARE
# ==============================================================================
cat("Merging datasets for comparison...\n")

# Prepare original aggregate data
setDT(agg_df)
agg_df_renamed <- agg_df[, .(year = as.integer(year), cbsa = as.integer(GEOID), mean_age_birth_agg = mean_age_at_birth)]

# Merge for comparison
validation_data <- merge(msa_agg_from_micro, agg_df_renamed, by = c("year", "cbsa"))

cat(paste("Validation data has", nrow(validation_data), "MSA-year observations for comparison.\n\n"))


# 6. VISUAL VALIDATION
# ==============================================================================
cat("Generating validation plots...\n")

# Scatter plot of mean age at birth
p_scatter <- ggplot(validation_data, aes(x = mean_age_birth_agg, y = mean_age_birth_micro)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(
        title = "Validation: Mean Age at Birth (Microdata vs. Aggregate)",
        subtitle = "Each point is an MSA-year observation",
        x = "Mean Age (Original Aggregate Data)",
        y = "Mean Age (Calculated from Microdata)"
    ) +
    theme_minimal()
ggsave(here(VALIDATION_OUTPUT_DIR, "validation_scatter_mean_age.png"), p_scatter, width = 8, height = 6)

# Time series plot of the national average
national_trends <- validation_data[, .(
    agg = mean(mean_age_birth_agg, na.rm = TRUE),
    micro = mean(mean_age_birth_micro, na.rm = TRUE)
), by = year]

p_timeseries <- ggplot(national_trends, aes(x = year)) +
    geom_line(aes(y = agg, color = "Original Aggregate"), size = 1) +
    geom_line(aes(y = micro, color = "From Microdata"), size = 1) +
    labs(
        title = "Validation: National Average Mean Age at Birth",
        subtitle = "Comparing trends from both data sources",
        x = "Year", y = "Average Mean Age at Birth", color = "Data Source"
    ) +
    theme_minimal()
ggsave(here(VALIDATION_OUTPUT_DIR, "validation_timeseries_mean_age.png"), p_timeseries, width = 8, height = 6)


cat("Validation plots saved. Script finished.\n") 