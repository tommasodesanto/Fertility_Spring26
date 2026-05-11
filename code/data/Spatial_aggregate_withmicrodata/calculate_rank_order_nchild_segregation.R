# ==============================================================================
# Utility Script: Calculate Rank-Order Family Size Segregation Index
#
# Description:
# This script implements the Rank-Order Information Theory Index to measure
# family size (nchild) segregation within each MSA. This is a sophisticated
# measure for continuous/ordered variables, as described in academic literature.
#
# The index H^R measures the extent to which the distribution of family sizes
# of smaller units (PUMAs) deviates from the distribution of the broader
# geographic unit (the MSA).
#
# - H^R = 0: No segregation. All PUMAs have the same distribution of family
#   sizes as the parent MSA.
# - H^R = 1: Extreme segregation. PUMAs are perfectly stratified by family size.
# ==============================================================================

# ---- 1. Setup ----
cat("--- Setting up environment... ---\n")
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}
install_and_load(c("data.table", "tidyverse", "spatstat.geom"))


# ---- 2. Configuration ----
cat("--- Configuring file paths... ---\n")
input_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
output_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "msa_rank_order_nchild_segregation.csv")
output_map_file <- file.path("Spatial_aggregate_withmicrodata", "analysis_output", "temp_run", "MSA_rank_order_nchild_segregation_map.png")
PERCENTILES <- seq(0.01, 0.99, by = 0.01) # Define the discrete percentiles for integration


# ---- 3. Data Loading ----
cat("--- Loading and preparing microdata... ---\n")
if (!file.exists(input_file)) { stop(paste("Input data file not found:", input_file)) }
acs_data <- readRDS(input_file)
setDT(acs_data)
# We need valid nchild values and household weights
acs_data <- acs_data[!is.na(nchild) & !is.na(hhwt) & hhwt > 0]


# ---- 4. Rank-Order Index Calculation Function ----
cat("--- Defining the Rank-Order Index calculation function... ---\n")

calculate_rank_order_index <- function(df, p_vector) {
  # Total population of the broader unit (MSA)
  total_msa_pop <- sum(df$hhwt)
  
  # Get the nchild values at each percentile for the whole MSA
  msa_nchild_quantiles <- weighted.quantile(df$nchild, df$hhwt, probs = p_vector)
  
  # Calculate omega_i (population share) for each PUMA
  puma_pop_shares <- df[, .(puma_pop = sum(hhwt)), by = puma][, omega_i := puma_pop / total_msa_pop]
  
  # Function to calculate F_bar_i,p
  # This is the CDF of a PUMA's nchild distribution, evaluated at a specific nchild value.
  get_puma_cdf <- function(nchild_val, puma_df) {
    weighted.mean(puma_df$nchild <= nchild_val, puma_df$hhwt)
  }
  
  # Calculate all F_bar values. This is the most intensive part.
  f_bar_results <- CJ(puma = unique(df$puma), p = p_vector)
  f_bar_results[, y_p := msa_nchild_quantiles[match(p, p_vector)]]
  f_bar_results[, F_bar_i_p := get_puma_cdf(y_p, df[puma == .BY$puma]), by = .(puma, p)]
  
  # Join with population shares
  f_bar_results <- merge(f_bar_results, puma_pop_shares, by = "puma")
  
  # Calculate the divergence term inside the summation
  f_bar_results[, term1 := ifelse(F_bar_i_p == 0, 0, F_bar_i_p * log(F_bar_i_p / p))]
  f_bar_results[, term2 := ifelse(F_bar_i_p == 1, 0, (1 - F_bar_i_p) * log((1 - F_bar_i_p) / (1 - p)))]
  f_bar_results[, divergence := omega_i * (term1 + term2)]
  
  # Sum over all PUMAs for each percentile p
  integral_components <- f_bar_results[, .(sum_divergence = sum(divergence)), by = p]
  
  # Approximate the integral by summing the components and multiplying by step width
  step_width <- p_vector[2] - p_vector[1]
  integral_approx <- sum(integral_components$sum_divergence) * step_width
  
  # Final H^R is 2 times the integral approximation
  rank_order_index <- 2 * integral_approx
  
  return(rank_order_index)
}


# ---- 5. Apply Calculation to all MSAs ----
cat("--- Applying calculation to all MSA-years (this may take a while)... ---\n")
# Filter for MSAs with more than one PUMA, as segregation is undefined for one unit.
msa_puma_counts <- acs_data[, .(num_pumas = uniqueN(puma)), by = .(met2013, year)]
msas_to_process <- msa_puma_counts[num_pumas > 1]

# Apply the function to each valid MSA-year group
segregation_results <- msas_to_process[, .(
  rank_order_seg_index = calculate_rank_order_index(acs_data[met2013 == .BY$met2013 & year == .BY$year], PERCENTILES)
), by = .(met2013, year)]

cat("--- Calculation complete. ---\n")


# ---- 6. Save Data and Generate Map ----
cat("--- Saving data and generating map... ---\n")
# Add names and save CSV
msa_geometries_for_names <- tigris::core_based_statistical_areas(cb = TRUE, year = 2021)
msa_name_lookup <- as.data.table(msa_geometries_for_names)[, .(met2013 = GEOID, msa_name = NAME)]
segregation_results[, met2013 := as.character(met2013)]
segregation_results <- merge(segregation_results, msa_name_lookup, by = "met2013", all.x = TRUE)
setcolorder(segregation_results, c("met2013", "msa_name", "year"))
fwrite(segregation_results, output_file)
cat(paste0("--- Segregation index data saved to ", output_file, " ---\n"))

# Map the results
states_geom <- tigris::states(cb = TRUE) %>% sf::st_transform(sf::st_crs(msa_geometries_for_names))
latest_year <- max(segregation_results$year, na.rm = TRUE)
map_data_sf <- merge(
  msa_geometries_for_names,
  segregation_results[year == latest_year],
  by.x = "GEOID",
  by.y = "met2013"
) %>% sf::st_as_sf()

segregation_map <- ggplot() +
  geom_sf(data = states_geom, fill = "gray90", color = "white", linewidth = 0.5) +
  geom_sf(data = map_data_sf, aes(fill = rank_order_seg_index), color = NA) +
  coord_sf(crs = sf::st_crs(5070), xlim = c(-2500000, 2500000), ylim = c(-2300000, 3000000)) +
  scale_fill_viridis_c(
    option = "magma", name = "Rank-Order Family Size\\nSegregation Index (H^R)",
    guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.position = "top")
  ) +
  labs(
    title = "Family Size Segregation Within U.S. Metropolitan Areas",
    subtitle = paste("Based on the Rank-Order Index of nchild for", latest_year),
    caption = "Source: ACS Microdata. High values (yellow) indicate PUMA family size distributions differ greatly from the overall MSA."
  ) +
  theme_void() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

ggsave(output_map_file, plot = segregation_map, width = 11, height = 8.5, dpi = 300, bg = "white")
cat(paste0("--- Map saved to ", output_map_file, " ---\n"))
cat("--- Script complete. ---\n") 