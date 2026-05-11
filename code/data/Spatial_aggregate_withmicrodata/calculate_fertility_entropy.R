# ==============================================================================
# Utility Script: Calculate MSA Fertility Entropy
#
# Description:
# This script calculates a fertility entropy index for each MSA and year based on
# the distribution of births across PUMAs within that MSA.
#
# Entropy measures the evenness of a distribution. In this context:
# - High Entropy (value near 1): Births are evenly distributed (high spatial dispersion)
#   across the MSA's many PUMAs. This indicates LOW concentration.
# - Low Entropy (value near 0): Births are concentrated in a small number
#   of PUMAs. This indicates HIGH concentration.
#
# The script also calculates a normalized entropy score, which adjusts for the
# number of PUMAs in an MSA, making scores comparable across different MSAs.
# ==============================================================================

# ---- 1. Setup ----
cat("--- Setting up environment... ---\n")
# This function ensures packages are installed and loaded
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}
install_and_load(c("data.table", "tidyverse"))


# ---- 2. Configuration ----
cat("--- Configuring file paths... ---\n")
# Relative path from the project root
input_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
output_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "msa_fertility_entropy.csv")


# ---- 3. Data Loading ----
cat("--- Loading microdata... ---\n")
if (!file.exists(input_file)) {
  stop(paste("Input data file not found:", input_file))
}
acs_data <- readRDS(input_file)
setDT(acs_data)


# ---- 4. Entropy Calculation ----
cat("--- Calculating fertility entropy for each MSA-year... ---\n")

# Step 1: Calculate weighted births per PUMA.
# We focus only on women who have recently had a birth.
puma_births <- acs_data[is_recent_mother == 1, .(
  weighted_births = sum(perwt, na.rm = TRUE)
), by = .(met2013, year, puma)]

# Step 2: Calculate total births per MSA and the proportion for each PUMA.
# This calculates the 'p_i' term for the standard Shannon entropy formula.
puma_births[, msa_total_births := sum(weighted_births), by = .(met2013, year)]
puma_births <- puma_births[msa_total_births > 0] # Avoid division by zero
puma_births[, proportion := weighted_births / msa_total_births]

# Step 3: Calculate the entropy term for each PUMA.
# The formula is p * log2(p). By convention, 0 * log(0) is treated as 0.
puma_births[proportion > 0, entropy_term := proportion * log2(proportion)]
puma_births[is.na(entropy_term), entropy_term := 0]

# Step 4: Sum the terms to get the final entropy for the MSA.
# We also count the number of PUMAs (N) in each MSA to calculate max and normalized entropy.
msa_entropy <- puma_births[, .(
  fertility_entropy = -sum(entropy_term, na.rm = TRUE),
  num_pumas = uniqueN(puma)
), by = .(met2013, year)]

# Step 5: Calculate normalized entropy.
# H_normalized = H / log2(N). This is useful for comparing MSAs of different sizes.
# If an MSA has only 1 PUMA, log2(1) = 0, and entropy is also 0.
# We define normalized entropy as 0 in this case to avoid division by zero.
msa_entropy[num_pumas > 1, normalized_fertility_entropy := fertility_entropy / log2(num_pumas)]
msa_entropy[num_pumas <= 1, normalized_fertility_entropy := 0]


# ---- 5. Add MSA Names and Save Output ----
cat("--- Adding MSA names to entropy data... ---\n")
# Fetch geometries just to get the names associated with the GEOIDs, with fallbacks
msa_geometries_for_names <- tryCatch({
  tigris::core_based_statistical_areas(cb = TRUE, year = 2022)
}, error = function(e) {
  cat("  -> Failed to get 2022 MSA geometries for names, trying 2021...\\n")
  tryCatch({
    tigris::core_based_statistical_areas(cb = TRUE, year = 2021)
  }, error = function(e2) {
    cat("  -> Failed to get 2021 MSA geometries for names, trying 2020...\\n")
    tigris::core_based_statistical_areas(cb = TRUE, year = 2020)
  })
})

msa_name_lookup <- as.data.table(msa_geometries_for_names)[, .(met2013 = GEOID, msa_name = NAME)]
msa_entropy[, met2013 := as.character(met2013)] # Ensure key is character for merge
msa_entropy <- merge(msa_entropy, msa_name_lookup, by = "met2013", all.x = TRUE)

# Re-order columns for better readability
setcolorder(msa_entropy, c("met2013", "msa_name", "year"))

cat(paste0("--- Saving entropy data to ", output_file, " ---\n"))
fwrite(msa_entropy, output_file)

cat("--- Script complete. ---\n")

# Display a sample of the results for verification
cat("Sample of the calculated entropy data (top 5 MSAs by normalized entropy):\n")
print(head(msa_entropy[order(normalized_fertility_entropy)], 5))

cat("\nSample of the calculated entropy data (bottom 5 MSAs by normalized entropy):\n")
print(head(msa_entropy[order(normalized_fertility_entropy)], 5))


# ---- 6. Generate and Save a Map ----
cat("\n--- Generating a map of normalized fertility entropy... ---\n")
output_map_file <- file.path("Spatial_aggregate_withmicrodata", "analysis_output", "temp_run", "MSA_normalized_fertility_entropy_map.png")

# Get MSA geometries using tigris, with fallbacks
msa_geometries <- tryCatch({
  tigris::core_based_statistical_areas(cb = TRUE, year = 2022)
}, error = function(e) {
  cat("  -> Failed to get 2022 MSA geometries, trying 2021...\\n")
  tryCatch({
    tigris::core_based_statistical_areas(cb = TRUE, year = 2021)
  }, error = function(e2) {
    cat("  -> Failed to get 2021 MSA geometries, trying 2020...\\n")
    tigris::core_based_statistical_areas(cb = TRUE, year = 2020)
  })
})

# Get state geometries for background context
states_geom <- tigris::states(cb = TRUE) %>% sf::st_transform(sf::st_crs(msa_geometries))

# Prepare data for mapping
latest_year <- max(msa_entropy$year, na.rm = TRUE)
map_data <- msa_entropy[year == latest_year]

# Join entropy data with geometries
map_data_sf <- merge(
  msa_geometries,
  map_data,
  by.x = "GEOID",
  by.y = "met2013",
  all.y = TRUE
) %>% sf::st_as_sf()


# Create the plot
entropy_map <- ggplot() +
  # Add state outlines as a base layer
  geom_sf(data = states_geom, fill = "gray90", color = "white", linewidth = 0.5) +
  # Add the MSA data on top
  geom_sf(data = map_data_sf, aes(fill = normalized_fertility_entropy), color = "white", linewidth = 0.1) +
  # Use a projection suitable for the continental US
  coord_sf(crs = sf::st_crs(5070), xlim = c(-2500000, 2500000), ylim = c(-2300000, 3000000)) + # Albers Equal Area
  scale_fill_viridis_c(
    option = "plasma",
    name = "Normalized Entropy",
    guide = guide_colorbar(
      barwidth = 15,
      barheight = 0.5,
      title.position = "top"
    )
  ) +
  labs(
    title = "Fertility Concentration Across U.S. Metropolitan Areas",
    subtitle = paste("Normalized Entropy of Births Across PUMAs for", latest_year),
    caption = "Source: ACS Microdata. High values (yellow) indicate births are evenly spread; low values (purple) indicate concentration."
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    plot.caption = element_text(hjust = 0.5)
  )

# Save the map
ggsave(
  output_map_file,
  plot = entropy_map,
  width = 11,
  height = 8.5,
  dpi = 300,
  bg = "white"
)

cat(paste0("--- Map saved to ", output_map_file, " ---\n")) 