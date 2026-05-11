# ==============================================================================
# Utility Script: Calculate MSA Fertility Segregation (Theil's H Index)
#
# Description:
# This script calculates the Theil's H entropy-based segregation index. It
# measures the extent to which fertility is segregated across PUMAs within an MSA.
#
# The index is calculated as H = (E_msa - E_pumas_avg) / E_msa, where:
# - E_msa is the overall entropy of the MSA.
# - E_pumas_avg is the population-weighted average of the PUMA entropies.
#
# The index ranges from 0 to 1:
# - H = 0 (No Segregation): Each PUMA's birth rate is identical to the MSA's.
# - H = 1 (Extreme Segregation): Each PUMA consists entirely of either
#   women who had a birth or women who did not.
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
install_and_load(c("data.table", "tidyverse", "tigris", "sf"))


# ---- 2. Configuration ----
cat("--- Configuring file paths... ---\n")
input_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
output_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "msa_fertility_segregation_index.csv")
output_map_file <- file.path("Spatial_aggregate_withmicrodata", "analysis_output", "temp_run", "MSA_fertility_segregation_map.png")


# ---- 3. Data Loading ----
cat("--- Loading microdata... ---\n")
if (!file.exists(input_file)) { stop(paste("Input data file not found:", input_file)) }
acs_data <- readRDS(input_file)
setDT(acs_data)


# ---- 4. Segregation Index Calculation ----
cat("--- Calculating Theil's H segregation index for each MSA-year... ---\n")

# Define a safe binary entropy function h(p) = -[p*log(p) + (1-p)*log(1-p)]
calc_entropy <- function(p) {
  p_safe <- p[p > 0 & p < 1] # Avoid log(0)
  entropy <- -sum(p_safe * log(p_safe) + (1 - p_safe) * log(1 - p_safe))
  return(ifelse(is.na(entropy), 0, entropy))
}

# Step 1: Calculate PUMA-level stats: birth rate (p_i) and population.
# We consider all women of child-bearing age for the rate calculation.
puma_stats <- acs_data[age %between% c(15, 50), .(
  puma_birth_rate = weighted.mean(is_recent_mother, perwt, na.rm = TRUE),
  puma_population = .N
), by = .(met2013, year, puma)]
puma_stats <- puma_stats[!is.na(puma_birth_rate)]

# Step 2: Calculate PUMA-level entropy (h_i).
puma_stats[, h_i := calc_entropy(puma_birth_rate)]

# Step 3: Calculate MSA-level components (H_hat and H_bar).
msa_stats <- puma_stats[, .(
  # H_bar: population-weighted average of PUMA entropies
  H_bar = weighted.mean(h_i, w = puma_population, na.rm = TRUE),
  # Also need total population and births to calculate MSA-level entropy
  msa_total_pop = sum(puma_population, na.rm = TRUE),
  msa_total_births = sum(puma_birth_rate * puma_population, na.rm = TRUE)
), by = .(met2013, year)]

# H_hat: MSA-level entropy
msa_stats[, msa_birth_rate := msa_total_births / msa_total_pop]
msa_stats[, H_hat := calc_entropy(msa_birth_rate)]

# Step 4: Calculate the final segregation index H.
msa_stats[, segregation_index_H := (H_hat - H_bar) / H_hat]
# Handle cases where H_hat is 0 (no diversity in MSA), segregation is 0.
msa_stats[H_hat == 0, segregation_index_H := 0]


# ---- 5. Add MSA Names and Save Output ----
cat("--- Adding MSA names and saving data... ---\n")
msa_geometries_for_names <- tryCatch({
  tigris::core_based_statistical_areas(cb = TRUE, year = 2022)
}, error = function(e) { tigris::core_based_statistical_areas(cb = TRUE, year = 2021) })
msa_name_lookup <- as.data.table(msa_geometries_for_names)[, .(met2013 = GEOID, msa_name = NAME)]
msa_stats[, met2013 := as.character(met2013)]
msa_stats <- merge(msa_stats, msa_name_lookup, by = "met2013", all.x = TRUE)
setcolorder(msa_stats, c("met2013", "msa_name", "year"))

fwrite(msa_stats, output_file)
cat(paste0("--- Segregation data saved to ", output_file, " ---\n"))


# ---- 6. Generate and Save a Map ----
cat("--- Generating a map of the segregation index... ---\n")
states_geom <- tigris::states(cb = TRUE) %>% sf::st_transform(sf::st_crs(msa_geometries_for_names))
latest_year <- max(msa_stats$year, na.rm = TRUE)
map_data_sf <- merge(
  msa_geometries_for_names,
  msa_stats[year == latest_year],
  by.x = "GEOID",
  by.y = "met2013"
) %>% sf::st_as_sf()

segregation_map <- ggplot() +
  geom_sf(data = states_geom, fill = "gray90", color = "white", linewidth = 0.5) +
  geom_sf(data = map_data_sf, aes(fill = segregation_index_H), color = "white", linewidth = 0.1) +
  coord_sf(crs = sf::st_crs(5070), xlim = c(-2500000, 2500000), ylim = c(-2300000, 3000000)) +
  scale_fill_viridis_c(
    option = "cividis", name = "Fertility Segregation (H)",
    guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.position = "top")
  ) +
  labs(
    title = "Fertility Segregation Within U.S. Metropolitan Areas",
    subtitle = paste("Theil's H Index for", latest_year),
    caption = "Source: ACS Microdata. High values (yellow) indicate high segregation of births across neighborhoods."
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    plot.caption = element_text(hjust = 0.5)
  )

ggsave(output_map_file, plot = segregation_map, width = 11, height = 8.5, dpi = 300, bg = "white")
cat(paste0("--- Map saved to ", output_map_file, " ---\n"))
cat("--- Script complete. ---\n") 