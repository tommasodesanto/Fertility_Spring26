# ==============================================================================
# Utility Script: Calculate MSA Fertility Rate Inequality
#
# Description:
# This script measures the inequality in fertility rates across the PUMAs
# within each MSA for a given year.
#
# It does this by calculating the population-weighted standard deviation of the
# PUMA-level birth rates.
#
# - A HIGH value indicates high inequality: The MSA has PUMAs with very
#   different fertility rates (e.g., some high, some low).
# - A LOW value indicates low inequality (homogeneity): All PUMAs within
#   the MSA have very similar fertility rates.
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
install_and_load(c("data.table", "tidyverse"))


# ---- 2. Configuration ----
cat("--- Configuring file paths... ---\n")
input_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
output_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "msa_fertility_rate_inequality.csv")
output_map_file <- file.path("Spatial_aggregate_withmicrodata", "analysis_output", "temp_run", "MSA_fertility_rate_inequality_map.png")


# ---- 3. Data Loading ----
cat("--- Loading microdata... ---\n")
if (!file.exists(input_file)) {
  stop(paste("Input data file not found:", input_file))
}
acs_data <- readRDS(input_file)
setDT(acs_data)


# ---- 4. Inequality Calculation ----
cat("--- Calculating fertility rate inequality for each MSA-year... ---\n")

# Step 1: Calculate birth rate and population for each PUMA-year.
# We consider all women of child-bearing age for the rate calculation.
puma_stats <- acs_data[age %between% c(15, 50), .(
  puma_birth_rate = weighted.mean(is_recent_mother, perwt, na.rm = TRUE),
  puma_population = sum(perwt, na.rm = TRUE)
), by = .(met2013, year, puma)]

# Remove any PUMAs where a rate couldn't be calculated
puma_stats <- puma_stats[!is.na(puma_birth_rate)]

# Step 2: Calculate the weighted standard deviation of birth rates for each MSA.
weighted_sd <- function(x, w, na.rm = TRUE) {
  if (na.rm) {
    w <- w[!is.na(x)]
    x <- x[!is.na(x)]
  }
  weighted_mean <- sum(w * x) / sum(w)
  sqrt(sum(w * (x - weighted_mean)^2) / sum(w))
}

# Apply the function to each MSA-year group.
# We also count PUMAs; inequality is meaningless with only one.
msa_inequality <- puma_stats[, .(
  fertility_rate_inequality = weighted_sd(puma_birth_rate, puma_population),
  num_pumas = uniqueN(puma)
), by = .(met2013, year)]

# Inequality is 0 if there's only one PUMA.
msa_inequality[num_pumas <= 1, fertility_rate_inequality := 0]


# ---- 5. Add MSA Names and Save Output ----
cat("--- Adding MSA names to inequality data... ---\n")
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
msa_inequality[, met2013 := as.character(met2013)]
msa_inequality <- merge(msa_inequality, msa_name_lookup, by = "met2013", all.x = TRUE)
setcolorder(msa_inequality, c("met2013", "msa_name", "year"))

fwrite(msa_inequality, output_file)
cat(paste0("--- Inequality data saved to ", output_file, " ---\n"))


# ---- 6. Generate and Save a Map ----
cat("--- Generating a map of fertility rate inequality... ---\n")
states_geom <- tigris::states(cb = TRUE) %>% sf::st_transform(sf::st_crs(msa_geometries_for_names))
latest_year <- max(msa_inequality$year, na.rm = TRUE)
map_data_sf <- merge(
  msa_geometries_for_names,
  msa_inequality[year == latest_year],
  by.x = "GEOID",
  by.y = "met2013"
) %>% sf::st_as_sf()

inequality_map <- ggplot() +
  geom_sf(data = states_geom, fill = "gray90", color = "white", linewidth = 0.5) +
  geom_sf(data = map_data_sf, aes(fill = fertility_rate_inequality), color = "white", linewidth = 0.1) +
  coord_sf(crs = sf::st_crs(5070), xlim = c(-2500000, 2500000), ylim = c(-2300000, 3000000)) +
  scale_fill_viridis_c(
    option = "rocket", name = "Fertility Rate Inequality\\n(Weighted St. Dev.)",
    guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.position = "top")
  ) +
  labs(
    title = "Fertility Rate Inequality Within U.S. Metropolitan Areas",
    subtitle = paste("Based on PUMA-level birth rates for", latest_year),
    caption = "Source: ACS Microdata. High values (yellow) indicate large variation in birth rates between a city's neighborhoods."
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    plot.caption = element_text(hjust = 0.5)
  )

ggsave(output_map_file, plot = inequality_map, width = 11, height = 8.5, dpi = 300, bg = "white")
cat(paste0("--- Map saved to ", output_map_file, " ---\n"))
cat("--- Script complete. ---\n") 