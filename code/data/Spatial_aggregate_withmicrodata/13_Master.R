# ======================================================================================
#
#    Microdata-Driven Fertility Analysis - V9.0 (LaTeX Export Refactor)
#    Author: Gemini AI Assistant (as programming partner)
#    Date: June 11, 2025
#
#    --- MODIFICATIONS ---
#    1. Refactored all LaTeX table exports to use robust helper functions.
#       - `save_table_latex` for fixest models, using etable and manual wrappers.
#       - `save_df_latex` for data frames (like tidy model summaries), using kableExtra.
#    2. Added the "Changes-on-Changes" analysis and its corresponding table.
#    3. Removed redundant/conflicting table export sections.
#
# ======================================================================================

# ======================================================================================
# 0. SETUP AND CONFIGURATION
# ======================================================================================
cat("--- 0. SETUP AND CONFIGURATION ---\n")

# ---- Function to install and load packages ----
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
  }
}

# ---- List of required packages ----
required_packages <- c(
  "tidyverse", "data.table", "fixest", "modelsummary", "marginaleffects",
  "sf", "tigris", "viridis", "patchwork", "scales", "spatstat.geom", "broom", "kableExtra"
)
install_and_load(required_packages)

# ---- Set global options for table exports ----
# Prevent etable from using siunitx for maximum LaTeX compatibility
options("fixest_etable_tex_siunitx" = FALSE)
# Force modelsummary to use the kableExtra backend instead of tabularray
options(modelsummary_factory_latex = "kableExtra")

# ---- Set up directories and file paths ----
# Use `here` for robust path management
if (!require(here)) { install.packages("here"); library(here) }
OUTPUT_DIR <- here("Spatial_aggregate_withmicrodata", "analysis_output_v7")
LATEX_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LATEX_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helper functions for saving outputs ----
save_plot <- function(plot_obj, filename, width = 10, height = 6.5, dpi = 300) {
  full_path <- file.path(OUTPUT_DIR, filename)
  dir.create(dirname(full_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(full_path, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("Saved plot: %s\n", full_path))
}
save_table_html <- function(table_obj, filename) {
  full_path <- file.path(OUTPUT_DIR, filename)
  if (is.data.frame(table_obj) || is.matrix(table_obj)) {
    write.csv(as.data.frame(table_obj), full_path, row.names = FALSE)
  } else {
    modelsummary::modelsummary(table_obj, output = full_path)
  }
  cat(sprintf("Saved HTML table: %s\n", full_path))
}

save_table_latex <- function(model_list, filename, title, label) {
  # This version uses capture.output to create a self-contained, robust table
  # that doesn't rely on complex LaTeX packages that might be missing.
  table_content <- capture.output(
    etable(model_list,
           signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
           fitstat = ~ n + r2 + ar2, # Using a robust set of stats
           tex = TRUE)
  )

  # Manually construct the full table float for maximum compatibility
  full_latex_code <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\caption{", title, "}"),
    paste0("\\label{", label, "}"),
    table_content,
    "\\end{table}"
  )

  full_path <- file.path(LATEX_DIR, filename)
  writeLines(full_latex_code, full_path)
  cat(sprintf("Saved robust LaTeX table: %s\\n", full_path))
}

save_df_latex <- function(df, filename, title, label) {
  # This uses kableExtra, which is more reliable than manually creating a table from a df.
  latex_table <- kableExtra::kbl(df,
                                 format = "latex",
                                 booktabs = TRUE,
                                 caption = title,
                                 label = label,
                                 digits = 3,
                                 align = "r") %>%
                 kableExtra::kable_styling(latex_options = "hold_position", font_size = 9)

  full_path <- file.path(LATEX_DIR, filename)
  writeLines(latex_table, full_path)
  cat(sprintf("Saved DataFrame as LaTeX table: %s\n", full_path))
}

# ---- Plotting Theme ----
theme_set(theme_minimal(base_size = 12) +
            theme(
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray30"),
              plot.caption = element_text(hjust = 0, color = "gray50"),
              axis.title = element_text(size = 11, face="bold"),
              legend.position = "bottom",
              strip.text = element_text(face="bold", size=11)
            ))

cat("Setup complete. Outputs ->", normalizePath(OUTPUT_DIR), "\n")


# ======================================================================================
# 1. DATA PREPARATION (Based on Cluster Script Logic)
# ======================================================================================
cat("\n--- 1. DATA LOADING AND PREPARATION ---\n")

# --- Configuration ---
# Use `here` for robust path management
if (!require(here)) { install.packages("here"); library(here) }

# --- Define File Paths ---
# Path to the final, cleaned, and sampled microdata file from script 03
INPUT_DATA_PATH <- here("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_10pct_stratified_sample_v7.rds")

# --- Analysis Parameters ---
MIN_MSA_OBS <- 100 # Minimum number of individual observations to include an MSA in a given year

if (!file.exists(INPUT_DATA_PATH)) { stop(paste("Data file not found:", INPUT_DATA_PATH)) }
acs_data <- readRDS(INPUT_DATA_PATH)
setDT(acs_data)
cat("  1.1 Engineering base features...\n")
acs_data[, educ_factor := factor(educ_factor, levels = c("Less than High School", "High School Grad", "Some College", "Bachelors", "Graduate"))]
acs_data[, age_sq := age^2]
acs_data[, met2013 := as.factor(met2013)]
acs_data[, age_group := cut(age, breaks = c(14, 24, 34, 44, 51), labels = c("15-24", "25-34", "35-44", "45-50"), right = FALSE)]
acs_data[, did_relocate := as.integer(migrate1d %in% c(21, 22, 23, 24))]
cat("  1.2 Calculating individual-level housing cost metrics...\n")
acs_data[hhincome <= 0, hhincome := NA]
acs_data[, individual_pti := valueh / hhincome]
acs_data[valueh == 9999999, individual_pti := NA]
acs_data[individual_pti <= 0 | individual_pti > 50, individual_pti := NA]
acs_data[, individual_rti := (rent * 12) / hhincome]
acs_data[rent <= 0, individual_rti := NA]
acs_data[individual_rti <= 0 | individual_rti > 1, individual_rti := NA]
cat("  1.3 Building the aggregate MSA-level panel...\n")
msa_panel <- acs_data[met2013 != "0" & hhincome > 0,
                      .(
                        mean_age_first_birth = weighted.mean(age_at_first_birth, perwt, na.rm = TRUE),
                        fertility_rate = weighted.mean(is_recent_mother, perwt, na.rm = TRUE),
                        pct_hh_with_children = weighted.mean(nchild > 0, hhwt, na.rm = TRUE),
                        mean_nchild = weighted.mean(nchild, hhwt, na.rm = TRUE),
                        msa_pti = tryCatch(weighted.median(individual_pti, hhwt, na.rm = TRUE), error = function(e) NA_real_),
                        msa_rti = tryCatch(weighted.median(individual_rti, hhwt, na.rm = TRUE), error = function(e) NA_real_),
                        msa_population = sum(perwt, na.rm = TRUE),
                        msa_median_hh_income = weighted.median(hhincome, hhwt, na.rm = TRUE),
                        pct_grad_plus = weighted.mean(educ_factor == "Graduate", perwt, na.rm = TRUE)
                      ), by = .(met2013, year)]

cat("  1.3b Calculating MSA-level age composition...\n")
# Define meaningful age brackets
acs_data[, age_bracket := cut(age,
                              breaks = c(0, 18, 25, 35, 45, 65, Inf),
                              labels = c("0-17", "18-24", "25-34", "35-44", "45-64", "65+"),
                              right = FALSE)]
# Aggregate by MSA, year, and age bracket
msa_age_counts <- acs_data[!is.na(age_bracket) & met2013 != "0",
                           .(pop_count = sum(perwt)), by = .(met2013, year, age_bracket)]
# Reshape from long to wide to get percentages
msa_age_comp <- dcast(msa_age_counts, met2013 + year ~ age_bracket, value.var = "pop_count", fill = 0)
# Calculate total population and percentages
msa_age_comp[, total_pop := `0-17` + `18-24` + `25-34` + `35-44` + `45-64` + `65+`]
setnames(msa_age_comp,
         c("0-17", "18-24", "25-34", "35-44", "45-64", "65+"),
         c("pct_pop_0_17", "pct_pop_18_24", "pct_pop_25_34", "pct_pop_35_44", "pct_pop_45_64", "pct_pop_65_plus"))
# Convert counts to percentages
cols_to_pct <- c("pct_pop_0_17", "pct_pop_18_24", "pct_pop_25_34", "pct_pop_35_44", "pct_pop_45_64", "pct_pop_65_plus")
msa_age_comp[, (cols_to_pct) := lapply(.SD, function(x) x / total_pop), .SDcols = cols_to_pct]
# Calculate median age and dependency ratios
msa_age_summary <- acs_data[met2013 != "0", .(
  median_age = weighted.median(age, perwt, na.rm=TRUE),
  youth_dependency = sum(perwt[age < 18], na.rm=T) / sum(perwt[age %between% c(18,64)], na.rm=T),
  aged_dependency = sum(perwt[age >= 65], na.rm=T) / sum(perwt[age %between% c(18,64)], na.rm=T)
), by = .(met2013, year)]
# Merge age composition and summary stats into the main MSA panel
msa_panel <- merge(msa_panel, msa_age_comp[, -c("total_pop")], by = c("met2013", "year"), all.x = TRUE)
msa_panel <- merge(msa_panel, msa_age_summary, by = c("met2013", "year"), all.x = TRUE)

cat("  1.4 Fetching MSA land area to calculate population density...\n")
options(tigris_use_cache = TRUE)
msa_geometries <- tryCatch({tigris::core_based_statistical_areas(cb = TRUE, year = 2022)}, error = function(e) {
  cat("  -> Failed to get 2022 MSA geometries, trying 2020...\n")
  tigris::core_based_statistical_areas(cb = TRUE, year = 2020)
})
land_area <- as.data.table(msa_geometries)[, .(GEOID, land_area_sq_miles = ALAND / 2589988.11)]
msa_panel <- merge(msa_panel, land_area, by.x = "met2013", by.y = "GEOID", all.x = TRUE)
msa_panel[, msa_pop_density := msa_population / land_area_sq_miles]
msa_panel_clean <- msa_panel[is.finite(msa_pti) & is.finite(msa_rti) & is.finite(mean_age_first_birth) & is.finite(msa_pop_density) & msa_pop_density > 0 & is.finite(msa_median_hh_income) & msa_median_hh_income > 0 & is.finite(mean_nchild)]
msa_panel_clean[, log_msa_pop_density := log(msa_pop_density)]
msa_panel_clean[, log_msa_median_hh_income := log(msa_median_hh_income)]
cat("  1.5 Merging aggregate metrics into microdata for individual analysis...\n")
analysis_data <- merge(acs_data, msa_panel_clean[, .(met2013, year, msa_pti, msa_rti, msa_pop_density, msa_population, log_msa_pop_density)], by = c("met2013", "year"), all.x = TRUE)
cat("Data preparation complete.\n")


# ======================================================================================
# 2. TRACK 1: AGGREGATE (MSA-LEVEL) ANALYSIS
# ======================================================================================
cat("\n--- 2. TRACK 1: AGGREGATE (MSA-LEVEL) ANALYSIS ---\n")
latest_year_data_agg <- msa_panel_clean[year == max(year)]

# ======================================================================================
# 2.1 AGGREGATE-LEVEL SCATTERPLOTS
# ======================================================================================
cat("  2.1 Generating aggregate-level plots...\n")
p_agg_pti <- ggplot(latest_year_data_agg, aes(x = msa_pti, y = fertility_rate)) +
  geom_point(alpha = 0.6, aes(color = log(msa_median_hh_income))) + geom_smooth(method = "lm", formula = y ~ x, color = "darkred") +
  scale_color_viridis_c(name = "Log Median\nHH Income") +
  labs(title = "MSA Fertility vs. Price-to-Income Ratio", subtitle = paste("Latest Year:", max(latest_year_data_agg$year)), x = "MSA Price-to-Income Ratio", y = "MSA Fertility Rate")
save_plot(p_agg_pti, "AGG_01_fertility_vs_pti.png")
p_agg_pti_scatter <- ggplot(latest_year_data_agg, aes(x = msa_pti, y = fertility_rate)) +
  geom_point(alpha = 0.7, aes(size = msa_population, color = log(msa_median_hh_income))) +
  geom_smooth(method = "lm", formula = y ~ x, color = "darkred", se = FALSE) +
  scale_color_viridis_c(name = "Log Median\nHH Income") +
  scale_size_continuous(name = "MSA Population", labels = scales::comma) +
  labs(title = "MSA Fertility vs. Price-to-Income Ratio (by Population)", subtitle = paste("Latest Year:", max(latest_year_data_agg$year)), x = "MSA Price-to-Income Ratio", y = "MSA Fertility Rate")
save_plot(p_agg_pti_scatter, "AGG_01b_fertility_vs_pti_scatter_pop.png")
p_agg_rti <- ggplot(latest_year_data_agg, aes(x = msa_rti, y = mean_age_first_birth)) +
  geom_point(alpha = 0.6, aes(color = pct_grad_plus)) + geom_smooth(method = "lm", formula = y ~ x, color = "darkblue") +
  scale_color_viridis_c(name = "Pct with\nGrad Degree", option="magma") +
  labs(title = "MSA Age at First Birth vs. Rent-to-Income Ratio", subtitle = paste("Latest Year:", max(latest_year_data_agg$year)), x = "MSA Rent-to-Income Ratio", y = "Mean Age at First Birth")
save_plot(p_agg_rti, "AGG_02_age_vs_rti.png")
p_agg_rti_scatter <- ggplot(latest_year_data_agg, aes(x = msa_rti, y = mean_age_first_birth)) +
  geom_point(alpha = 0.7, aes(size = msa_population, color = pct_grad_plus)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "darkblue", se = FALSE) +
  scale_color_viridis_c(name = "Pct with\nGrad Degree", option="magma") +
  scale_size_continuous(name = "MSA Population", labels = scales::comma) +
  labs(title = "MSA Age at First Birth vs. Rent-to-Income Ratio (by Population)", subtitle = paste("Latest Year:", max(latest_year_data_agg$year)), x = "MSA Rent-to-Income Ratio", y = "Mean Age at First Birth")
save_plot(p_agg_rti_scatter, "AGG_02b_age_vs_rti_scatter_pop.png")
p_agg_nchild_pti <- ggplot(latest_year_data_agg, aes(x = msa_pti, y = mean_nchild)) +
  geom_point(alpha = 0.6, aes(color = log(msa_median_hh_income))) + geom_smooth(method = "lm", formula = y ~ x, color = "darkgreen") +
  scale_color_viridis_c(name = "Log Median\nHH Income", option = "mako") +
  labs(title = "MSA Mean Children per Household vs. Price-to-Income Ratio", subtitle = paste("Latest Year:", max(latest_year_data_agg$year)), x = "MSA Price-to-Income Ratio", y = "Mean Children per Household")
save_plot(p_agg_nchild_pti, "AGG_02c_nchild_vs_pti.png")
p_agg_nchild_pti_scatter <- ggplot(latest_year_data_agg, aes(x = msa_pti, y = mean_nchild)) +
  geom_point(alpha = 0.7, aes(size = msa_population, color = log(msa_median_hh_income))) +
  geom_smooth(method = "lm", formula = y ~ x, color = "darkgreen", se = FALSE) +
  scale_color_viridis_c(name = "Log Median\nHH Income", option = "mako") +
  scale_size_continuous(name = "MSA Population", labels = scales::comma) +
  labs(title = "MSA Mean Children vs. Price-to-Income Ratio (by Population)", subtitle = paste("Latest Year:", max(latest_year_data_agg$year)), x = "MSA Price-to-Income Ratio", y = "Mean Children per Household")
save_plot(p_agg_nchild_pti_scatter, "AGG_02d_nchild_vs_pti_scatter_pop.png")

# ======================================================================================
# 2.2 ANALYSIS OF MSA AGE COMPOSITION
# ======================================================================================
cat("\n  2.2 Analyzing MSA Age Composition...\n")
age_comp_dir <- file.path(OUTPUT_DIR, "AGG_age_composition")
dir.create(age_comp_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2.2.1 National Trend: Stacked Area Chart ---
cat("    -> Plotting national age composition trends...\n")
national_age_trends <- msa_panel_clean[, lapply(.SD, weighted.mean, w = msa_population, na.rm = TRUE),
                                       .SDcols = c("median_age", cols_to_pct), by = year]

national_age_long <- melt(national_age_trends, id.vars = "year", measure.vars = cols_to_pct,
                          variable.name = "age_group", value.name = "percentage")

p_national_stacked_area <- ggplot(national_age_long, aes(x = year, y = percentage, fill = age_group)) +
  geom_area(alpha = 0.8, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(name = "Age Group", option = "D") +
  labs(
    title = "National Urban Age Composition Over Time",
    subtitle = "Population-weighted average across all MSAs",
    x = "Year", y = "Share of Population"
  )
save_plot(p_national_stacked_area, file.path("AGG_age_composition", "AGE_01_national_stacked_area.png"))

# --- 2.2.2 Spaghetti Plot for Top MSAs ---
cat("    -> Plotting median age trends for largest MSAs...\n")
top_msas <- msa_panel_clean[year == max(year)][order(-msa_population)][1:15, met2013]
p_top_msa_median_age <- ggplot(msa_panel_clean[met2013 %in% top_msas],
                               aes(x = year, y = median_age, group = met2013, color = met2013)) +
  geom_line(alpha = 0.8, linewidth = 1) +
  scale_color_viridis_d(guide = "none") + # Hide messy legend
  labs(
    title = "Median Age Trend in 15 Largest MSAs",
    subtitle = "Most large urban areas are aging, with some variation",
    x = "Year", y = "MSA Median Age"
  )
save_plot(p_top_msa_median_age, file.path("AGG_age_composition", "AGE_02_top_msa_median_age.png"))

# --- 2.2.3 Cross-sectional Relationship with Density and Cost (Latest Year) ---
cat("    -> Plotting cross-sectional relationships for latest year...\n")
p_age_vs_density <- ggplot(latest_year_data_agg, aes(x = log_msa_pop_density, y = median_age)) +
  geom_point(aes(color = log_msa_median_hh_income, size = msa_population), alpha = 0.7) +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE) +
  scale_color_viridis_c(name = "Log Median\nHH Income") +
  scale_size_continuous(name = "Population", labels = scales::comma) +
  labs(
    title = "Denser MSAs Tend to Have Younger Populations",
    subtitle = paste("Latest Year:", max(latest_year_data_agg$year)),
    x = "Log(Population Density)", y = "Median Age"
  )

p_fertility_vs_age <- ggplot(latest_year_data_agg, aes(x = median_age, y = fertility_rate)) +
  geom_point(aes(color = msa_pti, size = msa_population), alpha = 0.7) +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE) +
  scale_color_viridis_c(name = "Price-to-Income", option="plasma") +
  scale_size_continuous(name = "Population", labels = scales::comma) +
  labs(
    title = "Fertility is Lower in 'Older' MSAs",
    subtitle = paste("Latest Year:", max(latest_year_data_agg$year)),
    x = "Median Age", y = "Fertility Rate"
  )

combined_age_scatter <- p_age_vs_density + p_fertility_vs_age
save_plot(combined_age_scatter, file.path("AGG_age_composition", "AGE_03_age_scatters.png"), width = 14, height = 6)


# ======================================================================================
# 2.3 ANALYSIS OF DENSITY AND AGE TRENDS
# ======================================================================================
cat("\n  2.3 Analyzing if Denser Cities are Getting Younger Faster...\n")
trends_dir <- file.path(OUTPUT_DIR, "AGG_trend_analysis")
dir.create(trends_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2.3.1 Visualization: Age Trends by Density Quartile ---
# To avoid a city changing category over time, we define its density based on the first year of data.
cat("    -> Visualizing age trends by baseline density quartile...\n")
baseline_density <- msa_panel_clean[year == min(year), .(met2013, baseline_log_density = log_msa_pop_density)]

# Create density quartiles
density_breaks <- quantile(baseline_density$baseline_log_density, probs = 0:4/4, na.rm = TRUE)
baseline_density[, density_quartile := cut(baseline_log_density,
                                           breaks = density_breaks,
                                           labels = c("Q1 (Lowest Density)", "Q2", "Q3", "Q4 (Highest Density)"),
                                           include.lowest = TRUE)]

# Merge the static density quartile back into the main panel
msa_panel_with_trends <- merge(msa_panel_clean, baseline_density[, .(met2013, density_quartile)], by = "met2013")

# Calculate the population-weighted average median age for each quartile over time
age_trends_by_density <- msa_panel_with_trends[!is.na(density_quartile),
                                               .(avg_median_age = weighted.mean(median_age, msa_population, na.rm = TRUE)),
                                               by = .(year, density_quartile)]

# Plot the trends
p_age_trend_by_density <- ggplot(age_trends_by_density, aes(x = year, y = avg_median_age, color = density_quartile)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_viridis_d(name = "Baseline Density Quartile") +
  labs(
    title = "Age Trend by MSA Population Density",
    subtitle = "Highest-density cities have aged slower than lower-density cities",
    x = "Year",
    y = "Weighted Average of MSA Median Age",
    caption = "Density quartiles are fixed based on each MSA's density in the first year of the sample."
  )
save_plot(p_age_trend_by_density, file.path("AGG_trend_analysis", "TREND_01_age_by_density_quartile.png"))


# --- 2.3.2 Regression: Formal Test of the Interaction ---
cat("    -> Running regression model with a time-density interaction...\n")
msa_panel_for_reg <- merge(msa_panel_clean, baseline_density[, .(met2013, baseline_log_density)], by = "met2013")
m_age_trend_interaction <- feols(median_age ~ year * baseline_log_density, data = msa_panel_for_reg, fixef = "met2013")
m_age_trend_interaction_controlled <- feols(median_age ~ year * baseline_log_density + log_msa_median_hh_income + pct_grad_plus, data = msa_panel_for_reg, fixef = "met2013")
save_table_html(list("Simple" = m_age_trend_interaction, "Controlled" = m_age_trend_interaction_controlled), "AGG_TREND_02_age_trend_interaction.html")
save_table_latex(
    list("Simple"=m_age_trend_interaction, "Controlled"=m_age_trend_interaction_controlled),
    filename = "agg_age_trend_interaction_model.tex",
    title = "Interaction of Time and MSA Population Density on Median Age",
    label = "tab:age_density_interaction"
)


# ======================================================================================
# 2.4 AGGREGATE PANEL REGRESSIONS
# ======================================================================================
cat("\n  2.4 Running aggregate panel regressions...\n")
agg_controls <- "log_msa_median_hh_income + pct_grad_plus + log_msa_pop_density + median_age + pct_pop_25_34"
agg_models <- list()
for(outcome in c("mean_age_first_birth", "fertility_rate", "pct_hh_with_children", "mean_nchild")){
  for(metric in c("msa_pti", "msa_rti")){
    agg_models[[paste(outcome, metric, "simple", sep="_")]] <- feols(as.formula(paste(outcome, "~", metric)), data = msa_panel_clean, fixef = "year")
    agg_models[[paste(outcome, metric, "controlled", sep="_")]] <- feols(as.formula(paste(outcome, "~", metric, "+", agg_controls)), data = msa_panel_clean, fixef = c("met2013", "year"))
  }
}
save_table_html(agg_models, "AGG_03_panel_regressions.html")
save_table_latex(
    agg_models,
    filename = "agg_panel_regressions.tex",
    title = "Aggregate Panel Regression Models of Fertility and Housing Costs",
    label = "tab:panel_regressions"
)


# ======================================================================================
# 2.5 TIME-VARYING EFFECTS AT THE AGGREGATE LEVEL
# ======================================================================================
cat("  2.5 Analyzing time-varying effects at the aggregate level...\n")
agg_yearly_results <- list()
agg_outcomes <- c("fertility_rate", "mean_age_first_birth", "pct_hh_with_children", "mean_nchild")
agg_metrics <- c("msa_pti", "msa_rti", "log_msa_pop_density")
min_obs_required <- 5
for (outcome in agg_outcomes) {
  for (metric in agg_metrics) {
    for (yr in sort(unique(msa_panel_clean$year))) {
      yr_data <- msa_panel_clean[year == yr]
      if (sum(!is.na(yr_data[[metric]])) < min_obs_required) next
      m_simple <- lm(as.formula(paste(outcome, "~", metric)), data = yr_data)
      agg_yearly_results <- append(agg_yearly_results, list(tidy(m_simple) %>% filter(term==metric) %>% mutate(year=yr, outcome=outcome, metric=metric, controls="Simple", nobs=nobs(m_simple))))
      controls_to_use <- if (metric == "log_msa_pop_density") "log_msa_median_hh_income + pct_grad_plus + median_age" else "log_msa_pop_density + pct_grad_plus + median_age"
      m_controlled <- lm(as.formula(paste(outcome, "~", metric, "+", controls_to_use)), data = yr_data)
      agg_yearly_results <- append(agg_yearly_results, list(tidy(m_controlled) %>% filter(term==metric) %>% mutate(year=yr, outcome=outcome, metric=metric, controls="With Controls", nobs=nobs(m_controlled))))
    }
  }
}
agg_yearly_coefs <- rbindlist(agg_yearly_results)
save_table_html(agg_yearly_coefs, "AGG_04_yearly_coefficients_data.csv")
cat("  2.5.1 Saving individual time-varying plots...\n")
for (o in unique(agg_yearly_coefs$outcome)) {
  for (m in unique(agg_yearly_coefs$metric)) {
    plot_data <- agg_yearly_coefs[outcome == o & metric == m]
    if (nrow(plot_data) > 0) {
      p <- ggplot(plot_data, aes(x = year, y = estimate, color = controls, fill = controls)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_ribbon(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), alpha = 0.2, linetype = 0) +
        geom_line(linewidth = 1) +
        labs(title = paste("Time-Varying Effect of", m, "on", o), subtitle = "Coefficient from Yearly MSA-Level Cross-Sectional Regression", x = "Year", y = "Coefficient Estimate", color = "Model", fill = "Model")
      file_name <- paste0("AGG_05_", o, "_vs_", m, ".png")
      save_plot(p, file.path("AGG_05_individual_plots", file_name), width = 8, height = 6)
    }
  }
}
p_agg_time_varying <- ggplot(agg_yearly_coefs, aes(x = year, y = estimate, color = controls, fill = controls)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 1) +
  facet_grid(outcome~metric, scales = "free_y", labeller = "label_both") +
  labs(title = "Time-Varying Effects at the MSA-Level", subtitle = "Coefficients from Yearly Cross-Sectional Regressions", x = "Year", y = "Coefficient Estimate", color = "Model", fill = "Model")
save_plot(p_agg_time_varying, "AGG_05_time_varying_coefficients.png", width=12, height=10)


# ======================================================================================
# 2.6 LATEST-YEAR SCATTERPLOTS AND REGRESSION TABLE
# ======================================================================================
cat("\n  2.6 Generating latest-year cross-sectional scatters and table...\n")
latest_year_models <- list()
scatter_dir <- file.path(OUTPUT_DIR, "AGG_06_yearly_scatters")
dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)
for (outcome in agg_outcomes) {
  for (metric in agg_metrics) {
    p_scatter <- ggplot(latest_year_data_agg, aes(x = .data[[metric]], y = .data[[outcome]])) +
      geom_point(aes(size = msa_population), alpha = 0.6, color = "darkcyan") +
      geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
      scale_size_continuous(name = "MSA Population", labels = scales::comma) +
      labs(title = paste("Latest Year Relationship:", outcome, "vs.", metric), subtitle = paste("Data from", max(latest_year_data_agg$year), "| Each point is an MSA"), x = metric, y = outcome)
    scatter_filename <- paste0("scatter_", outcome, "_vs_", metric, ".png")
    save_plot(p_scatter, file.path("AGG_06_yearly_scatters", scatter_filename), width = 8, height = 6)
    model_formula <- as.formula(paste(outcome, "~", metric))
    latest_year_model <- lm(model_formula, data = latest_year_data_agg)
    model_name <- paste(outcome, "vs", metric)
    latest_year_models[[model_name]] <- latest_year_model
  }
}
cat("  -> Saving summary table for latest-year regressions...\\n")
modelsummary(
  latest_year_models,
  stars = TRUE,
  gof_map = c("nobs", "r.squared"),
  title = paste("MSA-Level Cross-Sectional Regressions for Latest Year (", max(latest_year_data_agg$year), ")", sep=""),
  output = file.path(LATEX_DIR, "agg_latest_year_regressions.tex")
)
save_table_html(latest_year_models, "AGG_07_latest_year_regressions.html")


# ======================================================================================
# 2.7 Changes-on-Changes Analysis
# ======================================================================================
cat("\n--- 2.7 Analyzing Long-Term Changes in Fertility vs. Age ---\\n")

start_year <- min(msa_panel_clean$year)
end_year <- max(msa_panel_clean$year)
change_data <- msa_panel_clean[year %in% c(start_year, end_year), .(met2013, year, fertility_rate, median_age, msa_population)]
change_data_wide <- dcast(change_data, met2013 ~ year, value.var = c("fertility_rate", "median_age", "msa_population"))
setnames(change_data_wide, old = paste0("fertility_rate_", start_year), new = "fertility_rate_start", skip_absent = TRUE)
setnames(change_data_wide, old = paste0("fertility_rate_", end_year), new = "fertility_rate_end", skip_absent = TRUE)
setnames(change_data_wide, old = paste0("median_age_", start_year), new = "median_age_start", skip_absent = TRUE)
setnames(change_data_wide, old = paste0("median_age_", end_year), new = "median_age_end", skip_absent = TRUE)
setnames(change_data_wide, old = paste0("msa_population_", start_year), new = "msa_population_start", skip_absent = TRUE)
setnames(change_data_wide, old = paste0("msa_population_", end_year), new = "msa_population_end", skip_absent = TRUE)
change_data_wide[, delta_fertility := fertility_rate_end - fertility_rate_start]
change_data_wide[, delta_median_age := median_age_end - median_age_start]
change_data_wide <- change_data_wide[is.finite(delta_fertility) & is.finite(delta_median_age)]
m_change_on_change <- lm(delta_fertility ~ delta_median_age, data = change_data_wide, weights = msa_population_start)
change_models <- list("Change-on-Change" = m_change_on_change)
change_title <- "Long-Term Change in Fertility vs. Change in Median Age"
tidy_change_model <- broom::tidy(m_change_on_change)
save_df_latex(
    df = tidy_change_model,
    filename = "CHANGES_01_regression.tex",
    title = change_title,
    label = "tab:changes_regression"
)
save_table_html(change_models, "CHANGES_01_regression.html")
p_change_on_change <- ggplot(change_data_wide, aes(x = delta_median_age, y = delta_fertility)) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray50") +
  geom_vline(xintercept = mean(change_data_wide$delta_median_age, na.rm=T), linetype="dashed", color="gray50") +
  geom_point(aes(size = msa_population_start, color = delta_fertility), alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, color = "firebrick", aes(weight = msa_population_start)) +
  scale_size_continuous(name = "MSA Population (Start Year)", labels = scales::comma) +
  scale_color_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0, name = "Change in GFR") +
  labs(title = "Relationship Between Long-Term Changes in Age and Fertility",
       subtitle = paste("Change from", start_year, "to", end_year, ". Each point is an MSA."),
       x = "Change in MSA Median Age",
       y = "Change in General Fertility Rate (GFR)",
       caption = "Regression line is weighted by starting year population.")
save_plot(p_change_on_change, "AGG_08_changes_on_changes_plot.png", width = 11, height = 7)
cat("  -> Changes-on-Changes analysis complete.\\n")


# ======================================================================================
# 3. TRACK 2: INDIVIDUAL-LEVEL ANALYSIS
# ======================================================================================
cat("\n--- 3. TRACK 2: INDIVIDUAL-LEVEL ANALYSIS ---\n")
indiv_data_clean <- analysis_data[age %between% c(15, 50) & met2013 != "0" & !is.na(msa_pti) & !is.na(msa_rti) & !is.na(educ_factor)]
cat("  3.1 Running main individual-level regressions...\n")
indiv_controls <- "age + age_sq + educ_factor"
indiv_models <- list(
  "Fertility ~ PTI (Simple)" = feglm(is_recent_mother ~ msa_pti | met2013 + year, data=indiv_data_clean, family="binomial", weights=~perwt),
  "Fertility ~ PTI (Controlled)" = feglm(as.formula(paste("is_recent_mother ~ msa_pti +", indiv_controls)), data=indiv_data_clean, family="binomial", fixef=c("met2013", "year"), weights=~perwt),
  "Fertility ~ RTI (Simple)" = feglm(is_recent_mother ~ msa_rti | met2013 + year, data=indiv_data_clean, family="binomial", weights=~perwt),
  "Fertility ~ RTI (Controlled)" = feglm(as.formula(paste("is_recent_mother ~ msa_rti +", indiv_controls)), data=indiv_data_clean, family="binomial", fixef=c("met2013", "year"), weights=~perwt)
)
save_table_html(indiv_models, "IND_01_main_logit_models.html")
save_table_latex(
    indiv_models,
    filename = "indiv_main_logit_models.tex",
    title = "Individual-Level Logit Models of Recent Birth",
    label = "tab:indiv_main_logit"
)

cat("  3.2 Analyzing heterogeneity with interaction models...\n")
int_model_educ_pti <- feglm(is_recent_mother ~ educ_factor * msa_pti + age + age_sq | met2013 + year, data = indiv_data_clean, family = "binomial", weights = ~perwt)
int_model_age_pti <- feglm(is_recent_mother ~ age_group * msa_pti + age + age_sq + educ_factor | met2013 + year, data = indiv_data_clean, family = "binomial", weights = ~perwt)
int_model_age_rti <- feglm(is_recent_mother ~ age_group * msa_rti + age + age_sq + educ_factor | met2013 + year, data = indiv_data_clean, family = "binomial", weights = ~perwt)
save_table_html(
    list("PTI × Education" = int_model_educ_pti, "PTI × Age Group"  = int_model_age_pti, "RTI × Age Group"  = int_model_age_rti),
    "IND_02_interaction_models_PTI_RTI.html"
)
save_table_latex(
    list("PTI × Education" = int_model_educ_pti, "PTI × Age Group"  = int_model_age_pti, "RTI × Age Group"  = int_model_age_rti),
    filename = "indiv_interaction_models_pti_rti.tex",
    title = "Individual-Level Logit Models with Interactions",
    label = "tab:indiv_interaction_pti_rti"
)
mfx_educ_plot <- plot_predictions(int_model_educ_pti, condition = c("msa_pti", "educ_factor")) + labs(title = "Effect of Housing Price on Fertility by Education", subtitle = paste("N =", nobs(int_model_educ_pti)), x = "MSA Price-to-Income Ratio")
save_plot(mfx_educ_plot, "IND_03_mfx_pti_by_education.png")
mfx_age_plot <- plot_predictions(int_model_age_pti, condition = c("msa_pti", "age_group")) + labs(title = "Effect of Housing Price on Fertility by Age Group", subtitle = paste("N =", nobs(int_model_age_pti)), x = "MSA Price-to-Income Ratio")
save_plot(mfx_age_plot, "IND_04_mfx_pti_by_age_group.png")
mfx_age_rti_plot <- plot_predictions(int_model_age_rti, condition = c("msa_rti", "age_group")) + labs(title = "Effect of Rent Burden on Fertility by Age Group", subtitle = paste("N =", nobs(int_model_age_rti)), x = "MSA Rent-to-Income Ratio")
save_plot(mfx_age_rti_plot, "IND_05_mfx_rti_by_age_group.png")


cat("  3.3 Analyzing time-varying effects at the individual level...\n")
indiv_yearly_results <- list()
for (metric in c("msa_pti", "msa_rti")) {
  for (yr in sort(unique(indiv_data_clean$year))) {
    yr_data <- indiv_data_clean[year == yr]
    m_simple <- safely(feglm)(as.formula(paste("is_recent_mother ~", metric)), data=yr_data, family="binomial", fixef="met2013", weights=~perwt)
    if(!is.null(m_simple$result)) indiv_yearly_results <- append(indiv_yearly_results, list(tidy(m_simple$result) %>% filter(grepl(metric, term)) %>% mutate(year=yr, metric=metric, controls="Simple", nobs=nobs(m_simple$result))))
    m_controlled <- safely(feglm)(as.formula(paste("is_recent_mother ~", metric, "+", indiv_controls)), data=yr_data, family="binomial", fixef="met2013", weights=~perwt)
    if(!is.null(m_controlled$result)) indiv_yearly_results <- append(indiv_yearly_results, list(tidy(m_controlled$result) %>% filter(grepl(metric, term)) %>% mutate(year=yr, metric=metric, controls="With Controls", nobs=nobs(m_controlled$result))))
  }
}
indiv_yearly_coefs <- rbindlist(indiv_yearly_results)
save_table_html(indiv_yearly_coefs, "IND_05_yearly_coefficients_data.csv")
if(nrow(indiv_yearly_coefs) > 0){
  p_indiv_time_varying <- ggplot(indiv_yearly_coefs, aes(x = year, y = estimate, color = controls, fill = controls)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), alpha = 0.2, linetype = 0) +
    geom_line(linewidth = 1) +
    facet_wrap(~metric, scales = "free_y", labeller = "label_both") +
    labs(title = "Time-Varying Effects on Individual Fertility Probability", subtitle = "Coefficients from Yearly Microdata Regressions (with MSA FE)", x = "Year", y = "Coefficient Estimate", color = "Model", fill = "Model")
  save_plot(p_indiv_time_varying, "IND_06_time_varying_coefficients.png", width = 12, height = 7)
} else {
  cat("Skipping individual time-varying plot as no valid coefficients were found.\n")
}


# ======================================================================================
# 4. TRACK 3: MIGRATION ANALYSIS
# ======================================================================================
cat("\n--- 4. MIGRATION ANALYSIS ---\n")
cat("  4.1 Generating migration plots...\n")
mig_data <- analysis_data[is_recent_mother == 1]
mig_agg_plot_data <- mig_data[, .(relocation_rate = weighted.mean(did_relocate, perwt, na.rm=TRUE)), by = .(met2013, year, msa_pti)]
p_mig_agg <- ggplot(mig_agg_plot_data[year == max(year)], aes(x = msa_pti, y = relocation_rate)) +
  geom_point(alpha=0.7, color="firebrick") + geom_smooth(method="lm", se=FALSE) +
  scale_y_continuous(labels=scales::percent) +
  labs(title="MSA Relocation Rate vs. Housing Cost for New Mothers", subtitle = paste("N =", nrow(mig_agg_plot_data[year==max(year)])), x="MSA Price-to-Income Ratio", y="Relocation Rate of New Mothers")
save_plot(p_mig_agg, "MIG_01_agg_relocation_vs_pti.png")
p_mig_indiv <- mig_data[, .(relocation_rate = weighted.mean(did_relocate, perwt, na.rm=TRUE)), by = .(educ_factor)] %>%
  ggplot(aes(x = fct_reorder(educ_factor, relocation_rate), y = relocation_rate, fill = educ_factor)) +
  geom_col(show.legend = FALSE) + scale_y_continuous(labels=scales::percent) + scale_fill_viridis_d(option="rocket") +
  labs(title="Relocation Rate of New Mothers by Education", subtitle = paste("N =", nrow(mig_data)), x="Educational Attainment", y="Relocation Rate (Moved from different PUMA)")
save_plot(p_mig_indiv, "MIG_02_indiv_relocation_by_educ.png")
cat("  4.2 Running migration models...\n")
mig_agg_data <- mig_data[, .(relocation_rate = weighted.mean(did_relocate, perwt, na.rm = TRUE), msa_pti = first(msa_pti), msa_rti = first(msa_rti)), by = .(met2013, year)]
m_mig_agg_pti <- feols(relocation_rate ~ msa_pti | year, data = mig_agg_data)
m_mig_agg_rti <- feols(relocation_rate ~ msa_rti | year, data = mig_agg_data)
mig_indiv_data <- mig_data[!is.na(msa_pti) & !is.na(msa_rti) & !is.na(educ_factor)]
m_mig_indiv_pti <- feglm(did_relocate ~ msa_pti + age + age_sq + educ_factor | year, data=mig_indiv_data, family="binomial", weights=~perwt)
m_mig_indiv_rti <- feglm(did_relocate ~ msa_rti + age + age_sq + educ_factor | year, data=mig_indiv_data, family="binomial", weights=~perwt)
save_table_html(list("Agg PTI"=m_mig_agg_pti, "Agg RTI"=m_mig_agg_rti, "Indiv PTI"=m_mig_indiv_pti, "Indiv RTI"=m_mig_indiv_rti), "MIG_03_migration_models.html")
save_table_latex(
    list("Agg PTI"=m_mig_agg_pti, "Agg RTI"=m_mig_agg_rti, "Indiv PTI"=m_mig_indiv_pti, "Indiv RTI"=m_mig_indiv_rti),
    filename = "mig_migration_models.tex",
    title = "Logit Models of Inter-MSA Relocation",
    label = "tab:migration_models"
)


# ======================================================================================
# 5. MIGRATION ORIGINS AND DESTINATIONS ANALYSIS
# ======================================================================================
cat("\n--- 5. MIGRATION ORIGINS AND DESTINATIONS ANALYSIS ---\n")

cat("  5.1 Preparing data with origin and destination housing costs...\n")
origin_msa_var <- "migmet131"
if (!origin_msa_var %in% names(analysis_data)) { stop(sprintf("FATAL ERROR: Migration MSA variable '%s' not found.", origin_msa_var)) }
analysis_data[, (origin_msa_var) := as.character(get(origin_msa_var))]
msa_metrics_lookup <- msa_panel_clean[, .(met2013, year, msa_pti_origin = msa_pti, msa_rti_origin = msa_rti)]
origin_analysis_data <- merge(analysis_data, msa_metrics_lookup, by.x = c(origin_msa_var, "year"), by.y = c("met2013", "year"), all.x = TRUE)
origin_analysis_data <- origin_analysis_data[age %between% c(18, 44) & get(origin_msa_var) != "0" & met2013 != "0" & !is.na(msa_pti_origin) & !is.na(msa_pti)]
cat(sprintf("  -> Prepared dataset with %s valid observations for origin analysis.\n", format(nrow(origin_analysis_data), big.mark = ",")))

cat("  5.2 Analyzing PUSH factor: Does a high-cost origin predict relocation?\n")
m_push_pti <- feglm(did_relocate ~ msa_pti_origin + age + age_sq + educ_factor | year, data = origin_analysis_data, family = "binomial", weights = ~perwt)
m_push_rti <- feglm(did_relocate ~ msa_rti_origin + age + age_sq + educ_factor | year, data = origin_analysis_data, family = "binomial", weights = ~perwt)
save_table_latex(
  list("Relocation ~ Origin PTI"=m_push_pti, "Relocation ~ Origin RTI"=m_push_rti),
  filename = "mig_origin_push_models.tex",
  title = "Migration PUSH Factor Analysis: Effect of Origin Housing Costs on Relocation",
  label = "tab:mig_origin_push"
)

cat("  5.3 Analyzing OUTCOME: Do people move to cheaper cities?\n")
inter_msa_movers <- origin_analysis_data[did_relocate == 1 & get(origin_msa_var) != met2013]
inter_msa_movers[, pti_change := msa_pti - msa_pti_origin]
inter_msa_movers[, rti_change := msa_rti - msa_rti_origin]
p_pti_change_dist <- ggplot(inter_msa_movers, aes(x = pti_change)) +
  geom_density(fill = "firebrick", alpha = 0.7) + geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(title = "Change in Price-to-Income Ratio for Movers", subtitle = "For inter-MSA movers only.", x = "PTI (Destination) - PTI (Origin)", y = "Density")
p_rti_change_dist <- ggplot(inter_msa_movers, aes(x = rti_change)) +
  geom_density(fill = "steelblue", alpha = 0.7) + geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(title = "Change in Rent-to-Income Ratio for Movers", subtitle = "For inter-MSA movers only.", x = "RTI (Destination) - RTI (Origin)", y = "Density")
combined_change_plot <- p_pti_change_dist + p_rti_change_dist
save_plot(combined_change_plot, "MIG_ORIGIN_02_cost_change_distribution.png", width = 12, height = 6)
m_dest_vs_origin_pti <- feols(msa_pti ~ msa_pti_origin | year, data = inter_msa_movers)
m_dest_vs_origin_rti <- feols(msa_rti ~ msa_rti_origin | year, data = inter_msa_movers)
save_table_latex(
  list("Destination PTI ~ Origin PTI"=m_dest_vs_origin_pti, "Destination RTI ~ Origin RTI"=m_dest_vs_origin_rti),
  filename = "mig_origin_destination_models.tex",
  title = "Migration OUTCOME Analysis: Housing Costs in Destination vs. Origin for Movers",
  label = "tab:mig_origin_destination"
)
cat("--- MIGRATION ORIGIN/DESTINATION ANALYSIS COMPLETE ---\n")

# ======================================================================================
# 6. COMPARATIVE MIGRATION ANALYSIS: New Mothers vs. General Population
# ======================================================================================
cat("\n--- 6. COMPARATIVE MIGRATION ANALYSIS ---\n")
run_comparative_migration_analysis <- function(data, group_name, file_prefix) {
  cat(sprintf("\n  -- Running Analysis for: %s --\n", group_name))
  cat("    -> Modeling PUSH factors (origin costs -> relocation)...\n")
  m_push_pti <- feglm(did_relocate ~ msa_pti_origin + age + age_sq + educ_factor | year, data = data, family = "binomial", weights = ~perwt)
  m_push_rti <- feglm(did_relocate ~ msa_rti_origin + age + age_sq + educ_factor | year, data = data, family = "binomial", weights = ~perwt)
  save_table_latex(
      list("PTI Push"=m_push_pti, "RTI Push"=m_push_rti),
      filename = paste0(file_prefix, "_push.tex"),
      title = paste("Push Factors for", group_name),
      label = paste0("tab:", file_prefix, "_push")
  )
  inter_msa_movers <- data[did_relocate == 1 & get(origin_msa_var) != met2013]
  if (nrow(inter_msa_movers) < 50) {
    cat(sprintf("    -> WARNING: Too few inter-MSA movers (%s) to analyze outcomes. Skipping.\n", nrow(inter_msa_movers)))
    return(list(push_pti = m_push_pti, push_rti = m_push_rti, outcome_pti = NULL, outcome_rti = NULL))
  }
  cat("    -> Analyzing OUTCOMES for inter-MSA movers...\n")
  p_pti_change <- ggplot(inter_msa_movers, aes(x = msa_pti - msa_pti_origin)) + geom_density(fill = "firebrick", alpha = 0.7) + geom_vline(xintercept = 0, linetype = "dashed") + labs(title = paste("Change in PTI for Movers:", group_name), x = "PTI (Destination) - PTI (Origin)")
  p_rti_change <- ggplot(inter_msa_movers, aes(x = msa_rti - msa_rti_origin)) + geom_density(fill = "steelblue", alpha = 0.7) + geom_vline(xintercept = 0, linetype = "dashed") + labs(title = paste("Change in RTI for Movers:", group_name), x = "RTI (Destination) - RTI (Origin)")
  save_plot(p_pti_change + p_rti_change, paste0(file_prefix, "_cost_change_dist.png"), width = 12, height = 6)
  m_outcome_pti <- feols(msa_pti ~ msa_pti_origin | year, data = inter_msa_movers)
  m_outcome_rti <- feols(msa_rti ~ msa_rti_origin | year, data = inter_msa_movers)
  save_table_latex(
      list("PTI Outcome"=m_outcome_pti, "RTI Outcome"=m_outcome_rti),
      filename = paste0(file_prefix, "_outcome.tex"),
      title = paste("Outcome Models for", group_name),
      label = paste0("tab:", file_prefix, "_outcome")
  )
  cat(sprintf("    -> Analysis for %s complete.\n", group_name))
  return(list(push_pti = m_push_pti, push_rti = m_push_rti, outcome_pti = m_outcome_pti, outcome_rti = m_outcome_rti))
}
cat("  6.1 Defining and filtering analysis groups...\n")
base_analysis_data <- origin_analysis_data[age %between% c(18, 50)]
general_women_data <- base_analysis_data
cat(sprintf("  -> General Women (18-50) group size: %s\n", format(nrow(general_women_data), big.mark=",")))
new_mothers_data <- base_analysis_data[is_recent_mother == 1]
cat(sprintf("  -> New Mothers (18-50) group size: %s\n", format(nrow(new_mothers_data), big.mark=",")))
cat("  6.2 Running analysis on both groups...\n")
results_general <- run_comparative_migration_analysis(data = general_women_data, group_name = "General Women (18-50)", file_prefix = "mig_compare_general")
results_mothers <- run_comparative_migration_analysis(data = new_mothers_data, group_name = "New Mothers (18-50)", file_prefix = "mig_compare_mothers")
cat("\n  6.3 Generating final comparative LaTeX table...\n")
comparative_models <- list(
  "Push PTI (General)" = results_general$push_pti, "Push PTI (Mothers)" = results_mothers$push_pti,
  "Push RTI (General)" = results_general$push_rti, "Push RTI (Mothers)" = results_mothers$push_rti,
  "Outcome PTI (General)" = results_general$outcome_pti, "Outcome PTI (Mothers)" = results_mothers$outcome_pti,
  "Outcome RTI (General)" = results_general$outcome_rti, "Outcome RTI (Mothers)" = results_mothers$outcome_rti
)
comparative_models <- Filter(Negate(is.null), comparative_models)
save_table_latex(
    comparative_models,
    filename = "mig_comparative_origin_models.tex",
    title = "Comparative Migration PUSH Factors: General Women vs. New Mothers",
    label = "tab:mig_comparative_origin"
)
cat("--- COMPARATIVE MIGRATION ANALYSIS COMPLETE ---\n")

# ======================================================================================
# 7. DEEP DIVE: MSA-Level Context vs. Individual-Level Burden
# ======================================================================================
cat("\n--- 7. DEEP DIVE: MSA-Level vs. Individual-Level Analysis ---\n")
cat("  7.1 Filtering data for valid individual PTI/RTI...\n")
pti_comparison_data <- base_analysis_data[!is.na(individual_pti)]
cat(sprintf("  -> Observations for PTI comparison: %s\n", format(nrow(pti_comparison_data), big.mark=",")))
rti_comparison_data <- base_analysis_data[!is.na(individual_rti)]
cat(sprintf("  -> Observations for RTI comparison: %s\n", format(nrow(rti_comparison_data), big.mark=",")))
cat("  7.2 Modeling for General Women (18-50)...\n")
m_gen_msa_pti <- feglm(did_relocate ~ msa_pti_origin + age + age_sq + educ_factor | year, data = pti_comparison_data, family = "binomial", weights = ~perwt)
m_gen_ind_pti <- feglm(did_relocate ~ individual_pti + age + age_sq + educ_factor | year, data = pti_comparison_data, family = "binomial", weights = ~perwt)
m_gen_msa_rti <- feglm(did_relocate ~ msa_rti_origin + age + age_sq + educ_factor | year, data = rti_comparison_data, family = "binomial", weights = ~perwt)
m_gen_ind_rti <- feglm(did_relocate ~ individual_rti + age + age_sq + educ_factor | year, data = rti_comparison_data, family = "binomial", weights = ~perwt)
cat("  7.3 Modeling for New Mothers (18-50)...\n")
m_mom_msa_pti <- feglm(did_relocate ~ msa_pti_origin + age + age_sq + educ_factor | year, data = pti_comparison_data[is_recent_mother==1], family = "binomial", weights = ~perwt)
m_mom_ind_pti <- feglm(did_relocate ~ individual_pti + age + age_sq + educ_factor | year, data = pti_comparison_data[is_recent_mother==1], family = "binomial", weights = ~perwt)
m_mom_msa_rti <- feglm(did_relocate ~ msa_rti_origin + age + age_sq + educ_factor | year, data = rti_comparison_data[is_recent_mother==1], family = "binomial", weights = ~perwt)
m_mom_ind_rti <- feglm(did_relocate ~ individual_rti + age + age_sq + educ_factor | year, data = rti_comparison_data[is_recent_mother==1], family = "binomial", weights = ~perwt)
cat("  7.4 Generating final comparative LaTeX tables...\n")
final_models_list <- list(
  "General: MSA PTI" = m_gen_msa_pti, "General: Indiv. PTI" = m_gen_ind_pti, "Mothers: MSA PTI" = m_mom_msa_pti, "Mothers: Indiv. PTI" = m_mom_ind_pti,
  "General: MSA RTI" = m_gen_msa_rti, "General: Indiv. RTI" = m_gen_ind_rti, "Mothers: MSA RTI" = m_mom_msa_rti, "Mothers: Indiv. RTI" = m_mom_ind_rti
)
save_table_latex(
    final_models_list,
    filename = "mig_final_msa_vs_individual_models.tex",
    title = "Comparison of Aggregate (MSA) vs. Individual-Level Cost Metrics in Predicting Relocation",
    label = "tab:mig_final_msa_vs_indiv"
)
cat("--- FINAL MIGRATION COMPARISON COMPLETE ---\n")

# ======================================================================================
# 8. TIME-VARYING ANALYSIS OF THE MIGRATION "PUSH" EFFECT
# ======================================================================================
cat("\n--- 8. TIME-VARYING MIGRATION PUSH ANALYSIS ---\n")
cat("  8.1 Defining time periods for analysis...\n")
base_analysis_data[, period := fcase(year %in% 2005:2007, "2005-2007 (Pre-GFC)", year %in% 2008:2012, "2008-2012 (GFC/Recovery)", year %in% 2013:2017, "2013-2017 (Mid-Cycle)", year >= 2018, "2018+ (Late-Cycle/Pandemic)")]
period_levels <- c("2005-2007 (Pre-GFC)", "2008-2012 (GFC/Recovery)", "2013-2017 (Mid-Cycle)", "2018+ (Late-Cycle/Pandemic)")
base_analysis_data[, period := factor(period, levels = period_levels)]

# Define the groups for repeated analyses
groups_to_run <- list(
  "General Women" = base_analysis_data,
  "New Mothers" = base_analysis_data[is_recent_mother == 1]
)
all_years <- sort(unique(base_analysis_data$year))

cat("  -> Time periods assigned to data.\n")
cat("  8.2 Running regressions for each time period and group...\n")
time_series_results <- list()
for (group_name in names(groups_to_run)) {
  for (p in levels(base_analysis_data$period)) {
    period_data <- groups_to_run[[group_name]][period == p]
    cat(sprintf("    -> Running model for %s, Period: %s (N=%s)\n", group_name, p, nrow(period_data)))
    model_result <- tryCatch({feglm(did_relocate ~ msa_rti_origin + age + age_sq + educ_factor | year, data = period_data, family = "binomial", weights = ~perwt)}, error = function(e) NULL)
    if (!is.null(model_result)) {
      tidy_result <- broom::tidy(model_result) %>% filter(term == "msa_rti_origin") %>% mutate(period = p, group = group_name)
      time_series_results <- append(time_series_results, list(tidy_result))
    }
  }
}
combined_coeffs_df <- rbindlist(time_series_results)
save_df_latex(
    df = combined_coeffs_df,
    filename = "mig_time_series_coefficients.tex",
    title = "Time-Varying PUSH Effect of Housing Costs on Relocation (Grouped by Era)",
    label = "tab:mig_time_series_coeffs"
)
cat("  8.3 Plotting the time-varying push effect...\n")
combined_coeffs_df[, `:=`(ci_lower = estimate - 1.96 * std.error, ci_upper = estimate + 1.96 * std.error)]
p_time_varying_push <- ggplot(combined_coeffs_df, aes(x = period, y = estimate, color = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(width = 0.3), width = 0.2, linewidth = 1) +
  scale_color_viridis_d(option = "plasma", end = 0.7, name = "Population Group") + theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "The Evolving 'Push' Effect of High Rent on Migration", subtitle = "Coefficient of 'MSA Origin RTI' on the probability of relocating, by period", x = "Time Period", y = "Coefficient Estimate (Log-Odds)")
save_plot(p_time_varying_push, "mig_time_series_push_effect_over_time.png", width = 10, height = 7)
cat("--- TIME-VARYING MIGRATION ANALYSIS COMPLETE ---\n")

# ======================================================================================
# 9. YEAR-BY-YEAR ANALYSIS OF THE MIGRATION "PUSH" EFFECT
# ======================================================================================
cat("\n--- 9. YEAR-BY-YEAR MIGRATION PUSH ANALYSIS ---\n")
cat("  9.1 Running regressions for each year and population group...\n")
yearly_results <- list()
for (group_name in names(groups_to_run)) {
  for (yr in all_years) {
    year_data <- groups_to_run[[group_name]][year == yr]
    if(nrow(year_data) < 100) next
    cat(sprintf("    -> Running model for %s, Year: %s (N=%s)\n", group_name, yr, nrow(year_data)))
    model_result <- tryCatch({feglm(did_relocate ~ msa_rti_origin + age + age_sq + educ_factor, data = year_data, family = "binomial", weights = ~perwt)}, error = function(e) NULL)
    if (!is.null(model_result)) {
      tidy_result <- broom::tidy(model_result) %>% filter(term == "msa_rti_origin") %>% mutate(year = yr, group = group_name, nobs = nobs(model_result))
      yearly_results <- append(yearly_results, list(tidy_result))
    }
  }
}
all_yearly_coeffs_df <- rbindlist(yearly_results)
save_df_latex(
    df = all_yearly_coeffs_df,
    filename = "mig_yearly_series_coefficients.tex",
    title = "Year-by-Year PUSH Effect of Housing Costs on Relocation",
    label = "tab:mig_yearly_series_coeffs"
)
cat("  9.2 Plotting the year-by-year push effect...\n")
all_yearly_coeffs_df[, `:=`(ci_lower = estimate - 1.96 * std.error, ci_upper = estimate + 1.96 * std.error)]
p_yearly_push <- ggplot(all_yearly_coeffs_df, aes(x = year, y = estimate, color = group, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 1) + geom_point(size = 2) + facet_wrap(~group, ncol = 1, scales = "free_y") +
  scale_color_viridis_d(option = "plasma", end = 0.7) + scale_fill_viridis_d(option = "plasma", end = 0.7) +
  theme(legend.position = "none") +
  labs(title = "Year-by-Year 'Push' Effect of High Rent on Migration", subtitle = "Coefficient of 'MSA Origin RTI' on the probability of relocating", x = "Year", y = "Coefficient Estimate (Log-Odds)")
save_plot(p_yearly_push, "mig_yearly_series_push_effect_by_year.png", width = 10, height = 8)
cat("--- YEAR-BY-YEAR MIGRATION ANALYSIS COMPLETE ---\n")

# ======================================================================================
# 10. YEAR-BY-YEAR ANALYSIS OF THE MIGRATION "PULL" EFFECT
# ======================================================================================
cat("\n--- 10. YEAR-BY-YEAR MIGRATION PULL ANALYSIS ---\n")
cat("  10.1 Preparing data of inter-MSA movers only...\n")
movers_only_data <- base_analysis_data[did_relocate == 1 & get(origin_msa_var) != met2013]
cat(sprintf("  -> Created dataset of %s inter-MSA movers for analysis.\n", format(nrow(movers_only_data), big.mark=",")))
cat("  10.2 Running 'pull' regressions for each year, group, and metric...\n")
yearly_pull_results <- list()
for (group_name in names(groups_to_run)) {
  for (yr in all_years) {
    for (metric in c("PTI", "RTI")) {
      year_data <- groups_to_run[[group_name]][year == yr]
      if(nrow(year_data) < 8) next
      cat(sprintf("    -> Running pull model for %s, %s, Year: %s (N=%s)\n", group_name, metric, yr, nrow(year_data)))
      outcome_var <- if (metric == "PTI") "msa_pti" else "msa_rti"; origin_var  <- if (metric == "PTI") "msa_pti_origin" else "msa_rti_origin"
      model_result <- tryCatch({feols(as.formula(paste(outcome_var, "~", origin_var)), data = year_data)}, error = function(e) NULL)
      if (!is.null(model_result)) {
        tidy_result <- broom::tidy(model_result) %>% filter(term == origin_var) %>% mutate(year = yr, group = group_name, metric = metric, nobs = nobs(model_result))
        yearly_pull_results <- append(yearly_pull_results, list(tidy_result))
      }
    }
  }
}
all_pull_coeffs_df <- rbindlist(yearly_pull_results)
save_df_latex(
    df = all_pull_coeffs_df,
    filename = "mig_yearly_pull_coefficients.tex",
    title = "Year-by-Year PULL Effect of Housing Costs on Movers' Destination Choice",
    label = "tab:mig_yearly_pull_coeffs"
)
cat("  10.3 Plotting the year-by-year pull effect...\n")
all_pull_coeffs_df[, `:=`(ci_lower = estimate - 1.96 * std.error, ci_upper = estimate + 1.96 * std.error)]
p_yearly_pull <- ggplot(all_pull_coeffs_df, aes(x = year, y = estimate, color = group, fill = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") + geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 1) + geom_point(size = 2) + facet_wrap(~metric, ncol = 1, scales = "free_y", labeller = label_both) +
  scale_color_viridis_d(option = "cividis") + scale_fill_viridis_d(option = "cividis") +
  labs(title = "Strength of the 'Pull' Toward Affordability for Movers, by Year", subtitle = "A value < 1 indicates a pull to cheaper cities.", x = "Year", y = "Coefficient Estimate")
save_plot(p_yearly_pull, "mig_yearly_pull_effect_by_year.png", width = 10, height = 8)
cat("--- YEAR-BY-YEAR MIGRATION PULL ANALYSIS COMPLETE ---\n")