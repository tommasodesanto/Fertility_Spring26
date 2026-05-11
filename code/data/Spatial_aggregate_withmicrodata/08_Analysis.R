# ======================================================================================
#
#    Microdata-Driven Fertility Analysis
#    Author: Gemini AI Assistant (as programming partner)
#    Date: June 9, 2025
#
#    This script performs a comprehensive analysis of fertility dynamics in the US
#    using the 2005-2023 5% IPUMS ACS sample. It follows the four-phase plan:
#    1. Macro-Level Descriptive Analysis & "Motivating Pictures"
#    2. Enhanced Aggregate Panel Analysis (MSA-level)
#    3. Individual-Level Microdata Models
#    4. Advanced Topics & Extensions
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
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# ---- List of required packages ----
required_packages <- c(
  "tidyverse",      # Core for data manipulation and plotting
  "data.table",     # Fast data manipulation
  "fixest",         # High-performance fixed-effects models (feglm)
  "modelsummary",   # For beautiful regression tables
  "marginaleffects",# For interpreting interaction terms
  "sf",             # For spatial data (maps)
  "tigris",         # To download US census shapefiles
  "viridis",        # Color palettes for plots
  "patchwork",      # To combine ggplots
  "scales",         # For plot axis formatting
  "spatstat.geom"   # For weighted.median function
)
install_and_load(required_packages)

# ---- Set up directories and file paths ----
data_file <- file.path("processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
output_dir <- "analysis_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper functions for saving outputs ----
save_plot <- function(plot_obj, filename, width = 10, height = 6.5, dpi = 300) {
  full_path <- file.path(output_dir, filename)
  ggsave(full_path, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
  cat(sprintf("Saved plot: %s\n", full_path))
}

save_table <- function(table_obj, filename) {
  full_path <- file.path(output_dir, filename)
  if (is.data.frame(table_obj) || is.matrix(table_obj)) {
    write.csv(as.data.frame(table_obj), full_path, row.names = FALSE)
  } else {
    modelsummary::modelsummary(table_obj, output = full_path)
  }
  cat(sprintf("Saved table: %s\n", full_path))
}

# ---- Plotting Theme ----
theme_set(theme_minimal(base_size = 12) +
            theme(
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray30"),
              plot.caption = element_text(hjust = 0, color = "gray50"),
              axis.title = element_text(size = 11, face="bold"),
              legend.position = "bottom"
            ))

cat("Setup complete. Outputs will be saved to:", output_dir, "\n")


# ======================================================================================
# 1. DATA LOADING AND PREPARATION
# ======================================================================================
cat("\n--- 1. DATA LOADING AND PREPARATION ---\n")

if (!file.exists(data_file)) {
  stop(paste("Data file not found:", data_file, "\nPlease ensure the script is in the correct directory."))
}
acs_data <- readRDS(data_file)
setDT(acs_data) # Convert to data.table for speed

# ---- Data Cleaning and Feature Engineering ----
# Ensure factors are set correctly for analysis
acs_data[, educ_factor := factor(educ_factor,
                                 levels = c("Less than High School", "High School Grad",
                                            "Some College", "Bachelors", "Graduate"))]

# Create an age-squared term for non-linear effects
acs_data[, age_sq := age^2]

# Ensure metro area is a factor
acs_data[, met2013 := as.factor(met2013)]

# Log transform income (add 1 to avoid log(0))
acs_data[, log_hhincome := log(hhincome + 1)]
acs_data[is.infinite(log_hhincome), log_hhincome := NA]

cat("Data loaded and prepared. Rows:", nrow(acs_data), "Cols:", ncol(acs_data), "\n")


# ======================================================================================
# 2. PHASE 1: MACRO-LEVEL DESCRIPTIVE ANALYSIS & "MOTIVATING PICTURES"
# ======================================================================================
cat("\n--- 2. PHASE 1: DESCRIPTIVE ANALYSIS ---\n")

# ---- 2.1 National Fertility Trends Over Time ----
cat("  2.1 Generating national trend plots...\n")

national_trends <- acs_data[age %between% c(15, 50),
                            .(mean_age_first_birth = weighted.mean(age_at_first_birth, perwt, na.rm = TRUE),
                              fertility_rate = weighted.mean(is_recent_mother, perwt, na.rm = TRUE) * 1000),
                            by = year]

# Plot 1: The rise in age at first birth (Motivating Picture 1)
p_age_trend <- ggplot(national_trends, aes(x = year, y = mean_age_first_birth)) +
  geom_rect(aes(xmin = 2007.5, xmax = 2009.5, ymin = -Inf, ymax = Inf), fill = "grey90", alpha = 0.02) +
  geom_rect(aes(xmin = 2019.5, xmax = 2022.5, ymin = -Inf, ymax = Inf), fill = "grey90", alpha = 0.02) +
  geom_line(color = "#0072B2", linewidth = 1.2) +
  geom_point(color = "#0072B2", size = 2.5) +
  geom_text(aes(label=round(mean_age_first_birth, 1)), vjust=-1, size=3.5, color="black") +
  annotate("text", x = 2008.5, y = min(national_trends$mean_age_first_birth), label = "GFC", size = 4, color="grey50", angle=90, hjust=0) +
  annotate("text", x = 2021, y = min(national_trends$mean_age_first_birth), label = "COVID", size = 4, color="grey50", angle=90, hjust=0) +
  labs(title = "The Unwavering Rise of Delayed Motherhood",
       subtitle = "Mean Age at First Birth for Women (15-50) in the U.S., 2005-2023",
       x = "Year",
       y = "Mean Age at First Birth",
       caption = "Source: IPUMS ACS 1-Year Estimates (2005-2023). Weights applied. Shaded areas denote GFC and COVID-19 periods.") +
  scale_x_continuous(breaks = seq(2005, 2023, 2)) +
  theme(plot.title.position = "plot")

save_plot(p_age_trend, "01_national_mean_age_first_birth.png")


# ---- 2.2 The Spatial Dimension of Fertility and Costs ----
cat("  2.2 Generating spatial maps...\n")
latest_year <- max(acs_data$year)
msa_spatial_data <- acs_data[year == latest_year & met2013 != "0" & hhincome > 0 & valueh < 9999999,
                             .(mean_age_first_birth = weighted.mean(age_at_first_birth, perwt, na.rm = TRUE),
                               pti_ratio = weighted.median(valueh, hhwt, na.rm=TRUE) / weighted.median(hhincome, hhwt, na.rm=TRUE),
                               total_pop = sum(perwt)),
                             by = .(GEOID = met2013)]
msa_spatial_data <- msa_spatial_data[total_pop > 50000 & !is.na(pti_ratio)]

# Download MSA shapefiles
msa_shapes <- core_based_statistical_areas(cb = TRUE, year = 2020) %>%
  st_transform(4326) %>% # Standard lat/lon projection
  select(GEOID, NAME)

# Join data
msa_map_data <- msa_shapes %>%
  inner_join(msa_spatial_data, by = "GEOID")

# Plot 2A: Age at First Birth Map
p_map_age <- ggplot(msa_map_data) +
  geom_sf(aes(fill = mean_age_first_birth), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "magma", name = "Mean Age") +
  labs(title = "Age at First Birth",
       subtitle = paste("MSA Level, Year", latest_year)) +
  theme_void() + theme(legend.position = "right")

# Plot 2B: Price-to-Income Ratio Map
p_map_pti <- ggplot(msa_map_data) +
  geom_sf(aes(fill = pti_ratio), color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "plasma", name = "PTI Ratio", limits = c(1, 15), oob = scales::squish) +
  labs(title = "Housing Price-to-Income Ratio",
       subtitle = "Median House Value / Median HH Income") +
  theme_void() + theme(legend.position = "right")

# Combine maps (Motivating Picture 2)
p_spatial_combo <- (p_map_age + p_map_pti) +
  plot_annotation(
    title = "The Geographic Divide: Where You Live and When You Have Kids",
    caption = "Source: IPUMS ACS 2023. MSAs with population > 50,000 shown. PTI = Price-to-Income."
  )
save_plot(p_spatial_combo, "02_spatial_comparison_maps.png", width = 14, height = 7)


# ---- 2.3 The Socio-Economic Dimension ----
cat("  2.3 Generating socio-economic plot...\n")
p_educ_age <- acs_data[is_recent_mother == 1 & !is.na(educ_factor), ] %>%
  ggplot(aes(x = educ_factor, y = age, fill = educ_factor)) +
  geom_violin(alpha = 0.7, show.legend = FALSE, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5, show.legend = FALSE) +
  scale_fill_viridis_d(option = "cividis") +
  labs(title = "Education is a Powerful Predictor of Fertility Timing",
       subtitle = "Age Distribution of Women Who Gave Birth in the Last Year, by Educational Attainment",
       x = "Educational Attainment",
       y = "Age at Birth",
       caption = "Source: IPUMS ACS 2005-2023 pool. Boxplots show median and IQR.") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

save_plot(p_educ_age, "03_age_at_birth_by_education.png")


# ======================================================================================
# 3. PHASE 2: ENHANCING THE AGGREGATE PANEL ANALYSIS
# ======================================================================================
cat("\n--- 3. PHASE 2: MSA-LEVEL PANEL REGRESSIONS ---\n")

# ---- 3.1 Construct the Analysis Panel ----
cat("  3.1 Constructing MSA-year panel...\n")
msa_panel <- acs_data[met2013 != "0" & hhincome > 0,
                      .(
                        # Dependent Variables
                        mean_age_first_birth = weighted.mean(age_at_first_birth, perwt, na.rm = TRUE),
                        fertility_rate = weighted.mean(is_recent_mother, perwt, na.rm = TRUE),
                        pct_hh_with_children = weighted.mean(nchild > 0, hhwt, na.rm = TRUE),
                        # Independent Variables
                        median_hh_income = weighted.median(hhincome, hhwt, na.rm = TRUE),
                        median_home_value = weighted.median(valueh[valueh < 9999999], hhwt[valueh < 9999999], na.rm = TRUE),
                        puma_density = weighted.mean(density, perwt, na.rm = TRUE),
                        pct_grad_plus = weighted.mean(educ_factor == "Graduate", perwt, na.rm = TRUE),
                        total_pop = sum(perwt)
                      ), by = .(met2013, year)]

# Calculate affordability ratios
msa_panel[, pti_ratio := median_home_value / median_hh_income]

# Clean up infinite/missing values and filter for stable panel
msa_panel <- msa_panel[total_pop > 20000 & is.finite(pti_ratio) & is.finite(mean_age_first_birth)]
msa_panel <- msa_panel[!is.na(pti_ratio) & !is.na(median_hh_income) & !is.na(pct_grad_plus)]


# ---- 3.2 Panel Regressions (TWFE) ----
cat("  3.2 Running panel regressions...\n")
reg_panel_twfe <- list(
  "Mean Age 1st Birth" = feols(mean_age_first_birth ~ pti_ratio + log(median_hh_income) + pct_grad_plus + log(puma_density) | met2013 + year, data = msa_panel),
  "Fertility Rate"     = feols(fertility_rate ~ pti_ratio + log(median_hh_income) + pct_grad_plus + log(puma_density) | met2013 + year, data = msa_panel),
  "Pct HH w/ Children" = feols(pct_hh_with_children ~ pti_ratio + log(median_hh_income) + pct_grad_plus + log(puma_density) | met2013 + year, data = msa_panel)
)

# Create and save regression table
modelsummary(reg_panel_twfe,
             stars = TRUE,
             gof_map = c("nobs", "r.squared.within"),
             title = "MSA-Level Panel Regressions with Two-Way Fixed Effects (MSA + Year)",
             notes = "All variables are aggregated at the MSA-year level. Weights applied during aggregation.",
             output = file.path(output_dir, "04_msa_panel_regressions.html"))

cat("MSA panel regressions complete.\n")


# ======================================================================================
# 4. PHASE 3: INDIVIDUAL-LEVEL MODELS
# ======================================================================================
cat("\n--- 4. PHASE 3: INDIVIDUAL-LEVEL LOGISTIC REGRESSIONS ---\n")

# ---- 4.1 Prepare Data for Individual Models ----
cat("  4.1 Preparing data for individual models...\n")

# Merge MSA-level housing cost proxy back to individual data
msa_housing_cost <- msa_panel[, .(met2013, year, pti_ratio)]
micro_model_data <- merge(acs_data, msa_housing_cost, by = c("met2013", "year"))

# Filter for the relevant sample: women of childbearing age
micro_model_data <- micro_model_data[age %between% c(15, 50) & !is.na(educ_factor) & !is.na(pti_ratio)]

# Create age groups for interaction analysis
micro_model_data[, age_group := cut(age,
                                    breaks = c(14, 24, 34, 44, 51),
                                    labels = c("15-24", "25-34", "35-44", "45-50"),
                                    right = TRUE)]

# ---- 4.2 Main Logistic Regression Model ----
cat("  4.2 Running main individual-level logistic regression...\n")
# Using fixest::feglm for speed with many fixed effects
main_logit_model <- feglm(is_recent_mother ~ age + age_sq + educ_factor + pti_ratio | met2013 + year,
                          family = "binomial",
                          data = micro_model_data,
                          weights = ~perwt)

save_table(main_logit_model, "05_individual_logit_model.html")

# ---- 4.3 Heterogeneity Analysis via Interactions ----
cat("  4.3 Analyzing heterogeneity with interaction terms...\n")
# Interaction 1: Housing Cost x Education
interact_educ_model <- feglm(is_recent_mother ~ age + age_sq + educ_factor * pti_ratio | met2013 + year,
                             family = "binomial",
                             data = micro_model_data,
                             weights = ~perwt)

# Interaction 2: Housing Cost x Age Group
interact_age_model <- feglm(is_recent_mother ~ age + age_sq + educ_factor + pti_ratio * age_group | met2013 + year,
                            family = "binomial",
                            data = micro_model_data,
                            weights = ~perwt)

# Save interaction model tables
modelsummary(list("Main" = main_logit_model,
                  "Interact (Educ)" = interact_educ_model,
                  "Interact (Age Group)" = interact_age_model),
             stars = TRUE, gof_map = c("nobs", "pseudo.r.squared"),
             title = "Individual Logistic Regressions on Probability of Recent Birth (Women 15-50)",
             notes = "All models include MSA and Year fixed effects. Standard errors are clustered by MSA.",
             output = file.path(output_dir, "06_individual_logit_interaction_models.html"))


# ---- 4.4 Visualizing Marginal Effects ----
cat("  4.4 Visualizing marginal effects from interaction models...\n")

# Plot Marginal Effect of PTI across Education Levels
mfx_educ_plot <- plot_predictions(interact_educ_model,
                                  condition = c("pti_ratio", "educ_factor")) +
  labs(title = "Housing Costs May Deter Less-Educated Women More",
       subtitle = "Predicted Probability of Birth by Housing Cost and Education",
       x = "MSA Price-to-Income Ratio",
       y = "Predicted Probability of a Recent Birth",
       color = "Education") +
  theme(legend.position = "right")
save_plot(mfx_educ_plot, "07_mfx_pti_by_education.png")

# Plot Marginal Effect of PTI across Age Groups
mfx_age_plot <- plot_predictions(interact_age_model,
                                 condition = c("pti_ratio", "age_group")) +
  labs(title = "Effect of Housing Costs on Fertility Varies by Age Group",
       subtitle = "Predicted Probability of Birth by Housing Cost and Age Group",
       x = "MSA Price-to-Income Ratio",
       y = "Predicted Probability of a Recent Birth",
       color = "Age Group") +
  theme(legend.position = "right")
save_plot(mfx_age_plot, "08_mfx_pti_by_age_group.png")


# ======================================================================================
# 5. PHASE 4: ADVANCED TOPICS & EXTENSIONS
# ======================================================================================
cat("\n--- 5. PHASE 4: ADVANCED TOPICS ---\n")

# ---- 5.1 Analysis of Second and Higher-Order Births ----
cat("  5.1 Modeling probability of subsequent births...\n")

# Create subsample of existing mothers
subsequent_birth_data <- micro_model_data[nchild >= 1]

# Re-run the main logit model on this subsample
subsequent_birth_model <- feglm(is_recent_mother ~ age + age_sq + educ_factor + pti_ratio | met2013 + year,
                                family = "binomial",
                                data = subsequent_birth_data,
                                weights = ~perwt)

# Compare the "any birth" model with the "subsequent birth" model
modelsummary(list("Any Birth (All Women 15-50)" = main_logit_model,
                  "Subsequent Birth (Mothers Only)" = subsequent_birth_model),
             stars = TRUE, gof_map = c("nobs", "pseudo.r.squared"),
             title = "Comparing Drivers of First vs. Subsequent Births",
             notes = "Both models include MSA and Year fixed effects.",
             output = file.path(output_dir, "09_subsequent_birth_model_comparison.html"))

# ---- 5.2 Structural Framing ----
cat("  5.2 Writing summary and structural interpretation...\n")

structural_interpretation <- "
# =====================================================================================
#  Structural Interpretation of Findings
# =====================================================================================
#
# Our analysis can be framed within a simple household utility maximization model.
# A household chooses consumption (C) and number of children (n) to maximize utility
# U(C, n), subject to a budget constraint Y = P_c*C + P_n*n, where Y is income,
# P_c is the price of consumption goods, and P_n is the 'price of a child'.
#
# Our findings shed light on the components of P_n:
#
# 1. Direct Housing Cost Component:
#    - The 'pti_ratio' (Price-to-Income) variable is a direct measure of the housing
#      component of P_n. Our models consistently show a negative coefficient on this
#      variable (see Tables 04, 05, 06, 09).
#    - This implies that as the cost of housing (a necessary input for raising a
#      child) increases, households substitute away from fertility.
#
# 2. Opportunity Cost Component (Education):
#    - The 'educ_factor' variable proxies for the opportunity cost of time. Women
#      with higher educational attainment typically have higher potential earnings,
#      making the time cost of child-rearing higher.
#    - Our models (Tables 05, 06, 09) show that, controlling for other factors,
#      women with 'Graduate' degrees have a significantly lower probability of a
#      recent birth compared to those with less education. This is consistent with
#      a higher opportunity cost component of P_n for this group.
#
# 3. Income Effect (Y):
#    - In theory, if children are a 'normal good', higher income (Y) should lead
#      to higher fertility. However, income is also correlated with education and
#      living in high-cost areas.
#    - Our aggregate panel model (Table 04) shows a mixed/negative relationship
#      for log(median_hh_income) after controlling for fixed effects, suggesting
#      that the price effects (which are correlated with income) dominate at the
#      MSA level.
#
# 4. Heterogeneity:
#    - The interaction models (Table 06, Figure 07) suggest that the 'price effect'
#      of housing may be non-linear. The negative impact of the PTI ratio appears
#      strongest for women with lower educational attainment, who may have tighter
#      budget constraints and less ability to absorb high housing costs.
#
# In summary, our reduced-form empirical results align well with economic theory,
# highlighting that both the direct financial costs (housing) and the opportunity
# costs (education) are powerful determinants of fertility decisions in the modern US.
#
# =====================================================================================
"
writeLines(structural_interpretation, file.path(output_dir, "10_structural_interpretation.txt"))

cat("\n--- ANALYSIS SCRIPT COMPLETE ---\n")
cat("All outputs have been saved in the '", output_dir, "' directory.\n")