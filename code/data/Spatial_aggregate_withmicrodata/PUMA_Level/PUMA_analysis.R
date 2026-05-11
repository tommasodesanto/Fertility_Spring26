# ======================================================================================
#
#    PUMA-Level Fertility Analysis
#    Author: Gemini AI Assistant (as programming partner)
#    Date: July 15, 2024
#
#    This script adapts the V3 master analysis to perform a full dual-track
#    analysis at the PUMA (Public Use Microdata Area) level instead of the MSA level.
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
  "sf", "viridis", "patchwork", "scales", "spatstat.geom", "broom"
)
install_and_load(required_packages)

# ---- Set modelsummary backend to kableExtra for simpler LaTeX ----
options(modelsummary_factory_latex = "kableExtra")

# ---- Set up directories and file paths ----
# Using the same 5% sample as the MSA analysis
data_file <- file.path("Spatial_aggregate_withmicrodata", "processed_data", "fertility_microdata_clean_5pct_stratified_sample.rds")
# New output directory for PUMA-level results
output_dir <- file.path("Spatial_aggregate_withmicrodata", "PUMA_Level", "analysis_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Clear the output directory for a clean run
cat("Clearing previous results from PUMA output directory...\n")
unlink(file.path(output_dir, "*"), recursive = TRUE)

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
              legend.position = "bottom",
              strip.text = element_text(face="bold", size=11)
            ))

cat("Setup complete. Outputs ->", output_dir, "\n")


# ======================================================================================
# 1. DATA PREPARATION (FOR PUMA-LEVEL ANALYSIS)
# ======================================================================================
cat("\n--- 1. DATA LOADING AND PREPARATION ---\n")

if (!file.exists(data_file)) { stop(paste("Data file not found:", data_file)) }
acs_data <- readRDS(data_file)
setDT(acs_data)

# ---- Feature Engineering (already done in script 03, just ensuring factors are correct) ----
cat("  1.1 Verifying base features...\n")
acs_data[, educ_factor := factor(educ_factor, levels = c("Less than High School", "High School Grad", "Some College", "Bachelors", "Graduate"))]
acs_data[, age_sq := age^2]
acs_data[, age_group := cut(age, breaks = c(14, 24, 34, 44, 51), labels = c("15-24", "25-34", "35-44", "45-50"), right = FALSE)]
acs_data[, did_relocate := as.integer(migrate1d %in% c(21, 22, 23, 24))]

# ---- 1.2 Individual-Level Housing Cost Metrics ----
cat("  1.2 Calculating individual-level housing cost metrics...\n")
acs_data[hhincome <= 0, hhincome := NA]
acs_data[, individual_pti := valueh / hhincome]
acs_data[valueh == 9999999, individual_pti := NA]
acs_data[individual_pti <= 0 | individual_pti > 50, individual_pti := NA]

acs_data[, individual_rti := (rent * 12) / hhincome]
acs_data[rent <= 0, individual_rti := NA]
acs_data[individual_rti <= 0 | individual_rti > 1, individual_rti := NA]

# ---- 1.3 Build Aggregate (PUMA-level) Panel ----
cat("  1.3 Building the aggregate PUMA-level panel...\n")
# CRITICAL CHANGE: Group by state and puma instead of met2013
# UPDATED: Also calculate population density directly from the 'density' variable.
puma_panel <- acs_data[hhincome > 0,
                      .(
                        # Fertility outcomes
                        mean_age_first_birth = weighted.mean(age_at_first_birth, perwt, na.rm = TRUE),
                        fertility_rate = weighted.mean(is_recent_mother, perwt, na.rm = TRUE),
                        pct_hh_with_children = weighted.mean(nchild > 0, hhwt, na.rm = TRUE),

                        # Housing Metrics: Median of individual ratios
                        puma_pti = tryCatch(weighted.median(individual_pti, hhwt, na.rm = TRUE), error = function(e) NA_real_),
                        puma_rti = tryCatch(weighted.median(individual_rti, hhwt, na.rm = TRUE), error = function(e) NA_real_),

                        # Other PUMA-level controls
                        puma_population = sum(perwt, na.rm = TRUE),
                        puma_median_hh_income = weighted.median(hhincome, hhwt, na.rm = TRUE),
                        # Use the existing 'density' variable, weighted by household weight
                        puma_pop_density = weighted.mean(density, hhwt, na.rm = TRUE),
                        pct_grad_plus = weighted.mean(educ_factor == "Graduate", perwt, na.rm = TRUE)
                      ), by = .(statefip, puma, year)]

# Create a unique PUMA identifier for merging
puma_panel[, puma_id := paste0(statefip, puma)]

# Clean and finalize panel
puma_panel_clean <- puma_panel[is.finite(puma_pti) & is.finite(puma_rti) & is.finite(mean_age_first_birth) & is.finite(puma_pop_density) & puma_pop_density > 0 & is.finite(puma_median_hh_income) & puma_median_hh_income > 0]
puma_panel_clean[, log_puma_pop_density := log(puma_pop_density)]
puma_panel_clean[, log_puma_median_hh_income := log(puma_median_hh_income)]

# ---- 1.4 Create Final Individual-Level Dataset ----
cat("  1.4 Merging PUMA aggregate metrics back into microdata...\n")
analysis_data <- merge(acs_data, puma_panel_clean, by = c("statefip", "puma", "year"), all.x = TRUE)
cat("Data preparation complete.\n")


# ======================================================================================
# 2. TRACK 1: AGGREGATE (PUMA-LEVEL) ANALYSIS
# ======================================================================================
cat("\n--- 2. TRACK 1: AGGREGATE (PUMA-LEVEL) ANALYSIS ---\n")
latest_year_data_agg <- puma_panel_clean[year == max(year)]

# ---- 2.1 Graphical EDA ----
cat("  2.1 Generating PUMA-level plots...\n")
p_agg_rti <- ggplot(latest_year_data_agg, aes(x = puma_rti, y = mean_age_first_birth)) +
  geom_point(alpha = 0.4, aes(color = pct_grad_plus)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "darkblue") +
  scale_color_viridis_c(name = "Pct with\nGrad Degree", option="magma") +
  labs(title = "PUMA Age at First Birth vs. Rent-to-Income Ratio", 
       subtitle = paste("Latest Year:", max(latest_year_data_agg$year), "| N =", nrow(latest_year_data_agg)),
       x = "PUMA Rent-to-Income Ratio", y = "Mean Age at First Birth")
save_plot(p_agg_rti, "PUMA_AGG_01_age_vs_rti.png")

# ---- 2.2 Aggregate Panel Regressions (FE Models) ----
cat("  2.2 Running PUMA aggregate panel regressions...\n")
agg_controls <- "log_puma_median_hh_income + pct_grad_plus + log_puma_pop_density"
agg_models <- list()
for(outcome in c("mean_age_first_birth", "fertility_rate")){
  for(metric in c("puma_pti", "puma_rti")){
    agg_models[[paste(outcome, metric, "simple", sep="_")]] <- feols(as.formula(paste(outcome, "~", metric)), data = puma_panel_clean, fixef = "year")
    agg_models[[paste0(outcome, "_controlled")]] <- feols(as.formula(paste(outcome, "~", metric, "+", agg_controls)), data = puma_panel_clean, fixef = c("statefip", "puma", "year"))
  }
}

# ---- 2.3 Time-Varying Effects (Aggregate Coefficients) ----
cat("  2.3 Analyzing time-varying effects at the PUMA level...\n")
agg_yearly_results <- list()
agg_outcomes <- c("fertility_rate", "mean_age_first_birth")
agg_metrics <- c("puma_pti", "puma_rti", "log_puma_pop_density")
min_obs_required <- 20

for (outcome in agg_outcomes) {
  for (metric in agg_metrics) {
    for (yr in sort(unique(puma_panel_clean$year))) {
        yr_data <- puma_panel_clean[year == yr]
        if (sum(!is.na(yr_data[[metric]])) < min_obs_required) {
          next
        }
        m_simple <- lm(as.formula(paste(outcome, "~", metric)), data = yr_data)
        n_obs <- nobs(m_simple)
        agg_yearly_results <- append(agg_yearly_results, list(broom::tidy(m_simple) %>% filter(term==metric) %>% mutate(year=yr, outcome=outcome, metric=metric, controls="Simple", nobs = n_obs)))
        
        controls_to_use <- if (metric == "log_puma_pop_density") {
          "log_puma_median_hh_income + pct_grad_plus"
        } else {
          "log_puma_pop_density + pct_grad_plus"
        }
        m_controlled <- lm(as.formula(paste(outcome, "~", metric, "+", controls_to_use)), data = yr_data)
        n_obs_ctrl <- nobs(m_controlled)
        agg_yearly_results <- append(agg_yearly_results, list(broom::tidy(m_controlled) %>% filter(term==metric) %>% mutate(year=yr, outcome=outcome, metric=metric, controls="With Controls", nobs = n_obs_ctrl)))
    }
  }
}
agg_yearly_coefs <- rbindlist(agg_yearly_results)
save_table(agg_yearly_coefs, "PUMA_AGG_03_yearly_coefficients_data.csv")

# ---- 2.4 Saving Individual Time-Varying Plots ----
cat("  2.4 Saving individual time-varying PUMA plots...\n")
individual_plots_dir <- file.path(output_dir, "PUMA_AGG_04_individual_plots")
dir.create(individual_plots_dir, showWarnings = FALSE, recursive = TRUE)

plot_combinations <- unique(agg_yearly_coefs[, .(outcome, metric)])

for (i in 1:nrow(plot_combinations)) {
    o <- plot_combinations$outcome[i]
    m <- plot_combinations$metric[i]
    plot_data <- agg_yearly_coefs[outcome == o & metric == m]

    if (nrow(plot_data) > 0) {
        # Create a summary data frame for N obs to prevent overplotting
        nobs_text_data <- plot_data[, .(nobs = first(nobs)), by = .(year)]

        p_individual <- ggplot(plot_data, aes(x = year, y = estimate, color = controls, fill = controls)) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
            geom_ribbon(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), alpha = 0.2, linetype = 0) +
            geom_line(linewidth = 1) +
            geom_text(data = nobs_text_data, aes(x = year, y = -Inf, label = paste0("N=", nobs)),
                      inherit.aes = FALSE, vjust = -0.5, size = 2.5, color = "gray40") +
            labs(
                title = paste("Time-Varying Effect of", m),
                subtitle = paste("Dependent Variable:", o),
                x = "Year", y = "Coefficient Estimate"
            ) +
            theme(plot.margin = margin(t = 5, r = 5, b = 25, l = 5)) # Add bottom margin

        plot_filename <- file.path(individual_plots_dir, sprintf("PUMA_plot_%s_on_%s.png", m, o))
        ggsave(plot_filename, plot = p_individual, width = 9, height = 6, dpi = 300, bg = "white")
        cat(sprintf("  -> Saved individual PUMA plot: %s\n", plot_filename))
    }
}

cat("\n--- PUMA analysis script finished successfully. ---\n")

# ======================================================================================
# 3. TRACK 2: INDIVIDUAL-LEVEL (MICRODATA) ANALYSIS
# ======================================================================================
cat("\n--- 3. TRACK 2: INDIVIDUAL-LEVEL PUMA ANALYSIS ---\n")
indiv_data_clean <- analysis_data[!is.na(puma_pti) & !is.na(puma_rti) & !is.na(educ_factor)]

# ---- 3.1 Main Individual Models (FE Models) ----
cat("  3.1 Running main individual-level regressions...\n")
indiv_controls <- "age + age_sq + educ_factor"
indiv_models <- list(
  "Fertility ~ PTI (Simple)" = feglm(is_recent_mother ~ puma_pti | puma_id + year, data=indiv_data_clean, family="binomial", weights=~perwt),
  "Fertility ~ PTI (Controlled)" = feglm(as.formula(paste("is_recent_mother ~ puma_pti +", indiv_controls)), data=indiv_data_clean, family="binomial", fixef=c("puma_id", "year"), weights=~perwt),
  "Fertility ~ RTI (Simple)" = feglm(is_recent_mother ~ puma_rti | puma_id + year, data=indiv_data_clean, family="binomial", weights=~perwt),
  "Fertility ~ RTI (Controlled)" = feglm(as.formula(paste("is_recent_mother ~ puma_rti +", indiv_controls)), data=indiv_data_clean, family="binomial", fixef=c("puma_id", "year"), weights=~perwt)
)

# ---- 3.2 Heterogeneity via Interactions ----
cat("  3.2 Analyzing heterogeneity with interaction models...\n")
int_model_educ_pti <- feglm(is_recent_mother ~ educ_factor * puma_pti + age + age_sq | statefip + puma + year, data = indiv_data_clean, family = "binomial", weights = ~perwt)
int_model_age_pti <- feglm(is_recent_mother ~ age_group * puma_pti + age + age_sq + educ_factor | statefip + puma + year, data = indiv_data_clean, family = "binomial", weights = ~perwt)

mfx_educ_plot <- plot_predictions(int_model_educ_pti, condition = c("puma_pti", "educ_factor")) + 
  labs(title = "Effect of Housing Price on Fertility by Education (PUMA)", 
       subtitle = paste("N =", nobs(int_model_educ_pti)),
       x = "PUMA Price-to-Income Ratio")
save_plot(mfx_educ_plot, "PUMA_IND_03_mfx_pti_by_education.png")

mfx_age_plot <- plot_predictions(int_model_age_pti, condition = c("puma_pti", "age_group")) + 
  labs(title = "Effect of Housing Price on Fertility by Age Group (PUMA)", 
       subtitle = paste("N =", nobs(int_model_age_pti)),
       x = "PUMA Price-to-Income Ratio")
save_plot(mfx_age_plot, "PUMA_IND_04_mfx_pti_by_age_group.png")

# ---- 3.3 Time-Varying Effects (Individual Coefficients) ----
cat("  3.3 Analyzing time-varying effects at the individual level...\n")
indiv_yearly_results <- list()
for (metric in c("puma_pti", "puma_rti")) {
    for (yr in sort(unique(indiv_data_clean$year))) {
        yr_data <- indiv_data_clean[year == yr]
        m_simple <- safely(feglm)(as.formula(paste("is_recent_mother ~", metric)), data=yr_data, family="binomial", fixef="puma_id", weights=~perwt)
        if(!is.null(m_simple$result)) {
          n_obs <- nobs(m_simple$result)
          indiv_yearly_results <- append(indiv_yearly_results, list(tidy(m_simple$result) %>% filter(grepl(metric, term)) %>% mutate(year=yr, metric=metric, controls="Simple", nobs = n_obs)))
        }
        
        m_controlled <- safely(feglm)(as.formula(paste("is_recent_mother ~", metric, "+", indiv_controls)), data=yr_data, family="binomial", fixef="puma_id", weights=~perwt)
        if(!is.null(m_controlled$result)) {
          n_obs <- nobs(m_controlled$result)
          indiv_yearly_results <- append(indiv_yearly_results, list(tidy(m_controlled$result) %>% filter(grepl(metric, term)) %>% mutate(year=yr, metric=metric, controls="With Controls", nobs = n_obs)))
        }
    }
}
indiv_yearly_coefs <- rbindlist(indiv_yearly_results)
save_table(indiv_yearly_coefs, "PUMA_IND_05_yearly_coefficients_data.csv")

if(nrow(indiv_yearly_coefs) > 0){
    nobs_text_data_indiv <- indiv_yearly_coefs[, .(nobs = first(nobs)), by = .(year, metric)]

    p_indiv_time_varying <- ggplot(indiv_yearly_coefs, aes(x = year, y = estimate, color = controls, fill = controls)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_ribbon(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), alpha = 0.2, linetype = 0) +
      geom_line(linewidth = 1) +
      geom_text(data = nobs_text_data_indiv, aes(x = year, y = -Inf, label = paste0("N=", nobs)),
                inherit.aes = FALSE, vjust = -0.5, size = 2.5, color = "gray40") +
      facet_wrap(~metric, scales = "free_y", labeller = "label_both") +
      labs(title = "Time-Varying Effects on Individual Fertility Probability (PUMA)", subtitle = "Coefficients from Yearly Microdata Regressions (with PUMA FE)",
           x = "Year", y = "Coefficient Estimate", color = "Model", fill = "Model") +
      theme(plot.margin = margin(t = 5, r = 5, b = 25, l = 5)) # Add bottom margin
    save_plot(p_indiv_time_varying, "PUMA_IND_06_time_varying_coefficients.png", width = 12, height = 7)
} else {
    cat("Skipping PUMA individual time-varying plot as no valid coefficients were found.\n")
}

cat("\n--- FULL PUMA analysis script finished successfully. ---\n") 