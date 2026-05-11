# =============================================================================
# MSA PANEL DATA ANALYSIS (Based on ACS/Decennial Data from 2005)
# Version 2.1 - Overhauled with Baseline Affordability & Pct HH Children
#               EXTENDED with New Plots & Output Formats (v2 - More Coef Plots - INTEGRAL)
# =============================================================================

# =============================================================================
# 0. SETUP AND CONFIGURATION
# =============================================================================
cat("Setting up environment...\n")

# Function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat(paste("Installing", pkg, "...\n"))
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# List of required packages
required_packages <- c(
  "tidyverse",   # Data manipulation and ggplot2
  "plm",         # Panel data models
  "lfe",         # High-dimensional fixed effects (optional)
  "modelsummary",# For regression tables
  "broom",       # For tidying model outputs
  "corrplot",    # Correlation plots
  "viridis",     # Color scales
  "scales",      # Formatting for plots
  "ggrepel",     # Label positioning in ggplot
  "patchwork",   # Combine ggplots
  "lubridate",   # Date handling
  "naniar",      # For visualizing missing data
  "purrr"        # For map functions
)
install_and_load(required_packages)

# BASE PROJECT DIRECTORY (USER: Set this to your main R project folder)
BASE_PROJECT_DIR <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/R" # Please update this path

# Set paths
data_file_path <- file.path(BASE_PROJECT_DIR, "msa_panel_output", "msa_panel_complete.csv")

# Create output subdirectory for analysis results
output_dir_analysis <- file.path(BASE_PROJECT_DIR, "msa_panel_analysis_output_v2.1_extended_v2") # New folder for this version
dir.create(output_dir_analysis, showWarnings = FALSE, recursive = TRUE)

# Set theme for plots
theme_set(theme_minimal(base_size = 11) +
            theme(
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 12, hjust = 0.5),
              axis.title = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 9),
              plot.caption = element_text(size = 8, hjust = 1, color = "gray40"),
              plot.background = element_rect(fill = "white", colour = "white"),
              panel.background = element_rect(fill = "white", colour = "grey90"),
              panel.grid.major = element_line(colour = "grey92"),
              panel.grid.minor = element_line(colour = "grey95", linetype = "dashed"),
              strip.background = element_rect(fill = "grey85", colour = "grey90"),
              strip.text = element_text(face = "bold")
            ))

# Helper function for saving plots
save_plot <- function(plot_obj, filename, width = 10, height = 6, dpi = 300) {
  full_path <- file.path(output_dir_analysis, filename)
  ggsave(
    filename = full_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  cat(sprintf("Saved plot: %s\n", full_path))
}

# Helper function for saving tables
save_table <- function(table_obj, filename) {
  full_path <- file.path(output_dir_analysis, filename)
  if (endsWith(tolower(filename), ".txt")) {
    sink(full_path)
    tryCatch({
      print(table_obj)
    }, finally = {
      sink()
    })
  } else if (is.data.frame(table_obj) || is.matrix(table_obj)) {
    write_csv(as.data.frame(table_obj), full_path)
  } else {
    capture.output(print(table_obj), file = full_path)
  }
  cat(sprintf("Saved table/output: %s\n", full_path))
}


cat("Setup complete. Output will be saved to:", output_dir_analysis, "\n")

# =============================================================================
# 1. DATA LOADING AND INITIAL INSPECTION
# =============================================================================
cat("\n--- 1. DATA LOADING AND INITIAL INSPECTION ---\n")

if (!file.exists(data_file_path)) {
  stop(paste("Data file not found:", data_file_path, 
             "\nPlease ensure BASE_PROJECT_DIR and data_file_path are correctly set."))
}
msa_data_raw <- read_csv(data_file_path, show_col_types = FALSE)

cat("Raw data loaded successfully.\n")
print(glimpse(msa_data_raw))

cat("\nDistribution of data_source:\n"); print(table(msa_data_raw$data_source, useNA = "ifany"))
cat("\nDistribution of fertility_data_source:\n"); print(table(msa_data_raw$fertility_data_source, useNA = "ifany"))
cat("\nYears with S1301 birth rate data:\n"); print(table(msa_data_raw$year[!is.na(msa_data_raw$birth_rate_s1301)], useNA = "ifany"))

if (nrow(msa_data_raw) > 0 && ncol(msa_data_raw) > 0) {
  p_miss_vis <- vis_miss(msa_data_raw, cluster = TRUE, sort_miss = TRUE) +
    labs(title = "Missing Data Pattern in Loaded Dataset",
         subtitle = paste0(nrow(msa_data_raw), " observations, ", ncol(msa_data_raw), " variables")) +
    theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_miss_vis, "missing_data_visualization.png", width = 12, height = 7)
}

# =============================================================================
# 2. DATA CLEANING AND PREPARATION
# =============================================================================
cat("\n--- 2. DATA CLEANING AND PREPARATION ---\n")

msa_panel <- msa_data_raw %>%
  mutate(
    year = as.numeric(year)
  ) %>%
  group_by(GEOID) %>% # Group by GEOID first
  mutate(
    NAME_latest = last(NAME[order(year)]), # Use latest MSA name for consistency
  ) %>%
  ungroup() %>% # Ungroup before further global mutations
  mutate(
    NAME = NAME_latest, # Overwrite original NAME
    log_total_population = ifelse(total_population > 0, log(total_population), NA_real_),
    log_median_hh_income = ifelse(median_hh_income > 0, log(median_hh_income), NA_real_),
    log_median_home_value = ifelse(median_home_value > 0, log(median_home_value), NA_real_),
    log_median_gross_rent = ifelse(median_gross_rent > 0, log(median_gross_rent), NA_real_),
    log_price_to_income = ifelse(price_to_income > 0, log(price_to_income), NA_real_),
    log_annual_rent_to_income = ifelse(annual_rent_to_income > 0, log(annual_rent_to_income), NA_real_),
    
    msa_size_category = factor(case_when(
      total_population >= 5000000 ~ "Very Large (5M+)",
      total_population >= 1000000 ~ "Large (1M-5M)",
      total_population >= 500000  ~ "Medium (500K-1M)",
      total_population >= 250000  ~ "Small (250K-500K)",
      TRUE                        ~ "Very Small (<250K)"
    ), levels = c("Very Small (<250K)", "Small (250K-500K)", "Medium (500K-1M)", "Large (1M-5M)", "Very Large (5M+)")),
    
    period = factor(case_when(
      year <= 2007 ~ "Pre-GFC (<=2007)",
      year <= 2012 ~ "GFC & Early Recovery (2008-2012)",
      year <= 2019 ~ "Mid-Cycle (2013-2019)",
      TRUE         ~ "COVID & Post (>=2020)"
    ), levels = c("Pre-GFC (<=2007)", "GFC & Early Recovery (2008-2012)", "Mid-Cycle (2013-2019)", "COVID & Post (>=2020)")),
    
    housing_affordability_quintile_dynamic = ifelse(!is.na(price_to_income), ntile(price_to_income, 5), NA_integer_), # Dynamic per year
    housing_affordability_category_dynamic = factor(case_when(
      housing_affordability_quintile_dynamic == 1 ~ "Most Affordable (Q1 Dynamic)",
      housing_affordability_quintile_dynamic == 2 ~ "Affordable (Q2 Dynamic)",
      housing_affordability_quintile_dynamic == 3 ~ "Moderate (Q3 Dynamic)",
      housing_affordability_quintile_dynamic == 4 ~ "Expensive (Q4 Dynamic)",
      housing_affordability_quintile_dynamic == 5 ~ "Least Affordable (Q5 Dynamic)",
      TRUE ~ NA_character_
    ), levels = c("Most Affordable (Q1 Dynamic)", "Affordable (Q2 Dynamic)", "Moderate (Q3 Dynamic)", "Expensive (Q4 Dynamic)", "Least Affordable (Q5 Dynamic)"))
  ) %>%
  group_by(GEOID) %>%
  mutate(
    price_to_income_baseline_val = first(price_to_income[!is.na(price_to_income)][order(year)]),
  ) %>%
  ungroup() %>%
  mutate(
    baseline_afford_quintile = ifelse(!is.na(price_to_income_baseline_val), ntile(price_to_income_baseline_val, 5), NA_integer_), # Baseline across all MSAs' first PTI
    baseline_afford_category = factor(case_when(
      baseline_afford_quintile == 1 ~ "Most Affordable (Q1 Baseline)",
      baseline_afford_quintile == 2 ~ "Affordable (Q2 Baseline)",
      baseline_afford_quintile == 3 ~ "Moderate (Q3 Baseline)",
      baseline_afford_quintile == 4 ~ "Expensive (Q4 Baseline)",
      baseline_afford_quintile == 5 ~ "Least Affordable (Q5 Baseline)",
      TRUE ~ NA_character_
    ), levels = c("Most Affordable (Q1 Baseline)", "Affordable (Q2 Baseline)", "Moderate (Q3 Baseline)", "Expensive (Q4 Baseline)", "Least Affordable (Q5 Baseline)"))
  ) %>%
  select(
    GEOID, NAME, year, data_source, collection_date, period,
    total_population, log_total_population, msa_size_category,
    children_under_18, pct_under_18, children_under_5,
    women_15_49, women_15_50, child_woman_ratio,
    women_gave_birth, birth_rate_per_1000, fertility_data_source, birth_rate_s1301,
    total_households, households_with_children, pct_households_with_children,
    median_hh_income, log_median_hh_income,
    median_home_value, log_median_home_value,
    median_gross_rent, log_median_gross_rent,
    price_to_income, log_price_to_income, housing_affordability_category_dynamic,
    price_to_income_baseline_val, baseline_afford_category,
    annual_rent_to_income, log_annual_rent_to_income,
    pct_college_plus,
    -NAME_latest 
  ) %>%
  arrange(GEOID, year)

cat("Data cleaning and preparation complete.\n")

if (n_distinct(msa_panel$GEOID) > 1 && n_distinct(msa_panel$year) > 1) {
  years_in_panel <- sort(unique(msa_panel$year))
  balanced_panel <- msa_panel %>%
    group_by(GEOID) %>%
    filter(n_distinct(year) == length(years_in_panel) && all(years_in_panel %in% year)) %>%
    ungroup()
  cat(sprintf("Full panel: %d observations, %d unique MSAs, %d unique years.\n",
              nrow(msa_panel), n_distinct(msa_panel$GEOID), length(years_in_panel)))
  cat(sprintf("Balanced panel: %d observations, %d unique MSAs.\n",
              nrow(balanced_panel), n_distinct(balanced_panel$GEOID)))
  if (nrow(balanced_panel) > 0) save_table(balanced_panel, "msa_panel_balanced.csv")
} else {
  balanced_panel <- msa_panel
  cat("Not enough diversity in GEOIDs or years for a meaningful balanced panel from this dataset.\n")
}

sanity_vars <- c("log_total_population", "log_median_hh_income", "log_price_to_income", "log_annual_rent_to_income", "birth_rate_per_1000", "pct_households_with_children")
for (v in sanity_vars) {
  if (v %in% names(msa_panel) && sum(!is.na(msa_panel[[v]])) > 0) {
    p_sanity <- ggplot(msa_panel, aes_string(x = v)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black", alpha = 0.7, na.rm = TRUE) +
      geom_density(color = "red", na.rm = TRUE) +
      labs(title = paste("Distribution of", v), x = v, y = "Density")
    save_plot(p_sanity, paste0("dist_", v, ".png"), width = 7, height = 5)
  }
}

# =============================================================================
# 3. EXPLORATORY DATA ANALYSIS (EDA)
# =============================================================================
cat("\n--- 3. EXPLORATORY DATA ANALYSIS ---\n")

numeric_cols_for_summary <- msa_panel %>% select(where(is.numeric), -year, -contains("quintile"), -price_to_income_baseline_val)
summary_stats_table <- numeric_cols_for_summary %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    N_Obs = sum(!is.na(value)),
    Mean = round(mean(value, na.rm = TRUE), 2), Median = round(median(value, na.rm = TRUE), 2),
    SD = round(sd(value, na.rm = TRUE), 2), Min = round(min(value, na.rm = TRUE), 2),
    Max = round(max(value, na.rm = TRUE), 2), Pct_Missing = round(mean(is.na(value)) * 100, 1),
    .groups = "drop"
  )
print(summary_stats_table)
save_table(summary_stats_table, "summary_statistics_all_vars.csv")

if (nrow(msa_panel) > 0 && "total_population" %in% names(msa_panel)) {
  national_trends <- msa_panel %>%
    group_by(year) %>%
    summarise(
      across(c(birth_rate_per_1000, child_woman_ratio, pct_households_with_children,
               median_hh_income, median_home_value, price_to_income, annual_rent_to_income,
               pct_college_plus, pct_under_18),
             ~weighted.mean(., w = total_population, na.rm = TRUE), .names = "avg_{.col}"),
      n_msas = n() ) %>% ungroup() %>% arrange(year)
  save_table(national_trends, "national_trends_summary.csv")
  
  p_fert_trends <- ggplot(national_trends, aes(x = year)) +
    geom_line(aes(y = avg_birth_rate_per_1000, color = "Birth Rate (B13002)"), linewidth = 1) +
    geom_point(aes(y = avg_birth_rate_per_1000, color = "Birth Rate (B13002)"), size = 2) +
    geom_line(aes(y = avg_child_woman_ratio / 5, color = "Child-Woman Ratio (scaled /5)"), linewidth = 1) + 
    geom_point(aes(y = avg_child_woman_ratio / 5, color = "Child-Woman Ratio (scaled /5)"), size = 2) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    labs(title = "National Avg Fertility Trends (Pop. Weighted)", x = "Year", y = "Rate / Scaled Rate", color = "Indicator") +
    theme(legend.position = "bottom")
  save_plot(p_fert_trends, "national_fertility_trends.png")
  
  p_hh_child_trends <- ggplot(national_trends, aes(x = year)) +
    geom_line(aes(y = avg_pct_households_with_children, color = "Pct Households w/ Children"), linewidth = 1) +
    geom_point(aes(y = avg_pct_households_with_children, color = "Pct Households w/ Children"), size = 2) +
    geom_line(aes(y = avg_pct_under_18, color = "Pct Population < 18"), linewidth = 1, linetype="dashed") +
    geom_point(aes(y = avg_pct_under_18, color = "Pct Population < 18"), size = 2) +
    scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 0.1)) +
    scale_color_viridis_d(option = "magma", end=0.8) +
    labs(title = "National Avg Child Presence Indicators (Pop. Weighted)", x = "Year", y = "Percentage", color = "Indicator") +
    theme(legend.position = "bottom")
  save_plot(p_hh_child_trends, "national_child_presence_trends.png")
  
  p_econ_trends <- ggplot(national_trends, aes(x = year)) +
    geom_line(aes(y = avg_median_hh_income, color = "Median HH Income"), linewidth = 1) +
    geom_point(aes(y = avg_median_hh_income, color = "Median HH Income"), size = 2) +
    scale_y_continuous(labels = scales::dollar_format()) + scale_color_viridis_d(option = "viridis") +
    labs(title = "National Avg Economic Trends (Pop. Weighted)", subtitle = "Nominal values", x = "Year", y = "USD", color = "Indicator") +
    theme(legend.position = "bottom")
  save_plot(p_econ_trends, "national_econ_trends.png")
  
  p_housing_ratio_trends <- ggplot(national_trends, aes(x = year)) +
    geom_line(aes(y = avg_price_to_income, color = "Price-to-Income Ratio"), linewidth = 1) +
    geom_point(aes(y = avg_price_to_income, color = "Price-to-Income Ratio"), size = 2) +
    geom_line(aes(y = avg_annual_rent_to_income * 10, color = "Annual Rent-to-Income Ratio (x10)"), linewidth = 1) + 
    geom_point(aes(y = avg_annual_rent_to_income * 10, color = "Annual Rent-to-Income Ratio (x10)"), size = 2) +
    scale_color_viridis_d(option = "cividis") +
    labs(title = "National Avg Housing Affordability Ratios (Pop. Weighted)", x = "Year", y = "Ratio / Scaled Ratio", color = "Indicator") +
    theme(legend.position = "bottom")
  save_plot(p_housing_ratio_trends, "national_housing_ratio_trends.png")
}

latest_year_val <- if (nrow(msa_panel) > 0) max(msa_panel$year, na.rm = TRUE) else NA
if (!is.na(latest_year_val)) {
  latest_year_data <- msa_panel %>% filter(year == latest_year_val)
  
  if (nrow(latest_year_data) > 0) {
    p_br_msasize <- ggplot(latest_year_data, aes(x = msa_size_category, y = birth_rate_per_1000, fill = msa_size_category)) +
      geom_boxplot(na.rm = TRUE, alpha=0.7, outlier.shape = NA) + geom_jitter(width=0.1, alpha=0.3, na.rm=TRUE, size=0.8) +
      scale_fill_viridis_d(guide = "none") +
      labs(title = paste("Birth Rates by MSA Size, Year", latest_year_val), x = "MSA Size Category", y = "Birth Rate per 1000") +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
    save_plot(p_br_msasize, paste0("latest_year_",latest_year_val,"_br_by_msasize.png"))
    
    p_inc_fert <- ggplot(latest_year_data, aes(x = median_hh_income, y = birth_rate_per_1000)) +
      geom_point(aes(size = total_population, color = msa_size_category), alpha = 0.6, na.rm = TRUE) +
      geom_smooth(method = "lm", se = FALSE, color = "darkred", na.rm = TRUE, formula = y ~ x) +
      scale_x_continuous(labels = scales::dollar_format()) + scale_size_continuous(name = "Population", labels = scales::comma_format()) +
      scale_color_viridis_d(name = "MSA Size") +
      labs(title = paste("Income vs. Birth Rate, Year", latest_year_val), x = "Median Household Income", y = "Birth Rate per 1000")
    save_plot(p_inc_fert, paste0("latest_year_",latest_year_val,"_income_vs_fertility.png"))
    
    p_pti_fert <- ggplot(latest_year_data, aes(x = price_to_income, y = birth_rate_per_1000)) +
      geom_point(aes(size = total_population, color = msa_size_category), alpha = 0.6, na.rm = TRUE) +
      geom_smooth(method = "lm", se = FALSE, color = "darkblue", na.rm = TRUE, formula = y ~ x) +
      scale_size_continuous(name = "Population", labels = scales::comma_format()) + scale_color_viridis_d(name = "MSA Size") +
      labs(title = paste("Price-to-Income Ratio vs. Birth Rate, Year", latest_year_val), x = "Median Home Value / Median HH Income", y = "Birth Rate per 1000")
    save_plot(p_pti_fert, paste0("latest_year_",latest_year_val,"_pti_vs_fertility.png"))
    
    if (nrow(latest_year_data) >= 10) {
      extreme_fert_msas <- latest_year_data %>% filter(!is.na(birth_rate_per_1000)) %>% arrange(birth_rate_per_1000) %>%
        slice(c(1:5, (n()-4):n())) %>% mutate(Rank_Category = rep(c("Lowest 5 Birth Rates", "Highest 5 Birth Rates"), each = 5)) %>%
        select(Rank_Category, NAME, year, birth_rate_per_1000, child_woman_ratio, median_hh_income, price_to_income, total_population)
      save_table(extreme_fert_msas, paste0("latest_year_",latest_year_val,"_extreme_fertility_msas.csv"))
    }
  }
} else {
  cat("Could not determine latest year for cross-sectional EDA.\n")
  latest_year_data <- data.frame() 
}


if (nrow(msa_panel) > 0 && all(c("birth_rate_per_1000", "child_woman_ratio") %in% names(msa_panel))) {
  p_comp_fert_metrics <- ggplot(msa_panel, aes(x = child_woman_ratio, y = birth_rate_per_1000)) +
    geom_point(aes(color = factor(year), size = total_population), alpha = 0.5, na.rm = TRUE) +
    geom_smooth(method="lm", color="black", na.rm=TRUE, formula = y ~ x) +
    scale_color_viridis_d(name = "Year", option = "turbo") + scale_size_continuous(name = "Population", labels = scales::comma_format()) +
    labs(title = "Comparison: CWR vs. Birth Rate", x = "Child-Woman Ratio", y = "Birth Rate per 1000") +
    theme(legend.position = "right")
  save_plot(p_comp_fert_metrics, "fertility_metrics_comparison.png", width = 11, height = 7)
  cor_fert_metrics <- cor(msa_panel$child_woman_ratio, msa_panel$birth_rate_per_1000, use = "pairwise.complete.obs")
  cat(sprintf("\nCorrelation between CWR and Birth Rate: %.3f\n", cor_fert_metrics))
}

# =============================================================================
# 4. CORRELATION ANALYSIS
# =============================================================================
cat("\n--- 4. CORRELATION ANALYSIS ---\n")
cor_vars <- c("birth_rate_per_1000", "child_woman_ratio", "pct_households_with_children", "pct_under_18",
              "median_hh_income", "median_home_value", "price_to_income", "annual_rent_to_income",
              "pct_college_plus", "log_total_population")
cor_vars_exist <- cor_vars[cor_vars %in% names(msa_panel)]

if (length(cor_vars_exist) > 1) {
  cor_data <- msa_panel %>% select(all_of(cor_vars_exist))
  cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")
  if (nrow(cor_matrix) > 1 && ncol(cor_matrix) > 1) {
    png(file.path(output_dir_analysis, "correlation_heatmap.png"), width=10, height=8, units="in", res=300, bg="white")
    corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", addCoef.col = "black",
             tl.col = "black", tl.srt = 45, diag = FALSE,
             col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(100),
             title = "Correlation Matrix of Key Variables", mar = c(0,0,1,0))
    dev.off(); cat("Saved correlation heatmap.\n")
    save_table(as.data.frame(cor_matrix), "correlation_matrix.csv")
  }
}

# =============================================================================
# 5. PANEL REGRESSION ANALYSIS
# =============================================================================
cat("\n--- 5. PANEL REGRESSION ANALYSIS ---\n")
ivs <- c("log_median_hh_income", "price_to_income", "annual_rent_to_income", "pct_college_plus", "log_total_population")
dvs <- c("birth_rate_per_1000", "child_woman_ratio", "pct_households_with_children")
lfe_available <- "lfe" %in% .packages()

missing_ivs <- setdiff(ivs, names(msa_panel))
if (length(missing_ivs) > 0) {
  stop(paste("Error: The following IVs are not found in msa_panel:", paste(missing_ivs, collapse=", ")))
}
missing_dvs <- setdiff(dvs, names(msa_panel))
if (length(missing_dvs) > 0) {
  stop(paste("Error: The following DVs are not found in msa_panel:", paste(missing_dvs, collapse=", ")))
}

reg_data_base <- msa_panel %>%
  mutate(across(all_of(ivs), ~ scale(.)[,1], .names = "std_{.col}")) %>%
  select(GEOID, NAME, year, all_of(dvs), starts_with("std_"))

models_collated <- list()
for (dv_name in dvs) {
  cat(paste0("\n--- Regressions for DV: ", dv_name, " ---\n"))
  current_ivs_std <- paste0("std_", ivs)
  
  missing_std_ivs <- setdiff(current_ivs_std, names(reg_data_base))
  if (length(missing_std_ivs) > 0) {
    cat(paste0("Skipping ", dv_name, ": Standardized IVs missing: ", paste(missing_std_ivs, collapse=", "), "\n"))
    next
  }
  
  model_data <- reg_data_base %>% select(GEOID, NAME, year, all_of(dv_name), all_of(current_ivs_std)) %>% drop_na()
  
  if (nrow(model_data) < (length(current_ivs_std) + 20) || n_distinct(model_data$GEOID) < 5) {
    cat(paste0("Skipping ", dv_name, ": Insufficient data (", nrow(model_data), " obs, ", n_distinct(model_data$GEOID), " MSAs) after NA removal for model estimation.\n")); next
  }
  
  p_model_data <- pdata.frame(model_data, index = c("GEOID", "year"))
  formula_str <- paste(dv_name, "~", paste(current_ivs_std, collapse = " + "))
  
  models_collated[[paste0(dv_name, "_OLS")]] <- lm(as.formula(formula_str), data = p_model_data)
  models_collated[[paste0(dv_name, "_YFE")]] <- lm(as.formula(paste(formula_str, "+ factor(year)")), data = p_model_data)
  
  fe_suffix <- if(lfe_available) "_lfe" else "_plm"
  if (lfe_available) {
    models_collated[[paste0(dv_name, "_MSA_FE", fe_suffix)]] <- tryCatch(felm(as.formula(paste(formula_str, "| GEOID")), data = p_model_data), error=function(e) {warning(e); NULL})
    models_collated[[paste0(dv_name, "_TWFE", fe_suffix)]] <- tryCatch(felm(as.formula(paste(formula_str, "| GEOID + year")), data = p_model_data), error=function(e) {warning(e); NULL})
  } else {
    models_collated[[paste0(dv_name, "_MSA_FE", fe_suffix)]] <- tryCatch(plm(as.formula(formula_str), data = p_model_data, model = "within", effect = "individual"), error=function(e) {warning(e); NULL})
    models_collated[[paste0(dv_name, "_TWFE", fe_suffix)]] <- tryCatch(plm(as.formula(formula_str), data = p_model_data, model = "within", effect = "twoways"), error=function(e) {warning(e); NULL})
  }
  models_collated <- Filter(Negate(is.null), models_collated) 
}

if (length(models_collated) > 0) {
  reg_table_main_obj <- modelsummary(
    models_collated, fmt = "%.3f", estimate = "{estimate}{stars}", statistic = "({std.error})",
    coef_map = paste0("std_", ivs), gof_map = c("nobs", "r.squared"),
    title = "Panel Regression: Fertility & Household Composition Metrics",
    notes = list("Std. errors in parentheses.", "IVs are standardized.")
  )
  save_table(reg_table_main_obj, "panel_regressions_main.txt")
  cat("Main panel regression table (console format) saved to panel_regressions_main.txt\n")
  print(reg_table_main_obj)
  
  birth_rate_models <- models_collated[grep("^birth_rate_per_1000_", names(models_collated))]
  if (length(birth_rate_models) > 0) {
    key_ivs_std_plot <- paste0("std_", ivs) 
    plot_data_list <- list()
    for(i in seq_along(birth_rate_models)){
      model_name <- names(birth_rate_models)[i]
      model_obj <- birth_rate_models[[i]]
      if(!is.null(model_obj)){
        tidy_df <- broom::tidy(model_obj, conf.int = TRUE) %>%
          filter(term %in% key_ivs_std_plot) %>%
          mutate(model = model_name)
        plot_data_list[[model_name]] <- tidy_df
      }
    }
    
    if(length(plot_data_list) > 0){
      plot_data_coefs <- bind_rows(plot_data_list) %>%
        mutate(term = factor(term, levels = rev(key_ivs_std_plot)), 
               model_type = case_when(
                 grepl("_OLS$", model) ~ "OLS",
                 grepl("_YFE$", model) ~ "Year FE",
                 grepl("_MSA_FE", model) ~ "MSA FE",
                 grepl("_TWFE", model) ~ "Two-Way FE",
                 TRUE ~ model 
               ),
               model_type = factor(model_type, levels = c("OLS", "Year FE", "MSA FE", "Two-Way FE")))
      
      plot_data_coefs$term_clean <- gsub("std_log_", "Std. Log ", plot_data_coefs$term)
      plot_data_coefs$term_clean <- gsub("std_", "Std. ", plot_data_coefs$term_clean)
      plot_data_coefs$term_clean <- gsub("_", " ", plot_data_coefs$term_clean)
      plot_data_coefs$term_clean <- tools::toTitleCase(plot_data_coefs$term_clean)
      
      term_levels_cleaned <- sapply(rev(key_ivs_std_plot), function(t) {
        clean_t <- gsub("std_log_", "Std. Log ", t)
        clean_t <- gsub("std_", "Std. ", clean_t)
        clean_t <- gsub("_", " ", clean_t)
        tools::toTitleCase(clean_t)
      })
      plot_data_coefs$term_clean <- factor(plot_data_coefs$term_clean, levels = term_levels_cleaned)
      
      p_birth_rate_coeffs <- ggplot(plot_data_coefs, aes(x = estimate, y = term_clean, color = model_type)) +
        geom_point(position = position_dodge(width = 0.5), size = 2.5) +
        geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                       height = 0, position = position_dodge(width = 0.5), linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        scale_color_viridis_d(name = "Model Specification") +
        labs(title = "Coefficient Plot for Birth Rate Regressions",
             subtitle = "Estimates and 95% Confidence Intervals for Standardized IVs",
             x = "Coefficient Estimate", y = "Independent Variable") +
        theme(legend.position = "bottom")
      save_plot(p_birth_rate_coeffs, "birth_rate_panel_model_coefficients.png", width = 12, height = 8)
    } else {
      cat("No valid model results to plot for birth rate coefficients.\n")
    }
  } else {
    cat("No birth rate models found in models_collated for coefficient plot.\n")
  }
} else {
  cat("No panel regression models were successfully estimated.\n")
}

# =============================================================================
# 6. HETEROGENEITY ANALYSIS (Example: Income Effect by MSA Size)
# =============================================================================
cat("\n--- 6. HETEROGENEITY ANALYSIS ---\n")
for (dv_hetero in c("birth_rate_per_1000", "pct_households_with_children")) {
  cat(paste0("\n--- Heterogeneity for DV: ", dv_hetero, " ---\n"))
  
  if (!dv_hetero %in% names(reg_data_base)) {
    cat(paste0("Skipping heterogeneity for ", dv_hetero, ": DV not found in reg_data_base.\n"))
    next
  }
  
  hetero_data <- reg_data_base %>%
    left_join(msa_panel %>% select(GEOID, year, msa_size_category) %>% distinct(), by = c("GEOID", "year")) %>%
    filter(!is.na(msa_size_category)) %>%
    drop_na(all_of(dv_hetero), std_log_median_hh_income, std_price_to_income) 
  
  required_std_ivs_hetero <- c("std_log_median_hh_income", "std_price_to_income", "std_annual_rent_to_income", "std_pct_college_plus", "std_log_total_population")
  if(!all(required_std_ivs_hetero %in% names(hetero_data))){
    cat(paste0("Skipping heterogeneity for ", dv_hetero, ": Not all required standardized IVs present in hetero_data.\n"))
    next
  }
  
  if (nrow(hetero_data) > 0 && n_distinct(hetero_data$msa_size_category) > 1) {
    formula_hetero_str <- paste(dv_hetero, "~ std_log_median_hh_income + std_price_to_income + std_annual_rent_to_income + std_pct_college_plus + std_log_total_population + factor(year)")
    
    models_by_size <- hetero_data %>% group_by(msa_size_category) %>%
      filter(n() > (length(required_std_ivs_hetero) + n_distinct(year) + 10) & n_distinct(GEOID) > 1) %>%
      do(model = tryCatch(lm(as.formula(formula_hetero_str), data = .), error = function(e) {warning(e); NULL})) %>%
      ungroup() %>% filter(!is.null(model))
    
    if (nrow(models_by_size) > 0) {
      coef_effects <- models_by_size %>% rowwise() %>% mutate(tidy_model = list(broom::tidy(model, conf.int = TRUE))) %>%
        unnest(tidy_model) %>% filter(term == "std_log_median_hh_income") %>%
        select(msa_size_category, estimate, conf.low, conf.high, p.value)
      save_table(coef_effects, paste0("hetero_income_effects_by_msasize_on_", dv_hetero, ".csv"))
      
      p_hetero <- ggplot(coef_effects, aes(x = msa_size_category, y = estimate, fill = msa_size_category)) +
        geom_col(alpha=0.8) + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
        geom_hline(yintercept = 0, linetype = "dashed") + scale_fill_viridis_d(guide="none") +
        labs(title = paste("Heterogeneous Effect of Std. Log Income on", dv_hetero),
             subtitle = "By MSA Size (Controls: Other IVs & Year FE)", x = "MSA Size", y = "Coefficient for Std. Log Income") +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
      save_plot(p_hetero, paste0("hetero_income_effect_on_", dv_hetero, ".png"), width=10, height=7)
    } else { cat("Not enough data for heterogeneity models for", dv_hetero, " by MSA size.\n") }
  } else { cat("Skipping heterogeneity for", dv_hetero, ": insufficient data or MSA size categories.\n") }
}

# =============================================================================
# 7. TIME SERIES PATTERNS FOR SELECT MSAS
# =============================================================================
cat("\n--- 7. TIME SERIES PATTERNS FOR SELECT MSAS ---\n")
if (nrow(msa_panel) > 0 && n_distinct(msa_panel$GEOID) > 0) {
  set.seed(42) 
  selected_geoids_ts <- sample(unique(msa_panel$GEOID), min(n_distinct(msa_panel$GEOID), 4))
  ts_data_sel <- msa_panel %>% filter(GEOID %in% selected_geoids_ts) %>% mutate(msa_label = str_extract(NAME, "^[^,]+"))
  
  if (nrow(ts_data_sel) > 0 && "birth_rate_per_1000" %in% names(ts_data_sel)) {
    p_ts_br <- ggplot(ts_data_sel, aes(x = year, y = birth_rate_per_1000, color = msa_label, group = msa_label)) +
      geom_line(linewidth = 1, na.rm = TRUE) + geom_point(size = 2, na.rm = TRUE) + scale_color_viridis_d(name = "MSA") +
      labs(title = "Birth Rate Trends in Selected MSAs", x = "Year", y = "Birth Rate per 1000") + theme(legend.position = "bottom")
    save_plot(p_ts_br, "birth_rate_trends_selected_msas.png", width = 10, height = 7)
  } else {
    cat("Not enough data or birth_rate_per_1000 missing for selected MSA time series plot.\n")
  }
}

# =============================================================================
# 8. HOUSING AFFORDABILITY ANALYSIS (DYNAMIC AND BASELINE)
# =============================================================================
cat("\n--- 8. HOUSING AFFORDABILITY IMPACTS ---\n")

if (nrow(msa_panel) > 0 && "housing_affordability_category_dynamic" %in% names(msa_panel) && sum(!is.na(msa_panel$housing_affordability_category_dynamic)) > 0) {
  aff_sum_dyn <- msa_panel %>% filter(!is.na(housing_affordability_category_dynamic) & !is.na(birth_rate_per_1000) & !is.na(total_population)) %>%
    group_by(year, housing_affordability_category_dynamic) %>%
    summarise(avg_birth_rate = weighted.mean(birth_rate_per_1000, total_population, na.rm = TRUE), n=n(), .groups = 'drop') %>% filter(n > 0)
  if(nrow(aff_sum_dyn) > 0) {
    p_br_aff_dyn <- ggplot(aff_sum_dyn, aes(x=year, y=avg_birth_rate, color=housing_affordability_category_dynamic, group=housing_affordability_category_dynamic)) +
      geom_line(linewidth=1.2, na.rm=TRUE) + geom_point(size=2, na.rm=TRUE) + scale_color_viridis_d(name="Affordability\n(Dynamic PTI Quintile)", option="mako") +
      labs(title="Birth Rates by DYNAMIC Housing Affordability", subtitle="Avg. birth rates within annually-defined PTI quintiles (Pop. Weighted)", x="Year", y="Avg Birth Rate") + theme(legend.position="bottom")
    save_plot(p_br_aff_dyn, "birth_rate_by_dynamic_housing_affordability.png", width=12, height=7)
  }
}

if (nrow(msa_panel) > 0 && "baseline_afford_category" %in% names(msa_panel) && sum(!is.na(msa_panel$baseline_afford_category)) > 0) {
  aff_sum_base <- msa_panel %>% filter(!is.na(baseline_afford_category) & !is.na(birth_rate_per_1000) & !is.na(total_population)) %>%
    group_by(year, baseline_afford_category) %>%
    summarise(avg_birth_rate = weighted.mean(birth_rate_per_1000, total_population, na.rm = TRUE), n=n(), .groups = 'drop') %>% filter(n > 0)
  if(nrow(aff_sum_base) > 0) {
    p_br_aff_base <- ggplot(aff_sum_base, aes(x=year, y=avg_birth_rate, color=baseline_afford_category, group=baseline_afford_category)) +
      geom_line(linewidth=1.2, na.rm=TRUE) + geom_point(size=2, na.rm=TRUE) + scale_color_viridis_d(name="Affordability\n(BASELINE PTI Quintile)", option="rocket") +
      labs(title="Birth Rates by BASELINE Housing Affordability", subtitle="Avg. birth rates for fixed MSA cohorts by baseline PTI (Pop. Weighted)", x="Year", y="Avg Birth Rate") + theme(legend.position="bottom")
    save_plot(p_br_aff_base, "birth_rate_by_baseline_housing_affordability.png", width=12, height=7)
  }
  
  reg_data_aff <- msa_panel %>% filter(!is.na(baseline_afford_category) & !is.na(birth_rate_per_1000)) %>%
    mutate(std_log_median_hh_income = scale(log_median_hh_income)[,1], 
           std_pct_college_plus = scale(pct_college_plus)[,1], 
           std_log_total_population = scale(log_total_population)[,1],
           baseline_afford_category = relevel(baseline_afford_category, ref = "Most Affordable (Q1 Baseline)")) %>%
    select(GEOID, year, birth_rate_per_1000, baseline_afford_category, std_log_median_hh_income, std_pct_college_plus, std_log_total_population) %>% drop_na()
  
  if (nrow(reg_data_aff) > 20 && n_distinct(reg_data_aff$GEOID) > 1 && n_distinct(reg_data_aff$baseline_afford_category) > 1) {
    p_reg_data_aff <- pdata.frame(reg_data_aff, index = c("GEOID", "year"))
    formula_aff_reg <- birth_rate_per_1000 ~ baseline_afford_category + std_log_median_hh_income + std_pct_college_plus + std_log_total_population
    
    aff_models <- list()
    aff_models[["Year FE"]] <- lm(update(formula_aff_reg, . ~ . + factor(year)), data = p_reg_data_aff)
    if (lfe_available) { aff_models[["TWFE_lfe"]] <- tryCatch(felm(update(formula_aff_reg, . ~ . | GEOID + year), data=p_reg_data_aff), error=function(e) {warning(e);NULL}) }
    else { aff_models[["TWFE_plm"]] <- tryCatch(plm(formula_aff_reg, data=p_reg_data_aff, model="within", effect="twoways"), error=function(e) {warning(e);NULL}) }
    aff_models <- Filter(Negate(is.null), aff_models)
    
    if(length(aff_models) > 0) {
      aff_reg_tbl_obj <- modelsummary(aff_models, fmt="%.3f", estimate="{estimate}{stars}", statistic="({std.error})",
                                      coef_rename=function(x) gsub("baseline_afford_category", "Afford Cat: ", x), gof_map=c("nobs", "r.squared"),
                                      title="Regression: Birth Rate vs. Baseline Housing Affordability",
                                      notes=list("Ref for Affordability: Most Affordable (Q1 Baseline).", "Controls: Std. Log Income, Pct College, Log Pop."))
      save_table(aff_reg_tbl_obj, "regression_baseline_affordability.txt")
      cat("Baseline affordability regression table (console format) saved to regression_baseline_affordability.txt\n")
      print(aff_reg_tbl_obj)
    }
  } else {
    cat("Insufficient data for baseline affordability regressions.\n")
  }
}

# =============================================================================
# 9. SUMMARY TABLES BY PERIOD
# =============================================================================
cat("\n--- 9. SUMMARY TABLES BY PERIOD ---\n")
if (nrow(msa_panel) > 0 && "period" %in% names(msa_panel)) {
  period_summary <- msa_panel %>% filter(!is.na(period)) %>% group_by(period) %>%
    summarise(n_msa_years = n(), 
              avg_birth_rate = mean(birth_rate_per_1000, na.rm = TRUE),
              avg_pct_hh_children = mean(pct_households_with_children, na.rm = TRUE),
              avg_median_income = mean(median_hh_income, na.rm = TRUE),
              avg_price_to_income = mean(price_to_income, na.rm = TRUE),
              .groups = 'drop') %>% mutate(across(where(is.numeric), ~round(., 2)))
  save_table(period_summary, "summary_by_period.csv"); print(period_summary)
}

# =============================================================================
# 10. COVID IMPACT ANALYSIS (Simplified: Pre/Post 2020 comparison)
# =============================================================================
cat("\n--- 10. COVID IMPACT ANALYSIS (SIMPLIFIED) ---\n")
if (nrow(msa_panel) > 0 && min(msa_panel$year, na.rm=T) < 2020 && max(msa_panel$year, na.rm=T) >= 2020) {
  covid_data <- msa_panel %>%
    mutate(covid_era = factor(ifelse(year >= 2020, "Post-2019", "Pre-2020"), levels = c("Pre-2020", "Post-2019"))) %>%
    filter(!is.na(birth_rate_per_1000) & !is.na(covid_era))
  
  paired_covid <- covid_data %>% group_by(GEOID, NAME, covid_era) %>%
    summarise(avg_val = mean(birth_rate_per_1000, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = covid_era, values_from = avg_val) %>% drop_na()
  
  if (nrow(paired_covid) >= 10) {
    tt_res <- t.test(paired_covid$`Post-2019`, paired_covid$`Pre-2020`, paired = TRUE)
    save_table(broom::tidy(tt_res), "covid_impact_ttest_birthrate.csv"); print(tt_res)
    
    paired_covid <- paired_covid %>% mutate(change = `Post-2019` - `Pre-2020`)
    p_covid_chg <- ggplot(paired_covid, aes(x = change)) +
      geom_histogram(bins = 20, fill = "tomato", color = "black", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
      geom_vline(xintercept = mean(paired_covid$change, na.rm=T), color = "darkred") +
      annotate("text", x=mean(paired_covid$change, na.rm=T)*1.1, y=Inf, vjust=2, label=paste("Mean:", round(mean(paired_covid$change, na.rm=T),2)), color="darkred") +
      labs(title = "Distribution of Change in Avg. Birth Rate per MSA", subtitle = "Post-2019 Avg minus Pre-2020 Avg", x = "Change in Birth Rate", y = "Number of MSAs")
    save_plot(p_covid_chg, "covid_impact_birth_rate_change_dist.png")
  } else {
    cat("Not enough paired data for COVID impact t-test on birth rates. Descriptive summary instead:\n")
    desc_covid_sum <- covid_data %>% group_by(covid_era) %>% summarise(mean_br = mean(birth_rate_per_1000, na.rm=T), n=n())
    save_table(desc_covid_sum, "covid_era_descriptive_summary_birthrate.csv"); print(desc_covid_sum)
  }
} else {
  cat("Insufficient year range for COVID impact analysis (need data pre and post 2020).\n")
}

# =============================================================================
# 11. NEW ANALYSIS SECTION: QUANTILE-BASED PLOTS 
# =============================================================================
cat("\n--- 11. QUANTILE-BASED PLOTS (Cross-sectional for latest year) ---\n")

outcome_vars_q <- c("pct_households_with_children", "birth_rate_per_1000")
quantile_def_vars_q <- c("total_population", "median_gross_rent", "median_hh_income")
num_quantiles <- 5 # Quintiles

if (!is.na(latest_year_val) && nrow(latest_year_data) > 0) {
  for (ov_q in outcome_vars_q) {
    if (!ov_q %in% names(latest_year_data)) {
      cat(sprintf("Outcome variable '%s' not found in latest_year_data. Skipping for quantile plots.\n", ov_q))
      next
    }
    for (qv_q in quantile_def_vars_q) {
      if (!qv_q %in% names(latest_year_data)) {
        cat(sprintf("Quantile-defining variable '%s' not found in latest_year_data. Skipping for quantile plots.\n", qv_q))
        next
      }
      
      cat(sprintf("Generating quantile plot for DV: %s by IV: %s quantiles (%s)\n", ov_q, qv_q, latest_year_val))
      
      quantile_data <- latest_year_data %>%
        filter(!is.na(.data[[ov_q]]) & !is.na(.data[[qv_q]]) & !is.na(total_population)) %>%
        mutate(quantile_group = ntile(.data[[qv_q]], num_quantiles)) %>%
        filter(!is.na(quantile_group))
      
      if (nrow(quantile_data) < num_quantiles * 5) { 
        cat(sprintf("Insufficient data for %s by %s quantiles.\n", ov_q, qv_q))
        next
      }
      
      summary_by_quantile <- quantile_data %>%
        group_by(quantile_group) %>%
        summarise(
          avg_outcome = weighted.mean(.data[[ov_q]], w = total_population, na.rm = TRUE),
          n_msas = n(),
          .groups = "drop"
        ) %>%
        mutate(quantile_group = factor(quantile_group, 
                                       labels = paste0("Q", 1:num_quantiles, " (", qv_q, ")"))) 
      
      if(nrow(summary_by_quantile) > 0){
        plot_title <- paste0("Average ", gsub("_", " ", ov_q), " by ", gsub("_", " ", qv_q), 
                             sprintf(" %s-tiles (%s %s)", ifelse(num_quantiles==5, "Quin", ifelse(num_quantiles==10, "Dec", as.character(num_quantiles))), "Year", latest_year_val))
        
        p_quantile <- ggplot(summary_by_quantile, aes(x = quantile_group, y = avg_outcome, fill = quantile_group)) +
          geom_col(alpha = 0.8, show.legend = FALSE) +
          geom_text(aes(label = sprintf("%.2f (N=%d)", avg_outcome, n_msas)), vjust = -0.5, size = 3) +
          scale_fill_viridis_d() +
          labs(title = plot_title,
               x = paste(tools::toTitleCase(gsub("_", " ", qv_q)), "Quantile Group"),
               y = paste("Avg.", tools::toTitleCase(gsub("_", " ", ov_q))),
               caption = "Averages are population-weighted. Quintiles based on MSAs in the latest year.") +
          theme(axis.text.x = element_text(angle = 15, hjust = 1))
        
        if(ov_q == "pct_households_with_children") {
          p_quantile <- p_quantile + scale_y_continuous(labels = scales::percent_format(scale = 1))
        }
        
        save_plot(p_quantile, paste0(ov_q, "_by_", qv_q, "_quantiles_", latest_year_val, ".png"), width = 10, height = 7)
      } else {
        cat(sprintf("No data to plot for %s by %s quantiles after summarization.\n", ov_q, qv_q))
      }
    }
  }
} else {
  cat("Skipping Quantile-Based Plots: latest_year_val is NA or latest_year_data is empty.\n")
}

# =============================================================================
# 12. MODIFIED ANALYSIS SECTION: PLOTTING REGRESSION COEFFICIENTS OVER TIME
# =============================================================================
cat("\n--- 12. PLOTTING REGRESSION COEFFICIENTS OVER TIME (EXTENDED) ---\n")

dvs_time_coef <- c("birth_rate_per_1000", "pct_households_with_children")
ivs_yearly_reg_std_full <- paste0("std_", c("log_median_hh_income", "price_to_income", "annual_rent_to_income", "pct_college_plus", "log_total_population"))
key_ivs_to_plot_ts <- c("std_log_median_hh_income", "std_price_to_income", "std_log_total_population", "std_annual_rent_to_income")

if (!exists("reg_data_base") || nrow(reg_data_base) == 0) {
  stop("reg_data_base not found or empty. Please ensure Section 5 runs correctly.")
}
missing_cols_in_reg_data_base <- setdiff(c(dvs_time_coef, ivs_yearly_reg_std_full), names(reg_data_base))
if (length(missing_cols_in_reg_data_base) > 0) {
  stop(paste("reg_data_base is missing required columns for Section 12:", paste(missing_cols_in_reg_data_base, collapse = ", ")))
}

all_yearly_coefs_list <- list() 
unique_years <- sort(unique(reg_data_base$year))

if (length(unique_years) > 1 && nrow(reg_data_base) > 0) {
  for (dv_tc in dvs_time_coef) {
    cat(sprintf("\nProcessing DV for Time-Varying Coefs: %s\n", dv_tc))
    yearly_coefs_for_dv_list <- list()
    
    for (yr in unique_years) {
      cat(sprintf("  Year: %d\n", yr))
      
      year_data_for_reg <- reg_data_base %>% 
        filter(year == yr) %>%
        select(all_of(dv_tc), all_of(ivs_yearly_reg_std_full))
      
      year_data_full_model <- year_data_for_reg %>% drop_na()
      
      if (nrow(year_data_full_model) > (length(ivs_yearly_reg_std_full) + 10)) {
        formula_full_str <- paste(dv_tc, "~", paste(ivs_yearly_reg_std_full, collapse = " + "))
        model_full <- tryCatch(lm(as.formula(formula_full_str), data = year_data_full_model), error = function(e) {
          cat(sprintf("    Error (Full Model, DV: %s, Year: %d): %s\n", dv_tc, yr, e$message)); NULL })
        
        if (!is.null(model_full)) {
          tidy_full <- broom::tidy(model_full, conf.int = TRUE) %>%
            filter(term %in% key_ivs_to_plot_ts) %>%
            mutate(
              year = yr, 
              dv = dv_tc, 
              model_type = "full_controls",
              controls = map_chr(term, ~ paste(setdiff(ivs_yearly_reg_std_full, .x), collapse = ", "))
            )
          if(nrow(tidy_full) > 0) yearly_coefs_for_dv_list[[paste0("Y",yr, "_DV",dv_tc, "_Full")]] <- tidy_full
        }
      } else {
        cat(sprintf("    Skipping Full Model (DV: %s, Year: %d): Insufficient observations (%d).\n", dv_tc, yr, nrow(year_data_full_model)))
      }
      
      for (key_iv_bivariate in key_ivs_to_plot_ts) {
        year_data_bivariate <- year_data_for_reg %>%
          select(all_of(dv_tc), all_of(key_iv_bivariate)) %>%
          drop_na()
        
        if (nrow(year_data_bivariate) > (1 + 10)) { 
          formula_bivariate_str <- paste(dv_tc, "~", key_iv_bivariate)
          model_bivariate <- tryCatch(lm(as.formula(formula_bivariate_str), data = year_data_bivariate), error = function(e) {
            cat(sprintf("    Error (Bivariate Model, DV: %s, IV: %s, Year: %d): %s\n", dv_tc, key_iv_bivariate, yr, e$message)); NULL})
          
          if (!is.null(model_bivariate)) {
            tidy_bivariate <- broom::tidy(model_bivariate, conf.int = TRUE) %>%
              filter(term == key_iv_bivariate) %>% 
              mutate(year = yr, dv = dv_tc, model_type = "no_controls", controls = "none")
            if(nrow(tidy_bivariate) > 0) yearly_coefs_for_dv_list[[paste0("Y",yr, "_DV",dv_tc, "_Bivar_", gsub("std_","",key_iv_bivariate))]] <- tidy_bivariate
          }
        } else {
          cat(sprintf("    Skipping Bivariate Model (DV: %s, IV: %s, Year: %d): Insufficient observations (%d).\n", dv_tc, key_iv_bivariate, yr, nrow(year_data_bivariate)))
        }
      }
    } 
    if (length(yearly_coefs_for_dv_list) > 0) {
      all_yearly_coefs_list[[dv_tc]] <- bind_rows(yearly_coefs_for_dv_list)
    }
  } 
  
  if (length(all_yearly_coefs_list) > 0) {
    final_yearly_coefs_table <- bind_rows(all_yearly_coefs_list)
    save_table(final_yearly_coefs_table, "yearly_regression_coefficients_ext.csv")
    cat("\nSaved extended yearly regression coefficients to yearly_regression_coefficients_ext.csv\n")
    
    cat("\nGenerating time-series coefficient plots...\n")
    if (!all(c("dv", "term", "model_type") %in% names(final_yearly_coefs_table))) {
      cat("Error: final_yearly_coefs_table is missing dv, term, or model_type columns. Cannot generate plots.\n")
    } else {
      plot_combinations <- final_yearly_coefs_table %>% 
        distinct(dv, term, model_type)
      
      for (i in 1:nrow(plot_combinations)) {
        dv_name_plot <- plot_combinations$dv[i]
        iv_name_plot <- plot_combinations$term[i]
        m_type <- plot_combinations$model_type[i]
        
        plot_data_segment <- final_yearly_coefs_table %>%
          filter(dv == dv_name_plot, term == iv_name_plot, model_type == m_type)
        
        if (nrow(plot_data_segment) > 1) { 
          iv_name_plot_clean <- tools::toTitleCase(gsub("_", " ", gsub("std_", "Std. ", iv_name_plot)))
          dv_name_plot_clean <- tools::toTitleCase(gsub("_", " ", dv_name_plot))
          model_type_clean <- ifelse(m_type == "full_controls", "Full Controls", "No Controls (Bivariate)")
          
          plot_title_str <- paste("Time Series of Coefficient:", iv_name_plot_clean)
          plot_subtitle_str <- paste("DV:", dv_name_plot_clean, "| Model:", model_type_clean, "| Yearly OLS")
          
          p_coef_time <- ggplot(plot_data_segment, aes(x = year, y = estimate)) +
            geom_line(color = "blue", linewidth = 1) +
            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") + 
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
            labs(title = plot_title_str,
                 subtitle = plot_subtitle_str,
                 x = "Year", y = "Coefficient Estimate (with 95% CI)") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + 
            theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))
          
          filename_short_iv <- gsub("std_","",iv_name_plot) 
          filename_short_iv <- gsub("log_","",filename_short_iv)
          filename_short_iv <- gsub("_","" ,filename_short_iv)
          
          filename_short_dv <- gsub("_per_1000","",dv_name_plot)
          filename_short_dv <- gsub("_","" ,filename_short_dv)
          filename_short_dv <- gsub("pcthouseholds","pcthh",filename_short_dv)
          
          plot_filename <- paste0("coef_ts_", filename_short_iv, "_on_", filename_short_dv, "_", m_type, ".png")
          save_plot(p_coef_time, plot_filename, width = 11, height = 6.5) 
        } else {
          cat(sprintf("  Skipping plot for DV: %s, IV: %s, Model: %s - Not enough data points (%d).\n", dv_name_plot, iv_name_plot, m_type, nrow(plot_data_segment)))
        }
      }
    }
  } else {
    cat("No yearly coefficients were successfully collected for plotting.\n")
  }
} else {
  cat("Skipping time-varying coefficient plots: Not enough unique years or data in reg_data_base.\n")
}


# =============================================================================
# 13. FINAL SUMMARY OF ANALYSIS 
# =============================================================================
cat("\n--- 13. ANALYSIS SCRIPT COMPLETE ---\n")
cat(paste0("All outputs saved in: ", output_dir_analysis, "\n"))

summary_text_parts <- list(
  "MSA PANEL DATA ANALYSIS SUMMARY (v2.1 - Extended v2 - INTEGRAL)",
  "===============================================================",
  paste0("Data Source: ", data_file_path),
  paste0("Total Observations in Processed Panel: ", if(exists("msa_panel")) nrow(msa_panel) else "N/A"),
  paste0("Unique MSAs: ", if(exists("msa_panel")) n_distinct(msa_panel$GEOID) else "N/A"),
  paste0("Year Range: ", if(exists("msa_panel") && nrow(msa_panel)>0) min(msa_panel$year, na.rm = TRUE) else "N/A", 
         " - ", if(exists("msa_panel") && nrow(msa_panel)>0) max(msa_panel$year, na.rm = TRUE) else "N/A"),
  "\nKey Variables Analyzed:",
  "- Fertility: birth_rate_per_1000, child_woman_ratio",
  "- Household Composition: pct_households_with_children",
  "- Demographics: total_population, pct_under_18, pct_college_plus",
  "- Housing/Economics: median_hh_income, price_to_income, annual_rent_to_income, median_gross_rent",
  "\nAnalysis Sections Performed:",
  "1. Data Loading & Initial Inspection", 
  "2. Data Cleaning & Preparation",
  "3. Exploratory Data Analysis",
  "4. Correlation Analysis", 
  "5. Panel Regression (Main table as .txt, birth rate coef plot)",
  "6. Heterogeneity Analysis",
  "7. Time Series for Selected MSAs",
  "8. Housing Affordability (Baseline affordability reg table as .txt)",
  "9. Summary Tables by Period", 
  "10. COVID Impact Analysis (Simplified)",
  "11. Quantile-Based Plots (Cross-sectional for latest year)",
  "12. Plotting Regression Coefficients Over Time (Extended):",
  "    - Yearly cross-sectional regressions for key DVs.",
  "    - Time-series plots for coefficients of: std_log_median_hh_income, std_price_to_income, std_log_total_population, std_annual_rent_to_income.",
  "    - Each IV coefficient plotted under two scenarios: 'full controls' and 'no controls (bivariate)'.",
  "    - Coefficients saved to yearly_regression_coefficients_ext.csv.",
  "13. Final Summary of Analysis",
  "\nNote: Some analyses/plots might be skipped if data is insufficient."
)
summary_final_text <- paste(summary_text_parts, collapse = "\n")
writeLines(summary_final_text, file.path(output_dir_analysis, "analysis_run_summary_v2.1_extended_v2.txt")) 
cat("\n", summary_final_text, "\n")