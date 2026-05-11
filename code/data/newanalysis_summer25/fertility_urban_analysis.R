# =============================================================================
# FERTILITY, HOUSING, AND URBAN AGGLOMERATION: COMPREHENSIVE ANALYSIS
# =============================================================================
# This script documents key patterns in fertility across the urban hierarchy,
# connecting to the agglomeration literature and creating publication-ready output

# library(tidycensus) # Make sure this is installed and working
# library(tidyverse)
# library(scales)
# library(ggrepel)
# library(viridis)
# library(patchwork)
# library(knitr)
# library(modelsummary) # If you want to use it for regression tables later

# For CRAN checks or non-interactive use, ensure libraries are explicitly called or attached
if (require(tidycensus) & require(tidyverse) & require(scales) & require(ggrepel) & 
    require(viridis) & require(patchwork) & require(knitr)) {
  
  # Set theme for all plots
  theme_set(theme_minimal(base_size = 12) + 
              theme(panel.grid.minor = element_blank(),
                    plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 12, color = "gray50")))
  
  # =============================================================================
  # PART 1: DATA COLLECTION - MULTIPLE YEARS FOR TIME TRENDS
  # =============================================================================
  
  # Function to get comprehensive MSA data for a given year
  get_msa_fertility_data <- function(year) {
    
    cat("Fetching data for year", year, "\n")
    
    # Skip if year is before ACS availability (DP tables are generally available from 2009/2010 for 5-year)
    if(year < 2010) { # Adjusting earliest year for 5-year ACS
      cat("Year", year, "is too early for reliable ACS 5-year DP tables. Skipping.\n")
      return(NULL)
    }
    
    survey <- "acs5" # Using 5-year estimates for more robust MSA data
    
    tryCatch({
      # Get household and family data (DP02 table)
      household_data <- get_acs(
        geography = "metropolitan statistical area/micropolitan statistical area",
        table = "DP02",
        year = year,
        survey = survey,
        cache_table = TRUE
      )
      
      # Process DP02 data
      dp02_processed <- household_data %>%
        filter(
          variable %in% c(
            "DP02_0001",  # Total households
            "DP02_0002",  # Family households  
            "DP02_0003",  # With own children under 18
            "DP02_0016"   # Average household size
          )
        ) %>%
        select(GEOID, NAME, variable, estimate) %>%
        pivot_wider(names_from = variable, values_from = estimate) %>%
        mutate(
          share_households_with_children = if("DP02_0003" %in% names(.) && "DP02_0001" %in% names(.)) DP02_0003 / DP02_0001 else NA,
          share_families_with_children = if("DP02_0003" %in% names(.) && "DP02_0002" %in% names(.)) DP02_0003 / DP02_0002 else NA,
          avg_household_size = if("DP02_0016" %in% names(.)) DP02_0016 else NA
        )
      
      # Get economic indicators
      econ_vars <- c(
        "B25077_001", # Median home value
        "B19013_001", # Median household income
        "B25064_001", # Median gross rent
        "B01003_001"  # Total population
      )
      
      econ_data <- get_acs(
        geography = "metropolitan statistical area/micropolitan statistical area",
        variables = econ_vars,
        year = year,
        survey = survey,
        output = "wide",
        cache_table = TRUE
      ) %>%
        rename(
          median_home_value = B25077_001E,
          median_household_income = B19013_001E,
          median_rent = B25064_001E,
          population = B01003_001E
        )
      
      # For education data
      edu_processed <- NULL # Initialize
      # Table B15003 is generally available, but let's be safe
      tryCatch({
        edu_data_raw <- get_acs(
          geography = "metropolitan statistical area/micropolitan statistical area",
          table = "B15003", # Educational Attainment for Pop 25+
          year = year,
          survey = survey,
          cache_table = TRUE
        )
        
        edu_processed <- edu_data_raw %>%
          filter(variable %in% c("B15003_001", "B15003_022", 
                                 "B15003_023", "B15003_024", "B15003_025")) %>%
          select(GEOID, variable, estimate) %>%
          pivot_wider(names_from = variable, values_from = estimate) %>%
          mutate(
            college_pct = (B15003_022 + B15003_023 + B15003_024 + B15003_025) / B15003_001 * 100
          ) %>%
          select(GEOID, college_pct)
      }, error = function(e) {
        cat("Could not get detailed education data (B15003) for", year, ":", e$message, "\n")
        # Fallback to DP02 education data if B15003 fails
        tryCatch({
          dp02_edu_vars <- c("DP02_0059", "DP02_0067") # Total pop 25+, Bachelor's degree or higher
          dp02_edu_data <- get_acs(
            geography = "metropolitan statistical area/micropolitan statistical area",
            variables = dp02_edu_vars,
            year = year,
            survey = survey,
            output = "wide",
            cache_table = TRUE
          )
          edu_processed <- dp02_edu_data %>%
            mutate(college_pct = DP02_0067E / DP02_0059E * 100) %>%
            select(GEOID, college_pct)
          cat("Using DP02 education data for", year, "\n")
        }, error = function(e_dp){
          cat("Could not get any education data for", year, ":", e_dp$message, "\n")
          # edu_processed remains NULL
        })
      })
      
      # Combine all data
      msa_data <- dp02_processed %>%
        left_join(econ_data, by = c("GEOID", "NAME"))
      
      # Add education if available
      if(!is.null(edu_processed)) {
        msa_data <- msa_data %>%
          left_join(edu_processed, by = "GEOID")
      } else {
        msa_data$college_pct <- NA
      }
      
      # Process final dataset
      msa_data <- msa_data %>%
        mutate(
          msa_name = str_extract(NAME, "^[^,]+"),
          cbsa_code = as.numeric(str_extract(GEOID, "\\d+$")), # GEOID is like "31000USXXXXX"
          year = as.integer(year),
          price_to_income = median_home_value / median_household_income,
          annual_rent = median_rent * 12,
          rent_to_income = annual_rent / median_household_income
        ) %>%
        filter(!str_detect(NAME, ", PR"), # Remove Puerto Rico
               population > 50000, # Keep only MSAs (or larger Micropolitans)
               !is.na(share_families_with_children),
               !is.na(population) 
        ) %>%
        select(GEOID, NAME, msa_name, cbsa_code, year, 
               starts_with("share_"), population, 
               starts_with("median_"), price_to_income, 
               rent_to_income, college_pct, avg_household_size)
      
      return(msa_data)
      
    }, error = function(e) {
      warning(paste("Error fetching or processing data for year", year, ":", e$message))
      return(NULL)
    })
  }
  
  # Fetch data for multiple years - recent years are more stable
  # ACS 5-year data ending in these years
  years <- c(2013, 2015, 2017, 2019, 2021, 2022) 
  msa_data_list <- lapply(years, get_msa_fertility_data)
  
  msa_data_list <- msa_data_list[!sapply(msa_data_list, is.null)]
  
  if(length(msa_data_list) == 0) {
    stop("No data could be fetched. Check API key, internet connection, or ACS availability for selected years.")
  }
  
  msa_panel <- bind_rows(msa_data_list)
  
  current_year <- max(msa_panel$year)
  msa_current <- msa_panel %>% filter(year == current_year)
  
  # =============================================================================
  # PART 2: CREATE MSA CATEGORIES AND AGGLOMERATION MEASURES
  # =============================================================================
  
  msa_panel <- msa_panel %>%
    mutate(
      log_population = log(population)
    ) %>%
    group_by(year) %>%
    mutate(
      size_category = case_when(
        population >= 5000000 ~ "Mega (>5M)",
        population >= 2500000 ~ "Very Large (2.5-5M)",
        population >= 1000000 ~ "Large (1-2.5M)",
        population >= 500000 ~ "Medium (0.5-1M)",
        population >= 250000 ~ "Small (250-500K)",
        TRUE ~ "Very Small (<250K)"
      ),
      size_category = factor(size_category, 
                             levels = c("Very Small (<250K)", "Small (250-500K)", 
                                        "Medium (0.5-1M)", "Large (1-2.5M)", 
                                        "Very Large (2.5-5M)", "Mega (>5M)")),
      pop_rank = rank(-population),
      is_top_10 = pop_rank <= 10,
      is_top_20 = pop_rank <= 20,
      is_top_50 = pop_rank <= 50,
      relative_pop = population / mean(population, na.rm = TRUE),
      relative_income = median_household_income / mean(median_household_income, na.rm = TRUE),
      relative_price = median_home_value / mean(median_home_value, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Re-filter msa_current with the new columns
  msa_current <- msa_panel %>% filter(year == current_year)
  
  # =============================================================================
  # PART 3: DOCUMENT CORE STYLIZED FACTS
  # =============================================================================
  
  fig1_data <- msa_current %>%
    filter(!is.na(share_families_with_children), !is.na(log_population))
  
  if(nrow(fig1_data) > 0) {
    p1 <- ggplot(fig1_data, aes(x = log_population, y = share_families_with_children)) +
      geom_point(aes(size = population, color = price_to_income), alpha = 0.7) +
      geom_smooth(method = "lm", aes(weight = population), color = "darkred", se = TRUE) +
      geom_text_repel(data = fig1_data %>% filter(is_top_20 | share_families_with_children < 
                                                    quantile(share_families_with_children, 0.1, na.rm=TRUE) |
                                                    price_to_income > quantile(price_to_income, 0.9, na.rm=TRUE)), # Adjusted filter for labeling
                      aes(label = msa_name), size = 3, max.overlaps = 15) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      scale_size_continuous(range = c(2, 10), guide = "none") +
      scale_color_viridis_c(name = "House Price\nto Income") +
      labs(
        title = "The Urban Fertility Penalty",
        subtitle = paste("Share of families with children vs. city size (", current_year, ")", sep = ""),
        x = "Log Population",
        y = "Share of Families with Children",
        caption = "Source: ACS 5-Year Estimates. Weighted regression line by MSA population."
      )
    print(p1)
  } else {
    cat("Not enough data for Figure 1 (Urban Fertility Penalty).\n")
    p1 <- NULL
  }
  
  time_trends_data <- msa_panel %>%
    filter(!is.na(share_families_with_children), !is.na(size_category)) %>%
    group_by(year, size_category) %>%
    summarise(
      avg_fertility = weighted.mean(share_families_with_children, w = population, na.rm = TRUE),
      .groups = "drop"
    )
  
  if(nrow(time_trends_data) > 0 && length(unique(time_trends_data$year)) > 1) {
    p2 <- ggplot(time_trends_data, aes(x = year, y = avg_fertility, color = size_category)) +
      geom_line(linewidth = 1.2) + # Use linewidth for ggplot2 3.4.0+
      geom_point(size = 3) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      scale_x_continuous(breaks = unique(time_trends_data$year)) +
      scale_color_viridis_d(name = "MSA Size") +
      labs(
        title = "Diverging Fertility Trends by City Size",
        subtitle = "The gap between large and small cities may be changing",
        x = "Year",
        y = "Average Share of Families with Children",
        caption = "Population-weighted averages within each size category"
      )
    print(p2)
  } else {
    cat("Not enough data for Figure 2 (Time Trends).\n")
    p2 <- NULL
  }
  
  # =============================================================================
  # PART 4: THE HOUSING COST CHANNEL
  # =============================================================================
  
  fig3_data <- msa_current %>%
    filter(!is.na(price_to_income), !is.na(share_families_with_children))
  
  if(nrow(fig3_data) > 0) {
    p3 <- ggplot(fig3_data, aes(x = price_to_income, y = share_families_with_children)) +
      geom_point(aes(size = population, color = log_population), alpha = 0.7) +
      geom_smooth(method = "lm", aes(weight = population), color = "darkred", se = TRUE) +
      geom_smooth(method = "loess", aes(weight = population), color = "darkblue", 
                  se = FALSE, linetype = "dashed", span = 0.75) + # Adjusted span
      geom_text_repel(data = fig3_data %>% 
                        filter(is_top_10 | price_to_income > quantile(price_to_income, 0.9, na.rm=TRUE) | 
                                 share_families_with_children < quantile(share_families_with_children, 0.1, na.rm=TRUE)), # Adjusted filter
                      aes(label = msa_name), size = 3, max.overlaps = 20) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      scale_size_continuous(range = c(2, 10), guide = "none") +
      scale_color_viridis_c(name = "Log Pop", option = "plasma") +
      labs(
        title = "Housing Affordability and Family Formation",
        subtitle = "Higher housing costs often correlate with lower family-with-children rates",
        x = "House Price to Income Ratio",
        y = "Share of Families with Children",
        caption = "Red line: linear fit; Blue dashed: local polynomial"
      )
    print(p3)
  } else {
    cat("Not enough data for Figure 3 (Housing Cost Channel).\n")
    p3 <- NULL
  }
  
  # =============================================================================
  # PART 5: AGGLOMERATION-FERTILITY TRADE-OFF
  # =============================================================================
  
  msa_current_for_agglom <- msa_current # Make a copy to add earnings
  tryCatch({
    earnings_data <- get_acs(
      geography = "metropolitan statistical area/micropolitan statistical area",
      variables = "B20002_001",  # Median earnings for workers
      year = current_year,
      survey = "acs5",
      output = "wide",
      cache_table = TRUE
    ) %>%
      rename(median_earnings = B20002_001E) %>%
      select(GEOID, median_earnings)
    
    msa_current_for_agglom <- msa_current_for_agglom %>%
      left_join(earnings_data, by = "GEOID")
    
  }, error = function(e) {
    cat("Could not get earnings data (B20002_001), using household income as proxy for earnings premium.\n")
    cat("Note: This is a rough proxy. Median earnings are preferred.\n")
    msa_current_for_agglom <- msa_current_for_agglom %>%
      mutate(median_earnings = median_household_income * 0.7) # A very rough adjustment factor
  })
  
  agglomeration_analysis <- msa_current_for_agglom %>%
    filter(!is.na(median_earnings), !is.na(size_category), !is.na(share_families_with_children)) %>%
    mutate(
      is_baseline = size_category == "Small (250-500K)"
    )
  
  if(any(agglomeration_analysis$is_baseline)) {
    baseline_values <- agglomeration_analysis %>%
      filter(is_baseline) %>%
      summarise(
        baseline_wage = weighted.mean(median_earnings, w = population, na.rm = TRUE),
        baseline_fertility = weighted.mean(share_families_with_children, w = population, na.rm = TRUE)
      )
    
    agglomeration_analysis <- agglomeration_analysis %>%
      mutate(
        wage_premium = (median_earnings / baseline_values$baseline_wage - 1) * 100,
        fertility_penalty = (share_families_with_children - baseline_values$baseline_fertility) * 100 # This is a difference, not a penalty in direction
      )
    
    agglomeration_summary <- agglomeration_analysis %>%
      group_by(size_category) %>%
      summarise(
        n = n(),
        total_pop = sum(population, na.rm=TRUE) / 1e6,
        avg_wage_premium = weighted.mean(wage_premium, w = population, na.rm = TRUE),
        avg_fertility_diff = weighted.mean(fertility_penalty, w = population, na.rm = TRUE), # Renamed for clarity
        avg_college = weighted.mean(college_pct, w = population, na.rm = TRUE),
        .groups = "drop"
      )
    
    if(nrow(agglomeration_summary) > 0) {
      p4_data <- agglomeration_summary %>%
        select(size_category, avg_wage_premium, avg_fertility_diff) %>%
        pivot_longer(cols = c(avg_wage_premium, avg_fertility_diff),
                     names_to = "measure", values_to = "value") %>%
        mutate(
          measure = factor(measure, 
                           levels = c("avg_wage_premium", "avg_fertility_diff"),
                           labels = c("Wage Premium", "Fertility Difference"))
        )
      
      p4 <- ggplot(p4_data, aes(x = size_category, y = value, fill = measure)) +
        geom_col(position = "dodge", width = 0.7) +
        geom_hline(yintercept = 0, linetype = "solid", color = "gray50") +
        scale_fill_manual(values = c("Wage Premium" = "#2E86AB", "Fertility Difference" = "#A23B72")) +
        scale_y_continuous(labels = function(x) paste0(x, "%")) +
        coord_flip() +
        labs(
          title = "The Urban Trade-off: Wages vs. Fertility",
          subtitle = "Percentage difference from small MSAs (250-500K population)",
          x = NULL,
          y = "Difference from Small MSAs (%)",
          fill = NULL
        ) +
        theme(legend.position = "bottom")
      print(p4)
    } else {
      cat("Not enough data for Figure 4 (Agglomeration Trade-off).\n")
      p4 <- NULL
      agglomeration_summary <- NULL # ensure it's NULL if not created
    }
  } else {
    cat("No baseline MSAs (Small 250-500K) found for agglomeration analysis. Skipping Figure 4.\n")
    p4 <- NULL
    agglomeration_summary <- NULL
  }
  
  # =============================================================================
  # PART 6: REGRESSION ANALYSIS
  # =============================================================================
  
  reg_data <- msa_current %>%
    filter(!is.na(share_families_with_children),
           !is.na(median_household_income),
           !is.na(price_to_income),
           !is.na(log_population)) %>%
    mutate(
      log_income = log(median_household_income),
      log_home_value = log(median_home_value),
      log_rent = log(median_rent)
    )
  
  if(nrow(reg_data) > 30) { # Need enough observations for regression
    m1 <- lm(share_families_with_children ~ log_population, 
             data = reg_data, weights = population)
    m2 <- lm(share_families_with_children ~ log_population + price_to_income, 
             data = reg_data, weights = population)
    m3 <- lm(share_families_with_children ~ log_population + price_to_income + 
               log_income, data = reg_data, weights = population)
    
    if(sum(!is.na(reg_data$college_pct)) > 30 && "college_pct" %in% names(reg_data)) { # Check if column exists and has enough data
      m4 <- lm(share_families_with_children ~ log_population + price_to_income + 
                 log_income + college_pct, data = reg_data, weights = population)
    } else {
      cat("Skipping college_pct in Model 4 due to insufficient data.\n")
      m4 <- m3 # Fallback
    }
    
    if (all(c("log_home_value", "log_rent") %in% names(reg_data))) {
      m5 <- lm(share_families_with_children ~ log_population + log_home_value + 
                 log_rent + log_income, 
               data = reg_data, weights = population)
    } else {
      cat("Skipping Model 5 due to missing log_home_value or log_rent.\n")
      m5 <- m4 # Fallback
    }
    
    cat("\n=== REGRESSION RESULTS (Coefficients and Std. Errors) ===\n")
    models_list <- list(m1=m1, m2=m2, m3=m3, m4=m4, m5=m5)
    for (i in seq_along(models_list)) {
      cat(paste0("\n--- ", names(models_list)[i], " ---\n"))
      print(summary(models_list[[i]])$coefficients[,1:2])
    }
    
    # If modelsummary is available:
    # if (require(modelsummary)) {
    #   modelsummary_table <- modelsummary(models_list,
    #                                      stars = TRUE,
    #                                      gof_map = c("nobs", "r.squared"))
    #   print(modelsummary_table)
    # }
    
  } else {
    cat("Insufficient data for regression analysis.\n")
    m1 <- m2 <- m3 <- m4 <- m5 <- NULL
  }
  
  # =============================================================================
  # PART 7: CREATE SUMMARY TABLES
  # =============================================================================
  
  if(nrow(msa_current) > 0 && "size_category" %in% names(msa_current)) {
    summary_table_data <- msa_current %>%
      group_by(size_category) %>%
      summarise(
        `Number of MSAs` = n(),
        `Total Population (millions)` = round(sum(population, na.rm=TRUE) / 1e6, 1),
        `Share with Children (%)` = round(weighted.mean(share_families_with_children, 
                                                        w = population, na.rm = TRUE) * 100, 1),
        `Median Income ($)` = dollar(round(weighted.mean(median_household_income, 
                                                         w = population, na.rm = TRUE), -2)),
        `House Price/Income` = round(weighted.mean(price_to_income, 
                                                   w = population, na.rm = TRUE), 2),
        `Rent/Income (%)` = round(weighted.mean(rent_to_income, 
                                                w = population, na.rm = TRUE) * 100, 1),
        .groups = "drop"
      )
    cat("\n=== SUMMARY STATISTICS BY CITY SIZE (Current Year: ", current_year, ") ===\n")
    print(kable(summary_table_data, format = "pipe")) # Use pipe for better markdown rendering
  } else {
    cat("Not enough data for summary table by city size.\n")
    summary_table_data <- NULL
  }
  
  if(length(unique(msa_panel$year)) >= 2 && "size_category" %in% names(msa_panel)) {
    min_yr <- min(msa_panel$year)
    max_yr <- max(msa_panel$year)
    
    time_change_table_data <- msa_panel %>%
      filter(year %in% c(min_yr, max_yr)) %>%
      group_by(year, size_category) %>%
      summarise(
        avg_fertility = weighted.mean(share_families_with_children, w = population, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_wider(names_from = year, values_from = avg_fertility, 
                  names_prefix = "year_") %>%
      mutate(
        change_ppt = (.[[paste0("year_", max_yr)]] - .[[paste0("year_", min_yr)]]) * 100
      )
    cat("\n=== CHANGE IN SHARE OF FAMILIES WITH CHILDREN (", min_yr, "to", max_yr, ") ===\n")
    print(kable(time_change_table_data, format = "pipe", digits=2))
  } else {
    cat("Not enough data for time change table.\n")
    time_change_table_data <- NULL
  }
  
  # =============================================================================
  # PART 8: SAVE ALL OUTPUTS
  # =============================================================================
  
  dir.create("fertility_analysis_outputs", showWarnings = FALSE)
  
  if(!is.null(p1)) ggsave("fertility_analysis_outputs/fig1_urban_fertility_penalty.png", p1, width = 10, height = 8, dpi = 300, bg="white")
  if(!is.null(p2)) ggsave("fertility_analysis_outputs/fig2_time_trends.png", p2, width = 10, height = 6, dpi = 300, bg="white")
  if(!is.null(p3)) ggsave("fertility_analysis_outputs/fig3_housing_fertility.png", p3, width = 10, height = 8, dpi = 300, bg="white")
  if(!is.null(p4)) ggsave("fertility_analysis_outputs/fig4_agglomeration_tradeoff.png", p4, width = 10, height = 6, dpi = 300, bg="white")
  
  write_csv(msa_panel, "fertility_analysis_outputs/msa_panel_data.csv")
  if(!is.null(summary_table_data)) write_csv(summary_table_data, "fertility_analysis_outputs/summary_by_size.csv")
  if(!is.null(agglomeration_summary)) write_csv(agglomeration_summary, "fertility_analysis_outputs/agglomeration_summary.csv")
  if(!is.null(time_change_table_data)) write_csv(time_change_table_data, "fertility_analysis_outputs/time_change_summary.csv")
  
  if(!is.null(m1)) { # Check if models were created
    save(m1, m2, m3, m4, m5, file = "fertility_analysis_outputs/regression_models.RData")
  }
  
  # =============================================================================
  # PART 9: KEY FINDINGS SUMMARY
  # =============================================================================
  
  cat("\n=== KEY FINDINGS (Based on current year: ", current_year, ") ===\n\n")
  
  if(!is.null(agglomeration_summary) && nrow(agglomeration_summary) > 0) {
    largest_vs_smallest <- agglomeration_summary %>%
      filter(size_category %in% c("Mega (>5M)", "Small (250-500K)")) %>%
      select(size_category, avg_fertility_diff) # Using renamed variable
    
    if(nrow(largest_vs_smallest) == 2 && "Mega (>5M)" %in% largest_vs_smallest$size_category) {
      penalty <- largest_vs_smallest$avg_fertility_diff[largest_vs_smallest$size_category == "Mega (>5M)"]
      cat("1. URBAN FERTILITY DIFFERENCE (vs. Small MSAs):\n")
      cat(sprintf("   - Mega cities (>5M) have %.1f percentage points %s families with children\n", 
                  abs(penalty), ifelse(penalty < 0, "fewer", "more")))
    }
  }
  
  if(!is.null(m1) && !is.null(m2)) { # Check if models exist
    housing_effect <- summary(m2)$coefficients["price_to_income", "Estimate"]
    cat("\n2. HOUSING COST CHANNEL (from regression m2):\n")
    cat(sprintf("   - A 1-unit increase in price/income ratio is associated with a %.1f ppt change in families w/ children\n", 
                housing_effect * 100)) # Sign indicates direction
    
    # Robust check for coefficient existence before division
    coef_m1_log_pop <- tryCatch(coef(m1)["log_population"], error = function(e) NA)
    coef_m2_log_pop <- tryCatch(coef(m2)["log_population"], error = function(e) NA)
    
    if (!is.na(coef_m1_log_pop) && !is.na(coef_m2_log_pop) && coef_m1_log_pop != 0) {
      explained_pct <- (1 - coef_m2_log_pop / coef_m1_log_pop) * 100
      cat(sprintf("   - Adding housing costs (price_to_income) changes the log_population coefficient by ~%.0f%%\n",
                  explained_pct))
    } else {
      cat("   - Could not calculate the percentage change in log_population coefficient.\n")
    }
  }
  
  cat("\n3. DATA SUMMARY:\n")
  cat(sprintf("   - Total MSAs analyzed in current year (%d): %d\n", current_year, nrow(msa_current)))
  cat(sprintf("   - Years covered in panel: %s\n", paste(sort(unique(msa_panel$year)), collapse = ", ")))
  cat(sprintf("   - Population covered in current year data: %.1f million\n", sum(msa_current$population, na.rm=T)/1e6))
  
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("Results saved to 'fertility_analysis_outputs' directory\n")
  
} else {
  stop("One or more required packages are not available. Please install them: 
       install.packages(c('tidycensus', 'tidyverse', 'scales', 'ggrepel', 'viridis', 'patchwork', 'knitr'))")
}