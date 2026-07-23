#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
})

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", file_arg[startsWith(file_arg, file_flag)])
  if (length(script_path) == 0) {
    stop("Run this file with Rscript.")
  }
  dirname(normalizePath(script_path[1]))
}

ratio_for_sample <- function(df) {
  wealth_rows <- df %>%
    filter(
      is.finite(weight),
      weight > 0,
      is.finite(net_worth)
    )
  earnings_rows <- wealth_rows %>%
    filter(
      age <= 62,
      is.finite(earnings),
      earnings >= 0
    )
  aggregate_wealth <- sum(wealth_rows$weight * wealth_rows$net_worth)
  aggregate_gross_earnings <- sum(earnings_rows$weight * earnings_rows$earnings)
  tibble(
    wealth_households = nrow(wealth_rows),
    earnings_households = nrow(earnings_rows),
    aggregate_wealth = aggregate_wealth,
    aggregate_gross_earnings = aggregate_gross_earnings,
    gross_wealth_to_gross_earnings = aggregate_wealth / aggregate_gross_earnings
  )
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
psid_path <- file.path(dirname(repo_root), "PSID", "PSIDSHELF_MOBILITY.dta")
out_dir <- file.path(
  repo_root,
  "code",
  "data",
  "psid_followup_mar2026",
  "output",
  "aggregate_wealth_earnings_audit"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(psid_path)) {
  stop("Missing PSID shelf file: ", psid_path)
}

df <- read_dta(
  psid_path,
  col_select = c(
    "year", "RELTOHEAD_", "AGEREP", "NETWORTHR", "EARNINDRRC", "IW"
  )
) %>%
  transmute(
    year = as.integer(year),
    relation_to_head = as.numeric(RELTOHEAD_),
    age = as.numeric(AGEREP),
    net_worth = as.numeric(NETWORTHR),
    earnings = as.numeric(EARNINDRRC),
    weight = as.numeric(IW)
  ) %>%
  filter(
    relation_to_head == 10,
    age >= 18,
    age <= 82,
    year >= 1984,
    year <= 2019
  )

yearly <- df %>%
  group_by(year) %>%
  group_modify(~ ratio_for_sample(.x)) %>%
  ungroup() %>%
  mutate(
    after_payroll_equivalent_0179 = gross_wealth_to_gross_earnings / (1 - 0.179)
  )

vintages <- bind_rows(
  ratio_for_sample(filter(df, year >= 1984, year <= 2003)) %>%
    mutate(vintage = "1984-2003"),
  ratio_for_sample(filter(df, year >= 2005, year <= 2019)) %>%
    mutate(vintage = "2005-2019"),
  ratio_for_sample(df) %>%
    mutate(vintage = "1984-2019")
) %>%
  relocate(vintage) %>%
  mutate(
    after_payroll_equivalent_0179 = gross_wealth_to_gross_earnings / (1 - 0.179),
    effective_tax_rate_implied_by_6_9 =
      1 - gross_wealth_to_gross_earnings / 6.9
  )

write_csv(yearly, file.path(out_dir, "yearly_ratios.csv"))
write_csv(vintages, file.path(out_dir, "vintage_ratios.csv"))

old <- vintages %>% filter(vintage == "1984-2003")
recent <- vintages %>% filter(vintage == "2005-2019")
readme <- sprintf(
  paste0(
    "# PSID aggregate wealth/earnings audit\n\n",
    "This is a diagnostic, not an adopted calibration target. It uses one ",
    "reference-person row per family, ages 18--82, `NETWORTHR` in the ",
    "numerator, and RP/spouse combined `EARNINDRRC` for households whose ",
    "reference person is at most 62 in the denominator. All aggregates use ",
    "the PSID longitudinal weight `IW`.\n\n",
    "- 1984--2003 gross wealth/gross earnings: `%.6f`.\n",
    "- 1984--2003 equivalent after the model's 17.9%% payroll wedge: `%.6f`.\n",
    "- 2005--2019 gross wealth/gross earnings: `%.6f`.\n",
    "- 2005--2019 equivalent after the model's 17.9%% payroll wedge: `%.6f`.\n\n",
    "The exercise does not reproduce Hendricks's TAXSIM-based after-tax ",
    "earnings, adjust the PSID for top-wealth undercoverage, or settle the ",
    "appropriate data vintage. Its purpose is to expose the gross/after-tax ",
    "normalization and vintage sensitivity before a hard target is adopted.\n"
  ),
  old$gross_wealth_to_gross_earnings,
  old$after_payroll_equivalent_0179,
  recent$gross_wealth_to_gross_earnings,
  recent$after_payroll_equivalent_0179
)
writeLines(readme, file.path(out_dir, "README.md"))

print(vintages)
