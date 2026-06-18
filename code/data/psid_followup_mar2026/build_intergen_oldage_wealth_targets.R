#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
})

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

weighted_quantile_safe <- function(x, w, prob = 0.5) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= prob)[1]]
}

post_moment <- function(rows, moment, value, sample, source, formula, status, n = NA_real_, weight = NA_real_) {
  rows[[length(rows) + 1]] <- tibble(
    moment = moment,
    value = as.numeric(value),
    sample = sample,
    n = as.numeric(n),
    weight = as.numeric(weight),
    source = source,
    formula = formula,
    status = status
  )
  rows
}

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", file_arg[startsWith(file_arg, file_flag)])
  if (length(script_path) == 0) {
    stop("Run this script with Rscript so --file= is available.")
  }
  dirname(normalizePath(script_path[1]))
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
psid_path <- file.path(dirname(repo_root), "PSID", "PSIDSHELF_MOBILITY.dta")
out_dir <- file.path(repo_root, "code", "data", "psid_followup_mar2026", "output", "intergen_oldage_wealth_targets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(psid_path)) {
  stop("Missing PSID shelf file: ", psid_path)
}

message("Reading PSID shelf for old-age wealth candidate targets...")
df <- read_dta(
  psid_path,
  col_select = c(
    "ID", "year", "AGEREP", "DEATHYEAR", "HOMEOWN", "RELCHIREP",
    "INCFAMR", "EARNINDR", "NETWORTHR", "NETWORTH2R", "NETWORTH3R",
    "HOMEEQUITYR", "WLTHSAVETOTR", "WLTHFUNDTOTR", "WLTHODEBTOTR", "IW"
  )
) %>%
  transmute(
    id = ID,
    year = as.integer(year),
    age = as.numeric(AGEREP),
    death_year = as.numeric(DEATHYEAR),
    own = case_when(HOMEOWN == 1 ~ 1, HOMEOWN == 2 ~ 0, TRUE ~ NA_real_),
    n_children = as.numeric(RELCHIREP),
    parent = n_children > 0,
    childless = n_children == 0,
    income = as.numeric(INCFAMR),
    earnings = as.numeric(EARNINDR),
    total_nw = as.numeric(NETWORTHR),
    nonhousing_nw = as.numeric(NETWORTH2R),
    broad_nw = as.numeric(NETWORTH3R),
    home_equity = as.numeric(HOMEEQUITYR),
    liquid_savings = as.numeric(WLTHSAVETOTR),
    liquid_funds = as.numeric(WLTHFUNDTOTR),
    other_debt = as.numeric(WLTHODEBTOTR),
    weight = as.numeric(IW)
  ) %>%
  filter(
    year >= 1984,
    year <= 2019,
    is.na(death_year) | year <= death_year,
    age >= 65,
    age <= 75,
    !is.na(weight),
    weight > 0,
    !is.na(n_children)
  ) %>%
  mutate(
    nonhousing_nw_to_income = if_else(income > 1000, nonhousing_nw / income, NA_real_),
    total_nw_to_income = if_else(income > 1000, total_nw / income, NA_real_),
    liquid_savings_to_income = if_else(income > 1000, liquid_savings / income, NA_real_)
  )

source <- "PSID PSIDSHELF_MOBILITY, 1984-2019, ages 65-75, completed children RELCHIREP"
status <- "candidate; re-audit wealth definitions/sample/formula before SMM use"
rows <- list()

for (sample_name in c("all_old_65_75", "parent_old_65_75", "childless_old_65_75")) {
  block <- switch(sample_name,
    all_old_65_75 = df,
    parent_old_65_75 = df %>% filter(parent),
    childless_old_65_75 = df %>% filter(childless)
  )
  sample_label <- switch(sample_name,
    all_old_65_75 = "age 65-75",
    parent_old_65_75 = "age 65-75, parent",
    childless_old_65_75 = "age 65-75, childless"
  )

  for (var in c("nonhousing_nw", "total_nw", "liquid_savings", "nonhousing_nw_to_income", "total_nw_to_income", "liquid_savings_to_income")) {
    rows <- post_moment(rows, paste0(sample_name, "_", var, "_mean"),
      weighted_mean_safe(block[[var]], block$weight),
      sample_label, source, paste0("weighted mean ", var), status,
      sum(!is.na(block[[var]])), sum(block$weight[!is.na(block[[var]])], na.rm = TRUE))
    rows <- post_moment(rows, paste0(sample_name, "_", var, "_median"),
      weighted_quantile_safe(block[[var]], block$weight, 0.5),
      sample_label, source, paste0("weighted median ", var), status,
      sum(!is.na(block[[var]])), sum(block$weight[!is.na(block[[var]])], na.rm = TRUE))
  }

  rows <- post_moment(rows, paste0(sample_name, "_own_rate"),
    weighted_mean_safe(block$own, block$weight),
    sample_label, source, "weighted mean homeownership", status,
    sum(!is.na(block$own)), sum(block$weight[!is.na(block$own)], na.rm = TRUE))
}

for (var in c("nonhousing_nw", "total_nw", "liquid_savings", "nonhousing_nw_to_income", "total_nw_to_income", "liquid_savings_to_income", "own")) {
  parent_block <- df %>% filter(parent)
  childless_block <- df %>% filter(childless)
  rows <- post_moment(rows, paste0("old_parent_minus_childless_", var, "_mean_gap"),
    weighted_mean_safe(parent_block[[var]], parent_block$weight) -
      weighted_mean_safe(childless_block[[var]], childless_block$weight),
    "age 65-75, parent minus childless", source,
    paste0("parent weighted mean ", var, " minus childless weighted mean ", var),
    status,
    sum(!is.na(df[[var]])), sum(df$weight[!is.na(df[[var]])], na.rm = TRUE))
}

targets <- bind_rows(rows)
write_csv(targets, file.path(out_dir, "intergen_oldage_wealth_targets.csv"))

by_child <- df %>%
  mutate(child_status = if_else(parent, "parent", "childless")) %>%
  group_by(child_status) %>%
  summarise(
    n = n(),
    total_weight = sum(weight, na.rm = TRUE),
    own_rate = weighted_mean_safe(own, weight),
    nonhousing_nw_mean = weighted_mean_safe(nonhousing_nw, weight),
    nonhousing_nw_median = weighted_quantile_safe(nonhousing_nw, weight, 0.5),
    nonhousing_nw_to_income_mean = weighted_mean_safe(nonhousing_nw_to_income, weight),
    nonhousing_nw_to_income_median = weighted_quantile_safe(nonhousing_nw_to_income, weight, 0.5),
    total_nw_mean = weighted_mean_safe(total_nw, weight),
    total_nw_median = weighted_quantile_safe(total_nw, weight, 0.5),
    .groups = "drop"
  )
write_csv(by_child, file.path(out_dir, "intergen_oldage_wealth_by_child_status.csv"))

readme <- c(
  "# Intergen Old-Age Wealth Candidate Targets",
  "",
  "These outputs are candidate replacement moments, not approved SMM targets.",
  "They must be reaudited for wealth definitions, completed-children timing, weights, and model-object match before use.",
  "",
  "- `intergen_oldage_wealth_targets.csv`: collapsed candidate old-age wealth and ownership moments.",
  "- `intergen_oldage_wealth_by_child_status.csv`: compact parent/childless old-age table.",
  "",
  paste0("Source: ", source),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote PSID old-age candidate targets to ", out_dir)
