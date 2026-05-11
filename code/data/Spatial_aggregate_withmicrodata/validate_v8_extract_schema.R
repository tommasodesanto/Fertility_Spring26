#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript validate_v8_extract_schema.R /path/to/new_extract.dta")
}

new_extract_path <- normalizePath(args[1], mustWork = TRUE)
baseline_path <- normalizePath("code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta", mustWork = TRUE)
v8_def_path <- normalizePath("code/data/Spatial_aggregate_withmicrodata/extract_definitions/fertility_extract_def_v8_origin_geo.rds", mustWork = TRUE)

baseline_vars <- toupper(names(read_dta(baseline_path, n_max = 1)))
new_vars <- toupper(names(read_dta(new_extract_path, n_max = 1)))
v8_vars <- toupper(names(readRDS(v8_def_path)$variables))

cat("New extract path:\n")
cat(new_extract_path, "\n\n")

cat("Missing relative to extract27 baseline:\n")
print(setdiff(baseline_vars, new_vars))
cat("\n")

cat("New variables relative to extract27 baseline:\n")
print(setdiff(new_vars, baseline_vars))
cat("\n")

expected_new <- c("MIGPLAC1", "MIGMETRO1")
cat("Expected new variables relative to extract27 baseline:\n")
print(expected_new)
cat("\n")

cat("Missing relative to v8 definition:\n")
print(setdiff(v8_vars, new_vars))
cat("\n")

cat("Unexpected relative to v8 definition:\n")
print(setdiff(new_vars, v8_vars))
cat("\n")
