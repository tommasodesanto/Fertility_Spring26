#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_flag <- "--file="
script_path <- sub(file_flag, "", args[startsWith(args, file_flag)])
if (length(script_path) == 0) {
  stop("Run this script with Rscript so --file= is available.")
}

script_dir <- dirname(normalizePath(script_path[1]))
rscript <- file.path(R.home("bin"), "Rscript")
root_dir <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)

scripts <- c(
  file.path(root_dir, "code", "data", "Spatial_aggregate_withmicrodata", "build_extract28_origin_geo_supplement.R"),
  file.path(script_dir, "build_mms_geography.R"),
  file.path(script_dir, "compute_mms_descriptives.R"),
  file.path(script_dir, "build_migpuma_origin_bridge.R"),
  file.path(script_dir, "analyze_mms_fertility_moves.R")
)

for (script in scripts) {
  status <- system2(rscript, script)
  if (!identical(status, 0L)) {
    stop("Failed while running: ", script)
  }
}
