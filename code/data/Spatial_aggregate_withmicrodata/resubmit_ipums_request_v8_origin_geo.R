#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ipumsr)
})

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
def_path <- file.path(script_dir, "extract_definitions", "fertility_extract_def_v8_origin_geo.rds")
download_dir <- file.path(script_dir, "raw_data")
legacy_key_script <- "/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/fertility_codes/resubmit_ipums_request.R"

if (!file.exists(def_path)) {
  stop("Missing v8 extract definition. Run create_extract_def_v8_origin_geo.R first.")
}

api_key <- Sys.getenv("IPUMS_API_KEY", unset = "")
if (api_key == "") {
  if (file.exists(legacy_key_script)) {
    key_line <- grep("set_ipums_api_key\\(", readLines(legacy_key_script, warn = FALSE), value = TRUE)
    if (length(key_line) > 0) {
      api_key <- sub('.*set_ipums_api_key\\("([^"]+)".*', "\\1", key_line[1])
    }
  }
}

if (api_key == "") {
  stop("Set IPUMS_API_KEY in the environment or keep the legacy resubmit script available.")
}

set_ipums_api_key(api_key)

extract_def <- readRDS(def_path)
submitted_extract <- submit_extract(extract_def)

cat("Submitted extract #", submitted_extract$number, "\n", sep = "")
cat("Waiting for IPUMS to prepare the extract...\n")

waited_extract <- wait_for_extract(submitted_extract, timeout = 3600)

dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
download_extract(waited_extract, download_dir = download_dir, overwrite = TRUE)

cat("Download complete in:\n")
cat(download_dir, "\n")
