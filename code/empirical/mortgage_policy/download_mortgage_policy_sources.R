#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", args[startsWith(args, file_flag)])
  if (length(script_path) == 0) {
    stop("Run this script with Rscript so --file= is available.")
  }
  dirname(normalizePath(script_path[1]))
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "1" else "0")
  value %in% c("1", "true", "TRUE", "yes", "YES")
}

script_dir <- get_script_dir()
manifest_path <- file.path(script_dir, "source_manifest.csv")
raw_dir <- file.path(script_dir, "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

manifest <- fread(manifest_path, na.strings = c("", "NA"))

download_all <- env_flag("DOWNLOAD_ALL")
download_hmda <- env_flag("DOWNLOAD_HMDA")
download_fhfa <- env_flag("DOWNLOAD_FHFA")
force <- env_flag("FORCE")
year_text <- Sys.getenv("YEARS", unset = "")

todo <- manifest[
  download_all |
    download_default == TRUE |
    (download_hmda & source_group %in% c("hmda_lar", "hmda_dictionary")) |
    (download_fhfa & source_group == "fhfa_cll")
]

if (nzchar(year_text)) {
  year_keep <- trimws(strsplit(year_text, ",")[[1]])
  todo <- todo[
    !source_group %in% c("hmda_lar", "fhfa_cll") |
      as.character(year) %in% year_keep
  ]
}

if (!nrow(todo)) {
  stop("No sources selected. Set DOWNLOAD_ALL=1, DOWNLOAD_HMDA=1, or DOWNLOAD_FHFA=1.")
}

for (ii in seq_len(nrow(todo))) {
  row <- todo[ii]
  dest <- file.path(script_dir, row$raw_file)
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(dest) && !force) {
    message(sprintf("Skipping existing file: %s", dest))
    next
  }
  message(sprintf("Downloading %s %s %s", row$source_group, row$product, row$year))
  message(sprintf("  %s", row$url))
  ok <- tryCatch({
    download.file(row$url, destfile = dest, mode = "wb", quiet = FALSE)
    TRUE
  }, error = function(e) {
    message(sprintf("FAILED: %s", conditionMessage(e)))
    FALSE
  })
  if (!ok) {
    stop(sprintf("Download failed for %s", row$url))
  }
}

message("Done.")
