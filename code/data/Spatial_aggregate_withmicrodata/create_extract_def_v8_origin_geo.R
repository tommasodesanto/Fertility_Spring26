#!/usr/bin/env Rscript

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
def_dir <- file.path(script_dir, "extract_definitions")
dir.create(def_dir, recursive = TRUE, showWarnings = FALSE)

v7_path <- file.path(def_dir, "fertility_extract_def_v7_final_vars.rds")
if (!file.exists(v7_path)) {
  stop("Missing base extract definition: ", v7_path)
}

v7_def <- readRDS(v7_path)

# Keep the current extract27 payload, then add origin-geography variables needed
# for valid migration-PUMA interpretation and origin-side location classification.
extra_vars <- c(
  "CBSERIAL",
  "GQ",
  "OWNERSHPD",
  "RELATED",
  "RACED",
  "BPLD",
  "EMPSTATD",
  "MIGRATE1D",
  "VETSTATD",
  "MIGPLAC1",
  "MIGMETRO1"
)

make_var_spec <- function(name) {
  structure(list(name = name), class = c("var_spec", "ipums_spec", "list"))
}

v8_vars <- unique(c(names(v7_def$variables), extra_vars))
v8_var_specs <- lapply(v8_vars, make_var_spec)
names(v8_var_specs) <- v8_vars

v8_def <- v7_def
v8_def$description <- "Fertility Analysis V8 - V7 vars plus origin geography and extract27 alignment vars"
v8_def$variables <- v8_var_specs

out_path <- file.path(def_dir, "fertility_extract_def_v8_origin_geo.rds")
saveRDS(v8_def, out_path)

cat("Saved:\n")
cat(out_path, "\n\n")

cat("Added relative to v7:\n")
cat(paste(sort(setdiff(v8_vars, names(v7_def$variables))), collapse = "\n"), "\n")
