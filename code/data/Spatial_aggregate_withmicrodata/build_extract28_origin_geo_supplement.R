#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
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
raw_dir <- file.path(script_dir, "raw_data")
ddi_path <- file.path(raw_dir, "usa_00028.xml")
out_path <- file.path(raw_dir, "extract28_origin_geo_2012_2023_age22_45_households.csv")
summary_path <- file.path(raw_dir, "extract28_origin_geo_2012_2023_age22_45_households_summary.txt")

if (!file.exists(ddi_path)) {
  stop("Missing raw IPUMS DDI: ", ddi_path)
}

if (file.exists(out_path)) {
  file.remove(out_path)
}

ddi <- read_ipums_ddi(ddi_path, lower_vars = TRUE)

rows_seen <- 0L
rows_kept <- 0L

callback <- IpumsSideEffectCallback$new(function(x, index) {
  dt <- as.data.table(x)

  dt[, `:=`(
    year = as.integer(year),
    sample = as.integer(sample),
    serial = as.numeric(serial),
    pernum = as.integer(pernum),
    age = as.integer(age),
    gq = as.integer(gq),
    migplac1 = as.integer(migplac1),
    migmetro1 = as.integer(migmetro1)
  )]

  rows_seen <<- rows_seen + nrow(dt)

  dt <- dt[
    year >= 2012 & year <= 2023 &
      age >= 22 & age <= 45 &
      gq %in% c(1L, 2L),
    .(year, sample, serial, pernum, migplac1, migmetro1)
  ]

  rows_kept <<- rows_kept + nrow(dt)

  fwrite(
    dt,
    out_path,
    append = file.exists(out_path),
    col.names = !file.exists(out_path)
  )

  message(
    sprintf(
      "Chunk %d processed: kept %s rows (%s seen so far).",
      index,
      format(nrow(dt), big.mark = ","),
      format(rows_seen, big.mark = ",")
    )
  )
})

read_ipums_micro_chunked(
  ddi,
  callback = callback,
  chunk_size = 200000,
  vars = c("year", "sample", "serial", "pernum", "age", "gq", "migplac1", "migmetro1"),
  verbose = TRUE
)

summary_lines <- c(
  sprintf("Source DDI: %s", ddi_path),
  sprintf("Output CSV: %s", out_path),
  sprintf("Rows seen: %s", format(rows_seen, big.mark = ",")),
  sprintf("Rows kept: %s", format(rows_kept, big.mark = ",")),
  "Filter: year 2012-2023, age 22-45, GQ in {1,2}",
  "Columns: year, sample, serial, pernum, migplac1, migmetro1"
)

writeLines(summary_lines, summary_path)
message("Saved supplement: ", out_path)
