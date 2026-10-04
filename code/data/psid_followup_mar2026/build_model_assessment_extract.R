#!/usr/bin/env Rscript

# One-pass, selected-column cache for downstream model-assessment plots.
# The PSID shelf is never modified.  A valid cache is reused without reopening it.
suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

SOURCE_PATH <- "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
EXPECTED_RATIO <- 6.92658379107299
AGE_MIN <- 18L
AGE_MAX <- 85L
WORKING_AGE_MAX <- 65L
MAX_SELECTED_BYTES <- 8 * 1024^3

script_path <- function() {
  x <- commandArgs(trailingOnly = FALSE)
  x <- sub("^--file=", "", x[startsWith(x, "--file=")])
  if (!length(x)) stop("Run with Rscript.")
  normalizePath(x[[1L]])
}
sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (!length(out) || !grepl("^[0-9a-fA-F]{64} ", out[[1L]])) stop("Unable to SHA-256: ", path)
  sub(" .*", "", out[[1L]])
}
file_identity <- function(path) {
  z <- file.info(path)
  if (is.na(z$size) || is.na(z$mtime)) stop("Cannot stat source: ", path)
  list(size_bytes = unname(z$size), mtime_utc = format(z$mtime, tz = "UTC", usetz = TRUE))
}
json_escape <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\\\"', x, fixed = TRUE)
  x <- gsub("\n", "\\\\n", x, fixed = TRUE)
  x
}
json_value <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return("null")
  if (is.logical(x)) return(if (x) "true" else "false")
  if (is.numeric(x)) return(format(x, digits = 17, scientific = FALSE, trim = TRUE))
  paste0('"', json_escape(x), '"')
}
as_number <- function(x) as.numeric(x) # tagged missing values remain missing; no imputation/recoding.

builder <- script_path()
script_dir <- dirname(builder)
out_dir <- file.path(script_dir, "output", "model_assessment")
cache_path <- file.path(out_dir, "psid_2005_2007.csv")
metadata_path <- file.path(out_dir, "metadata.json")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SOURCE_PATH)) stop("Missing PSID shelf: ", SOURCE_PATH)
source_before <- file_identity(SOURCE_PATH)
builder_sha <- sha256(builder)

# Regeneration guard: a cache is valid only for this exact source identity, builder,
# and selected-cache hash.  It avoids another 5.9GB shelf scan for plotting reruns.
if (file.exists(cache_path) && file.exists(metadata_path)) {
  old <- paste(readLines(metadata_path, warn = FALSE), collapse = "\n")
  cache_sha <- sha256(cache_path)
  same_source <- grepl(paste0('"size_bytes":', source_before$size_bytes), old, fixed = TRUE) &&
    grepl(paste0('"mtime_utc":"', source_before$mtime_utc, '"'), old, fixed = TRUE)
  same_builder <- grepl(paste0('"sourcebuilder_sha":"', builder_sha, '"'), old, fixed = TRUE)
  same_cache <- grepl(paste0('"selected_cache_sha256":"', cache_sha, '"'), old, fixed = TRUE)
  if (same_source && same_builder && same_cache) {
    message("Valid cache present; source was not reopened: ", cache_path)
    quit(status = 0L)
  }
  stop("Existing cache/metadata fails source identity or SHA validation; refusing a silent rescan. Remove or archive the owned cache explicitly before regeneration.")
}

started <- Sys.time()
message("progress: before raw read | ", format(started, tz = "UTC", usetz = TRUE))
selected <- c("ID", "year", "RELTOHEAD_", "AGEREP", "IW", "NETWORTHR", "EARNINDRRC",
              "HOMEOWN", "HOMEVALUER", "HOMEEQUITYR", "ACTUALROOMS_")

# This is the sole shelf read in this program. read_dta imports only selected columns.
raw <- as.data.table(read_dta(SOURCE_PATH, col_select = tidyselect::all_of(selected)))
if (!setequal(names(raw), selected)) {
  stop("STOPREPORT missing or unexpected selected fields: got ", paste(names(raw), collapse = ", "))
}
selected_bytes <- as.numeric(object.size(raw))
if (!is.finite(selected_bytes) || selected_bytes > MAX_SELECTED_BYTES) {
  stop("STOPREPORT selected-column memory estimate exceeds 8 GiB: ", selected_bytes)
}
message("progress: after raw read | rows=", nrow(raw), " selected_bytes=", selected_bytes)

labels <- vapply(selected, function(nm) {
  x <- attr(raw[[nm]], "label")
  if (is.null(x)) "" else as.character(x)
}, character(1))

dt <- raw[, .(
  id = as_number(ID), year = as_number(year), relation_to_head = as_number(RELTOHEAD_),
  age = as_number(AGEREP), weight = as_number(IW), total_net_wealth = as_number(NETWORTHR),
  annual_gross_labor_earnings = as_number(EARNINDRRC), homeown_code = as_number(HOMEOWN),
  home_value_raw = as_number(HOMEVALUER), home_equity = as_number(HOMEEQUITYR), rooms_raw = as_number(ACTUALROOMS_)
)]
rm(raw); invisible(gc())

retained <- dt[year %in% c(2005, 2007) & relation_to_head == 10 & age >= AGE_MIN & age <= AGE_MAX &
                 is.finite(weight) & weight > 0]
if (!nrow(retained)) stop("STOPREPORT retained universe is empty.")
if (retained[, anyDuplicated(paste(id, year, sep = "|"))] != 0L) {
  stop("STOPREPORT duplicate reference-person records by ID/year.")
}

# The authoritative broad-net-worth audit has a complete wealth numerator and a
# working-age earnings denominator. Retention intentionally does not impose either mask.
wealth_audit <- retained[is.finite(total_net_wealth) &
                           (age > WORKING_AGE_MAX | (is.finite(annual_gross_labor_earnings) & annual_gross_labor_earnings >= 0))]
numerator <- sum(wealth_audit$weight * wealth_audit$total_net_wealth)
working <- wealth_audit[age <= WORKING_AGE_MAX]
denominator <- sum(working$weight * working$annual_gross_labor_earnings)
ratio <- numerator / denominator
ratio_error <- ratio - EXPECTED_RATIO
if (!is.finite(ratio) || !isTRUE(all.equal(ratio, EXPECTED_RATIO, tolerance = 1e-10))) {
  stop(sprintf("STOPREPORT broad wealth/earnings identity mismatch: observed=%.17g expected=%.17g error=%.17g", ratio, EXPECTED_RATIO, ratio_error))
}

cache <- retained[, .(
  year = as.integer(year), age = age, weight = weight, total_net_wealth = total_net_wealth,
  annual_gross_labor_earnings = annual_gross_labor_earnings,
  owner = fifelse(homeown_code == 1, 1, fifelse(homeown_code == 2, 0, NA_real_)),
  gross_home_value = fifelse(homeown_code == 1 & is.finite(home_value_raw), home_value_raw,
                             fifelse(homeown_code == 2, 0, NA_real_)),
  home_equity = home_equity,
  rooms = fifelse(is.finite(rooms_raw) & rooms_raw > 0 & rooms_raw != 99, rooms_raw, NA_real_)
)]
fwrite(cache, cache_path, na = "")
cache_sha <- sha256(cache_path)
source_after <- file_identity(SOURCE_PATH)
if (!identical(source_before, source_after)) stop("STOPREPORT source changed during extraction.")
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
if (elapsed > 600) stop("STOPREPORT extraction exceeded 600 seconds: ", elapsed)

counts <- c(
  retained_rows = nrow(cache), retained_2005 = cache[year == 2005, .N], retained_2007 = cache[year == 2007, .N],
  finite_net_wealth = sum(is.finite(cache$total_net_wealth)), finite_earnings = sum(is.finite(cache$annual_gross_labor_earnings)),
  observed_owner = sum(is.finite(cache$owner)), observed_gross_home_value = sum(is.finite(cache$gross_home_value)),
  observed_home_equity = sum(is.finite(cache$home_equity)), observed_rooms = sum(is.finite(cache$rooms)),
  audit_rows = nrow(wealth_audit), audit_working_rows = nrow(working)
)
label_json <- paste(sprintf('"%s":"%s"', json_escape(selected), json_escape(labels)), collapse = ",")
count_json <- paste(sprintf('"%s":%s', names(counts), format(counts, scientific = FALSE, trim = TRUE)), collapse = ",")
metadata <- paste0(
  "{\n",
  '"source_path":"', json_escape(SOURCE_PATH), '",\n',
  '"source_before":{"size_bytes":', source_before$size_bytes, ',"mtime_utc":"', source_before$mtime_utc, '"},\n',
  '"source_after":{"size_bytes":', source_after$size_bytes, ',"mtime_utc":"', source_after$mtime_utc, '"},\n',
  '"sourcebuilder_sha":"', builder_sha, '",\n', '"selected_cache_sha256":"', cache_sha, '",\n',
  '"currency":"2022 USD (as stated in source variable labels)",\n', '"variable_labels":{', label_json, '},\n',
  '"filters":"year in {2005,2007}; RELTOHEAD_=10; AGEREP 18..85; finite IW>0; no outcome complete-case restriction",\n',
  '"field_definitions":"owner=1 if HOMEOWN=1, 0 if HOMEOWN=2, else NA; gross_home_value=HOMEVALUER for valid finite owners, 0 for renters, else NA; rooms=ACTUALROOMS_ only when finite, >0, and !=99; no recoding of 0/99; numeric weights unchanged; no winsorization",\n',
  '"broad_wealth_earnings_audit":{"definition":"sum(IW*NETWORTHR)/sum_{age<=65}(IW*EARNINDRRC), with finite NETWORTHR and (age>65 or finite EARNINDRRC>=0)","expected":', format(EXPECTED_RATIO, digits = 17), ',"observed":', format(ratio, digits = 17), ',"error":', format(ratio_error, digits = 17), '},\n',
  '"selected_column_memory_bytes":', format(selected_bytes, scientific = FALSE, trim = TRUE), ',"elapsed_seconds":', format(elapsed, digits = 17), ',\n',
  '"row_outcome_counts":{', count_json, '},\n',
  '"exact_command":"OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript ', json_escape(builder), '"\n', "}\n"
)
writeLines(metadata, metadata_path, useBytes = TRUE)
message("progress: complete | rows=", nrow(cache), " elapsed_seconds=", elapsed, " cache_sha256=", cache_sha)
