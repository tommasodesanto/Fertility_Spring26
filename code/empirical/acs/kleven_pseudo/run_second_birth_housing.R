#!/usr/bin/env Rscript

# Run the approved second-birth housing diagnostic from the completed exact
# support checkpoint.  This driver never reloads national data or rematches;
# it enriches the saved proxy input with source-key housing fields, then fits
# the three declared support specifications separately.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
script_path <- if (length(args)) args[[1L]] else "run_second_birth_housing.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
root <- Sys.getenv("KLEVEN_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
matched_dir <- Sys.getenv("SECOND_BIRTH_MATCHED_DIR", file.path(
  root, "output/kleven_acs_pilot/second_birth_transformed_support_18080547_corrected"))
proxy_file <- Sys.getenv("SECOND_BIRTH_PROXY_FILE", file.path(matched_dir, "proxy_checkpoint.rds"))
matched_file <- Sys.getenv("SECOND_BIRTH_MATCHED_FILE", file.path(matched_dir, "matched_checkpoint.rds"))
packet_file <- Sys.getenv("SECOND_BIRTH_PACKET_FILE", file.path(matched_dir, "packet_with_author_cells.rds"))
outdir <- Sys.getenv("SECOND_BIRTH_HOUSING_OUTDIR", file.path(
  root, "output/kleven_acs_pilot/second_birth_housing_20260920"))
estimator_file <- Sys.getenv("SECOND_BIRTH_HOUSING_FILE", file.path(script_dir, "second_birth_housing.R"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
progress_file <- file.path(outdir, "progress.log")
checkpoint <- function(stage, detail = "") {
  z <- sprintf("%s\t%s\t%s", format(Sys.time(), "%FT%T%z"), stage, detail)
  write(z, progress_file, append = TRUE)
  cat(z, "\n")
}
fail <- function(...) stop(paste0(...), call. = FALSE)
require_cols <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) fail(label, " missing: ", paste(miss, collapse = ", "))
}
key4 <- function(d) paste(d$YEAR, d$SAMPLE, d$SERIAL, d$PERNUM, sep = "\034")
resolve_raw <- function(d, candidates, label) {
  nms <- names(d)
  for (candidate in candidates) {
    hit <- nms[toupper(nms) == toupper(candidate)]
    if (length(hit) == 1L) return(hit[[1L]])
  }
  fail(label, " was not retained in the source packet; checked ", paste(candidates, collapse = ", "))
}

if (!file.exists(estimator_file)) fail("estimator file missing: ", estimator_file)
deps <- c(estimator_file, proxy_file, matched_file, packet_file)
if (any(!file.exists(deps))) fail("required input missing: ", paste(deps[!file.exists(deps)], collapse = ", "))
checkpoint("startup", paste0("matched_dir=", matched_dir))

source(estimator_file, local = TRUE)
required_fns <- c("prepare_second_birth_housing", "second_birth_housing_support", "fit_second_birth_housing")
if (!all(vapply(required_fns, exists, logical(1), envir = environment(), inherits = FALSE)))
  fail("incomplete second-birth housing estimator interface")
checkpoint("estimator_loaded")

proxy <- readRDS(proxy_file)
matched <- readRDS(matched_file)
packet <- as.data.table(readRDS(packet_file))
require_cols(proxy, c("input", "anchors", "post_rows"), "proxy checkpoint")
require_cols(matched, c("links", "target_support"), "matching checkpoint")
source_rows <- as.data.table(copy(proxy$input))
packet <- as.data.table(packet)
require_cols(source_rows, c("person_key", "YEAR", "SAMPLE", "SERIAL", "PERNUM"), "proxy input")
require_cols(packet, c("YEAR", "SAMPLE", "SERIAL", "PERNUM"), "source packet")
if (anyDuplicated(source_rows$person_key)) fail("proxy input person_key is not unique")
packet_key <- key4(packet)
if (anyDuplicated(packet_key)) fail("saved source packet key is not unique")
source_key <- key4(source_rows)
packet_idx <- match(source_key, packet_key)
if (anyNA(packet_idx)) fail("proxy input contains a source key absent from saved packet")

raw_map <- c(
  ROOMS = resolve_raw(packet, c("ROOMS_RAW", "ROOMS"), "ROOMS raw field"),
  BEDROOMS = resolve_raw(packet, c("BEDROOMS_RAW", "BEDROOMS"), "BEDROOMS raw field"),
  OWNERSHP = resolve_raw(packet, c("OWNERSHP_RAW", "OWNERSHP", "OWNERSHIP"), "OWNERSHP raw field")
)
for (nm in names(raw_map)) source_rows[[paste0("raw_", tolower(nm))]] <- packet[[raw_map[[nm]]]][packet_idx]

# These are the reviewed extract27 coding rules.  Unknown values stay missing;
# in particular ROOMS=28 and BEDROOMS=23 are retained in raw_* and flagged in
# the receipt rather than being recoded or dropped from the source packet.
rooms_raw <- suppressWarnings(as.numeric(source_rows$raw_rooms))
bedrooms_raw <- suppressWarnings(as.numeric(source_rows$raw_bedrooms))
own_raw <- suppressWarnings(as.numeric(source_rows$raw_ownershp))
source_rows[, rooms9 := ifelse(!is.na(rooms_raw) & rooms_raw %in% c(1:27, 30), pmin(rooms_raw, 9), NA_real_)]
source_rows[, bedrooms5 := ifelse(!is.na(bedrooms_raw) & bedrooms_raw %in% 1:22,
                                  pmin(ifelse(bedrooms_raw == 22, 21, bedrooms_raw - 1), 5), NA_real_)]
source_rows[, ownership_lw := ifelse(!is.na(own_raw) & own_raw %in% c(1, 2),
                                     ifelse(own_raw == 1, 1, 0), NA_real_)]
if (anyNA(source_rows$rooms9) && all(is.na(rooms_raw))) fail("ROOMS raw field was retained but entirely missing")
if (anyNA(source_rows$bedrooms5) && all(is.na(bedrooms_raw))) fail("BEDROOMS raw field was retained but entirely missing")
if (anyNA(source_rows$ownership_lw) && all(is.na(own_raw))) fail("OWNERSHP raw field was retained but entirely missing")
source_rows[, source_person_key := person_key]
proxy$input <- source_rows

coding_receipt <- data.table(
  variable = c("ROOMS", "BEDROOMS", "OWNERSHP"),
  raw_field = unname(raw_map),
  valid_rows = c(sum(!is.na(source_rows$rooms9)), sum(!is.na(source_rows$bedrooms5)), sum(!is.na(source_rows$ownership_lw))),
  missing_or_unknown_rows = c(sum(is.na(source_rows$rooms9)), sum(is.na(source_rows$bedrooms5)), sum(is.na(source_rows$ownership_lw))),
  rooms_unknown_28 = c(sum(rooms_raw == 28, na.rm = TRUE), NA_integer_, NA_integer_),
  bedrooms_unknown_23 = c(NA_integer_, sum(bedrooms_raw == 23, na.rm = TRUE), NA_integer_),
  ownership_unknown_3_9 = c(NA_integer_, NA_integer_, sum(own_raw %in% c(3, 9), na.rm = TRUE))
)
fwrite(coding_receipt, file.path(outdir, "housing_coding_receipt.csv"))
checkpoint("source_housing_enriched", sprintf("rows=%s rooms_raw=%s bedrooms_raw=%s ownership_raw=%s", nrow(source_rows), raw_map[[1L]], raw_map[[2L]], raw_map[[3L]]))

specs <- list(
  primary = list(support_spec = "all_anchors", fertyr_spec = "all", label = "all gap>=2 anchors"),
  joint_negative = list(support_spec = "joint_negative", fertyr_spec = "all", label = "anchors with donors at both -2 and -1"),
  fertyr_event0_yes = list(support_spec = "all_anchors", fertyr_spec = "event0_yes", label = "gap>=2 anchors with observed FERTYR yes at event 0")
)
outcomes <- c("rooms9", "bedrooms5", "ownership_lw")
all_status <- list()
for (spec_name in names(specs)) {
  spec <- specs[[spec_name]]
  spec_dir <- file.path(outdir, spec_name)
  fit_dir <- file.path(spec_dir, "fits")
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint("prepare_start", spec_name)
  prepared <- tryCatch(
    prepare_second_birth_housing(proxy, matched, source_rows = source_rows,
      outcomes = outcomes, support_spec = spec$support_spec, fertyr_spec = spec$fertyr_spec),
    error = function(e) e)
  if (inherits(prepared, "error")) {
    status <- data.table(outcome = outcomes, status = "PREPARE_FAILED", nobs = NA_integer_, error = conditionMessage(prepared))
    fwrite(status, file.path(spec_dir, "fit_status.csv"))
    all_status[[spec_name]] <- status
    checkpoint("prepare_failed", paste(spec_name, conditionMessage(prepared)))
    next
  }
  support <- second_birth_housing_support(prepared, outcomes = outcomes)
  fwrite(support, file.path(spec_dir, "support.csv"))
  write_json(list(specification = spec_name, label = spec$label, support = prepared$support,
                  contract = prepared$contract, generated = as.character(Sys.time())),
             file.path(spec_dir, "prepared_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
  checkpoint("prepare_complete", sprintf("%s rows=%s anchors=%s joint=%s", spec_name, nrow(prepared$panel), prepared$support$n_all_anchors, prepared$support$n_joint_negative))

  statuses <- list(); curves <- list(); contrasts <- list()
  for (y in outcomes) {
    checkpoint("fit_start", paste(spec_name, y))
    fit_result <- tryCatch(
      fit_second_birth_housing(prepared, outcomes = y, save_dir = fit_dir),
      error = function(e) e)
    if (inherits(fit_result, "error")) {
      statuses[[y]] <- data.table(outcome = y, status = "FIT_FAILED", nobs = NA_integer_, error = conditionMessage(fit_result))
      checkpoint("fit_failed", paste(spec_name, y, conditionMessage(fit_result)))
      next
    }
    rec <- fit_result$fits[[y]]
    statuses[[y]] <- data.table(outcome = y, status = "FIT_COMPLETE", nobs = rec$nobs, error = NA_character_)
    curves[[y]] <- fit_result$curves
    one_contrast <- as.data.table(fit_result$contrasts[[y]])
    one_contrast[, outcome := y]
    contrasts[[y]] <- one_contrast
    saveRDS(fit_result, file.path(spec_dir, paste0("fit_result_", y, ".rds")), compress = FALSE)
    checkpoint("fit_complete", paste(spec_name, y))
  }
  status <- rbindlist(statuses, fill = TRUE)
  fwrite(status, file.path(spec_dir, "fit_status.csv"))
  if (length(curves)) fwrite(rbindlist(curves, fill = TRUE), file.path(spec_dir, "curves.csv"))
  if (length(contrasts)) fwrite(rbindlist(contrasts, fill = TRUE), file.path(spec_dir, "contrasts.csv"))
  all_status[[spec_name]] <- status
}

status_all <- rbindlist(all_status, idcol = "specification", fill = TRUE)
fwrite(status_all, file.path(outdir, "fit_status_all.csv"))
write_json(list(
  status = if (all(status_all$status == "FIT_COMPLETE")) "COMPLETE" else "COMPLETED_WITH_FIT_FAILURES",
  source_proxy_dir = matched_dir, proxy_file = proxy_file, matched_file = matched_file,
  packet_file = packet_file, estimator_file = estimator_file,
  specifications = lapply(specs, function(x) x[c("support_spec", "fertyr_spec", "label")]),
  outcomes = outcomes, fit_status = status_all, causal_claim = FALSE,
  cluster = "YEAR:SAMPLE:SERIAL", generated = as.character(Sys.time())),
  file.path(outdir, "housing_run_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
checkpoint("complete", paste0("status=", if (all(status_all$status == "FIT_COMPLETE")) "COMPLETE" else "COMPLETED_WITH_FIT_FAILURES"))
