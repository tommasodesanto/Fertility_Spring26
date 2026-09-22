#!/usr/bin/env Rscript
# Estimator-only continuation for the completed national ACS pool.
# This driver never rematches, rereads raw ACS/CPS inputs, or overwrites the
# production output.  It validates the exact saved pool/production/snapshot
# contract before loading the pooled panel and reuses only its own checkpoint
# directory.
options(stringsAsFactors = FALSE, scipen = 999)

nac_stop <- function(...) stop(paste0("national_acs_continuation: ",
                                      paste0(..., collapse = "")), call. = FALSE)
nac_require <- function(ok, ...) if (!isTRUE(ok)) nac_stop(...)
nac_read_json <- function(path) {
  nac_require(file.exists(path), "receipt missing: ", path)
  if (!requireNamespace("jsonlite", quietly = TRUE)) nac_stop("jsonlite is required")
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}
nac_write_json <- function(object, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) nac_stop("jsonlite is required")
  jsonlite::write_json(object, path, auto_unbox = TRUE, pretty = TRUE, digits = 17)
}
nac_under <- function(path, root, label) {
  p <- normalizePath(path, mustWork = FALSE)
  r <- normalizePath(root, mustWork = FALSE)
  nac_require(identical(p, r) || startsWith(paste0(p, "/"), paste0(r, "/")),
              label, " must be under root: ", p)
  p
}
nac_file_identity <- function(path, label = path) {
  z <- file.info(path)
  nac_require(nrow(z) == 1L && isTRUE(z$isdir == FALSE) && is.finite(z$size),
              label, " is not a regular file: ", path)
  list(path = normalizePath(path), bytes = as.numeric(z$size),
       mtime = format(z$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"))
}
nac_same_identity <- function(a, b) identical(a$bytes, b$bytes) && identical(a$mtime, b$mtime)
nac_states <- c(1L, 2L, 4L, 5L, 6L, 8L, 9L, 10L, 11L, 12L, 13L, 15L, 16L, 17L,
                18L, 19L, 20L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L,
                30L, 31L, 32L, 33L, 34L, 35L, 36L, 37L, 38L, 39L, 40L, 41L,
                42L, 44L, 45L, 46L, 47L, 48L, 49L, 50L, 51L, 53L, 54L, 55L, 56L)

nac_validate_inputs <- function(root, production_dir, snapshot_root, panel_file,
                                reuse_dir, expected_panel_bytes = 6354971545,
                                expected_pool_rows = 47310973) {
  nac_under(production_dir, root, "production_dir")
  nac_under(snapshot_root, root, "snapshot_root")
  nac_under(panel_file, root, "panel_file")
  expected_panel <- file.path(production_dir, "national_cps_acs_pseudo-panel_housing.rds")
  nac_require(identical(normalizePath(panel_file, mustWork = FALSE),
                         normalizePath(expected_panel, mustWork = FALSE)),
              "panel_file must be the exact production pooled-panel path")
  expected_reuse <- normalizePath(file.path(production_dir, "national_first_birth_housing"), mustWork = FALSE)
  nac_require(identical(normalizePath(reuse_dir, mustWork = FALSE), expected_reuse),
              "reuse_dir must be the completed national fit directory")
  nac_require(dir.exists(reuse_dir), "reuse checkpoint directory absent: ", reuse_dir)
  prod_receipt <- file.path(production_dir, "production_start_receipt.json")
  pool_receipt <- file.path(production_dir, "national_pool_manifest.json")
  stage_receipt <- file.path(production_dir, "stage_receipt.json")
  source_contract <- file.path(snapshot_root, "code/empirical/acs/kleven_pseudo/source_contract.json")
  for (p in c(prod_receipt, pool_receipt, source_contract, panel_file))
    nac_require(file.exists(p), "required contract file absent: ", p)
  prod <- nac_read_json(prod_receipt); pool <- nac_read_json(pool_receipt)
  source <- nac_read_json(source_contract)
  nac_require(identical(as.character(prod$status), "NATIONAL_PRODUCTION_START") &&
                identical(as.character(prod$phase), "production"),
              "production receipt is not the national production contract")
  nac_require(identical(as.integer(unlist(prod$states)), nac_states),
              "production receipt state set mismatch")
  nac_require(identical(as.character(pool$status), "NATIONAL_ESTIMATOR_POOL_COMPLETE") &&
                identical(as.integer(pool$pooled_rows), as.integer(expected_pool_rows)) &&
                identical(as.integer(unlist(pool$states)), nac_states),
              "pool manifest mismatch")
  state_receipts <- file.path(production_dir, sprintf("statefip_%02d/state_receipt.json", nac_states))
  nac_require(all(file.exists(state_receipts)), "one or more terminal state receipts are absent")
  state_status <- vapply(state_receipts, function(p) as.character(nac_read_json(p)$status), character(1))
  nac_require(all(state_status == "STATE_MATCH_COMPLETE"), "one or more state receipts are not terminal")
  stage <- if (file.exists(stage_receipt)) nac_read_json(stage_receipt) else NULL
  if (!is.null(stage)) {
    nac_require(identical(as.character(stage$status), "MATCH_HOUSING_STAGE_COMPLETE") &&
                  identical(as.character(stage$phase), "production"),
                "present stage receipt is not a completed production receipt")
  }
  nac_require(identical(as.character(source$remote_root), normalizePath(root, mustWork = FALSE)),
              "snapshot source contract remote root mismatch")
  prod_vendor <- sort(unlist(prod$vendor_sha256)); source_vendor <- sort(unlist(source$vendor_sha256))
  nac_require(length(prod_vendor) > 0L && identical(prod_vendor, source_vendor),
              "production vendor SHA-256 receipt differs from snapshot contract")
  panel_id <- nac_file_identity(panel_file, "pooled panel")
  nac_require(identical(as.numeric(panel_id$bytes), as.numeric(expected_panel_bytes)),
              "pooled panel byte identity mismatch: got ", panel_id$bytes,
              ", expected ", expected_panel_bytes)
  reuse_files <- file.path(reuse_dir, paste0("checkpoint_rooms9_",
                                              c("full", "event_only", "age_only", "state_year"), ".rds"))
  nac_require(all(file.exists(reuse_files)), "completed rooms checkpoint set is incomplete")
  receipt_paths <- c(production = prod_receipt, pool = pool_receipt, snapshot = source_contract,
                     setNames(state_receipts, paste0("state_", nac_states)))
  if (!is.null(stage)) receipt_paths <- c(receipt_paths, stage = stage_receipt)
  list(panel = panel_id, receipts = lapply(receipt_paths, nac_file_identity),
       receipt_paths = receipt_paths, production = prod, pool = pool, stage = stage)
}

nac_main <- function() {
  root <- Sys.getenv("PROJECT_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
  production_dir <- Sys.getenv("PRODUCTION_DIR", file.path(root, "output/national_acs_production_20260921a"))
  snapshot_root <- Sys.getenv("SNAPSHOT_ROOT", file.path(root, "snapshots/national_acs_f81814dd"))
  panel_file <- Sys.getenv("PANEL_FILE", file.path(production_dir, "national_cps_acs_pseudo-panel_housing.rds"))
  reuse_dir <- Sys.getenv("REUSE_DIR", file.path(production_dir, "national_first_birth_housing"))
  outdir <- Sys.getenv("OUTDIR", file.path(root, "output/national_acs_estimator_continuation_20260921a"))
  code_root <- Sys.getenv("CONTINUATION_CODE_ROOT", snapshot_root)
  nac_require(!dir.exists(outdir) || !length(list.files(outdir, all.files = TRUE, no.. = TRUE)),
              "OUTDIR exists and is non-empty; refusing overwrite: ", outdir)
  nac_under(outdir, root, "OUTDIR")
  contract <- nac_validate_inputs(root, production_dir, snapshot_root, panel_file, reuse_dir)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  start_panel <- contract$panel
  nac_write_json(list(status = "CONTINUATION_INPUT_CONTRACT_PASS", root = root,
                      production_dir = production_dir, snapshot_root = snapshot_root,
                      reuse_dir = reuse_dir, panel = start_panel,
                      receipts = contract$receipts, receipt_paths = contract$receipt_paths,
                      expected_panel_bytes = 6354971545, expected_pool_rows = 47310973,
                      generated = format(Sys.time(), tz = "UTC")),
                 file.path(outdir, "continuation_input_identity_start.json"))
  estimator_file <- file.path(code_root, "code/empirical/acs/kleven_pseudo/estimate_national_first_birth_housing.R")
  nac_require(file.exists(estimator_file), "continuation estimator absent: ", estimator_file)
  source(estimator_file, local = TRUE)
  panel <- readRDS(panel_file)
  nac_require(is.data.frame(panel) && identical(as.integer(nrow(panel)), 47310973L),
              "pooled panel row count changed")
  required <- c("matching_sample", "t_es_lw", "wgt", "gender", "statefip", "age_factor",
                "doiy_factor", "event_time", "rooms_raw", "ownershp_raw", "bedrooms_raw",
                "source_origin", "from_cps")
  nac_require(all(required %in% names(panel)), "pooled panel missing estimator columns: ",
              paste(setdiff(required, names(panel)), collapse = ","))
  latest <- function(x) nac_write_json(list(event = x, generated = format(Sys.time(), tz = "UTC")),
                                       file.path(outdir, "latest_checkpoint.json"))
  fit <- tryCatch(estimate_national_first_birth_housing(
    panel, output_dir = outdir, checkpoint = latest, reuse_dir = reuse_dir,
    outcomes = c("ownership_lw", "bedrooms5", "rooms9"),
    source_origin_col = "source_origin", from_cps_col = "from_cps",
    geography_label = "National ACS"), error = function(e) {
      nac_write_json(list(status = "CONTINUATION_FAILED", error = conditionMessage(e),
                          generated = format(Sys.time(), tz = "UTC")),
                     file.path(outdir, "continuation_failure.json"))
      stop(e)
    })
  end_panel <- nac_file_identity(panel_file, "pooled panel")
  nac_require(nac_same_identity(start_panel, end_panel), "pooled panel size/mtime changed during continuation")
  end_receipts <- lapply(contract$receipt_paths, nac_file_identity)
  nac_require(length(end_receipts) == length(contract$receipts) &&
                all(vapply(seq_along(end_receipts), function(i)
                  nac_same_identity(contract$receipts[[i]], end_receipts[[i]]), logical(1))),
              "production/snapshot receipt size or mtime changed during continuation")
  nac_require(identical(fit$status, "ESTIMATION_COMPLETE_DIAGNOSTIC") && length(fit$fits) == 12L,
              "continuation did not produce all 12 fits")
  nac_write_json(list(status = "NATIONAL_ACS_ESTIMATOR_CONTINUATION_COMPLETE",
                      fit_count = length(fit$fits), outcomes = unique(fit$curves$outcome),
                      specifications = unique(fit$curves$specification),
                      input_identity_start = start_panel, input_identity_end = end_panel,
                      receipt_identity_start = contract$receipts, receipt_identity_end = end_receipts,
                      reuse_dir = reuse_dir, output_dir = outdir,
                      generated = format(Sys.time(), tz = "UTC")),
                 file.path(outdir, "continuation_receipt.json"))
  cat("NATIONAL_ACS_ESTIMATOR_CONTINUATION_PASS\n")
}

if (!identical(Sys.getenv("NATIONAL_ACS_CONTINUATION_TEST"), "1")) nac_main()
