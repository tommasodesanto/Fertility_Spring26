#!/usr/bin/env Rscript
# Contract and interrupted-run regression tests for the estimator-only driver.
old <- Sys.getenv("NATIONAL_ACS_CONTINUATION_TEST")
on.exit(Sys.setenv(NATIONAL_ACS_CONTINUATION_TEST = old), add = TRUE)
Sys.setenv(NATIONAL_ACS_CONTINUATION_TEST = "1")
script_arg <- commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1L]
here <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
driver <- file.path(here, "run_national_acs_estimator_continuation.R")
source(driver, local = TRUE)

expect_error <- function(expr, pattern) {
  hit <- FALSE
  tryCatch(force(expr), error = function(e) {
    hit <<- grepl(pattern, conditionMessage(e), fixed = TRUE)
  })
  if (!hit) stop("expected error containing: ", pattern, call. = FALSE)
}

root <- tempfile("national_acs_continuation_contract_")
dir.create(root, recursive = TRUE)
prod <- file.path(root, "output/national_acs_production_20260921a")
snap <- file.path(root, "snapshots/national_acs_f81814dd")
reuse <- file.path(prod, "national_first_birth_housing")
dir.create(reuse, recursive = TRUE)
dir.create(file.path(snap, "code/empirical/acs/kleven_pseudo"), recursive = TRUE)
panel <- file.path(prod, "national_cps_acs_pseudo-panel_housing.rds")
saveRDS(data.frame(x = 1), panel, compress = FALSE)
for (nm in c("production_start_receipt.json", "national_pool_manifest.json", "stage_receipt.json"))
  file.create(file.path(prod, nm))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
jsonlite::write_json(list(status = "NATIONAL_PRODUCTION_START", phase = "production",
                          states = as.list(nac_states)), file.path(prod, "production_start_receipt.json"), auto_unbox = TRUE)
jsonlite::write_json(list(status = "NATIONAL_ESTIMATOR_POOL_COMPLETE", pooled_rows = 0,
                          states = as.list(nac_states)), file.path(prod, "national_pool_manifest.json"), auto_unbox = TRUE)
jsonlite::write_json(list(status = "MATCH_HOUSING_STAGE_COMPLETE", phase = "production"),
                     file.path(prod, "stage_receipt.json"), auto_unbox = TRUE)
jsonlite::write_json(list(remote_root = normalizePath(root)),
                     file.path(snap, "code/empirical/acs/kleven_pseudo/source_contract.json"), auto_unbox = TRUE)
for (nm in c("full", "event_only", "age_only", "state_year"))
  file.create(file.path(reuse, paste0("checkpoint_rooms9_", nm, ".rds")))

ok <- nac_validate_inputs(root, prod, snap, panel, reuse,
                          expected_panel_bytes = file.info(panel)$size,
                          expected_pool_rows = 0)
stopifnot(identical(as.numeric(ok$panel$bytes), as.numeric(file.info(panel)$size)))

bad_panel <- file.path(prod, "different_panel.rds")
writeLines("different input", bad_panel)
expect_error(nac_validate_inputs(root, prod, snap, bad_panel, reuse,
                                 expected_panel_bytes = file.info(panel)$size,
                                 expected_pool_rows = 0), "byte identity mismatch")
expect_error(nac_validate_inputs(root, prod, snap, panel, file.path(prod, "other_fit_dir"),
                                 expected_panel_bytes = file.info(panel)$size,
                                 expected_pool_rows = 0), "completed national fit directory")

launcher <- readLines(file.path(here, "run_national_acs_estimator_continuation.sbatch"), warn = FALSE)
expect <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
expect(any(grepl("--cpus-per-task=4", launcher, fixed = TRUE)), "launcher CPU contract missing")
expect(any(grepl("--mem=128G", launcher, fixed = TRUE)), "launcher memory contract missing")
expect(any(grepl("--time=02:00:00", launcher, fixed = TRUE)), "launcher time contract missing")
expect(any(grepl("run_national_acs_estimator_continuation.R", launcher, fixed = TRUE)), "continuation driver missing")
expect(!any(grepl("run_national_acs_match_housing", launcher, fixed = TRUE)), "launcher invokes rematching")
cat("PASS: continuation receipt/path contract, mismatched-panel rejection, and 4CPU/128GB/2h estimator-only launcher\n")
