#!/usr/bin/env Rscript
# Recover a saved state panel after adapter-only metadata normalization.
# This intentionally skips raw-source reading, cleaning, matching, and housing.
options(stringsAsFactors = FALSE, scipen = 999)
root <- Sys.getenv("PROJECT_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
panel_file <- Sys.getenv("PANEL_FILE")
outdir <- Sys.getenv("OUTDIR")
if (!nzchar(panel_file) || !nzchar(outdir)) stop("PANEL_FILE and OUTDIR are required", call. = FALSE)
if (!file.exists(panel_file)) stop("saved panel absent: ", panel_file, call. = FALSE)
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  stop("OUTDIR exists and is non-empty; refusing overwrite", call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
matcher_file <- file.path(root, "code", "empirical", "acs", "kleven_pseudo", "run_national_acs_match_housing.R")
matcher_expr <- parse(file = matcher_file)
adapter_env <- new.env(parent = globalenv())
adapter_env$req <- function(ok, msg, stage = "recovery") if (!isTRUE(ok)) stop(msg, call. = FALSE)
for (e in matcher_expr) {
  is_assignment <- is.call(e) && identical(e[[1L]], as.name("<-"))
  if (is_assignment && (identical(e[[2L]], as.name("coalesce_field")) ||
                        identical(e[[2L]], as.name("normalize_lineage")) ||
                        identical(e[[2L]], as.name("narrow_estimator_panel")) ||
                        identical(e[[2L]], as.name("estimator_pool_columns"))))
    eval(e, envir = adapter_env)
}
if (!exists("normalize_lineage", envir = adapter_env, inherits = FALSE))
  stop("adapter normalize_lineage function absent", call. = FALSE)
estimator_file <- file.path(root, "code", "empirical", "acs", "kleven_pseudo", "estimate_national_first_birth_housing.R")
source(estimator_file, local = TRUE)
panel <- readRDS(panel_file)
panel <- adapter_env$narrow_estimator_panel(panel)
if (!all(panel$source_origin %in% c("ACS", "CPS"))) stop("normalized source origin invalid", call. = FALSE)
if (any(c("source_origin.x", "source_origin.y", "from_cps.x", "from_cps.y") %in% names(panel)))
  stop("stale adapter metadata remains after normalization", call. = FALSE)
checkpoint <- function(x) {
  jsonlite::write_json(x, file.path(outdir, "latest_checkpoint.json"), auto_unbox = TRUE, pretty = TRUE)
}
fit <- estimate_national_first_birth_housing(
  panel, output_dir = outdir, checkpoint = checkpoint,
  source_origin_col = "source_origin", from_cps_col = "from_cps",
  geography_label = "Vermont ACS")
if (!identical(fit$status, "ESTIMATION_COMPLETE_DIAGNOSTIC"))
  stop("estimator status invalid: ", fit$status, call. = FALSE)
jsonlite::write_json(list(status = fit$status, panel_file = panel_file,
                          output_dir = outdir, fit_count = length(fit$fits),
                          generated = format(Sys.time(), tz = "UTC")),
                     file.path(outdir, "estimator_recovery_receipt.json"),
                     auto_unbox = TRUE, pretty = TRUE)
cat("NATIONAL_ACS_ESTIMATOR_RECOVERY_PASS\n")
