# Computational stress probe, NOT a national estimate. Reads the VT
# recovery-smoke's saved slim analytic_frames.rds (already on compute; no
# raw source reload, no new science) and replicates complete household
# blocks (prefixing household_key/person_key uniquely per replicate, so
# within-household clustering is preserved and no cross-replicate
# collisions occur) up to at least TARGET_N Twin1 rows -- the actual
# national Twin1_eligible_N observed in job 18247804 (857,607) -- to
# rehearse the production fit_instrument_outcome()/ar_confidence_set()
# code path at realistic scale on real compute, using the SAME
# controls_fml/grid/weight/cluster as the driver. This measures memory and
# time; it does not compute or claim a national effect estimate.
suppressMessages({ library(data.table); library(fixest); library(jsonlite) })
options(warn = 1)
data.table::setDTthreads(1)
options(fixest_nthreads = 1)

#' A fit result passes only if: RF/FS/IV coefficients and SEs are all
#' finite; AR ran without any error grid point (ar_n_errors is a finite
#' length-1 numeric equal to 0, and there is no ar_error field -- missing
#' ar_n_errors is a FAIL, not a pass, since that would let a swallowed AR
#' failure through); AR reports valid finite grid endpoints and a
#' non-missing component count. A bounded/interior AR confidence interval
#' is NOT required -- boundary-touching, disconnected, or empty accepted
#' sets are legitimate honest outcomes. actual_nobs (rf/fs/iv) must equal
#' the independently computed expected complete-case count for that
#' design/outcome (same Y/D/Z/controls/weight/cluster columns), or the
#' mismatch is reported explicitly rather than assumed benign.
check_fit_gate <- function(r, expected_n) {
  fin1 <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  coefs_ok <- fin1(r$rf_coef) && fin1(r$rf_se) && fin1(r$fs_coef) && fin1(r$fs_se) &&
    fin1(r$iv_coef) && fin1(r$iv_se)
  ar_n_errors_ok <- fin1(r$ar_n_errors) && identical(as.numeric(r$ar_n_errors), 0)
  no_ar_error_field <- is.null(r$ar_error)
  ar_grid_ok <- fin1(r$ar_grid_min) && fin1(r$ar_grid_max) && fin1(r$ar_n_components)
  nobs_ok <- fin1(r$rf_nobs_fit) && fin1(r$fs_nobs_fit) && fin1(r$iv_nobs_fit) &&
    r$rf_nobs_fit == expected_n && r$fs_nobs_fit == expected_n && r$iv_nobs_fit == expected_n
  list(
    pass = identical(r$status, "full_fit") && coefs_ok && ar_n_errors_ok &&
      no_ar_error_field && ar_grid_ok && nobs_ok,
    status_full_fit = identical(r$status, "full_fit"), coefs_finite = coefs_ok,
    ar_n_errors_ok = ar_n_errors_ok, no_ar_error_field = no_ar_error_field, ar_grid_ok = ar_grid_ok,
    expected_n = expected_n, actual_rf_nobs = r$rf_nobs_fit %||% NA_integer_,
    actual_fs_nobs = r$fs_nobs_fit %||% NA_integer_, actual_iv_nobs = r$iv_nobs_fit %||% NA_integer_,
    nobs_match = nobs_ok
  )
}

args_env <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = ""); if (nzchar(v)) v else default
}
root <- args_env("ROOT")
source(file.path(root, "twins_samesex_iv_lib.R"))
vt_outdir <- args_env("VT_OUTDIR")
outdir <- args_env("PROBE_OUTDIR")
target_n <- as.integer(args_env("TARGET_TWIN1_N", "857607"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

af_path <- file.path(vt_outdir, "analytic_frames.rds")
if (!file.exists(af_path)) stop(sprintf("VT analytic_frames.rds not found: %s", af_path))
af <- readRDS(af_path)
id_path <- file.path(vt_outdir, "analytic_frames_identity.json")
identity <- jsonlite::fromJSON(id_path)
controls_fml <- identity$controls_fml
weight_var <- identity$weight_var
cluster_var <- identity$cluster_var

t0 <- Sys.time()
n_vt_t1 <- nrow(af$t1)
n_vt_ss <- nrow(af$ss)
if (n_vt_t1 == 0) stop("VT t1 frame has zero rows; cannot scale")
reps <- ceiling(target_n / n_vt_t1)

replicate_block <- function(dt, reps) {
  out <- vector("list", reps)
  for (k in seq_len(reps)) {
    d <- data.table::copy(dt)
    pfx <- sprintf("rep%03d_", k)
    d[, person_key := paste0(pfx, person_key)]
    d[, household_key := paste0(pfx, household_key)]  # preserves within-HH clustering per replicate
    out[[k]] <- d
  }
  data.table::rbindlist(out, use.names = TRUE)
}

t1_scaled <- replicate_block(af$t1, reps)
ss_scaled <- replicate_block(af$ss, reps)
n_t1_generated <- nrow(t1_scaled)
n_ss_generated <- nrow(ss_scaled)
n_t1_unique_hh <- uniqueN(t1_scaled$household_key)
n_ss_unique_hh <- uniqueN(ss_scaled$household_key)
elapsed_replicate <- as.numeric(Sys.time() - t0, units = "secs")

cat(sprintf("PROBE: VT t1=%d ss=%d -> reps=%d -> scaled t1=%d (target %d) ss=%d, replicate_sec=%.1f\n",
            n_vt_t1, n_vt_ss, reps, n_t1_generated, target_n, n_ss_generated, elapsed_replicate))

outcomes <- c("ROOMS_out", "OWNERSHP_out", "BEDROOMS_out")
ar_grid_rooms <- seq(-3, 3, by = 0.1)
ar_grid_own <- seq(-0.5, 0.5, by = 0.02)

receipt_dir <- file.path(outdir, "probe_primary_receipts")
dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
make_cb <- function(design, oc, expected_n = NA_integer_) {
  force(design); force(oc); force(expected_n)
  function(res) {
    res$design <- design; res$outcome <- oc; res$expected_complete_case_n <- expected_n
    fn <- file.path(receipt_dir, sprintf("%s__%s.json", design, oc))
    tmp <- paste0(fn, ".tmp")
    jsonlite::write_json(res, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 10, null = "null", na = "null")
    if (!file.rename(tmp, fn)) stop(sprintf("probe: atomic rename failed for %s", fn))
  }
}

ctrl_vars <- all.vars(stats::as.formula(paste("~", controls_fml)))
complete_case_n <- function(dt, outcome, treatment, instrument) {
  cols <- unique(c(outcome, treatment, instrument, weight_var, cluster_var, ctrl_vars))
  sum(stats::complete.cases(dt[, ..cols]))
}

results <- list()
gate_checks <- list()
t_fit0 <- Sys.time()
for (oc in outcomes) {
  ar_grid <- if (oc == "OWNERSHP_out") ar_grid_own else ar_grid_rooms
  exp_n1 <- complete_case_n(t1_scaled, oc, "treatment_2plus", "twin_like_proxy")
  r1 <- tryCatch(fit_instrument_outcome(t1_scaled, oc, "treatment_2plus", "twin_like_proxy",
                   controls_fml, weight_var, cluster_var, ar_grid,
                   primary_callback = make_cb("Twin1_probe", oc, exp_n1)),
                 error = function(e) list(status = "error", message = conditionMessage(e)))
  r1$design <- "Twin1_probe"; r1$outcome <- oc
  results[[length(results) + 1]] <- r1
  gate_checks[[length(gate_checks) + 1]] <- c(list(design = "Twin1_probe", outcome = oc),
                                                check_fit_gate(r1, exp_n1))
  if (n_ss_generated > 0) {
    exp_n2 <- complete_case_n(ss_scaled, oc, "treatment_3plus", "samesex")
    r2 <- tryCatch(fit_instrument_outcome(ss_scaled, oc, "treatment_3plus", "samesex",
                     controls_fml, weight_var, cluster_var, ar_grid,
                     primary_callback = make_cb("SameSex2_probe", oc, exp_n2)),
                   error = function(e) list(status = "error", message = conditionMessage(e)))
    r2$design <- "SameSex2_probe"; r2$outcome <- oc
    results[[length(results) + 1]] <- r2
    gate_checks[[length(gate_checks) + 1]] <- c(list(design = "SameSex2_probe", outcome = oc),
                                                  check_fit_gate(r2, exp_n2))
  }
  cat(sprintf("PROBE: outcome %s done, elapsed_sec=%.1f\n", oc, as.numeric(Sys.time() - t_fit0, units = "secs")))
}
elapsed_fit <- as.numeric(Sys.time() - t_fit0, units = "secs")

per_case_pass <- vapply(gate_checks, function(g) isTRUE(g$pass), logical(1))
counts_match <- (n_t1_generated >= target_n)  # honest outcome-level missingness is NOT required
                                               # to hit target_n; this only checks the GENERATED
                                               # replicated sample reached the target scale.

gate <- list(
  n_full_fit = sum(per_case_pass), n_total = length(gate_checks),
  all_cases_pass = all(per_case_pass),
  n_t1_generated = n_t1_generated, n_t1_unique_households = n_t1_unique_hh,
  n_ss_generated = n_ss_generated, n_ss_unique_households = n_ss_unique_hh,
  target_twin1_n = target_n, counts_match_target = counts_match,
  elapsed_replicate_sec = elapsed_replicate, elapsed_fit_sec = elapsed_fit,
  per_case = gate_checks,
  gate_pass = all(per_case_pass) && counts_match
)
jsonlite::write_json(gate, file.path(outdir, "probe_gate_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
utils::write.csv(data.table::rbindlist(lapply(results, function(r) data.table(
  design = r$design, outcome = r$outcome, status = r$status,
  rf_coef = r$rf_coef %||% NA_real_, iv_coef = r$iv_coef %||% NA_real_,
  ar_n_errors = r$ar_n_errors %||% NA_integer_)), fill = TRUE),
  file.path(outdir, "probe_results_table.csv"), row.names = FALSE)

cat(sprintf("PROBE_GATE: %s (pass=%d/%d cases, n_t1=%d>=%d)\n",
            if (gate$gate_pass) "PASS" else "FAIL",
            gate$n_full_fit, gate$n_total, n_t1_generated, target_n))
if (!gate$gate_pass) quit(status = 1)
cat("PROBE_COMPLETE\n")
