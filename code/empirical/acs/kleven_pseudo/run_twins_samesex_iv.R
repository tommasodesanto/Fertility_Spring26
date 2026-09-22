# National ACS Twin1 / SameSex2 housing IV driver. Reads state RDS
# partitions from the national_acs_source_stage_20260921 housing_narrow
# packet, builds the mother roster, constructs both instruments at each
# observed event age 0:5, runs RF/FS/2SLS/AR for ROOMS/OWNERSHP/BEDROOMS,
# and writes checkpoints + a compact status heartbeat.
#
# Env vars:
#   ROOT            repo root for source() (kleven_pseudo dir)
#   PARTITION_DIR    dir containing partitions/statefip_XX/housing_narrow.rds
#   STATEFIP_LIST    comma-separated state codes to load
#   OUTDIR           output directory (created if absent)
#   STATUS_PATH      path to status JSON heartbeat (optional)
#   SAMPLE_LABEL     "young" (21-35) or "wide" (25-45); default young

suppressMessages({
  library(data.table)
  library(fixest)
  library(jsonlite)
})

args_env <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

script_dir <- args_env("ROOT", dirname(sys.frame(1)$ofile))
source(file.path(script_dir, "twins_samesex_iv_lib.R"))

partition_dir <- args_env("PARTITION_DIR",
  "/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_source_stage_20260921/partitions")
statefip_list <- as.integer(strsplit(args_env("STATEFIP_LIST", "50"), ",")[[1]])
outdir <- args_env("OUTDIR", "output_acs_twins_samesex_smoke")
status_path <- args_env("STATUS_PATH", file.path(outdir, "status_heartbeat.json"))
sample_label <- args_env("SAMPLE_LABEL", "young")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Cleanup: a prior interrupted run in this OUTDIR can leave partial .tmp
# files from the atomic tmp+rename writes; remove them so stale partial
# state is never mistaken for a completed receipt.
unlink(list.files(outdir, pattern = "\\.tmp$", full.names = TRUE, recursive = TRUE))
last_heartbeat_time <- Sys.time()
heartbeat_interval_sec <- 300

write_status <- function(phase, extra = list()) {
  st <- c(list(phase = phase, generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
               statefip_list = statefip_list, sample_label = sample_label,
               outdir = normalizePath(outdir, mustWork = FALSE)), extra)
  tmp <- paste0(status_path, ".tmp")
  jsonlite::write_json(st, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 6)
  file.rename(tmp, status_path)
}

write_status("loading_partitions")

# Process one state partition at a time and keep only the much smaller
# per-mother roster in memory across states; never hold all 51 raw
# person-level partitions (59M rows nationally) simultaneously. This keeps
# the job inside the authorized 64GB envelope.
t0 <- Sys.time()
mr_list <- vector("list", length(statefip_list))
audit_list <- vector("list", length(statefip_list))
relate_audit_list <- vector("list", length(statefip_list))
total_rows <- 0L
sample_gate_totals <- list(n_input = 0L, n_excluded_year_out_of_range = 0L,
                            n_excluded_non_acs1yr_product = 0L, n_kept = 0L)
for (i in seq_along(statefip_list)) {
  sf <- statefip_list[i]
  p <- sprintf("%s/statefip_%02d/housing_narrow.rds", partition_dir, sf)
  if (!file.exists(p)) stop(sprintf("missing partition for statefip %d: %s", sf, p))
  dt_state_raw <- readRDS(p)
  total_rows <- total_rows + nrow(dt_state_raw)
  # Enforce the common national ACS 1-year product, 2005-2019, BEFORE any
  # roster fitting (SAMPLE == YEAR*100+1; source_audit_extract27_20260919.md
  # documented codes). Explicit exclusion counts accumulate across states.
  gate <- apply_source_sample_gate(dt_state_raw)
  dt_state <- gate$data
  sample_gate_totals$n_input <- sample_gate_totals$n_input + gate$counts$n_input
  sample_gate_totals$n_excluded_year_out_of_range <-
    sample_gate_totals$n_excluded_year_out_of_range + gate$counts$n_excluded_year_out_of_range
  sample_gate_totals$n_excluded_non_acs1yr_product <-
    sample_gate_totals$n_excluded_non_acs1yr_product + gate$counts$n_excluded_non_acs1yr_product
  sample_gate_totals$n_kept <- sample_gate_totals$n_kept + gate$counts$n_kept
  rm(dt_state_raw)
  built_state <- build_mother_roster(dt_state)
  mr_list[[i]] <- add_outcomes(built_state$mother_rows)
  audit_list[[i]] <- built_state$link_audit
  relate_audit_list[[i]] <- built_state$child_relate_audit
  rm(dt_state, built_state); gc(FALSE)
  write_status("loading_partitions", list(loaded = i, of = length(statefip_list),
                                           statefip = sf, rows_seen = total_rows,
                                           sample_gate_totals = sample_gate_totals,
                                           elapsed_sec = as.numeric(Sys.time() - t0, units = "secs")))
}
mr_pre_minor_gate <- data.table::rbindlist(mr_list, use.names = TRUE, fill = TRUE)
built <- list(
  link_audit = data.table::rbindlist(audit_list, use.names = TRUE, fill = TRUE)[
    , .(N = sum(N)), by = .(link_invalid_reason, child_relate)],
  child_relate_audit = data.table::rbindlist(relate_audit_list, use.names = TRUE, fill = TRUE)[
    , .(N = sum(N)), by = child_relate]
)
rm(mr_list, audit_list, relate_audit_list); gc(FALSE)
write_status("loaded", list(rows = total_rows, mother_rows = nrow(mr_pre_minor_gate),
                             sample_gate_totals = sample_gate_totals,
                             elapsed_sec = as.numeric(Sys.time() - t0, units = "secs")))

## ---- Sample-scope gates applied before roster fitting ---------------------
# Oldest-linked-child minor (<18) gate, applied before Twin1/SameSex2
# construction, with explicit exclusion counts (not a silent filter).
minor_gate <- apply_oldest_child_minor_gate(mr_pre_minor_gate)
mr <- minor_gate$data
rm(mr_pre_minor_gate); gc(FALSE)

mr[, weight_positive := PERWT > 0 & is.finite(PERWT)]
mr <- mr[weight_positive == TRUE]

age_lo <- if (sample_label == "young") 21 else 25
age_hi <- if (sample_label == "young") 35 else 45
mr_sample <- mr[AGE_norm >= age_lo & AGE_norm <= age_hi]

n_unique_mothers <- uniqueN(mr_sample$person_key)
n_households <- uniqueN(mr_sample$household_key)

write_status("roster_built", list(
  n_unique_mothers = n_unique_mothers, n_households = n_households,
  any_sex_missing = sum(mr_sample$any_sex_missing),
  mother_is_householder_share = mean(mr_sample$mother_is_householder, na.rm = TRUE),
  nchild_link_mismatch = sum(mr_sample$nchild_link_mismatch),
  sample_gate_totals = sample_gate_totals,
  n_excluded_oldest_child_adult = minor_gate$n_excluded_oldest_child_adult,
  elapsed_sec = as.numeric(Sys.time() - t0, units = "secs")))

jsonlite::write_json(built$link_audit, file.path(outdir, "link_audit.json"), auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(built$child_relate_audit, file.path(outdir, "child_relate_audit.json"),
                      auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(sample_gate_totals, file.path(outdir, "sample_gate_receipt.json"),
                      auto_unbox = TRUE, pretty = TRUE)

## ---- Instrument construction ---------------------------------------------
t1 <- build_twin1(mr_sample, age_grid = 0:5)
ss <- build_samesex2(mr_sample, age_grid = 0:5)

t1[, mother_weight := PERWT]
t1[, household_key := household_key]
ss[, mother_weight := PERWT]

counts <- list(
  sample_label = sample_label, age_lo = age_lo, age_hi = age_hi,
  n_unique_mothers_sample = n_unique_mothers,
  n_households_sample = n_households,
  Twin1_eligible_N = nrow(t1),
  Twin1_proxy_positive_N = sum(t1$twin_like_proxy),
  Twin1_contamination_risk_N = sum(t1$contamination_risk),
  Twin1_treatment_2plus_N = sum(t1$treatment_2plus),
  SameSex2_eligible_pool_N = nrow(ss),
  SameSex2_primary_age_tie_excluded_N = sum(ss$primary_age_tie),
  SameSex2_sex_missing_excluded_N = sum(!ss$primary_age_tie & (is.na(ss$sex1) | is.na(ss$sex2))),
  SameSex2_primary_eligible_N = sum(ss$eligible),
  SameSex2_positive_samesex_N = sum(ss$eligible & ss$samesex == 1),
  SameSex2_both_boys_N = sum(ss$both_boys, na.rm = TRUE),
  SameSex2_both_girls_N = sum(ss$both_girls, na.rm = TRUE),
  SameSex2_treatment_3plus_N = sum(ss$eligible & ss$treatment_3plus)
)
jsonlite::write_json(counts, file.path(outdir, "counts.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
write_status("instruments_built", counts)

## ---- Controls -------------------------------------------------------------
build_controls <- function(d, event_age_col) {
  d[, mat_age_at_event := round(AGE_norm - get(event_age_col))]
  d[, survey_year := YEAR]
  d
}
t1 <- build_controls(t1, "event_age")
ss <- build_controls(ss, "event_age")

controls_fml <- "i(mat_age_at_event) + i(RACE) + i(survey_year) + i(event_age)"

outcomes <- c("ROOMS_out", "OWNERSHP_out", "BEDROOMS_out")

## ---- Save slim analytic frames + identity receipt BEFORE any estimation --
## Persisted for a future reviewed estimator-only recovery loader (not
## implemented yet -- no automatic resume this pass). Counts/audits are
## already saved above, so SameSex2 keeps only the eligible==TRUE rows.
ctrl_vars <- all.vars(stats::as.formula(paste("~", controls_fml)))
t1_needed <- unique(c("person_key", "household_key", "mother_weight",
                       "twin_like_proxy", "treatment_2plus", "event_age", outcomes, ctrl_vars))
ss_needed <- unique(c("person_key", "household_key", "mother_weight",
                       "samesex", "treatment_3plus", "event_age", outcomes, ctrl_vars))
t1_slim <- t1[, ..t1_needed]
ss_slim <- ss[eligible == TRUE, ..ss_needed]
analytic_frames_path <- file.path(outdir, "analytic_frames.rds")
tmp_af <- paste0(analytic_frames_path, ".tmp")
saveRDS(list(t1 = t1_slim, ss = ss_slim), tmp_af)
file.rename(tmp_af, analytic_frames_path)
driver_file <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  normalizePath(sub("^--file=", "", a[grepl("^--file=", a)]), mustWork = FALSE)
}, error = function(e) NA_character_)
identity <- list(
  analytic_frames_path = analytic_frames_path,
  analytic_frames_size_bytes = file.info(analytic_frames_path)$size,
  analytic_frames_md5 = unname(tools::md5sum(analytic_frames_path)),
  driver_script = driver_file,
  driver_md5 = if (!is.na(driver_file)) unname(tools::md5sum(driver_file)) else NA_character_,
  lib_path = file.path(script_dir, "twins_samesex_iv_lib.R"),
  lib_md5 = unname(tools::md5sum(file.path(script_dir, "twins_samesex_iv_lib.R"))),
  sample_label = sample_label, age_lo = age_lo, age_hi = age_hi,
  controls_fml = controls_fml, outcomes = outcomes,
  weight_var = "mother_weight", cluster_var = "household_key",
  ar_grid_rooms = list(min = -3, max = 3, by = 0.1),
  ar_grid_own = list(min = -0.5, max = 0.5, by = 0.02),
  partition_dir = partition_dir, statefip_list = statefip_list,
  n_t1_rows = nrow(t1_slim), n_ss_rows = nrow(ss_slim),
  automatic_resume_implemented = FALSE,
  note = "Identity recorded for a future REVIEWED estimator-only loader; no automatic resume is implemented this pass.",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
jsonlite::write_json(identity, file.path(outdir, "analytic_frames_identity.json"), auto_unbox = TRUE, pretty = TRUE)
write_status("analytic_frames_saved", list(n_t1_rows = nrow(t1_slim), n_ss_rows = nrow(ss_slim)))
rm(mr, mr_sample, minor_gate); gc(FALSE)
ar_grid_rooms <- seq(-3, 3, by = 0.1)
ar_grid_own <- seq(-0.5, 0.5, by = 0.02)

receipt_dir <- file.path(outdir, "fit_receipts")
dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)

# Small atomic per-fit checkpoint: written immediately after each individual
# design/outcome fit completes, independent of the outcome-level loop, so
# completed RF/FS/IV work survives an interruption mid-outcome. Also serves
# as the driver's progress heartbeat (well under a 5-minute cadence for 18
# fits nationally).
checkpoint_fit <- function(res) {
  results[[length(results) + 1]] <<- res
  fname <- file.path(receipt_dir, sprintf("%s__%s.json", res$design %||% "NA", res$outcome %||% "NA"))
  tmp <- paste0(fname, ".tmp")
  jsonlite::write_json(res, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 10, null = "null", na = "null")
  file.rename(tmp, fname)
  saveRDS(results, paste0(file.path(outdir, "checkpoint_results.rds"), ".tmp"))
  file.rename(paste0(file.path(outdir, "checkpoint_results.rds"), ".tmp"),
              file.path(outdir, "checkpoint_results.rds"))
  write_status("estimating", list(last_completed = fname, n_fits_done = length(results),
                                   elapsed_sec = as.numeric(Sys.time() - t0, units = "secs")))
}

primary_receipt_dir <- file.path(outdir, "primary_receipts")
dir.create(primary_receipt_dir, recursive = TRUE, showWarnings = FALSE)
# Writes RF/FS/IV (no AR yet) atomically, invoked from inside
# fit_instrument_outcome BEFORE AR runs -- so a slow/failing AR cannot
# erase already-computed primary effects.
make_primary_cb <- function(design, oc) {
  force(design); force(oc)
  function(res) {
    res$design <- design; res$outcome <- oc
    fn <- file.path(primary_receipt_dir, sprintf("%s__%s.json", design, oc))
    tmp <- paste0(fn, ".tmp")
    jsonlite::write_json(res, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 10, null = "null", na = "null")
    file.rename(tmp, fn)
  }
}

results <- list()
for (oc in outcomes) {
  ar_grid <- if (oc == "OWNERSHP_out") ar_grid_own else ar_grid_rooms
  res_t1 <- tryCatch(fit_instrument_outcome(
    t1, outcome = oc, treatment = "treatment_2plus", instrument = "twin_like_proxy",
    controls_fml = controls_fml, weight_var = "mother_weight",
    cluster_var = "household_key", ar_grid = ar_grid,
    primary_callback = make_primary_cb("Twin1_pooled0_5", oc)),
    error = function(e) list(status = "error", message = conditionMessage(e)))
  res_t1$design <- "Twin1_pooled0_5"; res_t1$outcome <- oc
  checkpoint_fit(res_t1)

  res_ss <- tryCatch(fit_instrument_outcome(
    ss[eligible == TRUE], outcome = oc, treatment = "treatment_3plus", instrument = "samesex",
    controls_fml = controls_fml, weight_var = "mother_weight",
    cluster_var = "household_key", ar_grid = ar_grid,
    primary_callback = make_primary_cb("SameSex2_pooled0_5", oc)),
    error = function(e) list(status = "error", message = conditionMessage(e)))
  res_ss$design <- "SameSex2_pooled0_5"; res_ss$outcome <- oc
  checkpoint_fit(res_ss)

  for (ea in c(3, 5)) {
    t1_ea <- t1[event_age == ea]
    ss_ea <- ss[eligible == TRUE & event_age == ea]
    r1 <- tryCatch(fit_instrument_outcome(t1_ea, oc, "treatment_2plus", "twin_like_proxy",
                     controls_fml, "mother_weight", "household_key", ar_grid,
                     primary_callback = make_primary_cb(paste0("Twin1_event", ea), oc)),
                   error = function(e) list(status = "error", message = conditionMessage(e)))
    r1$design <- paste0("Twin1_event", ea); r1$outcome <- oc
    checkpoint_fit(r1)
    r2 <- tryCatch(fit_instrument_outcome(ss_ea, oc, "treatment_3plus", "samesex",
                     controls_fml, "mother_weight", "household_key", ar_grid,
                     primary_callback = make_primary_cb(paste0("SameSex2_event", ea), oc)),
                   error = function(e) list(status = "error", message = conditionMessage(e)))
    r2$design <- paste0("SameSex2_event", ea); r2$outcome <- oc
    checkpoint_fit(r2)
  }
}

results_dt <- data.table::rbindlist(lapply(results, function(r) {
  data.table(
    design = r$design %||% NA_character_, outcome = r$outcome %||% NA_character_,
    status = r$status %||% NA_character_,
    error_message = r$message %||% NA_character_,
    n_usable_prefit = r$n_usable %||% NA_integer_,
    n_households_prefit = r$n_households %||% NA_integer_,
    n_instrument_positive = r$n_instrument_positive %||% NA_integer_,
    rf_coef = r$rf_coef %||% NA_real_, rf_se = r$rf_se %||% NA_real_,
    rf_ci_lower = r$rf_ci_lower %||% NA_real_, rf_ci_upper = r$rf_ci_upper %||% NA_real_,
    rf_nobs_fit = r$rf_nobs_fit %||% NA_integer_, rf_error = r$rf_error %||% NA_character_,
    fs_coef = r$fs_coef %||% NA_real_, fs_se = r$fs_se %||% NA_real_,
    fs_ci_lower = r$fs_ci_lower %||% NA_real_, fs_ci_upper = r$fs_ci_upper %||% NA_real_,
    fs_nobs_fit = r$fs_nobs_fit %||% NA_integer_, fs_error = r$fs_error %||% NA_character_,
    first_stage_F = r$first_stage_F %||% NA_real_, first_stage_F_df = r$first_stage_F_df %||% NA_character_,
    iv_coef = r$iv_coef %||% NA_real_, iv_se = r$iv_se %||% NA_real_,
    iv_ci_lower = r$iv_ci_lower %||% NA_real_, iv_ci_upper = r$iv_ci_upper %||% NA_real_,
    iv_nobs_fit = r$iv_nobs_fit %||% NA_integer_, iv_error = r$iv_error %||% NA_character_,
    ar_summary_lower = r$ar_summary_lower %||% NA_real_, ar_summary_upper = r$ar_summary_upper %||% NA_real_,
    ar_fully_interior_bounded = r$ar_fully_interior_bounded %||% NA,
    ar_extent_unknown_due_to_errors = r$ar_extent_unknown_due_to_errors %||% NA,
    ar_n_components = r$ar_n_components %||% NA_integer_,
    ar_all_grid_points_accepted = r$ar_all_grid_points_accepted %||% NA,
    ar_empty_accepted_set = r$ar_empty_accepted_set %||% NA,
    ar_n_errors = r$ar_n_errors %||% NA_integer_,
    ar_grid_min = r$ar_grid_min %||% NA_real_, ar_grid_max = r$ar_grid_max %||% NA_real_,
    ar_components_json = as.character(r$ar_components_json %||% NA_character_))
}), fill = TRUE)

utils::write.csv(results_dt, file.path(outdir, "rf_fs_iv_ar_table.csv"), row.names = FALSE)
jsonlite::write_json(counts, file.path(outdir, "counts_final.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)

# Explicit fit-status tally so a run with error/insufficient_support rows
# cannot be reported as a clean success by omission.
status_tally <- results_dt[, .N, by = status]
jsonlite::write_json(status_tally, file.path(outdir, "fit_status_tally.json"),
                      auto_unbox = TRUE, pretty = TRUE)

write_status("complete", list(elapsed_sec = as.numeric(Sys.time() - t0, units = "secs"),
                               n_rows_table = nrow(results_dt),
                               fit_status_tally = as.list(setNames(status_tally$N, status_tally$status))))
cat("DRIVER_COMPLETE\n")
