# SameSex2 negative-rooms diagnostic driver. Compute-only (production
# mode) or DIAGNOSTIC_FIXTURE_MODE=1 (synthetic statefip 99, tiny fixture
# analytic_frames.rds + identity + saved receipt built fresh each run --
# for a genuine local end-to-end test of this exact driver, not just unit
# functions). Production default REQUIRES the real 51-state list and
# STATEFIP_SMOKE=50 (VT); DIAGNOSTIC_FIXTURE_MODE overrides both.
# Every case's receipt is written atomically before the next case runs.
suppressMessages({ library(data.table); library(fixest); library(jsonlite) })
options(warn = 1)
data.table::setDTthreads(1)
options(fixest_nthreads = 1)

args_env <- function(name, default = NULL) { v <- Sys.getenv(name, unset = ""); if (nzchar(v)) v else default }
ROOT <- args_env("ROOT")
# The national identity's lib_md5 was pinned against the ORIGINAL
# production twins_samesex_iv_lib.R (roster/instrument construction), NOT
# this diagnosis file -- verify_reproduction must check that hash, since
# that is the code whose reproduction on af$ss is actually being tested.
# The diagnosis lib's own hash is tracked separately (informational, not
# part of the reproduction gate).
ORIGINAL_LIB_PATH <- file.path(ROOT, "twins_samesex_iv_lib.R")
DIAG_LIB_PATH <- file.path(ROOT, "samesex_diagnosis_lib.R")
source(ORIGINAL_LIB_PATH)
source(DIAG_LIB_PATH)

FIXTURE_MODE <- isTRUE(as.logical(args_env("DIAGNOSTIC_FIXTURE_MODE", "0")))
NATIONAL_OUTDIR <- args_env("NATIONAL_OUTDIR")
PARTITION_DIR <- args_env("PARTITION_DIR",
  "/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_source_stage_20260921/partitions")
OUTDIR <- args_env("OUTDIR")
if (dir.exists(OUTDIR)) stop(sprintf("OUTDIR already exists, refusing to reuse: %s", OUTDIR))
dir.create(OUTDIR, recursive = TRUE)
STATUS_PATH <- args_env("STATUS_PATH", file.path(OUTDIR, "status_heartbeat.json"))
STATEFIP_SMOKE <- as.integer(args_env("STATEFIP_SMOKE", if (FIXTURE_MODE) "99" else "50"))
default_states <- if (FIXTURE_MODE) "99" else
  "1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56"
STATEFIP_LIST <- as.integer(strsplit(args_env("STATEFIP_LIST", default_states), ",")[[1]])
if (!FIXTURE_MODE && (length(STATEFIP_LIST) != 51 || STATEFIP_SMOKE != 50)) {
  stop("Production mode requires the real 51-state list and STATEFIP_SMOKE=50 (VT); set DIAGNOSTIC_FIXTURE_MODE=1 to override for local testing.")
}

t0 <- Sys.time()
atomic_rename <- function(tmp, final) if (!file.rename(tmp, final)) stop(sprintf("atomic rename failed: %s", final))
write_json_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp"); jsonlite::write_json(obj, tmp, auto_unbox = TRUE, pretty = TRUE, digits = 10, null = "null", na = "null")
  atomic_rename(tmp, path)
}
write_status <- function(phase, extra = list()) {
  write_json_atomic(c(list(phase = phase, generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                            elapsed_sec = as.numeric(Sys.time() - t0, units = "secs")), extra), STATUS_PATH)
}
n_stage_regressions <- 0L
bump <- function(k = 1L) {
  n_stage_regressions <<- n_stage_regressions + k
  if (n_stage_regressions > MAX_STAGE_REGRESSIONS)
    stop(sprintf("stage-regression budget exceeded: %d > %d", n_stage_regressions, MAX_STAGE_REGRESSIONS))
}
case_dir <- file.path(OUTDIR, "case_receipts"); dir.create(case_dir, recursive = TRUE)
save_case <- function(name, obj) write_json_atomic(obj, file.path(case_dir, paste0(name, ".json")))

## ---- Fixture-mode setup: build a tiny synthetic "national" cache in OUTDIR
if (FIXTURE_MODE) {
  fix_partition_dir <- file.path(OUTDIR, "fixture_partitions")
  dir.create(file.path(fix_partition_dir, "statefip_99"), recursive = TRUE)
  set.seed(77)
  n_fix <- 4000
  rows <- vector("list", n_fix); pid <- 0L
  for (i in seq_len(n_fix)) {
    yr <- sample(2005:2019, 1); ma <- sample(21:40, 1); nk <- sample(0:3, 1, prob = c(.15, .25, .35, .25))
    serial <- i
    mrow <- data.table(YEAR = yr, SAMPLE = yr * 100L + 1L, SERIAL = serial, PERNUM = 1L, SEX = 2L, AGE = ma,
                        RELATE = 1L, MOMLOC = 0L, NCHILD = nk, OWNERSHP = sample(1:2, 1),
                        ROOMS = sample(c(2:9, 28), 1), BEDROOMS = sample(1:6, 1), PERWT = sample(20:200, 1),
                        HHWT = sample(20:200, 1), STATEFIP = 99L, RACE = sample(1:4, 1), EDUC = sample(1:10, 1),
                        MARST = 1L, FERTYR = sample(c(1, 2, NA), 1))
    kid_rows <- list()
    if (nk > 0) {
      maxk <- min(ma - 15, 18)
      if (maxk >= 0) {
        ages <- sort(sample(0:maxk, nk, replace = TRUE), decreasing = TRUE)
        for (k in seq_len(nk)) kid_rows[[k]] <- data.table(YEAR = yr, SAMPLE = yr * 100L + 1L, SERIAL = serial,
          PERNUM = 1L + k, SEX = sample(1:2, 1), AGE = ages[k], RELATE = 3L, MOMLOC = 1L, NCHILD = NA_integer_,
          OWNERSHP = mrow$OWNERSHP, ROOMS = mrow$ROOMS, BEDROOMS = mrow$BEDROOMS, PERWT = mrow$PERWT,
          HHWT = mrow$HHWT, STATEFIP = 99L, RACE = mrow$RACE, EDUC = NA_integer_, MARST = NA_integer_, FERTYR = NA_integer_)
      }
    }
    rows[[i]] <- rbindlist(c(list(mrow), kid_rows), fill = TRUE)
  }
  dt_fix <- rbindlist(rows, fill = TRUE)
  dt_fix[, CBSERIAL := SERIAL]; dt_fix[, CLUSTER := SERIAL]; dt_fix[, PUMA := 100]; dt_fix[, STRATA := 1]
  dt_fix[, GQ := 1]; dt_fix[, OWNERSHPD := OWNERSHP]; dt_fix[, NCHLT5 := 0]; dt_fix[, ELDCH := 0]; dt_fix[, YNGCH := 0]
  dt_fix[, POPLOC := 0]; dt_fix[, SPLOC := 0]
  dt_fix[, ROOMS_RAW := ROOMS]; dt_fix[, BEDROOMS_RAW := BEDROOMS]; dt_fix[, OWNERSHP_RAW := OWNERSHP]
  saveRDS(dt_fix, file.path(fix_partition_dir, "statefip_99", "housing_narrow.rds"))
  PARTITION_DIR <- fix_partition_dir

  # Build the tiny "national" analytic frame + identity + saved receipt by
  # running the SAME extraction/roster pipeline once here, exactly
  # mirroring what the production driver does, so Stage A's reproduction
  # gate is a genuine (not fabricated) round-trip.
  meta_fix <- extract_samesex_roster_metadata_one_state(dt_fix)
  gate0 <- apply_source_sample_gate(dt_fix)
  built0 <- build_mother_roster(gate0$data)
  mr0 <- add_outcomes(built0$mother_rows)
  mg0 <- apply_oldest_child_minor_gate(mr0)
  mr1 <- mg0$data[PERWT > 0 & is.finite(PERWT) & AGE_norm >= 21 & AGE_norm <= 35]
  ss0 <- build_samesex2(mr1, age_grid = 0:5)
  ss0[, mother_weight := PERWT]
  outcomes_fix <- c("ROOMS_out", "OWNERSHP_out", "BEDROOMS_out")
  controls_fix <- "i(RACE)"
  ss_slim_fix <- ss0[eligible == TRUE, c("person_key", "household_key", "mother_weight", "samesex",
                                          "treatment_3plus", "event_age", outcomes_fix, "RACE"), with = FALSE]
  nat_outdir <- file.path(OUTDIR, "fixture_national_outdir")
  dir.create(file.path(nat_outdir, "fit_receipts"), recursive = TRUE)
  af_path <- file.path(nat_outdir, "analytic_frames.rds")
  saveRDS(list(t1 = data.table(), ss = ss_slim_fix), af_path)
  ref_fit <- rf_fs_diagnostic_fit(ss_slim_fix, "ROOMS_out", "treatment_3plus", "samesex", controls_fix,
                                   "mother_weight", "household_key")
  saved_receipt_fix <- list(outcome = "ROOMS_out", instrument = "samesex", treatment = "treatment_3plus",
                             rf_coef = ref_fit$rf_coef, rf_se = ref_fit$rf_se, fs_coef = ref_fit$fs_coef,
                             fs_se = ref_fit$fs_se, rf_nobs_fit = ref_fit$rf_nobs, fs_nobs_fit = ref_fit$fs_nobs)
  write_json_atomic(saved_receipt_fix, file.path(nat_outdir, "fit_receipts", "SameSex2_pooled0_5__ROOMS_out.json"))
  identity_fix <- list(analytic_frames_path = af_path, analytic_frames_size_bytes = file.info(af_path)$size,
                        analytic_frames_md5 = unname(tools::md5sum(af_path)), lib_md5 = unname(tools::md5sum(ORIGINAL_LIB_PATH)),
                        controls_fml = controls_fix, weight_var = "mother_weight", cluster_var = "household_key",
                        outcomes = outcomes_fix)
  write_json_atomic(identity_fix, file.path(nat_outdir, "analytic_frames_identity.json"))
  NATIONAL_OUTDIR <- nat_outdir
  write_status("fixture_mode_setup_complete", list(n_fixture_rows = nrow(dt_fix), n_ss_eligible = nrow(ss_slim_fix)))
}

## ---- Stage A: identity + numeric reproduction gate ------------------------
write_status("stage_A_reproduction_gate")
identity <- jsonlite::fromJSON(file.path(NATIONAL_OUTDIR, "analytic_frames_identity.json"))
saved_receipt_path <- file.path(NATIONAL_OUTDIR, "fit_receipts", "SameSex2_pooled0_5__ROOMS_out.json")
repro <- verify_reproduction(national_outdir = NATIONAL_OUTDIR, identity = identity, lib_path = ORIGINAL_LIB_PATH,
                              saved_receipt_path = saved_receipt_path)
bump(2)
save_case("A_reproduction_gate", repro[intersect(names(repro), c("pass", "gate", "stage", "checks"))])
write_status("stage_A_done", list(reproduction_pass = repro$pass))
if (!isTRUE(repro$pass)) stop(sprintf("Stage A reproduction gate FAILED at stage '%s' -- stopping before any diagnostic.", repro$stage))

af <- repro$af
outcomes <- identity$outcomes
controls_fml <- identity$controls_fml
weight_var <- identity$weight_var
cluster_var <- identity$cluster_var

## ---- Stage B: event age 0:5 x outcomes RF/FS, EACH case persisted --------
write_status("stage_B_event_age")
cases_B <- event_age_case_list(outcomes, ages = 0:5)
summary_rows <- vector("list", length(cases_B))
for (i in seq_along(cases_B)) {
  oc <- cases_B[[i]]$outcome; ea <- cases_B[[i]]$event_age
  d <- af$ss[event_age == ea]
  r <- rf_fs_diagnostic_fit(d, oc, "treatment_3plus", "samesex", controls_fml, weight_var, cluster_var)
  bump(2)
  save_case(sprintf("B_event_age_%d__%s", ea, oc), r)
  summary_rows[[i]] <- event_age_summary_row(r, oc, ea)
  write_status("stage_B_progress", list(case = i, of = length(cases_B), outcome = oc, event_age = ea))
}
tabB <- data.table::rbindlist(summary_rows, fill = TRUE)
utils::write.csv(tabB, file.path(OUTDIR, "B_event_age_rf_fs.csv"), row.names = FALSE)
write_status("stage_B_done", list(n_cases = nrow(tabB), n_full_fit = sum(tabB$status == "full_fit", na.rm = TRUE)))

## ---- Stage C-smoke: one state, mechanical loop check + strict partial gate
write_status("stage_C_smoke", list(statefip = STATEFIP_SMOKE))
p_smoke <- sprintf("%s/statefip_%02d/housing_narrow.rds", PARTITION_DIR, STATEFIP_SMOKE)
dt_smoke <- readRDS(p_smoke)
meta_smoke <- extract_samesex_roster_metadata_one_state(dt_smoke)
rm(dt_smoke); gc(FALSE)
audit_smoke <- sex_order_audit(meta_smoke)
overlap_smoke <- af$ss[person_key %in% meta_smoke[eligible == TRUE]$person_key]
n_overlap <- nrow(overlap_smoke)
smoke_gate <- if (n_overlap > 0) verify_key_coverage_and_recompute(meta_smoke[eligible == TRUE & person_key %in% overlap_smoke$person_key], overlap_smoke) else NULL
smoke_pass <- n_overlap > 0 && !is.null(smoke_gate) && isTRUE(smoke_gate$one_to_one_and_recomputed_equal)
save_case("C_smoke", list(statefip = STATEFIP_SMOKE, n_metadata_rows = nrow(meta_smoke), n_overlap = n_overlap,
                           audit = audit_smoke, gate = smoke_gate, smoke_pass = smoke_pass))
write_status("stage_C_smoke_done", list(n_metadata_rows = nrow(meta_smoke), n_overlap = n_overlap, smoke_pass = smoke_pass))
if (!smoke_pass) {
  stop(sprintf("Stage C smoke FAILED: n_overlap=%d smoke_pass=%s -- stopping before the full state pass (no overlap, mismatch, duplicates, or other gate failure).",
               n_overlap, smoke_pass))
}

## ---- Stage C-full: all states, memory-bounded, reuses VT smoke metadata --
write_status("stage_C_full_start", list(n_states = length(STATEFIP_LIST)))
meta_list <- vector("list", length(STATEFIP_LIST))
audit_list <- vector("list", length(STATEFIP_LIST))
for (i in seq_along(STATEFIP_LIST)) {
  sf <- STATEFIP_LIST[i]
  if (sf == STATEFIP_SMOKE) {
    m <- meta_smoke  # reuse already-read/extracted VT metadata, no second read
  } else {
    p <- sprintf("%s/statefip_%02d/housing_narrow.rds", PARTITION_DIR, sf)
    dt_state <- readRDS(p)
    m <- extract_samesex_roster_metadata_one_state(dt_state)
    rm(dt_state); gc(FALSE)
  }
  meta_list[[i]] <- m
  audit_list[[i]] <- c(list(statefip = sf, n_rows = nrow(m)), sex_order_audit(m))
  write_status("stage_C_full_progress", list(loaded = i, of = length(STATEFIP_LIST), statefip = sf))
}
metadata_all <- data.table::rbindlist(meta_list, use.names = TRUE, fill = TRUE)
rm(meta_list, meta_smoke); gc(FALSE)
audit_all <- sex_order_audit(metadata_all)
audit_all$n_full_pool <- nrow(metadata_all)
audit_all$n_eligible_primary <- sum(metadata_all$eligible)
save_case("C_full_audit", list(n_states = length(STATEFIP_LIST), per_state = audit_list, pooled_audit = audit_all))

coverage <- verify_key_coverage_and_recompute(metadata_all, af$ss)
save_case("C_full_coverage_gate", coverage)
write_status("stage_C_full_done", list(coverage_pass = coverage$one_to_one_and_recomputed_equal))
if (!isTRUE(coverage$one_to_one_and_recomputed_equal)) {
  stop("Stage C full coverage/recompute gate FAILED (exact 1:1 person-key match required) -- stopping before D/E (no sample silently changed).")
}

## ---- Persist metadata BEFORE D/E ------------------------------------------
meta_rds_path <- file.path(OUTDIR, "samesex_roster_metadata.rds")
tmp_meta <- paste0(meta_rds_path, ".tmp")
saveRDS(metadata_all, tmp_meta); atomic_rename(tmp_meta, meta_rds_path)
metadata_identity <- list(
  metadata_rds_path = meta_rds_path, metadata_rds_size_bytes = file.info(meta_rds_path)$size,
  metadata_rds_md5 = unname(tools::md5sum(meta_rds_path)),
  original_production_lib_md5 = unname(tools::md5sum(ORIGINAL_LIB_PATH)),
  diagnosis_lib_md5 = unname(tools::md5sum(DIAG_LIB_PATH)),
  national_analytic_frames_md5 = identity$analytic_frames_md5,
  n_states = length(STATEFIP_LIST), statefip_list = STATEFIP_LIST,
  n_full_pool = nrow(metadata_all), n_eligible_primary = sum(metadata_all$eligible),
  coverage_gate = coverage, generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
write_json_atomic(metadata_identity, file.path(OUTDIR, "samesex_roster_metadata_identity.json"))
af$t1 <- NULL; gc(FALSE)

## ---- Merge metadata onto the cached ss for D/E (exact person_key join) ---
ss_meta <- merge(af$ss, metadata_all[, .(person_key, sex1, sex2, a1, a2, a3, linked_child_count, NCHILD_norm)],
                  by = "person_key")
if (nrow(ss_meta) != nrow(af$ss)) stop("post-merge row count mismatch -- coverage gate should have prevented this")

## ---- Stage D: additive sex controls + joint BB/GG + descriptive cells ----
write_status("stage_D_sex_controls")
cellcounts <- sex_cell_counts(ss_meta, outcomes, weight_var)
utils::write.csv(cellcounts, file.path(OUTDIR, "D_sex_cell_counts.csv"), row.names = FALSE)
for (oc in outcomes) {
  add_fit <- additive_sex_control_rf_fs(ss_meta, oc, controls_fml, weight_var, cluster_var)
  bump(2); save_case(paste0("D_additive_sex_control__", oc), add_fit)
  write_status("stage_D_progress", list(outcome = oc, sub = "additive"))
  joint_fit <- joint_bb_gg_rf_fs(ss_meta, oc, controls_fml, weight_var, cluster_var)
  bump(2); save_case(paste0("D_joint_bb_gg__", oc), joint_fit)
  write_status("stage_D_progress", list(outcome = oc, sub = "joint"))
}
write_status("stage_D_done")

## ---- Stage E: ambiguity sensitivity excluding a2==a3 ----------------------
write_status("stage_E_ambiguity_sensitivity")
sens <- ambiguity_sensitivity_exclude_a2_eq_a3(ss_meta)
for (oc in outcomes) {
  r <- rf_fs_diagnostic_fit(sens$data, oc, "treatment_3plus", "samesex", controls_fml, weight_var, cluster_var)
  bump(2)
  r$n_excluded_a2_eq_a3 <- sens$n_excluded_a2_eq_a3; r$n_before <- sens$n_before
  save_case(paste0("E_ambiguity_sensitivity__", oc), r)
  write_status("stage_E_progress", list(outcome = oc))
}
write_status("stage_E_done", list(n_excluded_a2_eq_a3 = sens$n_excluded_a2_eq_a3))

all_case_files <- list.files(case_dir, pattern = "\\.json$")
write_status("complete", list(n_stage_regressions = n_stage_regressions, n_case_receipts = length(all_case_files)))
cat(sprintf("SAMESEX_DIAGNOSIS_COMPLETE n_stage_regressions=%d n_case_receipts=%d\n",
            n_stage_regressions, length(all_case_files)))
