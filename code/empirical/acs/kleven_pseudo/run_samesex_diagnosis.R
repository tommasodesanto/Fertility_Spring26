# SameSex2 negative-rooms diagnostic driver. Compute-only. Stages:
#   A) identity+numeric reproduction gate against the saved national
#      SameSex2 pooled ROOMS receipt (job 18274624) -- stop if it fails.
#   B) event age 0:5 x {rooms,bedrooms,ownership} RF/FS on the cached
#      national ss frame (analytic_frames.rds) -- no new extraction.
#   C-smoke) per-state metadata extraction on ONE state (STATEFIP_SMOKE,
#      default VT=50) as a mechanical smoke of the exact production
#      roster pipeline, with a partial (this-state-only) recompute check
#      against the cached ss -- NOT the full coverage gate yet.
#   C-full) the same per-state extraction across all 51 states, releasing
#      each wide state table immediately, then the STRICT 1:1 person-key
#      coverage + recompute-equality gate against the cached ss (this is
#      the real, exact-coverage gate; the smoke stage cannot achieve full
#      coverage on one state alone).
#   D) additive first/second-child-sex RF/FS + joint BB/GG RF/FS + BB/BG/
#      GB/GG descriptive cell counts, gated on C-full passing.
#   E) ambiguity sensitivity excluding a2==a3, gated on C-full passing.
# Every case's receipt is written atomically before the next case runs.
# Env vars: ROOT (code root, unchanged production+diagnosis libs),
# PARTITION_DIR, NATIONAL_OUTDIR (completed national run's small outputs,
# for identity/receipt/analytic_frames.rds), OUTDIR (fresh), STATUS_PATH,
# STATEFIP_SMOKE (default 50), STATEFIP_LIST (default all 51).
suppressMessages({ library(data.table); library(fixest); library(jsonlite) })
options(warn = 1)
data.table::setDTthreads(1)
options(fixest_nthreads = 1)

args_env <- function(name, default = NULL) { v <- Sys.getenv(name, unset = ""); if (nzchar(v)) v else default }
ROOT <- args_env("ROOT")
source(file.path(ROOT, "twins_samesex_iv_lib.R"))
source(file.path(ROOT, "samesex_diagnosis_lib.R"))

NATIONAL_OUTDIR <- args_env("NATIONAL_OUTDIR")
PARTITION_DIR <- args_env("PARTITION_DIR",
  "/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_source_stage_20260921/partitions")
OUTDIR <- args_env("OUTDIR")
if (dir.exists(OUTDIR)) stop(sprintf("OUTDIR already exists, refusing to reuse: %s", OUTDIR))
dir.create(OUTDIR, recursive = TRUE)
STATUS_PATH <- args_env("STATUS_PATH", file.path(OUTDIR, "status_heartbeat.json"))
STATEFIP_SMOKE <- as.integer(args_env("STATEFIP_SMOKE", "50"))
STATEFIP_LIST <- as.integer(strsplit(args_env(
  "STATEFIP_LIST",
  "1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56"
), ",")[[1]])

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
bump <- function(k = 1L) n_stage_regressions <<- n_stage_regressions + k

case_dir <- file.path(OUTDIR, "case_receipts"); dir.create(case_dir, recursive = TRUE)
save_case <- function(name, obj) write_json_atomic(obj, file.path(case_dir, paste0(name, ".json")))

## ---- Stage A: identity + numeric reproduction gate ------------------------
write_status("stage_A_reproduction_gate")
identity <- jsonlite::fromJSON(file.path(NATIONAL_OUTDIR, "analytic_frames_identity.json"))
saved_receipt_path <- file.path(NATIONAL_OUTDIR, "fit_receipts", "SameSex2_pooled0_5__ROOMS_out.json")
repro <- verify_reproduction(
  analytic_frames_path = identity$analytic_frames_path,
  expected_size = identity$analytic_frames_size_bytes, expected_md5 = identity$analytic_frames_md5,
  saved_receipt_path = saved_receipt_path,
  controls_fml = identity$controls_fml, weight_var = identity$weight_var, cluster_var = identity$cluster_var
)
bump(2)  # RF + FS
save_case("A_reproduction_gate", repro[c("pass", "gate", "stage", "checks")])
write_status("stage_A_done", list(reproduction_pass = repro$pass))
if (!isTRUE(repro$pass)) stop(sprintf("Stage A reproduction gate FAILED at stage '%s' -- stopping before any diagnostic.", repro$stage))

af <- readRDS(identity$analytic_frames_path)
outcomes <- identity$outcomes
controls_fml <- identity$controls_fml
weight_var <- identity$weight_var
cluster_var <- identity$cluster_var

## ---- Stage B: event age 0:5 x 3 outcomes RF/FS ----------------------------
write_status("stage_B_event_age")
tabB <- event_age_rf_fs_diagnostic(af$ss, outcomes, controls_fml, weight_var, cluster_var, ages = 0:5)
bump(nrow(tabB) * 2)
utils::write.csv(tabB, file.path(OUTDIR, "B_event_age_rf_fs.csv"), row.names = FALSE)
save_case("B_event_age_summary", list(n_rows = nrow(tabB), n_full_fit = sum(tabB$status == "full_fit", na.rm = TRUE)))
write_status("stage_B_done", list(n_rows = nrow(tabB)))

## ---- Stage C-smoke: one state, mechanical loop check + partial recompute -
write_status("stage_C_smoke", list(statefip = STATEFIP_SMOKE))
p_smoke <- sprintf("%s/statefip_%02d/housing_narrow.rds", PARTITION_DIR, STATEFIP_SMOKE)
dt_smoke <- readRDS(p_smoke)
meta_smoke <- extract_samesex_roster_metadata_one_state(dt_smoke)
rm(dt_smoke); gc(FALSE)
audit_smoke <- sex_order_audit(meta_smoke)
# Partial check: keys from this ONE state that also exist in the cached
# national ss must recompute-match; full 1:1 coverage is impossible from a
# single state and is NOT asserted here (that is Stage C-full's job).
overlap_smoke <- af$ss[person_key %in% meta_smoke[eligible == TRUE]$person_key]
partial_chk <- if (nrow(overlap_smoke) > 0) verify_key_coverage_and_recompute(meta_smoke, overlap_smoke) else
  list(note = "no overlapping keys found for this state in the cached ss (unexpected if STATEFIP_SMOKE is a real contributing state)")
save_case("C_smoke", list(statefip = STATEFIP_SMOKE, n_metadata_rows = nrow(meta_smoke),
                           audit = audit_smoke, partial_recompute_check = partial_chk))
write_status("stage_C_smoke_done", list(n_metadata_rows = nrow(meta_smoke)))
if (!is.null(partial_chk$n_recompute_mismatches) && partial_chk$n_recompute_mismatches > 0) {
  stop("Stage C smoke found recompute mismatches on the overlapping keys -- stopping before the full 51-state pass.")
}
rm(meta_smoke, overlap_smoke); gc(FALSE)

## ---- Stage C-full: all 51 states, memory-bounded, then STRICT gate -------
write_status("stage_C_full_start", list(n_states = length(STATEFIP_LIST)))
meta_list <- vector("list", length(STATEFIP_LIST))
audit_list <- vector("list", length(STATEFIP_LIST))
for (i in seq_along(STATEFIP_LIST)) {
  sf <- STATEFIP_LIST[i]
  p <- sprintf("%s/statefip_%02d/housing_narrow.rds", PARTITION_DIR, sf)
  dt_state <- readRDS(p)
  m <- extract_samesex_roster_metadata_one_state(dt_state)
  rm(dt_state); gc(FALSE)
  meta_list[[i]] <- m
  audit_list[[i]] <- c(list(statefip = sf, n_rows = nrow(m)), sex_order_audit(m))
  write_status("stage_C_full_progress", list(loaded = i, of = length(STATEFIP_LIST), statefip = sf))
}
metadata_all <- data.table::rbindlist(meta_list, use.names = TRUE, fill = TRUE)
rm(meta_list); gc(FALSE)
audit_all <- sex_order_audit(metadata_all)
save_case("C_full_audit", list(n_states = length(STATEFIP_LIST), per_state = audit_list, pooled_audit = audit_all))

coverage <- verify_key_coverage_and_recompute(metadata_all, af$ss)
save_case("C_full_coverage_gate", coverage)
write_status("stage_C_full_done", list(coverage_pass = coverage$one_to_one_and_recomputed_equal))
if (!isTRUE(coverage$one_to_one_and_recomputed_equal)) {
  stop("Stage C full 51-state coverage/recompute gate FAILED -- stopping before D/E (no sample silently changed).")
}

## ---- Merge metadata onto the cached ss for D/E (exact person_key join) ---
ss_meta <- merge(af$ss, metadata_all[, .(person_key, sex1, sex2, a3)], by = "person_key")
stopifnot(nrow(ss_meta) == nrow(af$ss))  # coverage gate already guarantees this; assert defensively

## ---- Stage D: additive sex controls + joint BB/GG + descriptive cells ----
write_status("stage_D_sex_controls")
cellcounts <- sex_cell_counts(ss_meta, outcomes, weight_var)
utils::write.csv(cellcounts, file.path(OUTDIR, "D_sex_cell_counts.csv"), row.names = FALSE)
for (oc in outcomes) {
  add_fit <- additive_sex_control_rf_fs(ss_meta, oc, controls_fml, weight_var, cluster_var)
  bump(2)
  save_case(paste0("D_additive_sex_control__", oc), add_fit)
  joint_fit <- joint_bb_gg_rf_fs(ss_meta, oc, controls_fml, weight_var, cluster_var)
  bump(2)
  save_case(paste0("D_joint_bb_gg__", oc), joint_fit)
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
}
write_status("stage_E_done", list(n_excluded_a2_eq_a3 = sens$n_excluded_a2_eq_a3))

write_status("complete", list(n_stage_regressions = n_stage_regressions))
cat(sprintf("SAMESEX_DIAGNOSIS_COMPLETE n_stage_regressions=%d\n", n_stage_regressions))
