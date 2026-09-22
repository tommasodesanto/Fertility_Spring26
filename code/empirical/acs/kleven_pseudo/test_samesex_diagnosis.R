# Local synthetic tests for samesex_diagnosis_lib.R. No cluster access; no
# national analytic_frames.rds locally, so verify_reproduction() is tested
# via a small in-memory equivalent (fixture DGP), not the real national
# frame (that check runs on compute per the diagnostic plan).
.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.script_path) == 1) dirname(normalizePath(.script_path)) else getwd()
source(file.path(.script_dir, "twins_samesex_iv_lib.R"))
source(file.path(.script_dir, "samesex_diagnosis_lib.R"))
suppressMessages(library(data.table))

stopifnot_msg <- function(cond, msg) if (!isTRUE(cond)) stop(paste("FAIL:", msg)) else cat("PASS:", msg, "\n")

## ---- sex_order_audit: valid sex, age order, a2==a3 tie, NCHILD mismatch, dup keys
meta <- data.table(
  person_key = c("p1", "p2", "p3", "p3"),  # p3 duplicated on purpose
  household_key = c("h1", "h2", "h3", "h4"),
  sex1 = c(1, 2, 9, 1),    # row3 invalid sex1
  sex2 = c(2, 1, 2, 2),
  a1 = c(5, 4, 3, 6), a2 = c(3, 2, 3, 4), a3 = c(1, NA, 3, 4),  # row3: a2==a3 tie; row4: a2==a3 tie
  linked_child_count = c(3, 2, 3, 3), NCHILD_norm = c(3, 2, 2, 3)  # row3: NCHILD mismatch (2 vs 3)
)
audit <- sex_order_audit(meta)
stopifnot_msg(audit$n_sex1_invalid == 1, "sex_order_audit: detects 1 invalid sex1 (code 9)")
stopifnot_msg(audit$n_third_child_age_tie_a2_eq_a3 == 2, "sex_order_audit: detects 2 rows with a2==a3")
stopifnot_msg(audit$n_nchild_link_mismatch == 1, "sex_order_audit: detects 1 NCHILD-vs-linked mismatch")
stopifnot_msg(audit$n_duplicate_person_key == 1, "sex_order_audit: detects 1 duplicated person_key")
stopifnot_msg(audit$n_age_order_violation_a1_lt_a2 == 0, "sex_order_audit: no a1<a2 violations in fixture (order preserved)")

## ---- ambiguity_sensitivity_exclude_a2_eq_a3 -------------------------------
ss_fix <- data.table(a2 = c(3, 2, NA, 5), a3 = c(3, NA, 1, 4), samesex = c(1, 0, 1, 0))
sens <- ambiguity_sensitivity_exclude_a2_eq_a3(ss_fix)
stopifnot_msg(sens$n_excluded_a2_eq_a3 == 1, "ambiguity sensitivity: excludes exactly the 1 row with a2==a3")
stopifnot_msg(nrow(sens$data) == 3, "ambiguity sensitivity: retains the other 3 rows (incl. NA a2 or a3)")

## ---- sex_cell_counts: BB/BG/GB/GG ----------------------------------------
cellfix <- data.table(sex1 = c(1, 1, 2, 2, 1), sex2 = c(1, 2, 1, 2, 1),
                       treatment_3plus = c(1, 0, 1, 0, 0), mother_weight = c(1, 2, 3, 4, 5),
                       ROOMS_out = c(6, 5, 7, 4, 6))
cells <- sex_cell_counts(cellfix, outcomes = "ROOMS_out", weight_var = "mother_weight")
stopifnot_msg(all(c("BB", "BG", "GB", "GG") %in% cells$cell), "sex_cell_counts: all 4 ordered cells present")
stopifnot_msg(cells[cell == "BB", n_unweighted] == 2, "sex_cell_counts: BB has 2 unweighted obs (rows 1,5)")
stopifnot_msg(cells[cell == "BB", n_weighted] == 6, "sex_cell_counts: BB weighted N = 1+5 = 6")

## ---- add_first_second_child_sex_controls ----------------------------------
ctrlfix <- data.table(sex1 = c(1, 2), sex2 = c(2, 2))
augmented <- add_first_second_child_sex_controls(ctrlfix)
stopifnot_msg(identical(augmented$firstchildboy, c(1, 0)), "first/second-child-sex controls: firstchildboy correct")
stopifnot_msg(identical(augmented$secondchildboy, c(0, 0)), "first/second-child-sex controls: secondchildboy correct")

## ---- verify_key_coverage_and_recompute: matched vs mismatched cases ------
metadata2 <- data.table(household_key = c("h1", "h2", "h3"), eligible = c(TRUE, TRUE, TRUE),
                         samesex = c(1, 0, 1), treatment_3plus = c(1, 0, 0), event_age = c(2, 3, 1))
cached_match <- data.table(household_key = c("h1", "h2", "h3"),
                            samesex = c(1, 0, 1), treatment_3plus = c(1, 0, 0), event_age = c(2, 3, 1))
chk_match <- verify_key_coverage_and_recompute(metadata2, cached_match)
stopifnot_msg(isTRUE(chk_match$one_to_one_and_recomputed_equal), "key coverage: matching metadata/cache passes")

cached_mismatch <- data.table::copy(cached_match)
cached_mismatch[household_key == "h3", event_age := 5]
chk_mismatch <- verify_key_coverage_and_recompute(metadata2, cached_mismatch)
stopifnot_msg(!isTRUE(chk_mismatch$one_to_one_and_recomputed_equal) && chk_mismatch$n_recompute_mismatches == 1,
              "key coverage: recompute mismatch (event_age) correctly detected, not silently accepted")

cached_extra <- rbindlist(list(cached_match, data.table(household_key = "h4", samesex = 1, treatment_3plus = 1, event_age = 0)))
chk_uncovered <- verify_key_coverage_and_recompute(metadata2, cached_extra)
stopifnot_msg(chk_uncovered$n_cached_unmatched == 1 && !isTRUE(chk_uncovered$one_to_one_and_recomputed_equal),
              "key coverage: a cached household missing from metadata is caught, not silently dropped")

## ---- event_age_rf_fs_diagnostic: all 6 ages kept on a synthetic DGP -------
set.seed(11)
n <- 3000
mat_age <- round(runif(n, 25, 40))
ea <- sample(0:5, n, replace = TRUE)
hh <- paste0("hh", seq_len(n))
Z <- rbinom(n, 1, 0.3)
D <- rbinom(n, 1, pmin(pmax(0.3 + 0.1 * Z, 0.01), 0.99))
Y <- 5 - 0.3 * D + rnorm(n)
ss_ea <- data.table(ROOMS_out = Y, treatment_3plus = D, samesex = Z, mat_age = mat_age,
                     mother_weight = 1, household_key = hh, event_age = ea)
tab <- event_age_rf_fs_diagnostic(ss_ea, outcomes = "ROOMS_out", controls_fml = "mat_age",
                                   weight_var = "mother_weight", cluster_var = "household_key", ages = 0:5)
stopifnot_msg(nrow(tab) == 6, "event-age diagnostic: all 6 ages (0:5) reported, none dropped")
stopifnot_msg(all(!is.na(tab$rf_coef)), "event-age diagnostic: every age produced a fitted RF coefficient")

## ---- extract_samesex_roster_metadata_one_state: runs on existing local fixture
fixture_path <- "/tmp/acs_iv_fixture/partitions/statefip_99/housing_narrow.rds"
if (file.exists(fixture_path)) {
  dt_state <- readRDS(fixture_path)
  meta_state <- extract_samesex_roster_metadata_one_state(dt_state)
  stopifnot_msg(nrow(meta_state) > 0, "metadata extraction: produces a non-empty metadata table on the local fixture")
  stopifnot_msg(all(c("sex1", "sex2", "a1", "a2", "a3", "linked_child_count", "NCHILD_norm",
                       "household_key", "person_key") %in% names(meta_state)),
                "metadata extraction: all required metadata columns present")
  audit_state <- sex_order_audit(meta_state)
  cat(sprintf("Local fixture metadata audit: n=%d sex1_invalid=%d a2eqa3=%d nchild_mismatch=%d dup_hh=%d\n",
              audit_state$n_rows, audit_state$n_sex1_invalid, audit_state$n_third_child_age_tie_a2_eq_a3,
              audit_state$n_nchild_link_mismatch, audit_state$n_duplicate_household_key))
} else {
  cat("SKIP: local fixture RDS not found at", fixture_path, "-- metadata-extraction-on-real-data check skipped\n")
}

cat("\nALL SAMESEX DIAGNOSIS TESTS PASSED\n")
