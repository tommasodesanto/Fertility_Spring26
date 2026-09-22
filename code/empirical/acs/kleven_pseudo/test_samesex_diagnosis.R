# Local synthetic tests for samesex_diagnosis_lib.R. No cluster access; no
# national analytic_frames.rds locally, so verify_reproduction()'s numeric
# stage is exercised via a small in-memory fixture DGP with matching
# identity fields, not the real national frame (that check runs on
# compute per the diagnostic plan). The fixture in this file is
# SYNTHETIC, not real ACS data -- any counts reported below (e.g. tie
# rates on the local fixture) describe fixture properties only, not a
# real-data finding.
.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.script_path) == 1) dirname(normalizePath(.script_path)) else getwd()
source(file.path(.script_dir, "twins_samesex_iv_lib.R"))
source(file.path(.script_dir, "samesex_diagnosis_lib.R"))
suppressMessages(library(data.table))

stopifnot_msg <- function(cond, msg) if (!isTRUE(cond)) stop(paste("FAIL:", msg)) else cat("PASS:", msg, "\n")

## ---- sex_order_audit: valid sex, order, a2==a3 tie, NCHILD mismatch, dup PERSON key, multi-mother HH
meta <- data.table(
  person_key = c("p1", "p2", "p3", "p4"),
  household_key = c("h1", "h1", "h2", "h3"),  # h1 has 2 distinct mothers -- legitimate, not an error
  sex1 = c(1, 2, 9, 1), sex2 = c(2, 1, 2, 2),  # row3 invalid sex1
  a1 = c(5, 4, 3, 6), a2 = c(3, 2, 3, 4), a3 = c(1, NA, 3, 4),  # row3,row4: a2==a3 tie
  linked_child_count = c(3, 2, 3, 3), NCHILD_norm = c(3, 2, 2, 3)  # row3: NCHILD mismatch
)
audit <- sex_order_audit(meta)
stopifnot_msg(audit$n_sex1_invalid == 1, "sex_order_audit: detects 1 invalid sex1 (code 9)")
stopifnot_msg(audit$n_third_child_age_tie_a2_eq_a3 == 2, "sex_order_audit: detects 2 rows with a2==a3")
stopifnot_msg(audit$n_nchild_link_mismatch == 1, "sex_order_audit: detects 1 NCHILD-vs-linked mismatch")
stopifnot_msg(audit$n_duplicate_person_key == 0, "sex_order_audit: no duplicate PERSON keys in fixture")
stopifnot_msg(audit$n_households_with_multiple_mothers == 1, "sex_order_audit: h1's 2 mothers flagged as multi-mother HH, not duplicates")

meta_dup <- rbindlist(list(meta, meta[1]))  # duplicate person_key p1, identical values
audit_dup <- sex_order_audit(meta_dup)
stopifnot_msg(audit_dup$n_duplicate_person_key == 1, "sex_order_audit: duplicate PERSON key flagged even with identical values")

## ---- verify_key_coverage_and_recompute: PERSON-KEY join, multi-mother HH allowed
metadata2 <- data.table(
  person_key = c("p1", "p2", "p3"), household_key = c("h1", "h1", "h2"),  # two mothers share h1
  eligible = c(TRUE, TRUE, TRUE), samesex = c(1, 0, 1), treatment_3plus = c(1, 0, 0), event_age = c(2, 3, 1)
)
cached_match <- data.table(
  person_key = c("p1", "p2", "p3"), household_key = c("h1", "h1", "h2"),
  samesex = c(1, 0, 1), treatment_3plus = c(1, 0, 0), event_age = c(2, 3, 1)
)
chk <- verify_key_coverage_and_recompute(metadata2, cached_match)
stopifnot_msg(isTRUE(chk$one_to_one_and_recomputed_equal), "key coverage: person-key join, shared household h1 for 2 mothers, passes (not flagged duplicate)")

## two distinct mothers same HH with DIFFERENT Z/D -- must pass (legitimate), not be treated as a conflict
metadata3 <- data.table(person_key = c("p1", "p2"), household_key = c("hA", "hA"),
                         eligible = c(TRUE, TRUE), samesex = c(1, 0), treatment_3plus = c(1, 0), event_age = c(2, 4))
cached3 <- data.table(person_key = c("p1", "p2"), household_key = c("hA", "hA"),
                       samesex = c(1, 0), treatment_3plus = c(1, 0), event_age = c(2, 4))
chk3 <- verify_key_coverage_and_recompute(metadata3, cached3)
stopifnot_msg(isTRUE(chk3$one_to_one_and_recomputed_equal),
              "key coverage: two distinct mothers, same HH, different Z/D values -- correctly passes (legitimate, not a conflict)")

## duplicate PERSON key in cache -- must FAIL even if values identical
cached_dup <- rbindlist(list(cached_match, cached_match[1]))
chk_dup <- verify_key_coverage_and_recompute(metadata2, cached_dup)
stopifnot_msg(chk_dup$n_duplicate_person_key_cache == 1 && !isTRUE(chk_dup$one_to_one_and_recomputed_equal),
              "key coverage: duplicate cache person_key fails even with identical values")

## missing key (in cache, not in metadata) -- must FAIL
cached_missing <- rbindlist(list(cached_match, data.table(person_key = "p9", household_key = "h9",
                                                            samesex = 1, treatment_3plus = 1, event_age = 0)))
chk_missing <- verify_key_coverage_and_recompute(metadata2, cached_missing)
stopifnot_msg(chk_missing$n_cached_keys_missing_from_metadata == 1 && !isTRUE(chk_missing$one_to_one_and_recomputed_equal),
              "key coverage: a cached person_key absent from metadata is caught, not silently dropped")

## recompute mismatch (event_age changed) -- must FAIL
cached_mismatch <- data.table::copy(cached_match); cached_mismatch[person_key == "p3", event_age := 5]
chk_mismatch <- verify_key_coverage_and_recompute(metadata2, cached_mismatch)
stopifnot_msg(chk_mismatch$n_recompute_mismatches == 1 && !isTRUE(chk_mismatch$one_to_one_and_recomputed_equal),
              "key coverage: recomputed event_age mismatch detected, not silently accepted")

## extra key (metadata has a person not in cache) -- must FAIL (strict: no
## legitimate reason for extra eligible mothers vs the original cache)
metadata_extra <- rbindlist(list(metadata2, data.table(person_key = "p9", household_key = "h9",
                                                         eligible = TRUE, samesex = 1, treatment_3plus = 1, event_age = 0)))
chk_extra <- verify_key_coverage_and_recompute(metadata_extra, cached_match)
stopifnot_msg(chk_extra$n_metadata_keys_not_in_cache == 1 && !isTRUE(chk_extra$one_to_one_and_recomputed_equal),
              "key coverage: an extra metadata person_key not in cache fails (strict, not just a warning)")

## ---- ambiguity_sensitivity_exclude_a2_eq_a3 -------------------------------
ss_fix <- data.table(a2 = c(3, 2, NA, 5), a3 = c(3, NA, 1, 4), samesex = c(1, 0, 1, 0))
sens <- ambiguity_sensitivity_exclude_a2_eq_a3(ss_fix)
stopifnot_msg(sens$n_excluded_a2_eq_a3 == 1, "ambiguity sensitivity: excludes exactly the 1 row with a2==a3")

## ---- sex_cell_counts: BB/BG/GB/GG ----------------------------------------
cellfix <- data.table(sex1 = c(1, 1, 2, 2, 1), sex2 = c(1, 2, 1, 2, 1),
                       treatment_3plus = c(1, 0, 1, 0, 0), mother_weight = c(1, 2, 3, 4, 5),
                       ROOMS_out = c(6, 5, 7, 4, 6))
cells <- sex_cell_counts(cellfix, outcomes = "ROOMS_out", weight_var = "mother_weight")
stopifnot_msg(all(c("BB", "BG", "GB", "GG") %in% cells$cell), "sex_cell_counts: all 4 ordered cells present")
stopifnot_msg(cells[cell == "BB", n_weighted] == 6, "sex_cell_counts: BB weighted N = 1+5 = 6")

## ---- rf_fs_diagnostic_fit: RF/FS only, no IV ever; unequal weights, shared HH, missing controls
set.seed(21)
n <- 1500
mat_age <- round(runif(n, 25, 40)); mat_age[sample(n, 30)] <- NA  # missing controls
hh <- paste0("hh", seq_len(n) %% (n / 2))  # shared HH across rows
w <- sample(1:9, n, replace = TRUE)        # unequal weights
Z <- rbinom(n, 1, 0.3)
D <- rbinom(n, 1, pmin(pmax(0.3 + 0.2 * Z, 0.01), 0.99))
Y <- 5 - 0.4 * D + rnorm(n)
dfit <- data.table(Y = Y, D = D, Z = Z, mat_age = mat_age, mother_weight = w, household_key = hh)
fit <- rf_fs_diagnostic_fit(dfit, "Y", "D", "Z", "mat_age", "mother_weight", "household_key")
stopifnot_msg(fit$status == "full_fit", "rf_fs_diagnostic_fit: reaches full_fit on unequal weights/shared HH/missing controls")
stopifnot_msg(is.null(fit[["iv_coef"]]) && is.null(fit[["ar_summary_lower"]]),
              "rf_fs_diagnostic_fit: NEVER returns an iv_coef or AR field -- IV/AR is not fit")
stopifnot_msg(!is.null(fit$rf_b) && !is.null(fit$rf_V), "rf_fs_diagnostic_fit: named b and full V returned")
stopifnot_msg(abs(sqrt(fit$rf_V[["Z"]][["Z"]]) - fit$rf_se) < 1e-9, "rf_fs_diagnostic_fit: V diagonal matches reported SE")
stopifnot_msg(fit$rf_nobs == fit$n_usable && isTRUE(fit$rf_nobs_matches_usable),
              "rf_fs_diagnostic_fit: nobs equals the ONE common complete-case sample size (controls included upfront, not left to feols to drop)")
stopifnot_msg(fit$n_usable < n, "rf_fs_diagnostic_fit: the common complete-case sample correctly excludes the missing-mat_age rows")

## ---- joint_bb_gg_rf_fs: four-cell toy DGP with DISTINCT BB/GG levels, catches wrong reference
set.seed(22)
n2 <- 2000
sex1 <- sample(1:2, n2, replace = TRUE); sex2 <- sample(1:2, n2, replace = TRUE)
bb <- as.numeric(sex1 == 1 & sex2 == 1); gg <- as.numeric(sex1 == 2 & sex2 == 2)
mat_age2 <- round(runif(n2, 25, 40))
D2 <- rbinom(n2, 1, 0.3)
Y2 <- 5 + 2.0 * bb - 1.0 * gg + rnorm(n2, sd = 0.5)  # BB and GG have DISTINCT, opposite-sign true effects vs mixed
dtoy <- data.table(ROOMS_out = Y2, treatment_3plus = D2, sex1 = sex1, sex2 = sex2, mat_age = mat_age2,
                    mother_weight = 1, household_key = paste0("h", seq_len(n2)))
jr <- joint_bb_gg_rf_fs(dtoy, "ROOMS_out", "mat_age", "mother_weight", "household_key")
stopifnot_msg(jr$rf$status == "full_fit", "joint_bb_gg_rf_fs: reaches full_fit")
stopifnot_msg(abs(jr$rf$bb_coef - 2.0) < 0.3, "joint_bb_gg_rf_fs: BB coefficient recovers its true +2.0 vs mixed reference")
stopifnot_msg(abs(jr$rf$gg_coef - (-1.0)) < 0.3, "joint_bb_gg_rf_fs: GG coefficient recovers its true -1.0 vs mixed reference (wrong reference would blend these)")
stopifnot_msg(jr$rf$bb_coef > 0 && jr$rf$gg_coef < 0,
              "joint_bb_gg_rf_fs: BB and GG have opposite signs as designed -- a BB-vs-(GG+mixed) reference would NOT show this contrast cleanly")

## ---- additive_sex_control_rf_fs: logs constant added control if collinear
dconst <- data.table::copy(dfit)
dconst[, sex1 := 1]; dconst[, sex2 := 1]  # both constant -> firstchildboy/secondchildboy constant
fit_const <- additive_sex_control_rf_fs(dconst, "Y", "mat_age", "mother_weight", "household_key")
stopifnot_msg(length(fit_const$constant_added_controls) == 2,
              "additive_sex_control_rf_fs: logs both added controls as constant when sex1/sex2 are constant")

## ---- event_age_case_list / event_age_summary_row: all 6 ages, no IV/AR ---
cases <- event_age_case_list(outcomes = "Y", ages = 0:5)
stopifnot_msg(length(cases) == 6, "event_age_case_list: all 6 ages (0:5) enumerated for 1 outcome")
ea <- sample(0:5, n, replace = TRUE)
dfit[, event_age := ea]
rows_ea <- lapply(cases, function(cs) {
  d <- dfit[event_age == cs$event_age]
  r <- rf_fs_diagnostic_fit(d, cs$outcome, "D", "Z", "mat_age", "mother_weight", "household_key")
  event_age_summary_row(r, cs$outcome, cs$event_age)
})
tab <- rbindlist(rows_ea)
stopifnot_msg(nrow(tab) == 6, "event-age diagnostic: all 6 ages (0:5) reported, none dropped")
stopifnot_msg(!("iv_coef" %in% names(tab)), "event-age diagnostic table has no iv_coef column (RF/FS only)")

## ---- extract_samesex_roster_metadata_one_state: runs on existing local fixture (SYNTHETIC)
fixture_path <- "/tmp/acs_iv_fixture/partitions/statefip_99/housing_narrow.rds"
if (file.exists(fixture_path)) {
  dt_state <- readRDS(fixture_path)
  meta_state <- extract_samesex_roster_metadata_one_state(dt_state)
  stopifnot_msg(nrow(meta_state) > 0, "metadata extraction: produces a non-empty metadata table on the local SYNTHETIC fixture")
  stopifnot_msg(all(c("sex1", "sex2", "a1", "a2", "a3", "linked_child_count", "NCHILD_norm",
                       "household_key", "person_key") %in% names(meta_state)),
                "metadata extraction: all required metadata columns present")
  audit_state <- sex_order_audit(meta_state)
  cat(sprintf("SYNTHETIC local fixture metadata audit (NOT a real-data finding): n=%d a2eqa3=%d\n",
              audit_state$n_rows, audit_state$n_third_child_age_tie_a2_eq_a3))
} else {
  cat("SKIP: local fixture RDS not found at", fixture_path, "\n")
}

cat("\nALL SAMESEX DIAGNOSIS TESTS PASSED\n")
