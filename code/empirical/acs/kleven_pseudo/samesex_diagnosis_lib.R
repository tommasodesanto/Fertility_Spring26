# SameSex2 negative-rooms diagnostic library. Investigates whether the
# national SameSex2 rooms RF (job 18274624: -0.0376 [-0.0476,-0.0276])
# reflects a coding/link/order problem, a direct sex-composition/room-
# sharing channel, or a robust-but-nonidentified association. Does NOT
# assume room-sharing is proven and does NOT search for a specification
# that flips the sign. All functions here operate on already-approved
# roster/instrument objects from twins_samesex_iv_lib.R -- no new sample,
# weight, clustering, or roster-construction rule is introduced.
suppressMessages({
  library(data.table)
  library(fixest)
  library(jsonlite)
})

## ---- A) Reproduction gate -------------------------------------------------
#' Refit SameSex2 pooled ROOMS RF/FS from the saved national analytic
#' frames and compare to the saved fit_receipts JSON within `tolerance`.
#' Must be run before any new diagnostic is trusted. Fails loudly (returns
#' pass=FALSE with the exact diffs) rather than silently proceeding.
verify_reproduction <- function(analytic_frames_path, identity_path, saved_receipt_path,
                                 tolerance = 1e-8) {
  af <- readRDS(analytic_frames_path)
  identity <- jsonlite::fromJSON(identity_path)
  saved <- jsonlite::fromJSON(saved_receipt_path)
  fit <- fit_instrument_outcome(
    af$ss, outcome = "ROOMS_out", treatment = "treatment_3plus", instrument = "samesex",
    controls_fml = identity$controls_fml, weight_var = identity$weight_var,
    cluster_var = identity$cluster_var, ar_grid = NULL
  )
  cmp <- function(a, b) isTRUE(abs(a - b) <= tolerance)
  checks <- list(
    rf_coef = cmp(fit$rf_coef, saved$rf_coef), rf_se = cmp(fit$rf_se, saved$rf_se),
    fs_coef = cmp(fit$fs_coef, saved$fs_coef), fs_se = cmp(fit$fs_se, saved$fs_se),
    rf_nobs = identical(fit$rf_nobs_fit, saved$rf_nobs_fit),
    fs_nobs = identical(fit$fs_nobs_fit, saved$fs_nobs_fit)
  )
  list(pass = all(unlist(checks)), checks = checks, fit = fit, saved = saved)
}

## ---- B) Event-age-by-age RF/FS (all ages kept, never dropped) ------------
#' event0 is post-second-birth, NOT a pre-treatment/placebo test: a third
#' child can already exist by event age 0 and interview timing/selection
#' matters. This function reports every age 0:5 for every outcome -- it
#' never selects ages by significance or drops an unfavorable age.
event_age_rf_fs_diagnostic <- function(ss, outcomes, controls_fml, weight_var, cluster_var,
                                        ages = 0:5) {
  out <- list()
  for (oc in outcomes) {
    for (ea in ages) {
      d <- ss[event_age == ea]
      r <- tryCatch(fit_instrument_outcome(d, oc, "treatment_3plus", "samesex",
                      controls_fml, weight_var, cluster_var, ar_grid = NULL),
                    error = function(e) list(status = "error", message = conditionMessage(e)))
      r$outcome <- oc; r$event_age <- ea
      out[[length(out) + 1]] <- r
    }
  }
  data.table::rbindlist(lapply(out, function(r) data.table(
    outcome = r$outcome, event_age = r$event_age, status = r$status %||% NA_character_,
    n_usable = r$n_usable %||% NA_integer_, n_households = r$n_households %||% NA_integer_,
    n_instrument_positive = r$n_instrument_positive %||% NA_integer_,
    rf_coef = r$rf_coef %||% NA_real_, rf_se = r$rf_se %||% NA_real_,
    rf_ci_lower = r$rf_ci_lower %||% NA_real_, rf_ci_upper = r$rf_ci_upper %||% NA_real_,
    fs_coef = r$fs_coef %||% NA_real_, fs_se = r$fs_se %||% NA_real_,
    first_stage_F = r$first_stage_F %||% NA_real_
  )), fill = TRUE)
}

## ---- C) Per-state roster metadata extraction (memory-bounded) ------------
#' Mirrors the production driver's exact per-state pipeline (source/sample
#' gate -> roster -> outcomes) but retains sex1/sex2/a1/a2/a3/
#' linked_child_count/NCHILD_norm/person_key/household_key instead of
#' projecting to the slim regression frame, restricted to the SameSex2
#' eligible_pool. Processes ONE state at a time; the caller is responsible
#' for rbind-ing the small metadata (not the wide per-state table) across
#' states and for releasing each wide state object immediately (matching
#' the production driver's memory-bounded pattern; see
#' run_twins_samesex_iv.R's per-state loop).
extract_samesex_roster_metadata_one_state <- function(dt_state, age_lo = 21, age_hi = 35) {
  gate <- apply_source_sample_gate(dt_state)
  built <- build_mother_roster(gate$data)
  mr <- add_outcomes(built$mother_rows)
  minor_gate <- apply_oldest_child_minor_gate(mr)
  mr2 <- minor_gate$data
  mr2[, weight_positive := PERWT > 0 & is.finite(PERWT)]
  mr2 <- mr2[weight_positive == TRUE & AGE_norm >= age_lo & AGE_norm <= age_hi]
  ss <- build_samesex2(mr2, age_grid = 0:5)
  ss[eligible_pool == TRUE, .(
    person_key, household_key, sex1, sex2, a1, a2, a3,
    linked_child_count, NCHILD_norm, samesex, primary_age_tie, eligible,
    treatment_3plus, event_age, PERWT
  )]
}

#' One-to-one key coverage and Z/D/event_age recomputation check between
#' freshly extracted metadata and the cached analytic-frame ss (household
#' key + samesex + treatment_3plus + event_age must match exactly). Fails
#' rather than silently accepting a mismatched sample.
verify_key_coverage_and_recompute <- function(metadata, cached_ss) {
  cached_primary <- cached_ss[, .(household_key, samesex_cached = samesex,
                                   treatment_3plus_cached = treatment_3plus,
                                   event_age_cached = event_age)]
  m <- merge(metadata[eligible == TRUE], cached_primary, by = "household_key", all = TRUE)
  n_cached <- nrow(cached_primary)
  n_metadata_eligible <- sum(metadata$eligible)
  n_matched <- sum(!is.na(m$samesex) & !is.na(m$samesex_cached))
  n_cached_unmatched <- sum(is.na(m$samesex))
  n_metadata_unmatched <- sum(is.na(m$samesex_cached))
  mismatched <- m[!is.na(samesex) & !is.na(samesex_cached) &
                     (samesex != samesex_cached | treatment_3plus != treatment_3plus_cached |
                        event_age != event_age_cached)]
  list(
    n_cached = n_cached, n_metadata_eligible = n_metadata_eligible, n_matched = n_matched,
    n_cached_unmatched = n_cached_unmatched, n_metadata_unmatched = n_metadata_unmatched,
    n_recompute_mismatches = nrow(mismatched),
    one_to_one_and_recomputed_equal = (n_cached_unmatched == 0 && n_metadata_unmatched == 0 &&
                                          nrow(mismatched) == 0)
  )
}

## ---- Sex/order/link audits (report, never silently filter) ---------------
#' Valid sex values, descending age order, and third-child age-tie
#' ambiguity (a2==a3, distinct from the already-excluded a1==a2 primary
#' tie). Also NCHILD-vs-linked-count and duplicate-key checks.
sex_order_audit <- function(metadata) {
  list(
    n_rows = nrow(metadata),
    n_sex1_invalid = sum(!(metadata$sex1 %in% c(1, 2)), na.rm = TRUE),
    n_sex2_invalid = sum(!(metadata$sex2 %in% c(1, 2)), na.rm = TRUE),
    n_age_order_violation_a1_lt_a2 = sum(!is.na(metadata$a1) & !is.na(metadata$a2) &
                                            metadata$a1 < metadata$a2, na.rm = TRUE),
    n_age_order_violation_a2_lt_a3 = sum(!is.na(metadata$a2) & !is.na(metadata$a3) &
                                            metadata$a2 < metadata$a3, na.rm = TRUE),
    n_third_child_age_tie_a2_eq_a3 = sum(!is.na(metadata$a2) & !is.na(metadata$a3) &
                                            metadata$a2 == metadata$a3, na.rm = TRUE),
    n_nchild_link_mismatch = sum(!is.na(metadata$NCHILD_norm) &
                                    metadata$NCHILD_norm != metadata$linked_child_count, na.rm = TRUE),
    n_duplicate_household_key = sum(duplicated(metadata$household_key)),
    n_duplicate_person_key = sum(duplicated(metadata$person_key))
  )
}

## ---- D) Ordered-sex cells and first/second-child-sex controls ------------
#' BB/BG/GB/GG cells from sex1 (oldest) x sex2 (second-oldest); sex==1 is
#' male, sex==2 is female (matches build_mother_roster's female_code=2
#' convention used throughout). Reports weighted/unweighted N and D/Y
#' means per cell -- descriptive, not a fit.
sex_cell_counts <- function(ss_with_sex, outcomes, weight_var = "mother_weight") {
  d <- data.table::copy(ss_with_sex)
  d[, cell := data.table::fcase(
    sex1 == 1 & sex2 == 1, "BB", sex1 == 1 & sex2 == 2, "BG",
    sex1 == 2 & sex2 == 1, "GB", sex1 == 2 & sex2 == 2, "GG",
    default = NA_character_
  )]
  agg <- d[!is.na(cell), .(
    n_unweighted = .N,
    n_weighted = sum(get(weight_var)),
    D_mean_unweighted = mean(treatment_3plus),
    D_mean_weighted = weighted.mean(treatment_3plus, get(weight_var))
  ), by = cell]
  for (oc in outcomes) {
    agg[[paste0(oc, "_Y_mean_unweighted")]] <- vapply(agg$cell, function(cc)
      mean(d[cell == cc][[oc]], na.rm = TRUE), numeric(1))
    agg[[paste0(oc, "_Y_mean_weighted")]] <- vapply(agg$cell, function(cc) {
      sub <- d[cell == cc]
      ok <- !is.na(sub[[oc]])
      if (!any(ok)) return(NA_real_)
      stats::weighted.mean(sub[[oc]][ok], sub[[weight_var]][ok])
    }, numeric(1))
  }
  agg
}

#' Adds firstchildboy/secondchildboy indicators and a controls_fml variant
#' that includes them, for the "same-sex RF/FS controlling for each
#' child's sex separately" diagnostic (D).
add_first_second_child_sex_controls <- function(ss) {
  d <- data.table::copy(ss)
  d[, firstchildboy := as.numeric(sex1 == 1)]
  d[, secondchildboy := as.numeric(sex2 == 1)]
  d
}

#' RF/FS for a binary group indicator (e.g. both-boys vs everyone else,
#' both-girls vs everyone else) with the SAME baseline controls. This is a
#' separate, clearly-labeled diagnostic -- not a replacement for the
#' primary same-sex instrument.
group_vs_mixed_rf_fs <- function(ss, group_label, group_indicator_col, outcome,
                                  controls_fml, weight_var, cluster_var) {
  fit_instrument_outcome(ss, outcome, "treatment_3plus", group_indicator_col,
                          controls_fml, weight_var, cluster_var, ar_grid = NULL)
}

## ---- E) Ambiguity sensitivity: exclude a2==a3 -----------------------------
#' Labeled sensitivity excluding third-child age ties among the baseline
#' SameSex2 sample. Reports the exact N excluded; never silently changes
#' the primary sample.
ambiguity_sensitivity_exclude_a2_eq_a3 <- function(ss_with_a3) {
  n_before <- nrow(ss_with_a3)
  keep <- is.na(ss_with_a3$a2) | is.na(ss_with_a3$a3) | ss_with_a3$a2 != ss_with_a3$a3
  list(data = ss_with_a3[keep], n_before = n_before, n_excluded_a2_eq_a3 = sum(!keep))
}
