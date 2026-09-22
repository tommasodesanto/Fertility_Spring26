# SameSex2 negative-rooms diagnostic library. Investigates whether the
# national SameSex2 rooms RF (job 18274624: -0.0376 rooms [-0.0476,-0.0276])
# reflects a coding/link/order problem, a direct sex-composition/room-
# sharing channel, or a robust-but-nonidentified association. Does NOT
# assume room-sharing is proven and does NOT search for a specification
# that flips the sign. All functions here operate on already-approved
# roster/instrument objects from twins_samesex_iv_lib.R -- no new sample,
# weight, clustering, or roster-construction rule is introduced.
#
# The additive first/second-child-sex control diagnostic (D) follows the
# convention in Angrist & Evans (1998, AER, "Children and Their Parents'
# Labor Supply"), Table 7/8 notes (pp. 465-466), which include Boy1st/
# Boy2nd alongside the same-sex instrument -- our ACS sample, controls,
# and outcome are NOT an exact replication of that design, only the same
# additive-sex-control idea.
# https://www.dpipe.tsukuba.ac.jp/~naito/teaching_web/public_economics_2013_web/Angrist_Evans_AER_1998.pdf
# MOMLOC links are social/coresident, confirmed by IPUMS documentation:
# https://usa.ipums.org/usa-action/variables/MOMLOC
suppressMessages({
  library(data.table)
  library(fixest)
  library(jsonlite)
})

MAX_STAGE_REGRESSIONS <- 100L

## ---- shared: build ONE explicit common complete-case sample --------------
#' Common sample over outcome, treatment, instrument, weight (positive,
#' finite), cluster (non-missing), and every control variable (finite/
#' non-missing) -- fit on these identical rows so actual post-fit nobs
#' equals nrow(the common sample) by construction, not by feols's own
#' listwise deletion after the fact.
build_common_complete_case <- function(data, needed_numeric_cols, weight_var, cluster_var) {
  missing_cols <- setdiff(c(needed_numeric_cols, weight_var, cluster_var), names(data))
  if (length(missing_cols) > 0) {
    return(list(ok = FALSE, error = sprintf("missing columns: %s", paste(missing_cols, collapse = ", "))))
  }
  d <- data.table::copy(data)
  ok_mask <- rep(TRUE, nrow(d))
  for (v in needed_numeric_cols) {
    x <- suppressWarnings(as.numeric(d[[v]]))
    ok_mask <- ok_mask & !is.na(x) & is.finite(x)
  }
  w <- suppressWarnings(as.numeric(d[[weight_var]]))
  ok_mask <- ok_mask & !is.na(w) & is.finite(w) & w > 0
  ok_mask <- ok_mask & !is.na(d[[cluster_var]])
  list(ok = TRUE, data = d[ok_mask])
}

## ---- A) Reproduction gate: exact path + identity BEFORE read -------------
#' Verifies (1) the analytic-frames path EXACTLY equals the canonical
#' `file.path(national_outdir, "analytic_frames.rds")`, (2) file size/MD5
#' match the pinned identity, (3) the library file's own MD5 matches the
#' identity's pinned lib_md5 (code identity, not raw source data), all
#' BEFORE reading the RDS; then (4) the saved receipt's outcome/
#' instrument/treatment match expected; THEN refits RF/FS on the ONE
#' common complete-case sample and compares to the saved receipt within
#' `tolerance` for coef/SE, but EXACT numeric equality for nobs. Fails
#' loudly at the first failed gate, before any scientific diagnostic runs.
verify_reproduction <- function(national_outdir, identity, lib_path, saved_receipt_path,
                                 expected_outcome = "ROOMS_out", expected_instrument = "samesex",
                                 expected_treatment = "treatment_3plus", tolerance = 1e-8) {
  gate <- list()
  canonical_path <- file.path(national_outdir, "analytic_frames.rds")
  gate$path_matches_canonical <- isTRUE(identical(normalizePath(identity$analytic_frames_path, mustWork = FALSE),
                                                   normalizePath(canonical_path, mustWork = FALSE)))
  if (!gate$path_matches_canonical) return(list(pass = FALSE, gate = gate, stage = "path_not_canonical"))

  gate$path_exists <- file.exists(canonical_path)
  if (!gate$path_exists) return(list(pass = FALSE, gate = gate, stage = "path_missing"))

  actual_size <- file.info(canonical_path)$size
  gate$size_match <- isTRUE(actual_size == identity$analytic_frames_size_bytes)
  if (!gate$size_match) return(list(pass = FALSE, gate = gate, stage = "size_mismatch", actual_size = actual_size))

  actual_md5 <- unname(tools::md5sum(canonical_path))
  gate$md5_match <- isTRUE(actual_md5 == identity$analytic_frames_md5)
  if (!gate$md5_match) return(list(pass = FALSE, gate = gate, stage = "md5_mismatch", actual_md5 = actual_md5))

  lib_md5_actual <- unname(tools::md5sum(lib_path))
  gate$lib_md5_match <- isTRUE(lib_md5_actual == identity$lib_md5)
  if (!gate$lib_md5_match) return(list(pass = FALSE, gate = gate, stage = "lib_md5_mismatch", actual = lib_md5_actual))

  saved <- jsonlite::fromJSON(saved_receipt_path)
  identity_ok <- identical(saved$outcome, expected_outcome) &&
    identical(saved$instrument, expected_instrument) && identical(saved$treatment, expected_treatment)
  gate$saved_receipt_identity_match <- identity_ok
  if (!identity_ok) {
    return(list(pass = FALSE, gate = gate, stage = "saved_receipt_identity_mismatch",
                saved_outcome = saved$outcome, saved_instrument = saved$instrument, saved_treatment = saved$treatment))
  }

  af <- readRDS(canonical_path)
  fit <- rf_fs_diagnostic_fit(af$ss, outcome = expected_outcome, treatment = expected_treatment,
                               instrument = expected_instrument, controls_fml = identity$controls_fml,
                               weight_var = identity$weight_var, cluster_var = identity$cluster_var)
  cmp_tol <- function(a, b) isTRUE(abs(as.numeric(a) - as.numeric(b)) <= tolerance)
  cmp_exact <- function(a, b) isTRUE(as.numeric(a) == as.numeric(b))
  checks <- list(
    rf_coef = cmp_tol(fit$rf_coef, saved$rf_coef), rf_se = cmp_tol(fit$rf_se, saved$rf_se),
    fs_coef = cmp_tol(fit$fs_coef, saved$fs_coef), fs_se = cmp_tol(fit$fs_se, saved$fs_se),
    rf_nobs_exact = cmp_exact(fit$rf_nobs, saved$rf_nobs_fit),
    fs_nobs_exact = cmp_exact(fit$fs_nobs, saved$fs_nobs_fit)
  )
  list(pass = all(unlist(checks)), gate = gate, stage = "numeric_comparison", checks = checks,
       fit = fit, saved = saved, af = af)
}

## ---- 3) RF/FS-only diagnostic helper: NEVER fits IV/AR --------------------
#' Diagnostics-only fitter on ONE common complete-case sample (see
#' build_common_complete_case): RF (Y~Z+controls) and FS (D~Z+controls)
#' only. Named b, full clustered V, nobs (== n_usable by construction),
#' HH count, Z-positive count computed FROM the common sample, SE, normal
#' 95% CI, single-instrument cluster-robust Wald F. Explicit errors and
#' captured warnings returned, never silently swallowed. No positive- or
#' any- sign gate.
rf_fs_diagnostic_fit <- function(data, outcome, treatment, instrument, controls_fml,
                                  weight_var, cluster_var, z = stats::qnorm(0.975)) {
  ctrl_vars <- if (nzchar(controls_fml)) all.vars(stats::as.formula(paste("~", controls_fml))) else character(0)
  needed_numeric <- unique(c(outcome, treatment, instrument, ctrl_vars))
  cc <- build_common_complete_case(data, needed_numeric, weight_var, cluster_var)
  if (!cc$ok) return(list(status = "error", error = cc$error))
  usable <- cc$data
  usable[[instrument]] <- as.numeric(usable[[instrument]])
  usable[[treatment]] <- as.numeric(usable[[treatment]])
  n_usable <- nrow(usable)
  n_hh <- length(unique(usable[[cluster_var]]))
  n_zpos <- sum(usable[[instrument]] == 1, na.rm = TRUE)
  if (n_usable < 10) {
    return(list(status = "insufficient_support", n_usable = n_usable, n_households = n_hh, n_instrument_positive = n_zpos))
  }
  rhs <- if (nzchar(controls_fml)) paste0(instrument, " + ", controls_fml) else instrument
  w_fml <- stats::as.formula(paste0("~", weight_var)); c_fml <- stats::as.formula(paste0("~", cluster_var))
  rf_fml <- stats::as.formula(paste0(outcome, " ~ ", rhs)); fs_fml <- stats::as.formula(paste0(treatment, " ~ ", rhs))

  warns <- character(0)
  wh <- withCallingHandlers(
    list(
      rf_fit = tryCatch(fixest::feols(rf_fml, data = usable, weights = w_fml, cluster = c_fml, notes = FALSE), error = function(e) e),
      fs_fit = tryCatch(fixest::feols(fs_fml, data = usable, weights = w_fml, cluster = c_fml, notes = FALSE), error = function(e) e)
    ),
    warning = function(w) { warns[[length(warns) + 1]] <<- conditionMessage(w); invokeRestart("muffleWarning") }
  )
  rf_fit <- wh$rf_fit; fs_fit <- wh$fs_fit
  rf_err <- if (inherits(rf_fit, "error")) conditionMessage(rf_fit) else NA_character_
  fs_err <- if (inherits(fs_fit, "error")) conditionMessage(fs_fit) else NA_character_
  if (!is.na(rf_err)) rf_fit <- NULL
  if (!is.na(fs_err)) fs_fit <- NULL

  extract <- function(fit, coefname) {
    if (is.null(fit) || !(coefname %in% names(stats::coef(fit))))
      return(list(coef = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                  nobs = NA_integer_, b = NULL, V = NULL, F = NA_real_, nobs_matches_usable = NA))
    b <- stats::coef(fit); V <- as.matrix(stats::vcov(fit))
    coefv <- unname(b[coefname]); se <- unname(sqrt(diag(V))[coefname])
    wt <- tryCatch(fixest::wald(fit, coefname, print = FALSE), error = function(e) NULL)
    nobs_fit <- stats::nobs(fit)
    list(coef = coefv, se = se, ci_lower = coefv - z * se, ci_upper = coefv + z * se,
         nobs = nobs_fit, b = as.list(b), V = apply(V, 1, as.list),
         F = if (!is.null(wt)) unname(wt$stat) else NA_real_,
         nobs_matches_usable = isTRUE(nobs_fit == n_usable))
  }
  rfx <- extract(rf_fit, instrument); fsx <- extract(fs_fit, instrument)
  list(
    status = if (is.na(rf_err) && is.na(fs_err) && !is.na(rfx$coef) && !is.na(fsx$coef)) "full_fit" else "error",
    n_usable = n_usable, n_households = n_hh, n_instrument_positive = n_zpos,
    rf_coef = rfx$coef, rf_se = rfx$se, rf_ci_lower = rfx$ci_lower, rf_ci_upper = rfx$ci_upper,
    rf_nobs = rfx$nobs, rf_nobs_matches_usable = rfx$nobs_matches_usable,
    rf_b = rfx$b, rf_V = rfx$V, rf_error = rf_err,
    fs_coef = fsx$coef, fs_se = fsx$se, fs_ci_lower = fsx$ci_lower, fs_ci_upper = fsx$ci_upper,
    fs_nobs = fsx$nobs, fs_nobs_matches_usable = fsx$nobs_matches_usable,
    fs_b = fsx$b, fs_V = fsx$V, fs_error = fs_err, first_stage_F = fsx$F,
    formula = list(rf = deparse(rf_fml), fs = deparse(fs_fml)),
    weight_var = weight_var, cluster_var = cluster_var, warnings = warns
  )
}

## ---- B) Event-age-by-age RF/FS: driver persists EACH case atomically -----
#' Returns the exact ordered list of (outcome, event_age) cases to fit;
#' the driver calls rf_fs_diagnostic_fit() once per case and persists the
#' FULL result (named b/V/nobs/warnings included) before moving to the
#' next case -- this function itself only enumerates cases and builds the
#' compact summary table from already-fit results, it does not collapse
#' or discard the per-case detail.
event_age_case_list <- function(outcomes, ages = 0:5) {
  cases <- list()
  for (oc in outcomes) for (ea in ages) cases[[length(cases) + 1]] <- list(outcome = oc, event_age = ea)
  cases
}

event_age_summary_row <- function(r, oc, ea) {
  data.table(
    outcome = oc, event_age = ea, status = r$status %||% NA_character_,
    n_usable = r$n_usable %||% NA_integer_, n_households = r$n_households %||% NA_integer_,
    n_instrument_positive = r$n_instrument_positive %||% NA_integer_,
    rf_coef = r$rf_coef %||% NA_real_, rf_se = r$rf_se %||% NA_real_,
    rf_ci_lower = r$rf_ci_lower %||% NA_real_, rf_ci_upper = r$rf_ci_upper %||% NA_real_,
    fs_coef = r$fs_coef %||% NA_real_, fs_se = r$fs_se %||% NA_real_, first_stage_F = r$first_stage_F %||% NA_real_,
    n_warnings = length(r$warnings %||% character(0)), rf_error = r$rf_error %||% NA_character_,
    fs_error = r$fs_error %||% NA_character_
  )
}

## ---- C) Per-state roster metadata extraction (memory-bounded) ------------
#' Mirrors the production driver's exact per-state pipeline. Retains
#' a1,a2,a3,sex1,sex2,linked_child_count,NCHILD_norm,person_key,
#' household_key,samesex,primary_age_tie,eligible,treatment_3plus,
#' event_age,PERWT for the SameSex2 eligible_pool. a2 must equal
#' event_age by SameSex2's own construction (event_age <- a2); this is
#' asserted, not merely assumed.
extract_samesex_roster_metadata_one_state <- function(dt_state, age_lo = 21, age_hi = 35) {
  gate <- apply_source_sample_gate(dt_state)
  built <- build_mother_roster(gate$data)
  mr <- add_outcomes(built$mother_rows)
  minor_gate <- apply_oldest_child_minor_gate(mr)
  mr2 <- minor_gate$data
  mr2[, weight_positive := PERWT > 0 & is.finite(PERWT)]
  mr2 <- mr2[weight_positive == TRUE & AGE_norm >= age_lo & AGE_norm <= age_hi]
  ss <- build_samesex2(mr2, age_grid = 0:5)
  meta <- ss[eligible_pool == TRUE, .(
    person_key, household_key, sex1, sex2, a1, a2, a3,
    linked_child_count, NCHILD_norm, samesex, primary_age_tie, eligible,
    treatment_3plus, event_age, PERWT
  )]
  bad_a2 <- sum(!is.na(meta$a2) & !is.na(meta$event_age) & meta$a2 != meta$event_age)
  if (bad_a2 > 0) stop(sprintf("extract_samesex_roster_metadata_one_state: %d rows have a2 != event_age -- SameSex2 construction invariant violated", bad_a2))
  meta
}

#' STRICT exact mother-level (PERSON_KEY) coverage: metadata (already
#' eligible==TRUE) and cached_ss must be in EXACT 1:1 correspondence --
#' no missing keys, no extra keys (there is no legitimate reason for
#' eligible-mother counts to differ from the original cache), no
#' duplicate/missing keys on either side, and post-join household_key/
#' samesex/treatment_3plus/event_age must match exactly on every row.
verify_key_coverage_and_recompute <- function(metadata, cached_ss) {
  meta_e <- metadata[eligible == TRUE]
  dup_meta <- sum(duplicated(meta_e$person_key)); dup_cache <- sum(duplicated(cached_ss$person_key))
  na_meta <- sum(is.na(meta_e$person_key)); na_cache <- sum(is.na(cached_ss$person_key))
  keys_meta <- unique(meta_e$person_key); keys_cache <- unique(cached_ss$person_key)
  missing_from_meta <- setdiff(keys_cache, keys_meta)
  extra_in_meta <- setdiff(keys_meta, keys_cache)

  keys_ok <- (dup_meta == 0 && dup_cache == 0 && na_meta == 0 && na_cache == 0 &&
                length(missing_from_meta) == 0 && length(extra_in_meta) == 0 &&
                nrow(meta_e) == nrow(cached_ss))

  m <- merge(meta_e[, .(person_key, household_key, samesex, treatment_3plus, event_age)],
             cached_ss[, .(person_key, household_key_cached = household_key,
                            samesex_cached = samesex, treatment_3plus_cached = treatment_3plus,
                            event_age_cached = event_age)],
             by = "person_key")
  mismatched <- m[household_key != household_key_cached | samesex != samesex_cached |
                    treatment_3plus != treatment_3plus_cached | event_age != event_age_cached]
  list(
    n_cached = nrow(cached_ss), n_metadata_eligible = nrow(meta_e),
    n_duplicate_person_key_metadata = dup_meta, n_duplicate_person_key_cache = dup_cache,
    n_missing_person_key_metadata = na_meta, n_missing_person_key_cache = na_cache,
    n_cached_keys_missing_from_metadata = length(missing_from_meta),
    n_metadata_keys_not_in_cache = length(extra_in_meta),
    n_matched_join_rows = nrow(m), n_recompute_mismatches = nrow(mismatched), keys_ok = keys_ok,
    one_to_one_and_recomputed_equal = (keys_ok && nrow(m) == nrow(cached_ss) && nrow(mismatched) == 0)
  )
}

## ---- Sex/order/link audits (report, never silently filter) ---------------
sex_order_audit <- function(metadata) {
  list(
    n_rows = nrow(metadata),
    n_sex1_invalid = sum(!(metadata$sex1 %in% c(1, 2)), na.rm = TRUE),
    n_sex2_invalid = sum(!(metadata$sex2 %in% c(1, 2)), na.rm = TRUE),
    n_age_order_violation_a1_lt_a2 = sum(!is.na(metadata$a1) & !is.na(metadata$a2) & metadata$a1 < metadata$a2, na.rm = TRUE),
    n_age_order_violation_a2_lt_a3 = sum(!is.na(metadata$a2) & !is.na(metadata$a3) & metadata$a2 < metadata$a3, na.rm = TRUE),
    n_third_child_age_tie_a2_eq_a3 = sum(!is.na(metadata$a2) & !is.na(metadata$a3) & metadata$a2 == metadata$a3, na.rm = TRUE),
    n_nchild_link_mismatch = sum(!is.na(metadata$NCHILD_norm) & metadata$NCHILD_norm != metadata$linked_child_count, na.rm = TRUE),
    n_duplicate_person_key = sum(duplicated(metadata$person_key)),
    n_unique_household_key = length(unique(metadata$household_key)),
    n_households_with_multiple_mothers = { tab <- table(metadata$household_key); sum(tab > 1) }
  )
}

## ---- D) Ordered-sex cells, additive controls, and JOINT BB/GG dummies ----
sex_cell_counts <- function(ss_with_sex, outcomes, weight_var = "mother_weight") {
  d <- data.table::copy(ss_with_sex)
  d[, cell := data.table::fcase(
    sex1 == 1 & sex2 == 1, "BB", sex1 == 1 & sex2 == 2, "BG",
    sex1 == 2 & sex2 == 1, "GB", sex1 == 2 & sex2 == 2, "GG", default = NA_character_
  )]
  agg <- d[!is.na(cell), .(
    n_unweighted = .N, n_weighted = sum(get(weight_var)),
    D_mean_unweighted = mean(treatment_3plus), D_mean_weighted = weighted.mean(treatment_3plus, get(weight_var))
  ), by = cell]
  for (oc in outcomes) {
    agg[[paste0(oc, "_Y_mean_unweighted")]] <- vapply(agg$cell, function(cc) mean(d[cell == cc][[oc]], na.rm = TRUE), numeric(1))
    agg[[paste0(oc, "_Y_mean_weighted")]] <- vapply(agg$cell, function(cc) {
      sub <- d[cell == cc]; ok <- !is.na(sub[[oc]])
      if (!any(ok)) return(NA_real_)
      stats::weighted.mean(sub[[oc]][ok], sub[[weight_var]][ok])
    }, numeric(1))
  }
  agg
}

add_first_second_child_sex_controls <- function(ss) {
  d <- data.table::copy(ss)
  d[, firstchildboy := as.numeric(sex1 == 1)]
  d[, secondchildboy := as.numeric(sex2 == 1)]
  d
}

#' RF for same-sex + additive firstchildboy/secondchildboy controls.
#' Uses rf_fs_diagnostic_fit (common complete-case sample); logs any
#' added control that is constant on that exact common sample.
additive_sex_control_rf_fs <- function(ss, outcome, controls_fml, weight_var, cluster_var) {
  d <- add_first_second_child_sex_controls(ss)
  extended_controls <- paste0(controls_fml, " + firstchildboy + secondchildboy")
  r <- rf_fs_diagnostic_fit(d, outcome, "treatment_3plus", "samesex", extended_controls, weight_var, cluster_var)
  cc <- build_common_complete_case(d, all.vars(stats::as.formula(paste("~", extended_controls))), weight_var, cluster_var)
  const_flags <- if (isTRUE(cc$ok)) vapply(c("firstchildboy", "secondchildboy"), function(v)
    length(unique(cc$data[[v]])) <= 1, logical(1)) else c(firstchildboy = NA, secondchildboy = NA)
  r$constant_added_controls <- names(const_flags)[!is.na(const_flags) & const_flags]
  r
}

#' JOINT both-boys/both-girls RF and FS in ONE regression per outcome,
#' mixed (BG+GB) as the implicit omitted reference. Built on the SAME
#' common-complete-case machinery (outcome/treatment/bothboys/bothgirls/
#' controls/weight/cluster all required non-missing+finite). Captures
#' warnings/errors; returns an error status (never aborts the whole
#' packet) if either fit fails; asserts BB/GG coefficients are
#' present/finite before reporting full_fit.
joint_bb_gg_rf_fs <- function(ss, outcome, controls_fml, weight_var, cluster_var) {
  d <- data.table::copy(ss)
  d[, bothboys := as.numeric(sex1 == 1 & sex2 == 1)]
  d[, bothgirls := as.numeric(sex1 == 2 & sex2 == 2)]
  ctrl_vars <- if (nzchar(controls_fml)) all.vars(stats::as.formula(paste("~", controls_fml))) else character(0)
  needed_numeric <- unique(c(outcome, "treatment_3plus", "bothboys", "bothgirls", ctrl_vars))
  cc <- build_common_complete_case(d, needed_numeric, weight_var, cluster_var)
  if (!cc$ok) return(list(outcome = outcome, status = "error", error = cc$error))
  usable <- cc$data
  n_usable <- nrow(usable); n_hh <- length(unique(usable[[cluster_var]]))
  rhs <- paste0("bothboys + bothgirls", if (nzchar(controls_fml)) paste0(" + ", controls_fml) else "")
  w_fml <- stats::as.formula(paste0("~", weight_var)); c_fml <- stats::as.formula(paste0("~", cluster_var))
  warns <- character(0)
  fit_one <- function(lhs) {
    fml <- stats::as.formula(paste0(lhs, " ~ ", rhs))
    fit <- withCallingHandlers(
      tryCatch(fixest::feols(fml, data = usable, weights = w_fml, cluster = c_fml, notes = FALSE), error = function(e) e),
      warning = function(w) { warns[[length(warns) + 1]] <<- conditionMessage(w); invokeRestart("muffleWarning") }
    )
    if (inherits(fit, "error")) return(list(status = "error", error = conditionMessage(fit), formula = deparse(fml)))
    b <- stats::coef(fit); V <- as.matrix(stats::vcov(fit))
    if (!all(c("bothboys", "bothgirls") %in% names(b)) || anyNA(b[c("bothboys", "bothgirls")])) {
      return(list(status = "error", error = "bothboys/bothgirls coefficient missing or NA", formula = deparse(fml)))
    }
    list(status = "full_fit", nobs = stats::nobs(fit), nobs_matches_usable = isTRUE(stats::nobs(fit) == n_usable),
         bb_coef = unname(b["bothboys"]), bb_se = unname(sqrt(diag(V))["bothboys"]),
         gg_coef = unname(b["bothgirls"]), gg_se = unname(sqrt(diag(V))["bothgirls"]),
         bb_gg_cov = unname(V["bothboys", "bothgirls"]), b = as.list(b), V = apply(V, 1, as.list),
         formula = deparse(fml))
  }
  rf <- fit_one(outcome); fs <- fit_one("treatment_3plus")
  list(outcome = outcome, status = if (identical(rf$status, "full_fit") && identical(fs$status, "full_fit")) "full_fit" else "partial_or_error",
       n_usable = n_usable, n_households = n_hh, n_bothboys = sum(usable$bothboys), n_bothgirls = sum(usable$bothgirls),
       n_mixed_reference = n_usable - sum(usable$bothboys) - sum(usable$bothgirls),
       weight_var = weight_var, cluster_var = cluster_var, rf = rf, fs = fs, warnings = warns)
}

## ---- E) Ambiguity sensitivity: exclude a2==a3 -----------------------------
#' Requires `a2` and `a3` present in ss_with_a3; asserts before use rather
#' than failing with an opaque NULL-column error. Labeled sensitivity;
#' reports the exact N excluded; never silently changes the primary sample.
ambiguity_sensitivity_exclude_a2_eq_a3 <- function(ss_with_a3) {
  required <- c("a2", "a3")
  missing_cols <- setdiff(required, names(ss_with_a3))
  if (length(missing_cols) > 0) {
    stop(sprintf("ambiguity_sensitivity_exclude_a2_eq_a3: missing required column(s): %s",
                  paste(missing_cols, collapse = ", ")))
  }
  n_before <- nrow(ss_with_a3)
  keep <- is.na(ss_with_a3$a2) | is.na(ss_with_a3$a3) | ss_with_a3$a2 != ss_with_a3$a3
  list(data = ss_with_a3[keep], n_before = n_before, n_excluded_a2_eq_a3 = sum(!keep))
}
