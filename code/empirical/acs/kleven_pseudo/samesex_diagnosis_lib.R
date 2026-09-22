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

## ---- A) Reproduction gate: identity BEFORE read, then numeric refit ------
#' Verifies the canonical analytic-frames path/size/MD5 match the expected
#' identity BEFORE reading the RDS, and that the saved receipt is for the
#' expected outcome/instrument/treatment, THEN refits RF/FS and compares
#' coefficients/SEs/nobs to the saved receipt within `tolerance` (numeric
#' comparison, not identical(), since JSON round-trips integers as
#' doubles). Fails loudly at the first mismatched gate, before any
#' scientific diagnostic runs.
verify_reproduction <- function(analytic_frames_path, expected_size, expected_md5,
                                 saved_receipt_path, expected_outcome = "ROOMS_out",
                                 expected_instrument = "samesex", expected_treatment = "treatment_3plus",
                                 controls_fml, weight_var, cluster_var, tolerance = 1e-8) {
  gate <- list(path_exists = file.exists(analytic_frames_path))
  if (!gate$path_exists) {
    return(list(pass = FALSE, gate = gate, stage = "path_missing"))
  }
  actual_size <- file.info(analytic_frames_path)$size
  gate$size_match <- isTRUE(actual_size == expected_size)
  gate$actual_size <- actual_size
  if (!gate$size_match) {
    return(list(pass = FALSE, gate = gate, stage = "size_mismatch"))
  }
  actual_md5 <- unname(tools::md5sum(analytic_frames_path))
  gate$md5_match <- isTRUE(actual_md5 == expected_md5)
  gate$actual_md5 <- actual_md5
  if (!gate$md5_match) {
    return(list(pass = FALSE, gate = gate, stage = "md5_mismatch"))
  }

  saved <- jsonlite::fromJSON(saved_receipt_path)
  identity_ok <- identical(saved$outcome, expected_outcome) &&
    identical(saved$instrument, expected_instrument) &&
    identical(saved$treatment, expected_treatment)
  gate$saved_receipt_identity_match <- identity_ok
  if (!identity_ok) {
    return(list(pass = FALSE, gate = gate, stage = "saved_receipt_identity_mismatch",
                saved_outcome = saved$outcome, saved_instrument = saved$instrument,
                saved_treatment = saved$treatment))
  }

  af <- readRDS(analytic_frames_path)
  fit <- rf_fs_diagnostic_fit(af$ss, outcome = expected_outcome, treatment = expected_treatment,
                               instrument = expected_instrument, controls_fml = controls_fml,
                               weight_var = weight_var, cluster_var = cluster_var)
  cmp <- function(a, b) isTRUE(abs(as.numeric(a) - as.numeric(b)) <= tolerance)
  checks <- list(
    rf_coef = cmp(fit$rf_coef, saved$rf_coef), rf_se = cmp(fit$rf_se, saved$rf_se),
    fs_coef = cmp(fit$fs_coef, saved$fs_coef), fs_se = cmp(fit$fs_se, saved$fs_se),
    rf_nobs = cmp(fit$rf_nobs, saved$rf_nobs_fit), fs_nobs = cmp(fit$fs_nobs, saved$fs_nobs_fit)
  )
  list(pass = all(unlist(checks)), gate = gate, stage = "numeric_comparison",
       checks = checks, fit = fit, saved = saved)
}

## ---- 3) RF/FS-only diagnostic helper: NEVER fits IV/AR --------------------
#' Small diagnostics-only fitter: full complete-case Y,D,Z,controls,
#' weight,cluster sample; RF (Y~Z+controls) and FS (D~Z+controls) only.
#' Returns named b, full clustered V, nobs, HH count, Z-positive count,
#' coefficient/SE, normal-approximation 95% CI, and the single-instrument
#' cluster-robust Wald F. Explicit errors/warnings are captured and
#' returned, never silently swallowed. This is the ONLY fitting entry
#' point used by the B/D/E diagnostics below -- no 2SLS/AR anywhere here.
rf_fs_diagnostic_fit <- function(data, outcome, treatment, instrument, controls_fml,
                                  weight_var, cluster_var, z = stats::qnorm(0.975)) {
  ctrl_vars <- if (nzchar(controls_fml)) all.vars(stats::as.formula(paste("~", controls_fml))) else character(0)
  needed <- unique(c(outcome, treatment, instrument, weight_var, cluster_var, ctrl_vars))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    return(list(status = "error", error = sprintf("missing columns: %s", paste(missing_cols, collapse = ", "))))
  }
  d <- data.table::copy(data[, ..needed])
  d[[instrument]] <- as.numeric(d[[instrument]])
  d[[treatment]] <- as.numeric(d[[treatment]])
  usable <- d[!is.na(d[[outcome]]) & !is.na(d[[treatment]]) & !is.na(d[[instrument]])]
  n_prefit <- nrow(usable)
  n_hh_prefit <- length(unique(usable[[cluster_var]]))
  n_zpos_prefit <- sum(usable[[instrument]] == 1, na.rm = TRUE)
  if (n_prefit < 10) {
    return(list(status = "insufficient_support", n_usable_prefit = n_prefit,
                n_households_prefit = n_hh_prefit, n_instrument_positive_prefit = n_zpos_prefit))
  }
  rhs <- if (nzchar(controls_fml)) paste0(instrument, " + ", controls_fml) else instrument
  w_fml <- stats::as.formula(paste0("~", weight_var))
  c_fml <- stats::as.formula(paste0("~", cluster_var))

  warns <- character(0)
  wh <- withCallingHandlers(
    {
      rf_fit <- tryCatch(fixest::feols(stats::as.formula(paste0(outcome, " ~ ", rhs)),
                                        data = usable, weights = w_fml, cluster = c_fml, notes = FALSE),
                          error = function(e) e)
      fs_fit <- tryCatch(fixest::feols(stats::as.formula(paste0(treatment, " ~ ", rhs)),
                                        data = usable, weights = w_fml, cluster = c_fml, notes = FALSE),
                          error = function(e) e)
      list(rf_fit = rf_fit, fs_fit = fs_fit)
    },
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
                  nobs = NA_integer_, b = NULL, V = NULL, F = NA_real_))
    b <- stats::coef(fit); V <- as.matrix(stats::vcov(fit))
    coefv <- unname(b[coefname]); se <- unname(sqrt(diag(V))[coefname])
    wt <- tryCatch(fixest::wald(fit, coefname, print = FALSE), error = function(e) NULL)
    list(coef = coefv, se = se, ci_lower = coefv - z * se, ci_upper = coefv + z * se,
         nobs = stats::nobs(fit), b = as.list(b), V = apply(V, 1, as.list),
         F = if (!is.null(wt)) unname(wt$stat) else NA_real_)
  }
  rfx <- extract(rf_fit, instrument); fsx <- extract(fs_fit, instrument)
  list(
    status = if (is.na(rf_err) && is.na(fs_err) && !is.na(rfx$coef) && !is.na(fsx$coef)) "full_fit" else "error",
    n_usable_prefit = n_prefit, n_households_prefit = n_hh_prefit, n_instrument_positive_prefit = n_zpos_prefit,
    rf_coef = rfx$coef, rf_se = rfx$se, rf_ci_lower = rfx$ci_lower, rf_ci_upper = rfx$ci_upper,
    rf_nobs = rfx$nobs, rf_b = rfx$b, rf_V = rfx$V, rf_error = rf_err,
    fs_coef = fsx$coef, fs_se = fsx$se, fs_ci_lower = fsx$ci_lower, fs_ci_upper = fsx$ci_upper,
    fs_nobs = fsx$nobs, fs_b = fsx$b, fs_V = fsx$V, fs_error = fs_err, first_stage_F = fsx$F,
    warnings = warns
  )
}

## ---- B) Event-age-by-age RF/FS (all ages kept, never dropped) ------------
#' event0 is post-second-birth, NOT a pre-treatment/placebo test: a third
#' child can already exist by event age 0 and interview timing/selection
#' matters. Reports every age 0:5 for every outcome -- never selects ages
#' by significance or drops an unfavorable age. RF/FS only (no IV/AR).
event_age_rf_fs_diagnostic <- function(ss, outcomes, controls_fml, weight_var, cluster_var,
                                        ages = 0:5) {
  out <- list()
  for (oc in outcomes) {
    for (ea in ages) {
      d <- ss[event_age == ea]
      r <- rf_fs_diagnostic_fit(d, oc, "treatment_3plus", "samesex", controls_fml, weight_var, cluster_var)
      r$outcome <- oc; r$event_age <- ea
      out[[length(out) + 1]] <- r
    }
  }
  data.table::rbindlist(lapply(out, function(r) data.table(
    outcome = r$outcome, event_age = r$event_age, status = r$status %||% NA_character_,
    n_usable_prefit = r$n_usable_prefit %||% NA_integer_, n_households_prefit = r$n_households_prefit %||% NA_integer_,
    n_instrument_positive_prefit = r$n_instrument_positive_prefit %||% NA_integer_,
    rf_coef = r$rf_coef %||% NA_real_, rf_se = r$rf_se %||% NA_real_,
    rf_ci_lower = r$rf_ci_lower %||% NA_real_, rf_ci_upper = r$rf_ci_upper %||% NA_real_,
    fs_coef = r$fs_coef %||% NA_real_, fs_se = r$fs_se %||% NA_real_,
    first_stage_F = r$first_stage_F %||% NA_real_,
    n_warnings = length(r$warnings %||% character(0))
  )), fill = TRUE)
}

## ---- C) Per-state roster metadata extraction (memory-bounded) ------------
#' Mirrors the production driver's exact per-state pipeline (source/sample
#' gate -> roster -> outcomes) but retains sex1/sex2/a1/a2/a3/
#' linked_child_count/NCHILD_norm/person_key/household_key instead of
#' projecting to the slim regression frame, restricted to the SameSex2
#' eligible_pool. Processes ONE state at a time; the caller rbinds the
#' small metadata (not the wide per-state table) across states and
#' releases each wide state object immediately.
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

#' Exact mother-level (PERSON_KEY = YEAR/SAMPLE/SERIAL/PERNUM) join between
#' freshly extracted metadata and the cached analytic-frame ss. Multiple
#' mothers sharing one household_key is LEGITIMATE (not an error) --
#' hence the join is on person_key, not household_key. Requires: (1)
#' person_key is non-missing and UNIQUE on both sides -- a duplicate
#' person_key is always an error, even with identical values on the
#' duplicate rows; (2) exact key-set equality (no metadata-only or
#' cache-only person_keys); (3) after the join, household_key/samesex/
#' treatment_3plus/event_age match exactly. Fails rather than silently
#' accepting a changed sample.
verify_key_coverage_and_recompute <- function(metadata, cached_ss) {
  meta_e <- metadata[eligible == TRUE]
  dup_meta <- sum(duplicated(meta_e$person_key))
  dup_cache <- sum(duplicated(cached_ss$person_key))
  na_meta <- sum(is.na(meta_e$person_key))
  na_cache <- sum(is.na(cached_ss$person_key))
  keys_meta <- unique(meta_e$person_key)
  keys_cache <- unique(cached_ss$person_key)
  missing_from_meta <- setdiff(keys_cache, keys_meta)  # in cache, not recomputed -> FAIL
  extra_in_meta <- setdiff(keys_meta, keys_cache)      # recomputed, not in cache -> report, not necessarily fail (eligible_pool is broader than primary cache)

  keys_ok <- (dup_meta == 0 && dup_cache == 0 && na_meta == 0 && na_cache == 0 &&
                length(missing_from_meta) == 0)

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
    n_matched_join_rows = nrow(m), n_recompute_mismatches = nrow(mismatched),
    keys_ok = keys_ok,
    one_to_one_and_recomputed_equal = (keys_ok && nrow(m) == nrow(cached_ss) && nrow(mismatched) == 0)
  )
}

## ---- Sex/order/link audits (report, never silently filter) ---------------
#' Valid sex values, descending age order, and third-child age-tie
#' ambiguity (a2==a3, distinct from the already-excluded a1==a2 primary
#' tie). Also NCHILD-vs-linked-count and duplicate-PERSON-key checks
#' (duplicate household_key is expected/legitimate, not flagged as an
#' error -- multiple mothers can share one household).
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
    n_duplicate_person_key = sum(duplicated(metadata$person_key)),
    n_unique_household_key = length(unique(metadata$household_key)),
    n_households_with_multiple_mothers = {
      tab <- table(metadata$household_key); sum(tab > 1)
    }
  )
}

## ---- D) Ordered-sex cells, additive controls, and JOINT BB/GG dummies ----
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
    n_unweighted = .N, n_weighted = sum(get(weight_var)),
    D_mean_unweighted = mean(treatment_3plus), D_mean_weighted = weighted.mean(treatment_3plus, get(weight_var))
  ), by = cell]
  for (oc in outcomes) {
    agg[[paste0(oc, "_Y_mean_unweighted")]] <- vapply(agg$cell, function(cc)
      mean(d[cell == cc][[oc]], na.rm = TRUE), numeric(1))
    agg[[paste0(oc, "_Y_mean_weighted")]] <- vapply(agg$cell, function(cc) {
      sub <- d[cell == cc]; ok <- !is.na(sub[[oc]])
      if (!any(ok)) return(NA_real_)
      stats::weighted.mean(sub[[oc]][ok], sub[[weight_var]][ok])
    }, numeric(1))
  }
  agg
}

#' Adds firstchildboy/secondchildboy indicators (Angrist-Evans 1998-style
#' additive sex controls, not an exact replication -- see file header).
add_first_second_child_sex_controls <- function(ss) {
  d <- data.table::copy(ss)
  d[, firstchildboy := as.numeric(sex1 == 1)]
  d[, secondchildboy := as.numeric(sex2 == 1)]
  d
}

#' RF for the additive-sex-control specification: same-sex instrument PLUS
#' firstchildboy/secondchildboy as extra controls (labeled separate
#' diagnostic; does not replace the primary same-sex RF). RF/FS only.
additive_sex_control_rf_fs <- function(ss, outcome, controls_fml, weight_var, cluster_var) {
  d <- add_first_second_child_sex_controls(ss)
  extended_controls <- paste0(controls_fml, " + firstchildboy + secondchildboy")
  # If either added control is constant/collinear within this outcome's
  # complete-case subset, log it explicitly rather than silently dropping
  # or crashing.
  ctrl_vars <- c("firstchildboy", "secondchildboy")
  const_flags <- vapply(ctrl_vars, function(v) length(unique(stats::na.omit(d[[v]]))) <= 1, logical(1))
  r <- rf_fs_diagnostic_fit(d, outcome, "treatment_3plus", "samesex", extended_controls, weight_var, cluster_var)
  r$constant_added_controls <- names(const_flags)[const_flags]
  r
}

#' JOINT both-boys / both-girls RF and FS in ONE regression per outcome,
#' with mixed (BG+GB) as the implicit omitted reference -- NOT a
#' bothboys-vs-everyone-else (which would wrongly pool in bothgirls) or a
#' bothboys-vs-(bothgirls+mixed) comparison. Reports both dummy
#' coefficients, their SEs, and their covariance (for a downstream
#' BB-vs-GG contrast test using the joint V, not a naive SE subtraction).
joint_bb_gg_rf_fs <- function(ss, outcome, controls_fml, weight_var, cluster_var) {
  d <- data.table::copy(ss)
  d[, bothboys := as.numeric(sex1 == 1 & sex2 == 1)]
  d[, bothgirls := as.numeric(sex1 == 2 & sex2 == 2)]
  ctrl_vars <- if (nzchar(controls_fml)) all.vars(stats::as.formula(paste("~", controls_fml))) else character(0)
  needed <- unique(c(outcome, "treatment_3plus", "bothboys", "bothgirls", weight_var, cluster_var, ctrl_vars))
  usable <- d[, ..needed]
  usable <- usable[!is.na(usable[[outcome]]) & !is.na(bothboys) & !is.na(bothgirls)]
  n_prefit <- nrow(usable); n_hh <- length(unique(usable[[cluster_var]]))
  rhs <- paste0("bothboys + bothgirls", if (nzchar(controls_fml)) paste0(" + ", controls_fml) else "")
  w_fml <- stats::as.formula(paste0("~", weight_var)); c_fml <- stats::as.formula(paste0("~", cluster_var))
  fit_one <- function(lhs) {
    fit <- tryCatch(fixest::feols(stats::as.formula(paste0(lhs, " ~ ", rhs)),
                                   data = usable, weights = w_fml, cluster = c_fml, notes = FALSE),
                     error = function(e) NULL)
    if (is.null(fit)) return(list(status = "error"))
    b <- stats::coef(fit); V <- as.matrix(stats::vcov(fit))
    list(status = "full_fit", nobs = stats::nobs(fit),
         bb_coef = unname(b["bothboys"]), bb_se = unname(sqrt(diag(V))["bothboys"]),
         gg_coef = unname(b["bothgirls"]), gg_se = unname(sqrt(diag(V))["bothgirls"]),
         bb_gg_cov = unname(V["bothboys", "bothgirls"]), b = as.list(b), V = apply(V, 1, as.list))
  }
  list(outcome = outcome, n_usable_prefit = n_prefit, n_households_prefit = n_hh,
       n_bothboys = sum(usable$bothboys), n_bothgirls = sum(usable$bothgirls),
       n_mixed_reference = n_prefit - sum(usable$bothboys) - sum(usable$bothgirls),
       rf = fit_one(outcome), fs = fit_one("treatment_3plus"))
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
