# Twin1 / SameSex2 roster, instrument, and estimation library.
#
# Twin1: among mothers with a valid linked child and oldest-child age 0:5,
# instrument = first two oldest linked children observed at the same age
# (age-only proxy; ACS extract27 has no BIRTHQTR), no older linked child.
# Treatment = >=2 linked children at interview. NOT a confirmed-twin
# indicator and NOT the completed first-birth pseudo-panel event object.
#
# SameSex2: among mothers with >=2 valid linked children, Z = same sex of
# the age-ordered oldest two (ties excluded from primary group), D = >=3
# linked children. This is a second-birth-to-third-child margin, not a
# first-birth or second-birth-fertility hazard object.
#
# Neither instrument is treated as identified: RF is primary, 2SLS/AR are
# reported as assumption-dependent diagnostics per the housing exclusion
# concern (same-sex composition can directly change room-sharing demand;
# twin status associates with maternal health/spacing).

suppressMessages({
  library(data.table)
  library(fixest)
})

.safe_key <- function(d, cols) {
  do.call(paste, c(lapply(cols, function(v) as.character(d[[v]])), sep = "_"))
}

age_missing_codes <- c(999)
rooms_valid_codes <- c(1:27, 30)
rooms_cap <- 9
bedrooms_cap <- 5

#' Build one row per eligible mother with the ordered linked-child age/sex
#' roster. Mirrors the MOMLOC join and validity rules in
#' build_second_birth_proxy.R:135-166, generalized to keep all linked
#' children (not just the first three) and to carry child SEX.
build_mother_roster <- function(dt, female_code = 2, max_linked = 10) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  person_cols <- c("YEAR", "SAMPLE", "SERIAL", "PERNUM")
  hh_cols <- c("YEAR", "SAMPLE", "SERIAL")
  bad_key <- vapply(person_cols, function(v) {
    z <- suppressWarnings(as.numeric(dt[[v]]))
    any(is.na(dt[[v]]) | is.na(z) | !is.finite(z))
  }, logical(1))
  if (any(bad_key)) {
    stop(sprintf("person key fields contain missing/nonfinite values: %s",
                  paste(person_cols[bad_key], collapse = ", ")))
  }
  dt[, person_key := .safe_key(.SD, person_cols), .SDcols = person_cols]
  dt[, household_key := .safe_key(.SD, hh_cols), .SDcols = hh_cols]
  if (anyDuplicated(dt$person_key) > 0L) {
    stop("person key (YEAR,SAMPLE,SERIAL,PERNUM) is not unique")
  }

  dt[, AGE_norm := as.numeric(AGE)]
  dt[AGE_norm %in% age_missing_codes, AGE_norm := NA_real_]
  dt[, NCHILD_norm := as.numeric(NCHILD)]
  dt[, MOMLOC_norm := as.numeric(MOMLOC)]
  dt[, weight_invalid := is.na(PERWT) | !is.finite(as.numeric(PERWT)) |
       as.numeric(PERWT) <= 0]
  dt[, RELATE_norm := as.numeric(RELATE)]

  mother <- dt[!is.na(SEX) & SEX == female_code & !is.na(AGE_norm)]
  if (nrow(mother) == 0L) stop("no female mother rows with valid AGE")

  child <- dt[!is.na(MOMLOC_norm) & MOMLOC_norm > 0,
              .(YEAR, SAMPLE, SERIAL, child_pernum = PERNUM,
                mother_pernum = MOMLOC_norm, child_age = AGE_norm,
                child_sex = SEX, child_relate = RELATE_norm)]

  links <- merge(
    child,
    mother[, .(YEAR, SAMPLE, SERIAL, mother_pernum = as.numeric(PERNUM),
               mother_person_key = person_key, mother_age = AGE_norm)],
    by = c("YEAR", "SAMPLE", "SERIAL", "mother_pernum"),
    all = TRUE, allow.cartesian = FALSE
  )
  links[, is_mother_placeholder := is.na(child_pernum) & !is.na(mother_person_key)]
  links[, is_unmatched_child := !is.na(child_pernum) & is.na(mother_person_key)]
  links[, momloc_self := !is.na(child_pernum) & child_pernum == mother_pernum]
  links[, link_valid := !is.na(child_age) & !is.na(mother_age) &
          !momloc_self & child_age >= 0 & child_age < mother_age]
  links[, link_invalid_reason := data.table::fcase(
    is_mother_placeholder, "mother_without_link",
    is_unmatched_child, "child_without_female_mother",
    is.na(child_age), "missing child age",
    is.na(mother_age), "missing mother age",
    momloc_self, "MOMLOC self-link",
    !is.na(child_age) & child_age < 0, "negative child age",
    !is.na(child_age) & !is.na(mother_age) & child_age >= mother_age,
      "child age not younger than mother",
    default = NA_character_
  )]
  link_audit <- links[, .N, by = .(link_invalid_reason, child_relate)]

  valid_links <- links[link_valid == TRUE]
  valid_links[, sex_missing := is.na(child_sex) | !(child_sex %in% c(1, 2))]
  valid_links[, non_biological_relate := !is.na(child_relate) &
                !(child_relate %in% c(3, 4))]  # 3 biological, 4 adopted (IPUMS RELATE)

  setorder(valid_links, mother_person_key, -child_age)
  roster <- valid_links[, .(
    linked_ages = list(child_age),
    linked_sexes = list(child_sex),
    linked_relate = list(child_relate),
    linked_child_count = .N,
    any_sex_missing = any(sex_missing),
    any_non_biological = any(non_biological_relate)
  ), by = mother_person_key]

  mother_rows <- merge(mother, roster, by.x = "person_key",
                        by.y = "mother_person_key", all.x = TRUE, sort = FALSE)
  mother_rows[is.na(linked_child_count), linked_child_count := 0L]
  mother_rows[is.na(linked_ages), linked_ages := list(numeric(0)), by = person_key]
  mother_rows[is.na(linked_sexes), linked_sexes := list(numeric(0)), by = person_key]

  for (k in seq_len(max_linked)) {
    mother_rows[, paste0("a", k) := vapply(linked_ages, function(z)
      if (length(z) >= k) z[k] else NA_real_, numeric(1))]
    mother_rows[, paste0("sex", k) := vapply(linked_sexes, function(z)
      if (length(z) >= k) z[k] else NA_real_, numeric(1))]
  }
  mother_rows[, nchild_link_mismatch := is.na(NCHILD_norm) |
                NCHILD_norm != linked_child_count]
  mother_rows[, any_sex_missing := replace(any_sex_missing, is.na(any_sex_missing), FALSE)]
  mother_rows[, any_non_biological := replace(any_non_biological, is.na(any_non_biological), FALSE)]

  list(mother_rows = mother_rows, link_audit = link_audit)
}

#' Recode housing outcomes under the reviewed coding contract.
add_outcomes <- function(d) {
  d <- data.table::copy(d)
  d[, ROOMS_num := suppressWarnings(as.numeric(ROOMS))]
  d[, rooms_valid := !is.na(ROOMS_num) & ROOMS_num %in% rooms_valid_codes]
  d[, ROOMS_out := ifelse(rooms_valid, pmin(ROOMS_num, rooms_cap), NA_real_)]

  d[, OWNERSHP_num := suppressWarnings(as.numeric(OWNERSHP))]
  d[, ownershp_valid := OWNERSHP_num %in% c(1, 2)]
  d[, OWNERSHP_out := ifelse(ownershp_valid, as.numeric(OWNERSHP_num == 1), NA_real_)]

  d[, BEDROOMS_num := suppressWarnings(as.numeric(BEDROOMS))]
  d[, bedrooms_valid := !is.na(BEDROOMS_num) & BEDROOMS_num >= 1]
  d[, BEDROOMS_out := ifelse(bedrooms_valid, pmin(BEDROOMS_num - 1, bedrooms_cap), NA_real_)]
  d
}

#' Twin1: constructed separately AT EACH observed oldest-child age 0:5
#' (repeated cross-section, not a followed cohort). age_grid defaults to
#' 0:5 per the approved design.
build_twin1 <- function(mother_rows, age_grid = 0:5) {
  d <- data.table::copy(mother_rows)
  d[, oldest_age := a1]
  d[, eligible := linked_child_count >= 1 & !is.na(oldest_age) &
       oldest_age %in% age_grid]
  d[, twin_like_proxy := eligible & linked_child_count >= 2 &
       !is.na(a2) & a1 == a2]
  d[, contamination_risk := eligible & linked_child_count >= 3 &
       !is.na(a3) & a1 == a3]  # >2 same-age children at the oldest age
  d[, treatment_2plus := linked_child_count >= 2]
  d[, event_age := oldest_age]
  d[, event_year := YEAR - event_age]
  d[eligible == TRUE]
}

#' SameSex2: among mothers with >=2 valid linked children, age-ordered
#' oldest two, ties excluded from primary group. Event age = second-oldest
#' child's age (a2), grid 0:5.
build_samesex2 <- function(mother_rows, age_grid = 0:5) {
  d <- data.table::copy(mother_rows)
  d[, eligible_pool := linked_child_count >= 2]
  d[, primary_age_tie := eligible_pool & !is.na(a1) & !is.na(a2) & a1 == a2]
  d[, event_age := a2]
  d[, in_event_window := eligible_pool & !is.na(event_age) & event_age %in% age_grid]
  d[, eligible := in_event_window & !primary_age_tie &
       !is.na(sex1) & !is.na(sex2) & sex1 %in% c(1, 2) & sex2 %in% c(1, 2)]
  d[, samesex := as.numeric(sex1 == sex2)]
  d[, both_boys := eligible & sex1 == 1 & sex2 == 1]
  d[, both_girls := eligible & sex1 == 2 & sex2 == 2]
  d[, treatment_3plus := linked_child_count >= 3]
  d[, event_year := YEAR - event_age]
  d[eligible_pool == TRUE]
}

#' Anderson-Rubin confidence set for a single-instrument single-endog IV
#' fit via grid inversion: for each candidate beta0, test H0: instrument
#' coefficient = 0 in the regression of (y - beta0*d) on instrument+controls
#' via fixest::wald(); beta0 is in the AR set if not rejected at `alpha`.
ar_confidence_set <- function(data, outcome, endog, instrument, controls_fml,
                               weight_var, cluster_var, grid, alpha = 0.05) {
  keep <- rep(NA, length(grid))
  for (i in seq_along(grid)) {
    b0 <- grid[i]
    yy <- data[[outcome]] - b0 * data[[endog]]
    dloc <- data.table::copy(data)
    dloc[, .ar_y := yy]
    fml <- stats::as.formula(paste0(".ar_y ~ ", instrument,
                                     if (nzchar(controls_fml)) paste0(" + ", controls_fml) else ""))
    fit <- tryCatch(
      fixest::feols(fml, data = dloc, weights = stats::as.formula(paste0("~", weight_var)),
                    cluster = stats::as.formula(paste0("~", cluster_var)), notes = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) { keep[i] <- NA; next }
    wt <- tryCatch(fixest::wald(fit, instrument, print = FALSE), error = function(e) NULL)
    keep[i] <- if (is.null(wt)) NA else (wt$p > alpha)
  }
  in_set <- grid[which(keep)]
  list(grid = grid, keep = keep,
       bounded = length(in_set) > 0 && length(in_set) < length(grid),
       lower = if (length(in_set) > 0) min(in_set) else NA_real_,
       upper = if (length(in_set) > 0) max(in_set) else NA_real_,
       disconnected = length(in_set) > 0 &&
         any(diff(sort(unique(c(which(keep), NA))), na.rm = TRUE) > 1, na.rm = TRUE))
}

#' Fit RF, FS, and 2SLS for one outcome/instrument/treatment triple, with an
#' AR grid-search diagnostic. controls_fml is a fixest-style RHS string
#' (e.g. "poly(mat_age,2) + i(race) + i(survey_year) + i(event_age)").
fit_instrument_outcome <- function(data, outcome, treatment, instrument,
                                    controls_fml, weight_var = "mother_weight",
                                    cluster_var = "household_key",
                                    ar_grid = NULL) {
  data <- data.table::copy(data)
  data[[instrument]] <- as.numeric(data[[instrument]])
  data[[treatment]] <- as.numeric(data[[treatment]])
  usable <- data[!is.na(data[[outcome]]) & !is.na(data[[treatment]]) &
                   !is.na(data[[instrument]])]
  n_usable <- nrow(usable)
  n_hh <- length(unique(usable[[cluster_var]]))
  n_instr_pos <- sum(usable[[instrument]] == 1, na.rm = TRUE)
  out <- list(outcome = outcome, treatment = treatment, instrument = instrument,
              n_usable = n_usable, n_households = n_hh,
              n_instrument_positive = n_instr_pos)
  if (n_usable < 30 || n_instr_pos < 5 || (n_instr_pos == n_usable)) {
    out$status <- "insufficient_support"
    return(out)
  }
  w_fml <- stats::as.formula(paste0("~", weight_var))
  c_fml <- stats::as.formula(paste0("~", cluster_var))
  rhs <- if (nzchar(controls_fml)) paste0(instrument, " + ", controls_fml) else instrument

  rf_fit <- tryCatch(fixest::feols(
    stats::as.formula(paste0(outcome, " ~ ", rhs)),
    data = usable, weights = w_fml, cluster = c_fml, notes = FALSE),
    error = function(e) NULL)
  fs_fit <- tryCatch(fixest::feols(
    stats::as.formula(paste0(treatment, " ~ ", rhs)),
    data = usable, weights = w_fml, cluster = c_fml, notes = FALSE),
    error = function(e) NULL)

  out$rf_coef <- if (!is.null(rf_fit)) unname(stats::coef(rf_fit)[instrument]) else NA_real_
  out$rf_se <- if (!is.null(rf_fit)) unname(sqrt(diag(stats::vcov(rf_fit)))[instrument]) else NA_real_
  out$fs_coef <- if (!is.null(fs_fit)) unname(stats::coef(fs_fit)[instrument]) else NA_real_
  out$fs_se <- if (!is.null(fs_fit)) unname(sqrt(diag(stats::vcov(fs_fit)))[instrument]) else NA_real_
  fs_wald <- if (!is.null(fs_fit)) tryCatch(fixest::wald(fs_fit, instrument, print = FALSE), error = function(e) NULL) else NULL
  out$first_stage_F <- if (!is.null(fs_wald)) unname(fs_wald$stat) else NA_real_
  out$first_stage_F_df <- if (!is.null(fs_wald))
    sprintf("%s,%s", fs_wald$df1, fs_wald$df2) else NA_character_

  iv_rhs <- paste0(controls_fml, " | ", treatment, " ~ ", instrument)
  iv_fit <- tryCatch(fixest::feols(
    stats::as.formula(paste0(outcome, " ~ ", iv_rhs)),
    data = usable, weights = w_fml, cluster = c_fml, notes = FALSE),
    error = function(e) NULL)
  iv_coefname <- paste0("fit_", treatment)
  out$iv_coef <- if (!is.null(iv_fit)) unname(stats::coef(iv_fit)[iv_coefname]) else NA_real_
  out$iv_se <- if (!is.null(iv_fit)) unname(sqrt(diag(stats::vcov(iv_fit)))[iv_coefname]) else NA_real_

  if (!is.null(ar_grid) && !is.null(iv_fit) && !is.na(out$iv_coef)) {
    ar <- tryCatch(ar_confidence_set(usable, outcome, treatment, instrument,
                                      controls_fml, weight_var, cluster_var, ar_grid),
                   error = function(e) NULL)
    if (!is.null(ar)) {
      out$ar_lower <- ar$lower; out$ar_upper <- ar$upper
      out$ar_bounded <- ar$bounded; out$ar_disconnected <- ar$disconnected
      out$ar_grid_min <- min(ar_grid); out$ar_grid_max <- max(ar_grid)
    }
  }
  out$status <- "fit"
  out
}
