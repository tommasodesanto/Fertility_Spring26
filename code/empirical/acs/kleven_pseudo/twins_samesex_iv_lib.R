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
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

.safe_key <- function(d, cols) {
  do.call(paste, c(lapply(cols, function(v) as.character(d[[v]])), sep = "_"))
}

age_missing_codes <- c(999)
rooms_valid_codes <- c(1:27, 30)
rooms_cap <- 9
bedrooms_valid_codes <- 1:22
bedrooms_cap <- 5
acs_year_lo <- 2005L
acs_year_hi <- 2019L

#' Restrict to the common national ACS 1-year product, 2005-2019, before any
#' roster fitting. SAMPLE == YEAR*100 + 1 is the documented IPUMS ACS 1-year
#' code (source_audit_extract27_20260919.md records observed codes 200701,
#' 201101, 201501, 201901, 202301, all matching this pattern; not invented
#' here). Returns the gated data.table plus exclusion counts, so a downstream
#' driver can report what the sample gate dropped rather than silently
#' filtering.
apply_source_sample_gate <- function(dt, year_lo = acs_year_lo, year_hi = acs_year_hi) {
  dt <- data.table::as.data.table(dt)
  n0 <- nrow(dt)
  year_num <- suppressWarnings(as.numeric(dt$YEAR))
  sample_num <- suppressWarnings(as.numeric(dt$SAMPLE))
  in_year <- !is.na(year_num) & year_num >= year_lo & year_num <= year_hi
  in_product <- !is.na(sample_num) & !is.na(year_num) & sample_num == (year_num * 100 + 1)
  keep <- in_year & in_product
  gated <- dt[keep]
  list(
    data = gated,
    counts = list(
      n_input = n0,
      n_excluded_year_out_of_range = sum(!in_year),
      n_excluded_non_acs1yr_product = sum(in_year & !in_product),
      n_kept = nrow(gated),
      year_lo = year_lo, year_hi = year_hi
    )
  )
}

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
  # RELATE is the child's relationship to the HOUSEHOLD HEAD/HOUSEHOLDER, not
  # to the MOMLOC-linked mother (who need not be the householder, e.g. a
  # young mother living with her own parents). It cannot establish or refute
  # biological/adoptive status of the mother link, so it is retained only as
  # a raw descriptive tabulation (relationship-to-householder), never used to
  # filter or reclassify a valid MOMLOC link. Biological relation to the
  # linked mother is unresolved by this source and is not claimed here.

  setorder(valid_links, mother_person_key, -child_age)
  roster <- valid_links[, .(
    linked_ages = list(child_age),
    linked_sexes = list(child_sex),
    linked_relate = list(child_relate),
    linked_child_count = .N,
    any_sex_missing = any(sex_missing)
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
  # Descriptive only: relationship of THIS mother row to her own household's
  # head. RELATE==1 (Head/Householder) is the standard IPUMS-USA harmonized
  # code; other values (e.g. she is a child, sibling, or other relative of
  # the householder) do not affect roster validity, only interpretation.
  mother_rows[, mother_is_householder := RELATE_norm == 1]

  # Relationship-to-householder distribution among linked children, raw
  # codes only, for the receipt -- not a biology claim.
  child_relate_audit <- valid_links[, .N, by = child_relate]

  list(mother_rows = mother_rows, link_audit = link_audit,
       child_relate_audit = child_relate_audit)
}

#' Restrict the mother roster to households where the oldest linked child is
#' a minor (<18), the approved-design sample scope for both Twin1 and
#' SameSex2. Mothers with zero linked children pass through unaffected (they
#' are excluded downstream by the Twin1/SameSex2 eligibility rules on other
#' grounds, not this gate). Returns exclusion counts alongside the gated
#' roster so a driver can report what this filter removed.
apply_oldest_child_minor_gate <- function(mother_rows, minor_age_max = 17) {
  d <- data.table::copy(mother_rows)
  d[, oldest_child_adult := linked_child_count >= 1 & !is.na(a1) & a1 > minor_age_max]
  n_excluded <- sum(d$oldest_child_adult, na.rm = TRUE)
  list(data = d[oldest_child_adult == FALSE],
       n_excluded_oldest_child_adult = n_excluded,
       n_input = nrow(d))
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

  # Valid raw codes are 1:22 (22 is a top code for 21+ bedrooms); codes 0/23+
  # and other out-of-range values are unknown/not-in-universe and must stay
  # missing, matching run_second_birth_housing_readout.R:29's reviewed rule.
  d[, BEDROOMS_num := suppressWarnings(as.numeric(BEDROOMS))]
  d[, bedrooms_valid := !is.na(BEDROOMS_num) & BEDROOMS_num %in% bedrooms_valid_codes]
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

#' Anderson-Rubin confidence set for a single-instrument single-endog IV fit
#' via grid inversion: for each candidate beta0, test H0: instrument
#' coefficient = 0 in the regression of (y - beta0*d) on instrument+controls
#' via fixest::wald(); beta0 is in the AR set if not rejected at `alpha`.
#'
#' A finite grid can only ever report a truncated view of the true
#' (possibly unbounded) AR set. This function returns every maximal
#' contiguous run of accepted grid points as a separate component, and
#' flags any component touching a grid edge as extent-unknown (the true AR
#' set may extend past what the grid covers). It never collapses a
#' boundary-touching accepted run into a "bounded" interval, and it keeps
#' regression failures (`error`) distinct from rejections (`reject`) so a
#' string of numerical errors cannot be mistaken for a rejection region.
ar_confidence_set_legacy <- function(data, outcome, endog, instrument, controls_fml,
                               weight_var, cluster_var, grid, alpha = 0.05) {
  n_grid <- length(grid)
  status <- character(n_grid)  # "accept" | "reject" | "error"
  for (i in seq_len(n_grid)) {
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
    if (is.null(fit)) { status[i] <- "error"; next }
    wt <- tryCatch(fixest::wald(fit, instrument, print = FALSE), error = function(e) NULL)
    if (is.null(wt) || is.na(wt$p)) { status[i] <- "error"; next }
    status[i] <- if (wt$p > alpha) "accept" else "reject"
  }

  accepted <- status == "accept"
  n_accepted <- sum(accepted)
  n_errors <- sum(status == "error")

  components <- list()
  if (n_accepted > 0) {
    idx <- which(accepted)
    breaks <- which(diff(idx) > 1)
    starts <- c(idx[1], idx[breaks + 1])
    ends <- c(idx[breaks], idx[length(idx)])
    for (k in seq_along(starts)) {
      lo_idx <- starts[k]; hi_idx <- ends[k]
      touches_lo <- lo_idx == 1L
      touches_hi <- hi_idx == n_grid
      components[[k]] <- list(
        lower = grid[lo_idx], upper = grid[hi_idx],
        touches_grid_boundary = touches_lo || touches_hi,
        extent = if (touches_lo || touches_hi) "grid_truncated_extent_unknown" else "interior_bounded"
      )
    }
  }

  any_boundary_touch <- length(components) > 0 &&
    any(vapply(components, function(c) c$touches_grid_boundary, logical(1)))
  all_grid_accepted <- n_accepted == n_grid

  # A grid point that errored is neither accepted nor rejected: its true
  # accept/reject status is unknown, so it could hide a wider (or
  # disconnected) true AR set. Any error anywhere on the grid means set
  # completeness/extent is unknown -- never certify "fully bounded" when
  # n_errors > 0, even if every non-error point looks interior.
  list(
    grid_min = min(grid), grid_max = max(grid), n_grid = n_grid,
    n_accepted = n_accepted, n_rejected = sum(status == "reject"),
    n_errors = n_errors,
    n_components = length(components), components = components,
    fully_interior_bounded = length(components) > 0 && !any_boundary_touch && n_errors == 0,
    extent_unknown_due_to_errors = n_errors > 0,
    all_grid_points_accepted = all_grid_accepted,
    empty_accepted_set = n_accepted == 0,
    # Convenience summary across ALL accepted points (may span >1 component
    # or touch the boundary); NOT a claim of a single bounded CI -- consult
    # fully_interior_bounded / components before interpreting as a CI.
    summary_lower = if (n_accepted > 0) min(grid[accepted]) else NA_real_,
    summary_upper = if (n_accepted > 0) max(grid[accepted]) else NA_real_
  )
}

#' Fast Anderson-Rubin grid inversion: computationally equivalent to
#' ar_confidence_set_legacy() for a single instrument/single endogenous
#' regressor with IDENTICAL X/sample/weights/clustering across the Y, D, and
#' Y+D regressions, but runs exactly 3 regressions total instead of
#' length(grid) regressions. For beta0, the residualized instrument
#' coefficient is pi(beta0) = pi_Y - beta0*pi_D with variance
#' V(beta0) = V_Y - 2*beta0*C + beta0^2*V_D, where C is the covariance of
#' the two instrument coefficients recovered from Var(pi_Y+pi_D) via a third
#' regression of (Y+D) on the same RHS: C = (V_sum - V_Y - V_D)/2. This is
#' the same finite-sample Wald test as the legacy per-grid-point regression
#' (same df, same clustered-SE convention), not a different or new
#' scientific estimator -- it is an algebraic shortcut. If V(beta0) <= 0 the
#' grid point is an explicit "error" (invalid variance), never silently
#' accepted or rejected.
ar_confidence_set <- function(data, outcome, endog, instrument, controls_fml,
                               weight_var, cluster_var, grid, alpha = 0.05) {
  w_fml <- stats::as.formula(paste0("~", weight_var))
  c_fml <- stats::as.formula(paste0("~", cluster_var))
  rhs <- if (nzchar(controls_fml)) paste0(instrument, " + ", controls_fml) else instrument
  fit_one <- function(lhs_expr) {
    dloc <- data.table::copy(data)
    dloc[, .ar_lhs := eval(parse(text = lhs_expr), envir = dloc)]
    fit <- tryCatch(fixest::feols(stats::as.formula(paste0(".ar_lhs ~ ", rhs)),
                                   data = dloc, weights = w_fml, cluster = c_fml, notes = FALSE),
                     error = function(e) NULL)
    if (is.null(fit)) return(list(pi = NA_real_, V = NA_real_, df2 = NA_real_))
    wt <- tryCatch(fixest::wald(fit, instrument, print = FALSE), error = function(e) NULL)
    if (is.null(wt)) return(list(pi = NA_real_, V = NA_real_, df2 = NA_real_))
    b <- unname(stats::coef(fit)[instrument])
    v <- unname(diag(stats::vcov(fit))[instrument])
    list(pi = b, V = v, df2 = unname(wt$df2))
  }
  fy <- fit_one(outcome)
  fd <- fit_one(endog)
  fs <- fit_one(paste0("(", outcome, ") + (", endog, ")"))
  base_ok <- !is.na(fy$pi) && !is.na(fd$pi) && !is.na(fs$pi) &&
    !is.na(fy$V) && !is.na(fd$V) && !is.na(fs$V)
  n_grid <- length(grid)
  if (!base_ok) {
    return(list(grid_min = suppressWarnings(min(grid, na.rm = TRUE)),
                grid_max = suppressWarnings(max(grid, na.rm = TRUE)), n_grid = n_grid,
                n_accepted = 0L, n_rejected = 0L, n_errors = n_grid,
                n_components = 0L, components = list(),
                fully_interior_bounded = FALSE, extent_unknown_due_to_errors = TRUE,
                all_grid_points_accepted = FALSE, empty_accepted_set = TRUE,
                summary_lower = NA_real_, summary_upper = NA_real_,
                base_regression_failure = TRUE))
  }
  C <- (fs$V - fy$V - fd$V) / 2
  df2 <- fy$df2
  status <- character(n_grid)
  for (i in seq_len(n_grid)) {
    b0 <- grid[i]
    if (is.na(b0)) { status[i] <- "error"; next }
    pi_b0 <- fy$pi - b0 * fd$pi
    V_b0 <- fy$V - 2 * b0 * C + b0^2 * fd$V
    if (is.na(V_b0) || V_b0 <= 0 || is.na(df2)) { status[i] <- "error"; next }
    wald_stat <- pi_b0^2 / V_b0
    p <- stats::pf(wald_stat, df1 = 1, df2 = df2, lower.tail = FALSE)
    status[i] <- if (is.na(p)) "error" else if (p > alpha) "accept" else "reject"
  }
  accepted <- status == "accept"
  n_accepted <- sum(accepted)
  n_errors <- sum(status == "error")
  components <- list()
  if (n_accepted > 0) {
    idx <- which(accepted)
    breaks <- which(diff(idx) > 1)
    starts <- c(idx[1], idx[breaks + 1]); ends <- c(idx[breaks], idx[length(idx)])
    for (k in seq_along(starts)) {
      lo_idx <- starts[k]; hi_idx <- ends[k]
      touches <- lo_idx == 1L || hi_idx == n_grid
      components[[k]] <- list(lower = grid[lo_idx], upper = grid[hi_idx],
                               touches_grid_boundary = touches,
                               extent = if (touches) "grid_truncated_extent_unknown" else "interior_bounded")
    }
  }
  any_boundary_touch <- length(components) > 0 &&
    any(vapply(components, function(c) c$touches_grid_boundary, logical(1)))
  list(grid_min = min(grid), grid_max = max(grid), n_grid = n_grid,
       n_accepted = n_accepted, n_rejected = sum(status == "reject"), n_errors = n_errors,
       n_components = length(components), components = components,
       fully_interior_bounded = length(components) > 0 && !any_boundary_touch && n_errors == 0,
       extent_unknown_due_to_errors = n_errors > 0,
       all_grid_points_accepted = n_accepted == n_grid,
       empty_accepted_set = n_accepted == 0,
       summary_lower = if (n_accepted > 0) min(grid[accepted]) else NA_real_,
       summary_upper = if (n_accepted > 0) max(grid[accepted]) else NA_real_,
       base_regression_failure = FALSE)
}

#' Fit RF, FS, and 2SLS for one outcome/instrument/treatment triple, with an
#' AR grid-search diagnostic. controls_fml is a fixest-style RHS string
#' (e.g. "poly(mat_age,2) + i(race) + i(survey_year) + i(event_age)").
fit_instrument_outcome <- function(data, outcome, treatment, instrument,
                                    controls_fml, weight_var = "mother_weight",
                                    cluster_var = "household_key",
                                    ar_grid = NULL) {
  # Project to only the columns the formulas actually use BEFORE any copy,
  # so the (possibly very wide, list-column-bearing) national mother roster
  # is never carried into feols() or the AR loop. all.vars() correctly
  # skips fixest's i()/poly() wrapper function names and returns only the
  # underlying variable symbols.
  ctrl_vars <- if (nzchar(controls_fml)) all.vars(stats::as.formula(paste("~", controls_fml))) else character(0)
  needed_vars <- unique(c(outcome, treatment, instrument, weight_var, cluster_var, ctrl_vars))
  missing_vars <- setdiff(needed_vars, names(data))
  if (length(missing_vars) > 0) stop(sprintf("fit_instrument_outcome: missing required columns: %s",
                                              paste(missing_vars, collapse = ", ")))
  data <- data.table::copy(data[, ..needed_vars])
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

  extract_coef <- function(fit, coefname, z = stats::qnorm(0.975)) {
    if (is.null(fit) || is.na(coefname) || !(coefname %in% names(stats::coef(fit))))
      return(list(coef = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                  nobs = NA_integer_, error = "coefficient_not_found"))
    b <- unname(stats::coef(fit)[coefname])
    se <- unname(sqrt(diag(stats::vcov(fit)))[coefname])
    list(coef = b, se = se, ci_lower = b - z * se, ci_upper = b + z * se,
         nobs = stats::nobs(fit), error = NA_character_)
  }
  # Compact full receipt (named b, full clustered V, formula, weights/cluster
  # definitions, actual nobs) -- extracted values only, never the fit object
  # itself, so a small per-fit JSON/CSV can be written without serializing
  # fixest internals.
  full_receipt <- function(fit, label) {
    if (is.null(fit)) return(list(stage = label, status = "fit_failed"))
    list(stage = label, status = "ok",
         formula = deparse(stats::formula(fit)),
         weight_var = weight_var, cluster_var = cluster_var,
         nobs = stats::nobs(fit),
         b = as.list(stats::coef(fit)),
         V = apply(as.matrix(stats::vcov(fit)), 1, as.list))
  }

  rf_err <- if (is.null(rf_fit)) "rf_fit_failed" else NA_character_
  fs_err <- if (is.null(fs_fit)) "fs_fit_failed" else NA_character_
  rfx <- extract_coef(rf_fit, instrument)
  fsx <- extract_coef(fs_fit, instrument)
  out$rf_coef <- rfx$coef; out$rf_se <- rfx$se
  out$rf_ci_lower <- rfx$ci_lower; out$rf_ci_upper <- rfx$ci_upper
  out$rf_nobs_fit <- rfx$nobs
  out$rf_error <- rf_err %||% rfx$error
  out$rf_receipt <- full_receipt(rf_fit, "rf")
  if (!is.null(rf_fit)) {
    fit_idx <- tryCatch(fixest::obs(rf_fit), error = function(e) NULL)
    if (!is.null(fit_idx)) {
      out$n_households_fit_rf <- length(unique(usable[[cluster_var]][fit_idx]))
      out$n_instrument_positive_fit_rf <- sum(usable[[instrument]][fit_idx] == 1, na.rm = TRUE)
    }
  }
  out$fs_coef <- fsx$coef; out$fs_se <- fsx$se
  out$fs_ci_lower <- fsx$ci_lower; out$fs_ci_upper <- fsx$ci_upper
  out$fs_nobs_fit <- fsx$nobs
  out$fs_error <- fs_err %||% fsx$error
  out$fs_receipt <- full_receipt(fs_fit, "fs")

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
  ivx <- extract_coef(iv_fit, iv_coefname)
  out$iv_coef <- ivx$coef; out$iv_se <- ivx$se
  out$iv_ci_lower <- ivx$ci_lower; out$iv_ci_upper <- ivx$ci_upper
  out$iv_nobs_fit <- ivx$nobs
  out$iv_error <- if (is.null(iv_fit)) "iv_fit_failed" else ivx$error
  out$iv_receipt <- full_receipt(iv_fit, "iv")

  if (!is.null(ar_grid) && !is.null(iv_fit) && !is.na(out$iv_coef)) {
    ar <- tryCatch(ar_confidence_set(usable, outcome, treatment, instrument,
                                      controls_fml, weight_var, cluster_var, ar_grid),
                   error = function(e) NULL)
    if (!is.null(ar)) {
      out$ar_summary_lower <- ar$summary_lower; out$ar_summary_upper <- ar$summary_upper
      out$ar_fully_interior_bounded <- ar$fully_interior_bounded
      out$ar_extent_unknown_due_to_errors <- ar$extent_unknown_due_to_errors
      out$ar_n_components <- ar$n_components
      out$ar_all_grid_points_accepted <- ar$all_grid_points_accepted
      out$ar_empty_accepted_set <- ar$empty_accepted_set
      out$ar_n_errors <- ar$n_errors
      out$ar_grid_min <- ar$grid_min; out$ar_grid_max <- ar$grid_max
      out$ar_components_json <- jsonlite::toJSON(ar$components, auto_unbox = TRUE)
    } else {
      out$ar_error <- "ar_confidence_set_failed"
    }
  }
  # RF is the primary estimand and is always reported when it fits, even if
  # FS/IV fail (weak/absent first stage does not invalidate the RF). Status
  # distinguishes a clean full fit from a partial failure so downstream
  # readers cannot mistake a broken IV/FS leg for a complete result.
  # full_fit requires the fit object AND its required coefficient to have
  # actually been extracted (not "coefficient_not_found"/NA) -- a fitted
  # model missing the coefficient of interest is a partial failure, not a
  # clean fit, even though feols() itself did not error.
  rf_ok <- !is.null(rf_fit) && is.na(rfx$error) && !is.na(rfx$coef)
  fs_ok <- !is.null(fs_fit) && is.na(fsx$error) && !is.na(fsx$coef)
  iv_ok <- !is.null(iv_fit) && is.na(ivx$error) && !is.na(ivx$coef)
  out$status <- if (!rf_ok) {
    "error"
  } else if (!fs_ok || !iv_ok) {
    "partial_failure"
  } else {
    "full_fit"
  }
  out
}
