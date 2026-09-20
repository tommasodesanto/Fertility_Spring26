# Build a strict, auditable coresident second-birth proxy from normalized ACS rows.
#
# This file deliberately stops before matching or estimation.  It consumes a
# small data.table with canonical uppercase names, preserves raw fields, and
# returns row identities, links, eligibility flags, donors, and donor targets.
# The local extract labels FERTYR 0=N/A, 1=No, 2=Yes, 8=Suppressed; callers must
# pass those codes explicitly through fertyr_codes.

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("build_second_birth_proxy.R requires the data.table package")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.required_second_birth_fields <- c(
  "YEAR", "SAMPLE", "SERIAL", "PERNUM", "MOMLOC", "AGE", "SEX",
  "NCHILD", "FERTYR", "PERWT"
)

.validate_fertyr_codes <- function(fertyr_codes) {
  if (!is.list(fertyr_codes) ||
      !all(c("yes", "no", "unknown") %in% names(fertyr_codes))) {
    stop("fertyr_codes must explicitly contain yes, no, and unknown")
  }
  if (length(fertyr_codes$yes) != 1L || length(fertyr_codes$no) != 1L) {
    stop("fertyr_codes yes and no must each contain one code")
  }
  all_codes <- c(fertyr_codes$yes, fertyr_codes$no, fertyr_codes$unknown)
  if (anyDuplicated(all_codes) > 0L) {
    stop("fertyr_codes categories must be disjoint")
  }
  invisible(TRUE)
}

.normalize_fertyr <- function(x, fertyr_codes) {
  .validate_fertyr_codes(fertyr_codes)
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x == fertyr_codes$yes] <- "yes"
  out[!is.na(x) & x == fertyr_codes$no] <- "no"
  out[!is.na(x) & x %in% fertyr_codes$unknown] <- "unknown"
  # Values outside the explicit contract remain NA and are counted as invalid.
  out
}

.check_required <- function(dt, required) {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(sprintf("missing required metadata fields: %s",
                 paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

.safe_key <- function(dt, cols) {
  do.call(paste, c(unname(dt[, ..cols]), sep = "/"))
}

.empty_table <- function() data.table::data.table()

#' Build strict anchors, post rows, one-child donors, and event-0 donor targets.
#'
#' @param dt normalized data.table. Raw columns are retained and never replaced.
#' @param fertyr_codes explicit list(yes=2L, no=1L, unknown=c(0L, 8L)).
#' @param match_covariates explicit author covariates to carry to donor targets.
#'   AGE is generated as target_mother_age; other columns are copied from anchors.
#' @param female_code source SEX code for mothers (extract27 uses 2).
#' @param age_missing_codes raw age codes to flag as missing (extract27 uses 999).
#' @param event_window event times used for target expansion; only k < 0 expand.
#' @param age_at_event_bounds strict proxy-birth age band.
#' @param full_pre_min_gap gap required for all -5,...,-1 targets.
#' @param reference_min_gap gap required for author reference t=-2.
build_second_birth_proxy <- function(
    dt,
    fertyr_codes,
    match_covariates,
    female_code = 2L,
    age_missing_codes = 999L,
    event_window = -5:10,
    age_at_event_bounds = c(25, 45),
    full_pre_min_gap = 5L,
    reference_min_gap = 2L) {
  if (!data.table::is.data.table(dt)) {
    stop("dt must be a normalized data.table")
  }
  .check_required(dt, .required_second_birth_fields)
  .validate_fertyr_codes(fertyr_codes)
  if (missing(match_covariates) || length(match_covariates) == 0L) {
    stop("match_covariates must be supplied explicitly")
  }
  missing_match <- setdiff(match_covariates, names(dt))
  if (length(missing_match) > 0L) {
    stop(sprintf("missing required match covariates: %s",
                 paste(missing_match, collapse = ", ")))
  }
  if (length(age_at_event_bounds) != 2L || age_at_event_bounds[1] > age_at_event_bounds[2]) {
    stop("age_at_event_bounds must be an increasing two-element vector")
  }
  if (any(!is.finite(event_window)) || any(event_window != as.integer(event_window))) {
    stop("event_window must contain integer event times")
  }
  event_window <- as.integer(event_window)
  if (!any(event_window < 0L)) stop("event_window must include negative targets")
  if (!any(event_window == 0L)) stop("event_window must include event time zero")

  dt <- data.table::copy(dt)
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
  dt[, FERTYR_status := .normalize_fertyr(FERTYR, fertyr_codes)]
  dt[, weight_invalid := is.na(PERWT) | !is.finite(as.numeric(PERWT)) |
       as.numeric(PERWT) <= 0]

  # Only nonmissing female rows are mother candidates. Missing/other SEX rows
  # remain in the returned input and are counted in the audit.
  mother <- dt[!is.na(SEX) & SEX == female_code & !is.na(AGE_norm)]
  if (nrow(mother) == 0L) stop("no female mother rows with valid AGE")
  mother[, mother_sex_valid := TRUE]

  # MOMLOC points from child to mother. Join on the full household key before
  # comparing MOMLOC with the mother's PERNUM.
  child <- dt[!is.na(MOMLOC_norm) & MOMLOC_norm > 0,
              .(YEAR, SAMPLE, SERIAL, child_pernum = PERNUM,
                mother_pernum = MOMLOC_norm, child_age = AGE_norm,
                child_age_raw = AGE)]
  links <- merge(
    child,
    mother[, .(YEAR, SAMPLE, SERIAL, mother_pernum = as.numeric(PERNUM),
               mother_person_key = person_key, mother_age = AGE_norm)],
    by = c("YEAR", "SAMPLE", "SERIAL", "mother_pernum"),
    all = TRUE,
    allow.cartesian = FALSE
  )
  links[, is_child_record := !is.na(child_pernum)]
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
    child_age < 0, "negative child age",
    child_age >= mother_age, "child age not younger than mother",
    default = NA_character_
  )]
  valid_links <- links[link_valid == TRUE]

  # Collapse valid linked ages to one row per mother while retaining raw IDs.
  age_map <- valid_links[, .(linked_ages = list(sort(child_age, decreasing = TRUE)),
                             linked_child_count = .N),
                         by = .(mother_person_key)]
  mother_rows <- merge(
    mother,
    age_map,
    by.x = "person_key", by.y = "mother_person_key", all.x = TRUE,
    sort = FALSE
  )
  mother_rows[is.na(linked_child_count), linked_child_count := 0L]
  mother_rows[is.na(linked_ages), linked_ages := list(numeric(0)), by = person_key]
  mother_rows[, a1 := vapply(linked_ages, function(z) if (length(z) >= 1L) z[1] else NA_real_, numeric(1))]
  mother_rows[, a2 := vapply(linked_ages, function(z) if (length(z) >= 2L) z[2] else NA_real_, numeric(1))]
  mother_rows[, a3 := vapply(linked_ages, function(z) if (length(z) >= 3L) z[3] else NA_real_, numeric(1))]
  mother_rows[, event_time := a2]
  mother_rows[, event_year := YEAR - event_time]
  mother_rows[, age_at_event := AGE_norm - event_time]
  mother_rows[, birth_gap := a1 - a2]
  mother_rows[, age_tie := !is.na(a1) & !is.na(a2) & a1 <= a2 |
                (!is.na(a2) & !is.na(a3) & a2 <= a3)]
  mother_rows[, nchild_link_mismatch := is.na(NCHILD_norm) |
                NCHILD_norm != linked_child_count]
  mother_rows[, fertyr_invalid_raw := !is.na(FERTYR) &
                is.na(FERTYR_status)]
  mother_rows[, fertyr_no_at_event0 := !is.na(event_time) & event_time == 0 &
                !is.na(FERTYR_status) & FERTYR_status == "no"]
  mother_rows[, fertyr_event0_complete := !is.na(event_time) & event_time == 0 &
                !is.na(FERTYR_status) & FERTYR_status == "yes"]
  mother_rows[, match_covariate_missing := Reduce(`|`, lapply(
    match_covariates,
    function(v) if (v == "AGE") is.na(AGE_norm) else is.na(get(v))
  ))]
  mother_rows[, donor_match_eligible := !match_covariate_missing & !weight_invalid]
  mother_rows[, strict_eligible :=
                !is.na(NCHILD_norm) & NCHILD_norm == 2 &
                !is.na(linked_child_count) & linked_child_count == 2L &
                !is.na(a1) & !is.na(a2) & a1 >= 1 & a1 > a2 &
                !age_tie & !is.na(age_at_event) &
                age_at_event >= age_at_event_bounds[1] &
                age_at_event <= age_at_event_bounds[2] &
                birth_gap >= 1 & !fertyr_no_at_event0 &
                !fertyr_invalid_raw & !weight_invalid &
                !match_covariate_missing]
  mother_rows[, gap_full_pre := strict_eligible & birth_gap >= full_pre_min_gap]
  mother_rows[, gap_reference := strict_eligible & birth_gap >= reference_min_gap]
  mother_rows[, event_time_outside_window := !is.na(event_time) &
                !(event_time %in% event_window)]
  mother_rows[, exclusion_flag := data.table::fcase(
    is.na(AGE_norm), "missing_mother_age",
    nchild_link_mismatch, "nchild_link_mismatch",
    NCHILD_norm < 2 | linked_child_count < 2L, "fewer_than_two_linked_children",
    is.na(a1) | is.na(a2) | a1 < 1, "no_older_child_or_missing_child_age",
    age_tie, "ambiguous_age_tie",
    is.na(age_at_event) | age_at_event < age_at_event_bounds[1] |
      age_at_event > age_at_event_bounds[2], "age_at_event_outside_band",
    fertyr_invalid_raw, "fertyr_unrecognized_raw_code",
    fertyr_no_at_event0, "fertyr_observed_no_at_event0",
    weight_invalid, "missing_or_nonpositive_perwt",
    match_covariate_missing, "missing_match_covariate",
    is.na(birth_gap) | birth_gap < 1, "nonpositive_birth_gap",
    default = NA_character_
  )]

  # Preserve raw fields and person weights. No generated wgt is assigned here.
  donor_cols <- c("person_key", "household_key", "YEAR", "SAMPLE", "SERIAL",
                  "PERNUM", "PERWT", "AGE", "AGE_norm", "NCHILD", "SEX",
                  "a1", "linked_child_count", "match_covariate_missing",
                  "weight_invalid", "donor_match_eligible", match_covariates)
  donor_cols <- unique(donor_cols[donor_cols %in% names(mother_rows)])
  donors <- mother_rows[NCHILD_norm == 1 & linked_child_count == 1L &
                          !is.na(a1) & a1 >= 0,
                        ..donor_cols]
  data.table::setnames(donors, "a1", "sole_child_age")

  anchors <- mother_rows[strict_eligible == TRUE & event_time == 0]
  post_rows <- mother_rows[strict_eligible == TRUE & event_time >= 0 &
                             event_time %in% event_window]

  # Targets are expanded from anchors only. They carry anchor covariates, with
  # AGE represented by the generated target_mother_age.
  targets <- .empty_table()
  if (nrow(anchors) > 0L) {
    neg_k <- event_window[event_window < 0L]
    targets <- anchors[, {
      ans <- data.table::rbindlist(lapply(neg_k, function(k) {
        child_age <- birth_gap + k
        if (is.na(child_age) || child_age < 0) return(NULL)
        out <- data.table::data.table(
          anchor_person_key = person_key,
          anchor_household_key = household_key,
          anchor_event_year = event_year,
          anchor_age_at_event = age_at_event,
          birth_gap = birth_gap,
          target_event_time = k,
          target_year = event_year + k,
          target_mother_age = age_at_event + k,
          target_child_age = child_age,
          gap_full_pre = gap_full_pre,
          gap_reference = gap_reference
        )
        for (v in match_covariates) {
          out[[paste0("target_", v)]] <- if (v == "AGE")
            age_at_event + k else get(v)
        }
        out
      }), fill = TRUE)
      ans
    }, by = seq_len(nrow(anchors))]
    targets[, seq_len := NULL]
  }

  audit <- data.table::data.table(
    category = c(
      "input_person_rows", "input_unique_person_keys", "missing_sex_rows",
      "missing_age_rows",
      "missing_or_nonpositive_perwt_rows", "missing_match_covariate_mothers",
      "female_valid_age_mothers", "raw_momloc_child_rows",
      "valid_momloc_links", "mother_without_link_rows",
      "unmatched_child_rows", "invalid_missing_child_age_links",
      "invalid_self_momloc_links", "invalid_child_age_not_younger_links",
      "invalid_or_unrecognized_fertyr_rows", "strict_eligible_rows",
      "strict_rows_outside_event_window", "event0_anchor_rows", "post_rows",
      "one_child_donors",
      "gap_full_pre_anchors", "gap_reference_anchors", "donor_targets"),
    count = c(
      nrow(dt), length(unique(dt$person_key)), sum(is.na(dt$SEX)),
      sum(is.na(dt$AGE_norm)), sum(dt$weight_invalid),
      sum(mother_rows$match_covariate_missing, na.rm = TRUE), nrow(mother),
      nrow(child), nrow(valid_links),
      sum(links$link_invalid_reason == "mother_without_link", na.rm = TRUE),
      sum(links$link_invalid_reason == "child_without_female_mother", na.rm = TRUE),
      sum(links$link_invalid_reason == "missing child age" & links$is_child_record,
          na.rm = TRUE),
      sum(links$link_invalid_reason == "MOMLOC self-link" & links$is_child_record,
          na.rm = TRUE),
      sum(links$link_invalid_reason == "child age not younger than mother" &
            links$is_child_record, na.rm = TRUE),
      sum(mother_rows$fertyr_invalid_raw, na.rm = TRUE),
      sum(mother_rows$strict_eligible, na.rm = TRUE),
      sum(mother_rows$strict_eligible & mother_rows$event_time_outside_window,
          na.rm = TRUE), nrow(anchors), nrow(post_rows), nrow(donors),
      sum(anchors$gap_full_pre, na.rm = TRUE),
      sum(anchors$gap_reference, na.rm = TRUE), nrow(targets)))

  list(
    input = dt,
    mother_rows = mother_rows,
    links = links,
    anchors = anchors,
    post_rows = post_rows,
    one_child_donors = donors,
    donor_targets = targets,
    audit = audit,
    config = list(
      fertyr_codes = fertyr_codes,
      female_code = female_code,
      age_missing_codes = age_missing_codes,
      event_window = event_window,
      age_at_event_bounds = age_at_event_bounds,
      match_covariates = match_covariates,
      full_pre_min_gap = full_pre_min_gap,
      reference_min_gap = reference_min_gap
    )
  )
}
