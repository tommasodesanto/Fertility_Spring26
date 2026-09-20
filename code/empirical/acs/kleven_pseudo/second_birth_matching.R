# Exact-cell donor matching for the coresident second-birth proxy.
#
# This is the smallest executable adapter to the vendor matching contract:
# year, mother age, transformed author covariates, and sole-child age are
# exact; matching is with replacement and all ties are retained. It does not
# coarsen cells, match on outcomes, or assign weights to target rows.

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("second_birth_matching.R requires the data.table package")
}

.sbm_cell <- function(d, cols) {
  vals <- lapply(cols, function(x) {
    z <- as.character(d[[x]])
    z[is.na(z)] <- "<missing>"
    z
  })
  do.call(paste, c(vals, sep = "\034"))
}

.sbm_require <- function(d, cols, label) {
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(sprintf("%s missing columns: %s", label,
                                paste(miss, collapse = ", ")))
}

#' Match event-0 donor targets to one-child donor rows in exact cells.
#'
#' `targets` is the `donor_targets` table from build_second_birth_proxy().
#' `donors` is its `one_child_donors` table. `match_covariates` must name the
#' transformed author cells (`gender.num`, `edlevel.num`, `marst.num`,
#' `race.num`, `statefip.num`) materialized in both tables. Raw extract fields
#' SEX/EDUC/MARST/RACE/STATEFIP are rejected because their source-to-author
#' recodes, especially race/ethnicity, are not interchangeable. AGE is
#' represented by target_mother_age and AGE_norm and is omitted if supplied.
#'
#' The returned donor weights follow matching.R: each tied donor receives the
#' reciprocal of the number of controls tied to that target row, summed only
#' within a donor-by-pseudo-event cell; `wgt = PERWT * wgt_match`. Target rows
#' retain their original survey weight outside this function.
second_birth_match_exact <- function(
    targets, donors, match_covariates,
    full_pre_times = -5:-1, reference_time = -2) {
  if (!data.table::is.data.table(targets) || !data.table::is.data.table(donors))
    stop("targets and donors must be data.table objects")
  if (missing(match_covariates) || !length(match_covariates))
    stop("match_covariates must be supplied explicitly")
  if (anyDuplicated(match_covariates)) stop("match_covariates must be unique")
  raw_names <- intersect(match_covariates, c("SEX", "EDUC", "MARST", "RACE", "STATEFIP"))
  if (length(raw_names))
    stop(sprintf("raw covariates require explicit author-cell recodes: %s",
                 paste(raw_names, collapse = ", ")))
  if ("AGE" %in% match_covariates) match_covariates <- setdiff(match_covariates, "AGE")
  if (any(!is.finite(full_pre_times)) || any(full_pre_times != as.integer(full_pre_times)) ||
      !all(full_pre_times < 0L) || !reference_time < 0L)
    stop("full_pre_times and reference_time must be negative integer event times")
  full_pre_times <- as.integer(full_pre_times); reference_time <- as.integer(reference_time)

  target_cov <- paste0("target_", match_covariates)
  tcols <- c("anchor_person_key", "anchor_household_key", "target_event_time", "target_year",
             "target_mother_age", "target_child_age", target_cov)
  dcols <- c("person_key", "household_key", "PERWT", "YEAR", "AGE_norm", "sole_child_age",
             match_covariates)
  .sbm_require(targets, tcols, "targets")
  .sbm_require(donors, dcols, "donors")
  if (any(!is.finite(as.numeric(targets$target_event_time))) ||
      any(as.numeric(targets$target_event_time) >= 0))
    stop("targets must contain negative event times from event-0 anchors only")
  if (anyDuplicated(paste(targets$anchor_person_key, targets$target_event_time, sep = "\034")))
    stop("duplicate anchor/event target rows")
  if (anyDuplicated(donors$person_key)) stop("donor person_key must be unique")

  targets <- data.table::copy(targets)
  donors <- data.table::copy(donors)
  target_key_cols <- c("target_year", "target_mother_age", "target_child_age", target_cov)
  donor_key_cols <- c("YEAR", "AGE_norm", "sole_child_age", match_covariates)
  targets[, .sbm_cell := .sbm_cell(.SD, target_key_cols), .SDcols = target_key_cols]
  donors[, .sbm_cell := .sbm_cell(.SD, donor_key_cols), .SDcols = donor_key_cols]

  donor_bad <- is.na(donors$PERWT) | !is.finite(as.numeric(donors$PERWT)) |
    as.numeric(donors$PERWT) <= 0 |
    Reduce(`|`, lapply(donor_key_cols, function(x) is.na(donors[[x]])))
  target_bad <- Reduce(`|`, lapply(target_key_cols, function(x) is.na(targets[[x]])))
  donors[, donor_match_excluded := donor_bad]
  targets[, target_match_excluded := target_bad]
  d_ok <- donors[donor_match_excluded == FALSE]
  t_ok <- targets[target_match_excluded == FALSE]

  # Merge retains every target-donor tie. No nearest-cell fallback is allowed.
  links <- merge(
    t_ok[, .(anchor_person_key, anchor_household_key, target_event_time, target_year,
             target_mother_age, target_child_age, .sbm_cell)],
    d_ok[, .(donor_person_key = person_key, donor_household_key = household_key,
             donor_PERWT = PERWT, .sbm_cell)],
    by = ".sbm_cell", allow.cartesian = TRUE, sort = FALSE
  )
  if (nrow(links)) {
    # Matching::Match normalizes all tied controls within each treated row so
    # their weights sum to one. This differs from a raw link count.
    links[, wgt_match := 1 / .N, by = .(anchor_person_key, target_event_time)]
    data.table::setorder(links, anchor_person_key, target_event_time, donor_person_key)
    donor_weights <- links[, .(
      PERWT = donor_PERWT[1L], wgt_match = sum(wgt_match),
      n_target_rows = .N),
      by = .(donor_person_key, donor_household_key, target_event_time,
             target_year, target_mother_age, target_child_age)]
    donor_weights[, `:=`(person_key = donor_person_key,
                         wgt_original = as.numeric(PERWT),
                         wgt = as.numeric(PERWT) * as.numeric(wgt_match))]
  } else {
    donor_weights <- data.table::data.table(
      donor_person_key = character(), donor_household_key = character(),
      target_event_time = integer(), target_year = integer(),
      target_mother_age = numeric(), target_child_age = numeric(),
      PERWT = numeric(), wgt_match = numeric(), n_target_rows = integer(),
      person_key = character(),
      wgt_original = numeric(), wgt = numeric())
  }

  # Keyed cell counts avoid scanning the donor table for every target row.
  donor_cell_counts <- d_ok[, .(candidate_donors = .N), by = .sbm_cell]
  support <- merge(targets[, .(anchor_person_key, target_event_time, .sbm_cell,
                               target_match_excluded)], donor_cell_counts,
                   by = ".sbm_cell", all.x = TRUE, sort = FALSE)
  support[is.na(candidate_donors) | target_match_excluded == TRUE, candidate_donors := 0L]
  support[, .sbm_cell := NULL]
  support[, has_donor := candidate_donors > 0L]
  anchor_support <- support[, .(
    full_pre_expected = length(full_pre_times),
    full_pre_observed = sum(target_event_time %in% full_pre_times & has_donor),
    reference_has_donor = any(target_event_time == reference_time & has_donor),
    gap_full_pre = any(target_event_time == min(full_pre_times)),
    gap_reference = any(target_event_time == reference_time)
  ), by = anchor_person_key]
  anchor_support[, `:=`(
    full_pre_supported = gap_full_pre & full_pre_observed == full_pre_expected,
    reference_supported = gap_reference & reference_has_donor)]

  list(
    links = links,
    donor_weights = donor_weights,
    target_support = support,
    anchor_support = anchor_support,
    audit = list(
      target_rows = nrow(targets), target_rows_eligible = nrow(t_ok),
      donor_rows = nrow(donors), donor_rows_eligible = nrow(d_ok),
      target_rows_excluded = sum(target_bad), donor_rows_excluded = sum(donor_bad),
      matched_link_rows = nrow(links),
      full_pre_supported_anchors = sum(anchor_support$full_pre_supported),
      reference_supported_anchors = sum(anchor_support$reference_supported),
      matching = "exact cells, with replacement, all ties",
      weight_formula = "donor wgt = PERWT * sum(1/n_tied_controls); target weights unchanged"
    ),
    config = list(match_covariates = match_covariates,
                  full_pre_times = full_pre_times,
                  reference_time = reference_time)
  )
}
