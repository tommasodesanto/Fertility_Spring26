# National pooled ACS first-birth housing estimator.
#
# The adapter (the caller) owns matching and source-key coalescing.  This file
# only consumes a compact matched panel.  Housing is observed only for true ACS
# rows; CPS rows remain in the audit but cannot receive an imputed outcome.

nfh_stop <- function(...) stop(paste0("national_first_birth_housing: ",
                                      paste0(..., collapse = "")), call. = FALSE)
nfh_require <- function(ok, msg) if (!isTRUE(ok)) nfh_stop(msg)
nfh_cols <- function(x, cols, label = "panel") {
  miss <- setdiff(cols, names(x))
  if (length(miss)) nfh_stop(label, " missing columns: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

# Matching adapters can leave both a merge-suffixed field and a plain field.
# The suffixed value is authoritative; an existing plain value is not silently
# substituted for a present .x value.
nfh_resolve <- function(x, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) {
    if (required) nfh_stop("no column for ", label, " (tried ", paste(candidates, collapse = ", "), ")")
    return(NA_character_)
  }
  hit[[1L]]
}

nfh_as_num <- function(x, label) {
  z <- suppressWarnings(as.numeric(as.character(x)))
  if (any(!is.na(x) & !is.finite(z))) nfh_stop(label, " contains nonnumeric values")
  z
}

nfh_json_line <- function(path, object) {
  if (!is.null(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!requireNamespace("jsonlite", quietly = TRUE)) nfh_stop("jsonlite is required for failure receipts")
    jsonlite::write_json(object, path, auto_unbox = TRUE, pretty = TRUE)
  }
}

nfh_safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

nfh_weighted_mean <- function(y, w) {
  ok <- is.finite(y) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(y[ok] * w[ok]) / sum(w[ok])
}

nfh_ess <- function(w) {
  ok <- is.finite(w) & w > 0
  if (!any(ok)) return(0)
  (sum(w[ok])^2) / sum(w[ok]^2)
}

nfh_code_housing <- function(panel, rooms_col, ownership_col, bedrooms_col = NA_character_) {
  rooms_raw <- nfh_as_num(panel[[rooms_col]], rooms_col)
  # IPUMS ACS: 0 and 28 are non-room/unknown; only 1:27 and 30 are valid.
  rooms_valid <- !is.na(rooms_raw) & rooms_raw %in% c(1:27, 30)
  rooms9 <- rep(NA_real_, nrow(panel)); rooms9[rooms_valid] <- pmin(rooms_raw[rooms_valid], 9)
  own_raw <- nfh_as_num(panel[[ownership_col]], ownership_col)
  own_valid <- !is.na(own_raw) & own_raw %in% c(1, 2)
  ownership_lw <- rep(NA_real_, nrow(panel)); ownership_lw[own_valid] <- ifelse(own_raw[own_valid] == 1, 1, 0)
  out <- panel
  out$rooms_raw <- rooms_raw; out$rooms9 <- rooms9
  out$ownership_raw <- own_raw; out$ownership_lw <- ownership_lw
  out$rooms_unknown_code <- !is.na(rooms_raw) & rooms_raw == 28
  out$rooms_missing_code <- !is.na(rooms_raw) & rooms_raw == 0
  out$rooms_invalid_code <- !is.na(rooms_raw) & !rooms_valid
  out$ownership_unknown_code <- !is.na(own_raw) & !(own_raw %in% c(0, 1, 2, 9))
  out$ownership_missing_code <- !is.na(own_raw) & own_raw %in% c(0, 9)
  out$ownership_invalid_code <- !is.na(own_raw) & !own_valid
  if (!is.na(bedrooms_col)) {
    br <- nfh_as_num(panel[[bedrooms_col]], bedrooms_col)
    bv <- !is.na(br) & br %in% 1:22
    bedrooms5 <- rep(NA_real_, nrow(panel)); bedrooms5[bv] <- pmin(pmax(br[bv] - 1, 0), 5)
    out$bedrooms_raw <- br; out$bedrooms5 <- bedrooms5
    out$bedrooms_unknown_code <- !is.na(br) & br == 23
    out$bedrooms_missing_code <- !is.na(br) & br == 0
    out$bedrooms_invalid_code <- !is.na(br) & !bv
  }
  out
}

nfh_event_term <- function(coef_names, event, ref = -2) {
  if (identical(as.character(event), as.character(ref))) return(NA_character_)
  term <- paste0("event_time::", as.character(event))
  hit <- which(coef_names == term)
  if (length(hit) != 1L) nfh_stop("missing or ambiguous event coefficient ", term,
                                    "; coefficient names were: ", paste(coef_names, collapse = ", "))
  term
}

nfh_contrast <- function(reg, event_times, ref, outcome, specification, variance_label) {
  b <- stats::coef(reg); V <- as.matrix(stats::vcov(reg))
  terms <- vapply(event_times, nfh_event_term, character(1), coef_names = names(b), ref = ref)
  # vapply returns NA_character_ for the reference; the contrast never uses it.
  get_term <- function(tt) nfh_event_term(names(b), tt, ref)
  t3 <- get_term(3); tm1 <- get_term(-1)
  if (is.na(t3) || is.na(tm1)) nfh_stop("+3/-1 contrast includes the reference event")
  L <- setNames(numeric(length(b)), names(b)); L[[t3]] <- 1; L[[tm1]] <- -1
  variance <- as.numeric(t(L) %*% V %*% L)
  tol <- 1e-10 * max(1, max(abs(V), na.rm = TRUE))
  if (!is.finite(variance) || variance < -tol)
    nfh_stop("negative or nonfinite contrast variance for ", outcome, "/", specification,
             ": ", format(variance, scientific = TRUE))
  se <- sqrt(max(0, variance))
  est <- unname(b[[t3]] - b[[tm1]])
  data.frame(outcome = outcome, specification = specification, contrast = "+3_minus_-1",
             estimate = est, std.error = se, conf.low = est - qnorm(.975) * se,
             conf.high = est + qnorm(.975) * se, variance = variance,
             variance_label = variance_label, stringsAsFactors = FALSE)
}

nfh_event_curve <- function(reg, data, outcome, specification, event_times, ref, variance_label) {
  b <- stats::coef(reg); V <- as.matrix(stats::vcov(reg)); out <- vector("list", length(event_times))
  for (i in seq_along(event_times)) {
    ev <- event_times[[i]]; term <- nfh_event_term(names(b), ev, ref)
    if (is.na(term)) { est <- 0; variance <- 0; term_out <- NA_character_ }
    else {
      est <- unname(b[[term]]); variance <- unname(V[term, term]); term_out <- term
      tol <- 1e-10 * max(1, max(abs(V), na.rm = TRUE))
      if (!is.finite(est) || !is.finite(variance) || variance < -tol)
        nfh_stop("invalid event coefficient/variance for ", outcome, "/", specification, "/", ev)
      variance <- max(0, variance)
    }
    se <- sqrt(variance)
    out[[i]] <- data.frame(outcome = outcome, specification = specification,
      event_time = ev, term = term_out, reference = identical(as.character(ev), as.character(ref)),
      estimate = est, std.error = se, conf.low = est - qnorm(.975) * se,
      conf.high = est + qnorm(.975) * se, variance = variance,
      variance_label = variance_label, nobs = stats::nobs(reg),
      clusters = length(unique(data$source_hh_cluster)), stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

estimate_national_first_birth_housing <- function(
    panel, output_dir = NULL, checkpoint = NULL, event_times = as.integer(-5:10),
    ref = -2L, women_only = TRUE, outcomes = NULL, weight_col = "wgt",
    rooms_col = NULL, ownership_col = NULL, bedrooms_col = NULL,
    source_origin_col = NULL, from_cps_col = NULL) {
  if (!is.data.frame(panel)) nfh_stop("panel must be a data.frame")
  if (!requireNamespace("fixest", quietly = TRUE)) nfh_stop("fixest is required")
  if (!is.null(output_dir) && !requireNamespace("jsonlite", quietly = TRUE)) nfh_stop("jsonlite is required when output_dir is used")
  event_times <- as.integer(event_times)
  if (!identical(as.integer(ref), -2L)) nfh_stop("ref must be -2")
  if (!length(event_times) || anyDuplicated(event_times) || !(ref %in% event_times) ||
      !all(c(-1L, 3L) %in% event_times)) nfh_stop("event_times must include unique -1, 3 and reference -2")
  if (is.null(rooms_col)) rooms_col <- nfh_resolve(panel, c("rooms_raw", "ROOMS_RAW", "raw_rooms"), "raw ROOMS")
  if (is.null(ownership_col)) ownership_col <- nfh_resolve(panel, c("ownershp_raw", "OWNERSHP_RAW", "ownership_raw", "raw_ownership"), "raw OWNERSHP")
  if (is.null(bedrooms_col)) bedrooms_col <- nfh_resolve(panel, c("bedrooms_raw", "BEDROOMS_RAW", "raw_bedrooms"), "raw BEDROOMS", required = FALSE)
  age_col <- nfh_resolve(panel, c("age_factor", "author_age_factor"), "author age factor")
  year_col <- nfh_resolve(panel, c("doiy_factor", "author_year_factor"), "author interview-year factor")
  state_col <- nfh_resolve(panel, c("statefip", "statefip.x"), "author state FIPS")
  gender_col <- nfh_resolve(panel, c("gender", "gender.x"), "gender")
  event_col <- nfh_resolve(panel, c("event_time", "t_es_lw", "t_es"), "event time")
  nfh_cols(panel, c(weight_col, age_col, year_col, state_col, gender_col, event_col), "panel")
  # Metadata resolution is explicitly .x first, then plain, to match the
  # national matcher adapter's coalesced fields.
  origin_col <- if (is.null(source_origin_col)) nfh_resolve(panel, c("source_origin.x", "source_origin"), "source origin") else source_origin_col
  cps_col <- if (is.null(from_cps_col)) nfh_resolve(panel, c("from_cps.x", "from_cps"), "from_cps") else from_cps_col
  year_source_col <- nfh_resolve(panel, c("YEAR.x", "source_year.x", "YEAR", "source_year", "source_YEAR", "source_doiy"), "source year")
  sample_source_col <- nfh_resolve(panel, c("SAMPLE.x", "source_sample.x", "SAMPLE", "source_sample", "sample"), "source SAMPLE")
  serial_source_col <- nfh_resolve(panel, c("SERIAL.x", "source_serial.x", "SERIAL", "source_serial", "serial"), "source SERIAL")
  source_state_col <- nfh_resolve(panel, c("STATEFIP.x", "source_statefip.x", "STATEFIP", "source_statefip", "source_state"), "source state", required = FALSE)
  source_age_col <- nfh_resolve(panel, c("AGE.x", "source_age.x", "AGE", "source_age", "source_AGE"), "source age", required = FALSE)
  if (any(c(age_col, year_col) %in% c("source_age", "source_year", "source_AGE", "source_YEAR")))
    nfh_stop("regression FE must use author age_factor/doiy_factor, not source age/year")
  weight <- nfh_as_num(panel[[weight_col]], weight_col)
  if (any(!is.na(weight) & weight <= 0)) nfh_stop(weight_col, " must be strictly positive when observed")
  event <- nfh_as_num(panel[[event_col]], event_col)
  if (any(!is.na(event) & !(event %in% event_times))) nfh_stop("event time outside requested window")
  origin <- as.character(panel[[origin_col]])
  cps <- if (is.logical(panel[[cps_col]])) as.integer(panel[[cps_col]]) else
    suppressWarnings(as.numeric(as.character(panel[[cps_col]])))
  if (any(is.na(origin) | !(origin %in% c("ACS", "CPS")))) nfh_stop("source origin must be ACS/CPS")
  if (any(is.na(cps) | !(cps %in% c(0, 1)))) nfh_stop("from_cps must contain only 0/1")
  true_acs <- origin == "ACS" & cps == 0
  if (!any(true_acs)) nfh_stop("no true ACS rows")
  dat <- panel
  dat$.row_id <- seq_len(nrow(dat)); dat$.weight <- weight; dat$.event <- event
  dat$source_year <- nfh_as_num(dat[[year_source_col]], year_source_col)
  dat$source_sample <- nfh_as_num(dat[[sample_source_col]], sample_source_col)
  dat$source_serial <- nfh_as_num(dat[[serial_source_col]], serial_source_col)
  if (any(!is.finite(dat$source_year[true_acs]) | !is.finite(dat$source_sample[true_acs]) |
          !is.finite(dat$source_serial[true_acs])))
    nfh_stop("true ACS rows require finite source YEAR/SAMPLE/SERIAL")
  dat$source_hh_cluster <- NA_character_
  dat$source_hh_cluster[true_acs] <- paste(origin[true_acs], dat$source_year[true_acs],
                                            dat$source_sample[true_acs], dat$source_serial[true_acs], sep = ":")
  dat <- nfh_code_housing(dat, rooms_col, ownership_col, bedrooms_col)
  # Preserve the author's canonical FE objects under their explicit names.
  # Source age/year fields above remain lineage metadata only.
  dat$age_factor <- dat[[age_col]]; dat$doiy_factor <- dat[[year_col]]
  dat$statefip <- dat[[state_col]]; dat$gender_value <- as.character(dat[[gender_col]])
  if (isTRUE(women_only)) {
    base_gender <- dat$gender_value == "Women"
    if (!any(base_gender, na.rm = TRUE)) nfh_stop("women_only=TRUE but no gender == Women rows")
  } else base_gender <- rep(TRUE, nrow(dat))
  base <- true_acs & base_gender & is.finite(dat$.weight) & dat$.weight > 0 &
    !is.na(dat$.event) & dat$.event %in% event_times
  base_n <- sum(base)
  if (!base_n) nfh_stop("no true ACS women in weighted event window")
  if (is.null(outcomes)) outcomes <- c("rooms9", "ownership_lw", if (!is.na(bedrooms_col)) "bedrooms5")
  outcomes <- unique(as.character(outcomes)); bad <- setdiff(outcomes, names(dat)); if (length(bad)) nfh_stop("unknown outcomes: ", paste(bad, collapse = ", "))
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  emit <- function(x) { if (!is.null(checkpoint)) checkpoint(x); invisible(NULL) }
  fit_specs <- c(full = "state + age + year", event_only = "none", age_only = "age", state_year = "state + year")
  formulas <- list(
    full = "OUTCOME ~ i(event_time, ref = -2) | statefip + age_factor + doiy_factor",
    event_only = "OUTCOME ~ i(event_time, ref = -2)",
    age_only = "OUTCOME ~ i(event_time, ref = -2) | age_factor",
    state_year = "OUTCOME ~ i(event_time, ref = -2) | statefip + doiy_factor")
  fits <- list(); curves <- list(); contrasts <- list(); status <- list(); kk <- 0L
  for (outcome in outcomes) {
    od <- dat[base & is.finite(dat[[outcome]]), , drop = FALSE]
    if (!nrow(od)) nfh_stop("no observed true ACS rows for ", outcome)
    od$event_time <- od$.event
    od$statefip <- factor(od$statefip)
    od$age_factor <- factor(od$age_factor)
    od$doiy_factor <- factor(od$doiy_factor)
    # Exact support check avoids accepting silently omitted event coefficients.
    if (any(table(factor(od$event_time, levels = event_times)) == 0))
      nfh_stop("missing event support for ", outcome)
    for (spec in names(formulas)) {
      nm <- paste(outcome, spec, sep = "::")
      emit(list(stage = "fit_start", outcome = outcome, specification = spec,
                nobs = nrow(od), base_nobs = base_n))
      one <- tryCatch({
        fml <- stats::as.formula(sub("OUTCOME", outcome, formulas[[spec]], fixed = TRUE))
        # All specifications use the same outcome-specific rows and author weight.
        reg <- fixest::feols(fml, data = od, weights = ~.weight,
                             vcov = ~source_hh_cluster)
        b <- stats::coef(reg); V <- as.matrix(stats::vcov(reg))
        expected <- setdiff(event_times, ref)
        for (ev in expected) if (!(paste0("event_time::", ev) %in% names(b)))
          nfh_stop("missing event coefficient for ", outcome, "/", spec, "/", ev)
        fit_data <- od[fixest::obs(reg), , drop = FALSE]
        cv <- nfh_event_curve(reg, fit_data, outcome, spec, event_times, ref, "source_hh_cluster")
        ct <- nfh_contrast(reg, event_times, ref, outcome, spec, "source_hh_cluster")
        fit <- list(outcome = outcome, specification = spec, formula = paste(deparse(fml), collapse = " "),
          coefficients = b, full_vcov = V, nobs = stats::nobs(reg), source_hh_clusters = length(unique(fit_data$source_hh_cluster)),
          event_curve = cv, contrast = ct, model = reg)
        if (!is.null(output_dir)) saveRDS(fit, file.path(output_dir, paste0("checkpoint_", nfh_safe_name(nm), ".rds")))
        emit(list(stage = "fit_complete", outcome = outcome, specification = spec,
                  nobs = fit$nobs, clusters = fit$source_hh_clusters, contrast = ct))
        fit
      }, error = function(e) {
        failure <- list(stage = "fit_failure", outcome = outcome, specification = spec, error = conditionMessage(e))
        if (!is.null(output_dir)) nfh_json_line(file.path(output_dir, paste0("failure_", nfh_safe_name(nm), ".json")), failure)
        emit(failure); stop(e)
      })
      fits[[nm]] <- one; kk <- kk + 1L; curves[[kk]] <- one$event_curve; contrasts[[kk]] <- one$contrast
      status[[kk]] <- data.frame(outcome = outcome, specification = spec, nobs = one$nobs,
        source_hh_clusters = one$source_hh_clusters, base_nobs = base_n,
        nobs_differs_from_fullspec = NA, stringsAsFactors = FALSE)
    }
  }
  status <- do.call(rbind, status); curves <- do.call(rbind, curves); contrasts <- do.call(rbind, contrasts)
  full_rows <- status$specification == "full"
  full_n <- setNames(status$nobs[full_rows], status$outcome[full_rows])
  status$nobs_differs_from_fullspec <- mapply(function(o, n) !identical(as.integer(n), as.integer(full_n[[o]])), status$outcome, status$nobs)
  raw_rows <- list(); counts <- list(); k <- 0L
  for (o in outcomes) for (ev in event_times) {
    z <- dat[base & dat$.event == ev & is.finite(dat[[o]]), , drop = FALSE]
    k <- k + 1L
    raw_rows[[k]] <- data.frame(outcome = o, event_time = ev, nobs = nrow(z),
      weighted_mean = nfh_weighted_mean(z[[o]], z$.weight), weight_sum = sum(z$.weight),
      source_hh_clusters = length(unique(z$source_hh_cluster)), stringsAsFactors = FALSE)
    counts[[k]] <- data.frame(outcome = o, event_time = ev, nobs = nrow(z),
      clusters = length(unique(z$source_hh_cluster)), event_ess = nfh_ess(z$.weight),
      weight_sum = sum(z$.weight), missing_outcome = sum(base & dat$.event == ev & is.na(dat[[o]])),
      stringsAsFactors = FALSE)
  }
  raw_baselines <- do.call(rbind, raw_rows); counts <- do.call(rbind, counts)
  if (!is.null(output_dir)) {
    utils::write.csv(curves, file.path(output_dir, "national_event_curves.csv"), row.names = FALSE)
    utils::write.csv(contrasts, file.path(output_dir, "national_contrasts.csv"), row.names = FALSE)
    utils::write.csv(raw_baselines, file.path(output_dir, "national_raw_baselines.csv"), row.names = FALSE)
    utils::write.csv(counts, file.path(output_dir, "national_counts_event_ess.csv"), row.names = FALSE)
    utils::write.csv(status, file.path(output_dir, "national_fit_status.csv"), row.names = FALSE)
    png(file.path(output_dir, "national_housing_event_curves.png"), width = 1800, height = 1100, res = 140)
    par(mfrow = c(length(outcomes), 1), mar = c(3, 4, 2, 1))
    for (o in outcomes) {
      z <- curves[curves$outcome == o & curves$specification == "full", ]
      display_scale <- if (identical(o, "ownership_lw")) 100 else 1
      display_label <- if (identical(o, "rooms9")) "Rooms (cap9)" else
        if (identical(o, "ownership_lw")) "Ownership (pp)" else "Bedrooms (cap5)"
      plot(z$event_time, display_scale * z$estimate, type = "b", pch = 16,
           ylim = display_scale * range(c(z$conf.low, z$conf.high), finite = TRUE),
           xlab = "Event time", ylab = display_label, main = paste("National ACS", display_label, "(full FE)"))
      zi <- is.finite(z$conf.low) & is.finite(z$conf.high) & (z$conf.high > z$conf.low)
      if (any(zi)) segments(z$event_time[zi], display_scale * z$conf.low[zi],
                             z$event_time[zi], display_scale * z$conf.high[zi])
      abline(v = ref, lty = 2); abline(h = 0, lty = 3)
    }
    dev.off()
  }
  list(status = "ESTIMATION_COMPLETE_DIAGNOSTIC", fits = fits, curves = curves,
       contrasts = contrasts, summary = contrasts, raw_baselines = raw_baselines,
       counts = counts, fit_status = status, data = dat, metadata = list(
         reference_event = ref, event_times = event_times, women_only = women_only,
         author_weight = weight_col, author_age_factor = age_col, author_year_factor = year_col,
         full_formula = formulas$full, source_cluster = "source_origin:source_YEAR:source_SAMPLE:source_SERIAL",
         source_fields = c(year_source_col, sample_source_col, serial_source_col),
         housing_validity = "ROOMS 1:27,30 capped at 9; OWNERSHP 1 owner/2 renter; no CPS imputation",
         base_nobs = base_n, actual_nobs_differences = status[status$nobs_differs_from_fullspec, , drop = FALSE]))
}
