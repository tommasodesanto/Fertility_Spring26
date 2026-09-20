# First-birth housing estimation driver (source-key-gated, diagnostic only).
#
# The future launcher must pass the already matched author panel, a verified
# NE housing key packet, its audit manifest, and an explicit coding_config.
# This file never loads raw ACS data, reruns matching, imputes CPS outcomes, or
# changes the vendor files.  Its estimator follows the vendor state-event
# specification at the LEVEL outcome: for each gender, it estimates
#
#   y_it ~ state_i + t_es_t:state_i | age_factor_i + doiy_factor_i,
#
# with author matching weights and state-household clustered standard errors.
# The reference event -2 is reported as zero; all other event coefficients must
# be present and finite.  No employment pre-earnings normalization is applied.

fbh_stop <- function(...) {
  stop(paste0("first_birth_housing: ", paste0(..., collapse = "")), call. = FALSE)
}

fbh_require <- function(ok, msg) if (!isTRUE(ok)) fbh_stop(msg)

fbh_require_columns <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) fbh_stop(label, " missing columns: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

fbh_key <- function(x, cols, label, allow_missing = FALSE) {
  fbh_require_columns(x, cols, label)
  vals <- lapply(cols, function(nm) {
    z <- x[[nm]]
    if (is.factor(z)) z <- as.character(z)
    if (is.numeric(z) && any(!is.finite(z) & !is.na(z)))
      fbh_stop(label, " has non-finite key component in ", nm)
    z <- as.character(z)
    if (!allow_missing && any(is.na(z) | !nzchar(z)))
      fbh_stop(label, " has missing/empty key component in ", nm)
    z
  })
  out <- Reduce(function(a, b) paste(a, b, sep = ":"), vals)
  if (any(!nzchar(out))) fbh_stop(label, " produced an empty key")
  out
}

fbh_scalar_true <- function(x, label) {
  if (!is.logical(x) || length(x) != 1L || is.na(x) || !x)
    fbh_stop(label, " must be TRUE")
  invisible(TRUE)
}

fbh_validate_manifest <- function(audit_manifest, source_key_cols) {
  if (!is.list(audit_manifest)) fbh_stop("audit_manifest must be a list")
  fbh_scalar_true(audit_manifest$verified, "audit_manifest$verified")
  if (!identical(as.character(audit_manifest$status), "PASS"))
    fbh_stop("audit_manifest$status must be PASS")
  if (!identical(as.character(audit_manifest$key_columns), as.character(source_key_cols)))
    fbh_stop("audit manifest key columns differ from the requested source key")
  fbh_scalar_true(audit_manifest$source_key_unique, "audit_manifest$source_key_unique")
  fbh_scalar_true(audit_manifest$overlap_verified, "audit_manifest$overlap_verified")
  invisible(TRUE)
}

# Join is intentionally many-panel-rows to one source row.  Repeated matched
# rows are retained; only the raw source key must be unique.
join_first_birth_housing <- function(panel, source_housing, audit_manifest,
                                      panel_key_cols = c("YEAR", "SAMPLE", "SERIAL", "PERNUM"),
                                      source_key_cols = panel_key_cols,
                                      source_columns = c(ROOMS = "ROOMS", BEDROOMS = "BEDROOMS",
                                                         OWNERSHP = "OWNERSHP"),
                                      labor_outcome = "emp_lw", weight_col = "wgt") {
  provenance <- c("src_key", "source_origin", "from_cps")
  fbh_require_columns(panel, c(panel_key_cols, provenance, labor_outcome, weight_col,
                               "t_es_lw"), "panel")
  if (!is.numeric(panel[[weight_col]]) && !is.integer(panel[[weight_col]]))
    fbh_stop("author weight ", weight_col, " must be numeric")
  fbh_require_columns(source_housing, c(source_key_cols, unname(source_columns)), "source_housing")
  if (length(source_columns) != 3L || is.null(names(source_columns)) ||
      !all(c("ROOMS", "BEDROOMS", "OWNERSHP") %in% names(source_columns)))
    fbh_stop("source_columns must name ROOMS, BEDROOMS, and OWNERSHP")
  fbh_validate_manifest(audit_manifest, source_key_cols)

  origin <- as.character(panel$source_origin)
  raw_cps <- panel$from_cps
  if (is.logical(raw_cps)) {
    cps <- raw_cps
  } else {
    cps_num <- suppressWarnings(as.numeric(raw_cps))
    if (any(is.na(cps_num) | !(cps_num %in% c(0, 1))))
      fbh_stop("from_cps must contain only 0/1 or logical values")
    cps <- cps_num == 1
  }
  if (any(is.na(origin) | is.na(cps) | !(origin %in% c("ACS", "CPS"))))
    fbh_stop("panel provenance has missing or unrecognized source_origin/from_cps")
  if (any((origin == "CPS") != cps))
    fbh_stop("panel source_origin and from_cps disagree")
  eligible <- origin == "ACS" & !cps

  skey <- fbh_key(source_housing, source_key_cols, "source_housing")
  if (anyDuplicated(skey)) fbh_stop("source housing key is duplicated; refusing ambiguous join")
  pkey <- rep(NA_character_, nrow(panel))
  if (any(eligible)) pkey[eligible] <- fbh_key(panel[eligible, , drop = FALSE],
                                               panel_key_cols, "ACS panel source key")
  m <- match(pkey, skey)

  out <- panel
  out$first_birth_row_id <- seq_len(nrow(out))
  canonical <- c(ROOMS = "rooms", BEDROOMS = "bedrooms", OWNERSHP = "ownership")
  for (nm in names(source_columns)) out[[paste0("raw_", canonical[[nm]])]] <- NA_real_
  for (nm in names(source_columns)) {
    vals <- source_housing[[unname(source_columns[[nm]])]]
    z <- rep(NA_real_, nrow(out))
    hit <- eligible & !is.na(m)
    z[hit] <- suppressWarnings(as.numeric(vals[m[hit]]))
    out[[paste0("raw_", canonical[[nm]])]] <- z
  }
  out$housing_join_status <- ifelse(cps, "cps_missing_outcome",
                                    ifelse(!eligible, "non_acs_excluded",
                                           ifelse(is.na(m), "acs_source_unmatched", "acs_source_matched")))
  out$housing_source_key <- pkey
  cluster_cols <- setdiff(source_key_cols, "PERNUM")
  if (!length(cluster_cols)) fbh_stop("source key must include a household cluster component")
  cluster_key <- rep(NA_character_, nrow(out))
  if (any(eligible & !is.na(m))) cluster_key[eligible & !is.na(m)] <-
    fbh_key(panel[eligible & !is.na(m), , drop = FALSE], cluster_cols, "ACS household cluster key")
  out$source_household_cluster <- ifelse(eligible & !is.na(m),
    paste(origin, cluster_key, sep = ":"), NA_character_)

  # A join may append housing fields only.  These fields must be byte-for-byte
  # unchanged, including author weights, event time, and labor outcome.
  if (!identical(out$first_birth_row_id, seq_len(nrow(panel))))
    fbh_stop("source join changed row order")
  for (nm in c(weight_col, "t_es_lw", labor_outcome)) {
    if (!isTRUE(all.equal(out[[nm]], panel[[nm]], check.attributes = FALSE)))
      fbh_stop("source join changed protected field ", nm)
  }
  status_tab <- table(out$housing_join_status, useNA = "ifany")
  source_clusters <- out$source_household_cluster[!is.na(out$source_household_cluster)]
  join_audit <- data.frame(
    panel_rows = nrow(panel), acs_rows = sum(eligible), cps_rows = sum(cps),
    matched_acs_rows = sum(out$housing_join_status == "acs_source_matched"),
    unmatched_acs_rows = sum(out$housing_join_status == "acs_source_unmatched"),
    cps_missing_outcome_rows = sum(out$housing_join_status == "cps_missing_outcome"),
    repeated_source_household_clusters = sum(duplicated(source_clusters)),
    source_household_clusters = length(unique(source_clusters)),
    weight_unchanged = TRUE, event_time_unchanged = TRUE,
    labor_outcome_unchanged = TRUE, stringsAsFactors = FALSE
  )
  attr(out, "housing_join_audit") <- join_audit
  out
}

fbh_validate_coding_config <- function(coding_config) {
  if (!is.list(coding_config)) fbh_stop("coding_config must be an explicit list")
  needed <- c("rooms_valid", "rooms_transform", "rooms_cap", "rooms_missing_codes",
              "rooms_unknown_codes", "bedrooms_valid", "bedrooms_transform", "bedrooms_cap",
              "bedrooms_missing_codes", "bedrooms_unknown_codes", "ownership_valid",
              "ownership_missing_codes", "ownership_unknown_codes", "allow_uncapped_sensitivity")
  if (!all(needed %in% names(coding_config)))
    fbh_stop("coding_config lacks explicit fields: ", paste(setdiff(needed, names(coding_config)), collapse = ", "))
  for (nm in c("rooms_valid", "rooms_transform", "bedrooms_valid", "bedrooms_transform", "ownership_valid"))
    if (!is.function(coding_config[[nm]])) fbh_stop("coding_config$", nm, " must be a function")
  if (!is.numeric(coding_config$rooms_cap) || length(coding_config$rooms_cap) != 1L ||
      !is.numeric(coding_config$bedrooms_cap) || length(coding_config$bedrooms_cap) != 1L)
    fbh_stop("housing caps must be scalar numeric values")
  if (!is.logical(coding_config$allow_uncapped_sensitivity) ||
      length(coding_config$allow_uncapped_sensitivity) != 1L ||
      is.na(coding_config$allow_uncapped_sensitivity))
    fbh_stop("coding_config$allow_uncapped_sensitivity must be one non-missing logical")
  invisible(TRUE)
}

fbh_apply_code <- function(raw, year, valid_fn, transform_fn, cap, missing_codes, unknown_codes,
                           prefix, allow_uncapped = FALSE) {
  valid <- suppressWarnings(as.logical(valid_fn(raw, year)))
  if (length(valid) != length(raw) || any(is.na(valid) & !is.na(raw)))
    fbh_stop(prefix, " validity function returned malformed values")
  valid[is.na(valid)] <- FALSE
  nonmissing <- !is.na(raw)
  invalid <- nonmissing & !valid
  missing_code <- nonmissing & raw %in% missing_codes
  unknown_code <- nonmissing & raw %in% unknown_codes
  transformed <- suppressWarnings(as.numeric(transform_fn(raw, year)))
  if (length(transformed) != length(raw)) fbh_stop(prefix, " transform returned wrong length")
  if (any(valid & !is.finite(transformed))) fbh_stop(prefix, " valid code transformed to non-finite value")
  clean <- rep(NA_real_, length(raw)); clean[valid] <- transformed[valid]
  primary <- clean
  primary[valid] <- pmin(clean[valid], cap)
  out <- list(clean = clean, primary = primary, invalid = invalid,
              missing_code = missing_code, unknown_code = unknown_code)
  if (allow_uncapped) out$uncapped <- clean
  out
}

code_first_birth_housing <- function(joined_panel, coding_config) {
  fbh_validate_coding_config(coding_config)
  fbh_require_columns(joined_panel, c("raw_rooms", "raw_bedrooms", "raw_ownership", "YEAR"), "joined_panel")
  yr <- suppressWarnings(as.numeric(joined_panel$YEAR))
  r <- fbh_apply_code(joined_panel$raw_rooms, yr, coding_config$rooms_valid,
                      coding_config$rooms_transform, coding_config$rooms_cap,
                      coding_config$rooms_missing_codes, coding_config$rooms_unknown_codes,
                      "ROOMS", coding_config$allow_uncapped_sensitivity)
  b <- fbh_apply_code(joined_panel$raw_bedrooms, yr, coding_config$bedrooms_valid,
                      coding_config$bedrooms_transform, coding_config$bedrooms_cap,
                      coding_config$bedrooms_missing_codes, coding_config$bedrooms_unknown_codes,
                      "BEDROOMS", coding_config$allow_uncapped_sensitivity)
  ov <- suppressWarnings(as.logical(coding_config$ownership_valid(joined_panel$raw_ownership, yr)))
  if (length(ov) != nrow(joined_panel) || any(is.na(ov) & !is.na(joined_panel$raw_ownership)))
    fbh_stop("OWNERSHP validity function returned malformed values")
  ov[is.na(ov)] <- FALSE
  own_nonmissing <- !is.na(joined_panel$raw_ownership)
  own_invalid <- own_nonmissing & !ov
  own_missing <- own_nonmissing & joined_panel$raw_ownership %in% coding_config$ownership_missing_codes
  own_unknown <- own_nonmissing & joined_panel$raw_ownership %in% coding_config$ownership_unknown_codes
  own <- rep(NA_real_, nrow(joined_panel))
  own[ov] <- ifelse(joined_panel$raw_ownership[ov] == 1, 1, 0)
  out <- joined_panel
  out$rooms_valid <- r$clean; out$rooms9 <- r$primary
  out$bedrooms_valid <- b$clean; out$bedrooms5 <- b$primary
  out$ownership_lw <- own
  out$rooms_invalid_code <- r$invalid; out$rooms_missing_code <- r$missing_code; out$rooms_unknown_code <- r$unknown_code
  out$bedrooms_invalid_code <- b$invalid; out$bedrooms_missing_code <- b$missing_code; out$bedrooms_unknown_code <- b$unknown_code
  out$ownership_invalid_code <- own_invalid; out$ownership_missing_code <- own_missing; out$ownership_unknown_code <- own_unknown
  if (coding_config$allow_uncapped_sensitivity) {
    out$rooms_uncapped <- r$uncapped; out$bedrooms_uncapped <- b$uncapped
  }
  code_audit <- data.frame(
    outcome = c("rooms", "bedrooms", "ownership"),
    raw_nonmissing = c(sum(!is.na(out$raw_rooms)), sum(!is.na(out$raw_bedrooms)), sum(!is.na(out$raw_ownership))),
    valid = c(sum(!is.na(out$rooms_valid)), sum(!is.na(out$bedrooms_valid)), sum(!is.na(out$ownership_lw))),
    invalid = c(sum(out$rooms_invalid_code), sum(out$bedrooms_invalid_code), sum(out$ownership_invalid_code)),
    missing_code = c(sum(out$rooms_missing_code), sum(out$bedrooms_missing_code), sum(out$ownership_missing_code)),
    unknown_code = c(sum(out$rooms_unknown_code), sum(out$bedrooms_unknown_code), sum(out$ownership_unknown_code)),
    primary_nonmissing = c(sum(!is.na(out$rooms9)), sum(!is.na(out$bedrooms5)), sum(!is.na(out$ownership_lw))),
    stringsAsFactors = FALSE
  )
  attr(out, "housing_code_audit") <- code_audit
  out
}

fbh_support <- function(d, outcomes, event_times, cohort_col) {
  base_cols <- c("statename", "gender", "event_time", cohort_col, "housing_join_status")
  support_keys <- d[, base_cols, drop = FALSE]
  for (nm in base_cols) {
    support_keys[[nm]] <- as.character(support_keys[[nm]])
    support_keys[[nm]][is.na(support_keys[[nm]])] <- "<missing>"
  }
  group_factor <- do.call(interaction, c(unname(support_keys),
                                         list(drop = TRUE, lex.order = TRUE)))
  groups <- split(seq_len(nrow(d)), group_factor, drop = TRUE)
  rows <- vector("list", length(groups) * length(outcomes)); k <- 0L
  for (ix in groups) for (o in outcomes) {
    k <- k + 1L
    one <- support_keys[ix[1], , drop = FALSE]
    cl <- d$source_household_cluster[ix]
    rows[[k]] <- data.frame(outcome = o, statename = as.character(one$statename),
      gender = as.character(one$gender), event_time = as.character(one$event_time),
      cohort = as.character(one[[cohort_col]]), housing_join_status = as.character(one$housing_join_status),
      rows = length(ix), valid_weight_rows = sum(is.finite(d$wgt[ix]) & d$wgt[ix] > 0),
      weight_sum = sum(d$wgt[ix][is.finite(d$wgt[ix]) & d$wgt[ix] > 0], na.rm = TRUE),
      missing_outcome = sum(is.na(d[[o]][ix])), source_household_clusters = length(unique(cl[!is.na(cl)])),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows[seq_len(k)])
}

fbh_term <- function(coef_names, state, event_time) {
  candidates <- c(paste0("statename", state, ":t_es", event_time),
                  paste0("t_es", event_time, ":statename", state))
  hit <- which(coef_names %in% candidates)
  if (length(hit) != 1L) fbh_stop("missing or ambiguous state/event coefficient for ", state, "/", event_time)
  coef_names[hit]
}

fbh_curve_and_contrast <- function(reg, Vh, data, outcome, event_times, ref,
                                   pre_times, post_times, state, gender, variance_label) {
  b <- stats::coef(reg); V <- as.matrix(stats::vcov(reg)); bh <- b
  rows <- vector("list", length(event_times))
  for (j in seq_along(event_times)) {
    tt <- event_times[j]
    if (tt == ref) {
      est <- 0; se <- 0; term <- NA_character_; se_h <- 0
    } else {
      term <- fbh_term(names(b), state, tt)
      if (!is.finite(b[[term]])) fbh_stop("non-finite coefficient for ", outcome, "/", state, "/", gender, "/", tt)
      est <- unname(b[[term]]); se <- sqrt(max(0, unname(V[term, term])))
      term_h <- fbh_term(names(bh), state, tt)
      se_h <- sqrt(max(0, unname(Vh[term_h, term_h])))
    }
    rows[[j]] <- data.frame(outcome = outcome, statename = state, gender = gender,
      event_time = tt, estimate = est, std.error = se,
      conf.low = est - qnorm(.975) * se, conf.high = est + qnorm(.975) * se,
      term = term, reference = tt == ref, variance = variance_label,
      nobs = stats::nobs(reg), source_household_clusters = length(unique(data$source_household_cluster)),
      stringsAsFactors = FALSE)
  }
  curve <- do.call(rbind, rows)
  contrast <- function(kind, times) {
    if (!all(times %in% event_times)) fbh_stop(kind, " contrast requests unsupported event time")
    get_result <- function(vec, mat) {
      L <- setNames(numeric(length(vec)), names(vec)); est <- 0
      for (tt in times) {
        sgn <- if (kind == "post_minus_pre") if (tt %in% post_times) 1 / length(post_times) else -1 / length(pre_times) else if (tt == "3") 1 else -1
        if (tt == ref) next
        tm <- fbh_term(names(vec), state, tt); L[[tm]] <- L[[tm]] + sgn; est <- est + sgn * vec[[tm]]
      }
      list(estimate = unname(est), se = sqrt(max(0, as.numeric(t(L) %*% mat %*% L))))
    }
    z <- get_result(b, V); zh <- get_result(bh, Vh)
    data.frame(outcome = outcome, statename = state, gender = gender, contrast = kind,
      estimate = z$estimate, std.error = z$se, conf.low = z$estimate - qnorm(.975) * z$se,
      conf.high = z$estimate + qnorm(.975) * z$se, std.error_hetero = zh$se,
      variance = variance_label, stringsAsFactors = FALSE)
  }
  list(curve = curve, summary = rbind(contrast("post_minus_pre", c(pre_times, post_times)),
                                      contrast("event_3_minus_event_neg1", c("3", "-1"))))
}

estimate_first_birth_housing <- function(panel, source_housing, audit_manifest, coding_config,
                                         event_times = as.character(c(-5:-1, 0:10)), ref = "-2",
                                         pre_times = as.character(-5:-1), post_times = as.character(0:10),
                                         division_value = "New England", state_max = 57L,
                                         outcomes = c("rooms9", "bedrooms5", "ownership_lw"),
                                         cohort_col = "cohort", labor_outcome = "emp_lw", weight_col = "wgt",
                                         checkpoint = NULL) {
  if (!(ref %in% event_times)) fbh_stop("reference event -2 is absent from event_times")
  fbh_require_columns(panel, c("census", "statefip", "statename", "gender", "age_factor",
                               "doiy_factor", cohort_col), "panel")
  if (!all(c(ref, pre_times, post_times, "3", "-1") %in% event_times))
    fbh_stop("event_times do not cover the requested reference/pre/post/+3/-1 support")
  joined <- join_first_birth_housing(panel, source_housing, audit_manifest,
                                     labor_outcome = labor_outcome, weight_col = weight_col)
  d <- code_first_birth_housing(joined, coding_config)
  d$event_time <- as.character(d$t_es_lw)
  d$in_event_window <- d$event_time %in% event_times
  d$analysis_geography <- as.character(d$census) == division_value &
    suppressWarnings(as.numeric(d$statefip)) < state_max
  if (!any(d$analysis_geography & d$in_event_window)) fbh_stop("no supported NE event rows")
  d$source_household_cluster <- ifelse(d$analysis_geography & !is.na(d$housing_source_key),
    paste(d$source_origin, d$YEAR, d$SAMPLE, d$SERIAL, sep = ":"), NA_character_)
  d$wgt <- d[[weight_col]]
  if (!is.null(checkpoint) && !is.function(checkpoint)) fbh_stop("checkpoint must be a function")
  support <- fbh_support(d[d$analysis_geography, , drop = FALSE], outcomes, event_times, cohort_col)
  analysis <- d[d$analysis_geography & d$in_event_window & is.finite(d$wgt) & d$wgt > 0, , drop = FALSE]
  if (!nrow(analysis)) fbh_stop("no positive-weight rows in analysis window")
  if (!requireNamespace("fixest", quietly = TRUE)) fbh_stop("fixest is required for estimation")
  analysis$statename <- droplevels(factor(analysis$statename))
  analysis$gender <- droplevels(factor(analysis$gender))
  analysis$t_es <- factor(analysis$event_time, levels = c(ref, setdiff(event_times, ref)))
  curves <- list(); summaries <- list(); fits <- list(); kk <- 0L
  for (outcome in outcomes) {
    if (!(outcome %in% names(analysis))) fbh_stop("unknown housing outcome ", outcome)
    od <- analysis[!is.na(analysis[[outcome]]) & !is.na(analysis$source_household_cluster), , drop = FALSE]
    if (!nrow(od)) fbh_stop("no observed housing outcome rows for ", outcome)
    if (length(unique(od$statename)) < 2L) fbh_stop("fewer than two supported states for ", outcome)
    for (gender in levels(od$gender)) {
      gd <- droplevels(od[od$gender == gender, , drop = FALSE])
      if (!nrow(gd)) fbh_stop("empty gender cell for ", outcome, "/", gender)
      counts <- table(gd$statename, gd$t_es)
      if (any(counts[, event_times, drop = FALSE] == 0))
        fbh_stop("missing state/event support for ", outcome, "/", gender)
      fml <- stats::as.formula(paste0(outcome, " ~ statename + t_es:statename | age_factor + doiy_factor"))
      reg <- fixest::feols(fml, data = gd, weights = stats::as.formula(paste0("~", weight_col)),
                           vcov = ~source_household_cluster)
      Vh <- as.matrix(stats::vcov(reg, vcov = "hetero"))
      fit_states <- levels(gd$statename)
      for (state in fit_states) {
        z <- fbh_curve_and_contrast(reg, Vh, gd, outcome, event_times, ref,
                                    pre_times, post_times, state, gender,
                                    "source_household_cluster")
        kk <- kk + 1L; curves[[kk]] <- z$curve; summaries[[kk]] <- z$summary
      }
      fit_name <- paste(outcome, gender, sep = "::")
      fits[[fit_name]] <- list(cluster = reg, hetero_vcov = Vh,
        nobs = stats::nobs(reg), source_household_clusters = length(unique(gd$source_household_cluster)))
      if (!is.null(checkpoint)) checkpoint(list(name = fit_name, outcome = outcome,
        gender = gender, curve = do.call(rbind, lapply(seq_along(fit_states), function(ii) {
          # The just-computed state curves are the last length(fit_states) entries.
          curves[[kk - length(fit_states) + ii]]
        })), summary = do.call(rbind, lapply(seq_along(fit_states), function(ii) {
          summaries[[kk - length(fit_states) + ii]]
        })), nobs = stats::nobs(reg)))
    }
  }
  curves <- do.call(rbind, curves); summaries <- do.call(rbind, summaries)
  list(status = "ESTIMATION_COMPLETE_DIAGNOSTIC", data = d, curves = curves,
       state_tables = curves, summary = summaries, support = support,
       fits = fits, join_audit = attr(joined, "housing_join_audit"),
       code_audit = attr(d, "housing_code_audit"),
       metadata = list(formula = "outcome ~ statename + t_es:statename | age_factor + doiy_factor",
         reference_event = ref, event_times = event_times, primary_variance = "source_household_cluster",
         sensitivity_variance = "heteroskedastic author sensitivity", level_coefficients = TRUE,
         outcomes = outcomes, causal_interpretation = "diagnostic matched pseudo-panel contrast"))
}
