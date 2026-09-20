# Mothers-only pooled-NE housing diagnostic for the coresident second-birth proxy.
# This file prepares two explicitly different estimands and fits only the
# descriptive event-factor model. It does not claim causal identification.

if (!requireNamespace("data.table", quietly = TRUE)) stop("second_birth_housing.R requires data.table")
.sb2_require <- function(d, cols, label) { miss <- setdiff(cols, names(d)); if (length(miss)) stop(sprintf("%s missing columns: %s", label, paste(miss, collapse = ", "))) }
.sb2_source_cluster <- function(d) { .sb2_require(d, c("YEAR", "SAMPLE", "SERIAL"), "source rows"); paste(d$YEAR, d$SAMPLE, d$SERIAL, sep = "\034") }
.sb2_attach_source_fields <- function(d, source, fields) {
  .sb2_require(d, "person_key", "proxy rows"); .sb2_require(source, c("person_key", fields), "source rows")
  if (anyDuplicated(source$person_key)) stop("source person_key must be unique")
  idx <- match(d$person_key, source$person_key)
  if (anyNA(idx)) stop("used proxy row is missing from the matched housing source")
  for (v in fields) if (!v %in% names(d)) d[[v]] <- source[[v]][idx]
  used_cluster <- d[, c("YEAR", "SAMPLE", "SERIAL"), with = FALSE]
  if (any(vapply(used_cluster, function(x) any(is.na(x) | !nzchar(as.character(x))), logical(1))))
    stop("used housing row has missing source household cluster component")
  d
}
.sb2_joint_negative_ids <- function(target_support, neg_times = c(-2L, -1L)) {
  .sb2_require(target_support, c("anchor_person_key", "target_event_time", "has_donor"), "target support")
  s <- target_support[target_event_time %in% neg_times, .(joint_negative_supported = length(unique(target_event_time)) == length(neg_times) && all(has_donor)), by = anchor_person_key]
  s[joint_negative_supported == TRUE, anchor_person_key]
}
.sb2_reaggregate_links <- function(matched, keep_ids = NULL, neg_times = c(-2L, -1L)) {
  links <- data.table::copy(matched$links)
  .sb2_require(links, c("anchor_person_key", "target_event_time", "target_year", "target_mother_age", "target_child_age", "donor_person_key", "donor_household_key", "donor_PERWT", "wgt_match"), "matching links")
  links <- links[target_event_time %in% neg_times]; if (!is.null(keep_ids)) links <- links[anchor_person_key %in% keep_ids]
  if (!nrow(links)) return(data.table::data.table(donor_person_key=character(), donor_household_key=character(), target_event_time=integer(), target_year=integer(), target_mother_age=numeric(), target_child_age=numeric(), PERWT=numeric(), wgt_match=numeric(), n_target_rows=integer(), person_key=character(), wgt_original=numeric(), wgt=numeric()))
  out <- links[, .(PERWT=donor_PERWT[1L], wgt_match=sum(wgt_match), n_target_rows=.N), by=.(donor_person_key, donor_household_key, target_event_time, target_year, target_mother_age, target_child_age)]
  out[, `:=`(person_key=donor_person_key, wgt_original=as.numeric(PERWT), wgt=as.numeric(PERWT)*as.numeric(wgt_match))]
  out
}

prepare_second_birth_housing <- function(proxy, matched, source_rows = NULL,
    outcomes = c("rooms9", "bedrooms5", "ownership_lw"), post_times = 0:3,
    neg_times = c(-2L, -1L), support_spec = c("all_anchors", "joint_negative"),
    fertyr_spec = c("all", "event0_yes")) {
  support_spec <- match.arg(support_spec); fertyr_spec <- match.arg(fertyr_spec)
  if (is.null(source_rows)) source_rows <- proxy$input
  .sb2_require(proxy, c("anchors", "post_rows"), "proxy checkpoint"); .sb2_require(matched, c("links", "target_support"), "matching checkpoint")
  source_rows <- data.table::as.data.table(data.table::copy(source_rows))
  .sb2_require(source_rows, c("person_key", "YEAR", "SAMPLE", "SERIAL", "PERWT", "AGE_norm", "STATEFIP", outcomes), "housing source")
  if (anyDuplicated(source_rows$person_key)) stop("housing source person_key must be unique")
  anchors <- data.table::as.data.table(data.table::copy(proxy$anchors)); post <- data.table::as.data.table(data.table::copy(proxy$post_rows))
  .sb2_require(anchors, c("person_key", "birth_gap"), "proxy anchors"); .sb2_require(post, c("person_key", "event_time", "birth_gap", "PERWT"), "proxy post rows")
  anchor_ids <- anchors[birth_gap >= 2 & is.finite(birth_gap), person_key]
  if (fertyr_spec == "event0_yes") { if (!"FERTYR_status" %in% names(anchors)) stop("event0_yes requires FERTYR_status"); anchor_ids <- anchors[birth_gap >= 2 & FERTYR_status == "yes", person_key] }
  joint_ids <- .sb2_joint_negative_ids(matched$target_support, neg_times)
  selected_ids <- if (support_spec == "joint_negative") intersect(anchor_ids, joint_ids) else anchor_ids
  post <- post[birth_gap >= 2 & event_time %in% post_times]
  if (fertyr_spec == "event0_yes") { .sb2_require(post, "FERTYR_status", "proxy post rows"); post <- post[event_time != 0 | FERTYR_status == "yes"] }
  used_fields <- c("YEAR","SAMPLE","SERIAL","PERWT","AGE_norm","STATEFIP",outcomes)
  if ("SEX" %in% names(source_rows)) used_fields <- c(used_fields, "SEX")
  post <- data.table::copy(.sb2_attach_source_fields(post, source_rows, used_fields))
  if (any(!is.finite(as.numeric(post$PERWT)) | as.numeric(post$PERWT) <= 0)) stop("used post row has missing or nonpositive PERWT")
  if ("SEX" %in% names(source_rows) && any(is.na(post$SEX) | as.numeric(post$SEX) != 2)) stop("used post row is not a mother (SEX must equal 2)")
  post[, `:=`(event_time=as.integer(event_time), weight=as.numeric(PERWT), pseudo_role="post", source_household_cluster=.sb2_source_cluster(.SD)), .SDcols=c("YEAR","SAMPLE","SERIAL")]
  post[, anchor_population := if (fertyr_spec == "event0_yes")
    "event0_yes_anchor_validation" else "all_gap_ge_2"]
  negative_ids <- if (support_spec == "joint_negative") selected_ids else anchor_ids
  donor_w <- .sb2_reaggregate_links(matched, keep_ids=negative_ids, neg_times=neg_times)
  donor_fields <- c("YEAR","SAMPLE","SERIAL","AGE_norm","STATEFIP",outcomes)
  if ("SEX" %in% names(source_rows)) donor_fields <- c(donor_fields, "SEX")
  donor <- data.table::copy(.sb2_attach_source_fields(donor_w, source_rows, donor_fields))
  if (any(!is.finite(as.numeric(donor$wgt)) | as.numeric(donor$wgt) <= 0)) stop("used donor row has missing or nonpositive weight")
  if ("SEX" %in% names(source_rows) && any(is.na(donor$SEX) | as.numeric(donor$SEX) != 2)) stop("used donor row is not a mother (SEX must equal 2)")
  donor[, `:=`(event_time=as.integer(target_event_time), weight=as.numeric(wgt), pseudo_role="negative_donor", source_household_cluster=.sb2_source_cluster(.SD), anchor_population=if (support_spec == "joint_negative") "joint_negative_supported" else if (fertyr_spec == "event0_yes") "event0_yes_anchor_validation" else "all_gap_ge_2"), .SDcols=c("YEAR","SAMPLE","SERIAL")]
  if (anyDuplicated(post[, .(person_key,event_time)])) stop("duplicate post person/event rows")
  panel <- data.table::rbindlist(list(post[, c("person_key","YEAR","SAMPLE","SERIAL","AGE_norm","STATEFIP","event_time","weight","pseudo_role","source_household_cluster","anchor_population",outcomes), with=FALSE], donor[, c("person_key","YEAR","SAMPLE","SERIAL","AGE_norm","STATEFIP","event_time","weight","pseudo_role","source_household_cluster","anchor_population",outcomes), with=FALSE]), fill=TRUE)
  panel[, event_factor := factor(event_time, levels=sort(unique(c(neg_times,post_times))))]
  list(panel=panel, post=post, negative_donor=donor, support=list(all_anchor_ids=anchor_ids, joint_negative_ids=joint_ids, selected_ids=selected_ids, n_all_anchors=length(anchor_ids), n_joint_negative=length(intersect(anchor_ids,joint_ids)), negative_times=neg_times, post_times=post_times, support_spec=support_spec, fertyr_spec=fertyr_spec), contract=list(population="mothers-only pooled NE strict gap>=2", post_weight="original PERWT once", donor_weight="donor PERWT times reaggregated Matching tie fraction", cluster="YEAR:SAMPLE:SERIAL", causal_claim=FALSE))
}

second_birth_housing_support <- function(prepared, outcomes=c("rooms9","bedrooms5","ownership_lw")) {
  d <- data.table::copy(prepared$panel)
  data.table::rbindlist(lapply(outcomes, function(y) {
    d[, .sb2_observed := is.finite(as.numeric(get(y))) & is.finite(weight) & weight > 0]
    d[, .(outcome = y, n_rows = .N, n_observed = sum(.sb2_observed),
          weighted_rows = sum(weight[.sb2_observed]),
          weight_ess = if (sum(.sb2_observed) > 0)
            sum(weight[.sb2_observed])^2 / sum(weight[.sb2_observed]^2) else NA_real_,
          source_household_clusters = data.table::uniqueN(source_household_cluster[.sb2_observed])),
      by = .(event_time, pseudo_role)]
  }), fill=TRUE)
}

.sb2_event_support <- function(d, y, expected_times = c(-2L, -1L, 0L, 1L, 2L, 3L)) {
  d <- data.table::copy(d)
  d[, .sb2_observed := is.finite(as.numeric(get(y))) & is.finite(weight) & weight > 0]
  out <- d[, .(n_rows = .N, n_observed = sum(.sb2_observed), weighted_rows = sum(weight[.sb2_observed]),
              weight_ess = if (sum(.sb2_observed)) sum(weight[.sb2_observed])^2 / sum(weight[.sb2_observed]^2) else NA_real_,
              source_household_clusters = data.table::uniqueN(source_household_cluster[.sb2_observed])),
          by = event_time]
  out <- merge(data.table::data.table(event_time = expected_times), out, by = "event_time", all.x = TRUE, sort = TRUE)
  out[is.na(n_rows), `:=`(n_rows = 0L, n_observed = 0L, weighted_rows = 0, source_household_clusters = 0L)]
  out
}

.sb2_event_terms <- function(fit, expected_times = c(-1L, 0L, 1L, 2L, 3L)) {
  b <- stats::coef(fit); v <- as.matrix(stats::vcov(fit))
  terms <- paste0("event_factor::", expected_times)
  if (any(!terms %in% names(b))) stop("estimation omitted an event coefficient; event cell is not identifiable")
  if (any(!is.finite(b[terms])) || any(!is.finite(v[terms, terms]))) stop("event coefficient or covariance is non-finite")
  setNames(terms, as.character(expected_times))
}

.sb2_curves <- function(fit, d, y, expected_times = c(-2L, -1L, 0L, 1L, 2L, 3L)) {
  s <- .sb2_event_support(d, y, expected_times)
  b <- stats::coef(fit); v <- as.matrix(stats::vcov(fit)); terms <- .sb2_event_terms(fit, expected_times[-1L])
  est <- setNames(rep(0, length(expected_times)), as.character(expected_times))
  se <- setNames(rep(0, length(expected_times)), as.character(expected_times))
  est[as.character(expected_times[-1L])] <- b[terms]
  se[as.character(expected_times[-1L])] <- sqrt(diag(v)[terms])
  crit <- stats::qnorm(.975)
  s[, `:=`(estimate = unname(est[as.character(event_time)]), std_error = unname(se[as.character(event_time)]),
           reference_event = event_time == -2L, outcome = y)]
  s[, `:=`(lower = estimate - crit * std_error, upper = estimate + crit * std_error)]
  raw <- d[is.finite(as.numeric(get(y))) & is.finite(weight) & weight > 0,
           .(raw_event_mean = sum(get(y) * weight) / sum(weight)), by = event_time]
  s <- merge(s, raw, by = "event_time", all.x = TRUE, sort = TRUE)
  baseline <- s[event_time == -2L, raw_event_mean][1L]
  s[, baseline_raw_mean := baseline]
  s[, pre_minus1_coefficient := ifelse(event_time == -1L, estimate, NA_real_)]
  s
}

fit_second_birth_housing <- function(prepared, outcomes=c("rooms9","bedrooms5","ownership_lw"), checkpoint=NULL, save_dir=NULL) {
  if (!requireNamespace("fixest", quietly=TRUE)) stop("fit requires fixest")
  d <- data.table::copy(prepared$panel); expected <- c(-2L,-1L,0L,1L,2L,3L)
  .sb2_require(d,c("event_time","AGE_norm","YEAR","STATEFIP","weight","source_household_cluster",outcomes),"estimation panel")
  if (!setequal(sort(unique(d$event_time)), expected)) stop("estimation panel must contain all six event-time cells")
  d[, event_factor := factor(event_time, levels=expected)]; if (!is.null(save_dir)) dir.create(save_dir,recursive=TRUE,showWarnings=FALSE)
  fits <- list(); curves <- list(); contrasts <- list()
  for (y in outcomes) {
    dd <- d[is.finite(as.numeric(get(y))) & is.finite(weight) & weight > 0]
    if (!nrow(dd)) stop(sprintf("outcome %s has no usable observations", y))
    event_counts <- dd[, .N, by = event_time]
    if (any(!expected %in% event_counts$event_time) || any(event_counts[match(expected, event_time), N] < 1L)) stop(sprintf("outcome %s is missing an event-time cell", y))
    fml <- stats::as.formula(sprintf("%s ~ i(event_factor, ref = '-2') | STATEFIP + AGE_norm + YEAR", y))
    nobs_prefit <- nrow(dd)
    fit <- fixest::feols(fml,data=dd,weights=~weight,cluster=~source_household_cluster)
    term_map <- .sb2_event_terms(fit, expected[-1L])
    cv <- as.matrix(stats::vcov(fit)); b <- stats::coef(fit)
    contrast_vec <- c(1, -1); contrast_terms <- term_map[c("3", "-1")]
    contrast_est <- sum(contrast_vec * b[contrast_terms]); contrast_var <- as.numeric(t(contrast_vec) %*% cv[contrast_terms, contrast_terms] %*% contrast_vec)
    if (!is.finite(contrast_var) || contrast_var < -1e-10) stop("contrast variance is materially negative or non-finite")
    contrast_var <- max(0, contrast_var); contrast_se <- sqrt(contrast_var)
    curve <- .sb2_curves(fit, dd, y, expected)
    curve[, `:=`(contrast_plus3_minus_minus1 = contrast_est, contrast_se = contrast_se,
                 contrast_lower = contrast_est - stats::qnorm(.975) * contrast_se,
                 contrast_upper = contrast_est + stats::qnorm(.975) * contrast_se,
                 nobs_prefit = nobs_prefit, nobs_fit = stats::nobs(fit))]
    rec <- list(outcome=y, formula=fml, fit=fit, coefficients=b[term_map], vcov=cv[term_map,term_map,drop=FALSE], nobs_prefit=nobs_prefit, nobs=stats::nobs(fit), curves=curve, contrast=list(estimate=contrast_est,variance=contrast_var,se=contrast_se,lower=contrast_est-stats::qnorm(.975)*contrast_se,upper=contrast_est+stats::qnorm(.975)*contrast_se,vcov=cv[contrast_terms,contrast_terms,drop=FALSE]), support=second_birth_housing_support(list(panel=dd),outcomes=y))
    fits[[y]] <- rec; curves[[y]] <- curve; contrasts[[y]] <- rec$contrast
    if (!is.null(save_dir)) saveRDS(rec,file.path(save_dir,paste0("fit_",y,"_checkpoint.rds")),compress=FALSE)
    if (!is.null(checkpoint)) checkpoint(rec)
  }
  list(fits=fits, curves=data.table::rbindlist(curves, fill=TRUE), contrasts=contrasts, support=second_birth_housing_support(prepared,outcomes), specification="pooled NE mothers-only; event factor ref -2; categorical STATEFIP + AGE_norm + YEAR fixed effects; source-household clustered")
}
