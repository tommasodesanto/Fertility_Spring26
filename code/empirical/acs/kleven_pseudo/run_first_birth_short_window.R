#!/usr/bin/env Rscript

# Prepared short-window sensitivity driver.  It reuses the verified source-key
# join and the approved estimator; it never rematches observations or submits a
# cluster job.  Normal execution is intentionally allocation-only because the
# persisted v5 panel is large.  Set SHORT_WINDOW_TEST=1 for the fixture test.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

SHORT_EVENTS <- c(-2L, -1L, 0L, 1L, 2L, 3L)
FIXED_COHORTS <- 2007:2016
SOURCE_YEARS <- 2005:2019
NE_STATES <- c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont")

positive_weight <- function(x) is.finite(x) & x > 0

short_window_support <- function(d, outcomes = c("rooms9", "bedrooms5", "ownership_lw"),
                                 events = SHORT_EVENTS, cohorts = FIXED_COHORTS,
                                 states = NULL, genders = NULL) {
  needed <- c("statename", "gender", "cohort", "event_time", "housing_join_status", "wgt")
  miss <- setdiff(needed, names(d))
  if (length(miss)) stop("support input missing: ", paste(miss, collapse = ", "))
  d <- as.data.table(d)
  d[, event_num := suppressWarnings(as.integer(as.character(event_time)))]
  d[, cohort_num := suppressWarnings(as.integer(as.character(cohort)))]
  d <- d[event_num %in% events & cohort_num %in% cohorts &
           housing_join_status == "acs_source_matched" & positive_weight(wgt)]
  if (!nrow(d)) stop("short-window support has no source-matched positive-weight rows")
  cells <- rbindlist(lapply(outcomes, function(o) {
    if (!o %in% names(d)) stop("support input missing outcome: ", o)
    z <- d[, .(source_rows = .N,
               outcome_valid_rows = sum(!is.na(get(o))),
               missing_outcome_rows = sum(is.na(get(o))),
               weight_sum = sum(wgt),
               sum_weight_sq = sum(wgt^2)),
           by = .(statename, gender, cohort = cohort_num, event_time = event_num)]
    z[, outcome := o]
    z
  }), fill = TRUE)
  if (is.null(states)) states <- sort(unique(as.character(d$statename)))
  if (is.null(genders)) genders <- sort(unique(as.character(d$gender)))
  keys <- CJ(outcome = outcomes, statename = states,
             gender = genders, cohort = cohorts,
             event_time = events, unique = TRUE)
  cells <- merge(keys, cells,
                 by = c("outcome", "statename", "gender", "cohort", "event_time"),
                 all.x = TRUE)
  for (j in c("source_rows", "outcome_valid_rows", "missing_outcome_rows",
              "weight_sum", "sum_weight_sq"))
    set(cells, which(is.na(cells[[j]])), j, 0)
  cells[, cell_support := outcome_valid_rows > 0]
  cells[]
}

short_window_gate <- function(cells, events = SHORT_EVENTS) {
  cells <- as.data.table(cells)
  cells[, .(event_cells = .N,
            supported_event_cells = sum(cell_support),
            all_six_supported = .N == length(events) && all(cell_support),
            source_rows = sum(source_rows),
            outcome_valid_rows = sum(outcome_valid_rows),
            missing_outcome_rows = sum(missing_outcome_rows),
            weight_sum = sum(weight_sum)),
        by = .(outcome, statename, gender, cohort)]
}

weighted_baseline_ess <- function(d, outcomes = c("rooms9", "bedrooms5", "ownership_lw"),
                                  baseline_event = -2L, cohorts = FIXED_COHORTS) {
  d <- as.data.table(d)
  d[, event_num := suppressWarnings(as.integer(as.character(event_time)))]
  d[, cohort_num := suppressWarnings(as.integer(as.character(cohort)))]
  d <- d[event_num == baseline_event & cohort_num %in% cohorts &
           housing_join_status == "acs_source_matched" & positive_weight(wgt)]
  if (!nrow(d)) stop("no source-matched baseline rows for weighted baseline/ESS")
  rbindlist(lapply(outcomes, function(o) {
    z <- d[!is.na(get(o)), .(
      n_rows = .N,
      weighted_baseline = sum(get(o) * wgt) / sum(wgt),
      weight_sum = sum(wgt),
      sum_weight_sq = sum(wgt^2),
      kish_ess = sum(wgt)^2 / sum(wgt^2),
      source_household_clusters = uniqueN(source_household_cluster)
    ), by = .(statename, gender, cohort = cohort_num)]
    z[, outcome := o]
    z
  }), fill = TRUE)
}

short_window_prepare <- function(d, states, genders,
                                 outcomes = c("rooms9", "bedrooms5", "ownership_lw"),
                                 events = SHORT_EVENTS, cohorts = FIXED_COHORTS) {
  d <- as.data.table(d)
  if (all(c("analysis_geography", "in_event_window") %in% names(d)))
    d <- d[analysis_geography & in_event_window]
  cells <- short_window_support(d, outcomes, events, cohorts, states, genders)
  gate <- short_window_gate(cells, events)
  baseline <- weighted_baseline_ess(d, outcomes, -2L, cohorts)
  list(data = d, cells = cells, gate = gate, baseline = baseline)
}

fixture <- function() {
  stopifnot(identical(SOURCE_YEARS, 2005:2019), identical(FIXED_COHORTS, 2007:2016))
  d <- CJ(statename = c("A", "B"), gender = "Men", cohort = 2007:2008,
          event_time = SHORT_EVENTS, unique = TRUE)
  d[, `:=`(housing_join_status = "acs_source_matched", wgt = 1: .N,
           rooms9 = 1, bedrooms5 = 2, ownership_lw = 1,
           source_household_cluster = paste0("h", seq_len(.N)))]
  # One unsupported outcome cell must fail its cohort gate; the other cells pass.
  d[statename == "B" & cohort == 2008 & event_time == 3, rooms9 := NA_real_]
  d[, `:=`(analysis_geography = TRUE, in_event_window = TRUE)]
  prepared <- short_window_prepare(d, states = c("A", "B"), genders = "Men",
                                   cohorts = 2007:2008)
  cells <- prepared$cells
  gate <- prepared$gate
  stopifnot(nrow(cells) == 3L * 2L * 1L * 2L * 6L)
  stopifnot(gate[outcome == "rooms9" & statename == "A" & cohort == 2007, all_six_supported])
  stopifnot(!gate[outcome == "rooms9" & statename == "B" & cohort == 2008, all_six_supported])
  b <- weighted_baseline_ess(d, cohorts = 2007:2008)
  dd <- d[statename == "A" & cohort == 2007 & event_time == -2]
  expected <- sum(dd$wgt)^2 / sum(dd$wgt^2)
  stopifnot(abs(b[outcome == "rooms9" & statename == "A" & cohort == 2007, kish_ess] - expected) < 1e-12)
  cat("short-window fixture PASS: six-cell gate and Kish ESS\n")
}

if (identical(Sys.getenv("SHORT_WINDOW_TEST", "0"), "1")) {
  fixture()
  quit(save = "no", status = 0, runLast = FALSE)
}

root <- Sys.getenv("KLEVEN_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
outdir <- Sys.getenv("OUTDIR", file.path(root, "output/kleven_acs_pilot/first_birth_short_window_20260920"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
progress_file <- file.path(outdir, "progress.log")
checkpoint <- function(stage, detail = "") {
  z <- sprintf("%s\t%s\t%s", format(Sys.time(), "%FT%T%z"), stage, detail)
  write(z, progress_file, append = TRUE)
  cat(z, "\n")
}
write_json(list(status = "STARTED_NO_SUBMIT", source_years = SOURCE_YEARS,
                implied_cohorts = FIXED_COHORTS, event_times = SHORT_EVENTS,
                reference_event = -2L, contrast = "+3 minus -1 using full covariance",
                matching = "reuse deterministic verified source-key join; no rematching",
                generated = as.character(Sys.time())),
           file.path(outdir, "short_window_start_receipt.json"), auto_unbox = TRUE, pretty = TRUE)

# The approved base driver performs the exact dependency, provenance, source
# digest, overlap, and real-key smoke checks.  RUN_ESTIMATION=0 leaves its
# loaded objects available in this private environment without fitting.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("cannot resolve short-window driver path")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
base_driver <- file.path(script_dir, "run_first_birth_housing.R")
ctx <- new.env(parent = globalenv())
old <- Sys.getenv(c("OUTDIR", "RUN_ESTIMATION", "TINY_DRIVER_SMOKE", "ESTIMATOR_FILE"), unset = NA_character_)
on.exit({
  if (is.na(old[["OUTDIR"]])) Sys.unsetenv("OUTDIR") else Sys.setenv(OUTDIR = old[["OUTDIR"]])
  if (is.na(old[["RUN_ESTIMATION"]])) Sys.unsetenv("RUN_ESTIMATION") else Sys.setenv(RUN_ESTIMATION = old[["RUN_ESTIMATION"]])
  if (is.na(old[["TINY_DRIVER_SMOKE"]])) Sys.unsetenv("TINY_DRIVER_SMOKE") else Sys.setenv(TINY_DRIVER_SMOKE = old[["TINY_DRIVER_SMOKE"]])
  if (is.na(old[["ESTIMATOR_FILE"]])) Sys.unsetenv("ESTIMATOR_FILE") else Sys.setenv(ESTIMATOR_FILE = old[["ESTIMATOR_FILE"]])
}, add = TRUE)
Sys.setenv(OUTDIR = outdir, RUN_ESTIMATION = "0", TINY_DRIVER_SMOKE = "0",
           ESTIMATOR_FILE = file.path(script_dir, "estimate_first_birth_housing.R"))
checkpoint("base_driver_start")
sys.source(base_driver, envir = ctx)
checkpoint("base_driver_ready")

panel <- ctx$panel
source_housing <- ctx$source_housing
manifest <- ctx$manifest
coding <- ctx$coding
estimator <- ctx$estimate_first_birth_housing
if (!is.function(estimator)) stop("base driver did not expose estimator")
source_year <- suppressWarnings(as.numeric(panel$YEAR))
cohort <- suppressWarnings(as.numeric(panel$cohort))
true_acs <- as.character(panel$source_origin) == "ACS" & as.integer(panel$from_cps) == 0L
keep <- true_acs & source_year %in% SOURCE_YEARS & cohort %in% FIXED_COHORTS
if (!any(keep)) stop("source-year/cohort restriction produced no panel rows")
panel_short <- panel[keep, , drop = FALSE]
checkpoint("panel_restricted", sprintf("rows=%s", nrow(panel_short)))

# Join and code once before fitting so support failure is a durable stop, not a
# post-fit discovery.  The estimator repeats this deterministic key lookup for
# its fit object; it never performs a new match.
joined_prefit <- ctx$join_first_birth_housing(panel_short, source_housing, manifest)
coded_prefit <- ctx$code_first_birth_housing(joined_prefit, coding)
coded_prefit$event_time <- as.character(coded_prefit$t_es_lw)
coded_prefit$in_event_window <- coded_prefit$event_time %in% as.character(SHORT_EVENTS)
coded_prefit$analysis_geography <- as.character(coded_prefit$census) == "New England" &
  suppressWarnings(as.numeric(coded_prefit$statefip)) < 57
coded_prefit$source_household_cluster <- ifelse(
  coded_prefit$analysis_geography & !is.na(coded_prefit$housing_source_key),
  paste(coded_prefit$source_origin, coded_prefit$YEAR, coded_prefit$SAMPLE,
        coded_prefit$SERIAL, sep = ":"), NA_character_)
coded_prefit$wgt <- coded_prefit$wgt
prefit <- short_window_prepare(coded_prefit, states = NE_STATES,
                               genders = c("Men", "Women"))
fwrite(prefit$cells, file.path(outdir, "short_window_support_cells_prefit.csv"))
fwrite(prefit$gate, file.path(outdir, "short_window_support_gate_prefit.csv"))
fwrite(prefit$baseline, file.path(outdir, "short_window_baseline_ess_prefit.csv"))
fwrite(attr(coded_prefit, "housing_code_audit"), file.path(outdir, "short_window_code_audit_prefit.csv"))
if (any(!prefit$gate$all_six_supported)) {
  bad <- prefit$gate[!all_six_supported]
  write_json(list(status = "SUPPORT_INCOMPLETE", source_years = SOURCE_YEARS,
                  implied_cohorts = FIXED_COHORTS, bad_groups = bad,
                  message = "at least one requested outcome/state/gender/cohort lacks positive outcome-valid support in one of six event cells; fitting stopped before estimation"),
             file.path(outdir, "short_window_failure_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
  stop("SUPPORT_INCOMPLETE: pre-fit six-cell support gate failed")
}
checkpoint("prefit_support_pass", sprintf("gate=%s/%s", sum(prefit$gate$all_six_supported), nrow(prefit$gate)))

result <- estimator(panel_short, source_housing, manifest, coding,
                    event_times = as.character(SHORT_EVENTS), ref = "-2",
                    pre_times = "-1", post_times = as.character(0:3),
                    cohort_col = "cohort",
                    checkpoint = function(x) {
                      nm <- gsub("[^A-Za-z0-9_.-]", "_", x$name)
                      saveRDS(x, file.path(outdir, paste0("fit_", nm, "_checkpoint.rds")))
                      checkpoint("fit_complete", x$name)
                    })
checkpoint("estimation_complete", result$status)

result_view <- result$data[result$data$analysis_geography & result$data$in_event_window, , drop = FALSE]
postfit <- short_window_prepare(result_view, states = NE_STATES,
                                genders = c("Men", "Women"))
cells <- postfit$cells
gate <- postfit$gate
baseline <- postfit$baseline
fwrite(cells, file.path(outdir, "short_window_support_cells.csv"))
fwrite(gate, file.path(outdir, "short_window_support_gate.csv"))
fwrite(baseline, file.path(outdir, "short_window_baseline_ess.csv"))
fwrite(result$curves, file.path(outdir, "short_window_curves.csv"))
fwrite(result$summary, file.path(outdir, "short_window_contrasts.csv"))
fwrite(result$join_audit, file.path(outdir, "short_window_join_audit.csv"))
write_json(list(
  status = "ESTIMATION_COMPLETE_DIAGNOSTIC",
  source_years = SOURCE_YEARS,
  implied_cohorts = FIXED_COHORTS,
  event_times = SHORT_EVENTS,
  reference_event = -2L,
  contrast = "+3 minus -1 using full covariance",
  matching = "reuse deterministic verified source-key join; no rematching",
  support_gate = "every outcome/state/gender/cohort cell has positive outcome-valid rows in all six events",
  gate_groups = nrow(gate), gate_passing = sum(gate$all_six_supported),
  weighted_baseline = "event -2, source-matched positive-weight rows, outcome-valid only",
  kish_ess = "(sum(w)^2)/sum(w^2), computed from joined data",
  joined_rows = nrow(result$data), join_audit = result$join_audit,
  generated = as.character(Sys.time())
), file.path(outdir, "short_window_result_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
checkpoint("receipts_written", sprintf("gate=%s/%s", sum(gate$all_six_supported), nrow(gate)))
