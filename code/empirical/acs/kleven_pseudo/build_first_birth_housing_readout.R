#!/usr/bin/env Rscript

# Build compact, durable readouts from a completed first-birth housing run.
# This script only reads saved fit/checkpoint and CSV receipts.  It does not
# load the ACS panel, rematch observations, or submit a cluster job.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(MASS)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[[1]] else Sys.getenv(
  "FIRST_BIRTH_OUTDIR",
  unset = file.path("code", "empirical", "acs", "kleven_pseudo",
                    "output", "first_birth_housing_20260920_job18079576")
)
if (!dir.exists(outdir)) stop("saved output directory is missing: ", outdir)
required <- c("curves.csv", "contrasts.csv", "support.csv", "join_audit.csv",
              "housing_code_audit.csv", "result_receipt.json", "readiness_receipt.json")
missing <- required[!file.exists(file.path(outdir, required))]
if (length(missing)) stop("saved readout inputs missing: ", paste(missing, collapse = ", "))

curves <- fread(file.path(outdir, "curves.csv"))
contrasts <- fread(file.path(outdir, "contrasts.csv"))
support <- fread(file.path(outdir, "support.csv"))
join_audit <- fread(file.path(outdir, "join_audit.csv"))
code_audit <- fread(file.path(outdir, "housing_code_audit.csv"))
receipt <- fromJSON(file.path(outdir, "result_receipt.json"), simplifyVector = TRUE)
readiness <- fromJSON(file.path(outdir, "readiness_receipt.json"), simplifyVector = TRUE)
audit_manifest <- readiness$audit_manifest
source_sha256 <- as.character(audit_manifest$source_sha256)
source_bytes <- as.numeric(audit_manifest$source_bytes)
if (!grepl("^[0-9a-fA-F]{64}$", source_sha256)) stop("readiness receipt has no valid 64-character source SHA-256")
if (!is.finite(source_bytes)) stop("readiness receipt has no finite source byte count")

num_event <- function(x) suppressWarnings(as.numeric(as.character(x)))
support[, event_num := num_event(event_time)]
matched <- support[housing_join_status == "acs_source_matched"]

# Source-matched support by outcome, state, and gender.  The compact support
# receipt does not retain sum(w^2), so weighted ESS is deliberately left NA.
state_gender_support <- matched[
  , .(source_matched_rows = sum(rows),
      outcome_valid_rows = sum(rows - missing_outcome),
      missing_outcome_rows = sum(missing_outcome),
      valid_weight_rows = sum(valid_weight_rows),
      weight_sum = sum(weight_sum)),
  by = .(outcome, statename, gender)
]
state_gender_support[, `:=`(
  outcome_valid_share = outcome_valid_rows / source_matched_rows,
  weighted_ess = NA_real_,
  weighted_ess_definition = "(sum(w)^2)/sum(w^2); sum(w^2) not retained in compact receipt"
)]
fwrite(state_gender_support, file.path(outdir, "state_gender_support.csv"))

# Short-window support is a gate over the six requested event cells.  It does
# not require zero missing outcomes: it requires positive observed outcome
# support in every source-matched cell and reports missingness separately.
short_events <- c(-2, -1, 0, 1, 2, 3)
short_cell <- matched[event_num %in% short_events,
  .(source_matched_rows = sum(rows),
    outcome_valid_rows = sum(rows - missing_outcome),
    missing_outcome_rows = sum(missing_outcome),
    valid_weight_rows = sum(valid_weight_rows),
    weight_sum = sum(weight_sum)),
  by = .(outcome, statename, gender, event_time = event_num)
]
keys <- CJ(outcome = unique(curves$outcome), statename = unique(curves$statename),
           gender = unique(curves$gender), event_time = short_events,
           unique = TRUE)
short_cell <- merge(keys, short_cell,
                    by = c("outcome", "statename", "gender", "event_time"),
                    all.x = TRUE)
for (j in c("source_matched_rows", "outcome_valid_rows", "missing_outcome_rows",
            "valid_weight_rows", "weight_sum")) set(short_cell, which(is.na(short_cell[[j]])), j, 0)
short_cell[, cell_support := outcome_valid_rows > 0]
short_gate <- short_cell[, .(
  event_cells = .N,
  supported_event_cells = sum(cell_support),
  all_six_supported = .N == length(short_events) && all(cell_support),
  source_matched_rows = sum(source_matched_rows),
  outcome_valid_rows = sum(outcome_valid_rows),
  missing_outcome_rows = sum(missing_outcome_rows),
  valid_weight_rows = sum(valid_weight_rows),
  weight_sum = sum(weight_sum)),
  by = .(outcome, statename, gender)
]
fwrite(short_cell, file.path(outdir, "short_window_support_cells.csv"))
fwrite(short_gate, file.path(outdir, "short_window_support_gate.csv"))

recipe <- list(
  status = "PREPARED_NOT_RUN",
  purpose = "same-estimator common implied-event-cohort short-window sensitivity",
  source_year_overlap = 2005:2019,
  source_year_filter = "future estimator invocation must restrict true ACS source YEAR to 2005:2019; the compact support receipt does not retain source YEAR by cell",
  event_times = short_events,
  reference_event = -2,
  pre_event_times = -1,
  post_event_times = 0:3,
  contrast = "+3 minus -1 using the full covariance matrix",
  matching = "reuse completed source-key matches; no rematching",
  cohort = "support-only implied_pseudo_event_year = doiy - numeric(t_es_lw), with endpoint bins retained only as labels",
  support_gate = "source-matched true ACS rows with outcome_valid_rows > 0 in every event cell; after applying source YEAR 2005:2019 in the future invocation; missingness reported, never imputed",
  variance = "source-household clustered primary, heteroskedastic sensitivity",
  runtime_note = "No allocation submitted. The full saved fit took 7:52 including serialization; a short-window fit is expected to be shorter but must be timed at launch.",
  causal_status = "diagnostic sensitivity conditional on constructed matches; no causal claim"
)
write_json(recipe, file.path(outdir, "short_window_sensitivity_recipe.json"),
           auto_unbox = TRUE, pretty = TRUE)
source_recipe <- Sys.getenv("FIRST_BIRTH_SOURCE_RECIPE", unset = "")
if (nzchar(source_recipe))
  write_json(recipe, source_recipe, auto_unbox = TRUE, pretty = TRUE)

# Joint pre-event Wald checks use the saved clustered covariance.  The reference
# event -2 is excluded by construction.  A generalized inverse is used only to
# report the estimable rank, with the rank/df recorded explicitly.
checkpoint_dir <- file.path(outdir, "checkpoints")
wald_rows <- list()
if (dir.exists(checkpoint_dir)) {
  cp <- list.files(checkpoint_dir, pattern = "_checkpoint\\.rds$", full.names = TRUE)
  for (f in cp) {
    fit <- readRDS(f)
    b <- fit$coefficients
    V <- fit$vcov$cluster
    states <- unique(as.character(fit$curve$statename))
    for (st in states) {
      ev <- c(-5, -4, -3, -1)
      terms <- paste0("statename", st, ":t_es", ev)
      idx <- match(terms, names(b))
      keep <- !is.na(idx)
      idx <- idx[keep]
      if (!length(idx)) next
      bb <- as.numeric(b[idx])
      VV <- V[idx, idx, drop = FALSE]
      ee <- eigen((VV + t(VV)) / 2, symmetric = TRUE, only.values = TRUE)$values
      tol <- max(dim(VV)) * max(ee, 0) * 1e-10
      rank <- sum(ee > tol)
      if (rank > 0) {
        stat <- as.numeric(t(bb) %*% ginv(VV) %*% bb)
        pval <- pchisq(stat, df = rank, lower.tail = FALSE)
      } else {
        stat <- NA_real_; pval <- NA_real_
      }
      cdat <- fit$curve[as.character(fit$curve$statename) == st &
                          num_event(fit$curve$event_time) %in% ev, ]
      wald_rows[[length(wald_rows) + 1L]] <- data.table(
        outcome = fit$outcome, gender = fit$gender, statename = st,
        pre_event_times = paste(ev[keep], collapse = ","),
        wald_chisq = stat, wald_rank = rank, wald_df = rank,
        wald_pvalue = pval,
        max_abs_pre_estimate = max(abs(cdat$estimate), na.rm = TRUE),
        mean_abs_pre_estimate = mean(abs(cdat$estimate), na.rm = TRUE),
        interpretation = "pattern diagnostic, not a validity proof"
      )
    }
  }
}
pretrend <- if (length(wald_rows)) rbindlist(wald_rows, fill = TRUE) else data.table()
fwrite(pretrend, file.path(outdir, "pretrend_wald_clustered.csv"))

# Save an explicit clustered-versus-heteroskedastic comparison for the two
# primary saved contrasts.  Estimates are unchanged; only uncertainty changes.
se_compare <- contrasts[contrast %in% c("post_minus_pre", "event_3_minus_event_neg1"),
  .(outcome, statename, gender, contrast, estimate,
    se_clustered = std.error, se_heteroskedastic = std.error_hetero,
    cluster_to_hetero_se = std.error / std.error_hetero,
    variance_primary = variance)]
fwrite(se_compare, file.path(outdir, "contrast_se_comparison.csv"))

# Re-render the saved curves with readable outcome labels.  State panels remain
# separate so lines cannot be mistaken for one another.
label_map <- c(rooms9 = "rooms", bedrooms5 = "bedrooms", ownership_lw = "ownership (pp)")
states <- sort(unique(as.character(curves$statename)))
gcols <- c(Men = "#1b6ca8", Women = "#c23b22")
png(file.path(outdir, "housing_event_curves_labeled.png"),
    width = max(1800, 500 * length(states)), height = 1300, res = 150)
par(mfrow = c(3, length(states)), mar = c(3.2, 3.2, 2.2, 0.8), oma = c(0, 0, 0, 0))
for (o in c("rooms9", "bedrooms5", "ownership_lw")) {
  zo <- curves[outcome == o]
  ylim <- range(c(zo$conf.low, zo$conf.high), finite = TRUE)
  ylim <- ylim + c(-1, 1) * max(diff(ylim) * 0.03, 0.01)
  for (st in states) {
    z <- zo[as.character(statename) == st]
    xx <- sort(unique(num_event(z$event_time)))
    plot(range(xx), ylim, type = "n", xlab = "Event time", ylab = label_map[[o]], main = st)
    for (g in names(gcols)) {
      q <- z[as.character(gender) == g][order(num_event(event_time))]
      if (nrow(q)) {
        polygon(c(num_event(q$event_time), rev(num_event(q$event_time))),
                c(q$conf.low, rev(q$conf.high)),
                col = adjustcolor(gcols[[g]], alpha.f = 0.15), border = NA)
        lines(num_event(q$event_time), q$estimate, type = "b", pch = 16, col = gcols[[g]])
      }
    }
    abline(v = -2, lty = 2)
    if (o == "rooms9" && st == states[[1]])
      legend("topleft", legend = names(gcols), col = gcols, lty = 1, pch = 16, bty = "n")
  }
}
dev.off()

# Write a concise durable narrative with the numbers that reviewers need.
ja <- join_audit[1]
room <- code_audit[outcome == "rooms"][1]
bed <- code_audit[outcome == "bedrooms"][1]
own <- code_audit[outcome == "ownership"][1]
support_total <- matched[, .(rows = sum(rows), valid = sum(rows - missing_outcome), missing = sum(missing_outcome))]
gate_n <- short_gate[, .(groups = .N, passing = sum(all_six_supported))]
v6_path <- "/scratch/td2248/projects/kleven_acs_pilot_20260917/overnight_ne_benchmark/ne_housing_v6_18047612/ne_housing_receipt.json"
candidate_dir <- file.path(dirname(outdir), "second_birth_candidate_support_20260920")
candidate_lines <- character()
if (file.exists(file.path(candidate_dir, "candidate_support_overall.csv"))) {
  cand <- fread(file.path(candidate_dir, "candidate_support_overall.csv"))
  cv <- setNames(as.character(cand$value), cand$metric)
  candidate_lines <- c(
    "",
    "## Second-birth candidate availability diagnostic",
    "",
    paste0("- The same allocation's bounded candidate-support step found ", cv[["target_rows"]], " target rows, ", cv[["target_rows_with_any_donor"]], " target rows with any donor, and ", cv[["target_rows_with_any_eligible_donor"]], " with any eligible donor under the declared exact demographic cells."),
    paste0("- It found ", cv[["full_pre_target_rows"]], " fixed full-pre-window target rows and ", cv[["full_pre_target_rows_with_any_donor"]], " with any donor; the event-time receipt is `second_birth_candidate_support_20260920/candidate_support_by_event_time.csv`. This is candidate availability only, with no matching assignment, coarsening, or causal interpretation."),
    ""
  )
}
report <- c(
  "# First-birth housing run readout",
  "",
  "This report is generated from the saved receipts and per-fit checkpoints for Torch job 18079576. It is a diagnostic matched pseudo-panel readout; the estimates are conditional on the constructed source-key matches and are not causal claims.",
  "",
  "## Run and failure evidence",
  "",
  paste0("- Job 18079496: failed after 33 seconds at the driver source-key assertion: `true ACS source key is not SAMPLE:YEAR:SERIAL:PERNUM`. The candidate-support step had already written its compact outputs. The failure receipt is `failure_receipt_18079496.json`.",
         " The source packet contained observed whitespace-delimited tokens (`SAMPLE YEAR SERIAL PERNUM`), while relabeled CPS keys retained the `CPS:` form."),
  "- Job 18079576: completed in 7:52 after the parser accepted the observed whitespace source-key form and retained the source-key guard. No new allocation is active.",
  "- The corrected run passed dependency smoke, real matched-key smoke, source overlap validation, and the estimator interface before fitting.",
  "",
  "## Saved-source and join evidence",
  "",
  paste0("- Verified source packet: ", format(readiness$source_packet_rows, big.mark = ","), " rows; SHA-256 `", source_sha256, "`; ", formatC(source_bytes, format = "f", digits = 0, big.mark = ","), " bytes.",
         " Verified unique key intersection is ", format(audit_manifest$unique_key_intersection, big.mark = ","), " and shared source years are 2005--2019."),
  paste0("- V5 panel rows: ", format(ja$panel_rows, big.mark = ","), "; true ACS rows: ", format(ja$acs_rows, big.mark = ","), "; original CPS rows: ", format(ja$original_cps_missing_outcome_rows, big.mark = ","), "; relabeled CPS rows: ", format(ja$relabeled_cps_missing_outcome_rows, big.mark = ","), "; total CPS rows in the join audit: ", format(ja$cps_rows, big.mark = ","), "."),
  paste0("- True ACS rows matched to the verified source packet: ", format(ja$matched_acs_rows, big.mark = ","), "; unmatched true ACS rows: ", format(ja$unmatched_acs_rows, big.mark = ","), "."),
  paste0("- Repeated source-household clusters: ", format(ja$repeated_source_household_clusters, big.mark = ","), "; distinct source-household clusters: ", format(ja$source_household_clusters, big.mark = ","), "."),
  "- Weights, event time, and labor outcomes are unchanged by the source join.",
  "",
  "## Housing coding and support",
  "",
  paste0("- ROOMS: ", format(room$valid, big.mark = ","), " valid; ", format(room$unknown_code, big.mark = ","), " unknown code 28; ", format(room$missing_code, big.mark = ","), " missing code. The primary outcome is capped at 9."),
  paste0("- BEDROOMS: ", format(bed$valid, big.mark = ","), " valid; codes 1--22 are transformed by x-1 and capped at 5; ", format(bed$missing_code, big.mark = ","), " missing code."),
  paste0("- OWNERSHP: ", format(own$valid, big.mark = ","), " valid; codes 1 and 2 are retained as 1/0; ", format(own$missing_code, big.mark = ","), " missing code."),
  paste0("- Across outcome-specific source-matched support cells (three outcomes are reported separately): ", format(support_total$rows, big.mark = ","), " rows, ", format(support_total$valid, big.mark = ","), " outcome-valid rows, and ", format(support_total$missing, big.mark = ","), " missing-outcome rows. State-by-gender detail is in `state_gender_support.csv`; short-window cell detail is in `short_window_support_cells.csv` and its compact gate is in `short_window_support_gate.csv` (", gate_n$passing, "/", gate_n$groups, " groups pass positive outcome support before the future source-year filter).") ,
  "- The compact support receipt does not retain sum(w^2), so weighted effective sample size is not fabricated. Its definition for a future retained-weight receipt is `(sum(w)^2)/sum(w^2)`.",
  "",
  "## Saved-fit inference diagnostics",
  "",
  "- Primary fit: level coefficients with source-household clustered covariance; heteroskedastic covariance is a sensitivity comparison. `contrast_se_comparison.csv` records both standard errors for post-minus-pre and +3-minus-(-1).",
  "- `pretrend_wald_clustered.csv` reports joint clustered-V Wald tests for event times -5, -4, -3, and -1, excluding the -2 reference event, with estimable rank and df. These are pretrend pattern diagnostics, not a validity proof.",
  "- The compact fit has one gender-specific regression across the six states; the curves retain the common fit-level nobs and source-household-cluster counts. The reference event is normalized to zero, so the saved curves do not identify raw housing level means at event -2.",
  "",
  "## V6 comparison",
  "",
  paste0("- The earlier V6 ownership continuation (receipt: `", v6_path, "`) used the same V5 panel lineage, a true-ACS observed-ownership sample with `n_ownership=962538`, percentage-point ownership levels, and heteroskedastic standard errors. It excluded 38,464 relabeled CPS rows and 22,375 unknown true-ACS ownership codes."),
  "- The new readout uses the verified extract27 source packet, retains the full source-key join audit, and uses source-household clustering as primary with heteroskedastic standard errors as sensitivity. The sample/source-year overlap and variance estimator differ from V6; this is a specification comparison, not a quantitative decomposition or evidence of a failed replication.",
  "",
  "## Prepared next sensitivity",
  "",
  "- `short_window_sensitivity_recipe.json` prepares, but does not run, the common implied-event-cohort window [-2,+3] over the verified 2005--2019 overlap, with -2 as reference, -1 as the pre-event, +3 minus -1 using the full covariance matrix, and no rematching.",
  "- The cohort is a support label `doiy - numeric(t_es_lw)`, not a biological birth year. The compact gate requires positive outcome-valid source-matched support in every event cell; the future invocation must reapply that gate after restricting true ACS source YEAR to 2005--2019 and reports missing outcomes separately.",
  candidate_lines,
  "",
  "## Files",
  "",
  "- `housing_event_curves_labeled.png` (saved-CSV figure with visible labels rooms, bedrooms, and ownership (pp), separated by state).",
  "- `run_first_birth_short_window.R` is the reviewed no-submit driver for fixed implied cohorts 2007--2016 and true ACS source years 2005--2019; `test_first_birth_short_window.R` checks the six-cell support gate and Kish ESS on a tiny fixture.",
  "- `state_gender_support.csv`, `short_window_support_cells.csv`, `short_window_support_gate.csv`, `pretrend_wald_clustered.csv`, and `contrast_se_comparison.csv` are compact machine-readable receipts.",
  "- This report and its generator are durable source artifacts; the large raw panel, full `fits.rds`, and cluster submission are not copied or recomputed locally."
)
writeLines(report, file.path(outdir, "first_birth_housing_run_report.md"))
source_report <- Sys.getenv("FIRST_BIRTH_SOURCE_REPORT", unset = "")
if (nzchar(source_report)) writeLines(report, source_report)
cat("WROTE", file.path(outdir, "first_birth_housing_run_report.md"), "\n")
cat("WROTE", file.path(outdir, "short_window_sensitivity_recipe.json"), "\n")
cat("WROTE", file.path(outdir, "housing_event_curves_labeled.png"), "\n")
