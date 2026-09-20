#!/usr/bin/env Rscript

# Compact readout from completed full- and short-window receipts.  No panel
# reload, rematching, estimation, or cluster submission occurs here.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
bundle <- if (length(args)) args[[1]] else file.path(
  "code", "empirical", "acs", "kleven_pseudo", "output",
  "first_birth_sensitivity_bundle_18080591")
full <- if (length(args) > 1) args[[2]] else file.path(
  "code", "empirical", "acs", "kleven_pseudo", "output",
  "first_birth_housing_20260920_job18079576")
dir.create(bundle, recursive = TRUE, showWarnings = FALSE)

second <- file.path(bundle, "second_birth")
short <- file.path(bundle, "short_window")
target <- fread(file.path(second, "target_support.csv"))
anchor <- fread(file.path(second, "anchor_support.csv"))
full_c <- fread(file.path(full, "contrasts.csv"))
short_c <- fread(file.path(short, "short_window_contrasts.csv"))
baseline <- fread(file.path(short, "short_window_baseline_ess.csv"))
curves <- fread(file.path(short, "short_window_curves.csv"))

# Joint -2/-1 donor support at the anchor level.
target[, event := as.integer(target_event_time)]
target <- target[target_match_excluded == FALSE]
joint <- target[, .(support_neg2 = any(event == -2L & has_donor),
                    support_neg1 = any(event == -1L & has_donor)), by = anchor_person_key]
joint_counts <- data.table(
  object = c("target_anchors", "target_neg2_supported", "target_neg1_supported",
             "target_joint_neg2_neg1_supported"),
  count = c(nrow(joint), sum(joint$support_neg2), sum(joint$support_neg1),
            sum(joint$support_neg2 & joint$support_neg1))
)
fwrite(joint_counts, file.path(bundle, "joint_neg2_neg1_support_counts.csv"))

# The anchor receipt contains source year in the person key and explicit gap
# indicators, but no state field.  Preserve the unavailable-state limitation.
anchor[, source_year := as.integer(tstrsplit(anchor_person_key, "/", fixed = TRUE)[[1]])]
year_support <- merge(
  anchor[, .(anchor_rows = .N,
            full_pre_supported = sum(full_pre_supported),
            reference_supported = sum(reference_supported),
            joint_fullpre_reference = sum(full_pre_supported & reference_supported)),
        by = source_year],
  target[, .(target_rows = uniqueN(anchor_person_key),
             target_neg2_supported = sum(event == -2L & has_donor),
             target_neg1_supported = sum(event == -1L & has_donor)),
        by = .(source_year = as.integer(tstrsplit(anchor_person_key, "/", fixed = TRUE)[[1]]))],
  by = "source_year", all = TRUE
)
fwrite(year_support, file.path(bundle, "joint_support_by_source_year.csv"))
gap_support <- anchor[, .(
  anchor_rows = .N,
  full_pre_supported = sum(full_pre_supported),
  reference_supported = sum(reference_supported),
  joint_fullpre_reference = sum(full_pre_supported & reference_supported)
), by = .(full_pre_gap = gap_full_pre, reference_gap = gap_reference)]
fwrite(gap_support, file.path(bundle, "joint_support_by_gap.csv"))

# Baseline/weight summary is computed from the joined short-window output.
baseline_summary <- baseline[, .(
  cells = .N,
  rows = sum(n_rows),
  weight_sum = sum(weight_sum),
  sum_weight_sq = sum(sum_weight_sq),
  weighted_baseline_min = min(weighted_baseline),
  weighted_baseline_max = max(weighted_baseline),
  kish_ess_min = min(kish_ess),
  kish_ess_max = max(kish_ess),
  source_household_clusters_min = min(source_household_clusters),
  source_household_clusters_max = max(source_household_clusters)
), by = outcome]
fwrite(baseline_summary, file.path(bundle, "short_baseline_weight_summary.csv"))

# Compare only the common +3 - (-1) contrast.  Post-minus-pre uses different
# windows in the two fits and is deliberately excluded from this table.
fc <- full_c[contrast == "event_3_minus_event_neg1",
             .(outcome, statename, gender,
               estimate_full = estimate, se_full_clustered = std.error,
               se_full_heteroskedastic = std.error_hetero)]
sc <- short_c[contrast == "event_3_minus_event_neg1",
              .(outcome, statename, gender,
                estimate_short = estimate, se_short_clustered = std.error,
                se_short_heteroskedastic = std.error_hetero)]
same_contrast <- merge(fc, sc, by = c("outcome", "statename", "gender"), all = TRUE)
same_contrast[, `:=`(
  estimate_short_minus_full = estimate_short - estimate_full,
  se_short_minus_full = se_short_clustered - se_full_clustered,
  contrast_definition = "+3 minus -1; full and short estimates use different event windows/cohort restrictions"
)]
fwrite(same_contrast, file.path(bundle, "full_vs_short_event3_minus_neg1.csv"))

# Render state-separated saved-CSV curves for the short window.  The saved
# ownership estimates remain proportions; convert them to percentage points
# only in this display layer.
label_map <- c(rooms9 = "Rooms (cap 9)", bedrooms5 = "Bedrooms (cap 5)", ownership_lw = "Ownership (pp)")
states <- sort(unique(as.character(curves$statename)))
cols <- c(Men = "#1b6ca8", Women = "#c23b22")
png(file.path(bundle, "short_window_housing_event_curves.png"),
    width = max(1800, 500 * length(states)), height = 1300, res = 150)
par(mfrow = c(3, length(states)), mar = c(3.4, 4.8, 2.6, 0.8))
for (o in c("rooms9", "bedrooms5", "ownership_lw")) {
  z0 <- copy(curves[outcome == o])
  if (o == "ownership_lw")
    z0[, `:=`(estimate = 100 * estimate,
              conf.low = 100 * conf.low,
              conf.high = 100 * conf.high)]
  ylim <- range(c(z0$conf.low, z0$conf.high), finite = TRUE)
  ylim <- ylim + c(-1, 1) * max(diff(ylim) * .03, .01)
  for (st in states) {
    z <- z0[as.character(statename) == st]
    xx <- as.numeric(as.character(z$event_time))
    plot(range(xx), ylim, type = "n", xlab = "Event time", ylab = label_map[[o]], main = st)
    for (g in names(cols)) {
      q <- z[as.character(gender) == g][order(as.numeric(as.character(event_time)))]
      polygon(c(as.numeric(q$event_time), rev(as.numeric(q$event_time))),
              c(q$conf.low, rev(q$conf.high)),
              col = adjustcolor(cols[[g]], alpha.f = .15), border = NA)
      lines(as.numeric(q$event_time), q$estimate, type = "b", pch = 16, col = cols[[g]])
    }
    abline(v = -2, lty = 2)
    if (o == "rooms9" && st == states[[1]])
      legend("topleft", legend = names(cols), col = cols, lty = 1, pch = 16, bty = "n")
  }
}
dev.off()

fmt <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
gate <- fread(file.path(short, "short_window_support_gate.csv"))
report <- c(
  "# Completed ACS sensitivity readout: job 18080591",
  "",
  "The corrected bundle completed both source-support and short-window housing stages. This readout uses saved receipts only and makes no causal claim.",
  "",
  paste0("- Joint target-anchor support at both event -2 and event -1: ", joint_counts[object == "target_joint_neg2_neg1_supported", count], " of ", joint_counts[object == "target_anchors", count], " target anchors; separate support counts are ", joint_counts[object == "target_neg2_supported", count], " at -2 and ", joint_counts[object == "target_neg1_supported", count], " at -1."),
  paste0("- Anchor receipt cross-check: ", sum(anchor$full_pre_supported & anchor$reference_supported), " anchors satisfy both full-pre and reference support; ", sum(anchor$full_pre_supported), " satisfy full-pre support and ", sum(anchor$reference_supported), " satisfy reference support."),
  "- Source-year support is reported in `joint_support_by_source_year.csv`; the compact anchor receipt has no state field, so a state breakdown is unavailable without reopening the joined source data.",
  paste0("- Gap flags in `joint_support_by_gap.csv`: reference_gap=TRUE for ", sum(anchor$gap_reference), " anchors and full_pre_gap=TRUE for ", sum(anchor$gap_full_pre), "; corresponding support flags are reference_supported=", sum(anchor$reference_supported), " and full_pre_supported=", sum(anchor$full_pre_supported), ". These are support flags, not estimates."),
  "- `full_vs_short_event3_minus_neg1.csv` retains state-by-gender cells; any cell summary across those rows is an unweighted descriptive average, not a pooled Northeast estimate.",
  "",
  "## Short-window housing support and weights",
  "",
  paste0("- The prefit gate passed ", sum(gate$all_six_supported), "/", nrow(gate), " requested outcome/state/gender/cohort groups across event times -2 through +3."),
  "- `short_baseline_weight_summary.csv` reports event -2 weighted housing baselines, source-household counts, weight sums, and Kish ESS computed from joined rows. The Kish definition is `(sum(w)^2)/sum(w^2)`.",
  "",
  "## Full versus short window",
  "",
  "- The only directly comparable contrast is +3 minus -1. The full-window fit uses the broader event window and its original cohort support; the short fit uses fixed implied cohorts 2007--2016 and true ACS source years 2005--2019. The post-minus-pre contrast is excluded because its post window differs.",
  "- `full_vs_short_event3_minus_neg1.csv` reports estimates and clustered/heteroskedastic standard errors side by side. Differences are descriptive specification comparisons conditional on constructed matches.",
  "- `short_window_housing_event_curves.png` renders the short-window saved curves with state-separated panels and readable outcome labels.",
  "",
  "## Durable next state",
  "",
  "- No new allocation was submitted after job 18080591. The next decision is second-birth estimator review; the short-window and source-support outputs are complete and reusable."
)
writeLines(report, file.path(bundle, "first_birth_sensitivity_readout_18080591.md"))
cat("WROTE", file.path(bundle, "first_birth_sensitivity_readout_18080591.md"), "\n")
