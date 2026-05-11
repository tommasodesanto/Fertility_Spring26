#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) == 0) {
  stop("Could not determine script path from commandArgs().")
}
script_dir <- dirname(normalizePath(script_path))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(repo_root, "output")
bin_mode <- Sys.getenv("ROOM_BIN_MODE", unset = "coarse")

acs_summary_path <- file.path(out_dir, "acs_2023_rooms_summary.csv")
acs_bins_path <- file.path(out_dir, "acs_2023_rooms_bins.csv")
model_summary_path <- file.path(out_dir, "model_rooms_summary.csv")
model_bins_path <- file.path(out_dir, "model_rooms_bins.csv")

for (path in c(acs_summary_path, acs_bins_path, model_summary_path, model_bins_path)) {
  if (!file.exists(path)) {
    stop(sprintf("Missing input file: %s", path))
  }
}

acs_summary <- fread(acs_summary_path)
acs_bins <- fread(acs_bins_path)
model_summary <- fread(model_summary_path)
model_bins <- fread(model_bins_path)

key_cols <- c("age_window", "tenure", "child_bin")
cmp_summary <- merge(
  acs_summary,
  model_summary,
  by = key_cols,
  suffixes = c("_acs", "_model")
)
cmp_summary[, `:=`(
  mean_gap = mean_model - mean_acs,
  p50_gap = p50_model - p50_acs,
  share_ge_8_gap = share_ge_8_model - share_ge_8_acs,
  share_ge_11_gap = share_ge_11_model - share_ge_11_acs
)]

cmp_bins <- merge(
  acs_bins,
  model_bins,
  by = c(key_cols, "bin"),
  suffixes = c("_acs", "_model")
)
cmp_bins[, share_gap := share_model - share_acs]

setorder(cmp_summary, age_window, tenure, child_bin)
setorder(cmp_bins, age_window, tenure, child_bin, bin)

fwrite(cmp_summary, file.path(out_dir, "room_distribution_compare_summary.csv"))
fwrite(cmp_bins, file.path(out_dir, "room_distribution_compare_bins.csv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_dt <- cmp_bins[age_window == "25_45" & child_bin %in% c("0", "1", "2+")]
  plot_dt <- rbindlist(list(
    plot_dt[, .(tenure, child_bin, bin, source = "ACS", share = share_acs)],
    plot_dt[, .(tenure, child_bin, bin, source = "Model", share = share_model)]
  ))
  bin_levels <- if (identical(bin_mode, "split_5_6")) {
    c("<=4", "5", "6", "7-8", "9-10", "11+")
  } else {
    c("<=4", "5-6", "7-8", "9-10", "11+")
  }
  plot_dt[, bin := factor(bin, levels = bin_levels)]
  plot_dt[, child_bin := factor(child_bin, levels = c("0", "1", "2+"))]
  plot_dt[, tenure := factor(tenure, levels = c("Renter", "Owner"))]
  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = bin, y = share, fill = source)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_grid(tenure ~ child_bin) +
    ggplot2::scale_fill_manual(values = c("ACS" = "#577590", "Model" = "#c8553d")) +
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%d%%", round(100 * x))) +
    ggplot2::labs(
      title = "Prime-Age Room Distribution: ACS vs Live Benchmark",
      subtitle = "Householder age 25-45, current children in household",
      x = "Rooms bin",
      y = "Share",
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  ggplot2::ggsave(
    filename = file.path(out_dir, "room_distribution_compare_25_45.png"),
    plot = p,
    width = 12,
    height = 7,
    dpi = 180
  )
}

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt_num <- function(x) sprintf("%.3f", x)

report_path <- file.path(out_dir, "room_distribution_compare_report.txt")
fid <- file(report_path, open = "wt")
on.exit(close(fid), add = TRUE)

cat("ACS vs model room-distribution comparison\n", file = fid)
cat(sprintf("date = %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = fid)

cat("What is and is not targeted in the live calibration\n", file = fid)
cat("- The live objective does not include any room-distribution moments.\n", file = fid)
cat("- These comparisons therefore audit the implied housing-size distribution rather than a targeted fit object.\n\n", file = fid)

cat("Top-line all-household tenure medians\n", file = fid)
topline <- cmp_summary[age_window == "all" & child_bin == "all", .(
  tenure,
  acs_p50 = p50_acs,
  model_p50 = p50_model,
  acs_ge_8 = share_ge_8_acs,
  model_ge_8 = share_ge_8_model,
  acs_ge_11 = share_ge_11_acs,
  model_ge_11 = share_ge_11_model
)]
cat(paste(capture.output(print(topline)), collapse = "\n"), file = fid)
cat("\n\n", file = fid)

cat("Prime-age (25-45) room summaries by tenure x current-child bin\n", file = fid)
prime <- cmp_summary[age_window == "25_45" & child_bin %in% c("0", "1", "2+"), .(
  tenure,
  child_bin,
  acs_p50 = p50_acs,
  model_p50 = p50_model,
  acs_ge_8 = share_ge_8_acs,
  model_ge_8 = share_ge_8_model,
  acs_ge_11 = share_ge_11_acs,
  model_ge_11 = share_ge_11_model
)]
cat(paste(capture.output(print(prime)), collapse = "\n"), file = fid)
cat("\n\n", file = fid)

cat("Prime-age (25-45) room-bin shares: ACS vs model\n", file = fid)
prime_bins <- cmp_bins[age_window == "25_45" & child_bin %in% c("0", "1", "2+"), .(
  tenure,
  child_bin,
  bin,
  share_acs,
  share_model,
  share_gap
)]
cat(paste(capture.output(print(prime_bins)), collapse = "\n"), file = fid)
cat("\n\n", file = fid)

cat("Headline interpretation\n", file = fid)

owner0 <- cmp_summary[age_window == "25_45" & tenure == "Owner" & child_bin == "0"]
owner1 <- cmp_summary[age_window == "25_45" & tenure == "Owner" & child_bin == "1"]
owner2 <- cmp_summary[age_window == "25_45" & tenure == "Owner" & child_bin == "2+"]
renter0 <- cmp_summary[age_window == "25_45" & tenure == "Renter" & child_bin == "0"]
renter1 <- cmp_summary[age_window == "25_45" & tenure == "Renter" & child_bin == "1"]
renter2 <- cmp_summary[age_window == "25_45" & tenure == "Renter" & child_bin == "2+"]
renter1_bin_78 <- cmp_bins[age_window == "25_45" & tenure == "Renter" & child_bin == "1" & bin == "7-8", share_acs]
renter1_model_bin_78 <- cmp_bins[age_window == "25_45" & tenure == "Renter" & child_bin == "1" & bin == "7-8", share_model]

cat(sprintf(
  "- Prime-age medians are much closer than in the older benchmark, but the model is too compressed. Owner medians are ACS %s/%s/%s for child bins 0/1/2+, versus model %s/%s/%s.\n",
  fmt_num(owner0$p50_acs), fmt_num(owner1$p50_acs), fmt_num(owner2$p50_acs),
  fmt_num(owner0$p50_model), fmt_num(owner1$p50_model), fmt_num(owner2$p50_model)
), file = fid)

cat(sprintf(
  "- Renter medians show the same compression. ACS is %s/%s/%s for child bins 0/1/2+, versus model %s/%s/%s.\n",
  fmt_num(renter0$p50_acs), fmt_num(renter1$p50_acs), fmt_num(renter2$p50_acs),
  fmt_num(renter0$p50_model), fmt_num(renter1$p50_model), fmt_num(renter2$p50_model)
), file = fid)

cat(sprintf(
  "- The distributional miss is shape, not just levels. For one-child renters, ACS has %s in 7-8 rooms and %s in 11+; the model has %s in 7-8 and %s in 11+.\n",
  fmt_pct(renter1_bin_78),
  fmt_pct(renter1$share_ge_11_acs),
  fmt_pct(renter1_model_bin_78),
  fmt_pct(renter1$share_ge_11_model)
), file = fid)

cat(sprintf(
  "- Owners are also missing the upper family-size tail. For one-child owners age 25-45, ACS has %s in 11+ rooms; the model has %s. For two-plus-child owners, ACS has %s in 11+; the model has %s.\n",
  fmt_pct(owner1$share_ge_11_acs), fmt_pct(owner1$share_ge_11_model),
  fmt_pct(owner2$share_ge_11_acs), fmt_pct(owner2$share_ge_11_model)
), file = fid)

cat("- On this evidence, the current starter-gap anchor is more credible on room levels than the older benchmark, but it still lacks the middle-to-upper family housing gradient visible in the ACS.\n", file = fid)

cat("\nroom_distribution_compare_report.txt\n")
