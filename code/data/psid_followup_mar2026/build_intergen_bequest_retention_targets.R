#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(haven)
})

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

weighted_quantile_safe <- function(x, w, prob = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= prob)[1L]]
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[startsWith(args, "--file=")]
  if (length(file_arg) == 0L) stop("Run this file with Rscript.")
  dirname(normalizePath(sub("--file=", "", file_arg[1L], fixed = TRUE)))
}

make_pairs <- function(dt, horizon, baseline_age_min, baseline_age_max,
                       sample_name, one_per_person = TRUE) {
  base <- copy(dt[
    age >= baseline_age_min & age <= baseline_age_max &
      is.finite(total_nw) & is.finite(nonhousing_nw) & is.finite(home_equity) &
      is.finite(income) & income > 1000 & is.finite(weight) & weight > 0 &
      is.finite(n_children)
  ])
  follow <- copy(dt[
    is.finite(total_nw) & is.finite(nonhousing_nw) & is.finite(home_equity) &
      is.finite(n_children)
  ])

  keep <- c(
    "id", "year", "age", "n_children", "income", "total_nw",
    "nonhousing_nw", "home_equity", "own", "weight", "married"
  )
  base <- base[, ..keep]
  follow <- follow[, ..keep]
  pairs <- merge(
    base, follow, by = "id", suffixes = c("0", "1"), allow.cartesian = TRUE
  )
  pairs <- pairs[
    year1 - year0 == horizon &
      age1 - age0 >= horizon - 1 & age1 - age0 <= horizon + 1
  ]
  pairs[, children_stable := n_children0 == n_children1]
  pairs <- pairs[children_stable == TRUE]

  if (one_per_person) {
    target_age <- (baseline_age_min + baseline_age_max) / 2
    pairs[, pair_score := abs(age0 - target_age)]
    setorder(pairs, id, pair_score, -year0)
    pairs <- pairs[, .SD[1L], by = id]
    pairs[, pair_score := NULL]
  }

  lag_candidates <- merge(
    pairs[, .(id, baseline_year = year0)],
    dt[, .(
      id, lag_year = year, lag_total_nw = total_nw,
      lag_income = income
    )],
    by = "id", allow.cartesian = TRUE
  )
  lag_candidates[, lag_gap := baseline_year - lag_year]
  lag_candidates <- lag_candidates[
    lag_gap >= 2 & lag_gap <= 6 & is.finite(lag_total_nw) &
      is.finite(lag_income) & lag_income > 1000
  ]
  lag_candidates[, lag_score := abs(lag_gap - 4)]
  setorder(lag_candidates, id, lag_score, -lag_year)
  lag <- lag_candidates[, .SD[1L], by = id][, .(
    id, lag_year, lag_gap, lag_total_nw, lag_income
  )]
  pairs <- merge(pairs, lag, by = "id", all.x = TRUE)

  pairs[, `:=`(
    horizon = horizon,
    sample = sample_name,
    n_children = n_children0,
    children_capped = pmin(n_children0, 3),
    fertility_bin = factor(fcase(
      n_children0 == 0, "0",
      n_children0 == 1, "1",
      n_children0 == 2, "2",
      n_children0 >= 3, "3+",
      default = NA_character_
    ), levels = c("0", "1", "2", "3+")),
    log_income0 = log(income0),
    total_initial_ihs = asinh(total_nw0 / income0)
  )]

  q1 <- weighted_quantile_safe(pairs$total_nw0, pairs$weight0, 1 / 3)
  q2 <- weighted_quantile_safe(pairs$total_nw0, pairs$weight0, 2 / 3)
  pairs[, initial_total_wealth_tertile := factor(fcase(
    total_nw0 <= q1, "bottom",
    total_nw0 <= q2, "middle",
    default = "top"
  ), levels = c("bottom", "middle", "top"))]
  pairs[, `:=`(wealth_tertile_cut1 = q1, wealth_tertile_cut2 = q2)]

  lag_ok <- is.finite(pairs$lag_total_nw)
  lag_q1 <- weighted_quantile_safe(
    pairs$lag_total_nw[lag_ok], pairs$weight0[lag_ok], 1 / 3
  )
  lag_q2 <- weighted_quantile_safe(
    pairs$lag_total_nw[lag_ok], pairs$weight0[lag_ok], 2 / 3
  )
  pairs[, lagged_total_wealth_tertile := factor(fcase(
    !is.finite(lag_total_nw), NA_character_,
    lag_total_nw <= lag_q1, "bottom",
    lag_total_nw <= lag_q2, "middle",
    default = "top"
  ), levels = c("bottom", "middle", "top"))]
  pairs[, `:=`(
    lagged_wealth_tertile_cut1 = lag_q1,
    lagged_wealth_tertile_cut2 = lag_q2
  )]

  for (asset in c("total", "nonhousing", "home_equity")) {
    w0 <- pairs[[paste0(asset, "_nw0")]]
    w1 <- pairs[[paste0(asset, "_nw1")]]
    if (asset == "home_equity") {
      w0 <- pairs$home_equity0
      w1 <- pairs$home_equity1
    }
    pairs[, (paste0(asset, "_initial_ihs")) := asinh(w0 / income0)]
    pairs[, (paste0(asset, "_ihs_change_annual")) :=
            (asinh(w1 / income0) - asinh(w0 / income0)) / horizon]
    pairs[, (paste0(asset, "_change_income_annual")) :=
            ((w1 - w0) / income0) / horizon]
    pairs[, (paste0(asset, "_retention_ratio")) :=
            fifelse(w0 > 10000, w1 / w0, NA_real_)]
  }
  pairs[]
}

summarize_group <- function(block, sample_name, asset, group_type, group_value) {
  w0 <- block[[paste0(asset, "0")]]
  w1 <- block[[paste0(asset, "1")]]
  if (asset == "total_nw") prefix <- "total"
  if (asset == "nonhousing_nw") prefix <- "nonhousing"
  if (asset == "home_equity") prefix <- "home_equity"
  ihs <- block[[paste0(prefix, "_ihs_change_annual")]]
  change_income <- block[[paste0(prefix, "_change_income_annual")]]
  retention <- block[[paste0(prefix, "_retention_ratio")]]
  data.table(
    sample = sample_name,
    asset = prefix,
    group_type = group_type,
    group_value = group_value,
    n = nrow(block),
    weight_sum = sum(block$weight0, na.rm = TRUE),
    baseline_wealth_mean = weighted_mean_safe(w0, block$weight0),
    baseline_wealth_median = weighted_quantile_safe(w0, block$weight0),
    followup_wealth_mean = weighted_mean_safe(w1, block$weight0),
    followup_wealth_median = weighted_quantile_safe(w1, block$weight0),
    ihs_change_annual_mean = weighted_mean_safe(ihs, block$weight0),
    ihs_change_annual_median = weighted_quantile_safe(ihs, block$weight0),
    change_income_annual_mean = weighted_mean_safe(change_income, block$weight0),
    change_income_annual_median = weighted_quantile_safe(change_income, block$weight0),
    retention_ratio_eligible_n = sum(is.finite(retention)),
    retention_ratio_median = weighted_quantile_safe(retention, block$weight0)
  )
}

summarize_pairs <- function(pairs) {
  rows <- list()
  assets <- c("total_nw", "nonhousing_nw", "home_equity")
  for (asset in assets) {
    rows[[length(rows) + 1L]] <- summarize_group(
      pairs, pairs$sample[1L], asset, "all", "all"
    )
    for (g in levels(pairs$fertility_bin)) {
      rows[[length(rows) + 1L]] <- summarize_group(
        pairs[fertility_bin == g], pairs$sample[1L], asset,
        "completed_fertility", g
      )
    }
    for (g in levels(pairs$initial_total_wealth_tertile)) {
      rows[[length(rows) + 1L]] <- summarize_group(
        pairs[initial_total_wealth_tertile == g], pairs$sample[1L], asset,
        "initial_total_wealth_tertile", g
      )
    }
    for (g in levels(pairs$lagged_total_wealth_tertile)) {
      rows[[length(rows) + 1L]] <- summarize_group(
        pairs[lagged_total_wealth_tertile == g], pairs$sample[1L], asset,
        "lagged_total_wealth_tertile", g
      )
    }
  }
  rbindlist(rows, fill = TRUE)
}

run_regressions <- function(pairs) {
  rows <- list()
  for (asset in c("total", "nonhousing", "home_equity")) {
    outcome <- paste0(asset, "_ihs_change_annual")
    initial <- paste0(asset, "_initial_ihs")
    formulas <- list(
      child_slope = as.formula(paste0(
        outcome, " ~ children_capped + ", initial, " + I(", initial,
        "^2) + log_income0 + age0 + own0 | year0"
      )),
      fertility_bins = as.formula(paste0(
        outcome, " ~ i(fertility_bin, ref = '2') + ", initial,
        " + I(", initial, "^2) + log_income0 + age0 + own0 | year0"
      )),
      child_slope_marriage_sensitivity = as.formula(paste0(
        outcome, " ~ children_capped + ", initial, " + I(", initial,
        "^2) + log_income0 + age0 + own0 + married0 + married1 | year0"
      ))
    )
    for (spec in names(formulas)) {
      fit <- feols(
        formulas[[spec]], data = pairs, weights = ~weight0,
        vcov = "hetero", notes = FALSE
      )
      ct <- as.data.table(coeftable(fit), keep.rownames = "term")
      setnames(
        ct, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
        c("estimate", "std_error", "t_value", "p_value")
      )
      ct <- ct[term == "children_capped" | grepl("^fertility_bin::", term)]
      ct[, `:=`(
        sample = pairs$sample[1L], asset = asset, specification = spec,
        outcome = "annualized change in asinh(asset / baseline family income)",
        n_obs = nobs(fit)
      )]
      rows[[length(rows) + 1L]] <- ct
    }
  }
  rbindlist(rows, fill = TRUE)
}

candidate_values <- function(pairs) {
  total_change <- pairs$total_ihs_change_annual
  avg <- weighted_mean_safe(total_change, pairs$weight0)
  mid <- weighted_mean_safe(
    total_change[pairs$lagged_total_wealth_tertile == "middle"],
    pairs$weight0[pairs$lagged_total_wealth_tertile == "middle"]
  )
  top <- weighted_mean_safe(
    total_change[pairs$lagged_total_wealth_tertile == "top"],
    pairs$weight0[pairs$lagged_total_wealth_tertile == "top"]
  )
  fit_data <- pairs[
    is.finite(total_ihs_change_annual) & is.finite(children_capped) &
      is.finite(total_initial_ihs) & is.finite(log_income0) &
      is.finite(age0) & is.finite(own0)
  ]
  fit <- lm(
    total_ihs_change_annual ~ children_capped + total_initial_ihs +
      I(total_initial_ihs^2) + log_income0 + age0 + own0 + factor(year0),
    data = fit_data, weights = weight0
  )
  child_slope <- unname(coef(fit)["children_capped"])
  c(theta0_average_retention = avg,
    theta1_lagged_top_minus_middle_retention = top - mid,
    theta_n_conditional_per_child_slope = child_slope)
}

bootstrap_candidates <- function(pairs, reps = 499L, seed = 20260715L) {
  set.seed(seed)
  point <- candidate_values(pairs)
  boot <- matrix(NA_real_, nrow = reps, ncol = length(point))
  colnames(boot) <- names(point)
  n <- nrow(pairs)
  for (r in seq_len(reps)) {
    idx <- sample.int(n, n, replace = TRUE)
    boot[r, ] <- tryCatch(candidate_values(pairs[idx]), error = function(e) {
      rep(NA_real_, length(point))
    })
  }
  data.table(
    moment = names(point),
    value = as.numeric(point),
    bootstrap_se = apply(boot, 2L, sd, na.rm = TRUE),
    bootstrap_p025 = apply(boot, 2L, quantile, probs = 0.025, na.rm = TRUE),
    bootstrap_p975 = apply(boot, 2L, quantile, probs = 0.975, na.rm = TRUE),
    bootstrap_reps = reps
  )
}

root <- normalizePath(file.path(script_dir(), "..", "..", ".."))
psid_path <- file.path(dirname(root), "PSID", "PSIDSHELF_MOBILITY.dta")
outdir <- file.path(
  root, "code", "data", "psid_followup_mar2026", "output",
  "intergen_bequest_retention_targets"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(psid_path)) stop("Missing PSID shelf: ", psid_path)

message("Reading selected PSID columns...")
raw <- as.data.table(read_dta(
  psid_path,
  col_select = c(
    "ID", "year", "AGEREP", "DEATHYEAR", "RELTOHEAD_", "RELCHINUM",
    "INCFAMR", "NETWORTHR", "NETWORTH2R", "HOMEEQUITYR", "HOMEOWN",
    "IW", "FAMMARRIED"
  )
))
dt <- raw[, .(
  id = as.numeric(ID), year = as.integer(year), age = as.numeric(AGEREP),
  death_year = as.numeric(DEATHYEAR), relation_to_head = as.numeric(RELTOHEAD_),
  n_children = as.numeric(RELCHINUM), income = as.numeric(INCFAMR),
  total_nw = as.numeric(NETWORTHR), nonhousing_nw = as.numeric(NETWORTH2R),
  home_equity = as.numeric(HOMEEQUITYR),
  own = fifelse(as.numeric(HOMEOWN) == 1, 1,
                fifelse(as.numeric(HOMEOWN) == 2, 0, NA_real_)),
  weight = as.numeric(IW), married = as.numeric(FAMMARRIED)
)]
rm(raw)
dt <- dt[
  relation_to_head == 10 & age >= 60 & age <= 82 &
    (is.na(death_year) | year <= death_year)
]

primary <- make_pairs(
  dt, 10L, 62, 68, "primary_10y_one_per_person", one_per_person = TRUE
)
sensitivity_early <- make_pairs(
  dt, 4L, 62, 68, "sensitivity_4y_early_one_per_person", one_per_person = TRUE
)
sensitivity_late <- make_pairs(
  dt, 4L, 68, 72, "sensitivity_4y_late_one_per_person", one_per_person = TRUE
)
if (nrow(primary) == 0L) stop("No primary ten-year transitions found.")
pair_samples <- list(primary, sensitivity_early, sensitivity_late)

pair_counts <- rbindlist(lapply(pair_samples, function(x) {
  x[, .(
    n = .N,
    weighted_mean_age0 = weighted_mean_safe(age0, weight0),
    weighted_mean_age1 = weighted_mean_safe(age1, weight0),
    weighted_mean_children = weighted_mean_safe(n_children, weight0),
    weighted_owner_rate0 = weighted_mean_safe(own0, weight0),
    baseline_year_min = min(year0), baseline_year_max = max(year0),
    wealth_tertile_cut1 = wealth_tertile_cut1[1L],
    wealth_tertile_cut2 = wealth_tertile_cut2[1L],
    lagged_wealth_n = sum(is.finite(lag_total_nw)),
    lagged_wealth_tertile_cut1 = lagged_wealth_tertile_cut1[1L],
    lagged_wealth_tertile_cut2 = lagged_wealth_tertile_cut2[1L]
  ), by = .(sample, fertility_bin)]
}), fill = TRUE)
fwrite(pair_counts, file.path(outdir, "transition_sample_counts.csv"))

summaries <- rbindlist(lapply(pair_samples, summarize_pairs))
fwrite(summaries, file.path(outdir, "retention_summaries.csv"))

regressions <- rbindlist(lapply(pair_samples, run_regressions))
setcolorder(regressions, c(
  "sample", "asset", "specification", "outcome", "term", "estimate",
  "std_error", "t_value", "p_value", "n_obs"
))
fwrite(regressions, file.path(outdir, "controlled_retention_gradients.csv"))

candidates <- bootstrap_candidates(primary)
candidates[, n := c(
  nrow(primary),
  primary[lagged_total_wealth_tertile %in% c("middle", "top"), .N],
  primary[
    is.finite(total_ihs_change_annual) & is.finite(children_capped) &
      is.finite(total_initial_ihs) & is.finite(log_income0) &
      is.finite(age0) & is.finite(own0), .N
  ]
)]
candidates[, `:=`(
  sample = primary$sample[1L],
  outcome = "annualized change in asinh(total wealth / baseline family income)",
  status = "candidate_not_approved"
)]
setcolorder(candidates, c(
  "sample", "moment", "n", "value", "bootstrap_se", "bootstrap_p025",
  "bootstrap_p975", "bootstrap_reps", "outcome", "status"
))
fwrite(candidates, file.path(outdir, "candidate_bequest_moments.csv"))

plot_data <- summaries[
  sample == "primary_10y_one_per_person" & asset == "total" &
    group_type %in% c("completed_fertility", "lagged_total_wealth_tertile")
]
plot_data[, panel := fifelse(
  group_type == "completed_fertility", "Completed fertility",
  "Lagged initial-wealth tertile"
)]
p <- ggplot(plot_data, aes(group_value, ihs_change_annual_mean, group = 1)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.8, color = "#1f5a94") +
  geom_point(size = 2.2, color = "#1f5a94") +
  facet_wrap(~panel, scales = "free_x") +
  labs(
    x = NULL,
    y = "Annualized change in asinh(total wealth / baseline income)"
  ) +
  theme_minimal(base_size = 11)
ggsave(
  file.path(outdir, "total_wealth_retention_gradients.png"),
  p, width = 8.5, height = 4.5, dpi = 180
)

candidate_line <- function(moment) {
  key <- moment
  row <- candidates[candidates$moment == key]
  sprintf("%.5f (bootstrap SE %.5f; 95%% interval [%.5f, %.5f])",
          row$value, row$bootstrap_se, row$bootstrap_p025, row$bootstrap_p975)
}
readme <- c(
  "# PSID late-life wealth-retention targets for the bequest block",
  "",
  "## Primary construction",
  "",
  "- One PSID reference person per family and one transition per person.",
  "- Baseline age 62--68; follow-up exactly ten survey years later.",
  "- Completed-child count must be unchanged across the two observations.",
  "- Wealth and family income are already measured in real 2022 dollars.",
  "- Total wealth is `NETWORTHR`; nonhousing wealth is `NETWORTH2R`; home",
  "  equity is `HOMEEQUITYR`.",
  "- The primary outcome is the annualized change in",
  "  `asinh(wealth / baseline family income)`. This permits zero or negative",
  "  wealth and avoids unstable individual percentage changes.",
  "- Raw wealth changes and the median follow-up/baseline ratio among households",
  "  with more than $10,000 in baseline asset wealth are reported as diagnostics.",
  "",
  "## Candidate identifying moments",
  "",
  paste0("- Overall total-wealth retention (`theta0`): ",
         candidate_line("theta0_average_retention"), "."),
  paste0("- Lagged top-minus-middle initial-wealth retention (`theta1`): ",
         candidate_line("theta1_lagged_top_minus_middle_retention"), "."),
  paste0("- Conditional per-child retention slope (`theta_n`): ",
         candidate_line("theta_n_conditional_per_child_slope"), "."),
  "",
  "The child slope conditions on baseline total wealth (quadratically after an",
  "inverse-hyperbolic-sine transform), baseline income, age, ownership, and",
  "baseline-year fixed effects. This asks whether families with more completed",
  "children retain more wealth from the same approximate starting position; it",
  "does not interpret lower raw wealth among large families as weaker bequests.",
  "Initial-wealth tertiles for the `theta1` candidate are formed using wealth",
  "observed two to six years before baseline. This avoids selecting groups on",
  "the same noisy baseline wealth realization used to measure subsequent change.",
  "",
  "## Interpretation boundary",
  "",
  "These are candidate targets, not approved calibration moments. They are",
  "conditional on surviving and remaining observable as a PSID reference person.",
  "They do not observe estates at death, so an HRS exit-estate moment remains the",
  "preferred direct check on the bequest mechanism. The early and late four-year",
  "results are model-period sensitivities, not the primary lifecycle target.",
  "",
  "## Files",
  "",
  "- `candidate_bequest_moments.csv`: the three proposed total-wealth targets.",
  "- `retention_summaries.csv`: all assets by fertility and initial-wealth bin.",
  "- `controlled_retention_gradients.csv`: conditional fertility estimates.",
  "- `transition_sample_counts.csv`: sample composition and wealth cutoffs.",
  "- `total_wealth_retention_gradients.png`: compact primary-sample plot."
)
writeLines(readme, file.path(outdir, "README.md"))

message("Wrote retention packet to ", outdir)
