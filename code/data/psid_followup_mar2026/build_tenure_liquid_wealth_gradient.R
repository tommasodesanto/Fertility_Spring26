#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

weighted_quantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= p)[1]]
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  sum(x[ok] * w[ok]) / sum(w[ok])
}

script_dir <- function() {
  flag <- commandArgs(trailingOnly = FALSE)
  flag <- flag[startsWith(flag, "--file=")]
  if (length(flag) == 0L) stop("Run with Rscript.")
  dirname(normalizePath(sub("--file=", "", flag[1L], fixed = TRUE)))
}

root <- normalizePath(file.path(script_dir(), "..", "..", ".."))
psid_path <- file.path(dirname(root), "PSID", "PSIDSHELF_MOBILITY.dta")
out_dir <- file.path(root, "code", "data", "psid_followup_mar2026", "output", "tenure_liquid_wealth_gradient")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading selected PSID columns...")
raw <- as.data.table(read_dta(
  psid_path,
  col_select = c(ID, year, AGEREP, DEATHYEAR, RELTOHEAD_, INCFAMR, HOMEOWN, IW,
                 NETWORTH2R, WLTHSAVETOTR, WLTHFUNDTOTR, WLTHODEBTOTR)
))

dt <- raw[, .(
  id = as.numeric(ID), year = as.integer(year), age = as.numeric(AGEREP),
  death_year = as.numeric(DEATHYEAR), reference_person = as.numeric(RELTOHEAD_) == 10,
  income = as.numeric(INCFAMR), own = fifelse(as.numeric(HOMEOWN) == 1, 1,
    fifelse(as.numeric(HOMEOWN) == 2, 0, NA_real_)), weight = as.numeric(IW),
  nonhousing_nw = as.numeric(NETWORTH2R),
  liquid_financial = as.numeric(WLTHSAVETOTR) + as.numeric(WLTHFUNDTOTR) - as.numeric(WLTHODEBTOTR)
)]
rm(raw)
invisible(gc())

dt <- dt[
  year >= 1984 & year <= 2019 & reference_person &
    (is.na(death_year) | year <= death_year) &
    is.finite(weight) & weight > 0 & !is.na(own)
]

entry <- dt[age >= 25 & age <= 30 & income > 1000,
  .(
    n_entry_years = .N,
    weight = mean(weight, na.rm = TRUE),
    own_rate_entry = mean(own, na.rm = TRUE),
    income = mean(income, na.rm = TRUE),
    liquid_financial = mean(liquid_financial, na.rm = TRUE),
    nonhousing_nw = mean(nonhousing_nw, na.rm = TRUE)
  ), by = id]
entry <- entry[own_rate_entry < 0.5 & is.finite(income) & income > 1000]

later <- dt[age >= 31 & age <= 35,
  .(own_rate_3135 = mean(own, na.rm = TRUE)), by = id]
later[, owner_by_35 := as.numeric(own_rate_3135 >= 0.5)]
sample <- merge(entry, later[, .(id, owner_by_35)], by = "id")
sample[, liquid_financial_to_income := liquid_financial / income]
sample[, nonhousing_nw_to_income := nonhousing_nw / income]

make_bins <- function(d, wealth_name) {
  x <- d[[wealth_name]]
  lo <- weighted_quantile(x, d$weight, 0.01)
  hi <- weighted_quantile(x, d$weight, 0.99)
  z <- copy(d[is.finite(get(wealth_name)) & between(get(wealth_name), lo, hi)])
  cuts <- vapply(1:4, function(q) weighted_quantile(z[[wealth_name]], z$weight, q / 5), numeric(1))
  z[, quintile := findInterval(get(wealth_name), cuts) + 1L]
  cells <- z[, .(
    n = .N,
    weight = sum(weight),
    wealth_mean = weighted_mean(get(wealth_name), weight),
    wealth_median = weighted_quantile(get(wealth_name), weight, 0.5),
    owner_by_35_rate = weighted_mean(owner_by_35, weight)
  ), by = quintile]
  setorder(cells, quintile)
  list(sample = z, cells = cells, gradient = cells[quintile == 5, owner_by_35_rate] - cells[quintile == 1, owner_by_35_rate])
}

financial <- make_bins(sample, "liquid_financial_to_income")
nonhousing <- make_bins(sample, "nonhousing_nw_to_income")

bootstrap_gradient <- function(z, reps = 500L, seed = 20260722L) {
  set.seed(seed)
  ids <- z$id
  draws <- numeric(reps)
  for (bb in seq_len(reps)) {
    draw_ids <- sample(ids, length(ids), replace = TRUE)
    counts <- table(draw_ids)
    b <- copy(z[match(names(counts), as.character(id))])
    b[, boot_weight := weight * as.numeric(counts)]
    cuts <- vapply(1:4, function(q) weighted_quantile(b$liquid_financial_to_income, b$boot_weight, q / 5), numeric(1))
    b[, quintile_boot := findInterval(liquid_financial_to_income, cuts) + 1L]
    q1 <- weighted_mean(b[quintile_boot == 1, owner_by_35], b[quintile_boot == 1, boot_weight])
    q5 <- weighted_mean(b[quintile_boot == 5, owner_by_35], b[quintile_boot == 5, boot_weight])
    draws[bb] <- q5 - q1
  }
  c(
    se = sd(draws),
    p025 = unname(quantile(draws, 0.025)),
    p975 = unname(quantile(draws, 0.975))
  )
}

boot <- bootstrap_gradient(financial$sample)
cells <- rbind(
  cbind(wealth_concept = "liquid_financial_to_income", financial$cells),
  cbind(wealth_concept = "nonhousing_nw_to_income_sensitivity", nonhousing$cells),
  fill = TRUE
)

moments <- data.table(
  moment = c(
    "initial_renter_owner_by35_gradient_liquid_financial_q5_q1",
    "initial_renter_owner_by35_gradient_nonhousing_nw_q5_q1"
  ),
  value = c(financial$gradient, nonhousing$gradient),
  standard_error = c(unname(boot["se"]), NA_real_),
  ci_lower = c(unname(boot["p025"]), NA_real_),
  ci_upper = c(unname(boot["p975"]), NA_real_),
  n = c(nrow(financial$sample), nrow(nonhousing$sample)),
  definition = c(
    "Among reference persons who are renters at ages 25-30: weighted ownership-rate difference by ages 31-35 between top and bottom quintiles of (cash/checking + bonds/CDs + stocks/IRAs - other debt)/family income; p1-p99 trimmed",
    "Same transition gradient using PSID NETWORTH2R/family income; sensitivity matching the older files' broad nonhousing-wealth convention"
  ),
  status = c("candidate primary target for tenure_choice_kappa", "sensitivity only; not genuinely liquid wealth")
)

fwrite(cells, file.path(out_dir, "psid_tenure_wealth_quintile_cells.csv"))
fwrite(moments, file.path(out_dir, "psid_tenure_wealth_gradient_moments.csv"))

writeLines(c(
  "# PSID tenure transition by initial liquid wealth",
  "",
  "Primary object: ownership by ages 31--35 among households initially renting",
  "at ages 25--30, split by initial liquid-financial-wealth/income quintile.",
  "This avoids the mechanical contemporaneous effect of home purchase and mortgage",
  "debt on measured wealth. NETWORTH2R is reported only as a sensitivity because",
  "it includes broad nonhousing assets and was misleadingly called liquid wealth",
  "in several older project files.",
  "",
  "The primary standard error uses a 500-draw person bootstrap."
), file.path(out_dir, "README.md"))

message("Wrote PSID tenure-gradient targets to ", out_dir)
