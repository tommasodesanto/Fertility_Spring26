#!/usr/bin/env Rscript

# PSID decomposition of RP/spouse real gross labor earnings into permanent
# (between-household) and transitory/measurement (within-household) variation.
# Survey measurement error inflates the WITHIN component, so the between share
# is a lower bound on permanent earnings dispersion.

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

YEAR_MIN <- 2005L
YEAR_MAX <- 2019L
AGE_MIN <- 25L
AGE_MAX <- 60L

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- args[startsWith(args, "--file=")]
  if (length(flag) == 0L) stop("Run this file with Rscript.")
  dirname(normalizePath(sub("--file=", "", flag[1L], fixed = TRUE)))
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_var_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  mean_x <- sum(w * x) / sum(w)
  sum(w * (x - mean_x)^2) / sum(w)
}

weighted_quantile_safe <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  order_index <- order(x)
  x <- x[order_index]; w <- w[order_index]
  x[which(cumsum(w) / sum(w) >= p)[1L]]
}

weighted_cov_safe <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  mx <- sum(w * x) / sum(w)
  my <- sum(w * y) / sum(w)
  sum(w * (x - mx) * (y - my)) / sum(w)
}

lognormal_p90_p50 <- function(variance) {
  if (!is.finite(variance) || variance < 0) return(NA_real_)
  exp(1.2816 * sqrt(variance))
}

script_path <- script_dir()
repo_root <- normalizePath(file.path(script_path, "..", "..", ".."))
psid_path <- file.path(dirname(repo_root), "PSID", "PSIDSHELF_MOBILITY.dta")
out_dir <- file.path(repo_root, "code", "data", "psid_followup_mar2026", "output")
out_path <- file.path(out_dir, "psid_income_between_within.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(psid_path)) stop("Missing PSID shelf file: ", psid_path)

message("Reading selected PSID columns...")
raw <- as.data.table(read_dta(
  psid_path,
  col_select = c("ID", "year", "RELTOHEAD_", "AGEREP", "EARNINDRRC", "IW")
))
dt <- raw[, .(
  id = as.numeric(ID),
  year = as.integer(year),
  relation_to_head = as.numeric(RELTOHEAD_),
  age = as.numeric(AGEREP),
  # EARNINDRRC is already RP/spouse combined real USD 2022 gross earnings.
  earnings_real_2022 = as.numeric(EARNINDRRC),
  weight = as.numeric(IW)
)]
rm(raw)

dt <- dt[
  relation_to_head == 10 &
    year >= YEAR_MIN & year <= YEAR_MAX &
    age >= AGE_MIN & age <= AGE_MAX &
    is.finite(earnings_real_2022) & earnings_real_2022 > 0 &
    is.finite(weight) & weight > 0
]
if (nrow(dt) == 0L) stop("The requested PSID earnings sample is empty.")
if (anyDuplicated(dt, by = c("id", "year"))) {
  stop("Duplicate reference-person (id, year) observations in the PSID sample.")
}
dt[, log_earnings := log(earnings_real_2022)]

by_wave <- dt[, .(
  statistic = c("p90_p50", "log_variance"),
  value = c(
    weighted_quantile_safe(earnings_real_2022, weight, 0.90) /
      weighted_quantile_safe(earnings_real_2022, weight, 0.50),
    weighted_var_safe(log_earnings, weight)
  ),
  unweighted_n = .N,
  weighted_n = sum(weight)
), by = year]
by_wave[, `:=`(section = "cross_section", component = "all", note =
  "IW-weighted; EARNINDRRC is RP/spouse combined real USD 2022 gross earnings")]

# This is the requested pooled cross-sectional log variance net of wave effects.
pooled_year_fit <- lm(log_earnings ~ factor(year), data = dt, weights = weight)
pooled_year_residual <- resid(pooled_year_fit)
pooled <- data.table(
  section = "cross_section",
  statistic = c("p90_p50", "log_variance_year_effects_removed"),
  component = "all",
  year = NA_integer_,
  value = c(
    weighted_quantile_safe(dt$earnings_real_2022, dt$weight, 0.90) /
      weighted_quantile_safe(dt$earnings_real_2022, dt$weight, 0.50),
    weighted_var_safe(pooled_year_residual, dt$weight)
  ),
  unweighted_n = nrow(dt),
  weighted_n = sum(dt$weight),
  note = c(
    "IW-weighted pooled distribution; nominal scale is common because EARNINDRRC is real USD 2022",
    "IW-weighted variance after weighted OLS removal of wave dummies"
  )
)

# Residualize log earnings on a quartic in age and wave dummies with IW weights.
residual_fit <- lm(
  log_earnings ~ age + I(age^2) + I(age^3) + I(age^4) + factor(year),
  data = dt,
  weights = weight
)
dt[, residual := as.numeric(resid(residual_fit))]

# A household's residual mean uses its IW weights.  Repeating that mean on each
# family-year makes total variance decompose exactly into weighted between plus
# weighted within variance over the restricted balanced-enough panel.
panel <- dt[, .(
  n_waves = .N,
  household_residual_mean = weighted_mean_safe(residual, weight),
  household_within_variance = weighted_var_safe(residual, weight),
  household_weight = sum(weight)
), by = id][n_waves >= 3L]
restricted <- merge(dt, panel, by = "id", all = FALSE)
restricted[, within_deviation := residual - household_residual_mean]

total_residual_variance <- weighted_var_safe(restricted$residual, restricted$weight)
naive_between_variance <- weighted_var_safe(
  restricted$household_residual_mean, restricted$weight
)
within_variance <- weighted_var_safe(restricted$within_deviation, restricted$weight)
mean_within_variance_over_n <- weighted_mean_safe(
  panel$household_within_variance / panel$n_waves, panel$household_weight
)
corrected_between_variance <- naive_between_variance - mean_within_variance_over_n

component_table <- data.table(
  section = "between_within",
  statistic = "residual_log_variance",
  component = c(
    "total_restricted_panel", "between_naive", "between_corrected",
    "within", "mean_within_variance_over_n", "between_share_naive",
    "between_share_corrected"
  ),
  year = NA_integer_,
  value = c(
    total_residual_variance, naive_between_variance, corrected_between_variance,
    within_variance, mean_within_variance_over_n,
    naive_between_variance / total_residual_variance,
    corrected_between_variance / total_residual_variance
  ),
  unweighted_n = nrow(restricted),
  weighted_n = sum(restricted$weight),
  note = c(
    rep("Residuals from IW-weighted quartic-age and wave-dummy OLS; households observed in at least 3 waves", 5L),
    "Naive between variance divided by total restricted-panel residual variance",
    "Corrected between variance divided by total restricted-panel residual variance; measurement error inflates WITHIN, so this share is a lower bound on permanent dispersion"
  )
)

# Adjacent PSID observations in this period are two years apart. Pair weights
# are the geometric mean of the two IW weights, matching the companion income
# process builder's pair-weight convention.
left <- restricted[, .(id, year, within_deviation, weight)]
right <- restricted[, .(
  id, year = year - 2L, within_deviation_lag = within_deviation,
  weight_lag = weight
)]
pairs <- merge(left, right, by = c("id", "year"), all = FALSE)
pairs[, pair_weight := sqrt(weight * weight_lag)]
rho_within_lag2 <- weighted_cov_safe(
  pairs$within_deviation, pairs$within_deviation_lag, pairs$pair_weight
) / sqrt(
  weighted_var_safe(pairs$within_deviation, pairs$pair_weight) *
    weighted_var_safe(pairs$within_deviation_lag, pairs$pair_weight)
)
autocorrelation <- data.table(
  section = "autocorrelation",
  statistic = "within_deviation_autocorrelation",
  component = "within",
  year = NA_integer_,
  value = rho_within_lag2,
  unweighted_n = nrow(pairs),
  weighted_n = sum(pairs$pair_weight),
  note = "Correlation of household-mean-demeaned residuals across exact two-year adjacent waves; pair weight is sqrt(IW_t * IW_t-2)"
)

lognormal_table <- data.table(
  section = "lognormal_implication",
  statistic = "implied_p90_p50",
  component = c("total_restricted_panel", "between_naive", "between_corrected", "within"),
  year = NA_integer_,
  value = vapply(
    c(total_residual_variance, naive_between_variance, corrected_between_variance, within_variance),
    lognormal_p90_p50,
    numeric(1)
  ),
  unweighted_n = nrow(restricted),
  weighted_n = sum(restricted$weight),
  note = "exp(1.2816 * sd) under lognormality"
)

setcolorder(by_wave, c("section", "statistic", "component", "year", "value", "unweighted_n", "weighted_n", "note"))
out <- rbindlist(list(by_wave, pooled, component_table, autocorrelation, lognormal_table),
                 use.names = TRUE, fill = TRUE)
fwrite(out, out_path)

summary_block <- data.table(
  metric = c(
    "total_residual_variance", "between_naive_variance", "between_corrected_variance",
    "within_variance", "between_share_corrected", "within_autocorrelation_lag2",
    "implied_p90_p50_total", "implied_p90_p50_between_naive",
    "implied_p90_p50_between_corrected", "implied_p90_p50_within"
  ),
  value = c(
    total_residual_variance, naive_between_variance, corrected_between_variance,
    within_variance, corrected_between_variance / total_residual_variance,
    rho_within_lag2,
    lognormal_p90_p50(total_residual_variance), lognormal_p90_p50(naive_between_variance),
    lognormal_p90_p50(corrected_between_variance), lognormal_p90_p50(within_variance)
  )
)
message("\nPSID between/within earnings summary")
print(summary_block)
message("Wrote ", out_path)
