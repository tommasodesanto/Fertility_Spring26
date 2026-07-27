#!/usr/bin/env Rscript

# PSID decomposition of RP/spouse real gross labor earnings into permanent
# (between-household) and transitory/measurement (within-household) variation.
# Survey measurement error inflates the WITHIN component, so the between share
# is a lower bound on permanent earnings dispersion.

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(haven)
})

YEAR_MIN <- 2005L
YEAR_MAX <- 2019L
PROCESS_YEAR_MIN <- 1984L
AGE_MIN <- 25L
AGE_MAX <- 60L
MD_LAGS <- c(0L, 1L, 2L, 4L, 6L, 8L, 10L, 12L, 16L, 20L, 24L, 28L, 32L)
MD_BOOT_REPS <- as.integer(Sys.getenv("PSID_MD_BOOT_REPS", "199"))
MD_BOOT_SEED <- 20260727L
MD_CHECKPOINT_EVERY <- 5L

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
md_out_dir <- file.path(out_dir, "psid_income_fixed_effect_md_20260727")
extract_path <- file.path(md_out_dir, "psid_income_md_extract.dta")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(md_out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(psid_path)) stop("Missing PSID shelf file: ", psid_path)
if (!is.finite(MD_BOOT_REPS) || MD_BOOT_REPS < 5L) {
  stop("PSID_MD_BOOT_REPS must be an integer of at least 5.")
}

if (file.exists(extract_path)) {
  message("Reading the narrow PSID income-process extract...")
  raw <- as.data.table(read_dta(extract_path))
} else {
  message("Reading selected PSID columns from the full shelf...")
  raw <- as.data.table(read_dta(
    psid_path,
    col_select = c(
      "ID", "year", "RELTOHEAD_", "AGEREP", "DEATHYEAR", "EARNINDRRC", "IW"
    )
  ))
}
dt_all <- raw[, .(
  id = as.numeric(ID),
  year = as.integer(year),
  relation_to_head = as.numeric(RELTOHEAD_),
  age = as.numeric(AGEREP),
  death_year = as.numeric(DEATHYEAR),
  # EARNINDRRC is already RP/spouse combined real USD 2022 gross earnings.
  earnings_real_2022 = as.numeric(EARNINDRRC),
  weight = as.numeric(IW)
)]
rm(raw)

dt_all <- dt_all[
  relation_to_head == 10 &
    year >= PROCESS_YEAR_MIN & year <= YEAR_MAX &
    (is.na(death_year) | year <= death_year) &
    age >= AGE_MIN & age <= AGE_MAX &
    is.finite(earnings_real_2022) & earnings_real_2022 > 0 &
    is.finite(weight) & weight > 0
]
if (nrow(dt_all) == 0L) stop("The requested PSID earnings sample is empty.")
if (anyDuplicated(dt_all, by = c("id", "year"))) {
  stop("Duplicate reference-person (id, year) observations in the PSID sample.")
}
dt_all[, log_earnings := log(earnings_real_2022)]
dt <- copy(dt_all[year >= YEAR_MIN])
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

# ---------------------------------------------------------------------------
# Fixed effect + persistent AR(1) + transitory minimum-distance decomposition
# ---------------------------------------------------------------------------
#
# For residual log earnings u_it,
#
#   u_it = a_i + p_it + e_it,    p_it = rho p_i,t-1 + eta_it,
#
# the autocovariance restrictions are
#
#   gamma(0) = var(a) + var(p) + var(e),
#   gamma(k) = var(a) + rho^k var(p),  k > 0.
#
# The long lags identify var(a): a short-lag persistent process can otherwise
# be mistaken for a fixed effect. The fixed-effect variance used by E6b is an
# externally estimated empirical object, not a new calibrated parameter.

message("\nPSID fixed-effect / AR(1) / transitory minimum-distance block")
inc <- copy(dt_all)
inc[, age_i := as.integer(round(age))]
inc <- inc[age_i >= AGE_MIN & age_i <= AGE_MAX]
setorder(inc, id, year)
inc[, rid := .I]
inc[, pid := match(id, unique(id))]

moment_pair_tables <- lapply(MD_LAGS[MD_LAGS > 0L], function(lag_years) {
  merge(
    inc[, .(id, pair_year = year, r1 = rid)],
    inc[, .(id, pair_year = year - lag_years, r2 = rid)],
    by = c("id", "pair_year"),
    all = FALSE
  )[, .(
    r1,
    r2,
    pair_weight_base = sqrt(inc$weight[r1] * inc$weight[r2]),
    pid = inc$pid[r1]
  )]
})
names(moment_pair_tables) <- as.character(MD_LAGS[MD_LAGS > 0L])

income_autocovariances <- function(person_frequency = NULL) {
  if (is.null(person_frequency)) {
    row_frequency <- rep(1, nrow(inc))
  } else {
    row_frequency <- person_frequency[inc$pid]
  }
  analysis_weight <- inc$weight * row_frequency
  keep <- analysis_weight > 0
  fit <- feols(
    log_earnings ~ 1 | age_i + year,
    data = inc[keep],
    weights = analysis_weight[keep],
    notes = FALSE
  )
  residual <- rep(NA_real_, nrow(inc))
  residual[keep] <- as.numeric(resid(fit))

  moments <- rep(NA_real_, length(MD_LAGS))
  pair_counts <- rep(NA_real_, length(MD_LAGS))
  moments[1L] <- weighted_var_safe(residual, analysis_weight)
  pair_counts[1L] <- sum(row_frequency)
  for (jj in seq_along(moment_pair_tables)) {
    pairs_j <- moment_pair_tables[[jj]]
    frequency_j <- if (is.null(person_frequency)) {
      rep(1, nrow(pairs_j))
    } else {
      person_frequency[pairs_j$pid]
    }
    pair_weight <- pairs_j$pair_weight_base * frequency_j
    moments[jj + 1L] <- weighted_cov_safe(
      residual[pairs_j$r1],
      residual[pairs_j$r2],
      pair_weight
    )
    pair_counts[jj + 1L] <- sum(frequency_j)
  }
  list(moments = moments, pair_counts = pair_counts)
}

predict_income_autocovariances <- function(parameters, lags = MD_LAGS) {
  fixed_variance <- parameters[["fixed_effect_variance"]]
  persistent_variance <- parameters[["persistent_variance"]]
  transitory_variance <- parameters[["transitory_variance"]]
  rho <- parameters[["rho_annual"]]
  predicted <- fixed_variance + persistent_variance * rho^lags
  predicted[lags == 0L] <- predicted[lags == 0L] + transitory_variance
  predicted
}

fit_income_md <- function(
    empirical_moments,
    lags,
    weight_matrix,
    restrict_fixed_effect = FALSE,
    start = NULL) {
  objective <- function(raw_parameters) {
    if (restrict_fixed_effect) {
      parameters <- c(
        fixed_effect_variance = 0,
        persistent_variance = raw_parameters[1L],
        transitory_variance = raw_parameters[2L],
        rho_annual = raw_parameters[3L]
      )
    } else {
      parameters <- setNames(
        raw_parameters,
        c(
          "fixed_effect_variance", "persistent_variance",
          "transitory_variance", "rho_annual"
        )
      )
    }
    gap <- empirical_moments -
      predict_income_autocovariances(parameters, lags)
    as.numeric(crossprod(gap, weight_matrix %*% gap))
  }
  objective_gradient <- function(raw_parameters) {
    if (restrict_fixed_effect) {
      parameters <- c(
        fixed_effect_variance = 0,
        persistent_variance = raw_parameters[1L],
        transitory_variance = raw_parameters[2L],
        rho_annual = raw_parameters[3L]
      )
    } else {
      parameters <- setNames(
        raw_parameters,
        c(
          "fixed_effect_variance", "persistent_variance",
          "transitory_variance", "rho_annual"
        )
      )
    }
    rho <- parameters[["rho_annual"]]
    persistent_variance <- parameters[["persistent_variance"]]
    gap <- empirical_moments -
      predict_income_autocovariances(parameters, lags)
    prediction_jacobian <- cbind(
      fixed_effect_variance = rep(1, length(lags)),
      persistent_variance = rho^lags,
      transitory_variance = as.numeric(lags == 0L),
      rho_annual = ifelse(
        lags == 0L,
        0,
        persistent_variance * lags * rho^(lags - 1L)
      )
    )
    if (restrict_fixed_effect) {
      prediction_jacobian <- prediction_jacobian[, -1L, drop = FALSE]
    }
    as.numeric(-2 * crossprod(
      prediction_jacobian,
      weight_matrix %*% gap
    ))
  }

  total_variance <- max(empirical_moments[lags == 0L], 0.05)
  long_covariance <- max(tail(empirical_moments, 1L), 0)
  candidate_starts <- list(
    c(long_covariance, max(empirical_moments[2L] - long_covariance, 0.05),
      max(total_variance - empirical_moments[2L], 0.05), 0.97),
    c(0, max(empirical_moments[2L], 0.05),
      max(total_variance - empirical_moments[2L], 0.05), 0.97),
    c(0.10 * total_variance, 0.65 * total_variance,
      0.25 * total_variance, 0.93),
    c(0.35 * total_variance, 0.45 * total_variance,
      0.20 * total_variance, 0.85)
  )
  if (!is.null(start)) candidate_starts <- c(list(start), candidate_starts)

  if (restrict_fixed_effect) {
    candidate_starts <- lapply(candidate_starts, function(x) x[c(2L, 3L, 4L)])
    lower <- c(0, 0, 0.50)
    upper <- c(3 * total_variance, 3 * total_variance, 0.9995)
  } else {
    lower <- c(0, 0, 0, 0.50)
    upper <- c(
      3 * total_variance, 3 * total_variance, 3 * total_variance, 0.9995
    )
  }

  best <- NULL
  for (candidate_start in candidate_starts) {
    candidate_start <- pmin(pmax(candidate_start, lower + 1e-10), upper - 1e-10)
    candidate <- optim(
      candidate_start,
      objective,
      gr = objective_gradient,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 5000, factr = 1e4, pgtol = 1e-10)
    )
    if (is.null(best) || candidate$value < best$value) best <- candidate
  }
  refined <- nlminb(
    start = best$par,
    objective = objective,
    gradient = objective_gradient,
    lower = lower,
    upper = upper,
    control = list(
      iter.max = 10000,
      eval.max = 20000,
      rel.tol = 1e-12,
      x.tol = 1e-10
    )
  )
  if (
      refined$convergence == 0L ||
      refined$objective < best$value - 1e-8) {
    best <- list(
      par = refined$par,
      value = refined$objective,
      convergence = refined$convergence
    )
  }
  if (best$convergence != 0L) {
    rho_lower <- lower[length(lower)]
    rho_upper <- upper[length(upper)]
    to_unconstrained <- function(raw_parameters) {
      variance_parameters <- pmax(
        raw_parameters[-length(raw_parameters)],
        1e-10
      )
      rho_share <- (
        raw_parameters[length(raw_parameters)] - rho_lower
      ) / (rho_upper - rho_lower)
      rho_share <- pmin(pmax(rho_share, 1e-10), 1 - 1e-10)
      c(log(variance_parameters), qlogis(rho_share))
    }
    from_unconstrained <- function(unconstrained_parameters) {
      rho_share <- plogis(tail(unconstrained_parameters, 1L))
      c(
        exp(head(unconstrained_parameters, -1L)),
        rho_lower + (rho_upper - rho_lower) * rho_share
      )
    }
    transformed_objective <- function(unconstrained_parameters) {
      objective(from_unconstrained(unconstrained_parameters))
    }
    transformed_gradient <- function(unconstrained_parameters) {
      raw_parameters <- from_unconstrained(unconstrained_parameters)
      rho_share <- plogis(tail(unconstrained_parameters, 1L))
      transform_derivative <- c(
        head(raw_parameters, -1L),
        (rho_upper - rho_lower) * rho_share * (1 - rho_share)
      )
      objective_gradient(raw_parameters) * transform_derivative
    }
    transformed <- optim(
      to_unconstrained(best$par),
      transformed_objective,
      gr = transformed_gradient,
      method = "BFGS",
      control = list(maxit = 10000, reltol = 1e-12)
    )
    transformed_raw <- from_unconstrained(transformed$par)
    if (
        transformed$convergence == 0L ||
        transformed$value < best$value - 1e-8) {
      best <- list(
        par = transformed_raw,
        value = transformed$value,
        convergence = transformed$convergence
      )
    }
  }

  if (restrict_fixed_effect) {
    parameters <- c(
      fixed_effect_variance = 0,
      persistent_variance = best$par[1L],
      transitory_variance = best$par[2L],
      rho_annual = best$par[3L]
    )
  } else {
    parameters <- setNames(
      best$par,
      c(
        "fixed_effect_variance", "persistent_variance",
        "transitory_variance", "rho_annual"
      )
    )
  }
  list(parameters = parameters, objective = best$value, convergence = best$convergence)
}

make_md_weight <- function(moment_draws, mode = "shrink_full") {
  draw_covariance <- cov(moment_draws, use = "pairwise.complete.obs")
  diagonal_variance <- pmax(diag(draw_covariance), 1e-10)
  if (mode == "diagonal" || nrow(moment_draws) < ncol(moment_draws) + 5L) {
    return(diag(1 / diagonal_variance))
  }
  scale_matrix <- diag(sqrt(diagonal_variance))
  inverse_scale <- diag(1 / sqrt(diagonal_variance))
  correlation <- inverse_scale %*% draw_covariance %*% inverse_scale
  # Ten-percent shrinkage makes the bootstrap covariance invertible while
  # retaining the joint information in the long-lag moment vector.
  regularized_correlation <- 0.90 * correlation + 0.10 * diag(ncol(correlation))
  inverse_scale %*% solve(regularized_correlation) %*% inverse_scale
}

point_moments <- income_autocovariances()
bootstrap_checkpoint_path <- file.path(
  md_out_dir, "bootstrap_autocovariance_checkpoint.csv"
)
moment_column_names <- paste0("lag_", MD_LAGS)
if (file.exists(bootstrap_checkpoint_path)) {
  bootstrap_checkpoint <- fread(bootstrap_checkpoint_path)
  expected_names <- c("replication", moment_column_names)
  if (!identical(names(bootstrap_checkpoint), expected_names)) {
    stop("Existing MD bootstrap checkpoint has an incompatible schema.")
  }
  bootstrap_checkpoint <- bootstrap_checkpoint[replication <= MD_BOOT_REPS]
} else {
  bootstrap_checkpoint <- data.table()
}

n_persons_md <- max(inc$pid)
completed_reps <- if (nrow(bootstrap_checkpoint)) {
  bootstrap_checkpoint$replication
} else {
  integer(0)
}
for (replication in setdiff(seq_len(MD_BOOT_REPS), completed_reps)) {
  set.seed(MD_BOOT_SEED + replication)
  person_frequency <- tabulate(
    sample.int(n_persons_md, n_persons_md, replace = TRUE),
    nbins = n_persons_md
  )
  draw <- income_autocovariances(person_frequency)$moments
  new_checkpoint_row <- as.data.table(as.list(
    c(replication = replication, draw)
  ))
  setnames(new_checkpoint_row, c("replication", moment_column_names))
  bootstrap_checkpoint <- rbind(
    bootstrap_checkpoint,
    new_checkpoint_row,
    use.names = TRUE
  )
  setorder(bootstrap_checkpoint, replication)
  if (
      replication %% MD_CHECKPOINT_EVERY == 0L ||
      replication == MD_BOOT_REPS) {
    fwrite(bootstrap_checkpoint, bootstrap_checkpoint_path)
    message("  MD bootstrap checkpoint ", replication, " / ", MD_BOOT_REPS)
  }
}

bootstrap_moment_matrix <- as.matrix(
  bootstrap_checkpoint[order(replication), ..moment_column_names]
)
md_weight <- make_md_weight(bootstrap_moment_matrix, "shrink_full")
md_weight_diagonal <- make_md_weight(bootstrap_moment_matrix, "diagonal")

md_fit <- fit_income_md(
  point_moments$moments, MD_LAGS, md_weight, restrict_fixed_effect = FALSE
)
md_fit_no_fixed <- fit_income_md(
  point_moments$moments, MD_LAGS, md_weight, restrict_fixed_effect = TRUE
)

bootstrap_parameter_draws <- matrix(
  NA_real_,
  nrow = MD_BOOT_REPS,
  ncol = 4L,
  dimnames = list(NULL, names(md_fit$parameters))
)
bootstrap_parameter_convergence <- rep(NA_integer_, MD_BOOT_REPS)
for (replication in seq_len(MD_BOOT_REPS)) {
  draw_fit <- fit_income_md(
    bootstrap_moment_matrix[replication, ],
    MD_LAGS,
    md_weight,
    restrict_fixed_effect = FALSE,
    start = unname(md_fit$parameters)
  )
  bootstrap_parameter_draws[replication, ] <- draw_fit$parameters
  bootstrap_parameter_convergence[replication] <- draw_fit$convergence
}

parameter_transform <- function(parameters) {
  fixed_variance <- parameters[["fixed_effect_variance"]]
  persistent_variance <- parameters[["persistent_variance"]]
  transitory_variance <- parameters[["transitory_variance"]]
  rho <- parameters[["rho_annual"]]
  total_variance <- fixed_variance + persistent_variance + transitory_variance
  c(
    parameters,
    fixed_effect_sd = sqrt(fixed_variance),
    persistent_sd = sqrt(persistent_variance),
    transitory_sd = sqrt(transitory_variance),
    sigma_eta_annual = sqrt(persistent_variance * (1 - rho^2)),
    rho_4yr = rho^4,
    fixed_share_of_residual_variance = fixed_variance / total_variance
  )
}

point_parameter_vector <- parameter_transform(md_fit$parameters)
bootstrap_parameter_transforms <- t(apply(
  bootstrap_parameter_draws,
  1L,
  function(x) parameter_transform(setNames(x, names(md_fit$parameters)))
))
parameter_table <- data.table(
  parameter = names(point_parameter_vector),
  estimate = as.numeric(point_parameter_vector),
  bootstrap_se = apply(bootstrap_parameter_transforms, 2L, sd, na.rm = TRUE),
  bootstrap_p025 = apply(
    bootstrap_parameter_transforms, 2L, quantile, 0.025,
    na.rm = TRUE, names = FALSE
  ),
  bootstrap_p975 = apply(
    bootstrap_parameter_transforms, 2L, quantile, 0.975,
    na.rm = TRUE, names = FALSE
  ),
  bootstrap_reps = MD_BOOT_REPS,
  near_lower_bound = c(
    md_fit$parameters <= c(1e-6, 1e-6, 1e-6, 0.5005),
    rep(FALSE, length(point_parameter_vector) - 4L)
  ),
  estimator = if (MD_BOOT_REPS >= length(MD_LAGS) + 5L) {
    "minimum_distance_10pct_shrunk_bootstrap_covariance"
  } else {
    "smoke_diagonal_bootstrap_covariance"
  }
)
fwrite(parameter_table, file.path(md_out_dir, "md_parameter_estimates.csv"))

fitted_moments <- predict_income_autocovariances(md_fit$parameters, MD_LAGS)
fitted_no_fixed <- predict_income_autocovariances(
  md_fit_no_fixed$parameters, MD_LAGS
)
moment_table <- data.table(
  lag_years = MD_LAGS,
  empirical_autocovariance = point_moments$moments,
  bootstrap_se = apply(bootstrap_moment_matrix, 2L, sd, na.rm = TRUE),
  fitted_fixed_ar1_transitory = fitted_moments,
  fitted_ar1_transitory_no_fixed = fitted_no_fixed,
  gap_full = point_moments$moments - fitted_moments,
  n_pairs = point_moments$pair_counts
)
fwrite(moment_table, file.path(md_out_dir, "md_autocovariance_fit.csv"))

sensitivity_definitions <- list(
  full_lags_shrink_covariance = list(max_lag = 32L, weight = md_weight),
  full_lags_diagonal_covariance = list(max_lag = 32L, weight = md_weight_diagonal),
  through_16_years = list(max_lag = 16L, weight = md_weight),
  through_24_years = list(max_lag = 24L, weight = md_weight)
)
sensitivity_table <- rbindlist(lapply(
  names(sensitivity_definitions),
  function(label) {
    definition <- sensitivity_definitions[[label]]
    selected <- which(MD_LAGS <= definition$max_lag)
    selected_weight <- definition$weight[selected, selected, drop = FALSE]
    fit <- fit_income_md(
      point_moments$moments[selected],
      MD_LAGS[selected],
      selected_weight,
      restrict_fixed_effect = FALSE
    )
    transformed <- parameter_transform(fit$parameters)
    data.table(
      specification = label,
      max_lag = definition$max_lag,
      parameter = names(transformed),
      estimate = as.numeric(transformed),
      objective = fit$objective,
      convergence = fit$convergence
    )
  }
))
fwrite(sensitivity_table, file.path(md_out_dir, "md_sensitivity.csv"))

bootstrap_parameter_table <- as.data.table(bootstrap_parameter_transforms)
bootstrap_parameter_table[, replication := seq_len(.N)]
bootstrap_parameter_table[, convergence := bootstrap_parameter_convergence]
setcolorder(bootstrap_parameter_table, "replication")
fwrite(
  bootstrap_parameter_table,
  file.path(md_out_dir, "bootstrap_parameter_draws.csv")
)

fit_summary <- data.table(
  specification = c("fixed_ar1_transitory", "ar1_transitory_no_fixed"),
  objective = c(md_fit$objective, md_fit_no_fixed$objective),
  objective_gap_vs_unrestricted = c(
    0, md_fit_no_fixed$objective - md_fit$objective
  ),
  moments = length(MD_LAGS),
  free_parameters = c(4L, 3L),
  overidentifying_degrees_of_freedom = c(
    length(MD_LAGS) - 4L, length(MD_LAGS) - 3L
  ),
  n_person_years = nrow(inc),
  n_persons = uniqueN(inc$id),
  year_min = min(inc$year),
  year_max = max(inc$year),
  converged = c(md_fit$convergence == 0L, md_fit_no_fixed$convergence == 0L),
  bootstrap_converged_share = c(
    mean(bootstrap_parameter_convergence == 0L),
    NA_real_
  )
)
fwrite(fit_summary, file.path(md_out_dir, "md_fit_summary.csv"))

png(
  file.path(md_out_dir, "md_autocovariance_fit.png"),
  width = 1800,
  height = 1100,
  res = 180
)
plot(
  moment_table$lag_years,
  moment_table$empirical_autocovariance,
  type = "b",
  pch = 19,
  xlab = "Exact calendar-year lag",
  ylab = "Residual log-earnings autocovariance",
  main = "PSID gross labor earnings: long-lag minimum-distance fit"
)
arrows(
  moment_table$lag_years,
  moment_table$empirical_autocovariance - 1.96 * moment_table$bootstrap_se,
  moment_table$lag_years,
  moment_table$empirical_autocovariance + 1.96 * moment_table$bootstrap_se,
  angle = 90,
  code = 3,
  length = 0.03
)
lines(
  moment_table$lag_years,
  moment_table$fitted_fixed_ar1_transitory,
  type = "b",
  pch = 17,
  col = "#1b6ca8"
)
lines(
  moment_table$lag_years,
  moment_table$fitted_ar1_transitory_no_fixed,
  type = "b",
  pch = 15,
  col = "#b23a48"
)
legend(
  "topright",
  legend = c("PSID moments", "fixed + AR(1) + transitory", "AR(1) + transitory"),
  col = c("black", "#1b6ca8", "#b23a48"),
  pch = c(19, 17, 15),
  lty = 1,
  bty = "n"
)
dev.off()

readme_lines <- c(
  "# PSID fixed-effect / AR(1) / transitory earnings decomposition",
  "",
  "This packet extends `build_psid_income_between_within.R`. The sample is",
  "PSID reference persons ages 25--60, alive in the observation year, with",
  "positive `IW` and positive RP/spouse combined real gross labor earnings",
  "(`EARNINDRRC`), 1984--2019. Log earnings are residualized on integer-age",
  "and year fixed effects with `IW` weights.",
  "",
  "The model is `u_it = a_i + p_it + e_it`, with",
  "`p_it = rho p_i,t-1 + eta_it`. It fits residual autocovariances at exact",
  "calendar-year lags 0, 1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 28, and 32.",
  "The long lags distinguish a fixed effect from a highly persistent AR(1).",
  "Pair weights are the geometric mean of the two survey weights.",
  "",
  sprintf(
    "Inference uses %d deterministic person-cluster bootstrap replications",
    MD_BOOT_REPS
  ),
  "(seed 20260727), re-running residualization and moment construction in",
  "each replication. The minimum-distance weight is the inverse bootstrap",
  "moment covariance after shrinking its correlation matrix 10 percent",
  "toward the identity. The checkpoint is resume-safe at five-rep intervals.",
  "",
  "The estimated fixed-effect variance is an external empirical restriction",
  "for E6b, not a freely calibrated structural parameter. Sensitivity rows",
  "change the maximum lag and the covariance weighting; the no-fixed-effect",
  "fit is reported as a nested diagnostic.",
  "",
  "The narrow extract can be rebuilt without loading the full shelf into R:",
  "`stata-mp -b do extract_psid_income_md.do <shelf.dta> <extract.dta>`."
)
writeLines(readme_lines, file.path(md_out_dir, "README.md"))

message("\nMinimum-distance parameter estimates")
print(parameter_table)
message("Wrote fixed-effect MD packet to ", md_out_dir)
