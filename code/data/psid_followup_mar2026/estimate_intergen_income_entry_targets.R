#!/usr/bin/env Rscript

# Income-process, entry-wealth, composition-share, and old-age ownership
# targets for the intergen one-market model, from PSIDSHELF_MOBILITY.
#
# Blocks
#   1. Family log-income process (persistent + transitory) for reference
#      persons aged 25-60, all waves 1984-2019.
#   2. Entry nonhousing wealth / income quintile nodes for reference persons
#      aged 18-24 (and 18-26 fallback), replicating the construction behind
#      PSID_ENTRY_WEALTH_RATIO_NODES_2535 (childless via RELCHIREP == 0,
#      renter via HOMEOWN == 2, INCFAMR > 1000, IW weights).
#   3. IW-weighted share of reference-person family-years aged 65-75 with
#      NETWORTH2R / INCFAMR >= 1.0 (thresholds 0.5 and 2.0 as robustness).
#   4. IW-weighted homeownership rate of reference persons aged 65-75.
#
# Conventions follow audit_intergen_bequest_family_size_targets.R:
# reference persons RELTOHEAD_ == 10, age AGEREP, weight IW, income INCFAMR,
# nonhousing wealth NETWORTH2R, years 1984-2019, alive filter on DEATHYEAR,
# person-cluster bootstrap (resample unique person ids, reweight family-year
# rows by draw frequency), seed 20260715.

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(haven)
})

BOOT_REPS_DEFAULT <- 499L
BOOT_SEED <- 20260715L

MODEL_RHO_ANNUAL <- 0.9601845894041878
MODEL_SIGMA_ETA_ANNUAL <- 0.06453733259357768
MODEL_STATIONARY_VAR <- 0.05337
LEGACY_OLD_AGE_OWN_RATE <- 0.76426097

STORED_ENTRY_NODES_2535 <- c(-2.51940697, -0.07907025, 0.10228762, 0.35287169, 3.03955200)
STORED_ENTRY_WEIGHTS_2535 <- c(0.2000027, 0.2000966, 0.2000024, 0.1998877, 0.2000107)
STORED_ENTRY_MEAN_2535 <- 0.17922556
STORED_ENTRY_MEDIAN_2535 <- 0.09996729

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
  x[which(cumsum(w) / sum(w) >= prob)[1]]
}

weighted_var_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  m <- sum(w * x) / sum(w)
  sum(w * (x - m)^2) / sum(w)
}

weighted_cov_pairs <- function(u1, u2, pw) {
  ok <- is.finite(u1) & is.finite(u2) & is.finite(pw) & pw > 0
  if (!any(ok)) return(NA_real_)
  u1 <- u1[ok]
  u2 <- u2[ok]
  pw <- pw[ok]
  s <- sum(pw)
  m1 <- sum(pw * u1) / s
  m2 <- sum(pw * u2) / s
  sum(pw * (u1 - m1) * (u2 - m2)) / s
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[startsWith(args, "--file=")]
  if (length(file_arg) == 0L) stop("Run this file with Rscript.")
  dirname(normalizePath(sub("--file=", "", file_arg[1L], fixed = TRUE)))
}

root <- normalizePath(file.path(script_dir(), "..", "..", ".."))
psid_path <- file.path(dirname(root), "PSID", "PSIDSHELF_MOBILITY.dta")
outdir <- file.path(
  root, "code", "data", "psid_followup_mar2026", "output",
  "intergen_income_entry_targets_20260716"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(psid_path)) stop("Missing PSID shelf: ", psid_path)

message("Reading selected PSID columns...")
raw <- as.data.table(read_dta(
  psid_path,
  col_select = c(
    "ID", "year", "AGEREP", "DEATHYEAR", "RELTOHEAD_",
    "RELCHIREP", "RELCHINUM", "INCFAMR", "NETWORTH2R", "HOMEOWN", "IW"
  )
))

dt <- raw[, .(
  id = as.numeric(ID),
  year = as.integer(year),
  age = as.numeric(AGEREP),
  death_year = as.numeric(DEATHYEAR),
  reference_person = as.numeric(RELTOHEAD_) == 10,
  n_children_rep = as.numeric(RELCHIREP),
  n_children_num = as.numeric(RELCHINUM),
  income = as.numeric(INCFAMR),
  nonhousing_nw = as.numeric(NETWORTH2R),
  own = fifelse(
    as.numeric(HOMEOWN) == 1, 1,
    fifelse(as.numeric(HOMEOWN) == 2, 0, NA_real_)
  ),
  weight = as.numeric(IW)
)]
rm(raw)
invisible(gc())

# Common base filters shared by every block: wave range, alive in the
# observation year, positive longitudinal weight. Non-reference rows are kept
# only for the replication check of the original 2535 construction.
dt <- dt[
  year >= 1984 & year <= 2019 &
    (is.na(death_year) | year <= death_year) &
    is.finite(weight) & weight > 0
]

data_problems <- character(0)

person_bootstrap <- function(sample_dt, stat_fun, reps, seed = BOOT_SEED) {
  # Person-cluster bootstrap: resample unique person ids with replacement and
  # reweight family-year rows by the draw frequency of their person id.
  ids <- unique(sample_dt$id)
  pid <- match(sample_dt$id, ids)
  point <- stat_fun(sample_dt$weight)
  draws <- matrix(NA_real_, nrow = reps, ncol = length(point),
                  dimnames = list(NULL, names(point)))
  set.seed(seed)
  for (bb in seq_len(reps)) {
    freq <- tabulate(
      sample.int(length(ids), length(ids), replace = TRUE),
      nbins = length(ids)
    )
    draws[bb, ] <- stat_fun(sample_dt$weight * freq[pid])
  }
  list(point = point, draws = draws)
}

# ---------------------------------------------------------------------------
# BLOCK 1: family log-income process, reference persons 25-60, 1984-2019.
# ---------------------------------------------------------------------------
message("BLOCK 1: income process...")

inc <- dt[
  reference_person &
    is.finite(age) & age >= 25 & age <= 60 &
    is.finite(income) & income > 1000
]
inc[, age_i := as.integer(round(age))]
inc <- inc[age_i >= 25 & age_i <= 60]
inc[, y := log(income)]

n_dup <- nrow(inc) - uniqueN(inc[, .(id, year)])
if (n_dup > 0L) {
  data_problems <- c(data_problems, sprintf(
    "Block 1: %d duplicate reference-person (id, year) rows dropped (kept first).",
    n_dup
  ))
  setorder(inc, id, year)
  inc <- inc[, .SD[1L], by = .(id, year)]
}
inc[, rid := .I]

LAGS <- c(1L, 2L, 4L, 6L, 8L)

# Pair row-index tables, built once. Both legs of a pair belong to the same
# person, so one person frequency rescales the whole pair in the bootstrap.
res_key <- inc[, .(rid, id, year)]
pair_tables <- lapply(LAGS, function(k) {
  p <- merge(
    res_key[, .(id, jy = year, r1 = rid)],
    res_key[, .(id, jy = year - k, r2 = rid)],
    by = c("id", "jy")
  )
  data.table(
    r1 = p$r1, r2 = p$r2,
    pw0 = sqrt(inc$weight[p$r1] * inc$weight[p$r2]),
    pid_pair = p$id,
    annual_era = (p$jy + k) <= 1997L
  )
})
names(pair_tables) <- paste0("lag", LAGS)

income_moments <- function(row_weights) {
  # Residualize log income on age dummies and year dummies, weighted, then
  # compute the weighted residual variance and pooled pair autocovariances.
  keep <- which(row_weights > 0)
  fit <- feols(
    y ~ 1 | age_i + year,
    data = inc[keep],
    weights = row_weights[keep],
    notes = FALSE
  )
  u <- rep(NA_real_, nrow(inc))
  u[keep] <- as.numeric(resid(fit))
  v0 <- weighted_var_safe(u, row_weights)
  n0 <- sum(row_weights > 0)
  covs <- numeric(length(LAGS))
  npairs <- numeric(length(LAGS))
  for (jj in seq_along(LAGS)) {
    pt <- pair_tables[[jj]]
    fpair <- row_weights[pt$r1] / inc$weight[pt$r1] # person draw frequency
    pw <- pt$pw0 * fpair
    covs[jj] <- weighted_cov_pairs(u[pt$r1], u[pt$r2], pw)
    npairs[jj] <- sum(fpair)
  }
  list(
    moments = c(v0, covs),
    lags = c(0L, LAGS),
    nobs = c(n0, npairs)
  )
}

fit_income_process <- function(moments, lags, wts, restrict_sigma_e = FALSE) {
  # u = p + e, p_t = rho p_{t-1} + eta:
  #   cov(k) = (sigma_eta^2 / (1 - rho^2)) rho^k  for k >= 1
  #   var(0) = sigma_eta^2 / (1 - rho^2) + sigma_e^2
  # Weighted nonlinear least squares on the moment vector, weights = N pairs.
  ok <- is.finite(moments) & is.finite(wts) & wts > 0
  moments <- moments[ok]
  lags <- lags[ok]
  wts <- wts[ok] / sum(wts[ok])
  obj <- function(par) {
    rho <- 1 / (1 + exp(-par[1]))
    s_eta <- exp(par[2])
    s_e <- if (restrict_sigma_e) 0 else exp(par[3])
    sp2 <- s_eta^2 / (1 - rho^2)
    pred <- ifelse(lags == 0L, sp2 + s_e^2, sp2 * rho^lags)
    sum(wts * (moments - pred)^2)
  }
  best <- NULL
  v0 <- moments[lags == 0L][1]
  for (rho0 in c(0.85, 0.92, 0.96, 0.99)) {
    kpos <- lags > 0L
    sp2_0 <- mean(moments[kpos] / rho0^lags[kpos])
    sp2_0 <- max(sp2_0, 1e-6)
    s_eta0 <- sqrt(sp2_0 * (1 - rho0^2))
    s_e0 <- sqrt(max(v0 - sp2_0, 1e-6))
    par0 <- c(log(rho0 / (1 - rho0)), log(s_eta0))
    if (!restrict_sigma_e) par0 <- c(par0, log(s_e0))
    cand <- optim(par0, obj, method = "Nelder-Mead",
                  control = list(maxit = 20000, reltol = 1e-14))
    cand <- optim(cand$par, obj, method = "BFGS",
                  control = list(maxit = 2000, reltol = 1e-14))
    if (is.null(best) || cand$value < best$value) best <- cand
  }
  rho <- 1 / (1 + exp(-best$par[1]))
  s_eta <- exp(best$par[2])
  s_e <- if (restrict_sigma_e) 0 else exp(best$par[3])
  c(
    rho_annual = rho,
    sigma_eta_annual = s_eta,
    sigma_e = s_e,
    rho_4yr = rho^4,
    stationary_var_persistent = s_eta^2 / (1 - rho^2),
    wnls_objective = best$value
  )
}

base_mom <- income_moments(inc$weight)
fit_full <- fit_income_process(base_mom$moments, base_mom$lags, base_mom$nobs,
                               restrict_sigma_e = FALSE)
fit_ar1 <- fit_income_process(base_mom$moments, base_mom$lags, base_mom$nobs,
                              restrict_sigma_e = TRUE)

# Person-cluster bootstrap: rerun residualization + moment construction + both
# fits inside every replication. Time the first replications and fall back to
# 199 reps if 499 would be too slow.
b1_ids <- unique(inc$id)
b1_pid <- match(inc$id, b1_ids)
b1_stat <- function(freq) {
  mm <- income_moments(inc$weight * freq[b1_pid])
  c(
    full = fit_income_process(mm$moments, mm$lags, mm$nobs, FALSE)[1:5],
    ar1 = fit_income_process(mm$moments, mm$lags, mm$nobs, TRUE)[1:5],
    empirical = setNames(mm$moments, paste0("mom_lag", mm$lags))
  )
}

b1_reps <- BOOT_REPS_DEFAULT
set.seed(BOOT_SEED)
t0 <- Sys.time()
first_draws <- list()
for (bb in 1:5) {
  freq <- tabulate(sample.int(length(b1_ids), length(b1_ids), replace = TRUE),
                   nbins = length(b1_ids))
  first_draws[[bb]] <- b1_stat(freq)
}
per_rep <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5
projected_minutes <- per_rep * BOOT_REPS_DEFAULT / 60
message(sprintf("Block 1 bootstrap: %.2f s/rep, projected %.1f min for %d reps.",
                per_rep, projected_minutes, BOOT_REPS_DEFAULT))
if (projected_minutes > 45) {
  b1_reps <- 199L
  data_problems <- c(data_problems, sprintf(
    "Block 1 bootstrap reduced to 199 reps (projected %.0f min for 499 at %.2f s/rep).",
    projected_minutes, per_rep
  ))
}
b1_draws <- matrix(NA_real_, nrow = b1_reps, ncol = length(first_draws[[1]]),
                   dimnames = list(NULL, names(first_draws[[1]])))
for (bb in 1:5) b1_draws[bb, ] <- first_draws[[bb]]
for (bb in 6:b1_reps) {
  freq <- tabulate(sample.int(length(b1_ids), length(b1_ids), replace = TRUE),
                   nbins = length(b1_ids))
  b1_draws[bb, ] <- b1_stat(freq)
  if (bb %% 50L == 0L) message("  Block 1 bootstrap rep ", bb, " / ", b1_reps)
}
b1_se <- apply(b1_draws, 2L, sd, na.rm = TRUE)

pred_cov <- function(fit, lag) {
  sp2 <- fit[["stationary_var_persistent"]]
  ifelse(lag == 0L, sp2 + fit[["sigma_e"]]^2, sp2 * fit[["rho_annual"]]^lag)
}

autocov_table <- data.table(
  lag_years = base_mom$lags,
  n_pairs = base_mom$nobs,
  n_pairs_annual_era = c(
    nrow(inc),
    vapply(pair_tables, function(pt) sum(pt$annual_era), numeric(1))
  ),
  n_pairs_biennial_era = c(
    0L,
    vapply(pair_tables, function(pt) sum(!pt$annual_era), numeric(1))
  ),
  n_persons = c(
    uniqueN(inc$id),
    vapply(pair_tables, function(pt) uniqueN(pt$pid_pair), numeric(1))
  ),
  empirical_moment = base_mom$moments,
  bootstrap_se = b1_se[paste0("empirical.mom_lag", base_mom$lags)],
  fitted_persistent_transitory = pred_cov(fit_full, base_mom$lags),
  fitted_ar1_only = pred_cov(fit_ar1, base_mom$lags)
)
fwrite(autocov_table, file.path(outdir, "block1_income_process_autocovariances.csv"))

b1_params <- rbindlist(list(
  data.table(
    specification = "persistent_plus_transitory",
    parameter = names(fit_full)[1:5],
    estimate = as.numeric(fit_full[1:5]),
    bootstrap_se = as.numeric(b1_se[paste0("full.", names(fit_full)[1:5])])
  ),
  data.table(
    specification = "ar1_only_sigma_e_0",
    parameter = names(fit_ar1)[1:5],
    estimate = as.numeric(fit_ar1[1:5]),
    bootstrap_se = as.numeric(b1_se[paste0("ar1.", names(fit_ar1)[1:5])])
  ),
  data.table(
    specification = "current_model_values",
    parameter = c("rho_annual", "sigma_eta_annual", "sigma_e",
                  "rho_4yr", "stationary_var_persistent"),
    estimate = c(MODEL_RHO_ANNUAL, MODEL_SIGMA_ETA_ANNUAL, NA_real_,
                 MODEL_RHO_ANNUAL^4, MODEL_STATIONARY_VAR),
    bootstrap_se = NA_real_
  )
))
b1_params[, `:=`(
  n_person_years = nrow(inc),
  n_persons = uniqueN(inc$id),
  bootstrap_reps = fifelse(specification == "current_model_values", NA_integer_, b1_reps),
  wnls_objective = fcase(
    specification == "persistent_plus_transitory", fit_full[["wnls_objective"]],
    specification == "ar1_only_sigma_e_0", fit_ar1[["wnls_objective"]],
    default = NA_real_
  )
)]
fwrite(b1_params, file.path(outdir, "block1_income_process_estimates.csv"))

# ---------------------------------------------------------------------------
# BLOCK 2: entry nonhousing wealth / income for reference persons 18-24.
# Construction mirrors PSID_ENTRY_WEALTH_RATIO_NODES_2535: childless via
# RELCHIREP == 0, renter via HOMEOWN == 2, INCFAMR > 1000, ratio =
# NETWORTH2R / INCFAMR, IW weights, weighted quintile-bin means.
# ---------------------------------------------------------------------------
message("BLOCK 2: entry wealth 18-24...")

quintile_bin_table <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  bin <- cut(cw, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1 + 1e-12),
             labels = FALSE, include.lowest = TRUE)
  bt <- data.table(bin = bin, x = x, w = w)
  bt[, .(
    bin_mean_ratio = sum(w * x) / sum(w),
    bin_weight = sum(w) / sum(bt$w)
  ), by = bin][order(bin)]
}

make_entry_sample <- function(source, age_lo, age_hi, reference_only = TRUE) {
  out <- source[
    is.finite(age) & age >= age_lo & age <= age_hi &
      n_children_rep == 0 & !is.na(own) & own == 0 &
      is.finite(income) & income > 1000 & is.finite(nonhousing_nw)
  ]
  if (reference_only) out <- out[reference_person == TRUE]
  out[, ratio := nonhousing_nw / income]
  out[is.finite(ratio)]
}

entry_block_rows <- function(sample_label, es, reps) {
  boot <- person_bootstrap(
    es,
    function(w) c(mean_ratio = weighted_mean_safe(es$ratio, w)),
    reps
  )
  bins <- quintile_bin_table(es$ratio, es$weight)
  rows <- data.table(
    sample = sample_label,
    statistic = c(
      paste0("quintile_bin_mean_ratio_", bins$bin),
      paste0("quintile_bin_weight_", bins$bin),
      "overall_weighted_mean_ratio",
      "overall_weighted_median_ratio"
    ),
    value = c(
      bins$bin_mean_ratio, bins$bin_weight,
      weighted_mean_safe(es$ratio, es$weight),
      weighted_quantile_safe(es$ratio, es$weight, 0.5)
    ),
    bootstrap_se = NA_real_,
    n_family_years = nrow(es),
    n_persons = uniqueN(es$id),
    bootstrap_reps = NA_integer_
  )
  rows[statistic == "overall_weighted_mean_ratio",
       `:=`(bootstrap_se = sd(boot$draws[, "mean_ratio"], na.rm = TRUE),
            bootstrap_reps = reps)]
  rows
}

entry_1824 <- make_entry_sample(dt, 18, 24)
entry_1826 <- make_entry_sample(dt, 18, 26)
entry_2535_replica <- make_entry_sample(dt, 25, 35, reference_only = FALSE)

message(sprintf(
  "Entry samples: 18-24 %d family-years / %d persons; 18-26 %d / %d; 25-35 replication %d / %d.",
  nrow(entry_1824), uniqueN(entry_1824$id),
  nrow(entry_1826), uniqueN(entry_1826$id),
  nrow(entry_2535_replica), uniqueN(entry_2535_replica$id)
))
usable_1824 <- nrow(entry_1824) >= 500L
if (!usable_1824) {
  data_problems <- c(data_problems, sprintf(
    "Block 2: only %d family-years at ages 18-24 (< 500); the 18-26 version is the usable object.",
    nrow(entry_1824)
  ))
}

block2 <- rbindlist(list(
  entry_block_rows("reference_childless_renter_18_24", entry_1824, BOOT_REPS_DEFAULT),
  entry_block_rows("reference_childless_renter_18_26", entry_1826, BOOT_REPS_DEFAULT),
  entry_block_rows("replication_2535_original_filters_all_persons",
                   entry_2535_replica, BOOT_REPS_DEFAULT)
))
block2[, usable_flag := fcase(
  sample == "reference_childless_renter_18_24" & usable_1824, "primary",
  sample == "reference_childless_renter_18_24" & !usable_1824, "underpowered_n_below_500",
  sample == "reference_childless_renter_18_26" & usable_1824, "robustness",
  sample == "reference_childless_renter_18_26" & !usable_1824, "primary_fallback",
  default = "replication_check"
)]
stored_check <- data.table(
  sample = "stored_2535_calibration_constants",
  statistic = c(
    paste0("quintile_bin_mean_ratio_", 1:5),
    paste0("quintile_bin_weight_", 1:5),
    "overall_weighted_mean_ratio",
    "overall_weighted_median_ratio"
  ),
  value = c(STORED_ENTRY_NODES_2535, STORED_ENTRY_WEIGHTS_2535,
            STORED_ENTRY_MEAN_2535, STORED_ENTRY_MEDIAN_2535),
  bootstrap_se = NA_real_,
  n_family_years = NA_integer_,
  n_persons = NA_integer_,
  bootstrap_reps = NA_integer_,
  usable_flag = "stored_constant"
)
block2 <- rbindlist(list(block2, stored_check), use.names = TRUE)
fwrite(block2, file.path(outdir, "block2_entry_wealth_18_24.csv"))

# ---------------------------------------------------------------------------
# BLOCK 3: composition share, reference persons 65-75, ratio >= threshold.
# Same filters as the estate targets: finite NETWORTH2R / INCFAMR ratio and
# INCFAMR > 1000.
# ---------------------------------------------------------------------------
message("BLOCK 3: composition share 65-75...")

comp <- dt[
  reference_person & is.finite(age) & age >= 65 & age <= 75 &
    is.finite(income) & income > 1000 & is.finite(nonhousing_nw)
]
comp[, ratio := nonhousing_nw / income]
comp <- comp[is.finite(ratio)]

comp_stat <- function(w) c(
  share_ratio_ge_0p5 = weighted_mean_safe(comp$ratio >= 0.5, w),
  share_ratio_ge_1p0 = weighted_mean_safe(comp$ratio >= 1.0, w),
  share_ratio_ge_2p0 = weighted_mean_safe(comp$ratio >= 2.0, w)
)
comp_boot <- person_bootstrap(comp, comp_stat, BOOT_REPS_DEFAULT)
comp_se <- apply(comp_boot$draws, 2L, sd, na.rm = TRUE)
block3 <- data.table(
  moment = names(comp_boot$point),
  threshold = c(0.5, 1.0, 2.0),
  role = c("robustness", "headline", "robustness"),
  value = as.numeric(comp_boot$point),
  bootstrap_se = as.numeric(comp_se),
  implied_weight_inv_se2 = 1 / as.numeric(comp_se)^2,
  n_family_years = nrow(comp),
  n_persons = uniqueN(comp$id),
  bootstrap_reps = BOOT_REPS_DEFAULT
)
fwrite(block3, file.path(outdir, "block3_composition_share_6575.csv"))

# ---------------------------------------------------------------------------
# BLOCK 4: old-age homeownership, reference persons 65-75, HOMEOWN tenure.
# ---------------------------------------------------------------------------
message("BLOCK 4: old-age ownership 65-75...")

own_old <- dt[
  reference_person & is.finite(age) & age >= 65 & age <= 75 & !is.na(own)
]
own_stat <- function(w) c(own_rate_6575 = weighted_mean_safe(own_old$own, w))
own_boot <- person_bootstrap(own_old, own_stat, BOOT_REPS_DEFAULT)
own_se <- sd(own_boot$draws[, "own_rate_6575"], na.rm = TRUE)
block4 <- data.table(
  moment = c("psid_reference_own_rate_6575", "legacy_stored_old_age_own_rate"),
  value = c(as.numeric(own_boot$point), LEGACY_OLD_AGE_OWN_RATE),
  bootstrap_se = c(own_se, NA_real_),
  implied_weight_inv_se2 = c(1 / own_se^2, NA_real_),
  gap_psid_minus_legacy = c(as.numeric(own_boot$point) - LEGACY_OLD_AGE_OWN_RATE, NA_real_),
  n_family_years = c(nrow(own_old), NA_integer_),
  n_persons = c(uniqueN(own_old$id), NA_integer_),
  bootstrap_reps = c(BOOT_REPS_DEFAULT, NA_integer_)
)
fwrite(block4, file.path(outdir, "block4_oldage_ownership_6575.csv"))

# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
if (length(data_problems) == 0L) data_problems <- "None."
readme <- c(
  "# Intergen income-process, entry-wealth, composition, and ownership targets",
  "",
  "Built by `estimate_intergen_income_entry_targets.R` from",
  "`PSID/PSIDSHELF_MOBILITY.dta`. Shared conventions: reference persons",
  "(`RELTOHEAD_ == 10`), age `AGEREP`, weight `IW > 0`, income `INCFAMR`,",
  "nonhousing wealth `NETWORTH2R`, waves 1984-2019, alive in the observation",
  "year (`year <= DEATHYEAR` or missing), person-cluster bootstrap (resample",
  "unique person ids, reweight family-year rows by draw frequency, seed",
  "20260715).",
  "",
  "## Block 1 (`block1_income_process_*.csv`)",
  "",
  "Reference persons aged 25-60 with `INCFAMR > 1000`, all waves. `y = log(INCFAMR)`",
  "is residualized on age dummies (25-60) and year dummies in an IW-weighted",
  "fixed-effects regression. Autocovariances pool every valid same-person pair at",
  "each calendar-year lag k in {1, 2, 4, 6, 8}; lag-1 pairs exist only in the",
  "annual era (second year <= 1997); other lags pool annual-era and biennial-era",
  "pairs, with the era split reported. Pair weight = sqrt(IW_t * IW_{t+k}).",
  "The persistent+transitory model u = p + e, p_t = rho p_{t-1} + eta implies",
  "cov(k) = (sigma_eta^2/(1-rho^2)) rho^k for k >= 1 and",
  "var(0) = sigma_eta^2/(1-rho^2) + sigma_e^2; (rho, sigma_eta, sigma_e) is",
  "estimated by nonlinear least squares on the moment vector, weighting each",
  "moment by its pair count. The `ar1_only_sigma_e_0` rows restrict sigma_e = 0.",
  "Bootstrap SEs re-run residualization, moments, and both fits inside every",
  "replication. `current_model_values` is the comparison row for the live model",
  "parameters (rho, sigma_eta annual; stationary persistent variance 0.05337).",
  "",
  "## Block 2 (`block2_entry_wealth_18_24.csv`)",
  "",
  "Replicates the construction behind `PSID_ENTRY_WEALTH_RATIO_NODES_2535`",
  "(calibration.py): childless via `RELCHIREP == 0`, renter via `HOMEOWN == 2`",
  "(`HOMEOWN == 3` excluded), `INCFAMR > 1000`, ratio `NETWORTH2R/INCFAMR`,",
  "IW-weighted quintile-bin means (bins cut at cumulative weight 0.2, 0.4, 0.6,",
  "0.8). New samples restrict to reference persons aged 18-24 (primary if at",
  "least 500 family-years) and 18-26 (fallback/robustness). NETWORTH2R exists",
  "only in wealth-supplement waves, which limits N. The",
  "`replication_2535_original_filters_all_persons` rows rebuild the original",
  "25-35 object (which did NOT restrict to reference persons) as a validation",
  "check against the stored constants (`stored_2535_calibration_constants`).",
  "",
  "## Block 3 (`block3_composition_share_6575.csv`)",
  "",
  "IW-weighted share of reference-person family-years aged 65-75 (finite",
  "`NETWORTH2R/INCFAMR`, `INCFAMR > 1000`) with ratio >= 1.0 (headline);",
  "thresholds 0.5 and 2.0 are robustness rows. Implied SMM weight = 1/SE^2.",
  "",
  "## Block 4 (`block4_oldage_ownership_6575.csv`)",
  "",
  "IW-weighted homeownership rate (`HOMEOWN == 1` vs `== 2`) of reference",
  "persons aged 65-75, compared with the legacy stored `old_age_own_rate`",
  "0.76426097 (which the moment-SE audit attributes to an ACS/MMS builder; the",
  "known PSID check value there is 0.863730454).",
  "",
  "## Data problems",
  "",
  paste0("- ", data_problems),
  "",
  "## Files",
  "",
  "- `block1_income_process_estimates.csv`: fitted income-process parameters + model comparison row.",
  "- `block1_income_process_autocovariances.csv`: empirical vs fitted autocovariance table with pair counts.",
  "- `block2_entry_wealth_18_24.csv`: entry wealth quintile nodes, weights, means, medians, replication check.",
  "- `block3_composition_share_6575.csv`: old-age wealth/income composition shares.",
  "- `block4_oldage_ownership_6575.csv`: old-age ownership rate vs legacy constant."
)
writeLines(readme, file.path(outdir, "README.md"))

message("Wrote income/entry target packet to ", outdir)
