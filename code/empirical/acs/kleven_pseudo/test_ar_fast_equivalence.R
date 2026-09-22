# Equivalence check: fast analytic AR (ar_confidence_set) vs the
# per-grid-point legacy implementation (ar_confidence_set_legacy), across
# strong/weak/negative first-stage, unequal weights, and household
# clustering. Not exhaustive; bounded local check.
.args <- commandArgs(trailingOnly = FALSE)
.sd <- dirname(sub("^--file=", "", .args[grepl("^--file=", .args)]))
source(file.path(.sd, "twins_samesex_iv_lib.R"))
suppressMessages(library(data.table))

check_case <- function(label, n, fs_coef, sd_y = 1.2, seed = 1, grid = seq(-2, 2, by = 0.2)) {
  set.seed(seed)
  mat_age <- round(runif(n, 25, 40))
  hh <- paste0("hh", seq_len(n) %% (n / 3))  # multi-row clusters
  w <- sample(1:5, n, replace = TRUE)
  Z <- rbinom(n, 1, 0.2)
  D <- rbinom(n, 1, pmin(pmax(0.3 + fs_coef * Z + 0.01 * (mat_age - 30), 0.01), 0.99))
  Y <- 5 + 0.6 * D + 0.02 * (mat_age - 30) + rnorm(n, sd = sd_y)
  d <- data.table(Y = Y, D = D, Z = Z, mat_age = mat_age, household_key = hh, w = w)
  fast <- ar_confidence_set(d, "Y", "D", "Z", "mat_age", "w", "household_key", grid)
  legacy <- ar_confidence_set_legacy(d, "Y", "D", "Z", "mat_age", "w", "household_key", grid)
  bounds_ok <- isTRUE(all.equal(fast$summary_lower, legacy$summary_lower, tolerance = 1e-6)) &&
    isTRUE(all.equal(fast$summary_upper, legacy$summary_upper, tolerance = 1e-6)) &&
    fast$n_accepted == legacy$n_accepted && fast$n_components == legacy$n_components
  # Numerical equivalence proper: per-grid-point stat and p-value, not just
  # aggregated bounds/counts.
  pergrid_ok <- isTRUE(all.equal(fast$pvalues, legacy$pvalues, tolerance = 1e-6)) &&
    isTRUE(all.equal(fast$stats, legacy$stats, tolerance = 1e-6)) &&
    identical(fast$status_vec, legacy$status_vec)
  ok <- bounds_ok && pergrid_ok
  cat(sprintf("[%s] bounds=[%s,%s]/[%s,%s] pergrid_max_abs_p_diff=%.2e pergrid_max_abs_stat_diff=%.2e status_identical=%s -> %s\n",
              label, fast$summary_lower, fast$summary_upper, legacy$summary_lower, legacy$summary_upper,
              max(abs(fast$pvalues - legacy$pvalues), na.rm = TRUE),
              max(abs(fast$stats - legacy$stats), na.rm = TRUE),
              identical(fast$status_vec, legacy$status_vec),
              if (ok) "MATCH" else "MISMATCH"))
  ok
}

results <- c(
  check_case("strong_FS", 1200, 0.5, seed = 1),
  check_case("weak_FS", 1200, 0.02, seed = 2),
  check_case("negative_FS", 1200, -0.4, seed = 3),
  check_case("unequal_weights_multirow_HH", 900, 0.35, seed = 4)
)

## Missing-control fixture: some rows have NA mat_age, so the complete-case
## sample must match across fast's 3 base regressions and legacy's per-grid
## regressions. This exercises the identical-sample guard directly.
set.seed(5)
n <- 800
mat_age <- round(runif(n, 25, 40)); mat_age[sample(n, 40)] <- NA
hh <- paste0("hh", seq_len(n) %% 300)
w <- sample(1:5, n, replace = TRUE)
Z <- rbinom(n, 1, 0.25)
D <- rbinom(n, 1, pmin(pmax(0.3 + 0.4 * Z, 0.01), 0.99))
Y <- 5 + 0.5 * D + rnorm(n)
d_na <- data.table(Y = Y, D = D, Z = Z, mat_age = mat_age, household_key = hh, w = w)
fast_na <- ar_confidence_set(d_na, "Y", "D", "Z", "mat_age", "w", "household_key", seq(-1, 1, 0.25))
legacy_na <- ar_confidence_set_legacy(d_na, "Y", "D", "Z", "mat_age", "w", "household_key", seq(-1, 1, 0.25))
na_ok <- isTRUE(all.equal(fast_na$pvalues, legacy_na$pvalues, tolerance = 1e-6))
cat(sprintf("[missing_controls] identical_sample_guard_triggered=%s pergrid_match=%s -> %s\n",
            isTRUE(fast_na$base_regression_failure), na_ok, if (na_ok) "MATCH" else "MISMATCH"))
results <- c(results, na_ok)
if (all(results)) cat("AR_FAST_EQUIVALENCE_ALL_MATCH\n") else { cat("AR_FAST_EQUIVALENCE_MISMATCH_DETECTED\n"); quit(status = 1) }
