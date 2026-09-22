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
  ok <- isTRUE(all.equal(fast$summary_lower, legacy$summary_lower, tolerance = 1e-6)) &&
    isTRUE(all.equal(fast$summary_upper, legacy$summary_upper, tolerance = 1e-6)) &&
    fast$n_accepted == legacy$n_accepted && fast$n_components == legacy$n_components
  cat(sprintf("[%s] fast=[%s,%s] n=%d legacy=[%s,%s] n=%d -> %s\n",
              label, fast$summary_lower, fast$summary_upper, fast$n_accepted,
              legacy$summary_lower, legacy$summary_upper, legacy$n_accepted,
              if (ok) "MATCH" else "MISMATCH"))
  ok
}

results <- c(
  check_case("strong_FS", 1200, 0.5, seed = 1),
  check_case("weak_FS", 1200, 0.02, seed = 2),
  check_case("negative_FS", 1200, -0.4, seed = 3)
)
if (all(results)) cat("AR_FAST_EQUIVALENCE_ALL_MATCH\n") else { cat("AR_FAST_EQUIVALENCE_MISMATCH_DETECTED\n"); quit(status = 1) }
