#!/usr/bin/env Rscript

# Reproducibility receipt for the two ACS room targets in the active E5 profile.
#
# This intentionally writes to a separate receipt directory.  It does not edit
# the target constants, calibration weights, or the older 14-moment SE packet.

suppressPackageStartupMessages({
  library(data.table)
})

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  hit <- file_arg[startsWith(file_arg, "--file=")]
  if (length(hit) != 1L) stop("Run with Rscript so --file= is available.")
  dirname(normalizePath(sub("^--file=", "", hit)))
}

parse_args <- function() {
  ans <- list(B = 1000L, seed = 20260705L, output_dir = NULL)
  for (arg in commandArgs(trailingOnly = TRUE)) {
    if (startsWith(arg, "--B=")) ans$B <- as.integer(sub("^--B=", "", arg))
    if (startsWith(arg, "--seed=")) ans$seed <- as.integer(sub("^--seed=", "", arg))
    if (startsWith(arg, "--output-dir=")) {
      ans$output_dir <- sub("^--output-dir=", "", arg)
    }
  }
  if (!is.finite(ans$B) || ans$B < 2L) stop("--B must be an integer of at least 2.")
  if (!is.finite(ans$seed)) stop("--seed must be a finite integer.")
  ans
}

sha256_file <- function(path) {
  openssl <- Sys.which("openssl")
  if (!nzchar(openssl)) stop("The openssl executable is required for SHA-256 verification.")
  output <- system2(
    openssl, c("dgst", "-sha256", shQuote(path)),
    stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(output, "status")) && attr(output, "status") != 0L) {
    stop("openssl SHA-256 failed for ", path, ": ", paste(output, collapse = " "))
  }
  hash <- trimws(sub("^.*= *", "", tail(output, 1L)))
  if (!grepl("^[0-9a-f]{64}$", hash)) stop("Could not parse SHA-256 for ", path)
  hash
}

assert_close <- function(actual, expected, tolerance, label) {
  gap <- abs(actual - expected)
  if (!is.finite(gap) || gap > tolerance) {
    stop(sprintf(
      "%s failed: actual=%.17g expected=%.17g abs_diff=%.17g tolerance=%.17g",
      label, actual, expected, gap, tolerance
    ))
  }
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) stop("Weighted mean has no valid positive-weight observations.")
  sum(x[ok] * w[ok]) / sum(w[ok])
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
args <- parse_args()
if (is.null(args$output_dir)) {
  args$output_dir <- file.path(script_dir, "output_active_acs_room_target_receipt_20260817")
} else if (!grepl("^/", args$output_dir)) {
  args$output_dir <- file.path(repo_root, args$output_dir)
}
out_dir <- normalizePath(args$output_dir, mustWork = FALSE)
if (file.exists(out_dir)) {
  stop("Output directory already exists; choose a new --output-dir to preserve the prior receipt: ", out_dir)
}

paths <- c(
  raw_extract = file.path(repo_root, "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"),
  lookup_2010 = file.path(repo_root, "code/data/mms_center_periphery/data/puma_mms_lookup_2010.csv"),
  lookup_2020 = file.path(repo_root, "code/data/mms_center_periphery/data/puma_mms_lookup_2020.csv"),
  lookup_build_summary = file.path(repo_root, "code/data/mms_center_periphery/data/build_summary.txt"),
  point_builder = file.path(repo_root, "code/data/mms_center_periphery/build_intergen_one_market_housing_targets.R"),
  bootstrap_builder = file.path(repo_root, "code/data/moment_standard_errors/build_moment_bootstrap_se.R"),
  analysis_cache = file.path(repo_root, "code/data/moment_standard_errors/cache/acs_analysis_samples.rds"),
  receipt_builder = normalizePath(file.path(script_dir, "build_active_acs_room_target_receipt.R"))
)
missing <- names(paths)[!file.exists(paths)]
if (length(missing) > 0L) stop("Missing required receipt inputs: ", paste(missing, collapse = ", "))

# These fingerprints freeze the exact empirical inputs and builders audited on
# 2026-08-17.  A changed input must be reviewed and deliberately re-receipted.
expected_sha256 <- c(
  raw_extract = "edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e",
  lookup_2010 = "7d2f6927ce9d5f17eb9d1fac7bfcedc0d2702f294ceaffb1a280e25dd178fbd4",
  lookup_2020 = "464236a42dc2691a08718be8c449e4f3b1399c9017985f20fd2ca825cbe79a5d",
  lookup_build_summary = "4b6cee19a5fa1427234e1fde1dea62413549cd1cd31237cb38d5ec052378523e",
  point_builder = "a08ef1a07b7db004b66596d3d126c6515975e6a7eee713e6cde7b423058ad099",
  bootstrap_builder = "7731b89c94dcdc2c1dab2408e10bd209e8ee2992d155cbefb782c4fa940c9248",
  analysis_cache = "0eae9aaed4e1d3b9655235be967378c1d389f29ec9b64f005a811c2f0ced7df0"
)

message("Hashing pinned ACS receipt inputs (the 9.2 GB raw extract is the slow step)...")
observed_sha256 <- vapply(paths, sha256_file, character(1))
for (nm in names(expected_sha256)) {
  if (!identical(observed_sha256[[nm]], expected_sha256[[nm]])) {
    stop(sprintf(
      "SHA-256 mismatch for %s: observed=%s expected=%s",
      nm, observed_sha256[[nm]], expected_sha256[[nm]]
    ))
  }
}

message("Loading the pinned ACS household-head analysis cache...")
cache <- readRDS(paths[["analysis_cache"]])
if (!is.list(cache) || !all(c("targets", "women") %in% names(cache))) {
  stop("Pinned cache does not have the expected targets/women structure.")
}
dt <- as.data.table(cache$targets)
required_columns <- c(
  "year", "met2013", "age", "pernum", "hhwt", "ownershp", "rooms",
  "nchild", "yngch", "parent_u18", "owner", "renter", "mms_location"
)
missing_columns <- setdiff(required_columns, names(dt))
if (length(missing_columns) > 0L) {
  stop("Pinned cache is missing columns: ", paste(missing_columns, collapse = ", "))
}
if (nrow(dt) != 4103889L) stop("Unexpected pinned-cache target row count: ", nrow(dt))
if (min(dt$year) != 2012L || max(dt$year) != 2023L) stop("Unexpected year range in ACS cache.")
if (uniqueN(dt$met2013) != 42L) stop("Unexpected number of observed MMS metros in ACS cache.")

# The cache builder has already imposed year 2012--2023, identifiable metro,
# household head, positive HHWT, owner/renter, positive literal ROOMS, exact
# period-specific MMS PUMA match, and middle-MMS -> center.  The active mean
# rooms moment additionally imposes the stated upper age bound of 85.
mean_sample <- dt[age >= 18L & age <= 85L]
gap_sample <- mean_sample[
  age >= 30L & age <= 55L & parent_u18 == TRUE & nchild >= 1
]
gap_3plus <- gap_sample[nchild >= 3]
gap_1to2 <- gap_sample[nchild >= 1 & nchild <= 2]

target_values <- c(
  prime30_55_parent_3plus_minus_1to2_mean_rooms = 0.36769955881,
  aggregate_mean_occupied_rooms_18_85 = 5.779970481941968
)
point_values <- c(
  prime30_55_parent_3plus_minus_1to2_mean_rooms =
    weighted_mean(gap_3plus$rooms, gap_3plus$hhwt) -
    weighted_mean(gap_1to2$rooms, gap_1to2$hhwt),
  aggregate_mean_occupied_rooms_18_85 = weighted_mean(mean_sample$rooms, mean_sample$hhwt)
)
point_tolerance <- 1e-10
for (nm in names(target_values)) {
  assert_close(point_values[[nm]], target_values[[nm]], point_tolerance, paste("target reproduction", nm))
}

expected_counts <- c(mean_n = 3986589, gap_n = 918595, gap_3plus_n = 235597, gap_1to2_n = 682998)
observed_counts <- c(
  mean_n = nrow(mean_sample), gap_n = nrow(gap_sample),
  gap_3plus_n = nrow(gap_3plus), gap_1to2_n = nrow(gap_1to2)
)
if (!identical(as.numeric(observed_counts), as.numeric(expected_counts))) {
  stop("ACS receipt sample counts changed: ", paste(observed_counts, collapse = ", "))
}
expected_weights <- c(mean_weight = 419883313, gap_weight = 105193582, gap_3plus_weight = 28441511, gap_1to2_weight = 76752071)
observed_weights <- c(
  mean_weight = sum(mean_sample$hhwt), gap_weight = sum(gap_sample$hhwt),
  gap_3plus_weight = sum(gap_3plus$hhwt), gap_1to2_weight = sum(gap_1to2$hhwt)
)
for (nm in names(expected_weights)) {
  assert_close(observed_weights[[nm]], expected_weights[[nm]], 1e-8, paste("sample weight", nm))
}

# Collapse sufficient statistics by metro, then resample the 42 metros with
# replacement.  The same draw is used for both targets, preserving covariance.
metros <- sort(unique(mean_sample$met2013))
components <- data.table(met2013 = metros)
add_component <- function(base, sample, prefix) {
  cell <- sample[, .(
    numerator = sum(hhwt * rooms),
    denominator = sum(hhwt)
  ), by = met2013]
  setnames(cell, c("numerator", "denominator"), paste0(prefix, c("_num", "_den")))
  out <- merge(base, cell, by = "met2013", all.x = TRUE, sort = FALSE)
  for (column in paste0(prefix, c("_num", "_den"))) out[is.na(get(column)), (column) := 0]
  out
}
components <- add_component(components, mean_sample, "all")
components <- add_component(components, gap_3plus, "g3")
components <- add_component(components, gap_1to2, "g12")

draw_ratio <- function(num_name, den_name, multiplicity) {
  denominator <- sum(components[[den_name]] * multiplicity)
  if (!is.finite(denominator) || denominator <= 0) stop("Bootstrap draw has an empty denominator.")
  sum(components[[num_name]] * multiplicity) / denominator
}

set.seed(args$seed)
boot <- matrix(
  NA_real_, nrow = args$B, ncol = 2L,
  dimnames = list(NULL, names(target_values))
)
for (b in seq_len(args$B)) {
  multiplicity <- tabulate(
    sample.int(nrow(components), nrow(components), replace = TRUE),
    nbins = nrow(components)
  )
  boot[b, "prime30_55_parent_3plus_minus_1to2_mean_rooms"] <-
    draw_ratio("g3_num", "g3_den", multiplicity) -
    draw_ratio("g12_num", "g12_den", multiplicity)
  boot[b, "aggregate_mean_occupied_rooms_18_85"] <-
    draw_ratio("all_num", "all_den", multiplicity)
}
if (any(!is.finite(boot))) stop("Bootstrap contains nonfinite draws.")
boot_cov <- cov(boot)
boot_se <- sqrt(diag(boot_cov))
if (any(!is.finite(boot_se) | boot_se <= 0)) stop("Bootstrap produced a nonpositive or nonfinite SE.")
if (max(abs(diag(boot_cov) - boot_se^2)) > 1e-14) stop("Covariance diagonal does not equal SE squared.")
boot_cor <- cov2cor(boot_cov)
if (any(!is.finite(boot_cor)) || any(abs(boot_cor) > 1 + 1e-12)) stop("Invalid bootstrap correlation matrix.")

target_receipt <- data.table(
  moment_key = names(target_values),
  target = as.numeric(target_values),
  reproduced_point = as.numeric(point_values[names(target_values)]),
  absolute_difference = abs(point_values[names(target_values)] - target_values),
  reproduction_tolerance = point_tolerance,
  reproduction_pass = TRUE,
  n_unweighted = c(nrow(gap_sample), nrow(mean_sample)),
  hhwt_sum = c(sum(gap_sample$hhwt), sum(mean_sample$hhwt)),
  n_metro_clusters = length(metros),
  bootstrap_B = args$B,
  bootstrap_seed = args$seed,
  bootstrap_mean = colMeans(boot)[names(target_values)],
  metro_bootstrap_se = boot_se[names(target_values)],
  declared_synthetic_se = 0.05 * abs(as.numeric(target_values)),
  active_calibration_weight = 1 / (0.05 * abs(as.numeric(target_values)))^2,
  metro_bootstrap_used_in_active_calibration = FALSE,
  formula = c(
    "HHWT-mean ROOMS for NCHILD>=3 minus HHWT-mean ROOMS for NCHILD in 1:2",
    "sum(HHWT*ROOMS)/sum(HHWT)"
  )
)
setnames(target_receipt, "moment_key", "key")

covariance <- as.data.table(boot_cov, keep.rownames = "key")
correlation <- as.data.table(boot_cor, keep.rownames = "key")
file_info <- file.info(paths)
provenance <- data.table(
  artifact = names(paths),
  path = unname(paths),
  sha256 = unname(observed_sha256[names(paths)]),
  bytes = as.numeric(file_info$size),
  modified_at = format(file_info$mtime, "%Y-%m-%dT%H:%M:%S%z")
)

readme <- c(
  "# Active ACS Room-Target Receipt",
  "",
  sprintf("Generated: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Bootstrap: `B=%d`, seed `%d`, 42 `MET2013` clusters resampled with replacement.", args$B, args$seed),
  "",
  "## Result",
  "",
  sprintf(
    "- `prime30_55_parent_3plus_minus_1to2_mean_rooms`: target `%.14f`, reproduced `%.15f`, metro-bootstrap SE `%.9f`; active declared synthetic SE `%.9f` (weight `%.6f`).",
    target_values[[1]], point_values[[1]], boot_se[[1]],
    0.05 * target_values[[1]], 1 / (0.05 * target_values[[1]])^2
  ),
  sprintf(
    "- `aggregate_mean_occupied_rooms_18_85`: target `%.15f`, reproduced `%.15f`, metro-bootstrap SE `%.9f`; active declared synthetic SE `%.9f` (weight `%.6f`).",
    target_values[[2]], point_values[[2]], boot_se[[2]],
    0.05 * target_values[[2]], 1 / (0.05 * target_values[[2]])^2
  ),
  sprintf("- Same-draw covariance: `%.9f`; correlation: `%.9f`.", boot_cov[1, 2], boot_cor[1, 2]),
  "",
  "These metro-bootstrap SEs are a project-design uncertainty measure, not an official ACS replicate-weight variance estimate. The extract does not carry ACS replicate weights. They are reported for audit and are not substituted for the active calibration's declared synthetic five-percent SEs.",
  "",
  "## Exact empirical objects",
  "",
  "The source is the pooled 2012--2023 IPUMS ACS `extract27.dta`. The base sample keeps identifiable `MET2013` observations in households (`GQ` 1 or 2), the first person in the household (`PERNUM=1`), positive `HHWT`, owner or renter tenure, and positive literal `ROOMS`. It joins 2012--2021 observations to the 2010-PUMA lookup and 2022--2023 observations to the 2020-PUMA lookup on state, PUMA, and CBSA. Matched center and periphery PUMAs are retained and the lookup's middle category is collapsed to center.",
  "",
  "Mean rooms additionally keeps ages 18--85. The family-size contrast keeps ages 30--55 with at least one own child in the household and a youngest own child younger than 18; it compares `NCHILD>=3` with `NCHILD` 1--2. `NCHILD` is current own children in the household, not completed parity. Both estimators use `HHWT`.",
  "",
  "The pinned geography build summary records a 0.30 core-tract population target, 0.50 center-PUMA cutoff, 0.10 periphery-PUMA cutoff, 51 retained cities, and 1,361 PUMA-city rows in each period lookup. The realized target sample contains 42 distinct `MET2013` codes.",
  "",
  "## Why the older moment_se.csv is NA",
  "",
  "The committed `code/data/moment_standard_errors/output/moment_se.csv` was last regenerated with the PSID-only source branch. Its README and timing log contain no ACS load or bootstrap entry, and its reproduction log says `ACS source not run`. Thus the NA is a run-selection artifact, not a failed ACS point-estimate gate. The older 14-moment harness also does not contain the active aggregate mean-rooms row.",
  "",
  "## Provenance limitation",
  "",
  "This receipt independently pins the raw extract, lookup files, relevant builders, and the analysis cache, and reproduces the target formulas from the pinned cache. The legacy cache predates a source-hash manifest, so its relationship to the raw extract is reconstructed rather than cryptographically recorded at cache creation. A future raw rebuild should write these hashes when the cache is created; changing any pinned input makes this script fail before emitting a receipt.",
  "",
  "## Files",
  "",
  "- `target_receipt.csv`: point estimates, samples, weights, and bootstrap SEs.",
  "- `bootstrap_covariance.csv` and `bootstrap_correlation.csv`: same-draw two-moment matrices.",
  "- `provenance.csv`: SHA-256, byte size, and modification time for every input and builder."
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(target_receipt, file.path(out_dir, "target_receipt.csv"))
fwrite(covariance, file.path(out_dir, "bootstrap_covariance.csv"))
fwrite(correlation, file.path(out_dir, "bootstrap_correlation.csv"))
fwrite(provenance, file.path(out_dir, "provenance.csv"))
writeLines(readme, file.path(out_dir, "README.md"))

message("Receipt complete: ", out_dir)
print(target_receipt[, .(key, target, reproduced_point, absolute_difference, metro_bootstrap_se)])
print(boot_cov)
