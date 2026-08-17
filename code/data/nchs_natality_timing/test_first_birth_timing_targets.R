#!/usr/bin/env Rscript

# Deterministic checks for the NCHS timing-target measurement contract.

options(stringsAsFactors = FALSE)

script_file <- function() {
  arg <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", arg, value = TRUE)
  if (length(hit) == 1L) return(normalizePath(sub("^--file=", "", hit)))
  normalizePath("test_first_birth_timing_targets.R")
}

root <- dirname(script_file())
counts <- read.csv(file.path(root, "first_birth_counts_year_age.csv"), stringsAsFactors = FALSE)
manifest <- read.csv(file.path(root, "first_birth_counts_manifest.csv"), stringsAsFactors = FALSE)
exclusions <- read.csv(file.path(root, "order_exclusion_shares_by_year.csv"), stringsAsFactors = FALSE)
legacy <- read.csv(file.path(root, "timing_targets.csv"), stringsAsFactors = FALSE)
exact <- read.csv(file.path(root, "timing_targets_exact_age.csv"), stringsAsFactors = FALSE)
model <- read.csv(file.path(root, "timing_targets_model_comparable.csv"), stringsAsFactors = FALSE)
bins <- read.csv(file.path(root, "timing_age_bin_counts.csv"), stringsAsFactors = FALSE)
contract <- read.csv(file.path(root, "timing_target_contract.csv"), stringsAsFactors = FALSE)
metadata <- jsonlite::read_json(file.path(root, "timing_target_metadata.json"), simplifyVector = FALSE)

assert_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
assert_close <- function(actual, expected, tolerance = 5e-13, message = "values differ") {
  if (length(actual) != length(expected) || any(!is.finite(actual)) ||
      any(abs(actual - expected) > tolerance)) {
    stop(sprintf("%s: actual=%s expected=%s", message,
                 paste(actual, collapse = ","), paste(expected, collapse = ",")),
         call. = FALSE)
  }
}

assert_true(identical(sort(unique(counts$year)), 1987L:2023L), "NCHS data-vintage years changed")
assert_true(nrow(counts) == 37L * 38L, "year-by-age cache is not a complete 1987--2023 by 12--49 grid")
assert_true(identical(manifest$year, 1987L:2023L), "raw-file manifest years changed")
assert_true(identical(names(manifest), c("year", "file", "size_bytes", "sha256", "cache_contract")),
            "raw-file manifest no longer pins content hashes")
assert_true(all(grepl("^[0-9a-f]{64}$", manifest$sha256)), "raw-file SHA-256 receipt is malformed")
assert_true(all(manifest$cache_contract == "nchs_first_birth_counts_v4_lbo_unknown_code_fixed_raw_sha256"),
            "raw-file manifest has the wrong cache contract")
assert_true(exclusions$n_unknown_live_birth_order[exclusions$year == 2003L] > 0,
            "2003 LBO_REC code 9 is not recorded as unknown")
assert_true(all(exclusions$order_variable[exclusions$year >= 2003L] == "lbo_rec"),
            "modern live-birth-order field receipt changed")
assert_close(counts$n_first_births[counts$year == 2003L & counts$age == 25L], 80773,
             tolerance = 0, message = "2003 MAGER41 code 12 was not decoded to age 25")
assert_close(counts$n_first_births[counts$year == 2003L & counts$age == 12L], 0,
             tolerance = 0, message = "2003 recode was again treated as literal age")
assert_true(identical(names(legacy), c(
  "window", "cohort_range", "mean_age_first_birth", "share_30plus",
  "n_first_births", "truncation_estimated_share"
)), "legacy exact-age schema changed")
assert_true(identical(legacy, exact[, names(legacy)]), "legacy and explicit exact-age outputs differ")

counts$cohort <- counts$year - counts$age
windows <- list(primary_1979_1984 = 1979:1984, sensitivity_1975_1980 = 1975:1980)
for (window_name in names(windows)) {
  x <- subset(counts, cohort %in% windows[[window_name]])
  boundary_start <- ifelse(
    x$age < 22L, 18L,
    ifelse(x$age >= 42L, 42L, 18L + 4L * floor((x$age - 18L) / 4L))
  )
  common_keep <- x$age >= 18L & x$age <= 45L
  common_start <- 18L + 4L * floor((x$age[common_keep] - 18L) / 4L)

  boundary_row <- subset(model, window == window_name & timing_operator == "boundary_collapsed")
  common_row <- subset(model, window == window_name & timing_operator == "common_support_18_45")
  assert_true(nrow(boundary_row) == 1L && nrow(common_row) == 1L, "operator/window row is not unique")

  boundary_mean_start <- weighted.mean(boundary_start, x$n_first_births)
  boundary_mean_midpoint <- weighted.mean(boundary_start + 2, x$n_first_births)
  boundary_share_30 <- sum(x$n_first_births[boundary_start >= 30L]) / sum(x$n_first_births)
  common_mean_start <- weighted.mean(common_start, x$n_first_births[common_keep])
  common_mean_midpoint <- weighted.mean(common_start + 2, x$n_first_births[common_keep])
  common_share_30 <- sum(x$n_first_births[common_keep][common_start >= 30L]) /
    sum(x$n_first_births[common_keep])

  assert_close(boundary_row$mean_age_first_birth_period_start_labels, boundary_mean_start,
               message = "boundary period-start mean mismatch")
  assert_close(boundary_row$mean_age_first_birth, boundary_mean_midpoint,
               message = "boundary midpoint mean mismatch")
  assert_close(boundary_row$share_30plus, boundary_share_30,
               message = "boundary age-30 share mismatch")
  assert_close(common_row$mean_age_first_birth_period_start_labels, common_mean_start,
               message = "common-support period-start mean mismatch")
  assert_close(common_row$mean_age_first_birth, common_mean_midpoint,
               message = "common-support midpoint mean mismatch")
  assert_close(common_row$share_30plus, common_share_30,
               message = "common-support age-30 share mismatch")

  # Thirty is exactly a cell boundary under both operators.
  assert_close(
    sum(x$n_first_births[boundary_start >= 30L]),
    sum(x$n_first_births[x$age >= 30L]),
    tolerance = 0,
    message = "boundary collapse changed age-30 classification"
  )
  assert_close(
    sum(x$n_first_births[common_keep][common_start >= 30L]),
    sum(x$n_first_births[common_keep & x$age >= 30L]),
    tolerance = 0,
    message = "common-support operator changed age-30 classification"
  )

  assert_close(boundary_row$n_first_births_used, sum(x$n_first_births), tolerance = 0,
               message = "boundary collapse did not retain all observed births")
  assert_close(common_row$n_first_births_used, sum(x$n_first_births[common_keep]), tolerance = 0,
               message = "common-support used count mismatch")
  assert_close(common_row$n_first_births_excluded,
               sum(x$n_first_births[!common_keep]), tolerance = 0,
               message = "common-support excluded count mismatch")

  for (operator in c("boundary_collapsed", "common_support_18_45")) {
    b <- subset(bins, window == window_name & timing_operator == operator)
    assert_true(identical(b$bin_start, seq(18L, 42L, by = 4L)), "bin-start grid changed")
    assert_true(identical(b$bin_midpoint, seq(20L, 44L, by = 4L)), "bin-midpoint grid changed")
    expected_n <- if (operator == "boundary_collapsed") sum(x$n_first_births) else sum(x$n_first_births[common_keep])
    assert_close(sum(b$n_first_births), expected_n, tolerance = 0, message = "bin counts do not exhaust used births")
    assert_close(sum(b$share_first_births), 1, tolerance = 5e-15, message = "bin shares do not sum to one")
  }
}

primary_boundary <- subset(model, window == "primary_1979_1984" & timing_operator == "boundary_collapsed")
sensitivity_boundary <- subset(model, window == "sensitivity_1975_1980" & timing_operator == "boundary_collapsed")
primary_common <- subset(model, window == "primary_1979_1984" & timing_operator == "common_support_18_45")
sensitivity_common <- subset(model, window == "sensitivity_1975_1980" & timing_operator == "common_support_18_45")

assert_close(primary_boundary$mean_age_first_birth, 26.0446272574833, message = "primary boundary mean changed")
assert_close(primary_boundary$share_30plus, 0.260327401666964, message = "primary boundary share changed")
assert_close(sensitivity_boundary$mean_age_first_birth, 25.9190402686846, message = "sensitivity boundary mean changed")
assert_close(sensitivity_boundary$share_30plus, 0.250397460163499, message = "sensitivity boundary share changed")
assert_close(primary_common$mean_age_first_birth, 26.6650697126504, message = "primary common-support mean changed")
assert_close(primary_common$share_30plus, 0.287048349933470, message = "primary common-support share changed")
assert_close(sensitivity_common$mean_age_first_birth, 26.6050674798020, message = "sensitivity common-support mean changed")
assert_close(sensitivity_common$share_30plus, 0.279424130172943, message = "sensitivity common-support share changed")

assert_close(primary_boundary$primary_sensitivity_spread_mean_age, 0.125586988798773,
             message = "boundary mean window spread changed")
assert_close(primary_boundary$primary_sensitivity_spread_share_30plus, 0.00992994150346543,
             message = "boundary share window spread changed")
assert_close(primary_boundary$declared_se_mean_age, 0.15, message = "boundary mean SE changed")
assert_close(primary_boundary$declared_se_share_30plus, 0.01, message = "boundary share SE changed")
assert_close(primary_boundary$inverse_variance_weight_mean_age, 1 / 0.15^2, message = "boundary mean weight changed")
assert_close(primary_boundary$inverse_variance_weight_share_30plus, 1 / 0.01^2,
             tolerance = 1e-9,
             message = "boundary share weight changed")
assert_close(primary_common$declared_se_mean_age, 0.1, message = "common-support mean SE changed")
assert_close(primary_common$declared_se_share_30plus, 0.008, message = "common-support share SE changed")

# Midpoint versus period-start labeling must move data and any model moment by
# the same constant, leaving their gap unchanged.
illustrative_model_start_mean <- 24.4104
gap_start <- illustrative_model_start_mean - primary_boundary$mean_age_first_birth_period_start_labels
gap_midpoint <- (illustrative_model_start_mean + 2) - primary_boundary$mean_age_first_birth
assert_close(gap_start, gap_midpoint, message = "midpoint relabeling changed the model-data gap")

recommended_contract <- subset(contract, target_status == "recommended_model_comparable" & window == "primary_1979_1984")
assert_true(nrow(recommended_contract) == 2L, "recommended primary activation contract is incomplete")
assert_close(
  recommended_contract$target_value[recommended_contract$moment == "mean_age_first_birth"],
  primary_boundary$mean_age_first_birth,
  message = "contract mean does not match receipt"
)
assert_close(
  recommended_contract$target_value[recommended_contract$moment == "share_30plus"],
  primary_boundary$share_30plus,
  message = "contract share does not match receipt"
)

assert_true(metadata$contract_id == "nchs_first_birth_timing_v4_boundary_collapsed_midpoint_1987_2023_lbo_unknown_fixed",
            "metadata contract id changed")
assert_true(metadata$recommended_operator$id == "boundary_collapsed_midpoint_v1",
            "metadata operator id changed")
assert_close(unlist(metadata$recommended_operator$reported_midpoint_labels), seq(20L, 44L, by = 4L),
             tolerance = 0, message = "metadata midpoint labels changed")
assert_true(isFALSE(metadata$contains_five_dated_period_targets),
            "cohort contract incorrectly claims to contain five dated targets")
assert_true(metadata$maternal_age_source_2003$official_document_sha256 ==
              "82246197e30d54c56a69314bfdcb8f553e6ca0a0d509f8ce8112c2e996b5b2f5",
            "official 2003 NCHS documentation hash changed")
assert_true(metadata$pre_2003_measurement$official_1987_document_sha256 ==
              "145d8496d764ee66583ef1520ae56804e1c3373bea35d3a570ecbc1426bd5c79",
            "official 1987 NCHS documentation hash changed")
assert_true(metadata$pre_2003_measurement$official_1989_document_sha256 ==
              "92dab8115baec71eec3633239cbd042b2079ad6b80bd1b3a3a43c3276ac3a7cb",
            "official 1989 NCHS documentation hash changed")
assert_true(metadata$live_birth_order_codes$legacy_dlivord_unknown == 99L &&
              metadata$live_birth_order_codes$lbo_rec_unknown == 9L,
            "vintage-specific unknown live-birth-order codes changed")

for (relative_path in names(metadata$file_sha256)) {
  actual_hash <- digest::digest(file.path(root, relative_path), algo = "sha256", file = TRUE)
  assert_true(actual_hash == metadata$file_sha256[[relative_path]],
              paste("metadata file hash mismatch:", relative_path))
}
source_members <- unlist(metadata$source_bundle_members)
source_bundle_hash <- digest::digest(
  paste(names(source_members), source_members, sep = "=", collapse = "\n"),
  algo = "sha256", serialize = FALSE
)
assert_true(source_bundle_hash == metadata$source_bundle_sha256, "source-bundle hash mismatch")
contract_members <- unlist(metadata$contract_bundle_members)
contract_bundle_hash <- digest::digest(
  paste(names(contract_members), contract_members, sep = "=", collapse = "\n"),
  algo = "sha256", serialize = FALSE
)
assert_true(contract_bundle_hash == metadata$contract_bundle_sha256, "contract-bundle hash mismatch")

cat("NCHS_TIMING_TARGET_TESTS_PASS tests=28\n")
