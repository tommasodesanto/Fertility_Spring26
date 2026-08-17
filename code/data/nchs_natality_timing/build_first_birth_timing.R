#!/usr/bin/env Rscript

# Build cohort first-birth timing targets from the local NCHS natality archive.
# The only microdata read is the age/order pair needed for this construction.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(haven))

script_file <- function() {
  arg <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", arg, value = TRUE)
  if (length(hit) == 1L) return(normalizePath(sub("^--file=", "", hit)))
  normalizePath("build_first_birth_timing.R")
}

out_dir <- dirname(script_file())
archive_root <- Sys.getenv(
  "NCHS_NATALITY_ROOT",
  "/Users/tommasodesanto/Desktop/Projects/Datasets/natality_data"
)
force_rebuild <- "--force" %in% commandArgs(trailingOnly = TRUE)
year_min <- 1987L
year_max <- 2023L
ages <- 12:49
model_age_min <- 18L
model_age_max <- 45L
model_bin_width <- 4L
model_bin_starts <- seq(model_age_min, 42L, by = model_bin_width)
cache_contract <- "nchs_first_birth_counts_v4_lbo_unknown_code_fixed_raw_sha256"

cache_path <- file.path(out_dir, "first_birth_counts_year_age.csv")
manifest_path <- file.path(out_dir, "first_birth_counts_manifest.csv")
exclusions_path <- file.path(out_dir, "order_exclusion_shares_by_year.csv")
targets_path <- file.path(out_dir, "timing_targets.csv")
exact_targets_path <- file.path(out_dir, "timing_targets_exact_age.csv")
model_targets_path <- file.path(out_dir, "timing_targets_model_comparable.csv")
contract_path <- file.path(out_dir, "timing_target_contract.csv")
bin_counts_path <- file.path(out_dir, "timing_age_bin_counts.csv")
metadata_path <- file.path(out_dir, "timing_target_metadata.json")
truncation_path <- file.path(out_dir, "cohort_truncation.csv")
readme_path <- file.path(out_dir, "README.md")

find_us_files <- function(root) {
  candidates <- list.files(
    root, recursive = TRUE, full.names = TRUE,
    pattern = "^natality([0-9]{4}us|us[0-9]{4})\\.dta$"
  )
  year <- as.integer(sub(
    "^natality(?:us)?([0-9]{4})(?:us)?\\.dta$", "\\1",
    basename(candidates), perl = TRUE
  ))
  # Freeze the empirical vintage. A newly downloaded file must not silently
  # change a calibration target.
  keep <- !is.na(year) & year >= year_min & year <= year_max
  ans <- data.frame(year = year[keep], path = normalizePath(candidates[keep]),
                    stringsAsFactors = FALSE)
  ans <- ans[order(ans$year, ans$path), , drop = FALSE]
  if (anyDuplicated(ans$year)) stop("More than one US natality file found for a year.")
  if (!nrow(ans)) stop("No US natality files found for 1987--2023 under: ", root)
  missing_years <- setdiff(year_min:year_max, ans$year)
  if (length(missing_years)) stop("Missing required US natality years: ", paste(missing_years, collapse = ", "))
  ans
}

make_manifest <- function(files) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package digest is required to fingerprint the raw natality files.")
  }
  data.frame(
    year = files$year,
    file = files$path,
    # Store byte sizes as text so a CSV round trip cannot change their type or
    # precision (and thereby spuriously invalidate the cache).
    size_bytes = format(file.info(files$path)$size, scientific = FALSE, trim = TRUE),
    sha256 = vapply(
      files$path,
      function(path) digest::digest(path, algo = "sha256", file = TRUE),
      character(1)
    ),
    cache_contract = cache_contract,
    stringsAsFactors = FALSE
  )
}

same_manifest <- function(current) {
  if (!file.exists(manifest_path) || !file.exists(cache_path) || !file.exists(exclusions_path)) return(FALSE)
  old <- tryCatch(
    read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(old) || !identical(names(old), names(current))) return(FALSE)
  old$year <- as.integer(old$year)
  old$size_bytes <- as.character(old$size_bytes)
  old$sha256 <- as.character(old$sha256)
  identical(old, current)
}

variables_for_file <- function(path, year) {
  header <- read_dta(path, n_max = 0)
  available <- names(header)
  if (year <= 2002L) {
    age <- "dmage"
    order <- "dlivord"
  } else if (year == 2003L) {
    age <- "mager41"
    order <- "lbo_rec"
  } else {
    age <- "mager"
    order <- "lbo_rec"
  }
  if (!(age %in% available) || !(order %in% available)) {
    stop("Missing contract-required ", age, "/", order, " in ", path, " (year ", year, ").")
  }
  c(age = age, order = order)
}

decode_maternal_age <- function(raw_age, age_variable) {
  age <- as.numeric(raw_age)
  if (age_variable != "mager41") return(age)
  # Official 2003 NCHS detail-file documentation defines MAGER41 code 01 as
  # under 15, code 02 as age 15, ..., code 37 as age 50. Code 01 is assigned
  # 14 solely to place it below the model's first cutoff; the target cohorts
  # studied here are age 19--28 in 2003, so this grouped tail cannot affect
  # their exact-age or model-bin moments.
  decoded <- rep(NA_real_, length(age))
  decoded[age == 1] <- 14
  decoded[age >= 2 & age <= 41] <- age[age >= 2 & age <= 41] + 13
  decoded
}

build_cache <- function(files, manifest) {
  counts <- list()
  exclusions <- list()
  for (i in seq_len(nrow(files))) {
    year <- files$year[i]
    path <- files$path[i]
    vars <- variables_for_file(path, year)
    message(sprintf("Reading %d (%s, %s) ...", year, vars["age"], vars["order"]))
    # `col_select` means haven reads only these two columns, not the full record.
    d <- read_dta(path, col_select = tidyselect::all_of(unname(vars)))
    age <- decode_maternal_age(d[[vars["age"]]], vars["age"])
    order <- as.numeric(d[[vars["order"]]])
    valid_age <- !is.na(age) & age %in% ages
    # The legacy detail field uses 99 for unknown live-birth order, whereas
    # LBO_REC uses 9.  The codes are vintage-specific even though both fields
    # use 1 for a first birth.
    unknown_order <- is.na(order) | if (vars["order"] == "lbo_rec") order == 9 else order == 99
    first <- valid_age & !unknown_order & order == 1
    by_age <- tabulate(match(age[first], ages), nbins = length(ages))
    counts[[i]] <- data.frame(year = year, age = ages, n_first_births = by_age)
    exclusions[[i]] <- data.frame(
      year = year,
      n_birth_records_age_12_49 = sum(valid_age),
      n_unknown_live_birth_order = sum(valid_age & unknown_order),
      unknown_live_birth_order_share = mean(unknown_order[valid_age]),
      age_variable = vars["age"],
      age_transformation = ifelse(
        vars["age"] == "mager41",
        "NCHS 41-category recode decoded: 01->14 proxy for under-15; 02:41->code+13",
        "literal single-year age"
      ),
      order_variable = vars["order"],
      stringsAsFactors = FALSE
    )
    rm(d, age, order)
    gc(verbose = FALSE)
  }
  write.csv(do.call(rbind, counts), cache_path, row.names = FALSE)
  write.csv(do.call(rbind, exclusions), exclusions_path, row.names = FALSE)
  write.csv(manifest, manifest_path, row.names = FALSE)
}

weighted_mean_or_na <- function(x, w) if (sum(w) == 0) NA_real_ else sum(x * w) / sum(w)

window_definitions <- function() data.frame(
  window = c("primary_1979_1984", "sensitivity_1975_1980"),
  cohort_min = c(1979L, 1975L), cohort_max = c(1984L, 1980L)
)

common_support_bin <- function(age) {
  ans <- rep(NA_integer_, length(age))
  keep <- age >= model_age_min & age <= model_age_max
  ans[keep] <- model_age_min + model_bin_width * floor((age[keep] - model_age_min) / model_bin_width)
  ans
}

boundary_collapsed_bin <- function(age) {
  # The model has no pre-18 parity state and no post-45 fertility decision.
  # Preserve all observed first births by assigning the two tails to its first
  # and last cells. Interior cutoffs exactly match the model's four-year cells.
  ans <- common_support_bin(pmin(pmax(age, model_age_min), model_age_max))
  ans[age < 22L] <- 18L
  ans[age >= 42L] <- 42L
  ans
}

round_up_to <- function(x, step) {
  # The tolerance prevents a value already on the grid from rounding up only
  # because of binary floating-point representation.
  step * ceiling(x / step - 1e-12)
}

window_tail_estimate <- function(x, reference, operator) {
  cohorts <- sort(unique(x$cohort))
  rows <- lapply(cohorts, function(cohort_value) {
    z <- subset(x, cohort == cohort_value)
    max_age <- if (any(z$n_first_births > 0)) max(z$age[z$n_first_births > 0]) else NA_integer_
    if (operator == "common_support_18_45") {
      ref_denominator <- sum(reference$n_first_births[
        reference$age >= model_age_min & reference$age <= model_age_max
      ])
      ref_missed <- if (is.na(max_age)) NA_real_ else sum(reference$n_first_births[
        reference$age > max_age & reference$age <= model_age_max &
          reference$age >= model_age_min
      ])
      observed_weight <- sum(z$n_first_births[
        z$age >= model_age_min & z$age <= model_age_max
      ])
    } else {
      ref_denominator <- sum(reference$n_first_births)
      ref_missed <- if (is.na(max_age)) NA_real_ else sum(reference$n_first_births[
        reference$age > max_age
      ])
      observed_weight <- sum(z$n_first_births)
    }
    data.frame(
      cohort = cohort_value,
      max_observed_age = max_age,
      observed_weight = observed_weight,
      reference_denominator = ref_denominator,
      reference_missed = ref_missed,
      estimated_share = ref_missed / ref_denominator
    )
  })
  rows <- do.call(rbind, rows)
  weighted_mean_or_na(rows$estimated_share, rows$observed_weight)
}

measure_model_operator <- function(x, reference, operator) {
  if (operator == "boundary_collapsed") {
    used <- x
    used$bin_start <- boundary_collapsed_bin(used$age)
    support_description <- "all observed ages 12--49; exact ages <22 collapse to 18 and exact ages >=42 collapse to 42"
    tail_operator <- "boundary_collapsed"
  } else if (operator == "common_support_18_45") {
    used <- subset(x, age >= model_age_min & age <= model_age_max)
    used$bin_start <- common_support_bin(used$age)
    support_description <- "exact ages 18--45 only; exact ages below 18 and above 45 excluded"
    tail_operator <- "common_support_18_45"
  } else {
    stop("Unknown timing operator: ", operator)
  }
  if (anyNA(used$bin_start) || !all(used$bin_start %in% model_bin_starts)) {
    stop("Timing operator did not map every used age to a model bin.")
  }
  observed_total <- sum(x$n_first_births)
  used_total <- sum(used$n_first_births)
  below <- sum(x$n_first_births[x$age < model_age_min])
  above <- sum(x$n_first_births[x$age > model_age_max])
  share_30_bin <- sum(used$n_first_births[used$bin_start >= 30L]) / used_total
  share_30_exact <- sum(used$n_first_births[used$age >= 30L]) / used_total
  if (!isTRUE(all.equal(share_30_bin, share_30_exact, tolerance = 1e-15))) {
    stop("The age-30 threshold is not invariant under the timing operator.")
  }
  list(
    mean_age_first_birth = weighted_mean_or_na(used$bin_start, used$n_first_births),
    # Exact integer age a represents the continuous interval [a,a+1), so the
    # four-year cells [18,22),... have calendar-age midpoints 20,24,...,44.
    mean_age_first_birth_midpoint_labels = weighted_mean_or_na(used$bin_start + 2, used$n_first_births),
    share_30plus = share_30_bin,
    n_first_births_used = used_total,
    n_first_births_observed_all_ages = observed_total,
    n_first_births_below_18 = below,
    n_first_births_above_45 = above,
    n_first_births_excluded = observed_total - used_total,
    observed_excluded_share = (observed_total - used_total) / observed_total,
    support_description = support_description,
    estimated_unobserved_right_tail_share = window_tail_estimate(x, reference, tail_operator),
    used = used
  )
}

make_targets <- function(counts) {
  cohort_counts <- transform(counts, cohort = year - age)
  # The integer assignment is conventional: births in year t to an age-a mother
  # are assigned to t-a. Dates within a calendar year create an unavoidable +/-1 blur.
  reference <- subset(cohort_counts, cohort >= 1970 & cohort <= 1974)
  windows <- window_definitions()
  exact_targets <- list()
  model_targets <- list()
  bin_counts <- list()
  truncation <- list()
  for (i in seq_len(nrow(windows))) {
    w <- windows[i, ]
    x <- subset(cohort_counts, cohort >= w$cohort_min & cohort <= w$cohort_max)
    cohort_rows <- lapply(w$cohort_min:w$cohort_max, function(c) {
      z <- subset(x, cohort == c)
      max_age <- if (any(z$n_first_births > 0)) max(z$age[z$n_first_births > 0]) else NA_integer_
      above <- if (is.na(max_age)) NA_real_ else sum(reference$n_first_births[reference$age > max_age])
      ref_total <- sum(reference$n_first_births)
      data.frame(
        window = w$window, cohort = c, max_observed_age = max_age,
        n_observed_first_births = sum(z$n_first_births),
        reference_n_first_births = ref_total,
        reference_n_first_births_above_max_age = above,
        estimated_share_plausibly_missed = above / ref_total
      )
    })
    tr <- do.call(rbind, cohort_rows)
    truncation[[i]] <- tr
    exact_targets[[i]] <- data.frame(
      window = w$window,
      cohort_range = sprintf("%d-%d", w$cohort_min, w$cohort_max),
      mean_age_first_birth = weighted_mean_or_na(x$age, x$n_first_births),
      share_30plus = if (sum(x$n_first_births) == 0) NA_real_ else sum(x$n_first_births[x$age >= 30]) / sum(x$n_first_births),
      n_first_births = sum(x$n_first_births),
      truncation_estimated_share = weighted_mean_or_na(
        tr$estimated_share_plausibly_missed, tr$n_observed_first_births
      )
    )

    for (operator in c("boundary_collapsed", "common_support_18_45")) {
      measured <- measure_model_operator(x, reference, operator)
      key <- paste(w$window, operator, sep = "__")
      model_targets[[key]] <- data.frame(
        window = w$window,
        cohort_range = sprintf("%d-%d", w$cohort_min, w$cohort_max),
        timing_operator = operator,
        bin_label_convention = "midpoint",
        support_definition = measured$support_description,
        mean_age_first_birth = measured$mean_age_first_birth_midpoint_labels,
        mean_age_first_birth_period_start_labels = measured$mean_age_first_birth,
        mean_age_first_birth_midpoint_labels = measured$mean_age_first_birth_midpoint_labels,
        share_30plus = measured$share_30plus,
        n_first_births_used = measured$n_first_births_used,
        n_first_births_observed_all_ages = measured$n_first_births_observed_all_ages,
        n_first_births_below_18 = measured$n_first_births_below_18,
        n_first_births_above_45 = measured$n_first_births_above_45,
        n_first_births_excluded = measured$n_first_births_excluded,
        observed_excluded_share = measured$observed_excluded_share,
        estimated_unobserved_right_tail_share = measured$estimated_unobserved_right_tail_share,
        stringsAsFactors = FALSE
      )
      counts_by_bin <- aggregate(n_first_births ~ bin_start, measured$used, sum)
      counts_by_bin <- merge(
        data.frame(bin_start = model_bin_starts), counts_by_bin,
        by = "bin_start", all.x = TRUE, sort = TRUE
      )
      counts_by_bin$n_first_births[is.na(counts_by_bin$n_first_births)] <- 0
      counts_by_bin$share_first_births <- counts_by_bin$n_first_births / measured$n_first_births_used
      counts_by_bin$window <- w$window
      counts_by_bin$cohort_range <- sprintf("%d-%d", w$cohort_min, w$cohort_max)
      counts_by_bin$timing_operator <- operator
      counts_by_bin$bin_end <- counts_by_bin$bin_start + 3L
      counts_by_bin$bin_midpoint <- counts_by_bin$bin_start + 2
      counts_by_bin$bin_label_convention <- "midpoint"
      counts_by_bin$exact_age_rule <- ifelse(
        counts_by_bin$bin_start == 18L,
        ifelse(operator == "boundary_collapsed", "age <= 21", "18 <= age <= 21"),
        ifelse(
          counts_by_bin$bin_start == 42L,
          ifelse(operator == "boundary_collapsed", "age >= 42", "42 <= age <= 45"),
          sprintf("%d <= age <= %d", counts_by_bin$bin_start, counts_by_bin$bin_end)
        )
      )
      bin_counts[[key]] <- counts_by_bin[, c(
        "window", "cohort_range", "timing_operator", "bin_start", "bin_end", "bin_midpoint",
        "bin_label_convention", "exact_age_rule", "n_first_births", "share_first_births"
      )]
    }
  }

  exact <- do.call(rbind, exact_targets)
  model <- do.call(rbind, model_targets)
  for (operator in unique(model$timing_operator)) {
    idx <- model$timing_operator == operator
    primary <- model[idx & model$window == "primary_1979_1984", ]
    sensitivity <- model[idx & model$window == "sensitivity_1975_1980", ]
    mean_spread <- abs(primary$mean_age_first_birth - sensitivity$mean_age_first_birth)
    share_spread <- abs(primary$share_30plus - sensitivity$share_30plus)
    model$primary_sensitivity_spread_mean_age[idx] <- mean_spread
    model$primary_sensitivity_spread_share_30plus[idx] <- share_spread
    # Preserve the pre-existing conservative rounding convention: mean-age
    # uncertainty rounds up to 0.05 year and the share uncertainty to 0.001.
    model$declared_se_mean_age[idx] <- round_up_to(mean_spread, 0.05)
    model$declared_se_share_30plus[idx] <- round_up_to(share_spread, 0.001)
    model$inverse_variance_weight_mean_age[idx] <- 1 / model$declared_se_mean_age[idx]^2
    model$inverse_variance_weight_share_30plus[idx] <- 1 / model$declared_se_share_30plus[idx]^2
  }

  exact_mean_spread <- abs(diff(exact$mean_age_first_birth))
  exact_share_spread <- abs(diff(exact$share_30plus))
  exact$primary_sensitivity_spread_mean_age <- exact_mean_spread
  exact$primary_sensitivity_spread_share_30plus <- exact_share_spread
  exact$declared_se_mean_age <- round_up_to(exact_mean_spread, 0.05)
  exact$declared_se_share_30plus <- round_up_to(exact_share_spread, 0.001)
  exact$inverse_variance_weight_mean_age <- 1 / exact$declared_se_mean_age^2
  exact$inverse_variance_weight_share_30plus <- 1 / exact$declared_se_share_30plus^2

  contract_rows <- list()
  append_contract <- function(source, operator, status, support, label) {
    rows <- list()
    for (i in seq_len(nrow(source))) {
      for (moment in c("mean_age_first_birth", "share_30plus")) {
        is_mean <- moment == "mean_age_first_birth"
        rows[[length(rows) + 1L]] <- data.frame(
          target_status = status,
          timing_operator = operator,
          window = source$window[i],
          cohort_range = source$cohort_range[i],
          moment = moment,
          target_value = source[[moment]][i],
          primary_sensitivity_window_spread = if (is_mean) source$primary_sensitivity_spread_mean_age[i] else source$primary_sensitivity_spread_share_30plus[i],
          declared_standard_error = if (is_mean) source$declared_se_mean_age[i] else source$declared_se_share_30plus[i],
          inverse_variance_weight = if (is_mean) source$inverse_variance_weight_mean_age[i] else source$inverse_variance_weight_share_30plus[i],
          support_definition = support[i],
          age_label_convention = label,
          age_30plus_rule = if (moment != "share_30plus") {
            "not applicable"
          } else if (operator == "exact_age_observed") {
            "exact age >= 30"
          } else {
            "bin_start >= 30; identical to exact age >= 30 under this operator"
          },
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  }
  contract_rows[[1L]] <- append_contract(
    subset(model, timing_operator == "boundary_collapsed"),
    "boundary_collapsed", "recommended_model_comparable",
    subset(model, timing_operator == "boundary_collapsed")$support_definition,
    "four-year model-cell midpoints 20,24,...,44; data and live model use these labels directly"
  )
  contract_rows[[2L]] <- append_contract(
    subset(model, timing_operator == "common_support_18_45"),
    "common_support_18_45", "support_sensitivity",
    subset(model, timing_operator == "common_support_18_45")$support_definition,
    "four-year model-cell midpoints 20,24,...,44; data and live model use these labels directly"
  )
  contract_rows[[3L]] <- append_contract(
    exact, "exact_age_observed", "reference_not_model_comparable",
    rep("observed exact ages 12--49; cohort right edge is truncated by the 2023 data endpoint", nrow(exact)),
    "mother's exact integer age"
  )

  # Keep the original six-column file byte-compatible in schema for downstream
  # readers, but also write an explicitly named exact-age reference with its
  # uncertainty receipt.
  write.csv(exact[, c(
    "window", "cohort_range", "mean_age_first_birth", "share_30plus",
    "n_first_births", "truncation_estimated_share"
  )], targets_path, row.names = FALSE)
  write.csv(exact, exact_targets_path, row.names = FALSE)
  write.csv(model, model_targets_path, row.names = FALSE)
  write.csv(do.call(rbind, contract_rows), contract_path, row.names = FALSE)
  write.csv(do.call(rbind, bin_counts), bin_counts_path, row.names = FALSE)
  write.csv(do.call(rbind, truncation), truncation_path, row.names = FALSE)
}

write_readme <- function(files) {
  exclusions <- read.csv(exclusions_path, stringsAsFactors = FALSE)
  exact <- read.csv(exact_targets_path, stringsAsFactors = FALSE)
  model <- read.csv(model_targets_path, stringsAsFactors = FALSE)
  recommended <- subset(model, timing_operator == "boundary_collapsed")
  common <- subset(model, timing_operator == "common_support_18_45")
  rec_primary <- subset(recommended, window == "primary_1979_1984")
  rec_sensitivity <- subset(recommended, window == "sensitivity_1975_1980")
  common_primary <- subset(common, window == "primary_1979_1984")
  common_sensitivity <- subset(common, window == "sensitivity_1975_1980")
  exact_primary <- subset(exact, window == "primary_1979_1984")
  exact_sensitivity <- subset(exact, window == "sensitivity_1975_1980")
  fmt <- function(x) format(x, scientific = FALSE, trim = TRUE)
  rows <- vapply(seq_len(nrow(exclusions)), function(i) sprintf(
    "| %d | %s | %s | %s | %.4f%% |",
    exclusions$year[i], exclusions$age_variable[i], exclusions$order_variable[i],
    fmt(exclusions$n_unknown_live_birth_order[i]),
    100 * exclusions$unknown_live_birth_order_share[i]
  ), character(1))
  lines <- c(
    "# NCHS natality first-birth timing targets",
    "",
    "This directory builds first-birth timing targets from the local NCHS natality microdata archive. The archive is not part of this repository and is read-only:",
    "`/Users/tommasodesanto/Desktop/Projects/Datasets/natality_data/YYYY/`.",
    "Only national `us` files are used; `ps` territory and `_cwed` files are ignored. The target contract is deliberately frozen to data years 1987--2023. This start date observes every assigned-cohort age cell from age 12 through the age attainable in 2023 for both the 1975--1980 comparison cohorts and the 1979--1984 target cohorts; it does not remove right censoring for the younger cohorts. Later files cannot silently change the contract.",
    "",
    "## Reproduce",
    "",
    "From the repository root, run:",
    "",
    "```sh",
    "Rscript code/data/nchs_natality_timing/build_first_birth_timing.R",
    "```",
    "",
    "Set `NCHS_NATALITY_ROOT` to use another local archive. Add `--force` to rebuild the cache. Otherwise the builder skips the heavy microdata pass only when `first_birth_counts_manifest.csv` exactly matches the current US-file list, byte sizes, and raw-file SHA-256 hashes.",
    "",
    "## Construction",
    "",
    "A first birth is a record with live-birth order equal to 1. The builder retains mothers aged 12--49, excludes missing/unknown live-birth order (code 99 in legacy `dlivord`; code 9 in `lbo_rec`), and counts first births by calendar year and single-year age. It reads only the two required columns from each `.dta` file through `haven::read_dta(..., col_select = ...)`; it never loads a full natality record file.",
    "",
    "| Data years | Mother's single-year age | Live-birth order |",
    "| --- | --- | --- |",
    "| 1987--2002 | `dmage` | `dlivord` |",
    "| 2003 | `mager41`, decoded from the official 41-category recode | `lbo_rec` |",
    "| 2004--2023 | `mager` | `lbo_rec` |",
    "",
    "The 1989 certificate revision changes the source of maternal age toward completed age derived from date of birth when available; it does not change the completed-age concept or the definition of detail live-birth order. The annual first-birth distributions and unknown-order shares are smooth across this source-item transition, so no adjustment is applied. The relevant official documentation receipts are `Nat1987doc.pdf` (SHA-256 `145d8496d764ee66583ef1520ae56804e1c3373bea35d3a570ecbc1426bd5c79`) and `Nat1989doc.pdf` (SHA-256 `92dab8115baec71eec3633239cbd042b2079ad6b80bd1b3a3a43c3276ac3a7cb`).",
    "",
    "The cohort convention is `cohort = data_year - age`. Because both date of birth and age are measured within calendar years, this is an integer approximation with a mechanical +/-1-year cohort-assignment blur.",
    "",
    "The inherited cache previously treated 2003 `MAGER41` codes as literal ages. That was wrong: the official NCHS 2003 detail-file documentation defines code 01 as under age 15, code 02 as age 15, and subsequent codes in single years. The corrected cache decodes codes 02--41 as `code + 13`; code 01 is placed at 14 only to preserve its below-first-cutoff status. The target cohorts are ages 19--28 in 2003, where the recode is single-year exact. The official file is [Nat2003doc.pdf](https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/Nat2003doc.pdf), SHA-256 `82246197e30d54c56a69314bfdcb8f553e6ca0a0d509f8ce8112c2e996b5b2f5`.",
    "",
    "`timing_targets.csv` preserves the original exact-age output for backward compatibility, and `timing_targets_exact_age.csv` gives the same reference moments plus their uncertainty receipt. Exact-age moments are not directly comparable with a model that records births only in four-year cells.",
    "",
    "## Model-comparable timing contract",
    "",
    "The recommended minimal-change operator uses every observed first birth. Exact ages below 22 map to the model's first cell (start age 18); ages 22--25 map to 22; ages 26--29 to 26; ages 30--33 to 30; ages 34--37 to 34; ages 38--41 to 38; and exact ages 42 or older to the terminal cell (start age 42). This is called `boundary_collapsed` in the machine-readable files.",
    "",
    "Both data and the live model measurement are labeled by the continuous-age interval midpoints 20, 24, ..., 44. Relative to an archived period-start report, midpoint coding adds exactly 2 to both sides and therefore changes neither the model--data gap nor the loss. The age-30 threshold is aligned with a bin boundary: `bin_start >= 30` is identically the same classification as exact age at least 30 under both reported operators.",
    "",
    "| Contract | Cohorts | Mean first-birth age | Share age 30+ | First births used |",
    "| --- | --- | ---: | ---: | ---: |",
    sprintf("| Recommended boundary collapse | 1979--1984 | %.12f | %.12f | %s |", rec_primary$mean_age_first_birth, rec_primary$share_30plus, fmt(rec_primary$n_first_births_used)),
    sprintf("| Boundary-collapse comparison window | 1975--1980 | %.12f | %.12f | %s |", rec_sensitivity$mean_age_first_birth, rec_sensitivity$share_30plus, fmt(rec_sensitivity$n_first_births_used)),
    sprintf("| Common support, ages 18--45 | 1979--1984 | %.12f | %.12f | %s |", common_primary$mean_age_first_birth, common_primary$share_30plus, fmt(common_primary$n_first_births_used)),
    sprintf("| Common-support comparison window | 1975--1980 | %.12f | %.12f | %s |", common_sensitivity$mean_age_first_birth, common_sensitivity$share_30plus, fmt(common_sensitivity$n_first_births_used)),
    sprintf("| Exact-age reference | 1979--1984 | %.12f | %.12f | %s |", exact_primary$mean_age_first_birth, exact_primary$share_30plus, fmt(exact_primary$n_first_births)),
    sprintf("| Exact-age comparison window | 1975--1980 | %.12f | %.12f | %s |", exact_sensitivity$mean_age_first_birth, exact_sensitivity$share_30plus, fmt(exact_sensitivity$n_first_births)),
    "",
    sprintf("Conditioning on the literal common support 18--45 is retained as a sensitivity, not the recommendation. It excludes %s observed first births (%.4f percent) from the primary window, almost entirely teen births, and therefore changes the age-30 share's denominator from %.6f to %.6f. Boundary collapse instead preserves the unconditional observed-mothers denominator used alongside completed fertility and childlessness. Its cost is transparent coarsening: teen births are dated in the first model cell and births at 42 or older in the last. A future model with parity at age-18 entry could treat the teen-birth stock directly; the current target does not claim to do so.", fmt(common_primary$n_first_births_excluded), 100 * common_primary$observed_excluded_share, rec_primary$share_30plus, common_primary$share_30plus),
    "",
    sprintf("The primary-versus-comparison-window spread is recomputed after applying each operator. Both cohort windows are observed from age 12 onward in the 1987--2023 cache, so this comparison is not contaminated by differential left censoring. Under the recommended boundary collapse the spread is %.6f year for the mean and %.6f for the age-30 share. Following the existing conservative convention, the mean spread is rounded up to the nearest 0.05 year and the share spread to the nearest 0.001, giving standard errors %.3f and %.3f and inverse-variance weights %.6f and %.6f. Under common support the corresponding standard errors are %.3f and %.3f, with weights %.6f and %.6f. These are window-stability measures, not sampling standard errors.", rec_primary$primary_sensitivity_spread_mean_age, rec_primary$primary_sensitivity_spread_share_30plus, rec_primary$declared_se_mean_age, rec_primary$declared_se_share_30plus, rec_primary$inverse_variance_weight_mean_age, rec_primary$inverse_variance_weight_share_30plus, common_primary$declared_se_mean_age, common_primary$declared_se_share_30plus, common_primary$inverse_variance_weight_mean_age, common_primary$inverse_variance_weight_share_30plus),
    "",
    "Machine-readable artifacts:",
    "",
    "- `timing_target_contract.csv`: long-form activation contract; the recommended rows are labeled `recommended_model_comparable`.",
    "- `timing_targets_model_comparable.csv`: wide target, support, tail, uncertainty, and weight receipt for both model operators.",
    "- `timing_age_bin_counts.csv`: every window-by-operator-by-bin count and share.",
    "- `timing_targets_exact_age.csv`: exact-age reference and uncertainty receipt.",
    "- `timing_targets.csv`: backward-compatible six-column exact-age reference.",
    "- `timing_target_metadata.json`: fail-closed operator, source-document, file-hash, source-bundle, and contract-bundle provenance.",
    "",
    sprintf("Activation requires changing the target contract and its fingerprint: use the recommended primary target %.12f for the midpoint-labeled mean and %.12f for the age-30 share, with standard errors %.3f and %.3f and weights %.12g and %.12g. The live model measurement is already midpoint-coded and uses the aligned age-30 boundary. No structural model or fertility profile is changed by this builder.", rec_primary$mean_age_first_birth, rec_primary$share_30plus, rec_primary$declared_se_mean_age, rec_primary$declared_se_share_30plus, rec_primary$inverse_variance_weight_mean_age, rec_primary$inverse_variance_weight_share_30plus),
    "",
    "## Order exclusions",
    "",
    "The table below reports unknown-order exclusions among all age-12--49 birth records. The same values are available in machine-readable form in `order_exclusion_shares_by_year.csv`.",
    "",
    "| Year | Age variable | Order variable | Unknown order records | Share of age-12--49 records |",
    "| --- | --- | --- | ---: | ---: |",
    rows,
    "",
    "## Truncation",
    "",
    sprintf("For each target cohort, `cohort_truncation.csv` reports the maximum observed first-birth age. Its estimated missed tail is the share of first births above that age in the pooled 1970--1974 reference cohorts. Under the minimal 1987 start, that reference pool lacks some early-age cells, so this is an uncorrected right-tail diagnostic rather than an input to either target or weight. The primary boundary-collapse estimate is %.4f percent; the corresponding common-support diagnostic is %.4f percent. Neither target is corrected or imputed. An observed count of zero above age 45 in the primary window does not mean the true tail is zero: later ages for these cohorts occur after the fixed 2023 data endpoint.", 100 * rec_primary$estimated_unobserved_right_tail_share, 100 * common_primary$estimated_unobserved_right_tail_share)
  )
  writeLines(lines, readme_path)
}

write_metadata <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE) ||
      !requireNamespace("digest", quietly = TRUE)) {
    stop("Packages jsonlite and digest are required to write target provenance.")
  }
  relative_files <- c(
    "first_birth_counts_year_age.csv",
    "first_birth_counts_manifest.csv",
    "order_exclusion_shares_by_year.csv",
    "cohort_truncation.csv",
    "timing_targets.csv",
    "timing_targets_exact_age.csv",
    "timing_targets_model_comparable.csv",
    "timing_target_contract.csv",
    "timing_age_bin_counts.csv",
    "build_first_birth_timing.R",
    "test_first_birth_timing_targets.R"
  )
  absolute_files <- file.path(out_dir, relative_files)
  if (any(!file.exists(absolute_files))) {
    stop("Cannot write provenance; missing files: ",
         paste(relative_files[!file.exists(absolute_files)], collapse = ", "))
  }
  file_sha256 <- vapply(
    absolute_files,
    function(path) digest::digest(path, algo = "sha256", file = TRUE),
    character(1)
  )
  names(file_sha256) <- relative_files

  official_1987_document_sha256 <- "145d8496d764ee66583ef1520ae56804e1c3373bea35d3a570ecbc1426bd5c79"
  official_1989_document_sha256 <- "92dab8115baec71eec3633239cbd042b2079ad6b80bd1b3a3a43c3276ac3a7cb"
  official_document_sha256 <- "82246197e30d54c56a69314bfdcb8f553e6ca0a0d509f8ce8112c2e996b5b2f5"
  source_members <- c(
    first_birth_counts_year_age.csv = file_sha256[["first_birth_counts_year_age.csv"]],
    first_birth_counts_manifest.csv = file_sha256[["first_birth_counts_manifest.csv"]],
    official_1987_nchs_documentation_pdf = official_1987_document_sha256,
    official_1989_nchs_documentation_pdf = official_1989_document_sha256,
    official_2003_nchs_documentation_pdf = official_document_sha256,
    cache_contract = cache_contract
  )
  source_bundle_sha256 <- digest::digest(
    paste(names(source_members), source_members, sep = "=", collapse = "\n"),
    algo = "sha256", serialize = FALSE
  )
  contract_member_names <- c(
    "timing_target_contract.csv", "timing_targets_model_comparable.csv",
    "timing_age_bin_counts.csv", "timing_targets_exact_age.csv",
    "build_first_birth_timing.R", "test_first_birth_timing_targets.R"
  )
  contract_members <- file_sha256[contract_member_names]
  contract_bundle_sha256 <- digest::digest(
    paste(names(contract_members), contract_members, sep = "=", collapse = "\n"),
    algo = "sha256", serialize = FALSE
  )

  model_targets <- read.csv(model_targets_path, stringsAsFactors = FALSE)
  primary <- subset(
    model_targets,
    window == "primary_1979_1984" & timing_operator == "boundary_collapsed"
  )
  sensitivity <- subset(
    model_targets,
    window == "sensitivity_1975_1980" & timing_operator == "boundary_collapsed"
  )
  metadata <- list(
    schema_version = 1L,
    contract_id = "nchs_first_birth_timing_v4_boundary_collapsed_midpoint_1987_2023_lbo_unknown_fixed",
    cache_contract = cache_contract,
    data_vintage = list(year_min = year_min, year_max = year_max),
    contains_five_dated_period_targets = FALSE,
    target_interpretation = "pooled cohort timing target; annual period comparisons must apply the declared operator directly to the year-age cache",
    cohort_windows = list(
      primary = "1979-1984",
      sensitivity = "1975-1980",
      cohort_assignment = "calendar data year minus exact/decoded maternal age; unavoidable plus-or-minus-one calendar blur",
      left_support = "both primary and sensitivity windows are observed from age 12 through each cohort's age attainable in 2023"
    ),
    right_tail_diagnostic = list(
      reference_cohorts = "1970-1974",
      status = "diagnostic only; the 1987 start omits some early-age cells in this reference pool",
      use_in_target_or_weight = FALSE
    ),
    recommended_operator = list(
      id = "boundary_collapsed_midpoint_v1",
      boundary_treatment = "exact/decoded ages <=21 map to cell start 18; 22-25 to 22; 26-29 to 26; 30-33 to 30; 34-37 to 34; 38-41 to 38; ages >=42 map to cell start 42",
      bin_starts = unname(model_bin_starts),
      reported_midpoint_labels = unname(model_bin_starts + 2L),
      mean_formula = "first-birth-count-weighted mean of midpoint labels 20,24,...,44",
      share_30plus_formula = "sum first-birth counts with bin_start >= 30 divided by all used first-birth counts",
      share_30plus_classification_identity = "bin_start >= 30 is identical to exact/decoded age >= 30",
      observed_age_support = "all observed ages in the cache (nominally 12-49; 2003 MAGER41 code 01 is an under-15 group placed below the first cutoff)",
      boundary_exclusions = "none",
      primary_target = list(
        mean_age_first_birth = primary$mean_age_first_birth,
        share_30plus = primary$share_30plus,
        n_first_births = primary$n_first_births_used,
        mean_standard_error = primary$declared_se_mean_age,
        share_30plus_standard_error = primary$declared_se_share_30plus,
        mean_inverse_variance_weight = primary$inverse_variance_weight_mean_age,
        share_30plus_inverse_variance_weight = primary$inverse_variance_weight_share_30plus
      ),
      sensitivity_target = list(
        mean_age_first_birth = sensitivity$mean_age_first_birth,
        share_30plus = sensitivity$share_30plus,
        n_first_births = sensitivity$n_first_births_used
      )
    ),
    uncertainty = list(
      definition = "absolute primary-versus-sensitivity cohort-window spread under the same operator; both windows are observed from age 12 onward; not a sampling standard error",
      mean_rounding = "round upward to the nearest 0.05 year",
      share_rounding = "round upward to the nearest 0.001"
    ),
    raw_microdata_receipt = list(
      status = "every 1987--2023 national raw file is content-hashed before cache reuse",
      manifest = "first_birth_counts_manifest.csv",
      manifest_columns = c("year", "file", "size_bytes", "sha256", "cache_contract")
    ),
    live_birth_order_codes = list(
      first_birth = 1L,
      legacy_dlivord_unknown = 99L,
      lbo_rec_unknown = 9L,
      missing_values_also_excluded = TRUE
    ),
    pre_2003_measurement = list(
      maternal_age_field = "dmage",
      maternal_age_field_label = "Age of Mother",
      live_birth_order_field = "dlivord",
      live_birth_order_field_label = "Detail Live Birth Order",
      years = "1987--2002",
      source_item_transition = "1989 revision moves maternal age toward completed age derived from date of birth when available; no concept or coding adjustment applied",
      official_1987_document_url = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/Nat1987doc.pdf",
      official_1987_document_sha256 = official_1987_document_sha256,
      official_1989_document_url = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/Nat1989doc.pdf",
      official_1989_document_sha256 = official_1989_document_sha256
    ),
    maternal_age_source_2003 = list(
      field = "MAGER41",
      transformation = "code 01 assigned 14 only to preserve its under-15 boundary status; codes 02:41 decoded as code+13",
      relevance = "primary and sensitivity cohorts are ages 19-28 in 2003, where MAGER41 is single-year exact after decoding",
      official_document_url = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/DVS/natality/Nat2003doc.pdf",
      official_document_sha256 = official_document_sha256,
      official_document_size_bytes = 16116001
    ),
    file_sha256 = as.list(file_sha256),
    source_bundle_members = as.list(source_members),
    source_bundle_sha256 = source_bundle_sha256,
    contract_bundle_members = as.list(contract_members),
    contract_bundle_sha256 = contract_bundle_sha256
  )
  jsonlite::write_json(
    metadata, metadata_path, pretty = TRUE, auto_unbox = TRUE,
    digits = 16, null = "null"
  )
}

files <- find_us_files(archive_root)
manifest <- make_manifest(files)
if (force_rebuild || !same_manifest(manifest)) {
  message("Cache missing or archive manifest changed: building year-by-age counts.")
  build_cache(files, manifest)
} else {
  message("Archive manifest matches cache: skipping microdata pass.")
}
counts <- read.csv(cache_path, stringsAsFactors = FALSE)
make_targets(counts)
write_readme(files)
write_metadata()

targets <- read.csv(targets_path, stringsAsFactors = FALSE)
primary <- subset(targets, window == "primary_1979_1984")
cat("\nPRIMARY 1979-1984 COHORT TARGETS\n")
cat(sprintf("Mean age at first birth: %.6f\n", primary$mean_age_first_birth))
cat(sprintf("Share of first births at age >= 30: %.6f\n", primary$share_30plus))
