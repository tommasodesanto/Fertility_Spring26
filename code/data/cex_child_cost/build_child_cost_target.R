#!/usr/bin/env Rscript

# Pilot CEX estimate of the nonhousing expenditure increment associated with
# one additional dependent child.  The model period is four years, while CEX
# expenditure summaries are quarterly and income is annual.

suppressPackageStartupMessages({
  library(haven)
  library(sandwich)
})

weighted_quantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  x[which(cumsum(w) / sum(w) >= p)[1L]]
}

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[1]), "../../.."),
                      mustWork = FALSE)
if (!dir.exists(file.path(root, "code", "data", "cex_child_cost"))) {
  root <- normalizePath(".")
}

years <- as.integer(strsplit(Sys.getenv("CEX_YEARS", "2023"), ",", fixed = TRUE)[[1]])
stopifnot(length(years) >= 1L, all(is.finite(years)))
output_label <- Sys.getenv(
  "CEX_OUTPUT_LABEL",
  if (length(years) == 1L) paste0("child_cost_", years) else
    paste0("child_cost_", min(years), "_", max(years))
)
input_dirs <- file.path(root, "code", "data", "cex_child_cost", "raw", years, "extracted")
output_dir <- file.path(root, "code", "data", "cex_child_cost", "output", output_label)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(all(dir.exists(input_dirs)))
files_by_year <- lapply(input_dirs, list.files, pattern = "^fmli.*[.]dta$", full.names = TRUE)
names(files_by_year) <- years
stopifnot(all(lengths(files_by_year) >= 4L))

keep <- c(
  "CUID", "INTERI", "QINTRVYR", "QINTRVMO", "AGE_REF", "FAM_SIZE",
  "PERSLT18", "CUTENURE", "ROOMSQ", "FINCBTAX", "FINCBTXM", "FINCBTXI",
  "INCLASS2", "FINLWT21",
  "TOTEXPCQ", "HOUSCQ", "CARTKNCQ", "CARTKUCQ", "OTHVEHCQ",
  "CASHCOCQ", "PERINSCQ"
)

read_fmli <- function(path, package_year) {
  actual_names <- names(read_dta(path, n_max = 0))
  selected <- actual_names[match(keep, toupper(actual_names))]
  if (anyNA(selected)) {
    stop("Missing required columns in ", path, ": ",
         paste(keep[is.na(selected)], collapse = ", "))
  }
  x <- as.data.frame(read_dta(path, col_select = all_of(selected)))
  names(x) <- toupper(names(x))
  x$PACKAGE_YEAR <- package_year
  x
}

d <- do.call(rbind, unlist(Map(
  function(paths, package_year) lapply(paths, read_fmli, package_year = package_year),
  files_by_year, years
), recursive = FALSE))

cpi_u <- c(
  `2019` = 255.657, `2020` = 258.811, `2021` = 270.970,
  `2022` = 292.655, `2023` = 304.702, `2024` = 313.689
)
if (any(!as.character(d$QINTRVYR) %in% names(cpi_u))) {
  stop("Missing CPI-U value for an interview year.")
}
deflator <- unname(cpi_u[["2023"]] / cpi_u[as.character(d$QINTRVYR)])
monetary <- c(
  "FINCBTAX", "FINCBTXM", "TOTEXPCQ", "HOUSCQ", "CARTKNCQ", "CARTKUCQ",
  "OTHVEHCQ", "CASHCOCQ", "PERINSCQ"
)
d[monetary] <- lapply(d[monetary], function(x) x * deflator)
income_measure <- toupper(Sys.getenv("CEX_INCOME_MEASURE", "FINCBTXM"))
if (!income_measure %in% c("FINCBTXM", "FINCBTAX")) {
  stop("CEX_INCOME_MEASURE must be FINCBTXM or FINCBTAX.")
}
d$FINCBTAX_COLLECTED <- d$FINCBTAX
d$FINCBTAX <- d[[income_measure]]
complete_reporters_only <- identical(Sys.getenv("CEX_COMPLETE_REPORTERS_ONLY", "0"), "1")
if (complete_reporters_only) d <- d[!is.na(d$INCLASS2) & d$INCLASS2 != 7, ]

d$adults <- d$FAM_SIZE - d$PERSLT18
d$quarter <- interaction(d$QINTRVYR, d$QINTRVMO, drop = TRUE)

# The model's c composite is nonhousing consumption.  Starting from total CEX
# expenditure, remove housing, vehicle purchases, cash contributions, and
# personal insurance/pension contributions.  The remaining composite includes
# food, apparel, vehicle operating costs, health, entertainment, personal care,
# reading, education, tobacco, and miscellaneous consumption.
d$nonhousing_nondurable_q <- with(
  d,
  TOTEXPCQ - HOUSCQ - CARTKNCQ - CARTKUCQ - OTHVEHCQ - CASHCOCQ - PERINSCQ
)
d$nonhousing_nondurable_annual <- 4 * d$nonhousing_nondurable_q

s <- subset(
  d,
  AGE_REF >= 25 & AGE_REF <= 55 &
    PERSLT18 >= 0 & PERSLT18 <= 3 &
    adults >= 1 & FINCBTAX > 1000 &
    nonhousing_nondurable_annual > 0 &
    !is.na(ROOMSQ) & ROOMSQ > 0
)

# Limit sensitivity to expenditure outliers without trimming on the outcome by
# family type.  These are unweighted cutoffs, reported below.
s_pretrim <- s
trim <- quantile(s$nonhousing_nondurable_annual, c(0.01, 0.99), na.rm = TRUE)
s <- subset(
  s,
  nonhousing_nondurable_annual >= trim[[1]] &
    nonhousing_nondurable_annual <= trim[[2]]
)
matched_trim <- c(
  weighted_quantile(
    s_pretrim$nonhousing_nondurable_annual, s_pretrim$FINLWT21, 0.01
  ),
  weighted_quantile(
    s_pretrim$nonhousing_nondurable_annual, s_pretrim$FINLWT21, 0.99
  )
)
s_one_shot <- subset(
  s_pretrim,
  nonhousing_nondurable_annual >= matched_trim[[1]] &
    nonhousing_nondurable_annual <= matched_trim[[2]]
)

s$log_income <- log(s$FINCBTAX)
s$log_income_sq <- s$log_income^2
s$age_sq <- s$AGE_REF^2
s$any_child <- as.integer(s$PERSLT18 > 0)
s$additional_children <- pmax(s$PERSLT18 - 1, 0)
s$one_shot_parity <- ifelse(s$PERSLT18 == 0, 0,
                            ifelse(s$PERSLT18 <= 2, 1, 2))
s_one_shot$log_income <- log(s_one_shot$FINCBTAX)
s_one_shot$log_income_sq <- s_one_shot$log_income^2
s_one_shot$age_sq <- s_one_shot$AGE_REF^2
s_one_shot$one_shot_parity <- ifelse(
  s_one_shot$PERSLT18 == 0, 0, ifelse(s_one_shot$PERSLT18 <= 2, 1, 2)
)
s_one_shot$owner <- as.integer(s_one_shot$CUTENURE %in% c(1, 2))

common_controls <- paste(
  "log_income + log_income_sq + AGE_REF + age_sq + adults +",
  "factor(CUTENURE) + factor(quarter)"
)

fit_model <- function(rhs) {
  formula <- as.formula(paste("nonhousing_nondurable_annual ~", rhs))
  lm(formula, data = s, weights = FINLWT21)
}

models <- list(
  linear_no_rooms = fit_model(paste("PERSLT18 +", common_controls)),
  linear_rooms = fit_model(paste("PERSLT18 + ROOMSQ +", common_controls)),
  first_vs_additional = fit_model(paste(
    "any_child + additional_children + ROOMSQ +", common_controls
  )),
  one_shot_model_feasible = lm(
    nonhousing_nondurable_annual ~ one_shot_parity + ROOMSQ +
      log_income + log_income_sq + AGE_REF + age_sq + owner,
    data = s_one_shot, weights = FINLWT21
  )
)

extract_term <- function(model_name, term) {
  model <- models[[model_name]]
  vc <- vcovCL(model, cluster = ~CUID, type = "HC1")
  estimate <- unname(coef(model)[[term]])
  se <- sqrt(diag(vc))[[term]]
  data.frame(
    specification = model_name,
    term = term,
    estimate_annual_dollars = estimate,
    cluster_se = se,
    ci_low = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se,
    stringsAsFactors = FALSE
  )
}

results <- rbind(
  extract_term("linear_no_rooms", "PERSLT18"),
  extract_term("linear_rooms", "PERSLT18"),
  extract_term("first_vs_additional", "any_child"),
  extract_term("first_vs_additional", "additional_children"),
  extract_term("one_shot_model_feasible", "one_shot_parity")
)

weighted_mean_income <- weighted.mean(s$FINCBTAX, s$FINLWT21)
weighted_mean_consumption <- weighted.mean(
  s$nonhousing_nondurable_annual, s$FINLWT21
)
period_years <- 4
results$share_of_mean_income <- results$estimate_annual_dollars / weighted_mean_income
results$model_period_units <- period_years * results$share_of_mean_income
one_shot_row <- results$specification == "one_shot_model_feasible"
one_shot_mean_income <- weighted.mean(s_one_shot$FINCBTAX, s_one_shot$FINLWT21)
results$share_of_mean_income[one_shot_row] <-
  results$estimate_annual_dollars[one_shot_row] / one_shot_mean_income
results$model_period_units[one_shot_row] <-
  period_years * results$share_of_mean_income[one_shot_row]

sample_summary <- data.frame(
  statistic = c(
    "interview_observations", "unique_consumer_units", "weighted_mean_income",
    "weighted_mean_nonhousing_nondurable", "p01_expenditure_trim",
    "p99_expenditure_trim"
  ),
  value = c(
    nrow(s), length(unique(s$CUID)), weighted_mean_income,
    weighted_mean_consumption, trim[[1]], trim[[2]]
  )
)

child_cells <- do.call(rbind, lapply(0:3, function(n) {
  x <- s[s$PERSLT18 == n, ]
  data.frame(
    children_under_18 = n,
    observations = nrow(x),
    unique_consumer_units = length(unique(x$CUID)),
    weighted_mean_income = weighted.mean(x$FINCBTAX, x$FINLWT21),
    weighted_mean_nonhousing_nondurable = weighted.mean(
      x$nonhousing_nondurable_annual, x$FINLWT21
    )
  )
}))

write.csv(results, file.path(output_dir, "cex_child_cost_regressions.csv"), row.names = FALSE)
write.csv(sample_summary, file.path(output_dir, "cex_child_cost_sample.csv"), row.names = FALSE)
write.csv(child_cells, file.path(output_dir, "cex_child_cost_cells.csv"), row.names = FALSE)

readme <- c(
  "# CEX child-cost pilot target",
  "",
  paste0("Source: ", min(years), "--", max(years),
         " Consumer Expenditure Survey Interview PUMD; monetary values are 2023 dollars."),
  paste0("Income measure: ", income_measure,
         if (complete_reporters_only) "; complete income reporters only." else
           "; all reporters using BLS collected-or-imputed income."),
  "Sample: reference person age 25--55, 0--3 persons under 18, positive annual",
  "before-tax income, positive nonhousing nondurable expenditure, and valid rooms.",
  "Observations are interview-quarter consumer units and regressions use FINLWT21.",
  "Standard errors are clustered by CUID.",
  "",
  "Nonhousing nondurable expenditure equals four times quarterly total expenditure",
  "less housing, new/used/other vehicle purchases, cash contributions, and personal",
  "insurance/pension contributions. The primary pilot is `linear_rooms`.",
  "The `one_shot_model_feasible` row instead uses the model's 0 / 1--2 / 3+",
  "parity index and only covariates available in the one-shot model: income, age,",
  "owner/renter tenure, and rooms. It uses weighted 1st/99th-percentile expenditure trims,",
  "matching the deterministic model-distribution operation. This is the candidate",
  "matched auxiliary for c_bar_n.",
  "",
  "`model_period_units` divides the annual-dollar coefficient by weighted mean annual",
  "income and multiplies by the model's four-year period length.",
  "",
  "Important limitation: the coefficient is the net change in total household spending",
  "conditional on resources. It is not a structural equivalence scale: parents can fund",
  "child goods by reducing their own consumption. Treat it as a diagnostic until an",
  "Engel/Rothbarth child-cost target is selected or this interpretation is defended."
)
writeLines(readme, file.path(output_dir, "README.md"))

cat("Primary annual coefficient:",
    results$estimate_annual_dollars[results$specification == "linear_rooms"], "\n")
cat("Primary model-period units:",
    results$model_period_units[results$specification == "linear_rooms"], "\n")
