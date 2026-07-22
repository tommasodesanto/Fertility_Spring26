#!/usr/bin/env Rscript

# Pilot CEX estimate of the nonhousing expenditure increment associated with
# one additional dependent child.  The model period is four years, while CEX
# expenditure summaries are quarterly and income is annual.

suppressPackageStartupMessages({
  library(haven)
  library(sandwich)
})

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[1]), "../../.."),
                      mustWork = FALSE)
if (!dir.exists(file.path(root, "code", "data", "cex_child_cost"))) {
  root <- normalizePath(".")
}

input_dir <- file.path(root, "code", "data", "cex_child_cost", "raw", "2023", "extracted")
output_dir <- file.path(root, "code", "data", "cex_child_cost", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "^fmli.*[.]dta$", full.names = TRUE)
stopifnot(length(files) == 4L)

keep <- c(
  "CUID", "INTERI", "QINTRVYR", "QINTRVMO", "AGE_REF", "FAM_SIZE",
  "PERSLT18", "CUTENURE", "ROOMSQ", "FINCBTAX", "FINLWT21",
  "TOTEXPCQ", "HOUSCQ", "CARTKNCQ", "CARTKUCQ", "OTHVEHCQ",
  "CASHCOCQ", "PERINSCQ"
)

d <- do.call(rbind, lapply(files, function(path) {
  as.data.frame(read_dta(path, col_select = all_of(keep)))
}))

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
trim <- quantile(s$nonhousing_nondurable_annual, c(0.01, 0.99), na.rm = TRUE)
s <- subset(
  s,
  nonhousing_nondurable_annual >= trim[[1]] &
    nonhousing_nondurable_annual <= trim[[2]]
)

s$log_income <- log(s$FINCBTAX)
s$log_income_sq <- s$log_income^2
s$age_sq <- s$AGE_REF^2
s$any_child <- as.integer(s$PERSLT18 > 0)
s$additional_children <- pmax(s$PERSLT18 - 1, 0)

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
  ))
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
  extract_term("first_vs_additional", "additional_children")
)

weighted_mean_income <- weighted.mean(s$FINCBTAX, s$FINLWT21)
weighted_mean_consumption <- weighted.mean(
  s$nonhousing_nondurable_annual, s$FINLWT21
)
period_years <- 4
results$share_of_mean_income <- results$estimate_annual_dollars / weighted_mean_income
results$model_period_units <- period_years * results$share_of_mean_income

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
  "Source: 2023 Consumer Expenditure Survey Interview PUMD, interviews 2023Q2--2024Q1.",
  "Sample: reference person age 25--55, 0--3 persons under 18, positive annual",
  "before-tax income, positive nonhousing nondurable expenditure, and valid rooms.",
  "Observations are interview-quarter consumer units and regressions use FINLWT21.",
  "Standard errors are clustered by CUID.",
  "",
  "Nonhousing nondurable expenditure equals four times quarterly total expenditure",
  "less housing, new/used/other vehicle purchases, cash contributions, and personal",
  "insurance/pension contributions. The primary pilot is `linear_rooms`.",
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
