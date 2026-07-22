#!/usr/bin/env Rscript

# PSID auxiliary moment for the dispersion of tenure choices.  Predict ownership
# four years ahead from observables available at the initial date, using
# person-level cross-fitting, then measure the weighted mean squared residual.
# The identical auxiliary prediction exercise must be run on simulated model
# data before this moment is used to discipline tenure_choice_kappa.

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_quantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  x[which(cumsum(w) / sum(w) >= p)[1L]]
}

script_dir <- function() {
  flag <- commandArgs(trailingOnly = FALSE)
  flag <- flag[startsWith(flag, "--file=")]
  if (length(flag) == 0L) stop("Run with Rscript.")
  dirname(normalizePath(sub("--file=", "", flag[1L], fixed = TRUE)))
}

root <- normalizePath(file.path(script_dir(), "..", "..", ".."))
psid_path <- file.path(dirname(root), "PSID", "PSIDSHELF_MOBILITY.dta")
out_dir <- file.path(root, "code", "data", "psid_followup_mar2026", "output",
                     "tenure_residual_variance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading selected PSID columns...")
raw <- as.data.table(read_dta(
  psid_path,
  col_select = c(
    ID, year, AGEREP, DEATHYEAR, RELTOHEAD_, INCFAMR, HOMEOWN, IW,
    FAMCHILD, FAMMARRIED, WLTHSAVETOTR, WLTHFUNDTOTR, WLTHODEBTOTR
  )
))

dt <- raw[, .(
  id = as.numeric(ID),
  year = as.integer(year),
  age = as.numeric(AGEREP),
  death_year = as.numeric(DEATHYEAR),
  reference_person = as.numeric(RELTOHEAD_) == 10,
  income = as.numeric(INCFAMR),
  own = fifelse(as.numeric(HOMEOWN) == 1, 1,
                fifelse(as.numeric(HOMEOWN) == 2, 0, NA_real_)),
  weight = as.numeric(IW),
  children = as.numeric(FAMCHILD),
  married = as.numeric(FAMMARRIED) == 1,
  liquid_financial = as.numeric(WLTHSAVETOTR) + as.numeric(WLTHFUNDTOTR) -
    as.numeric(WLTHODEBTOTR)
)]
rm(raw)
invisible(gc())

dt <- dt[year >= 1984 & year <= 2019 &
           (is.na(death_year) | year <= death_year)]
setkey(dt, id, year)

# HOMEOWN is a family-level variable repeated on person records.  Initial
# observations must be reference persons; the same person's family tenure four
# years later supplies the outcome even if reference-person status changes.
future <- dt[, .(id, year = year - 4L, own4 = own)]
panel <- merge(dt[reference_person == TRUE], future, by = c("id", "year"), all = FALSE)
panel <- panel[
  age >= 25 & age <= 55 & income > 1000 & weight > 0 &
    is.finite(own) & is.finite(own4) & is.finite(children) &
    is.finite(married) & is.finite(liquid_financial)
]

panel[, liquid_to_income := liquid_financial / income]
lo <- weighted_quantile(panel$liquid_to_income, panel$weight, 0.01)
hi <- weighted_quantile(panel$liquid_to_income, panel$weight, 0.99)
panel[, liquid_to_income_w := pmin(pmax(liquid_to_income, lo), hi)]
panel[, log_income := log(income)]
panel[, children_capped := pmin(as.integer(children), 3L)]
panel[, fold := as.integer(id %% 5L)]

formula_aux <- own4 ~ own + log_income + I(log_income^2) + age + I(age^2) +
  liquid_to_income_w + I(liquid_to_income_w^2) + factor(children_capped) +
  married + factor(year)

panel[, predicted_own4 := NA_real_]
for (ff in 0:4) {
  fit <- glm(
    formula_aux,
    data = panel[fold != ff],
    weights = weight,
    family = quasibinomial(link = "logit")
  )
  panel[fold == ff, predicted_own4 := predict(
    fit, newdata = panel[fold == ff], type = "response"
  )]
}

stopifnot(all(is.finite(panel$predicted_own4)))
panel[, squared_residual := (own4 - predicted_own4)^2]
panel[, raw_change := (own4 - own)^2]

brier <- weighted_mean_safe(panel$squared_residual, panel$weight)
raw_change_mse <- weighted_mean_safe(panel$raw_change, panel$weight)
ownership_rate <- weighted_mean_safe(panel$own4, panel$weight)
prediction_variance <- weighted_mean_safe(
  (panel$predicted_own4 - weighted_mean_safe(panel$predicted_own4, panel$weight))^2,
  panel$weight
)

# Person bootstrap of the already cross-fitted losses.  This treats the fitted
# prediction rule as fixed; the README therefore labels the standard error as
# conditional rather than as a full generated-regressor bootstrap.
set.seed(20260722)
id_loss <- panel[, .(
  weighted_loss = sum(weight * squared_residual),
  total_weight = sum(weight)
), by = id]
boot <- replicate(500L, {
  ii <- sample.int(nrow(id_loss), nrow(id_loss), replace = TRUE)
  sum(id_loss$weighted_loss[ii]) / sum(id_loss$total_weight[ii])
})

moments <- data.table(
  moment = c(
    "crossfitted_four_year_tenure_residual_variance",
    "four_year_current_tenure_prediction_mse",
    "ownership_rate_four_year_ahead",
    "variance_crossfitted_predicted_ownership"
  ),
  value = c(brier, raw_change_mse, ownership_rate, prediction_variance),
  standard_error = c(sd(boot), NA_real_, NA_real_, NA_real_),
  ci_lower = c(unname(quantile(boot, 0.025)), NA_real_, NA_real_, NA_real_),
  ci_upper = c(unname(quantile(boot, 0.975)), NA_real_, NA_real_, NA_real_),
  status = c(
    "candidate primary target for tenure_choice_kappa",
    "benchmark: prediction using current tenure alone",
    "diagnostic; ownership level targets chi",
    "diagnostic: observable sorting component"
  )
)

cells <- panel[, .(
  observations = .N,
  unique_people = uniqueN(id),
  ownership_rate = weighted_mean_safe(own4, weight),
  predicted_ownership = weighted_mean_safe(predicted_own4, weight),
  residual_variance = weighted_mean_safe(squared_residual, weight)
), by = .(initial_tenure = own, age_group = cut(
  age, breaks = c(24, 34, 44, 55), labels = c("25-34", "35-44", "45-55")
))]

sample_summary <- data.table(
  statistic = c(
    "person_year_observations", "unique_people", "first_year", "last_initial_year",
    "weighted_initial_ownership", "weighted_future_ownership",
    "liquid_wealth_ratio_p01", "liquid_wealth_ratio_p99"
  ),
  value = c(
    nrow(panel), uniqueN(panel$id), min(panel$year), max(panel$year),
    weighted_mean_safe(panel$own, panel$weight), ownership_rate, lo, hi
  )
)

fwrite(moments, file.path(out_dir, "psid_tenure_residual_variance_moments.csv"))
fwrite(cells, file.path(out_dir, "psid_tenure_residual_variance_cells.csv"))
fwrite(sample_summary, file.path(out_dir, "psid_tenure_residual_variance_sample.csv"))

writeLines(c(
  "# PSID four-year residual tenure variance",
  "",
  "Primary candidate: the weighted Brier score from a five-fold cross-fitted logit",
  "for ownership four years ahead. Initial covariates are current tenure, liquid",
  "financial wealth/income, family income, age, number of children, marital status,",
  "and survey-year effects. The sample contains PSID reference persons ages 25--55.",
  "",
  "Economic interpretation: chi is assigned the ownership level; tenure_choice_kappa",
  "is assigned the residual dispersion in tenure after conditioning on observed states.",
  "The identical auxiliary regression and Brier score must be computed on simulated",
  "model histories. The wealth-gradient and switching-rate moments remain validation",
  "outcomes rather than substitutes for this residual-variation target.",
  "",
  "The reported person-bootstrap standard error conditions on the cross-fitted",
  "prediction rule. A final target should use a full bootstrap that re-estimates the",
  "five auxiliary logits in every draw."
), file.path(out_dir, "README.md"))

print(moments)
print(sample_summary)
