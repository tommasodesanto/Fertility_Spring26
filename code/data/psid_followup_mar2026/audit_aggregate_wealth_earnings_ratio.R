#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
})

RECENT_YEAR_MIN <- 2005L
RECENT_YEAR_MAX <- 2019L
MODEL_AGE_MIN <- 18L
MODEL_AGE_MAX <- 85L
MODEL_WORKING_AGE_MAX <- 65L
BOOTSTRAP_REPS <- 999L
BOOTSTRAP_SEED <- 20260723L

AGE_BINS <- data.table(
  age_bin = c("26-35", "36-45", "46-55", "56-65"),
  age_lo = c(26L, 36L, 46L, 56L),
  age_hi = c(35L, 45L, 55L, 65L)
)

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", file_arg[startsWith(file_arg, file_flag)])
  if (length(script_path) == 0L) {
    stop("Run this file with Rscript.")
  }
  dirname(normalizePath(script_path[1L]))
}

safe_ratio <- function(numerator, denominator) {
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  numerator / denominator
}

weighted_components <- function(sample_dt, denominator_name) {
  denominator <- sample_dt[[denominator_name]]
  c(
    numerator = sum(sample_dt$weight * sample_dt$net_worth),
    denominator = sum(sample_dt$weight * denominator)
  )
}

bootstrap_summary <- function(point, draws, reps) {
  data.table(
    estimate = point,
    bootstrap_se = sd(draws, na.rm = TRUE),
    bootstrap_p025 = quantile(draws, 0.025, na.rm = TRUE, names = FALSE),
    bootstrap_p975 = quantile(draws, 0.975, na.rm = TRUE, names = FALSE),
    bootstrap_reps = reps
  )
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
psid_path <- file.path(dirname(repo_root), "PSID", "PSIDSHELF_MOBILITY.dta")
out_dir <- file.path(
  repo_root,
  "code",
  "data",
  "psid_followup_mar2026",
  "output",
  "aggregate_wealth_earnings_audit"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(psid_path)) {
  stop("Missing PSID shelf file: ", psid_path)
}

dt <- as.data.table(read_dta(
  psid_path,
  col_select = c(
    "ID", "year", "RELTOHEAD_", "AGEREP", "NETWORTHR",
    "EARNINDRRC", "INCFAMR", "IW"
  )
))
setnames(
  dt,
  c(
    "ID", "year", "RELTOHEAD_", "AGEREP", "NETWORTHR",
    "EARNINDRRC", "INCFAMR", "IW"
  ),
  c(
    "id", "year", "relation_to_head", "age", "net_worth",
    "gross_labor_earnings", "total_family_income", "weight"
  )
)
for (name in names(dt)) {
  set(dt, j = name, value = as.numeric(dt[[name]]))
}

# A model state labelled 62 represents the four-year working period 62--65,
# and the terminal state labelled 82 represents ages 82--85. The data universe
# therefore uses reference persons 18--85 and labor earnings through age 65.
dt <- dt[
  relation_to_head == 10 &
    year >= RECENT_YEAR_MIN & year <= RECENT_YEAR_MAX &
    age >= MODEL_AGE_MIN & age <= MODEL_AGE_MAX &
    is.finite(weight) & weight > 0 &
    is.finite(net_worth) &
    (
      age > MODEL_WORKING_AGE_MAX |
        (is.finite(gross_labor_earnings) & gross_labor_earnings >= 0)
    )
]

if (nrow(dt) == 0L || uniqueN(dt$id) == 0L) {
  stop("The recent-vintage PSID analysis sample is empty.")
}

aggregate_components <- dt[, .(
  numerator = sum(weight * net_worth),
  denominator = sum(
    weight[age <= MODEL_WORKING_AGE_MAX] *
      gross_labor_earnings[age <= MODEL_WORKING_AGE_MAX]
  ),
  family_income_denominator = sum(
    weight[
      age <= MODEL_WORKING_AGE_MAX & is.finite(total_family_income)
    ] *
      total_family_income[
        age <= MODEL_WORKING_AGE_MAX & is.finite(total_family_income)
      ]
  ),
  family_years = .N,
  persons = uniqueN(id)
)]
aggregate_point <- safe_ratio(
  aggregate_components$numerator,
  aggregate_components$denominator
)

yearly <- dt[, {
  worker <- age <= MODEL_WORKING_AGE_MAX
  .(
    wealth_households = .N,
    earnings_households = sum(worker),
    aggregate_wealth = sum(weight * net_worth),
    aggregate_gross_labor_earnings =
      sum(weight[worker] * gross_labor_earnings[worker]),
    gross_wealth_to_gross_labor_earnings = safe_ratio(
      sum(weight * net_worth),
      sum(weight[worker] * gross_labor_earnings[worker])
    )
  )
}, by = year]

age_profile_point <- rbindlist(lapply(seq_len(nrow(AGE_BINS)), function(index) {
  definition <- AGE_BINS[index]
  sample_dt <- dt[age >= definition$age_lo & age <= definition$age_hi]
  components <- weighted_components(sample_dt, "gross_labor_earnings")
  data.table(
    age_bin = definition$age_bin,
    age_lo = definition$age_lo,
    age_hi = definition$age_hi,
    family_years = nrow(sample_dt),
    persons = uniqueN(sample_dt$id),
    aggregate_wealth = components[["numerator"]],
    aggregate_gross_labor_earnings = components[["denominator"]],
    estimate = safe_ratio(
      components[["numerator"]],
      components[["denominator"]]
    )
  )
}))

# Person-cluster bootstrap: each unique reference-person ID is resampled, and
# all of that person's family-year observations retain their survey weights.
component_names <- c("aggregate", AGE_BINS$age_bin)
by_id <- dt[, {
  worker <- age <= MODEL_WORKING_AGE_MAX
  values <- list(
    aggregate_num = sum(weight * net_worth),
    aggregate_den = sum(weight[worker] * gross_labor_earnings[worker])
  )
  for (index in seq_len(nrow(AGE_BINS))) {
    definition <- AGE_BINS[index]
    in_bin <- age >= definition$age_lo & age <= definition$age_hi
    stem <- gsub("-", "_", definition$age_bin)
    values[[paste0(stem, "_num")]] <- sum(weight[in_bin] * net_worth[in_bin])
    values[[paste0(stem, "_den")]] <-
      sum(weight[in_bin] * gross_labor_earnings[in_bin])
  }
  values
}, by = id]

bootstrap_draws <- matrix(
  NA_real_,
  nrow = BOOTSTRAP_REPS,
  ncol = length(component_names),
  dimnames = list(NULL, component_names)
)
set.seed(BOOTSTRAP_SEED)
n_ids <- nrow(by_id)
for (rep in seq_len(BOOTSTRAP_REPS)) {
  frequency <- tabulate(
    sample.int(n_ids, n_ids, replace = TRUE),
    nbins = n_ids
  )
  bootstrap_draws[rep, "aggregate"] <- safe_ratio(
    sum(frequency * by_id$aggregate_num),
    sum(frequency * by_id$aggregate_den)
  )
  for (index in seq_len(nrow(AGE_BINS))) {
    label <- AGE_BINS$age_bin[index]
    stem <- gsub("-", "_", label)
    bootstrap_draws[rep, label] <- safe_ratio(
      sum(frequency * by_id[[paste0(stem, "_num")]]),
      sum(frequency * by_id[[paste0(stem, "_den")]])
    )
  }
}

aggregate_target <- cbind(
  data.table(
    moment = "aggregate_net_worth_to_gross_labor_earnings",
    vintage = sprintf("%d-%d", RECENT_YEAR_MIN, RECENT_YEAR_MAX),
    age_universe = sprintf("%d-%d", MODEL_AGE_MIN, MODEL_AGE_MAX),
    working_age_universe =
      sprintf("%d-%d", MODEL_AGE_MIN, MODEL_WORKING_AGE_MAX),
    denominator = "RP/SP combined gross labor earnings (EARNINDRRC)",
    family_years = aggregate_components$family_years,
    persons = aggregate_components$persons,
    aggregate_wealth = aggregate_components$numerator,
    aggregate_gross_labor_earnings = aggregate_components$denominator
  ),
  bootstrap_summary(
    aggregate_point,
    bootstrap_draws[, "aggregate"],
    BOOTSTRAP_REPS
  ),
  data.table(
    bootstrap_cluster = "reference-person ID",
    bootstrap_seed = BOOTSTRAP_SEED,
    status = "hard_target"
  )
)

age_profile <- merge(
  age_profile_point,
  rbindlist(lapply(AGE_BINS$age_bin, function(label) {
    cbind(
      data.table(age_bin = label),
      bootstrap_summary(
        age_profile_point[age_bin == label]$estimate,
        bootstrap_draws[, label],
        BOOTSTRAP_REPS
      )[, -"estimate"]
    )
  })),
  by = "age_bin",
  all.x = TRUE,
  sort = FALSE
)
age_profile[, `:=`(
  vintage = sprintf("%d-%d", RECENT_YEAR_MIN, RECENT_YEAR_MAX),
  denominator = "RP/SP combined gross labor earnings (EARNINDRRC)",
  bootstrap_cluster = "reference-person ID",
  bootstrap_seed = BOOTSTRAP_SEED,
  status = "robustness_only"
)]

definition_sensitivity <- data.table(
  definition = c(
    "aggregate net worth / gross labor earnings",
    "aggregate net worth / total family income"
  ),
  denominator_variable = c("EARNINDRRC", "INCFAMR"),
  estimate = c(
    aggregate_point,
    safe_ratio(
      aggregate_components$numerator,
      aggregate_components$family_income_denominator
    )
  ),
  status = c("hard_target", "definitional_sensitivity_only")
)

fwrite(yearly, file.path(out_dir, "yearly_ratios.csv"))
fwrite(aggregate_target, file.path(out_dir, "aggregate_target.csv"))
fwrite(age_profile, file.path(out_dir, "age_profile_robustness.csv"))
fwrite(
  definition_sensitivity,
  file.path(out_dir, "definition_sensitivity.csv")
)
fwrite(
  as.data.table(bootstrap_draws),
  file.path(out_dir, "bootstrap_draws.csv")
)

readme <- sprintf(
  paste0(
    "# Corrected PSID aggregate wealth/earnings target\n\n",
    "The hard target is aggregate family net worth divided by aggregate ",
    "RP/spouse gross labor earnings. Gross labor earnings are used instead ",
    "of total family income because they map directly to the model's labor-",
    "income process; `INCFAMR` additionally contains pensions, transfers, ",
    "asset income, and other components.\n\n",
    "- Vintage: %d--%d PSID waves.\n",
    "- Wealth universe: reference persons ages %d--%d.\n",
    "- Earnings universe: reference persons ages %d--%d.\n",
    "- Point target: `%.6f`.\n",
    "- Person-cluster bootstrap SE: `%.6f` ",
    "(%d draws; 95%% interval [`%.6f`, `%.6f`]).\n\n",
    "The model states are four-year periods: the working state labelled 62 ",
    "represents ages 62--65, and the terminal state labelled 82 represents ",
    "ages 82--85. Those intervals determine the empirical age restrictions.\n\n",
    "The 26--35 / 36--45 / 46--55 / 56--65 profile is robustness-only for ",
    "both the one-shot and sequential strands. It is not part of either ",
    "target system. The PSID remains weak in the extreme top wealth tail; ",
    "that limitation should be reported when interpreting this aggregate ",
    "ratio.\n"
  ),
  RECENT_YEAR_MIN,
  RECENT_YEAR_MAX,
  MODEL_AGE_MIN,
  MODEL_AGE_MAX,
  MODEL_AGE_MIN,
  MODEL_WORKING_AGE_MAX,
  aggregate_target$estimate,
  aggregate_target$bootstrap_se,
  BOOTSTRAP_REPS,
  aggregate_target$bootstrap_p025,
  aggregate_target$bootstrap_p975
)
writeLines(readme, file.path(out_dir, "README.md"))

print(aggregate_target)
print(age_profile)
