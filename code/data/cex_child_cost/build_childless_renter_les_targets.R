#!/usr/bin/env Rscript

# Estimate a childless-renter linear expenditure system (LES) from CEX PUMD.
# The primary regression is
#
#   rent = a + b * allocatable_expenditure + d * rent_per_room + error,
#
# which maps into the Stone-Geary parameters as
#   alpha_cons = 1 - b,
#   h_bar_0    = d / (1 - b),
#   c_bar_0    = -a / b.
#
# Rent per room is a deterministic two-fold cross-fitted geography mean.  Thus
# the price assigned to a household never uses interviews from that household's
# fold, avoiding a direct own-rent/own-rooms mechanical relationship.

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
  if (length(years) == 1L) paste0("les_", years) else
    paste0("les_", min(years), "_", max(years))
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
  "PERSLT18", "CUTENURE", "ROOMSQ", "REGION", "BLS_URBN", "POPSIZE",
  "PSU", "FINCBTAX", "FINCBTXM", "FINCBTXI", "INCLASS2", "FINLWT21",
  "TOTEXPCQ", "HOUSCQ", "SHELTCQ", "RNTXRPCQ",
  "RENDWECQ", "RNTAPYCQ", "CARTKNCQ", "CARTKUCQ", "OTHVEHCQ",
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

# Convert every monetary flow to 2023 dollars using the annual-average CPI-U.
# Interview packages span Q2 of the named year through Q1 of the next year from
# 2020 onward, so deflation follows QINTRVYR rather than the package label.
cpi_u <- c(
  `2019` = 255.657, `2020` = 258.811, `2021` = 270.970,
  `2022` = 292.655, `2023` = 304.702, `2024` = 313.689
)
if (any(!as.character(d$QINTRVYR) %in% names(cpi_u))) {
  stop("Missing CPI-U value for an interview year.")
}
deflator <- unname(cpi_u[["2023"]] / cpi_u[as.character(d$QINTRVYR)])
monetary <- c(
  "FINCBTAX", "FINCBTXM", "TOTEXPCQ", "HOUSCQ", "SHELTCQ", "RNTXRPCQ", "RENDWECQ",
  "RNTAPYCQ", "CARTKNCQ", "CARTKUCQ", "OTHVEHCQ", "CASHCOCQ", "PERINSCQ"
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

to_character <- function(x) {
  y <- as.character(x)
  y[is.na(x)] <- ""
  y
}

# Public CEX identifies selected metropolitan PSUs and suppresses PSU elsewhere.
# For suppressed observations, use the public region x population-size x urban
# cell.  The price cell is annual (not interview-quarter) to retain adequate
# opposite-fold cell size in this one-year pilot.
d$psu_chr <- to_character(d$PSU)
d$coarse_geo <- interaction(d$REGION, d$POPSIZE, d$BLS_URBN, drop = TRUE)
d$geo <- ifelse(nchar(d$psu_chr) > 0, paste0("PSU_", d$psu_chr),
                paste0("COARSE_", d$coarse_geo))

# Keep every interview from a consumer unit in the same fold.  CUID is numeric
# in the public file; using its final digit makes this reproducible.
cuid_digits <- gsub("[^0-9]", "", to_character(d$CUID))
d$fold <- as.integer(substr(cuid_digits, nchar(cuid_digits), nchar(cuid_digits))) %% 2L

# RNTXRPCQ is cash rent excluding rent-as-pay, measured over the quarter.
d$monthly_rent <- d$RNTXRPCQ / 3
d$monthly_rent_per_room <- d$monthly_rent / d$ROOMSQ

price_pool <- subset(
  d,
  CUTENURE == 4 & ROOMSQ > 0 & RNTXRPCQ > 0 & FINLWT21 > 0 &
    is.finite(monthly_rent_per_room)
)

# Trim only the price-building sample to prevent miscoded room/rent ratios from
# moving a small public geography cell.  The cross-fitted mean is used rather
# than a median because it maps directly to a local unit-rent average.
price_cut <- quantile(price_pool$monthly_rent_per_room, c(0.01, 0.99), na.rm = TRUE)
price_pool <- subset(
  price_pool,
  monthly_rent_per_room >= price_cut[[1]] &
    monthly_rent_per_room <= price_cut[[2]]
)

weighted_cell_mean <- function(data, group_vars, min_n = 10L) {
  key <- interaction(data[group_vars], drop = TRUE, lex.order = TRUE)
  pieces <- split(seq_len(nrow(data)), key)
  rows <- lapply(pieces, function(ii) {
    x <- data[ii, ]
    data.frame(
      cell_key = as.character(key[ii[[1]]]),
      cell_n = nrow(x),
      price = weighted.mean(x$monthly_rent_per_room, x$FINLWT21),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[out$cell_n >= min_n & is.finite(out$price), ]
}

make_lookup <- function(pool, group_vars, prefix, min_n) {
  z <- weighted_cell_mean(pool, group_vars, min_n)
  names(z)[names(z) == "price"] <- paste0(prefix, "_price")
  names(z)[names(z) == "cell_n"] <- paste0(prefix, "_n")
  z
}

# Build each fold's price from the opposite fold, then cascade to broader cells
# only when the detailed cell has fewer than ten contributing interviews.
assign_prices <- function(target_fold) {
  train <- price_pool[price_pool$fold != target_fold, ]
  target <- d[d$fold == target_fold, ]

  target$geo_key <- as.character(interaction(target[c("geo")], drop = TRUE,
                                               lex.order = TRUE))
  geo_lu <- make_lookup(train, c("geo"), "geo", 10L)
  target <- merge(target, geo_lu, by.x = "geo_key", by.y = "cell_key", all.x = TRUE,
                  sort = FALSE)

  target$coarse_key <- as.character(interaction(target[c("coarse_geo")], drop = TRUE,
                                                  lex.order = TRUE))
  coarse_lu <- make_lookup(train, c("coarse_geo"), "coarse", 10L)
  target <- merge(target, coarse_lu, by.x = "coarse_key", by.y = "cell_key", all.x = TRUE,
                  sort = FALSE)

  region_lu <- make_lookup(train, c("REGION"), "region", 20L)
  target$region_key <- as.character(interaction(target[c("REGION")], drop = TRUE,
                                                  lex.order = TRUE))
  target <- merge(target, region_lu, by.x = "region_key", by.y = "cell_key", all.x = TRUE,
                  sort = FALSE)

  national_price <- weighted.mean(train$monthly_rent_per_room, train$FINLWT21)
  target$price_source <- ifelse(!is.na(target$geo_price), "geo",
                         ifelse(!is.na(target$coarse_price), "coarse",
                           ifelse(!is.na(target$region_price), "region", "national")))
  target$monthly_price_per_room <- ifelse(
    !is.na(target$geo_price), target$geo_price,
    ifelse(!is.na(target$coarse_price), target$coarse_price,
      ifelse(!is.na(target$region_price), target$region_price, national_price))
  )
  target
}

d_cf <- rbind(assign_prices(0L), assign_prices(1L))

# For renters, subtract shelter and excluded asset/transfer items from total
# expenditure, then add cash rent back.  Utilities and household operations stay
# in the consumption composite; vehicle purchases, cash contributions, and
# personal insurance/pensions do not.
d_cf$consumption_q <- with(
  d_cf,
  TOTEXPCQ - SHELTCQ - CARTKNCQ - CARTKUCQ - OTHVEHCQ - CASHCOCQ - PERINSCQ
)
d_cf$rent_annual <- 4 * d_cf$RNTXRPCQ
d_cf$consumption_annual <- 4 * d_cf$consumption_q
d_cf$allocatable_annual <- d_cf$consumption_annual + d_cf$rent_annual
d_cf$annual_price_per_room <- 12 * d_cf$monthly_price_per_room

s <- subset(
  d_cf,
  CUTENURE == 4 & AGE_REF >= 30 & AGE_REF <= 55 & PERSLT18 == 0 &
    FAM_SIZE >= 1 & FAM_SIZE <= 2 & ROOMSQ > 0 & RNTXRPCQ > 0 &
    FINCBTAX > 1000 & FINLWT21 > 0 & consumption_annual > 0 &
    allocatable_annual > 0 & annual_price_per_room > 0
)

# Symmetric expenditure and rent trimming is fixed before estimating all
# specifications and is reported in the output.
s_pretrim <- s
trim_x <- quantile(s$allocatable_annual, c(0.01, 0.99), na.rm = TRUE)
trim_rent <- quantile(s$rent_annual, c(0.01, 0.99), na.rm = TRUE)
s <- subset(
  s,
  allocatable_annual >= trim_x[[1]] & allocatable_annual <= trim_x[[2]] &
    rent_annual >= trim_rent[[1]] & rent_annual <= trim_rent[[2]]
)
matched_trim_x <- c(
  weighted_quantile(s_pretrim$allocatable_annual, s_pretrim$FINLWT21, 0.01),
  weighted_quantile(s_pretrim$allocatable_annual, s_pretrim$FINLWT21, 0.99)
)
matched_trim_rent <- c(
  weighted_quantile(s_pretrim$rent_annual, s_pretrim$FINLWT21, 0.01),
  weighted_quantile(s_pretrim$rent_annual, s_pretrim$FINLWT21, 0.99)
)
s_one_market <- subset(
  s_pretrim,
  allocatable_annual >= matched_trim_x[[1]] &
    allocatable_annual <= matched_trim_x[[2]] &
    rent_annual >= matched_trim_rent[[1]] & rent_annual <= matched_trim_rent[[2]]
)

# Center demographic controls so the intercept in the controlled regression is
# evaluated at the weighted sample mean household rather than at age zero.
s$age_centered <- s$AGE_REF - weighted.mean(s$AGE_REF, s$FINLWT21)
s$two_adults_centered <- as.integer(s$FAM_SIZE == 2) -
  weighted.mean(as.integer(s$FAM_SIZE == 2), s$FINLWT21)

models <- list(
  primary = lm(
    rent_annual ~ allocatable_annual + annual_price_per_room,
    data = s, weights = FINLWT21
  ),
  demographic_controls = lm(
    rent_annual ~ allocatable_annual + annual_price_per_room +
      age_centered + I(age_centered^2) + two_adults_centered,
    data = s, weights = FINLWT21
  ),
  own_price_diagnostic = lm(
    rent_annual ~ allocatable_annual + I(12 * monthly_rent_per_room),
    data = s, weights = FINLWT21
  )
)

extract_model <- function(name, model) {
  vc <- vcovCL(model, cluster = ~CUID, type = "HC1")
  cf <- coef(model)
  se <- sqrt(diag(vc))
  price_term <- if (name == "own_price_diagnostic") {
    "I(12 * monthly_rent_per_room)"
  } else {
    "annual_price_per_room"
  }
  a <- unname(cf[["(Intercept)"]])
  b <- unname(cf[["allocatable_annual"]])
  d_price <- unname(cf[[price_term]])
  alpha_cons <- 1 - b
  h_bar_0 <- d_price / alpha_cons
  c_bar_0_annual_dollars <- -a / b
  mean_income <- weighted.mean(s$FINCBTAX, s$FINLWT21)
  g_h <- rep(0, length(cf)); names(g_h) <- names(cf)
  g_h[["allocatable_annual"]] <- d_price / alpha_cons^2
  g_h[[price_term]] <- 1 / alpha_cons
  h_bar_0_se <- sqrt(drop(t(g_h) %*% vc %*% g_h))
  g_c <- rep(0, length(cf)); names(g_c) <- names(cf)
  g_c[["(Intercept)"]] <- -1 / b
  g_c[["allocatable_annual"]] <- a / b^2
  c_bar_0_annual_se <- sqrt(drop(t(g_c) %*% vc %*% g_c))
  data.frame(
    specification = name,
    observations = nobs(model),
    unique_consumer_units = length(unique(s$CUID)),
    intercept = a,
    intercept_se = unname(se[["(Intercept)"]]),
    expenditure_coefficient = b,
    expenditure_coefficient_se = unname(se[["allocatable_annual"]]),
    price_coefficient = d_price,
    price_coefficient_se = unname(se[[price_term]]),
    alpha_cons = alpha_cons,
    alpha_cons_se = unname(se[["allocatable_annual"]]),
    h_bar_0_rooms = h_bar_0,
    h_bar_0_rooms_se = h_bar_0_se,
    c_bar_0_annual_dollars = c_bar_0_annual_dollars,
    c_bar_0_annual_dollars_se = c_bar_0_annual_se,
    c_bar_0_model_units = 4 * c_bar_0_annual_dollars / mean_income,
    c_bar_0_model_units_se = 4 * c_bar_0_annual_se / mean_income,
    admissible = alpha_cons > 0 & alpha_cons < 1 & h_bar_0 > 0 &
      c_bar_0_annual_dollars > 0,
    r_squared = summary(model)$r.squared,
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, Map(extract_model, names(models), models))

# Model-feasible one-market auxiliary statistics. The one-market model cannot
# reproduce the cross-geography price coefficient. We therefore retain the
# expenditure slope, evaluate the fitted intercept at the weighted mean unit
# rent, and replace the price coefficient with mean rooms in the bottom
# weighted quintile of allocatable expenditure. All three are observables that
# can be computed from model policies rather than parameter identities.
one_market_model <- lm(
  rent_annual ~ allocatable_annual + annual_price_per_room,
  data = s_one_market, weights = FINLWT21
)
primary_cf <- coef(one_market_model)
mean_income_primary <- weighted.mean(s_one_market$FINCBTAX, s_one_market$FINLWT21)
mean_price_primary <- weighted.mean(
  s_one_market$annual_price_per_room, s_one_market$FINLWT21
)
low_allocatable_cutoff <- weighted_quantile(
  s_one_market$allocatable_annual, s_one_market$FINLWT21, 0.20
)
low_allocatable_sample <- s_one_market[
  s_one_market$allocatable_annual <= low_allocatable_cutoff,
]
one_market_auxiliary <- data.frame(
  moment = c(
    "childless_renter_rent_expenditure_slope",
    "childless_renter_intercept_at_mean_price_model_units",
    "bottom_quintile_childless_renter_mean_rooms"
  ),
  value = c(
    unname(primary_cf[["allocatable_annual"]]),
    4 * (
      unname(primary_cf[["(Intercept)"]]) +
        unname(primary_cf[["annual_price_per_room"]]) * mean_price_primary
    ) / mean_income_primary,
    weighted.mean(low_allocatable_sample$ROOMSQ, low_allocatable_sample$FINLWT21)
  ),
  sample = c(
    "childless cash renters ages 30-55; weighted-trim matched sample",
    "childless cash renters ages 30-55; weighted-trim matched sample",
    "bottom weighted allocatable-expenditure quintile of matched sample"
  ),
  observations = c(nobs(one_market_model), nobs(one_market_model),
                   nrow(low_allocatable_sample)),
  stringsAsFactors = FALSE
)

estimate_sensitivity <- function(data, name, age_lo = 30, age_hi = 55,
                                 family_rule = "one_two", consumption_q = NULL) {
  z <- data
  if (is.null(consumption_q)) consumption_q <- z$consumption_q
  z$consumption_alt_annual <- 4 * consumption_q
  z$rent_alt_annual <- 4 * z$RNTXRPCQ
  z$allocatable_alt_annual <- z$consumption_alt_annual + z$rent_alt_annual
  z$annual_price_alt <- 12 * z$monthly_price_per_room
  valid_family <- if (family_rule == "single") z$FAM_SIZE == 1 else
    z$FAM_SIZE >= 1 & z$FAM_SIZE <= 2
  valid <- complete.cases(z[c(
    "CUTENURE", "AGE_REF", "PERSLT18", "FAM_SIZE", "ROOMSQ", "RNTXRPCQ",
    "FINCBTAX", "FINLWT21", "consumption_alt_annual",
    "allocatable_alt_annual", "annual_price_alt"
  )]) &
    z$CUTENURE == 4 & z$AGE_REF >= age_lo & z$AGE_REF <= age_hi &
      z$PERSLT18 == 0 & valid_family & z$ROOMSQ > 0 & z$RNTXRPCQ > 0 &
      z$FINCBTAX > 1000 & z$FINLWT21 > 0 & z$consumption_alt_annual > 0 &
      z$allocatable_alt_annual > 0 & z$annual_price_alt > 0
  q <- z[valid, ]
  tx <- quantile(q$allocatable_alt_annual, c(0.01, 0.99), na.rm = TRUE)
  tr <- quantile(q$rent_alt_annual, c(0.01, 0.99), na.rm = TRUE)
  q <- q[
    q$allocatable_alt_annual >= tx[[1]] & q$allocatable_alt_annual <= tx[[2]] &
      q$rent_alt_annual >= tr[[1]] & q$rent_alt_annual <= tr[[2]],
  ]
  m <- lm(rent_alt_annual ~ allocatable_alt_annual + annual_price_alt,
          data = q, weights = FINLWT21)
  cf <- coef(m); a <- cf[["(Intercept)"]]; b <- cf[["allocatable_alt_annual"]]
  dp <- cf[["annual_price_alt"]]; alpha_alt <- 1 - b
  data.frame(
    case = name,
    observations = nrow(q),
    unique_consumer_units = length(unique(q$CUID)),
    alpha_cons = alpha_alt,
    h_bar_0_rooms = dp / alpha_alt,
    c_bar_0_model_units = 4 * (-a / b) / weighted.mean(q$FINCBTAX, q$FINLWT21),
    r_squared = summary(m)$r.squared,
    stringsAsFactors = FALSE
  )
}

sample_sensitivities <- rbind(
  estimate_sensitivity(d_cf, "age30_55_one_two"),
  estimate_sensitivity(d_cf, "age25_55_one_two", age_lo = 25),
  estimate_sensitivity(d_cf, "age30_55_single", family_rule = "single"),
  estimate_sensitivity(d_cf, "age25_55_single", age_lo = 25, family_rule = "single")
)

excluded_q <- with(d_cf, CARTKNCQ + CARTKUCQ + OTHVEHCQ + CASHCOCQ + PERINSCQ)
composite_sensitivities <- rbind(
  estimate_sensitivity(d_cf, "exclude_shelter",
                       consumption_q = with(d_cf, TOTEXPCQ - SHELTCQ - excluded_q)),
  estimate_sensitivity(d_cf, "exclude_all_housing",
                       consumption_q = with(d_cf, TOTEXPCQ - HOUSCQ - excluded_q)),
  estimate_sensitivity(d_cf, "include_all_except_cash_rent",
                       consumption_q = with(d_cf, TOTEXPCQ - RNTXRPCQ - excluded_q))
)

sample_summary <- data.frame(
  statistic = c(
    "interview_observations", "unique_consumer_units", "weighted_mean_income",
    "weighted_mean_allocatable_expenditure", "weighted_mean_annual_rent",
    "weighted_mean_rooms", "weighted_mean_price_per_room", "p01_allocatable",
    "p99_allocatable", "p01_annual_rent", "p99_annual_rent",
    "p01_price_building_rent_per_room", "p99_price_building_rent_per_room"
  ),
  value = c(
    nrow(s), length(unique(s$CUID)), weighted.mean(s$FINCBTAX, s$FINLWT21),
    weighted.mean(s$allocatable_annual, s$FINLWT21),
    weighted.mean(s$rent_annual, s$FINLWT21),
    weighted.mean(s$ROOMSQ, s$FINLWT21),
    weighted.mean(s$annual_price_per_room, s$FINLWT21),
    trim_x[[1]], trim_x[[2]], trim_rent[[1]], trim_rent[[2]],
    price_cut[[1]], price_cut[[2]]
  )
)

price_sources <- aggregate(
  FINLWT21 ~ price_source,
  data = s,
  FUN = sum
)
price_sources$weighted_share <- price_sources$FINLWT21 / sum(price_sources$FINLWT21)
price_sources$FINLWT21 <- NULL

price_cells <- aggregate(
  cbind(monthly_price_per_room, FINLWT21) ~ geo + fold,
  data = s,
  FUN = function(x) c(mean = mean(x), n = length(x))
)

write.csv(results, file.path(output_dir, "les_parameter_estimates.csv"), row.names = FALSE)
write.csv(sample_summary, file.path(output_dir, "les_sample_summary.csv"), row.names = FALSE)
write.csv(price_sources, file.path(output_dir, "les_price_sources.csv"), row.names = FALSE)
write.csv(price_cells, file.path(output_dir, "les_price_cells.csv"), row.names = FALSE)
write.csv(sample_sensitivities,
          file.path(output_dir, "les_sample_sensitivities.csv"), row.names = FALSE)
write.csv(composite_sensitivities,
          file.path(output_dir, "les_composite_sensitivities.csv"), row.names = FALSE)
write.csv(one_market_auxiliary,
          file.path(output_dir, "one_market_auxiliary_targets.csv"), row.names = FALSE)

readme <- c(
  "# CEX childless-renter LES pilot",
  "",
  paste0("Source: ", min(years), "--", max(years),
         " Consumer Expenditure Survey Interview PUMD; monetary values are 2023 dollars."),
  paste0("Income measure: ", income_measure,
         if (complete_reporters_only) "; complete income reporters only." else
           "; all reporters using BLS collected-or-imputed income."),
  "Sample: cash renters, reference person age 30--55, no persons under 18, one or two",
  "household members, positive income/consumption/rent, and valid rooms. The final",
  "sample trims allocatable expenditure and annual cash rent at the 1st/99th percentiles.",
  "Regressions use FINLWT21 and cluster standard errors by consumer unit.",
  "",
  "The primary local unit-rent measure is a deterministic two-fold cross-fitted mean",
  "cash rent per room. Identified metropolitan PSUs are used when available; suppressed",
  "PSUs use region x population-size x urban cells. Cells with fewer than ten opposite-fold",
  "interviews cascade to broader public geography. No household's own fold enters its price.",
  "",
  "The exact LES mapping is alpha_cons = 1-b, h_bar_0 = d/(1-b), and c_bar_0 = -a/b.",
  "c_bar_0_model_units converts the annual-dollar subsistence consumption estimate into",
  "four-year model units using weighted mean annual before-tax income.",
  "For the one-market model, one_market_auxiliary_targets.csv replaces the",
  "unavailable cross-sectional price-coefficient moment with mean rooms in the",
  "bottom weighted quintile of allocatable expenditure. It also reports the",
  "rent-expenditure slope and fitted intercept at the weighted mean unit rent.",
  "The matched auxiliary applies weighted 1st/99th percentile trims so the same",
  "deterministic operation can be applied to the model distribution;",
  "these must be recomputed from model policies rather than copied from parameters.",
  "",
  "This is a pooled-year feasibility estimate, not yet the final target. A credible final",
  "estimate should pool years, deflate expenditures, test alternative consumption composites,",
  "and verify stability of the price coefficient. The own-price specification is explicitly",
  "mechanical and is diagnostic only."
)
writeLines(readme, file.path(output_dir, "README.md"))

print(results)
print(sample_summary)
print(price_sources)
