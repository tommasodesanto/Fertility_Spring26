# Local synthetic-fixture tests for twins_samesex_iv_lib.R. Exercises key
# uniqueness, MOMLOC/RELATE link categories, sex availability, age-tie/order
# flags, NCHILD>=2/3 support, outcome coding, weight positivity, one row per
# maternal household anchor, and a manual RF/FS/2SLS sanity check against a
# known synthetic DGP. Run: Rscript test_twins_samesex_iv.R

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.script_path) == 1) dirname(normalizePath(.script_path)) else getwd()
source(file.path(.script_dir, "twins_samesex_iv_lib.R"))
suppressMessages(library(data.table))

stopifnot_msg <- function(cond, msg) if (!isTRUE(cond)) stop(paste("FAIL:", msg)) else cat("PASS:", msg, "\n")

set.seed(42)

## ---- Fixture 1: hand-built roster covering edge cases ----------------
# Households:
#  H1: mother + twin-like pair at age 0 (SS: both boys) -> Twin1 pos, SameSex2 same-sex
#  H2: mother + single child age 2 -> Twin1 eligible, not positive; SameSex2 ineligible (NCHILD<2)
#  H3: mother + 2 children ages 3,1 (diff sex) -> SameSex2 eligible, not same-sex
#  H4: mother + 3 children ages 4,2,2 (age tie at 2nd/3rd, oldest two ages 4,2 distinct) -> SameSex2 eligible on oldest two; contamination flag from 3rd tie only if a1==a3
#  H5: mother + child with missing sex -> excluded from SameSex2 primary
#  H6: mother + 2 children tied at oldest age (a1==a2) -> Twin1 positive AND SameSex2 primary_age_tie exclusion
#  H7: non-mother-linked child (MOMLOC points nowhere valid) -> unmatched, excluded
#  H8: self-link (MOMLOC == own PERNUM) -> excluded
#  H9: mother age 30, child age 35 (older than mother) -> invalid link excluded
#  H10: zero/negative weight mother -> flagged but not dropped by roster builder itself

rows <- rbindlist(list(
  data.table(YEAR=2015, SAMPLE=1, SERIAL=1, PERNUM=1, SEX=2, AGE=30, RELATE=1, MOMLOC=0, NCHILD=2, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=100, HHWT=100, STATEFIP=9, RACE=1, EDUC=6, MARST=1, FERTYR=2),
  data.table(YEAR=2015, SAMPLE=1, SERIAL=1, PERNUM=2, SEX=1, AGE=0, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=100, HHWT=100, STATEFIP=9, RACE=1, EDUC=NA, MARST=NA, FERTYR=NA),
  data.table(YEAR=2015, SAMPLE=1, SERIAL=1, PERNUM=3, SEX=1, AGE=0, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=100, HHWT=100, STATEFIP=9, RACE=1, EDUC=NA, MARST=NA, FERTYR=NA),

  data.table(YEAR=2016, SAMPLE=1, SERIAL=2, PERNUM=1, SEX=2, AGE=28, RELATE=1, MOMLOC=0, NCHILD=1, OWNERSHP=2, ROOMS=4, BEDROOMS=2, PERWT=80, HHWT=80, STATEFIP=9, RACE=2, EDUC=5, MARST=1, FERTYR=1),
  data.table(YEAR=2016, SAMPLE=1, SERIAL=2, PERNUM=2, SEX=2, AGE=2, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=2, ROOMS=4, BEDROOMS=2, PERWT=80, HHWT=80, STATEFIP=9, RACE=2, EDUC=NA, MARST=NA, FERTYR=NA),

  data.table(YEAR=2017, SAMPLE=1, SERIAL=3, PERNUM=1, SEX=2, AGE=33, RELATE=1, MOMLOC=0, NCHILD=2, OWNERSHP=1, ROOMS=8, BEDROOMS=4, PERWT=90, HHWT=90, STATEFIP=25, RACE=1, EDUC=8, MARST=1, FERTYR=2),
  data.table(YEAR=2017, SAMPLE=1, SERIAL=3, PERNUM=2, SEX=1, AGE=3, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=1, ROOMS=8, BEDROOMS=4, PERWT=90, HHWT=90, STATEFIP=25, RACE=1, EDUC=NA, MARST=NA, FERTYR=NA),
  data.table(YEAR=2017, SAMPLE=1, SERIAL=3, PERNUM=3, SEX=2, AGE=1, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=1, ROOMS=8, BEDROOMS=4, PERWT=90, HHWT=90, STATEFIP=25, RACE=1, EDUC=NA, MARST=NA, FERTYR=NA),

  data.table(YEAR=2018, SAMPLE=1, SERIAL=4, PERNUM=1, SEX=2, AGE=35, RELATE=1, MOMLOC=0, NCHILD=3, OWNERSHP=2, ROOMS=5, BEDROOMS=2, PERWT=70, HHWT=70, STATEFIP=25, RACE=3, EDUC=4, MARST=1, FERTYR=2),
  data.table(YEAR=2018, SAMPLE=1, SERIAL=4, PERNUM=2, SEX=1, AGE=4, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=2, ROOMS=5, BEDROOMS=2, PERWT=70, HHWT=70, STATEFIP=25, RACE=3, EDUC=NA, MARST=NA, FERTYR=NA),
  data.table(YEAR=2018, SAMPLE=1, SERIAL=4, PERNUM=3, SEX=2, AGE=2, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=2, ROOMS=5, BEDROOMS=2, PERWT=70, HHWT=70, STATEFIP=25, RACE=3, EDUC=NA, MARST=NA, FERTYR=NA),
  data.table(YEAR=2018, SAMPLE=1, SERIAL=4, PERNUM=4, SEX=1, AGE=2, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=2, ROOMS=5, BEDROOMS=2, PERWT=70, HHWT=70, STATEFIP=25, RACE=3, EDUC=NA, MARST=NA, FERTYR=NA),

  data.table(YEAR=2019, SAMPLE=1, SERIAL=5, PERNUM=1, SEX=2, AGE=27, RELATE=1, MOMLOC=0, NCHILD=1, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=60, HHWT=60, STATEFIP=6, RACE=1, EDUC=7, MARST=1, FERTYR=2),
  data.table(YEAR=2019, SAMPLE=1, SERIAL=5, PERNUM=2, SEX=NA, AGE=1, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=60, HHWT=60, STATEFIP=6, RACE=1, EDUC=NA, MARST=NA, FERTYR=NA),

  data.table(YEAR=2015, SAMPLE=1, SERIAL=6, PERNUM=1, SEX=2, AGE=31, RELATE=1, MOMLOC=0, NCHILD=2, OWNERSHP=2, ROOMS=3, BEDROOMS=1, PERWT=110, HHWT=110, STATEFIP=6, RACE=2, EDUC=6, MARST=1, FERTYR=2),
  data.table(YEAR=2015, SAMPLE=1, SERIAL=6, PERNUM=2, SEX=1, AGE=0, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=2, ROOMS=3, BEDROOMS=1, PERWT=110, HHWT=110, STATEFIP=6, RACE=2, EDUC=NA, MARST=NA, FERTYR=NA),
  data.table(YEAR=2015, SAMPLE=1, SERIAL=6, PERNUM=3, SEX=2, AGE=0, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=2, ROOMS=3, BEDROOMS=1, PERWT=110, HHWT=110, STATEFIP=6, RACE=2, EDUC=NA, MARST=NA, FERTYR=NA),

  data.table(YEAR=2016, SAMPLE=1, SERIAL=7, PERNUM=1, SEX=2, AGE=40, RELATE=1, MOMLOC=0, NCHILD=0, OWNERSHP=1, ROOMS=5, BEDROOMS=2, PERWT=50, HHWT=50, STATEFIP=6, RACE=1, EDUC=6, MARST=1, FERTYR=2),
  data.table(YEAR=2016, SAMPLE=1, SERIAL=7, PERNUM=2, SEX=1, AGE=5, RELATE=3, MOMLOC=9, NCHILD=NA, OWNERSHP=1, ROOMS=5, BEDROOMS=2, PERWT=50, HHWT=50, STATEFIP=6, RACE=1, EDUC=NA, MARST=NA, FERTYR=NA),

  data.table(YEAR=2017, SAMPLE=1, SERIAL=8, PERNUM=1, SEX=2, AGE=25, RELATE=1, MOMLOC=1, NCHILD=0, OWNERSHP=2, ROOMS=4, BEDROOMS=2, PERWT=45, HHWT=45, STATEFIP=6, RACE=1, EDUC=5, MARST=1, FERTYR=1),

  data.table(YEAR=2018, SAMPLE=1, SERIAL=9, PERNUM=1, SEX=2, AGE=30, RELATE=1, MOMLOC=0, NCHILD=1, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=65, HHWT=65, STATEFIP=6, RACE=1, EDUC=6, MARST=1, FERTYR=2),
  data.table(YEAR=2018, SAMPLE=1, SERIAL=9, PERNUM=2, SEX=1, AGE=35, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=65, HHWT=65, STATEFIP=6, RACE=1, EDUC=NA, MARST=NA, FERTYR=NA),

  data.table(YEAR=2019, SAMPLE=1, SERIAL=10, PERNUM=1, SEX=2, AGE=32, RELATE=1, MOMLOC=0, NCHILD=1, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=0, HHWT=0, STATEFIP=6, RACE=1, EDUC=6, MARST=1, FERTYR=2),
  data.table(YEAR=2019, SAMPLE=1, SERIAL=10, PERNUM=2, SEX=1, AGE=2, RELATE=3, MOMLOC=1, NCHILD=NA, OWNERSHP=1, ROOMS=6, BEDROOMS=3, PERWT=0, HHWT=0, STATEFIP=6, RACE=1, EDUC=NA, MARST=NA, FERTYR=NA)
), fill = TRUE)

## ---- Key uniqueness ----------------------------------------------------
stopifnot_msg(anyDuplicated(rows[, .(YEAR, SAMPLE, SERIAL, PERNUM)]) == 0L,
              "person key (YEAR,SAMPLE,SERIAL,PERNUM) unique in fixture")

built <- build_mother_roster(rows)
mr <- built$mother_rows
stopifnot_msg(anyDuplicated(mr$person_key) == 0L, "roster: one row per maternal household anchor")

## ---- Link-category audit ------------------------------------------------
stopifnot_msg(any(built$link_audit$link_invalid_reason == "child_without_female_mother", na.rm = TRUE),
              "unmatched-child link (H7) captured in audit")
h8 <- mr[household_key == "2017_1_8"]
stopifnot_msg(nrow(h8) == 1 && h8$linked_child_count == 0,
              "self-link household (H8) yields zero linked children")
h9 <- mr[household_key == "2018_1_9"]
stopifnot_msg(nrow(h9) == 1 && h9$linked_child_count == 0,
              "older-than-mother link (H9) excluded")

## ---- Sex availability / missingness -------------------------------------
h5 <- mr[household_key == "2019_1_5"]
stopifnot_msg(h5$linked_child_count == 1 && isTRUE(h5$any_sex_missing),
              "missing-sex child (H5) linked but flagged any_sex_missing")

## ---- Weight positivity ---------------------------------------------------
stopifnot_msg(any(rows$PERWT <= 0), "fixture includes a zero-weight mother (H10) for the weight-positivity check")
mr[, weight_positive := PERWT > 0 & is.finite(PERWT)]
stopifnot_msg(!mr[household_key == "2019_1_10", weight_positive],
              "zero-weight mother (H10) flagged weight_invalid / not weight_positive")

mr <- add_outcomes(mr)

## ---- Twin1 construction ---------------------------------------------------
t1 <- build_twin1(mr, age_grid = 0:5)
h1 <- t1[household_key == "2015_1_1"]
stopifnot_msg(nrow(h1) == 1 && isTRUE(h1$twin_like_proxy) && h1$event_age == 0,
              "Twin1: H1 (two age-0 linked children) is twin-like-positive at event age 0")
h2 <- t1[household_key == "2016_1_2"]
stopifnot_msg(nrow(h2) == 1 && !isTRUE(h2$twin_like_proxy) && !h2$treatment_2plus,
              "Twin1: H2 (single child) eligible but not twin-positive, treatment=0")
h6 <- t1[household_key == "2015_1_6"]
stopifnot_msg(nrow(h6) == 1 && isTRUE(h6$twin_like_proxy),
              "Twin1: H6 (age-tied oldest two) is twin-like-positive")
stopifnot_msg(all(t1$oldest_age %in% 0:5), "Twin1: age grid restricted to 0:5")

## ---- SameSex2 construction --------------------------------------------
ss <- build_samesex2(mr, age_grid = 0:5)
h3 <- ss[household_key == "2017_1_3"]
stopifnot_msg(nrow(h3) == 1 && isTRUE(h3$eligible) && h3$samesex == 0,
              "SameSex2: H3 (boy age3, girl age1) eligible, opposite-sex")
h4 <- ss[household_key == "2018_1_4"]
stopifnot_msg(nrow(h4) == 1 && isTRUE(h4$eligible) && h4$samesex == 0 && !h4$primary_age_tie,
              "SameSex2: H4 oldest-two (ages 4,2) not tied even though 2nd/3rd tie at age 2")
stopifnot_msg(isTRUE(h4$treatment_3plus), "SameSex2: H4 has 3 linked children -> treatment_3plus TRUE")
h6ss <- ss[household_key == "2015_1_6"]
stopifnot_msg(nrow(h6ss) == 1 && isTRUE(h6ss$primary_age_tie) && !isTRUE(h6ss$eligible),
              "SameSex2: H6 (oldest-two age-tied) excluded from primary group via primary_age_tie")
h2ss <- ss[household_key == "2016_1_2"]
stopifnot_msg(nrow(h2ss) == 0,
              "SameSex2: H2 (only 1 linked child) dropped, not in eligible_pool (NCHILD>=2 required)")
h5ss <- ss[household_key == "2019_1_5"]
stopifnot_msg(nrow(h5ss) == 0,
              "SameSex2: H5 (only 1 linked child, sex missing) dropped, not in eligible_pool")

## ---- Outcome coding -----------------------------------------------------
stopifnot_msg(all(mr[!is.na(ROOMS_out), ROOMS_out] <= rooms_cap), "ROOMS capped at 9")
stopifnot_msg(all(mr[!is.na(OWNERSHP_out), OWNERSHP_out] %in% c(0, 1)), "OWNERSHP recoded to 0/1")

cat("\n--- Fixture 1 structural checks: ALL PASS ---\n\n")

## ---- Fixture 2: synthetic DGP for RF/FS/2SLS/AR numerical sanity --------
# beta_true relates D (treatment) to Y (outcome) with a known instrument Z
# that shifts D with first-stage pi=0.5 and has no direct effect on Y
# (valid-IV synthetic case, purely to check the estimator code, NOT a claim
# about the real Twin1/SameSex2 exclusion restriction).
n <- 4000
set.seed(7)
mat_age <- round(runif(n, 25, 40))
race <- sample(1:3, n, replace = TRUE)
Z <- rbinom(n, 1, 0.15)
p_d <- pmin(pmax(0.3 + 0.5 * Z + 0.005 * (mat_age - 30), 0.01), 0.99)
D <- rbinom(n, 1, p_d)
beta_true <- 0.8
Y <- 5 + beta_true * D + 0.03 * (mat_age - 30) + rnorm(n, sd = 1.2)
sim <- data.table(Y = Y, D = D, Z = Z, mat_age = mat_age, race = factor(race),
                   household_key = paste0("hh", seq_len(n)),
                   mother_weight = 1)

controls_fml <- "mat_age + i(race)"
fit <- fit_instrument_outcome(sim, outcome = "Y", treatment = "D", instrument = "Z",
                               controls_fml = controls_fml, weight_var = "mother_weight",
                               cluster_var = "household_key",
                               ar_grid = seq(0, 2, by = 0.05))

cat(sprintf("Synthetic DGP: true beta=%.3f | RF=%.3f FS=%.3f IV=%.3f (SE=%.3f) FS-F=%.1f\n",
            beta_true, fit$rf_coef, fit$fs_coef, fit$iv_coef, fit$iv_se, fit$first_stage_F))

stopifnot_msg(fit$status == "fit", "synthetic DGP: estimator reaches status=fit")
stopifnot_msg(abs(fit$fs_coef - 0.5) < 0.1, "synthetic DGP: first-stage coef recovers ~0.5 within tolerance")
stopifnot_msg(fit$first_stage_F > 10, "synthetic DGP: first-stage F is strong (>10) as designed")
stopifnot_msg(abs(fit$iv_coef - beta_true) < 0.35,
              "synthetic DGP: 2SLS coefficient recovers true beta within tolerance")
stopifnot_msg(!is.na(fit$ar_lower) && !is.na(fit$ar_upper) && fit$ar_bounded,
              "synthetic DGP: AR confidence set is bounded on a strong instrument")
stopifnot_msg(fit$ar_lower <= beta_true && beta_true <= fit$ar_upper,
              "synthetic DGP: AR set covers the true beta")

## ---- Weak-instrument case: AR should be unbounded/wide, not fabricated --
Zw <- rbinom(n, 1, 0.15)
Dw <- rbinom(n, 1, pmin(pmax(0.3 + 0.01 * Zw + 0.005 * (mat_age - 30), 0.01), 0.99))  # ~no first stage
simw <- data.table(Y = Y, D = Dw, Z = Zw, mat_age = mat_age, race = factor(race),
                    household_key = paste0("hh", seq_len(n)), mother_weight = 1)
fitw <- fit_instrument_outcome(simw, outcome = "Y", treatment = "D", instrument = "Z",
                                controls_fml = controls_fml, weight_var = "mother_weight",
                                cluster_var = "household_key",
                                ar_grid = seq(-2, 4, by = 0.1))
cat(sprintf("Weak-IV DGP: FS=%.3f FS-F=%.2f AR bounded=%s [%.2f, %.2f] (grid %.1f..%.1f)\n",
            fitw$fs_coef, fitw$first_stage_F, fitw$ar_bounded, fitw$ar_lower, fitw$ar_upper,
            min(seq(-2,4,by=0.1)), max(seq(-2,4,by=0.1))))
stopifnot_msg(fitw$first_stage_F < 10, "weak-IV DGP: first-stage F is weak as designed")

cat("\n--- Fixture 2 numerical sanity checks: ALL PASS ---\n\n")
cat("ALL TESTS PASSED\n")
