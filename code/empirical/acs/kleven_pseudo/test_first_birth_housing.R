#!/usr/bin/env Rscript
# Tiny deterministic contract tests. No ACS/PSID files are loaded.
source("estimate_first_birth_housing.R", local = TRUE)

expect <- function(ok, msg) if (!isTRUE(ok)) stop(paste("FAIL:", msg), call. = FALSE)
expect_error <- function(expr, msg) {
  hit <- FALSE
  tryCatch(force(expr), error = function(e) hit <<- TRUE)
  expect(hit, msg)
}

key_cols <- c("YEAR", "SAMPLE", "SERIAL", "PERNUM")
events <- c("-2", "-1", "0", "1", "2", "3")

# Explicit post-audit coding object. ROOMS=28 is retained as an unknown flag
# and coded missing; ROOMS=0 is a configured missing code. The bedrooms
# transformation follows the extract27 coding: code 1 means zero bedrooms,
# codes 1:22 map to 0:21; primary bedrooms are capped at 5.
coding <- list(
  rooms_valid = function(x, year) !is.na(x) & x %in% c(1:27, 30),
  rooms_transform = function(x, year) x,
  rooms_cap = 9,
  rooms_missing_codes = 0,
  rooms_unknown_codes = 28,
  bedrooms_valid = function(x, year) !is.na(x) & x %in% 1:22,
  bedrooms_transform = function(x, year) ifelse(x == 22, 21, x - 1),
  bedrooms_cap = 5,
  bedrooms_missing_codes = 0,
  bedrooms_unknown_codes = 23,
  ownership_valid = function(x, year) !is.na(x) & x %in% c(1, 2),
  ownership_missing_codes = 0,
  ownership_unknown_codes = c(3, 9),
  allow_uncapped_sensitivity = FALSE
)

audit_manifest <- list(status = "PASS", verified = TRUE, key_columns = key_cols,
  source_key_unique = TRUE, overlap_verified = TRUE,
  source_packet = "verified NE housing key packet supplied by launcher")

# Two states x two genders x six event times x four age cells x four survey
# cells x two persons per source household. The second person repeats the
# source household cluster while retaining a unique source key through PERNUM.
g <- expand.grid(statename = c("Maine", "Vermont"), gender = c("Men", "Women"),
                 event_time = events, age_factor = paste0("a", 1:4),
                 doiy_factor = paste0("y", 1:4), person = 1:2,
                 KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
g$YEAR <- 2005L + as.integer(sub("y", "", g$doiy_factor)) - 1L
g$SAMPLE <- 200501L
g$SERIAL <- as.integer(rep(seq_len(nrow(g) / 2), times = 2))
g$PERNUM <- g$person
g$statefip <- ifelse(g$statename == "Maine", 23L, 50L)
g$census <- "New England"
g$src_key <- paste("ACS", g$YEAR, g$SAMPLE, g$SERIAL, g$PERNUM, sep = ":")
g$source_origin <- "ACS"
g$from_cps <- 0L
g$wgt <- 10 + (seq_len(nrow(g)) %% 7)
g$emp_lw <- as.numeric((seq_len(nrow(g)) %% 3) > 0)
g$cohort <- ifelse(g$statename == "Maine", 2010L, 2011L)
g$age_factor <- factor(g$age_factor)
g$doiy_factor <- factor(g$doiy_factor)
g$t_es_lw <- g$event_time
g$rooms_raw_expected <- 4 + as.integer(g$event_time) + (g$statename == "Vermont")
g$rooms_raw_expected[g$rooms_raw_expected < 1] <- 1
g$bedrooms_raw_expected <- ifelse(as.integer(g$event_time) >= 0, 4, 3)
g$ownership_raw_expected <- ifelse(g$statename == "Vermont" & as.integer(g$event_time) >= 0, 1,
                                  ifelse(g$statename == "Vermont", 2, 1))

source_housing <- g[, key_cols, drop = FALSE]
source_housing$ROOMS <- g$rooms_raw_expected
source_housing$BEDROOMS <- g$bedrooms_raw_expected
source_housing$OWNERSHP <- g$ownership_raw_expected

panel <- g[, c(key_cols, "src_key", "source_origin", "from_cps", "wgt", "emp_lw",
               "t_es_lw", "cohort", "census", "statefip", "statename", "gender",
               "age_factor", "doiy_factor"), drop = FALSE]
# Add a relabeled CPS donor and an original CPS row with no source key: neither
# has ACS housing data, but both provenance classes must survive the join.
cps <- panel[1, , drop = FALSE]
cps[key_cols] <- list(NA_integer_, NA_integer_, NA_integer_, NA_integer_)
cps$src_key <- "CPS:donor:1"; cps$source_origin <- "CPS"; cps$from_cps <- 1L
cps$wgt <- 17; cps$emp_lw <- 1; cps$t_es_lw <- "-2"; cps$cohort <- 2010L
original_cps <- panel[2, , drop = FALSE]
original_cps[key_cols] <- list(NA_integer_, NA_integer_, NA_integer_, NA_integer_)
original_cps$src_key <- "CPS:original:2"; original_cps$source_origin <- "CPS"; original_cps$from_cps <- 0L
original_cps$wgt <- 18; original_cps$emp_lw <- 1; original_cps$t_es_lw <- "-2"; original_cps$cohort <- 2010L
acs_2000_unavailable <- panel[3, , drop = FALSE]
acs_2000_unavailable[key_cols] <- list(2000L, 200001L, 999999L, 1L)
acs_2000_unavailable$src_key <- "2000:200001:999999:1"; acs_2000_unavailable$source_origin <- "ACS"; acs_2000_unavailable$from_cps <- 0L
acs_2000_unavailable$wgt <- 19; acs_2000_unavailable$emp_lw <- 1; acs_2000_unavailable$t_es_lw <- "-2"; acs_2000_unavailable$cohort <- 2000L
panel <- rbind(panel, cps, original_cps, acs_2000_unavailable)

joined <- join_first_birth_housing(panel, source_housing, audit_manifest)
expect(nrow(joined) == nrow(panel), "source join changes row count")
expect(identical(joined$first_birth_row_id, seq_len(nrow(panel))), "source join changes row order")
expect(isTRUE(all.equal(joined$wgt, panel$wgt, check.attributes = FALSE)), "weights changed by join")
expect(isTRUE(all.equal(joined$t_es_lw, panel$t_es_lw, check.attributes = FALSE)), "event time changed by join")
expect(isTRUE(all.equal(joined$emp_lw, panel$emp_lw, check.attributes = FALSE)), "labor outcome changed by join")
nn <- nrow(joined)
expect(joined$housing_join_status[nn - 2L] == "relabeled_cps_missing_outcome", "relabeled CPS row not retained as missing")
expect(joined$housing_join_status[nn - 1L] == "original_cps_missing_outcome", "original CPS row not retained as missing")
expect(joined$housing_join_status[nn] == "acs_source_unmatched", "unavailable 2000 ACS row not audited")
expect(all(is.na(joined$raw_rooms[(nn - 2L):nn])), "unavailable-source housing was imputed")
join_audit <- attr(joined, "housing_join_audit")
expect(join_audit$relabeled_cps_missing_outcome_rows == 1L && join_audit$original_cps_missing_outcome_rows == 1L,
       "CPS provenance audit counts failed")
expect(join_audit$repeated_source_household_clusters > 0,
       "repeated source household cluster was not recorded")

expect_error(join_first_birth_housing(panel, rbind(source_housing, source_housing[1, ]), audit_manifest),
             "duplicate source key accepted")
panel_bad <- panel; panel_bad$SAMPLE[1] <- NA_integer_
expect_error(join_first_birth_housing(panel_bad, source_housing, audit_manifest),
             "missing SAMPLE accepted")

# Code exact values across years: ROOMS 0/1/6/7/28/30 and BEDROOMS 0/1/6/7/22/23.
code_panel <- joined[seq_len(6), , drop = FALSE]
code_panel$YEAR <- c(2007L, 2007L, 2008L, 2008L, 2009L, 2009L)
code_panel$raw_rooms <- c(0, 1, 6, 7, 28, 30)
code_panel$raw_bedrooms <- c(0, 1, 6, 7, 22, 23)
code_panel$raw_ownership <- c(0, 1, 2, 9, 1, 2)
coded <- code_first_birth_housing(code_panel, coding)
expect(is.na(coded$rooms9[1]) && coded$rooms9[2] == 1 && coded$rooms9[3] == 6 &&
       coded$rooms9[4] == 7 && is.na(coded$rooms9[5]) && coded$rooms9[6] == 9,
       "ROOMS coding/cap contract failed")
expect(is.na(coded$bedrooms5[1]) && coded$bedrooms5[2] == 0 && coded$bedrooms5[3] == 5 &&
       coded$bedrooms5[4] == 5 && coded$bedrooms5[5] == 5 && is.na(coded$bedrooms5[6]),
       "BEDROOMS coding/cap contract failed")
expect(is.na(coded$ownership_lw[1]) && coded$ownership_lw[2] == 1 && coded$ownership_lw[3] == 0 &&
       is.na(coded$ownership_lw[4]), "OWNERSHP coding contract failed")
expect(coded$rooms_unknown_code[5] && coded$rooms_invalid_code[5], "ROOMS=28 flag missing")
expect(coded$bedrooms_unknown_code[6] && coded$bedrooms_invalid_code[6], "BEDROOMS=23 flag missing")
expect(!("rooms_uncapped" %in% names(coded)), "uncapped sensitivity produced without config permission")
expect_error(fbh_term("statenameMaine:t_es10", "Maine", "1"),
             "event 1 matched event 10")

# Make one duplicate row an explicitly configured unknown rooms code while a
# same-cell copy remains valid, so estimation support remains complete.
source_housing$ROOMS[1] <- 28
checkpoints <- list()
result <- estimate_first_birth_housing(panel, source_housing, audit_manifest, coding,
  event_times = events, pre_times = c("-2", "-1"), post_times = c("0", "1", "2", "3"),
  checkpoint = function(x) checkpoints[[length(checkpoints) + 1L]] <<- x)
expect(result$status == "ESTIMATION_COMPLETE_DIAGNOSTIC", "estimator did not complete tiny diagnostic")
expect(length(checkpoints) == 6L, "per-outcome/gender checkpoints were not emitted")
expect(is.numeric(checkpoints[[1L]]$coefficients) && is.matrix(checkpoints[[1L]]$vcov$cluster) &&
       is.matrix(checkpoints[[1L]]$vcov$heteroskedastic),
       "checkpoint omitted recoverable coefficients/covariances")
expect(all(result$curves$estimate[result$curves$reference] == 0), "reference event is not zero")
expect(all(is.finite(result$curves$estimate[!result$curves$reference])), "non-reference coefficient missing")
expect(all(result$curves$source_household_clusters < result$curves$nobs),
       "primary inference did not reuse source household clusters")
expect(all(result$summary$contrast %in% c("post_minus_pre", "event_3_minus_event_neg1")),
       "contrast labels malformed")
expect(all(is.finite(result$summary$std.error)) && all(result$summary$std.error >= 0),
       "cluster contrast covariance failed")
vt3 <- result$summary$statename == "Vermont" & result$summary$gender == "Men" &
  result$summary$outcome == "ownership_lw" & result$summary$contrast == "event_3_minus_event_neg1"
expect(sum(vt3) == 1L, "known Vermont ownership contrast is not unique")
expect(all(abs(result$summary$estimate[vt3] - 1) < 1e-8),
       "known Vermont ownership +3-minus--1 contrast changed")
expect(max(abs(result$curves$estimate[result$curves$outcome == "ownership_lw"])) < 2,
       "ownership coefficients were incorrectly scaled")
expect(max(abs(result$curves$estimate[result$curves$outcome == "rooms9"])) < 10,
       "rooms coefficients were incorrectly normalized")
expect(any(result$code_audit$unknown_code > 0), "unknown code audit was not returned")

cat("PASS: first-birth housing join, coding, level event curves, source clustering, and full-covariance contrasts\n")
