#!/usr/bin/env Rscript
# End-to-end deterministic national estimator fixture. It fits the actual
# four fixest specifications and compares the pooled full fit and covariance
# to a separately constructed direct fixest call.
suppressPackageStartupMessages(library(fixest))
script_file <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_file[grepl("^--file=", script_file)][1L])
source(file.path(dirname(normalizePath(script_file)), "estimate_national_first_birth_housing.R"), local = TRUE)

expect <- function(ok, msg) if (!isTRUE(ok)) stop(paste("FAIL:", msg), call. = FALSE)
expect_error <- function(expr, msg) {
  hit <- FALSE; tryCatch(force(expr), error = function(e) hit <<- TRUE)
  expect(hit, msg)
}

set.seed(170921)
ev <- -5:10
g <- expand.grid(statefip = c(6L, 17L, 50L), age = 25:28, year = 2005:2007,
                 event_time = ev, rep = 1:2, KEEP.OUT.ATTRS = FALSE)
g <- g[order(g$statefip, g$year, g$age, g$event_time, g$rep), ]
g$gender <- "Women"
g$age_factor <- factor(g$age)
g$doiy_factor <- factor(g$year)
g$wgt <- 1 + (seq_len(nrow(g)) %% 5)
g$wgt[2] <- NA_real_ # missing author weight remains missing and is excluded
g$event_time <- as.integer(g$event_time)
g$SAMPLE <- 200501L
g$SERIAL <- as.integer((seq_len(nrow(g)) - 1L) %% 80L + 1L) # reused clusters
g$source_year.x <- g$year
g$source_year <- g$year + 1000L # proves .x precedence
g$source_sample.x <- g$SAMPLE
g$source_sample <- g$SAMPLE + 1L
g$source_serial.x <- g$SERIAL
g$source_serial <- g$SERIAL + 1000L
g$source_origin.x <- "ACS"
g$source_origin <- "CPS"
g$from_cps.x <- 0L
g$from_cps <- 1L
g$source_statefip <- g$statefip
g$source_age <- g$age
g$rooms_raw <- pmax(1, pmin(30, 6 + (g$event_time >= 0) + (g$statefip == 50)))
g$ownershp_raw <- ifelse(g$statefip == 50 & g$event_time >= 0, 1, 2)
g$bedrooms_raw <- ifelse(g$event_time >= 0, 5L, 4L)

# Unknown/missing housing codes and a valid CPS housing value are preserved but
# cannot enter their outcome fits. The CPS row also checks no CPS imputation.
g$rooms_raw[3] <- 28
g$rooms_raw[4] <- 0
g$ownershp_raw[5] <- 0
g$source_origin.x[nrow(g)] <- "CPS"
g$from_cps.x[nrow(g)] <- 1L
g$source_year.x[nrow(g)] <- NA_real_
g$source_sample.x[nrow(g)] <- NA_real_
g$source_serial.x[nrow(g)] <- NA_real_
g$rooms_raw[nrow(g)] <- 30

outdir <- tempfile("national_first_birth_housing_")
dir.create(outdir)
checkpoints <- list()
res <- estimate_national_first_birth_housing(
  g, output_dir = outdir, event_times = ev,
  checkpoint = function(x) checkpoints[[length(checkpoints) + 1L]] <<- x)

expect(res$status == "ESTIMATION_COMPLETE_DIAGNOSTIC", "national estimator did not complete")
expect(nrow(res$fit_status) == 12L, "four specifications were not fit for each of three outcomes")
expect(all(c("rooms9", "ownership_lw", "bedrooms5") %in% unique(res$curves$outcome)), "outcomes missing")
expect(all(res$curves$event_time %in% ev), "curve event support changed")
expect(all(res$curves$estimate[res$curves$reference] == 0), "reference event is not zero")
expect(all(is.finite(res$contrasts$std.error) & res$contrasts$std.error >= 0), "contrast SE invalid")
expect(all(res$fit_status$nobs_differs_from_fullspec == FALSE), "specifications changed the outcome sample")
expect(all(res$fit_status$nobs < nrow(g)), "missing weights/CPS rows were not excluded")
expect(any(res$counts$missing_outcome > 0), "housing missingness was not reported")
expect(any(res$counts$event_ess < res$counts$nobs), "event ESS did not use weights")
expect(any(res$data$rooms_unknown_code) && any(res$data$rooms_missing_code), "ROOMS unknown/missing flags lost")
expect(any(res$data$ownership_missing_code), "ownership unknown code was not retained")
expect(identical(res$data$.row_id, seq_len(nrow(g))), "row identity changed")
expect(all(res$data$source_year[seq_len(nrow(g) - 1L)] == g$source_year.x[seq_len(nrow(g) - 1L)]), ".x source year was not authoritative")
expect(all(res$data$source_sample[seq_len(nrow(g) - 1L)] == g$source_sample.x[seq_len(nrow(g) - 1L)]), ".x source sample was not authoritative")
expect(all(res$data$source_serial[seq_len(nrow(g) - 1L)] == g$source_serial.x[seq_len(nrow(g) - 1L)]), ".x source serial was not authoritative")
expect(length(unique(res$data$source_hh_cluster)) < nrow(g), "household cluster did not omit PERNUM/reuse households")

# Direct full fixest comparison on the exact same women/true-ACS/positive-weight
# rooms sample used by the estimator.
z <- res$data
z <- z[res$data$source_origin.x == "ACS" & res$data$from_cps.x == 0 &
       is.finite(res$data$.weight) & res$data$.weight > 0 & !is.na(res$data$rooms9) &
       res$data$.event %in% ev, ]
z$event_time <- z$.event
z$statefip <- factor(z$statefip)
z$age_factor <- factor(z$age_factor)
z$doiy_factor <- factor(z$doiy_factor)
direct <- feols(rooms9 ~ i(event_time, ref = -2) |
                  statefip + age_factor + doiy_factor,
                data = z, weights = ~.weight, vcov = ~source_hh_cluster)
pooled <- res$fits[["rooms9::full"]]
expect(isTRUE(all.equal(unname(coef(direct)), unname(pooled$coefficients), tolerance = 1e-10)),
       "full coefficients differ from direct fixest")
expect(isTRUE(all.equal(unname(vcov(direct)), unname(pooled$full_vcov), tolerance = 1e-9)),
       "full covariance differs from direct fixest")
direct_contrast <- unname(coef(direct)[["event_time::3"]] - coef(direct)[["event_time::-1"]])
got_contrast <- res$contrasts$estimate[res$contrasts$outcome == "rooms9" & res$contrasts$specification == "full"]
expect(abs(direct_contrast - got_contrast) < 1e-10, "+3 minus -1 contrast differs from direct fixest")

expect(file.exists(file.path(outdir, "national_housing_event_curves.png")), "PNG diagnostic was not written")
expect(file.exists(file.path(outdir, "checkpoint_rooms9_full.rds")), "per-fit checkpoint was not written")
expect(length(checkpoints) == 24L, "per-fit start/completion checkpoints missing")

expect(is.na(res$data$ownership_lw[5]), "raw OWNERSHP=0 was not retained as missing")

cat("PASS: national pooled ACS full/event-only/age-only/state-year fits, full covariance contrast, source precedence, missingness, weights, clusters, checkpoints, and PNG\n")
