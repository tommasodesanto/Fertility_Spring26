#!/usr/bin/env Rscript
# Small deterministic tests for the national ACS matcher adapter.  The actual
# vendor cleaner/matcher is exercised by the VT Slurm smoke; these tests cover
# the source-identity and housing-bridge contracts without loading 59m rows.

options(stringsAsFactors = FALSE)
matcher_file <- "code/empirical/acs/kleven_pseudo/run_national_acs_match_housing.R"
source_file <- "code/empirical/acs/kleven_pseudo/run_national_acs_source_stage.R"
txt <- paste(readLines(matcher_file, warn = FALSE), collapse = "\n")
stxt <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
stopifnot(grepl("subset(select=c(doiy, month, serial, pernum, wgt, id,", txt, fixed = TRUE))
stopifnot(grepl("source_month", txt, fixed = TRUE), grepl("source_hhcluster", txt, fixed = TRUE))
stopifnot(grepl("sploc", txt, fixed = TRUE), grepl("fertyr", txt, fixed = TRUE))
stopifnot(grepl("verified_overlap", txt, fixed = TRUE), grepl("lhs_only_missing", txt, fixed = TRUE),
          grepl("rhs_only_missing", txt, fixed = TRUE), grepl("packet_states <- valid_states", txt, fixed = TRUE))
stopifnot(grepl("functions.R", stxt, fixed = TRUE), grepl("vendor_expected", stxt, fixed = TRUE))

# The unchanged vendor CPS cleaner carries the author weight as
# ifelse(is.na(asecwt), wtfinl, asecwt) (clean_cps.R:591).  Exercise the
# fallback explicitly so a missing ASECWT is not mistaken for a changed weight.
cps_weight_fixture <- data.frame(asecwt = c(NA_real_, 12), wtfinl = c(7, 8))
cps_weight_fixture$wgt <- with(cps_weight_fixture,
                                ifelse(is.na(asecwt), wtfinl, asecwt))
stopifnot(identical(as.numeric(cps_weight_fixture$wgt), c(7, 12)))

adapter_env <- new.env(parent = globalenv())
adapter_env$req <- function(ok, msg, stage = "test") if (!isTRUE(ok)) stop(msg, call. = FALSE)
adapter_env$estimator_pool_columns <- c("wgt","age_factor","doiy_factor","statefip","gender","event_time",
                                        "rooms_raw","ownershp_raw","bedrooms_raw","source_origin","from_cps",
                                        "source_year","source_sample","source_serial","source_pernum","source_hh_cluster")
adapter_env$resolve <- function(nms, want, required = TRUE) {
  hit <- nms[tolower(nms) == tolower(want)]
  if (!length(hit) && !required) return(NA_character_)
  if (length(hit) != 1L) stop("test resolver failed", call. = FALSE)
  hit[[1L]]
}
exprs <- parse(file = matcher_file)
extract_function <- function(name) {
  for (e in exprs) if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
      identical(e[[2L]], as.name(name)) && is.call(e[[3L]]) &&
      identical(e[[3L]][[1L]], as.name("function"))) {
    eval(e, envir = adapter_env)
    return(get(name, envir = adapter_env, inherits = FALSE))
  }
  stop("adapter function absent: ", name, call. = FALSE)
}
coalesce_field <- extract_function("coalesce_field")
normalize <- extract_function("normalize_lineage")
add_lineage <- extract_function("add_lineage")
narrow_pool <- extract_function("narrow_estimator_panel")

x <- data.frame(
  src_key.x = c("11:2005:10:1", NA), src_key = c(NA, "CPS:2005:3:20:1"),
  source_origin.x = c("ACS", NA), source_origin = c(NA, "CPS"),
  source_year.x = c(2005L, NA), source_year = c(NA, 2005L),
  source_doiy.x = c(2005L, NA), source_doiy = c(NA, 2005L),
  source_month.x = c(NA, NA), source_month = c(NA, 3L),
  source_sample.x = c(11L, NA), source_sample = c(NA, NA),
  source_serial.x = c(10L, NA), source_serial = c(NA, 20L),
  source_pernum.x = c(1L, NA), source_pernum = c(NA, 1L),
  source_sex.x = c(2L, NA), source_sex = c(NA, 2L),
  source_age.x = c(30L, NA), source_age = c(NA, 30L),
  source_hhcluster.x = c("11:2005:10", NA), source_hhcluster = c(NA, "CPS:2005:3:20"),
  source_ownershp_raw.x = c(1L, NA), source_ownershp_raw = c(NA, NA),
  ownershp_raw.x = c(NA, NA), ownershp_raw = c(NA, NA),
  rooms_raw.x = c(NA, NA), rooms_raw = c(NA, NA),
  bedrooms_raw.x = c(NA, NA), bedrooms_raw = c(NA, NA),
  from_cps.x = c(0L, NA), from_cps = c(NA, 1L)
)
n <- normalize(x)
stopifnot(all(as.integer(n$from_cps) == c(0L, 1L)))
stopifnot(identical(n$source_origin, c("ACS", "CPS")),
          identical(n$source_hhcluster, c("11:2005:10", "CPS:2005:3:20")),
          !any(c("source_origin.x", "source_origin.y", "from_cps.x",
                 "source_hhcluster.x") %in% names(n)))

pool_fixture <- data.frame(
  src_key = c("11:2005:10:1", "CPS:2005:3:20:1"), source_origin = c("ACS", "CPS"),
  source_doiy = c(2005L, 2005L), source_year = c(2005L, 2005L), source_month = c(NA, 3L),
  source_sample = c(11L, NA), source_serial = c(10L, 20L), source_pernum = c(1L, 1L),
  source_sex = c(2L, 2L), source_age = c(30L, 30L), source_hhcluster = c("11:2005:10", NA),
  source_ownershp_raw = c(1L, NA), from_cps = c(0L, 1L), ownershp_raw = c(1L, NA),
  rooms_raw = c(5, NA), bedrooms_raw = c(2, NA), wgt = c(1, 2),
  age_factor = factor(c(25, 26)), doiy_factor = factor(c(2005, 2005)),
  statefip = c(50L, 50L), gender = c("Women", "Women"), event_time = c(-1L, 0L),
  author_control = c("keep_a", "keep_b"), stringsAsFactors = FALSE)
pool_one <- narrow_pool(pool_fixture)
pool_two <- dplyr::bind_rows(narrow_pool(pool_fixture[1, , drop = FALSE]),
                             narrow_pool(pool_fixture[2, , drop = FALSE]))
pool_wide_then_narrow <- narrow_pool(dplyr::bind_rows(pool_fixture, pool_fixture))
pool_cols <- adapter_env$estimator_pool_columns
stopifnot(identical(names(pool_one), pool_cols),
          isTRUE(all.equal(pool_two, pool_one, check.attributes = TRUE)),
          nrow(pool_wide_then_narrow) == 4L,
          !"author_control" %in% names(pool_one))

acs_raw <- data.frame(YEAR = 2005L, SAMPLE = 11L, SERIAL = 10L, PERNUM = 1L,
                      SEX = 2L, AGE = 30L, OWNERSHP = 1L)
acs_lineage <- add_lineage(acs_raw, "ACS")
stopifnot(acs_lineage$src_key == "11:2005:10:1",
          acs_lineage$source_hhcluster == "11:2005:10",
          acs_lineage$source_ownershp_raw == 1L,
          acs_lineage$from_cps == 0L)
cps_raw <- data.frame(YEAR = c(2005L, 2005L), MONTH = c(3L, 4L), SERIAL = c(20L, 20L),
                      PERNUM = c(1L, 1L), SEX = c(2L, 2L), AGE = c(30L, 30L))
cps_lineage <- add_lineage(cps_raw, "CPS")
stopifnot(!anyDuplicated(cps_lineage$src_key),
          cps_lineage$src_key[1] == "CPS:2005:3:20:1",
          cps_lineage$src_key[2] == "CPS:2005:4:20:1")

h <- data.frame(SAMPLE = 11L, YEAR = 2005L, SERIAL = 10L, PERNUM = 1L,
                SEX = 2L, AGE = 30L, OWNERSHP_RAW = 1L,
                ROOMS_RAW = 5L, BEDROOMS_RAW = 2L)
panel_key <- paste(n$source_sample[1], n$source_year[1], n$source_serial[1], n$source_pernum[1], sep = ":")
housing_key <- paste(h$SAMPLE, h$YEAR, h$SERIAL, h$PERNUM, sep = ":")
stopifnot(identical(match(panel_key, housing_key), 1L),
          n$source_origin[1] == "ACS", n$from_cps[1] == 0L,
          n$source_hhcluster[1] == "11:2005:10",
          h$SEX == 2L, h$AGE == 30L, h$OWNERSHP_RAW == 1L)
cat("NATIONAL_ACS_MATCH_HOUSING_TEST_PASS\n")
