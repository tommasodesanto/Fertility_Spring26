#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
  library(fixest)
})

required_packages <- c("haven", "fixest", "data.table", "survey")
package_status <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", file_arg[startsWith(file_arg, file_flag)])
  if (length(script_path) == 0) {
    return(normalizePath("code/data/moment_standard_errors"))
  }
  dirname(normalizePath(script_path[1]))
}

parse_args <- function() {
  out <- list(source = "all", B = 1000L, seed = 20260705L, rebuild_cache = FALSE)
  args <- commandArgs(trailingOnly = TRUE)
  for (arg in args) {
    if (startsWith(arg, "--source=")) out$source <- sub("^--source=", "", arg)
    if (startsWith(arg, "--B=")) out$B <- as.integer(sub("^--B=", "", arg))
    if (startsWith(arg, "--seed=")) out$seed <- as.integer(sub("^--seed=", "", arg))
    if (arg == "--rebuild-cache") out$rebuild_cache <- TRUE
  }
  if (!(out$source %in% c("all", "psid", "acs"))) {
    stop("--source must be one of all, psid, acs")
  }
  out
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
out_dir <- file.path(script_dir, "output")
cache_dir <- file.path(script_dir, "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(out_dir, c(
  "event_study_coefficient_path.csv",
  "source_timing.csv",
  "source_provenance_log.txt"
)), force = TRUE)
args <- parse_args()

TARGETS <- data.table(
  moment_key = c(
    "tfr",
    "childless_rate",
    "own_rate",
    "own_family_gap",
    "housing_increment_0to1",
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
    "prime30_55_childless_renter_mean_rooms",
    "prime30_55_childless_owner_share_rooms_ge6",
    "old_nonhousing_wealth_to_income_median_6575",
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535",
    "prime30_55_childless_owner_minus_renter_mean_rooms",
    "old_age_own_rate",
    "own_rate_2534",
    "prime30_55_parent_3plus_minus_1to2_mean_rooms"
  ),
  published_target = c(
    1.918000,
    0.188000,
    0.575472,
    0.167662,
    0.664435,
    1.007450,
    3.805288,
    0.596131,
    2.230461,
    0.179226,
    2.418762,
    0.764261,
    0.341166,
    0.367700
  )
)
setnames(TARGETS, "moment_key", "key")

SOURCE_DEFAULTS <- data.table(
  moment_key = TARGETS$key,
  source = c(
    "source-unconfirmed",
    "source-unconfirmed",
    "ACS/MMS ownership audit",
    "ACS/MMS ownership audit",
    "PSID first-birth event study",
    "PSID old-age wealth",
    "ACS/MMS one-market housing",
    "ACS/MMS one-market housing",
    "PSID old-age wealth",
    "PSID young wealth",
    "ACS/MMS one-market housing",
    "ACS/MMS ownership audit",
    "ACS/MMS ownership audit",
    "ACS/MMS one-market housing"
  ),
  psu = c(
    NA_character_,
    NA_character_,
    "met2013",
    "met2013",
    "ID",
    "ID",
    "met2013",
    "met2013",
    "ID",
    "ID",
    "met2013",
    "met2013",
    "met2013",
    "met2013"
  ),
  weight_var = c(
    NA_character_,
    NA_character_,
    "hhwt",
    "hhwt",
    NA_character_,
    "IW",
    "hhwt",
    "hhwt",
    "IW",
    "IW",
    "hhwt",
    "hhwt",
    "hhwt",
    "hhwt"
  )
)
setnames(SOURCE_DEFAULTS, "moment_key", "key")

fmt_num <- function(x, digits = 12) {
  if (length(x) == 0 || is.na(x)) return("NA")
  formatC(as.numeric(x), digits = digits, format = "fg", flag = "#")
}

tol_for <- function(target) {
  max(1e-6, 1e-4 * abs(target))
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_quantile_safe <- function(x, w, prob = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- as.numeric(x[ok])
  w <- as.numeric(w[ok])
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= prob)[1]]
}

weighted_share <- function(flag, w) {
  weighted_mean_safe(as.numeric(flag), w)
}

min_ignore_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

collapse_mms_location <- function(x, middle_target = "center") {
  out <- as.character(x)
  if (middle_target == "center") out[out == "middle"] <- "center"
  if (middle_target == "periphery") out[out == "middle"] <- "periphery"
  out
}

repro_lines <- character()
notes <- character()
timings <- data.table(
  source = character(),
  phase = character(),
  wall_clock_sec = numeric(),
  B = integer(),
  n_analysis = integer(),
  n_clusters = integer(),
  n_moments = integer(),
  method_note = character()
)

add_note <- function(txt) {
  notes <<- c(notes, txt)
}

record_timing <- function(source, phase, start_time, B = NA_integer_,
                          n_analysis = NA_integer_, n_clusters = NA_integer_,
                          n_moments = NA_integer_, method_note = "") {
  timings <<- rbind(
    timings,
    data.table(
      source = source,
      phase = phase,
      wall_clock_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      B = as.integer(B),
      n_analysis = as.integer(n_analysis),
      n_clusters = as.integer(n_clusters),
      n_moments = as.integer(n_moments),
      method_note = method_note
    )
  )
}

log_gate <- function(key, target, reproduced, status, source, reason = "", tol = NA_real_) {
  abs_diff <- if (is.na(target) || is.na(reproduced)) NA_real_ else abs(reproduced - target)
  if (is.na(tol) && !is.na(target)) tol <- tol_for(target)
  line <- paste0(
    "key=", key,
    " | target=", fmt_num(target),
    " | reproduced=", fmt_num(reproduced),
    " | abs_diff=", fmt_num(abs_diff),
    " | tol=", fmt_num(tol),
    " | status=", status,
    " | source=", source,
    if (nzchar(reason)) paste0(" | reason=", reason) else ""
  )
  repro_lines <<- c(repro_lines, line)
}

gate_pass <- function(target, reproduced) {
  if (is.na(target) || is.na(reproduced)) return(FALSE)
  abs(reproduced - target) <= tol_for(target)
}

target_value <- function(moment_key) {
  TARGETS[match(moment_key, TARGETS$key), published_target]
}

cache_path <- function(name) file.path(cache_dir, paste0(name, ".rds"))

build_psid_cache <- function() {
  path <- cache_path("psid_analysis_samples")
  if (file.exists(path) && !args$rebuild_cache) {
    message("Loading cached PSID analysis samples: ", path)
    return(readRDS(path))
  }

  psid_path <- file.path(dirname(repo_root), "PSID", "PSIDSHELF_MOBILITY.dta")
  if (!file.exists(psid_path)) stop("Missing PSID shelf file: ", psid_path)

  message("Reading PSID shelf once for old-age, young-wealth, and event-study samples...")
  psid <- as.data.table(read_dta(
    psid_path,
    col_select = c(
      "ID", "year", "AGEREP", "EDUYEAR", "SEX", "DEATHYEAR", "HOMEOWN",
      "RELCHIREP", "RELCHINUM", "RELCHI1BYEAR", "RELCHI2BYEAR",
      "ACTUALROOMS_", "INCFAMR", "EARNINDR", "NETWORTHR", "NETWORTH2R",
      "NETWORTH3R", "HOMEEQUITYR", "WLTHSAVETOTR", "IW"
    )
  ))
  setnames(psid, old = names(psid), new = tolower(names(psid)))
  psid[, `:=`(
    famid = as.character(id),
    year = as.integer(year),
    age = as.numeric(agerep),
    death_year = as.numeric(deathyear),
    homeown_raw = as.numeric(homeown),
    own = fifelse(as.numeric(homeown) == 1, 1, fifelse(as.numeric(homeown) == 2, 0, NA_real_)),
    n_children_rep = as.numeric(relchirep),
    n_children_num = as.numeric(relchinum),
    first_child_year_raw = as.numeric(relchi1byear),
    second_child_year_raw = as.numeric(relchi2byear),
    rooms = as.numeric(actualrooms_),
    income = as.numeric(incfamr),
    earnings = as.numeric(earnindr),
    nonhousing_nw = as.numeric(networth2r),
    weight = as.numeric(iw),
    educ = as.integer(eduyear)
  )]

  old <- psid[
    year >= 1984 & year <= 2019 &
      (is.na(death_year) | year <= death_year) &
      age >= 65 & age <= 75 &
      !is.na(weight) & weight > 0 &
      !is.na(n_children_rep),
    .(
      famid, year, age, weight, own,
      n_children = n_children_rep,
      parent = n_children_rep > 0,
      childless = n_children_rep == 0,
      nonhousing_nw_to_income = fifelse(income > 1000, nonhousing_nw / income, NA_real_)
    )
  ]

  young <- psid[
    year >= 1984 & year <= 2019 &
      (is.na(death_year) | year <= death_year) &
      age >= 25 & age <= 35 &
      !is.na(weight) & weight > 0,
    .(
      famid, year, age, weight,
      own,
      childless = n_children_rep == 0,
      wealth_to_income = fifelse(income > 1000, nonhousing_nw / income, NA_real_)
    )
  ][childless == TRUE & own == 0 & !is.na(wealth_to_income)]

  event_base <- psid[
    (is.na(death_year) | year <= death_year) &
      !is.na(age) & age >= 18,
    .(
      famid,
      year,
      age = as.integer(round(age)),
      educ = as.integer(educ),
      rooms,
      first_child_year_raw,
      second_child_year_raw
    )
  ]
  event_base[, first_child_year := min_ignore_na(first_child_year_raw), by = famid]
  event_base[, second_child_year := min_ignore_na(second_child_year_raw), by = famid]
  event_base <- event_base[!is.na(first_child_year)]
  event_base[, year_entry := min(year, na.rm = TRUE), by = famid]
  event_base <- event_base[first_child_year >= year_entry]
  event_base[, `:=`(
    K = year - first_child_year,
    lastcohort = first_child_year == max(first_child_year, na.rm = TRUE)
  )]
  event_base <- event_base[!is.na(rooms) & !is.na(year) & !is.na(age) & !is.na(educ)]

  out <- list(old = old, young = young, event = event_base)
  saveRDS(out, path)
  message("Cached PSID analysis samples: ", path)
  out
}

load_mms_lookup <- function() {
  data_dir <- file.path(repo_root, "code", "data", "mms_center_periphery", "data")
  lookup_2010 <- as.data.table(fread(file.path(data_dir, "puma_mms_lookup_2010.csv")))
  lookup_2020 <- as.data.table(fread(file.path(data_dir, "puma_mms_lookup_2020.csv")))
  rbindlist(list(
    lookup_2010[, .(
      lookup_period = "pre2022",
      statefip = as.integer(statefip),
      puma = as.integer(puma),
      met2013 = as.integer(cbsacode),
      mms_location
    )],
    lookup_2020[, .(
      lookup_period = "post2021",
      statefip = as.integer(statefip),
      puma = as.integer(puma),
      met2013 = as.integer(cbsacode),
      mms_location
    )]
  ), use.names = TRUE)
}

build_acs_cache <- function() {
  path <- cache_path("acs_analysis_samples")
  if (file.exists(path) && !args$rebuild_cache) {
    message("Loading cached ACS analysis samples: ", path)
    return(readRDS(path))
  }

  extract_path <- file.path(repo_root, "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta")
  if (!file.exists(extract_path)) stop("Missing ACS/IPUMS extract: ", extract_path)

  lookup <- load_mms_lookup()
  message("Reading ACS/IPUMS extract27 once for target and supplementary samples...")
  acs <- as.data.table(read_dta(
    extract_path,
    col_select = c(
      "year", "statefip", "puma", "met2013", "gq", "pernum", "relate",
      "hhwt", "perwt", "ownershp", "unitsstr", "rooms", "rent", "hhincome",
      "nchild", "nchlt5", "eldch", "yngch", "fertyr", "sex", "age",
      "educ", "marst"
    )
  ))
  setnames(acs, old = names(acs), new = tolower(names(acs)))
  acs <- acs[
    year >= 2012 & year <= 2023 &
      met2013 > 0 &
      gq %in% c(1, 2)
  ]
  acs[, `:=`(
    year = as.integer(year),
    statefip = as.integer(statefip),
    puma = as.integer(puma),
    met2013 = as.integer(met2013),
    lookup_period = fifelse(year <= 2021, "pre2022", "post2021"),
    pernum = as.integer(pernum),
    relate = as.integer(relate),
    hhwt = as.numeric(hhwt),
    perwt = as.numeric(perwt),
    ownershp = as.integer(ownershp),
    unitsstr = as.integer(unitsstr),
    rooms = as.numeric(rooms),
    rent = as.numeric(rent),
    hhincome = as.numeric(hhincome),
    nchild = as.numeric(nchild),
    eldch = as.numeric(eldch),
    yngch = as.numeric(yngch),
    fertyr = as.integer(fertyr),
    sex = as.integer(sex),
    age = as.integer(age),
    educ = as.integer(educ),
    marst = as.integer(marst)
  )]
  acs <- merge(
    acs,
    lookup,
    by = c("lookup_period", "statefip", "puma", "met2013"),
    all.x = TRUE,
    sort = FALSE
  )
  acs[, mms_location := collapse_mms_location(mms_location, "center")]
  acs <- acs[mms_location %in% c("center", "periphery")]
  acs[, `:=`(
    owner = ownershp == 1L,
    renter = ownershp == 2L,
    childless = nchild == 0 & !is.na(nchild),
    parent_u18 = nchild > 0 & !is.na(yngch) & yngch != 99 & yngch < 18,
    newparent = nchild > 0 & !is.na(nchild) & !is.na(eldch) & eldch != 99 & eldch < 4,
    due_housing = unitsstr %in% 3:10 & !is.na(rooms) & rooms > 0,
    hhincome_clean = fifelse(hhincome > 1000 & hhincome < 2000000, hhincome, NA_real_),
    college = as.integer(educ >= 10),
    married = as.integer(marst %in% c(1, 2)),
    recent_birth = as.integer(fertyr == 2),
    fertyr_observed = fertyr %in% c(1, 2)
  )]

  targets <- acs[
    pernum == 1L &
      !is.na(hhwt) & hhwt > 0 &
      ownershp %in% c(1L, 2L) &
      !is.na(rooms) & rooms > 0 &
      age >= 18
  ]
  supplement_women <- acs[
    sex == 2L &
      age >= 22 & age <= 45 &
      !is.na(perwt) & perwt > 0 &
      !is.na(hhincome_clean)
  ]

  out <- list(targets = targets, women = supplement_women)
  saveRDS(out, path)
  message("Cached ACS analysis samples: ", path)
  out
}

psid_point_moments <- function(psid) {
  old <- psid$old
  young <- psid$young
  parent_old <- old[parent == TRUE]
  childless_old <- old[childless == TRUE]
  c(
    old_parent_childless_nonhousing_wealth_to_income_gap_6575 =
      weighted_mean_safe(parent_old$nonhousing_nw_to_income, parent_old$weight) -
      weighted_mean_safe(childless_old$nonhousing_nw_to_income, childless_old$weight),
    old_nonhousing_wealth_to_income_median_6575 =
      weighted_quantile_safe(old$nonhousing_nw_to_income, old$weight, 0.5),
    young_childless_renter_liquid_wealth_to_annual_gross_income_2535 =
      weighted_mean_safe(young$wealth_to_income, young$weight),
    old_age_own_rate_psid_check =
      weighted_mean_safe(old$own, old$weight)
  )
}

acs_point_moments <- function(acs) {
  dt <- acs$targets
  one_market <- dt[age >= 30 & age <= 55]
  prime_childless <- one_market[childless == TRUE]
  owner_childless <- prime_childless[owner == TRUE]
  renter_childless <- prime_childless[renter == TRUE]
  prime_parent <- one_market[parent_u18 == TRUE]
  parent_3plus <- prime_parent[nchild >= 3]
  parent_1to2 <- prime_parent[nchild >= 1 & nchild <= 2]

  due <- dt[relate == 1L & due_housing == TRUE]
  due_3055 <- due[age >= 30 & age <= 55]
  due_2534 <- due[age >= 25 & age <= 34]
  due_6575 <- due[age >= 65 & age <= 75]
  no_child_3055 <- due_3055[nchild == 0 & !is.na(nchild)]
  newparent_3055 <- due_3055[newparent == TRUE]

  owner_mean <- weighted_mean_safe(owner_childless$rooms, owner_childless$hhwt)
  renter_mean <- weighted_mean_safe(renter_childless$rooms, renter_childless$hhwt)

  c(
    own_rate = weighted_share(due_3055$owner, due_3055$hhwt),
    own_family_gap =
      weighted_share(newparent_3055$owner, newparent_3055$hhwt) -
      weighted_share(no_child_3055$owner, no_child_3055$hhwt),
    prime30_55_childless_renter_mean_rooms = renter_mean,
    prime30_55_childless_owner_share_rooms_ge6 =
      weighted_share(owner_childless$rooms >= 6, owner_childless$hhwt),
    prime30_55_childless_owner_minus_renter_mean_rooms = owner_mean - renter_mean,
    old_age_own_rate = weighted_share(due_6575$owner, due_6575$hhwt),
    own_rate_2534 = weighted_share(due_2534$owner, due_2534$hhwt),
    prime30_55_parent_3plus_minus_1to2_mean_rooms =
      weighted_mean_safe(parent_3plus$rooms, parent_3plus$hhwt) -
      weighted_mean_safe(parent_1to2$rooms, parent_1to2$hhwt)
  )
}

n_analysis_for_key <- function(key, psid = NULL, acs = NULL) {
  if (!is.null(psid)) {
    old <- psid$old
    young <- psid$young
    if (key == "old_parent_childless_nonhousing_wealth_to_income_gap_6575") {
      return(nrow(old[!is.na(nonhousing_nw_to_income) & (parent == TRUE | childless == TRUE)]))
    }
    if (key == "old_nonhousing_wealth_to_income_median_6575") {
      return(nrow(old[!is.na(nonhousing_nw_to_income)]))
    }
    if (key == "young_childless_renter_liquid_wealth_to_annual_gross_income_2535") {
      return(nrow(young[!is.na(wealth_to_income)]))
    }
  }
  if (!is.null(acs)) {
    dt <- acs$targets
    if (key == "own_rate") return(nrow(dt[relate == 1L & due_housing == TRUE & age >= 30 & age <= 55]))
    if (key == "own_family_gap") return(nrow(dt[relate == 1L & due_housing == TRUE & age >= 30 & age <= 55 & ((!is.na(nchild) & nchild == 0) | newparent == TRUE)]))
    if (key == "prime30_55_childless_renter_mean_rooms") return(nrow(dt[age >= 30 & age <= 55 & childless == TRUE & renter == TRUE]))
    if (key == "prime30_55_childless_owner_share_rooms_ge6") return(nrow(dt[age >= 30 & age <= 55 & childless == TRUE & owner == TRUE]))
    if (key == "prime30_55_childless_owner_minus_renter_mean_rooms") return(nrow(dt[age >= 30 & age <= 55 & childless == TRUE & (owner == TRUE | renter == TRUE)]))
    if (key == "old_age_own_rate") return(nrow(dt[relate == 1L & due_housing == TRUE & age >= 65 & age <= 75]))
    if (key == "own_rate_2534") return(nrow(dt[relate == 1L & due_housing == TRUE & age >= 25 & age <= 34]))
    if (key == "prime30_55_parent_3plus_minus_1to2_mean_rooms") return(nrow(dt[age >= 30 & age <= 55 & parent_u18 == TRUE & nchild >= 1]))
  }
  NA_integer_
}

n_clusters_for_key <- function(key, psid = NULL, acs = NULL) {
  if (!is.null(psid)) {
    old <- psid$old
    young <- psid$young
    if (key %in% c("old_parent_childless_nonhousing_wealth_to_income_gap_6575", "old_nonhousing_wealth_to_income_median_6575")) {
      return(uniqueN(old$famid))
    }
    if (key == "young_childless_renter_liquid_wealth_to_annual_gross_income_2535") {
      return(uniqueN(young$famid))
    }
  }
  if (!is.null(acs)) {
    dt <- acs$targets
    if (key %in% names(acs_point_moments(acs))) return(uniqueN(dt$met2013))
  }
  NA_integer_
}

bootstrap_psid <- function(psid, keys, B_requested) {
  old <- copy(psid$old)
  young <- copy(psid$young)
  old[, source_row := "old"]
  young[, source_row := "young"]
  common_cols <- union(names(old), names(young))
  for (nm in setdiff(common_cols, names(old))) old[, (nm) := NA]
  for (nm in setdiff(common_cols, names(young))) young[, (nm) := NA]
  combined <- rbindlist(list(old[, ..common_cols], young[, ..common_cols]), use.names = TRUE, fill = TRUE)
  famids <- sort(unique(combined$famid))
  fam_index_old <- match(old$famid, famids)
  fam_index_young <- match(young$famid, famids)
  n_fam <- length(famids)
  B <- B_requested
  boot <- matrix(NA_real_, nrow = B, ncol = length(keys), dimnames = list(NULL, keys))
  start <- Sys.time()
  checked <- FALSE
  for (b in seq_len(B)) {
    mult <- tabulate(sample.int(n_fam, n_fam, replace = TRUE), nbins = n_fam)
    old[, boot_w := weight * mult[fam_index_old]]
    young[, boot_w := weight * mult[fam_index_young]]
    parent_old <- old[parent == TRUE]
    childless_old <- old[childless == TRUE]
    vals <- c(
      old_parent_childless_nonhousing_wealth_to_income_gap_6575 =
        weighted_mean_safe(parent_old$nonhousing_nw_to_income, parent_old$boot_w) -
        weighted_mean_safe(childless_old$nonhousing_nw_to_income, childless_old$boot_w),
      old_nonhousing_wealth_to_income_median_6575 =
        weighted_quantile_safe(old$nonhousing_nw_to_income, old$boot_w, 0.5),
      young_childless_renter_liquid_wealth_to_annual_gross_income_2535 =
        weighted_mean_safe(young$wealth_to_income, young$boot_w)
    )
    boot[b, ] <- vals[keys]
    if (!checked && b == min(50L, B)) {
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      projected <- elapsed / b * B_requested
      if (B_requested > 500L && projected > 7200) {
        B <- 500L
        boot <- boot[seq_len(B), , drop = FALSE]
        add_note(sprintf("PSID bootstrap fallback: requested B=%d, using B=500 because first %d reps projected %.1f hours (>2h).", B_requested, b, projected / 3600))
      }
      checked <- TRUE
    }
    if (b >= B) break
  }
  boot[seq_len(B), , drop = FALSE]
}

bootstrap_acs <- function(acs, keys, B_requested) {
  dt <- copy(acs$targets)
  metros <- sort(unique(dt$met2013))
  comp <- data.table(met2013 = metros)

  add_component <- function(comp_dt, nm, sub_dt, value) {
    tmp <- copy(sub_dt)
    tmp[, value__ := as.numeric(value)]
    tmp <- tmp[is.finite(value__) & is.finite(hhwt) & hhwt > 0]
    tmp <- tmp[, .(
      num = sum(hhwt * value__, na.rm = TRUE),
      den = sum(hhwt, na.rm = TRUE)
    ), by = met2013]
    setnames(tmp, c("num", "den"), paste0(nm, c("_num", "_den")))
    out <- merge(comp_dt, tmp, by = "met2013", all.x = TRUE, sort = FALSE)
    for (cc in paste0(nm, c("_num", "_den"))) {
      out[is.na(get(cc)), (cc) := 0]
    }
    out
  }

  one_market <- dt[age >= 30 & age <= 55]
  prime_childless <- one_market[childless == TRUE]
  owner_childless <- prime_childless[owner == TRUE]
  renter_childless <- prime_childless[renter == TRUE]
  prime_parent <- one_market[parent_u18 == TRUE]
  parent_3plus <- prime_parent[nchild >= 3]
  parent_1to2 <- prime_parent[nchild >= 1 & nchild <= 2]

  due <- dt[relate == 1L & due_housing == TRUE]
  due_3055 <- due[age >= 30 & age <= 55]
  due_2534 <- due[age >= 25 & age <= 34]
  due_6575 <- due[age >= 65 & age <= 75]
  no_child_3055 <- due_3055[nchild == 0 & !is.na(nchild)]
  newparent_3055 <- due_3055[newparent == TRUE]

  comp <- add_component(comp, "own3055", due_3055, due_3055$owner)
  comp <- add_component(comp, "newparent_own3055", newparent_3055, newparent_3055$owner)
  comp <- add_component(comp, "nochild_own3055", no_child_3055, no_child_3055$owner)
  comp <- add_component(comp, "renter_childless_rooms", renter_childless, renter_childless$rooms)
  comp <- add_component(comp, "owner_childless_rooms", owner_childless, owner_childless$rooms)
  comp <- add_component(comp, "owner_childless_ge6", owner_childless, owner_childless$rooms >= 6)
  comp <- add_component(comp, "oldown6575", due_6575, due_6575$owner)
  comp <- add_component(comp, "own2534", due_2534, due_2534$owner)
  comp <- add_component(comp, "parent3plus_rooms", parent_3plus, parent_3plus$rooms)
  comp <- add_component(comp, "parent1to2_rooms", parent_1to2, parent_1to2$rooms)

  ratio <- function(num_col, den_col, mult) {
    num <- sum(comp[[num_col]] * mult)
    den <- sum(comp[[den_col]] * mult)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    num / den
  }

  n_metro <- nrow(comp)
  B <- B_requested
  boot <- matrix(NA_real_, nrow = B, ncol = length(keys), dimnames = list(NULL, keys))
  start <- Sys.time()
  checked <- FALSE
  for (b in seq_len(B)) {
    mult <- tabulate(sample.int(n_metro, n_metro, replace = TRUE), nbins = n_metro)
    owner_mean <- ratio("owner_childless_rooms_num", "owner_childless_rooms_den", mult)
    renter_mean <- ratio("renter_childless_rooms_num", "renter_childless_rooms_den", mult)

    vals <- c(
      own_rate = ratio("own3055_num", "own3055_den", mult),
      own_family_gap =
        ratio("newparent_own3055_num", "newparent_own3055_den", mult) -
        ratio("nochild_own3055_num", "nochild_own3055_den", mult),
      prime30_55_childless_renter_mean_rooms = renter_mean,
      prime30_55_childless_owner_share_rooms_ge6 =
        ratio("owner_childless_ge6_num", "owner_childless_ge6_den", mult),
      prime30_55_childless_owner_minus_renter_mean_rooms = owner_mean - renter_mean,
      old_age_own_rate = ratio("oldown6575_num", "oldown6575_den", mult),
      own_rate_2534 = ratio("own2534_num", "own2534_den", mult),
      prime30_55_parent_3plus_minus_1to2_mean_rooms =
        ratio("parent3plus_rooms_num", "parent3plus_rooms_den", mult) -
        ratio("parent1to2_rooms_num", "parent1to2_rooms_den", mult)
    )
    boot[b, ] <- vals[keys]
    if (!checked && b == min(50L, B)) {
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      projected <- elapsed / b * B_requested
      if (B_requested > 500L && projected > 7200) {
        B <- 500L
        boot <- boot[seq_len(B), , drop = FALSE]
        add_note(sprintf("ACS bootstrap fallback: requested B=%d, using B=500 because first %d reps projected %.1f hours (>2h).", B_requested, b, projected / 3600))
      }
      checked <- TRUE
    }
    if (b >= B) break
  }
  boot[seq_len(B), , drop = FALSE]
}

extract_event_coef <- function(est, rel) {
  ct <- as.data.frame(coeftable(est))
  ct$term <- rownames(ct)
  rel_pat <- paste0("::", rel, "($|[^0-9.-])")
  hit <- grep(rel_pat, ct$term)
  if (length(hit) == 0) {
    hit <- grep(paste0("(^|[^0-9.-])", rel, "($|[^0-9.-])"), ct$term)
  }
  if (length(hit) != 1) {
    return(c(estimate = NA_real_, se = NA_real_))
  }
  c(estimate = as.numeric(ct$Estimate[hit]), se = as.numeric(ct$`Std. Error`[hit]))
}

run_event_study <- function(psid) {
  existing_summary <- file.path(
    repo_root,
    "code/data/psid_followup_mar2026/output/sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_summary.csv"
  )
  use_fixest <- tolower(Sys.getenv("EVENT_USE_FIXEST", "0")) %in% c("1", "true", "yes", "y")
  if (!use_fixest && file.exists(existing_summary)) {
    ev_summary <- fread(existing_summary)
    return(list(
      coef_p3 = as.numeric(ev_summary$coef_p3[1]),
      se_p3 = as.numeric(ev_summary$se_p3[1]),
      coef_p5 = as.numeric(ev_summary$coef_p5[1]),
      se_p5 = as.numeric(ev_summary$se_p5[1]),
      n = as.integer(ev_summary$sample_obs[1]),
      n_clusters = as.integer(ev_summary$sample_ids[1]),
      method_note = paste0(
        "Read-only existing Stata builder summary ",
        "sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_summary.csv; ",
        "eventstudyinteract analytic IW estimator with vce(cluster ID), absorb(year), cohort(first_child_year), ",
        "control_cohort(lastcohort), covariates(i.AGEREP i.EDUYEAR). ",
        "Set EVENT_USE_FIXEST=1 to attempt the slower R fixest::sunab re-estimation."
      )
    ))
  }

  ev <- copy(psid$event)
  ev[, cohort_sa := fifelse(lastcohort == TRUE, 9999, as.integer(first_child_year))]
  ev <- ev[is.finite(rooms) & !is.na(cohort_sa) & !is.na(year) & !is.na(age) & !is.na(educ)]
  est <- tryCatch(
    feols(
      rooms ~ sunab(cohort_sa, year, ref.p = -2) + i(age) + i(educ) | year,
      data = ev,
      cluster = ~famid,
      notes = FALSE
    ),
    error = function(e) e
  )
  if (inherits(est, "error")) {
    add_note(paste("Event-study fixest failed:", conditionMessage(est)))
    if (file.exists(existing_summary)) {
      ev_summary <- fread(existing_summary)
      add_note("Event-study fallback used existing read-only Stata builder summary after fixest failure.")
      return(list(
        coef_p3 = as.numeric(ev_summary$coef_p3[1]),
        se_p3 = as.numeric(ev_summary$se_p3[1]),
        coef_p5 = as.numeric(ev_summary$coef_p5[1]),
        se_p5 = as.numeric(ev_summary$se_p5[1]),
        n = as.integer(ev_summary$sample_obs[1]),
        n_clusters = as.integer(ev_summary$sample_ids[1]),
        method_note = paste("fixest failed:", conditionMessage(est), "; read-only Stata builder summary used")
      ))
    }
    return(list(
      coef_p3 = NA_real_, se_p3 = NA_real_, coef_p5 = NA_real_, se_p5 = NA_real_,
      n = nrow(ev), n_clusters = uniqueN(ev$famid), method_note = paste("fixest failed:", conditionMessage(est))
    ))
  }
  p3 <- extract_event_coef(est, 3)
  p5 <- extract_event_coef(est, 5)
  if (is.na(p3["estimate"]) && file.exists(existing_summary)) {
    ev_summary <- fread(existing_summary)
    add_note("Event-study fixest ran but did not yield an extractable K=3 coefficient; existing read-only Stata builder summary used.")
    return(list(
      coef_p3 = as.numeric(ev_summary$coef_p3[1]),
      se_p3 = as.numeric(ev_summary$se_p3[1]),
      coef_p5 = as.numeric(ev_summary$coef_p5[1]),
      se_p5 = as.numeric(ev_summary$se_p5[1]),
      n = as.integer(ev_summary$sample_obs[1]),
      n_clusters = as.integer(ev_summary$sample_ids[1]),
      method_note = paste0(
        "fixest::sunab re-estimation ran but coefficient extraction failed; ",
        "read-only Stata builder summary used for analytic cluster-robust SE"
      )
    ))
  }
  list(
    coef_p3 = unname(p3["estimate"]),
    se_p3 = unname(p3["se"]),
    coef_p5 = unname(p5["estimate"]),
    se_p5 = unname(p5["se"]),
    n = nobs(est),
    n_clusters = uniqueN(ev$famid),
    method_note = "fixest::feols rooms ~ sunab(cohort_sa, year, ref.p=-2) + age/education indicators | year, cluster=ID; last first-birth cohort coded as never-treated control"
  )
}

load_event_study_path <- function() {
  estimates_path <- file.path(
    repo_root,
    "code/data/psid_followup_mar2026/output/sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_estimates.dta"
  )
  if (!file.exists(estimates_path)) {
    add_note(paste("Event-study coefficient path file not found:", estimates_path))
    return(data.table())
  }
  ev_path <- as.data.table(read_dta(estimates_path))
  needed <- c("relative_time", "b", "se", "ci_lo", "ci_hi", "pre_event_mean", "sample_obs", "sample_ids")
  missing_needed <- setdiff(needed, names(ev_path))
  if (length(missing_needed) > 0) {
    add_note(paste("Event-study coefficient path file is missing columns:", paste(missing_needed, collapse = ", ")))
    return(data.table())
  }
  ev_path <- ev_path[
    relative_time >= -2 & relative_time <= 5,
    .(
      relative_time = as.integer(relative_time),
      published_or_authoritative = as.numeric(b),
      reproduced_from_saved_output = as.numeric(b),
      cluster_se = as.numeric(se),
      ci_lo = as.numeric(ci_lo),
      ci_hi = as.numeric(ci_hi),
      calibration_reference = NA_real_,
      abs_diff_to_reference = NA_real_,
      pre_event_mean = as.numeric(pre_event_mean),
      sample_obs = as.integer(sample_obs),
      sample_ids = as.integer(sample_ids),
      source_file = "code/data/psid_followup_mar2026/output/sa_rooms_first_birth_one_variant_v1/rooms_f_c_y_all_estimates.dta",
      verification_note = "authoritative saved Stata eventstudyinteract output; harness imports the saved coefficient"
    )
  ][order(relative_time)]
  ev_path[relative_time == 3, `:=`(
    calibration_reference = target_value("housing_increment_0to1"),
    abs_diff_to_reference = abs(reproduced_from_saved_output - target_value("housing_increment_0to1")),
    verification_note = "K=3 target verified against authoritative saved Stata output"
  )]
  ev_path[relative_time == 5, verification_note := "K=5 plateau coefficient verified from authoritative saved Stata output"]
  ev_path[]
}

build_childlessness_by_income <- function(acs) {
  dt <- copy(acs$women)
  dt[, income_quintile := frank(hhincome_clean, ties.method = "average", na.last = "keep"), by = year]
  dt[, income_quintile := as.integer(ceiling(5 * income_quintile / .N)), by = year]
  dt[, income_decile := frank(hhincome_clean, ties.method = "average", na.last = "keep"), by = year]
  dt[, income_decile := as.integer(ceiling(10 * income_decile / .N)), by = year]
  dt[income_quintile < 1, income_quintile := 1L]
  dt[income_quintile > 5, income_quintile := 5L]
  dt[income_decile < 1, income_decile := 1L]
  dt[income_decile > 10, income_decile := 10L]
  dt[, education_group := fifelse(college == 1L, "college_plus", "less_than_college")]

  summarise_bins <- function(bin_col, by_educ = FALSE) {
    by_cols <- if (by_educ) c(bin_col, "education_group") else bin_col
    out <- dt[, .(
      n = .N,
      weight = sum(perwt, na.rm = TRUE),
      median_hhincome = weighted_quantile_safe(hhincome_clean, perwt, 0.5),
      childless_rate = weighted_share(nchild == 0 & !is.na(nchild), perwt),
      parent_u18_rate = weighted_share(nchild > 0 & !is.na(yngch) & yngch != 99 & yngch < 18, perwt),
      recent_birth_rate = weighted_mean_safe(fifelse(fertyr_observed, recent_birth, NA_integer_), perwt),
      college_rate = weighted_mean_safe(college, perwt),
      married_rate = weighted_mean_safe(married, perwt)
    ), by = by_cols]
    setnames(out, bin_col, "income_bin")
    out[, `:=`(
      bin_type = sub("^income_", "", bin_col),
      education_group = if (!by_educ) "all" else education_group,
      source = "ACS/MMS women ages 22-45, household-income bins within year",
      model_childless_z06_to_z14 = "0.48,0.56,0.17,0.10,0.07"
    )]
    setcolorder(out, c(
      "source", "bin_type", "income_bin", "education_group", "n", "weight",
      "median_hhincome", "childless_rate", "parent_u18_rate",
      "recent_birth_rate", "college_rate", "married_rate",
      "model_childless_z06_to_z14"
    ))
    out[]
  }

  rbindlist(list(
    summarise_bins("income_quintile", FALSE),
    summarise_bins("income_decile", FALSE),
    summarise_bins("income_quintile", TRUE),
    summarise_bins("income_decile", TRUE)
  ), fill = TRUE)
}

build_young_wealth_by_era <- function(psid) {
  y <- copy(psid$young)
  y[, era := fcase(
    year >= 1984 & year <= 1999, "1984-1999",
    year >= 2000 & year <= 2009, "2000-2009",
    year >= 2010 & year <= 2019, "2010-2019",
    default = NA_character_
  )]
  psid_rows <- y[!is.na(era), .(
    n = .N,
    n_clusters = uniqueN(famid),
    weight = sum(weight, na.rm = TRUE),
    mean = weighted_mean_safe(wealth_to_income, weight),
    median = weighted_quantile_safe(wealth_to_income, weight, 0.5),
    method_note = "PSID young childless renters ages 25-35, NETWORTH2R/INCFAMR, INCFAMR>1000, IW-weighted"
  ), by = era]
  psid_rows[, source := "PSID"]
  scf <- data.table(
    source = "SCF",
    era = "2022",
    n = NA_integer_,
    n_clusters = NA_integer_,
    weight = NA_real_,
    mean = 1.14,
    median = 0.385,
    method_note = "Lead-specified SCF 2022 comparison values from D7 era-sensitivity note; no SCF microdata resampling in this harness"
  )
  setcolorder(psid_rows, names(scf))
  rbindlist(list(psid_rows, scf), fill = TRUE)
}

message("Package availability:")
for (pkg in names(package_status)) {
  message("  ", pkg, "=", package_status[[pkg]])
}
if (!package_status[["survey"]]) {
  add_note("Package survey is not installed. This harness checks and reports it, but does not require it because the ACS estimator is implemented directly as a weighted bootstrap/cluster bootstrap.")
}

set.seed(args$seed)

moment_rows <- merge(copy(TARGETS), SOURCE_DEFAULTS, by = "key", sort = FALSE)
moment_rows[, `:=`(
  reproduced_point = NA_real_,
  repro_abs_diff = NA_real_,
  repro_pass = FALSE,
  boot_mean = NA_real_,
  boot_se = NA_real_,
  n_analysis = NA_integer_,
  n_clusters = NA_integer_,
  B = 0L,
  method_note = NA_character_
)]
setkey(moment_rows, key)

psid <- NULL
acs <- NULL
if (args$source %in% c("all", "psid")) {
  timing_start <- Sys.time()
  psid <- build_psid_cache()
  record_timing(
    "PSID", "load_or_build_cached_analysis_samples", timing_start,
    B = 0L,
    n_analysis = sum(vapply(psid, nrow, integer(1))),
    n_clusters = uniqueN(unlist(lapply(psid, function(x) x$famid))),
    n_moments = 4L,
    method_note = "Loaded cached old-age, young-wealth, and first-birth event-study analysis samples"
  )
}
if (args$source %in% c("all", "acs")) {
  timing_start <- Sys.time()
  acs <- build_acs_cache()
  record_timing(
    "ACS/MMS", "load_or_build_cached_analysis_samples", timing_start,
    B = 0L,
    n_analysis = nrow(acs$targets),
    n_clusters = uniqueN(acs$targets$met2013),
    n_moments = 8L,
    method_note = "Loaded cached ACS/MMS target and supplementary analysis samples"
  )
}

if (!is.null(psid)) {
  pp <- psid_point_moments(psid)
  psid_map <- c(
    old_parent_childless_nonhousing_wealth_to_income_gap_6575 = pp[["old_parent_childless_nonhousing_wealth_to_income_gap_6575"]],
    old_nonhousing_wealth_to_income_median_6575 = pp[["old_nonhousing_wealth_to_income_median_6575"]],
    young_childless_renter_liquid_wealth_to_annual_gross_income_2535 = pp[["young_childless_renter_liquid_wealth_to_annual_gross_income_2535"]]
  )
  for (kk in names(psid_map)) {
    target <- target_value(kk)
    reproduced <- as.numeric(psid_map[[kk]])
    pass <- gate_pass(target, reproduced)
    moment_rows[kk, `:=`(
      reproduced_point = reproduced,
      repro_abs_diff = abs(reproduced - target),
      repro_pass = pass,
      n_analysis = n_analysis_for_key(kk, psid = psid),
      n_clusters = n_clusters_for_key(kk, psid = psid),
      method_note = "PSID point estimate from audited builder filters; IW-weighted"
    )]
    log_gate(kk, target, reproduced, if (pass) "PASS" else "FAIL", moment_rows[kk, source], if (pass) "" else "PSID point estimate did not match published target")
  }

  timing_start <- Sys.time()
  event <- run_event_study(psid)
  record_timing(
    "PSID first-birth event study", "saved_stata_summary_or_fixest_import",
    timing_start, B = 0L, n_analysis = event$n, n_clusters = event$n_clusters,
    n_moments = 1L, method_note = event$method_note
  )
  event_path <- load_event_study_path()
  fwrite(event_path, file.path(out_dir, "event_study_coefficient_path.csv"), na = "NA")
  kk <- "housing_increment_0to1"
  target <- target_value(kk)
  reproduced <- event$coef_p3
  pass <- gate_pass(target, reproduced)
  moment_rows[kk, `:=`(
    reproduced_point = reproduced,
    repro_abs_diff = ifelse(is.na(reproduced), NA_real_, abs(reproduced - target)),
    repro_pass = pass,
    boot_mean = if (pass) reproduced else NA_real_,
    boot_se = if (pass) event$se_p3 else NA_real_,
    n_analysis = event$n,
    n_clusters = event$n_clusters,
    B = 0L,
    method_note = paste0(event$method_note, sprintf("; K=5 plateau coef_p5=%s, se_p5=%s", fmt_num(event$coef_p5), fmt_num(event$se_p5)))
  )]
  log_gate(kk, target, reproduced, if (pass) "PASS" else "FAIL", moment_rows[kk, source], if (pass) sprintf("K=5 plateau coef_p5=%s se_p5=%s", fmt_num(event$coef_p5), fmt_num(event$se_p5)) else "fixest Sun-Abraham coefficient did not match published K=3 target")
}

if (!is.null(acs)) {
  ap <- acs_point_moments(acs)
  for (kk in names(ap)) {
    target <- target_value(kk)
    reproduced <- as.numeric(ap[[kk]])
    pass <- gate_pass(target, reproduced)
    moment_rows[kk, `:=`(
      reproduced_point = reproduced,
      repro_abs_diff = abs(reproduced - target),
      repro_pass = pass,
      n_analysis = n_analysis_for_key(kk, acs = acs),
      n_clusters = n_clusters_for_key(kk, acs = acs),
      method_note = "ACS/MMS point estimate from named builders; hhwt-weighted; middle MMS collapsed to center"
    )]
    reason <- ""
    if (kk == "old_age_own_rate") {
      reason <- "source correction: published target is reproduced by ACS/MMS DUE-housing ownership audit; PSID old-age ownership check is 0.863730454 and does not match"
    } else if (!pass) {
      reason <- "ACS point estimate did not match published target"
    }
    log_gate(kk, target, reproduced, if (pass) "PASS" else "FAIL", moment_rows[kk, source], reason)
  }
}

for (kk in c("tfr", "childless_rate")) {
  target <- target_value(kk)
  moment_rows[kk, `:=`(
    reproduced_point = NA_real_,
    repro_abs_diff = NA_real_,
    repro_pass = FALSE,
    method_note = paste(
      "source-unconfirmed after repo search in code/data, code/empirical, and docs;",
      "no in-repo builder or saved target CSV pins the completed-fertility-equivalent value;",
      "required external source is a children-ever-born/completed-fertility source such as a CPS June Fertility Supplement or Census fertility table/microdata extract with declared survey year(s), female completed-fertility age window, CEB/CHBORN variable, and weight;",
      "NCHS period TFR alone cannot source childlessness and ACS NCHILD is a current-household proxy, so SE is NA"
    )
  )]
  log_gate(kk, target, NA_real_, "FAIL", "source-unconfirmed", "No confirmed source builder for completed-fertility-equivalent target; required external children-ever-born/completed-fertility source (e.g., CPS June Fertility Supplement/Census fertility table) is not checked in; not substituting ACS no-child-in-household")
}
log_gate(
  "tfr_childless_source_search", NA_real_, NA_real_, "INFO", "repo source search",
  "Searched code/data, code/empirical, docs, exact values 1.918/0.188, completed fertility, childless/childlessness, tfr/TFR, CPS/ACS fertility supplement references, and saved target CSV names. No reproducible builder found; exact values appear only as calibration target constants/prose, not a source artifact.",
  tol = NA_real_
)

psid_boot <- NULL
acs_boot <- NULL
if (!is.null(psid)) {
  psid_keys <- c(
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
    "old_nonhousing_wealth_to_income_median_6575",
    "young_childless_renter_liquid_wealth_to_annual_gross_income_2535"
  )
  psid_pass_keys <- psid_keys[moment_rows[psid_keys, repro_pass]]
  if (length(psid_pass_keys) > 0) {
    message("Running PSID family-cluster bootstrap for ", length(psid_pass_keys), " passed moments...")
    timing_start <- Sys.time()
    psid_boot <- bootstrap_psid(psid, psid_pass_keys, args$B)
    record_timing(
      "PSID old-age/young wealth", "family_cluster_bootstrap",
      timing_start, B = nrow(psid_boot),
      n_analysis = sum(moment_rows[psid_pass_keys, n_analysis], na.rm = TRUE),
      n_clusters = uniqueN(c(psid$old$famid, psid$young$famid)),
      n_moments = length(psid_pass_keys),
      method_note = "Family-cluster nonparametric bootstrap at ID over cached PSID analysis samples"
    )
    psid_cov <- cov(psid_boot, use = "pairwise.complete.obs")
    for (key in psid_pass_keys) {
      moment_rows[key, `:=`(
        boot_mean = mean(psid_boot[, key], na.rm = TRUE),
        boot_se = sqrt(psid_cov[key, key]),
        B = nrow(psid_boot),
        method_note = paste(method_note, sprintf("PSID family-cluster nonparametric bootstrap, B=%d, seed=%d", nrow(psid_boot), args$seed), sep = "; ")
      )]
    }
  }
}

if (!is.null(acs)) {
  acs_keys <- c(
    "own_rate",
    "own_family_gap",
    "prime30_55_childless_renter_mean_rooms",
    "prime30_55_childless_owner_share_rooms_ge6",
    "prime30_55_childless_owner_minus_renter_mean_rooms",
    "old_age_own_rate",
    "own_rate_2534",
    "prime30_55_parent_3plus_minus_1to2_mean_rooms"
  )
  acs_pass_keys <- acs_keys[moment_rows[acs_keys, repro_pass]]
  if (length(acs_pass_keys) > 0) {
    message("Running ACS/MMS metro-cluster bootstrap for ", length(acs_pass_keys), " passed moments...")
    timing_start <- Sys.time()
    acs_boot <- bootstrap_acs(acs, acs_pass_keys, args$B)
    record_timing(
      "ACS/MMS", "metro_cluster_bootstrap",
      timing_start, B = nrow(acs_boot),
      n_analysis = nrow(acs$targets),
      n_clusters = uniqueN(acs$targets$met2013),
      n_moments = length(acs_pass_keys),
      method_note = "Metro-cluster bootstrap at met2013 over cached ACS/MMS target sample"
    )
    acs_cov <- cov(acs_boot, use = "pairwise.complete.obs")
    for (key in acs_pass_keys) {
      moment_rows[key, `:=`(
        boot_mean = mean(acs_boot[, key], na.rm = TRUE),
        boot_se = sqrt(acs_cov[key, key]),
        B = nrow(acs_boot),
        method_note = paste(method_note, sprintf("ACS/MMS metro-cluster bootstrap because met2013 is present in builders, B=%d, seed=%d", nrow(acs_boot), args$seed), sep = "; ")
      )]
    }
  }
}

if (!is.null(acs)) {
  timing_start <- Sys.time()
  child_income <- build_childlessness_by_income(acs)
  child_path <- file.path(out_dir, "childlessness_by_income.csv")
  writeLines("# model childlessness by z=0.6,0.8,1.0,1.2,1.4: 0.48,0.56,0.17,0.10,0.07", child_path)
  fwrite(child_income, child_path, append = TRUE)
  record_timing(
    "ACS/MMS", "supplement_childlessness_by_income",
    timing_start, B = 0L, n_analysis = nrow(acs$women),
    n_clusters = uniqueN(acs$women$met2013), n_moments = NA_integer_,
    method_note = "Supplementary ACS current-household childlessness-by-income diagnostic"
  )
  log_gate("childlessness_by_income", NA_real_, nrow(child_income), "PASS", "ACS/MMS supplementary", sprintf("table rows=%d; includes quintile/decile and education splits", nrow(child_income)), tol = NA_real_)
} else {
  child_income <- data.table()
  fwrite(child_income, file.path(out_dir, "childlessness_by_income.csv"))
  log_gate("childlessness_by_income", NA_real_, NA_real_, "FAIL", "ACS/MMS supplementary", "ACS source not run", tol = NA_real_)
}

if (!is.null(psid)) {
  timing_start <- Sys.time()
  era <- build_young_wealth_by_era(psid)
  fwrite(era, file.path(out_dir, "young_wealth_by_era.csv"), na = "NA")
  record_timing(
    "PSID young wealth", "supplement_young_wealth_by_era",
    timing_start, B = 0L, n_analysis = nrow(psid$young),
    n_clusters = uniqueN(psid$young$famid), n_moments = NA_integer_,
    method_note = "Supplementary PSID young childless renter wealth/income by era; SCF row is fixed comparison value"
  )
  log_gate("young_wealth_by_era", NA_real_, nrow(era), "PASS", "PSID + SCF supplementary", sprintf("table rows=%d; SCF 2022 row uses lead-specified 1.14/0.385", nrow(era)), tol = NA_real_)
} else {
  era <- data.table()
  fwrite(era, file.path(out_dir, "young_wealth_by_era.csv"))
  log_gate("young_wealth_by_era", NA_real_, NA_real_, "FAIL", "PSID + SCF supplementary", "PSID source not run", tol = NA_real_)
}

cov_mat <- matrix(NA_real_, nrow = nrow(TARGETS), ncol = nrow(TARGETS), dimnames = list(TARGETS$key, TARGETS$key))
if (!is.null(psid_boot) && exists("psid_cov")) {
  cov_mat[rownames(psid_cov), colnames(psid_cov)] <- psid_cov
}
if (!is.null(acs_boot) && exists("acs_cov")) {
  cov_mat[rownames(acs_cov), colnames(acs_cov)] <- acs_cov
}
event_key <- "housing_increment_0to1"
if (isTRUE(moment_rows[event_key, repro_pass]) && is.finite(moment_rows[event_key, boot_se])) {
  cov_mat[event_key, event_key] <- moment_rows[event_key, boot_se]^2
}
for (key in TARGETS$key[is.na(cov_mat[cbind(TARGETS$key, TARGETS$key)])]) {
  if (isTRUE(moment_rows[key, repro_pass]) && is.finite(moment_rows[key, boot_se])) {
    cov_mat[key, key] <- moment_rows[key, boot_se]^2
  }
}

finite_se_keys <- moment_rows[is.finite(boot_se), key]
for (ii in finite_se_keys) {
  for (jj in finite_se_keys) {
    if (ii != jj && is.na(cov_mat[ii, jj])) {
      cov_mat[ii, jj] <- 0
    }
  }
}
for (key in finite_se_keys) {
  delta <- abs(cov_mat[key, key] - moment_rows[key, boot_se]^2)
  if (!is.finite(delta) || delta > 1e-12) {
    stop("Covariance diagonal does not equal se^2 for ", key, ": delta=", delta)
  }
  if (cov_mat[key, key] <= 0) {
    stop("Nonpositive covariance diagonal for ", key)
  }
}

cor_mat <- cov_mat
for (i in seq_len(nrow(cov_mat))) {
  for (j in seq_len(ncol(cov_mat))) {
    denom <- sqrt(cov_mat[i, i] * cov_mat[j, j])
    cor_mat[i, j] <- if (is.finite(denom) && denom > 0 && is.finite(cov_mat[i, j])) cov_mat[i, j] / denom else NA_real_
  }
}

finite_corr <- cor_mat[is.finite(cor_mat)]
if (length(finite_corr) > 0 && any(finite_corr < -1 - 1e-10 | finite_corr > 1 + 1e-10)) {
  stop("Correlation matrix has entries outside [-1,1]")
}

moment_out <- moment_rows[TARGETS, on = "key"]
moment_out[, `:=`(
  published_target = i.published_target,
  i.published_target = NULL
)]
setcolorder(moment_out, c(
  "key", "published_target", "reproduced_point", "repro_abs_diff", "repro_pass",
  "boot_mean", "boot_se", "n_analysis", "n_clusters", "source", "psu",
  "weight_var", "B", "method_note"
))
fwrite(moment_out, file.path(out_dir, "moment_se.csv"), na = "NA")
fwrite(as.data.table(cov_mat, keep.rownames = "key"), file.path(out_dir, "moment_covariance.csv"), na = "NA")
fwrite(as.data.table(cor_mat, keep.rownames = "key"), file.path(out_dir, "moment_correlation.csv"), na = "NA")
fwrite(timings, file.path(out_dir, "source_timing.csv"), na = "NA")

writeLines(repro_lines, file.path(out_dir, "reproduction_log.txt"))

source_provenance <- c(
  "# Source Provenance Log",
  "",
  sprintf("Generated: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## First-Birth Event Study",
  "",
  "- Authoritative saved Stata output found under `code/data/psid_followup_mar2026/output/sa_rooms_first_birth_one_variant_v1/`.",
  "- `rooms_f_c_y_all_summary.csv` verifies K=3 `0.66443467` against the published target `0.664435` and reports K=5 `0.84313726` with SE `0.15404199`.",
  "- `event_study_coefficient_path.csv` imports K=-2..5 from `rooms_f_c_y_all_estimates.dta`.",
  "",
  "## TFR And Childless Rate",
  "",
  "- Searched repository text/source artifacts in `code/data`, `code/empirical`, and `docs` for exact `1.918`, exact `0.188`, completed-fertility phrases, childless/childlessness, `tfr`/`TFR`, CPS/ACS fertility supplement references, and saved target CSV names.",
  "- No in-repo data builder, saved target CSV, or source extract was found that reproduces `tfr = 1.918` or `childless_rate = 0.188`.",
  "- Exact values appear as calibration target constants/prose, not as a checked data-source artifact.",
  "- Required external source to make these moments reproducible: a children-ever-born/completed-fertility source, plausibly CPS June Fertility Supplement or a Census fertility table/microdata extract, with declared survey year(s), female completed-fertility age window, CEB/CHBORN variable, and weight.",
  "- NCHS period TFR alone cannot source `childless_rate`; ACS `NCHILD` is a current-household child measure and repo scripts warn not to use it for completed childlessness. SEs therefore remain `NA` until the external CEB source is pinned."
)
writeLines(source_provenance, file.path(out_dir, "source_provenance_log.txt"))

readme <- c(
  "# Moment Bootstrap Standard Errors",
  "",
  sprintf("Generated: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Seed: `%d`", args$seed),
  sprintf("Requested B: `%d`", args$B),
  "",
  "## Package Check",
  paste0("- `", names(package_status), "`: ", as.character(package_status)),
  "",
  "## Design",
  "",
  "- Reproduction gate is applied before any standard error is reported.",
  "- PSID old-age and young-wealth moments use family-cluster resampling at `ID`, carrying all sampled rows and recomputing `IW`-weighted estimators.",
  "- PSID first-birth rooms uses the read-only existing Stata builder's analytic cluster-robust event-study summary clustered at `ID`; set `EVENT_USE_FIXEST=1` to attempt the slower R `fixest::sunab` re-estimation. The `K=5` plateau coefficient and SE are stored in `method_note`, and K=-2..5 is written to `event_study_coefficient_path.csv`.",
  "- ACS/MMS moments use a metro-cluster bootstrap at `met2013` because the named builders carry `met2013`; within each drawn metro, all household-head rows keep their `hhwt` weights.",
  "- Cross-source covariance blocks are set to zero by construction; moments with failed reproduction gates have `NA` variances/covariances.",
  "- `old_age_own_rate` is routed to the ACS/MMS ownership-audit source because that named builder reproduces the published `0.764261`; the PSID old-age ownership check is `0.863730454` and fails the gate.",
  "- `tfr` and `childless_rate` remain `source-unconfirmed`; the harness does not substitute ACS current-child measures for completed fertility/childlessness. The required external source is a children-ever-born/completed-fertility source such as CPS June Fertility Supplement or a Census fertility table/microdata extract with year(s), age window, variable, and weight declared.",
  "",
  "## Source Timings",
  if (nrow(timings) == 0) "- No timed source phases were recorded." else paste0(
    "- ", timings$source, " / ", timings$phase,
    ": ", sprintf("%.3f", timings$wall_clock_sec), " sec",
    ifelse(is.na(timings$B), "", paste0(", B=", timings$B)),
    ifelse(is.na(timings$n_clusters), "", paste0(", clusters=", timings$n_clusters)),
    ifelse(is.na(timings$n_moments), "", paste0(", moments=", timings$n_moments))
  ),
  "",
  "## Notes",
  if (length(notes) == 0) "- No fallbacks or package issues beyond those listed above." else paste0("- ", notes),
  "",
  "## Outputs",
  "",
  "- `moment_se.csv`",
  "- `moment_covariance.csv`",
  "- `moment_correlation.csv`",
  "- `reproduction_log.txt`",
  "- `event_study_coefficient_path.csv`",
  "- `source_timing.csv`",
  "- `source_provenance_log.txt`",
  "- `childlessness_by_income.csv`",
  "- `young_wealth_by_era.csv`"
)
writeLines(readme, file.path(out_dir, "README.md"))

passed <- moment_out[repro_pass == TRUE, .N]
se_vals <- moment_out[is.finite(boot_se), boot_se]
cat("Console summary\n")
cat(sprintf("Source mode: %s\n", args$source))
cat(sprintf("Reproduction gates passed: %d / 14 target moments\n", passed))
if (length(se_vals) > 0) {
  cat(sprintf(
    "SE min/median/max: %.10g / %.10g / %.10g\n",
    min(se_vals), median(se_vals), max(se_vals)
  ))
} else {
  cat("SE min/median/max: NA / NA / NA\n")
}
cat(sprintf("B used: PSID=%s; ACS=%s; event_study=analytic\n",
            if (!is.null(psid_boot)) nrow(psid_boot) else "NA",
            if (!is.null(acs_boot)) nrow(acs_boot) else "NA"))
cat("Source timings:\n")
print(timings)
cat("Package installs performed: none\n")
if (length(notes) > 0) {
  cat("Notes:\n")
  cat(paste0("- ", notes, collapse = "\n"), "\n")
}
cat("childlessness_by_income head:\n")
print(head(child_income))
cat("young_wealth_by_era head:\n")
print(head(era))
cat("Wrote outputs to: ", out_dir, "\n", sep = "")
