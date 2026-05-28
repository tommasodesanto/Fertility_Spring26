#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(haven)
})

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", file_arg[startsWith(file_arg, file_flag)])
  if (length(script_path) == 0) {
    stop("Run this script with Rscript so --file= is available.")
  }
  dirname(normalizePath(script_path[1]))
}

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(as.numeric(x[ok]), as.numeric(w[ok]))
}

collapse_mms_location <- function(x, middle_target) {
  out <- x
  if (middle_target == "center") {
    out[x == "middle"] <- "center"
  } else if (middle_target == "periphery") {
    out[x == "middle"] <- "periphery"
  } else if (middle_target != "drop") {
    stop("MMS_MIDDLE_TARGET must be one of: drop, center, periphery")
  }
  out
}

extract_python_target <- function(path, name) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  txt <- readLines(path, warn = FALSE)
  pat <- sprintf('"%s"\\s*:\\s*([0-9.]+)', name)
  hit <- regexec(pat, txt)
  vals <- regmatches(txt, hit)
  vals <- vals[lengths(vals) > 1]
  if (length(vals) == 0) {
    return(NA_real_)
  }
  as.numeric(vals[[1]][2])
}

fmt <- function(x) {
  ifelse(is.na(x), "", sprintf("%.3f", x))
}

window_table <- function(dt, sample_name, source_name, windows) {
  rbindlist(lapply(seq_len(nrow(windows)), function(ii) {
    label <- windows$window[ii]
    age_min <- windows$age_min[ii]
    age_max <- windows$age_max[ii]
    sub <- dt[age >= age_min & age <= age_max]
    by_loc <- sub[, .(owner_rate = weighted_mean_safe(owner, weight)), by = mms_location]
    center_rate <- by_loc[mms_location == "center", owner_rate]
    periphery_rate <- by_loc[mms_location == "periphery", owner_rate]
    no_child <- sub[nchild == 0 & !is.na(nchild)]
    parent <- sub[nchild > 0 & !is.na(nchild)]
    newparent <- sub[nchild > 0 & !is.na(nchild) & !is.na(eldch) & eldch != 99 & eldch < 4]

    data.table(
      source = source_name,
      sample = sample_name,
      window = label,
      age_min = age_min,
      age_max = age_max,
      n_records = nrow(sub),
      weight_sum = sum(sub$weight, na.rm = TRUE),
      overall = weighted_mean_safe(sub$owner, sub$weight),
      center = if (length(center_rate) == 0) NA_real_ else center_rate,
      periphery = if (length(periphery_rate) == 0) NA_real_ else periphery_rate,
      p_minus_c = if (length(center_rate) == 0 || length(periphery_rate) == 0) {
        NA_real_
      } else {
        periphery_rate - center_rate
      },
      parent_minus_nochildren = weighted_mean_safe(parent$owner, parent$weight) -
        weighted_mean_safe(no_child$owner, no_child$weight),
      newparent_minus_nochildren = weighted_mean_safe(newparent$owner, newparent$weight) -
        weighted_mean_safe(no_child$owner, no_child$weight)
    )
  }), fill = TRUE)
}

psid_window_table <- function(dt, sample_name, windows) {
  rbindlist(lapply(seq_len(nrow(windows)), function(ii) {
    label <- windows$window[ii]
    age_min <- windows$age_min[ii]
    age_max <- windows$age_max[ii]
    sub <- dt[age >= age_min & age <= age_max]
    no_child <- sub[children == 0 & !is.na(children)]
    parent <- sub[children > 0 & !is.na(children)]

    data.table(
      source = "PSID",
      sample = sample_name,
      window = label,
      age_min = age_min,
      age_max = age_max,
      n_records = nrow(sub),
      weight_sum = sum(sub$weight, na.rm = TRUE),
      overall = weighted_mean_safe(sub$owner, sub$weight),
      center = NA_real_,
      periphery = NA_real_,
      p_minus_c = NA_real_,
      parent_minus_nochildren = weighted_mean_safe(parent$owner, parent$weight) -
        weighted_mean_safe(no_child$owner, no_child$weight),
      newparent_minus_nochildren = NA_real_
    )
  }), fill = TRUE)
}

age_profile <- function(dt, sample_name, source_name) {
  dt[age >= 20 & age <= 84, .(
    source = source_name,
    sample = sample_name,
    n_records = .N,
    weight_sum = sum(weight, na.rm = TRUE),
    owner_rate = weighted_mean_safe(owner, weight)
  ), by = age][order(source, sample, age)]
}

age_location_profile <- function(dt, sample_name, source_name) {
  dt[age >= 20 & age <= 84, .(
    source = source_name,
    sample = sample_name,
    n_records = .N,
    weight_sum = sum(weight, na.rm = TRUE),
    owner_rate = weighted_mean_safe(owner, weight)
  ), by = .(age, mms_location)][order(source, sample, mms_location, age)]
}

append_metric <- function(rows, moment, label, table, sample_name, value_col, out_col) {
  val <- table[sample == sample_name & window == label, get(value_col)]
  if (length(val) == 0) {
    return(rows)
  }
  rows[moment == moment, (out_col) := val[1]]
  rows
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_dir <- file.path(script_dir, Sys.getenv("OWNERSHIP_AUDIT_OUTPUT_DIR", "output_ownership_audit"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

middle_target <- tolower(Sys.getenv("MMS_MIDDLE_TARGET", ""))
if (middle_target == "") {
  include_middle_in_center <- tolower(Sys.getenv("MMS_INCLUDE_MIDDLE_IN_CENTER", "true")) %in%
    c("1", "true", "yes", "y")
  middle_target <- if (include_middle_in_center) "center" else "drop"
}
if (!(middle_target %in% c("drop", "center", "periphery"))) {
  stop("MMS_MIDDLE_TARGET must be one of: drop, center, periphery")
}

lookup_2010 <- as.data.table(fread(file.path(data_dir, "puma_mms_lookup_2010.csv")))
lookup_2020 <- as.data.table(fread(file.path(data_dir, "puma_mms_lookup_2020.csv")))
dest_lookup <- rbindlist(list(
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

windows <- data.table(
  window = c("20_84", "20_24", "22_45", "25_34", "30_45", "30_55", "35_44", "65_75", "80_84"),
  age_min = c(20L, 20L, 22L, 25L, 30L, 30L, 35L, 65L, 80L),
  age_max = c(84L, 24L, 45L, 34L, 45L, 55L, 44L, 75L, 84L)
)

message("Loading ACS/IPUMS extract27 for ownership audit...")
extract_path <- file.path(repo_root, "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta")
acs <- as.data.table(read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "met2013", "gq", "pernum", "relate",
    "hhwt", "perwt", "ownershp", "unitsstr", "rooms", "nchild", "eldch", "age"
  )
))
acs <- acs[
  year >= 2012 & year <= 2023 &
    met2013 > 0 &
    gq %in% c(1, 2) &
    age >= 20 & age <= 84 &
    ownershp %in% c(1, 2)
]
acs[, `:=`(
  year = as.integer(year),
  statefip = as.integer(statefip),
  puma = as.integer(puma),
  met2013 = as.integer(met2013),
  gq = as.integer(gq),
  pernum = as.integer(pernum),
  relate = as.integer(relate),
  hhwt = as.numeric(hhwt),
  perwt = as.numeric(perwt),
  ownershp = as.integer(ownershp),
  unitsstr = as.integer(unitsstr),
  rooms = as.numeric(rooms),
  nchild = as.numeric(nchild),
  eldch = as.numeric(eldch),
  age = as.integer(age),
  lookup_period = fifelse(year <= 2021, "pre2022", "post2021"),
  owner = as.integer(ownershp == 1L)
)]

head_equivalence <- acs[, .(
  n_records = .N,
  perwt_sum = sum(perwt, na.rm = TRUE),
  hhwt_sum = sum(hhwt, na.rm = TRUE)
), by = .(
  pernum_eq_1 = pernum == 1L,
  relate_eq_1 = relate == 1L
)][order(-pernum_eq_1, -relate_eq_1)]

acs <- merge(
  acs,
  dest_lookup,
  by = c("lookup_period", "statefip", "puma", "met2013"),
  all.x = TRUE,
  sort = FALSE
)
acs[, mms_location := collapse_mms_location(mms_location, middle_target)]
acs <- acs[mms_location %in% c("center", "periphery")]

acs_person <- copy(acs[!is.na(perwt) & perwt > 0])
acs_person[, weight := perwt]

acs_head_all <- copy(acs[pernum == 1L & relate == 1L & !is.na(hhwt) & hhwt > 0])
acs_head_all[, weight := hhwt]

acs_head_due <- copy(acs_head_all[
  unitsstr %in% 3:10 &
    !is.na(rooms) &
    rooms > 0
])

acs_windows <- rbindlist(list(
  window_table(acs_person, "all_persons_perwt", "ACS", windows),
  window_table(acs_head_all, "household_heads_hhwt_all_housing", "ACS", windows),
  window_table(acs_head_due, "household_heads_hhwt_due_housing", "ACS", windows)
), fill = TRUE)

acs_profiles <- rbindlist(list(
  age_profile(acs_person, "all_persons_perwt", "ACS"),
  age_profile(acs_head_all, "household_heads_hhwt_all_housing", "ACS"),
  age_profile(acs_head_due, "household_heads_hhwt_due_housing", "ACS")
), fill = TRUE)

acs_location_profiles <- rbindlist(list(
  age_location_profile(acs_person, "all_persons_perwt", "ACS"),
  age_location_profile(acs_head_all, "household_heads_hhwt_all_housing", "ACS"),
  age_location_profile(acs_head_due, "household_heads_hhwt_due_housing", "ACS")
), fill = TRUE)

message("Loading PSID shelf for ownership cross-check...")
psid_path <- "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
psid_raw <- as.data.table(read_dta(
  psid_path,
  col_select = c(
    "ID", "year", "CURRENT", "AGEREP", "DEATHYEAR", "REL",
    "HOMEOWN", "IW", "RELCHINUM", "RELCHIREP"
  )
))
psid_raw <- psid_raw[
  year >= 1984 & year <= 2019 &
    CURRENT == 1 &
    (is.na(DEATHYEAR) | year <= DEATHYEAR) &
    !is.na(AGEREP) &
    AGEREP >= 20 & AGEREP <= 84 &
    !is.na(IW) & IW > 0
]
psid_raw[, `:=`(
  age = as.integer(AGEREP),
  weight = as.numeric(IW),
  owner = fifelse(HOMEOWN == 1, 1, fifelse(HOMEOWN == 2, 0, NA_real_)),
  children = fifelse(!is.na(RELCHINUM), as.numeric(RELCHINUM), as.numeric(RELCHIREP)),
  rel = as.integer(REL)
)]
psid_raw <- psid_raw[!is.na(owner)]

psid_period_tables <- function(dt, year_min, year_max, suffix) {
  sub <- dt[year >= year_min & year <= year_max]
  rbindlist(list(
    psid_window_table(sub, paste0("all_individuals_iw_", suffix), windows),
    psid_window_table(sub[rel == 1L], paste0("reference_person_iw_", suffix), windows),
    psid_window_table(sub[rel %in% c(1L, 2L)], paste0("reference_or_partner_iw_", suffix), windows)
  ), fill = TRUE)
}

psid_period_profiles <- function(dt, year_min, year_max, suffix) {
  sub <- dt[year >= year_min & year <= year_max]
  rbindlist(list(
    age_profile(sub, paste0("all_individuals_iw_", suffix), "PSID"),
    age_profile(sub[rel == 1L], paste0("reference_person_iw_", suffix), "PSID"),
    age_profile(sub[rel %in% c(1L, 2L)], paste0("reference_or_partner_iw_", suffix), "PSID")
  ), fill = TRUE)
}

psid_windows <- rbindlist(list(
  psid_period_tables(psid_raw, 1984, 2019, "1984_2019"),
  psid_period_tables(psid_raw, 2000, 2019, "2000_2019"),
  psid_period_tables(psid_raw, 2012, 2019, "2012_2019")
), fill = TRUE)

psid_profiles <- rbindlist(list(
  psid_period_profiles(psid_raw, 1984, 2019, "1984_2019"),
  psid_period_profiles(psid_raw, 2000, 2019, "2000_2019"),
  psid_period_profiles(psid_raw, 2012, 2019, "2012_2019")
), fill = TRUE)

all_windows <- rbindlist(list(acs_windows, psid_windows), fill = TRUE)
all_profiles <- rbindlist(list(acs_profiles, psid_profiles), fill = TRUE)

fwrite(acs_windows, file.path(out_dir, "acs_ownership_window_targets.csv"))
fwrite(acs_profiles, file.path(out_dir, "acs_ownership_age_profiles.csv"))
fwrite(acs_location_profiles, file.path(out_dir, "acs_ownership_age_location_profiles.csv"))
fwrite(psid_windows, file.path(out_dir, "psid_ownership_window_targets.csv"))
fwrite(psid_profiles, file.path(out_dir, "psid_ownership_age_profiles.csv"))
fwrite(all_windows, file.path(out_dir, "ownership_window_targets_all_sources.csv"))
fwrite(all_profiles, file.path(out_dir, "ownership_age_profiles_all_sources.csv"))
fwrite(head_equivalence, file.path(out_dir, "acs_head_equivalence_diagnostic.csv"))

plot_profiles <- all_profiles[
  (source == "ACS" & sample %in% c("all_persons_perwt", "household_heads_hhwt_due_housing")) |
    (source == "PSID" & sample %in% c("reference_person_iw_2012_2019", "reference_person_iw_1984_2019"))
]
plot_profiles[, series := fcase(
  source == "ACS" & sample == "all_persons_perwt",
  "ACS all persons",
  source == "ACS" & sample == "household_heads_hhwt_due_housing",
  "ACS household heads",
  source == "PSID" & sample == "reference_person_iw_2012_2019",
  "PSID ref. person, 2012-2019",
  source == "PSID" & sample == "reference_person_iw_1984_2019",
  "PSID ref. person, 1984-2019",
  default = NA_character_
)]
plot_profiles <- plot_profiles[!is.na(series)]
plot_profiles[, series := factor(
  series,
  levels = c(
    "ACS all persons",
    "ACS household heads",
    "PSID ref. person, 2012-2019",
    "PSID ref. person, 1984-2019"
  )
)]

lifecycle_plot <- ggplot(plot_profiles, aes(x = age, y = owner_rate, color = series, linetype = series)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = function(x) sprintf("%d%%", round(100 * x))
  ) +
  scale_x_continuous(breaks = seq(20, 85, 5), limits = c(20, 84)) +
  scale_color_manual(values = c(
    "ACS all persons" = "#8C2D04",
    "ACS household heads" = "#08519C",
    "PSID ref. person, 2012-2019" = "#238B45",
    "PSID ref. person, 1984-2019" = "#636363"
  )) +
  scale_linetype_manual(values = c(
    "ACS all persons" = "solid",
    "ACS household heads" = "solid",
    "PSID ref. person, 2012-2019" = "dashed",
    "PSID ref. person, 1984-2019" = "dotted"
  )) +
  labs(
    x = "Age",
    y = "Homeownership rate",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 14, 8, 10)
  )

ggsave(file.path(out_dir, "ownership_lifecycle_acs_psid.png"), lifecycle_plot, width = 8.5, height = 5.2, dpi = 300)
ggsave(file.path(out_dir, "ownership_lifecycle_acs_psid.pdf"), lifecycle_plot, width = 8.5, height = 5.2)

param_path <- file.path(repo_root, "code/model/dt_cp_model/parameters.py")
current_targets <- data.table(
  moment = c(
    "own_rate_30_55",
    "own_gradient_30_55",
    "newparent_minus_nochildren_30_55",
    "old_age_own_rate_65_75",
    "old_age_parent_childless_gap_65_75"
  ),
  current_calibration_target = c(
    extract_python_target(param_path, "own_rate"),
    extract_python_target(param_path, "own_gradient"),
    extract_python_target(param_path, "own_family_gap"),
    extract_python_target(param_path, "old_age_own_rate"),
    extract_python_target(param_path, "old_age_parent_childless_gap")
  )
)

comparison <- data.table(moment = c(
  "own_rate_20_84",
  "own_rate_20_24",
  "own_rate_22_45",
  "own_rate_25_34",
  "own_rate_30_55",
  "own_rate_35_44",
  "own_rate_65_75",
  "own_rate_80_84",
  "own_gradient_30_55",
  "parent_minus_nochildren_30_55",
  "newparent_minus_nochildren_30_55",
  "old_age_parent_childless_gap_65_75",
  "lifecycle_slope_80_84_minus_20_24",
  "lifecycle_slope_65_75_minus_25_34"
))
comparison <- merge(comparison, current_targets, by = "moment", all.x = TRUE, sort = FALSE)

add_source_value <- function(comp, out_col, table, sample_name) {
  get_val <- function(moment) {
    if (startsWith(moment, "own_rate_")) {
      label <- sub("^own_rate_", "", moment)
      return(table[sample == sample_name & window == label, overall][1])
    }
    if (moment == "own_gradient_30_55") {
      return(table[sample == sample_name & window == "30_55", p_minus_c][1])
    }
    if (moment == "parent_minus_nochildren_30_55") {
      return(table[sample == sample_name & window == "30_55", parent_minus_nochildren][1])
    }
    if (moment == "newparent_minus_nochildren_30_55") {
      return(table[sample == sample_name & window == "30_55", newparent_minus_nochildren][1])
    }
    if (moment == "old_age_parent_childless_gap_65_75") {
      return(table[sample == sample_name & window == "65_75", parent_minus_nochildren][1])
    }
    if (moment == "lifecycle_slope_80_84_minus_20_24") {
      old <- table[sample == sample_name & window == "80_84", overall][1]
      young <- table[sample == sample_name & window == "20_24", overall][1]
      return(old - young)
    }
    if (moment == "lifecycle_slope_65_75_minus_25_34") {
      old <- table[sample == sample_name & window == "65_75", overall][1]
      young <- table[sample == sample_name & window == "25_34", overall][1]
      return(old - young)
    }
    NA_real_
  }
  comp[, (out_col) := vapply(moment, get_val, numeric(1))]
  comp
}

comparison <- add_source_value(comparison, "acs_person_perwt", acs_windows, "all_persons_perwt")
comparison <- add_source_value(
  comparison,
  "acs_head_hhwt_all_housing",
  acs_windows,
  "household_heads_hhwt_all_housing"
)
comparison <- add_source_value(
  comparison,
  "acs_head_hhwt_due_housing",
  acs_windows,
  "household_heads_hhwt_due_housing"
)
comparison <- add_source_value(
  comparison,
  "psid_recent_all_individuals",
  psid_windows,
  "all_individuals_iw_2012_2019"
)
comparison <- add_source_value(
  comparison,
  "psid_recent_reference_person",
  psid_windows,
  "reference_person_iw_2012_2019"
)
comparison <- add_source_value(
  comparison,
  "psid_recent_reference_or_partner",
  psid_windows,
  "reference_or_partner_iw_2012_2019"
)
comparison <- add_source_value(
  comparison,
  "psid_long_reference_person",
  psid_windows,
  "reference_person_iw_1984_2019"
)

comparison[, due_reference := NA_real_]
comparison[moment == "own_rate_20_84", due_reference := 0.58]
comparison[moment == "lifecycle_slope_80_84_minus_20_24", due_reference := 0.64]

fwrite(comparison, file.path(out_dir, "ownership_target_comparison.csv"))

md_path <- file.path(out_dir, "ACS_MMS_OWNERSHIP_TARGET_AUDIT.md")
head_diag <- fread(file.path(out_dir, "acs_head_equivalence_diagnostic.csv"))
key <- comparison[moment %in% c(
  "own_rate_20_24",
  "own_rate_30_55",
  "own_gradient_30_55",
  "newparent_minus_nochildren_30_55",
  "old_age_own_rate_65_75",
  "lifecycle_slope_80_84_minus_20_24"
)]

md <- c(
  "# ACS/MMS Ownership Target Audit",
  "",
  sprintf("Generated: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Purpose: compare the current person-level ACS ownership construction with a household-head ownership object that matches the model and the DUE-style lifecycle target.",
  "",
  "## Sample Definitions",
  "",
  "- `ACS/all_persons_perwt`: person records, person weights, `OWNERSHP` repeated from the household record.",
  "- `ACS/household_heads_hhwt_all_housing`: `PERNUM == 1` and `RELATE == 1`, household weights.",
  "- `ACS/household_heads_hhwt_due_housing`: household heads, household weights, `UNITSSTR` in standard structures `3:10`, and positive rooms.",
  "- `PSID/reference_person_iw`: `REL == 1`, individual weights, `HOMEOWN == 1` as owner and `HOMEOWN == 2` as renter.",
  "",
  "## Key Moments",
  "",
  "| moment | current target | ACS person | ACS head | ACS head DUE-housing | PSID recent ref person | PSID long ref person | DUE reference |",
  "|---|---:|---:|---:|---:|---:|---:|---:|"
)
for (ii in seq_len(nrow(key))) {
  md <- c(md, sprintf(
    "| `%s` | %s | %s | %s | %s | %s | %s | %s |",
    key$moment[ii],
    fmt(key$current_calibration_target[ii]),
    fmt(key$acs_person_perwt[ii]),
    fmt(key$acs_head_hhwt_all_housing[ii]),
    fmt(key$acs_head_hhwt_due_housing[ii]),
    fmt(key$psid_recent_reference_person[ii]),
    fmt(key$psid_long_reference_person[ii]),
    fmt(key$due_reference[ii])
  ))
}
md <- c(
  md,
  "",
  "## Head Diagnostic",
  "",
  "This checks whether `PERNUM == 1` and `RELATE == 1` agree in the ACS extract before the MMS geography merge.",
  "",
  "| pernum == 1 | relate == 1 | records | person-weight sum | household-weight sum |",
  "|---:|---:|---:|---:|---:|"
)
for (ii in seq_len(nrow(head_diag))) {
  md <- c(md, sprintf(
    "| %s | %s | %s | %.0f | %.0f |",
    head_diag$pernum_eq_1[ii],
    head_diag$relate_eq_1[ii],
    format(head_diag$n_records[ii], big.mark = ",", scientific = FALSE),
    head_diag$perwt_sum[ii],
    head_diag$hhwt_sum[ii]
  ))
}
md <- c(
  md,
  "",
  "## Read",
  "",
  "The person-level ACS object is a share of people living in owner-occupied units. For young adults it includes adult children in parent-owned homes. The household-head object is the credible model target because the model has independent households and no co-residence state.",
  "",
  "The recent PSID columns use 2012--2019, the closest PSID period to the ACS target window. The long PSID column uses 1984--2019 and is useful for old-age checks but is not the cleanest cross-check for the current ACS ownership level.",
  "",
  "The DUE reference values are included for conceptual orientation only; their geography and period are not identical to the MMS target system.",
  "",
  "## Output Files",
  "",
  "- `ownership_target_comparison.csv`",
  "- `acs_ownership_window_targets.csv`",
  "- `acs_ownership_age_profiles.csv`",
  "- `acs_ownership_age_location_profiles.csv`",
  "- `psid_ownership_window_targets.csv`",
  "- `psid_ownership_age_profiles.csv`",
  "- `acs_head_equivalence_diagnostic.csv`",
  "- `ownership_lifecycle_acs_psid.png`",
  "- `ownership_lifecycle_acs_psid.pdf`"
)
writeLines(md, md_path)

cat("Wrote ownership audit packet to:\n")
cat(out_dir, "\n")
cat("\nKey comparison:\n")
print(key)
