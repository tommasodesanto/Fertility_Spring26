#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(haven)
})

SCF_SUMMARY_URL <- "https://www.federalreserve.gov/econres/files/scfp2022s.zip"

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

weighted_quantile_safe <- function(x, w, prob = 0.5) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  x <- as.numeric(x[ok])
  w <- as.numeric(w[ok])
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= prob)[1]]
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

clean_home_value <- function(x) {
  out <- as.numeric(x)
  out[!is.finite(out) | out <= 0 | out >= 9999998] <- NA_real_
  out
}

age_group <- function(age) {
  fcase(
    age >= 20 & age <= 34, "20_34",
    age >= 35 & age <= 44, "35_44",
    age >= 45 & age <= 54, "45_54",
    age >= 55 & age <= 64, "55_64",
    age >= 65 & age <= 74, "65_74",
    age >= 75 & age <= 84, "75_84",
    default = NA_character_
  )
}

summarise_home_value <- function(dt, by_cols) {
  dt[, .(
    n_records = .N,
    weight_sum = sum(weight, na.rm = TRUE),
    owner_rate = weighted_mean_safe(owner, weight),
    mean_owner_occupied_value = weighted_mean_safe(owner_occupied_value, weight),
    mean_home_value_owner = weighted_mean_safe(home_value_owner, weight),
    median_home_value_owner = weighted_quantile_safe(home_value_owner, weight, 0.5)
  ), by = by_cols]
}

read_scf_summary <- function(zip_path) {
  if (!file.exists(zip_path)) {
    allow_download <- tolower(Sys.getenv("SCF_ALLOW_DOWNLOAD", "0")) %in%
      c("1", "true", "yes", "y")
    if (!allow_download) {
      stop(
        "Missing SCF summary zip: ", zip_path, "\n",
        "Download it with:\n",
        "  curl -L ", SCF_SUMMARY_URL, " -o ", zip_path, "\n",
        "or rerun with SCF_ALLOW_DOWNLOAD=1."
      )
    }
    dir.create(dirname(zip_path), recursive = TRUE, showWarnings = FALSE)
    download.file(SCF_SUMMARY_URL, zip_path, mode = "wb", quiet = FALSE)
  }
  td <- tempfile("scf_summary_")
  dir.create(td)
  unzip(zip_path, exdir = td)
  dta <- list.files(td, pattern = "[.]dta$", full.names = TRUE, ignore.case = TRUE)
  if (length(dta) == 0) {
    stop("No Stata .dta file found inside ", zip_path)
  }
  as.data.table(read_dta(
    dta[1],
    col_select = c("yy1", "y1", "wgt", "age", "housecl", "houses", "networth")
  ))
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
out_dir <- file.path(script_dir, Sys.getenv(
  "HOME_VALUE_VALIDATION_OUTPUT_DIR",
  "output_housing_value_validation"
))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scf_zip <- Sys.getenv("SCF_SUMMARY_ZIP", file.path(script_dir, "data_scf", "scfp2022s.zip"))

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

message("Loading ACS/IPUMS extract27 with VALUEH...")
extract_path <- file.path(repo_root, "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta")
acs <- as.data.table(read_dta(
  extract_path,
  col_select = c(
    "year", "statefip", "puma", "met2013", "gq", "pernum", "relate",
    "hhwt", "ownershp", "unitsstr", "rooms", "age", "valueh"
  )
))
acs <- acs[
  year >= 2012 & year <= 2023 &
    gq %in% c(1, 2) &
    pernum == 1 &
    relate == 1 &
    age >= 20 & age <= 84 &
    ownershp %in% c(1, 2) &
    !is.na(hhwt) & hhwt > 0
]
acs[, `:=`(
  year = as.integer(year),
  statefip = as.integer(statefip),
  puma = as.integer(puma),
  met2013 = as.integer(met2013),
  gq = as.integer(gq),
  ownershp = as.integer(ownershp),
  unitsstr = as.integer(unitsstr),
  rooms = as.numeric(rooms),
  age = as.integer(age),
  weight = as.numeric(hhwt),
  lookup_period = fifelse(year <= 2021, "pre2022", "post2021"),
  owner = as.integer(ownershp == 1L),
  home_value_owner = fifelse(ownershp == 1L, clean_home_value(valueh), NA_real_)
)]
acs[, owner_occupied_value := fifelse(owner == 1L, home_value_owner, 0)]
acs[, age_group := age_group(age)]

acs_due <- copy(acs[
  unitsstr %in% 3:10 &
    !is.na(rooms) &
    rooms > 0
])

acs_mms <- merge(
  acs_due[met2013 > 0],
  dest_lookup,
  by = c("lookup_period", "statefip", "puma", "met2013"),
  all.x = TRUE,
  sort = FALSE
)
acs_mms[, mms_location := collapse_mms_location(mms_location, middle_target)]
acs_mms <- acs_mms[mms_location %in% c("center", "periphery")]

spatial_age_profiles <- summarise_home_value(
  acs_mms,
  c("age", "mms_location")
)[order(mms_location, age)]
spatial_age_profiles[, `:=`(
  source = "ACS",
  sample = "household_heads_hhwt_due_housing"
)]
setcolorder(spatial_age_profiles, c(
  "age", "mms_location", "source", "sample", "n_records", "weight_sum",
  "owner_rate", "mean_owner_occupied_value", "mean_home_value_owner",
  "median_home_value_owner"
))

acs_validation <- rbindlist(list(
  summarise_home_value(
    acs[year == 2022 & !is.na(age_group)],
    c("age_group")
  )[, `:=`(source = "ACS", sample = "national_heads_all_housing_2022")],
  summarise_home_value(
    acs_due[year == 2022 & !is.na(age_group)],
    c("age_group")
  )[, `:=`(source = "ACS", sample = "national_heads_due_housing_2022")],
  summarise_home_value(
    acs_mms[year == 2022 & !is.na(age_group)],
    c("age_group")
  )[, `:=`(source = "ACS", sample = "mms_heads_due_housing_2022")]
), fill = TRUE)

message("Loading SCF 2022 summary extract...")
scf <- read_scf_summary(scf_zip)
scf[, `:=`(
  implicate = as.integer(y1 %% 10),
  age = as.integer(age),
  weight = as.numeric(wgt),
  owner = as.integer(housecl == 1L),
  home_value_owner = fifelse(housecl == 1L & houses > 0, as.numeric(houses), NA_real_),
  owner_occupied_value = fifelse(housecl == 1L & houses > 0, as.numeric(houses), 0)
)]
scf <- scf[age >= 20 & age <= 84 & weight > 0]
scf[, age_group := age_group(age)]

scf_implicate <- summarise_home_value(
  scf[!is.na(age_group)],
  c("implicate", "age_group")
)
scf_validation <- scf_implicate[, .(
  n_records = mean(n_records, na.rm = TRUE),
  weight_sum = mean(weight_sum, na.rm = TRUE),
  owner_rate = mean(owner_rate, na.rm = TRUE),
  mean_owner_occupied_value = mean(mean_owner_occupied_value, na.rm = TRUE),
  mean_home_value_owner = mean(mean_home_value_owner, na.rm = TRUE),
  median_home_value_owner = mean(median_home_value_owner, na.rm = TRUE)
), by = age_group]
scf_validation[, `:=`(source = "SCF", sample = "summary_extract_2022")]

validation <- rbindlist(list(acs_validation, scf_validation), fill = TRUE)
setcolorder(validation, c(
  "source", "sample", "age_group", "n_records", "weight_sum", "owner_rate",
  "mean_owner_occupied_value", "mean_home_value_owner",
  "median_home_value_owner"
))
setorder(validation, age_group, source, sample)

wide_validation <- dcast(
  validation[sample %in% c("national_heads_due_housing_2022", "summary_extract_2022")],
  age_group ~ source,
  value.var = c("owner_rate", "mean_owner_occupied_value", "mean_home_value_owner", "median_home_value_owner")
)
for (metric in c("owner_rate", "mean_owner_occupied_value", "mean_home_value_owner", "median_home_value_owner")) {
  acs_col <- paste0(metric, "_ACS")
  scf_col <- paste0(metric, "_SCF")
  if (acs_col %in% names(wide_validation) && scf_col %in% names(wide_validation)) {
    wide_validation[, (paste0(metric, "_acs_minus_scf")) := get(acs_col) - get(scf_col)]
    wide_validation[, (paste0(metric, "_acs_over_scf")) := get(acs_col) / get(scf_col)]
  }
}

fwrite(spatial_age_profiles, file.path(out_dir, "acs_home_value_age_location_profiles.csv"))
fwrite(validation, file.path(out_dir, "acs_scf_home_value_validation.csv"))
fwrite(wide_validation, file.path(out_dir, "acs_scf_home_value_validation_wide.csv"))

validation_plot <- validation[
  sample %in% c("national_heads_due_housing_2022", "summary_extract_2022")
]
validation_plot[, series := fifelse(source == "ACS", "ACS heads, DUE housing", "SCF summary extract")]
validation_plot[, age_group := factor(
  age_group,
  levels = c("20_34", "35_44", "45_54", "55_64", "65_74", "75_84")
)]

plot_long <- melt(
  validation_plot,
  id.vars = c("series", "age_group"),
  measure.vars = c("owner_rate", "mean_owner_occupied_value", "mean_home_value_owner"),
  variable.name = "metric",
  value.name = "value"
)
plot_long[, metric_label := fcase(
  metric == "owner_rate", "Ownership rate",
  metric == "mean_owner_occupied_value", "E[owner x primary residence value]",
  metric == "mean_home_value_owner", "E[primary residence value | owner]",
  default = metric
)]
plot_long[metric != "owner_rate", value := value / 1000]

home_value_validation_plot <- ggplot(plot_long, aes(x = age_group, y = value, color = series, group = series)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  facet_wrap(
    ~metric_label,
    scales = "free_y",
    ncol = 1,
    labeller = label_wrap_gen(width = 36)
  ) +
  scale_color_manual(values = c(
    "ACS heads, DUE housing" = "#08519C",
    "SCF summary extract" = "#8C2D04"
  )) +
  labs(
    x = "Reference-person age group",
    y = "Rate, or 2022 dollars (thousands)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(
  file.path(out_dir, "acs_scf_home_value_validation.png"),
  home_value_validation_plot,
  width = 8.2,
  height = 8.6,
  dpi = 300
)
ggsave(
  file.path(out_dir, "acs_scf_home_value_validation.pdf"),
  home_value_validation_plot,
  width = 8.2,
  height = 8.6
)

spatial_plot <- copy(spatial_age_profiles[age >= 20 & age <= 55])
spatial_plot[, mms_location := factor(mms_location, levels = c("center", "periphery"))]
spatial_plot_long <- melt(
  spatial_plot,
  id.vars = c("age", "mms_location"),
  measure.vars = c("mean_owner_occupied_value", "mean_home_value_owner"),
  variable.name = "metric",
  value.name = "value"
)
spatial_plot_long[, `:=`(
  metric_label = fifelse(
    metric == "mean_owner_occupied_value",
    "E[owner x VALUEH]",
    "E[VALUEH | owner]"
  ),
  value = value / 1000
)]

spatial_value_plot <- ggplot(spatial_plot_long, aes(x = age, y = value, color = mms_location)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~metric_label, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("center" = "#08519C", "periphery" = "#8C2D04")) +
  labs(
    x = "Age",
    y = "2022 dollars (thousands)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(
  file.path(out_dir, "acs_home_value_spatial_lifecycle.png"),
  spatial_value_plot,
  width = 8.2,
  height = 6.4,
  dpi = 300
)
ggsave(
  file.path(out_dir, "acs_home_value_spatial_lifecycle.pdf"),
  spatial_value_plot,
  width = 8.2,
  height = 6.4
)

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", format(round(x, digits), big.mark = ",", scientific = FALSE))
}
key <- wide_validation[age_group %in% c("20_34", "35_44", "45_54", "65_74")]

md <- c(
  "# ACS Home-Value / SCF Validation",
  "",
  sprintf("Generated: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Purpose: build the ACS spatial housing-value lifecycle object and validate the ACS `VALUEH` level against the public 2022 SCF summary extract.",
  "",
  "## Definitions",
  "",
  "- ACS spatial object: MMS household heads, household weights, DUE-style housing restrictions, center/periphery geography, ages 20--84.",
  "- ACS validation object: national ACS 2022 household heads with the same DUE-style housing restrictions.",
  "- SCF validation object: 2022 SCF summary extract; point estimates average the five implicates.",
  "- Housing-value stock: `E[owner * VALUEH]` in ACS and `E[owner * HOUSES]` in SCF.",
  "- Conditional owner value: `E[VALUEH | owner]` in ACS and `E[HOUSES | owner]` in SCF.",
  "",
  "The SCF is national in the public data used here. It validates the ACS home-value measurement, but it is not a center/periphery counterpart.",
  "",
  "## National Validation",
  "",
  "| age group | ACS own | SCF own | ACS E[O x value] | SCF E[O x value] | ACS E[value|O] | SCF E[value|O] |",
  "|---|---:|---:|---:|---:|---:|---:|"
)
for (ii in seq_len(nrow(key))) {
  md <- c(md, sprintf(
    "| `%s` | %.3f | %.3f | %s | %s | %s | %s |",
    key$age_group[ii],
    key$owner_rate_ACS[ii],
    key$owner_rate_SCF[ii],
    fmt_num(key$mean_owner_occupied_value_ACS[ii], 0),
    fmt_num(key$mean_owner_occupied_value_SCF[ii], 0),
    fmt_num(key$mean_home_value_owner_ACS[ii], 0),
    fmt_num(key$mean_home_value_owner_SCF[ii], 0)
  ))
}
md <- c(
  md,
  "",
  "## Read",
  "",
  "ACS and SCF ownership levels should be close but not identical because ACS is a much larger housing survey sample while SCF is optimized for balance sheets and wealth tails.",
  "",
  "The conditional-owner value is the cleanest direct check on `VALUEH` versus SCF `HOUSES`. The unconditional housing-value stock also embeds ownership-rate differences.",
  "",
  "For the spatial model diagnostic, use `E[owner * VALUEH | age, location]`. It maps to the model object `E[1{owner} p_i H | age, location]`; it is not total net worth or home equity.",
  "",
  "## Output Files",
  "",
  "- `acs_home_value_age_location_profiles.csv`",
  "- `acs_scf_home_value_validation.csv`",
  "- `acs_scf_home_value_validation_wide.csv`",
  "- `acs_scf_home_value_validation.png`",
  "- `acs_scf_home_value_validation.pdf`",
  "- `acs_home_value_spatial_lifecycle.png`",
  "- `acs_home_value_spatial_lifecycle.pdf`"
)
writeLines(md, file.path(out_dir, "ACS_SCF_HOME_VALUE_VALIDATION.md"))

cat("Wrote ACS/SCF home-value validation packet to:\n")
cat(out_dir, "\n")
cat("\nKey national ACS-vs-SCF comparison:\n")
print(key)
