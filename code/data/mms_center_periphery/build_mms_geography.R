#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
  library(sf)
  library(stringr)
  library(tidyr)
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

is_true_env <- function(x) {
  tolower(x) %in% c("1", "true", "yes", "y")
}

script_dir <- get_script_dir()
data_suffix <- Sys.getenv("MMS_DATA_SUFFIX", "")
data_dir <- file.path(script_dir, paste0("data", data_suffix))
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

core_pop_share <- as.numeric(Sys.getenv("MMS_CORE_POP_SHARE", "0.10"))
center_puma_cutoff <- as.numeric(Sys.getenv("MMS_CENTER_PUMA_CORE_SHARE", "0.50"))
periphery_puma_cutoff <- as.numeric(Sys.getenv("MMS_PERIPHERY_PUMA_CORE_SHARE", "0.10"))

if (!is.finite(core_pop_share) || core_pop_share <= 0 || core_pop_share >= 1) {
  stop("MMS_CORE_POP_SHARE must be in (0,1).")
}
if (!is.finite(center_puma_cutoff) || !is.finite(periphery_puma_cutoff) ||
    center_puma_cutoff <= periphery_puma_cutoff ||
    periphery_puma_cutoff < 0 || center_puma_cutoff > 1) {
  stop("Invalid PUMA core-share cutoffs.")
}

paths <- list(
  tracts_regions = "/Users/tommasodesanto/Desktop/dataverse_files/data_construct/crosswalks/final_xwalks/tracts00_regions.dta",
  tract_coords_2000 = "/Users/tommasodesanto/Desktop/dataverse_files/data_construct/other_data/census2000/tract_coords_landwt_2000.csv",
  county_cbsa = "/Users/tommasodesanto/Desktop/Projects/Datasets/Geodata&Crosswalks/county_cbsa_crosswalk.dta",
  puma_shapes = "/Users/tommasodesanto/Desktop/Projects/Datasets/spatial_data/puma_fertility_contiguous.shp"
)

missing_paths <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0) {
  stop("Missing required inputs: ", paste(missing_paths, collapse = ", "))
}

message("Loading 2000 tract/CBD universe...")
tracts <- read_dta(paths$tracts_regions) %>%
  transmute(
    ctracts2000 = as.character(ctracts2000),
    region_all = as.integer(region_all),
    cbdname = as.character(cbdname),
    cbdx = as.numeric(cbdx),
    cbdy = as.numeric(cbdy),
    pop2000 = as.numeric(pop2000),
    cbddis = as.numeric(cbddis)
  ) %>%
  filter(!is.na(region_all), !is.na(pop2000), pop2000 > 0, !is.na(cbddis))

dup_ct <- tracts %>%
  count(ctracts2000, name = "n_regions") %>%
  filter(n_regions > 1)

if (nrow(dup_ct) > 0) {
  message("Tracts assigned to multiple regions: ", nrow(dup_ct), ". Keeping nearest CBD assignment.")
  tracts <- tracts %>%
    arrange(ctracts2000, cbddis, region_all) %>%
    group_by(ctracts2000) %>%
    slice(1) %>%
    ungroup()
}

region_pop <- tracts %>%
  group_by(region_all, cbdname) %>%
  summarise(region_pop2000 = sum(pop2000), .groups = "drop")

top_n_regions_env <- Sys.getenv("MMS_TOP_N_REGIONS", "")
top_n_regions <- suppressWarnings(as.integer(top_n_regions_env))
region_pop_cutoff <- as.numeric(Sys.getenv("MMS_REGION_POP_CUTOFF", "1e6"))

if (top_n_regions_env != "" && (!is.finite(top_n_regions) || top_n_regions <= 0)) {
  stop("MMS_TOP_N_REGIONS must be a positive integer when set.")
}

large_regions <- if (top_n_regions_env != "") {
  region_pop %>%
    arrange(desc(region_pop2000), region_all) %>%
    slice_head(n = top_n_regions)
} else {
  region_pop %>%
    filter(region_pop2000 >= region_pop_cutoff)
}

tracts <- tracts %>%
  semi_join(large_regions, by = c("region_all", "cbdname"))

region_pop <- region_pop %>%
  semi_join(large_regions, by = c("region_all", "cbdname"))

message("Mapping region_all markets to CBSAs through counties...")
county_cbsa <- read_dta(paths$county_cbsa) %>%
  transmute(
    cbsacode = suppressWarnings(as.integer(cbsacode)),
    cbsatitle = as.character(cbsatitle),
    metro_type = as.character(metropolitanmicropolitanstatis),
    county_fips = paste0(
      str_pad(as.character(fipsstatecode), width = 2, side = "left", pad = "0"),
      str_pad(as.character(fipscountycode), width = 3, side = "left", pad = "0")
    )
  ) %>%
  filter(!is.na(cbsacode), str_detect(metro_type, "Metropolitan"))

region_to_cbsa <- tracts %>%
  mutate(county_fips = substr(ctracts2000, 1, 5)) %>%
  left_join(county_cbsa, by = "county_fips") %>%
  filter(!is.na(cbsacode)) %>%
  group_by(region_all, cbdname, cbsacode, cbsatitle) %>%
  summarise(overlap_pop2000 = sum(pop2000), .groups = "drop") %>%
  left_join(region_pop, by = c("region_all", "cbdname")) %>%
  mutate(overlap_share = overlap_pop2000 / region_pop2000) %>%
  arrange(region_all, desc(overlap_pop2000), cbsacode) %>%
  group_by(region_all, cbdname) %>%
  slice(1) %>%
  ungroup()

write_csv(region_to_cbsa, file.path(data_dir, "region_to_cbsa.csv"))

tracts <- tracts %>%
  left_join(region_to_cbsa %>% select(region_all, cbdname, cbsacode, cbsatitle, overlap_share),
            by = c("region_all", "cbdname"))

message("Flagging downtown and outer-suburb tracts...")
downtown_flags <- tracts %>%
  arrange(region_all, cbddis, ctracts2000) %>%
  group_by(region_all) %>%
  mutate(
    total_pop2000 = sum(pop2000),
    cum_pop = cumsum(pop2000),
    cum_pop_before = lag(cum_pop, default = 0),
    downtown_tract = cum_pop_before < core_pop_share * total_pop2000
  ) %>%
  ungroup() %>%
  select(ctracts2000, downtown_tract)

outer_flags <- tracts %>%
  arrange(region_all, desc(cbddis), ctracts2000) %>%
  group_by(region_all) %>%
  mutate(
    total_pop2000 = sum(pop2000),
    outer_cum_pop = cumsum(pop2000),
    outer_cum_before = lag(outer_cum_pop, default = 0),
    outer50_tract = outer_cum_before < 0.50 * total_pop2000
  ) %>%
  ungroup() %>%
  select(ctracts2000, outer50_tract)

tracts <- tracts %>%
  left_join(downtown_flags, by = "ctracts2000") %>%
  left_join(outer_flags, by = "ctracts2000")

message("Loading 2000 tract centroid coordinates...")
coords_2000 <- read_csv(paths$tract_coords_2000, show_col_types = FALSE) %>%
  transmute(
    ctracts2000 = paste0(
      str_pad(as.character(county), width = 5, side = "left", pad = "0"),
      str_pad(as.character(as.integer(round(tract * 100))), width = 6, side = "left", pad = "0")
    ),
    lon = as.numeric(intptlon),
    lat = as.numeric(intptlat),
    land_sqmi_2000 = as.numeric(LandSQMI)
  ) %>%
  distinct(ctracts2000, .keep_all = TRUE)

tracts_geo <- tracts %>%
  left_join(coords_2000, by = "ctracts2000")

coord_match_rate <- mean(!is.na(tracts_geo$lon) & !is.na(tracts_geo$lat))
message("2000 tract centroid match rate: ", sprintf("%.1f%%", 100 * coord_match_rate))

tracts_geo <- tracts_geo %>%
  filter(!is.na(lon), !is.na(lat))

sf_use_s2(FALSE)

message("Reading local PUMA geometry...")
pumas <- st_read(paths$puma_shapes, quiet = TRUE) %>%
  transmute(
    GEOID20 = as.integer(GEOID20),
    pum2010 = as.integer(pum2010),
    STATEFP = as.integer(STATEFP)
  ) %>%
  distinct(GEOID20, pum2010, .keep_all = TRUE) %>%
  st_make_valid()

tract_points <- st_as_sf(
  tracts_geo,
  coords = c("lon", "lat"),
  crs = 4269,
  remove = FALSE
) %>%
  st_transform(st_crs(pumas))

message("Assigning 2000 tract centroids to PUMAs...")
tract_puma <- st_join(
  tract_points,
  pumas,
  left = FALSE,
  join = st_within
) %>%
  st_drop_geometry()

if (nrow(tract_puma) == 0) {
  stop("No tract centroids were assigned to a PUMA.")
}

message("Building PUMA-level downtown shares...")
puma_lookup <- tract_puma %>%
  group_by(region_all, cbdname, cbsacode, cbsatitle, pum2010, GEOID20, STATEFP) %>%
  summarise(
    tract_pop2000 = sum(pop2000),
    downtown_pop2000 = sum(pop2000[downtown_tract], na.rm = TRUE),
    outer50_pop2000 = sum(pop2000[outer50_tract], na.rm = TRUE),
    n_tracts = n(),
    .groups = "drop"
  ) %>%
  mutate(
    downtown_share = downtown_pop2000 / tract_pop2000,
    outer50_share = outer50_pop2000 / tract_pop2000,
    mms_location = case_when(
      downtown_share >= center_puma_cutoff ~ "center",
      downtown_share < periphery_puma_cutoff ~ "periphery",
      TRUE ~ "middle"
    )
  ) %>%
  arrange(cbsacode, pum2010, desc(tract_pop2000), region_all) %>%
  group_by(pum2010, GEOID20, cbsacode) %>%
  slice(1) %>%
  ungroup()

tract_puma_class <- tract_puma %>%
  left_join(
    puma_lookup %>%
      select(region_all, pum2010, GEOID20, mms_location),
    by = c("region_all", "pum2010", "GEOID20")
  )

city_coverage <- tract_puma_class %>%
  group_by(region_all, cbdname, cbsacode, cbsatitle) %>%
  summarise(
    region_pop2000 = sum(pop2000),
    downtown_tract_pop2000 = sum(pop2000[downtown_tract], na.rm = TRUE),
    downtown_pop_in_center_pumas = sum(pop2000[downtown_tract & mms_location == "center"], na.rm = TRUE),
    center_puma_pop2000 = sum(pop2000[mms_location == "center"], na.rm = TRUE),
    n_center_pumas = n_distinct(pum2010[mms_location == "center"]),
    downtown_capture_share = downtown_pop_in_center_pumas / downtown_tract_pop2000,
    .groups = "drop"
  ) %>%
  mutate(keep_city = downtown_capture_share >= 0.50)

write_csv(city_coverage, file.path(data_dir, "city_coverage.csv"))

kept_cities <- city_coverage %>%
  filter(keep_city) %>%
  select(region_all, cbsacode)

puma_lookup <- puma_lookup %>%
  semi_join(kept_cities, by = c("region_all", "cbsacode"))

tracts_out <- tracts_geo %>%
  semi_join(kept_cities, by = c("region_all", "cbsacode")) %>%
  select(
    ctracts2000, region_all, cbdname, cbsacode, cbsatitle, overlap_share,
    pop2000, cbddis, downtown_tract, outer50_tract, lon, lat
  )

write_csv(region_pop, file.path(data_dir, "region_pop2000.csv"))
write_csv(tracts_out, file.path(data_dir, "tract2000_mms_flags.csv.gz"))

puma_lookup_2010 <- puma_lookup %>%
  transmute(
    year_start = 2012L,
    year_end = 2021L,
    region_all,
    cbdname,
    cbsacode,
    cbsatitle,
    statefip = as.integer(substr(sprintf("%07d", pum2010), 1, 2)),
    puma = as.integer(substr(sprintf("%07d", pum2010), 3, 7)),
    puma_id = pum2010,
    downtown_share,
    outer50_share,
    tract_pop2000,
    mms_location
  ) %>%
  arrange(cbsacode, statefip, puma)

puma_lookup_2020 <- puma_lookup %>%
  transmute(
    year_start = 2022L,
    year_end = 2023L,
    region_all,
    cbdname,
    cbsacode,
    cbsatitle,
    statefip = as.integer(substr(sprintf("%07d", GEOID20), 1, 2)),
    puma = as.integer(substr(sprintf("%07d", GEOID20), 3, 7)),
    puma_id = GEOID20,
    downtown_share,
    outer50_share,
    tract_pop2000,
    mms_location
  ) %>%
  arrange(cbsacode, statefip, puma)

write_csv(puma_lookup_2010, file.path(data_dir, "puma_mms_lookup_2010.csv"))
write_csv(puma_lookup_2020, file.path(data_dir, "puma_mms_lookup_2020.csv"))

summary_lines <- c(
  sprintf("Core tract population target share: %.2f", core_pop_share),
  sprintf("Center PUMA cutoff on core share: %.2f", center_puma_cutoff),
  sprintf("Periphery PUMA cutoff on core share: %.2f", periphery_puma_cutoff),
  sprintf("Region selection rule: %s",
          if (top_n_regions_env != "") {
            sprintf("top %d by 2000 region population", top_n_regions)
          } else {
            sprintf("region_pop2000 >= %.0f", region_pop_cutoff)
          }),
  sprintf("Large regions kept before PUMA coverage screen: %d", nrow(region_pop)),
  sprintf("Cities passing downtown-coverage screen: %d", nrow(filter(city_coverage, keep_city))),
  sprintf("2000 tract centroid coordinate match rate: %.1f%%", 100 * coord_match_rate),
  sprintf("PUMA-city pairs in 2010 lookup: %d", nrow(puma_lookup_2010)),
  sprintf("PUMA-city pairs in 2020 lookup: %d", nrow(puma_lookup_2020))
)
writeLines(summary_lines, file.path(data_dir, "build_summary.txt"))

message("MMS geography build complete.")
