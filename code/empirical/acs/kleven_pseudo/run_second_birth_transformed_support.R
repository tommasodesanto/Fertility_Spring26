#!/usr/bin/env Rscript
# Source-join diagnostic for author-transformed second-birth matching cells.
# This script has no housing regression and is intended for the approved
# Torch launcher. It joins only the author-overlap years 2005--2019.

suppressPackageStartupMessages(library(data.table))
source("second_birth_matching.R", local = TRUE)
source("build_second_birth_proxy.R", local = TRUE)

root <- Sys.getenv("KLEVEN_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
packet_file <- Sys.getenv("SECOND_BIRTH_PACKET", file.path(root, "output/kleven_acs_pilot",
  "source_audit_extract27_20260919/ne_extract27_housing_key_packet.rds"))
prepared_file <- Sys.getenv("SECOND_BIRTH_PREPARED", file.path(root, "overnight_ne_benchmark",
  "ne_acs.RData"))
outdir <- Sys.getenv("SECOND_BIRTH_OUTDIR", file.path(root, "output/kleven_acs_pilot",
  "second_birth_transformed_support_20260920"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fail <- function(...) stop(paste0(...), call. = FALSE)
resolve_one <- function(nms, candidates, label) {
  hit <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(hit) != 1L) fail(label, " resolves to ", length(hit), " columns")
  hit[[1L]]
}
key3 <- function(d, y, s, p) paste(d[[y]], d[[s]], d[[p]], sep = "\034")
num_equal <- function(a, b) {
  aa <- suppressWarnings(as.numeric(a)); bb <- suppressWarnings(as.numeric(b))
  !is.na(aa) & !is.na(bb) & aa != bb
}

if (!file.exists(packet_file) || !file.exists(prepared_file))
  fail("packet or prepared ACS file is missing")
packet <- as.data.table(readRDS(packet_file))
setnames(packet, toupper(names(packet)))
required_packet <- c("YEAR", "SAMPLE", "SERIAL", "PERNUM", "SEX", "AGE", "MARST",
                     "RACE", "EDUC", "STATEFIP", "MOMLOC", "NCHILD", "FERTYR", "PERWT")
if (length(setdiff(required_packet, names(packet))))
  fail("packet missing: ", paste(setdiff(required_packet, names(packet)), collapse = ", "))

env <- new.env(parent = emptyenv())
load(prepared_file, envir = env)
if (!exists("data", envir = env, inherits = FALSE)) fail("prepared RData lacks object data")
prepared <- as.data.table(env$data)
setnames(prepared, tolower(names(prepared)))
overlap_years <- 2005:2019

py <- resolve_one(names(prepared), "year", "prepared year")
ps <- resolve_one(names(prepared), "serial", "prepared serial")
pp <- resolve_one(names(prepared), "pernum", "prepared pernum")
raw_shared <- c(sex = "sex", age = "age", race = "race", educ = "educ",
                marst = "marst", statefip = "statefip")
raw_cols <- vapply(raw_shared, function(x) resolve_one(names(prepared), x, paste("prepared", x)), character(1))
hisp_col <- resolve_one(names(prepared), c("hispanic", "hispan"), "prepared Hispanic")
clean_author_names <- c("gender", "edlevel", "marst", "race", "statefip")
use_clean_author_cells <- all(clean_author_names %in% names(prepared))
if (use_clean_author_cells) {
  author_cols <- setNames(clean_author_names, clean_author_names)
} else {
  author_cols <- NULL
  educd_col <- resolve_one(names(prepared), "educd", "prepared detailed education")
}

packet_overlap <- packet[YEAR %in% overlap_years]
prep_overlap <- prepared[get(py) %in% overlap_years]
pk <- key3(packet_overlap, "YEAR", "SERIAL", "PERNUM")
rk <- key3(prep_overlap, py, ps, pp)
if (anyDuplicated(pk)) fail("packet overlap key (YEAR,SERIAL,PERNUM) is not unique")
if (anyDuplicated(rk)) fail("prepared overlap key (year,serial,pernum) is not unique")
idx <- match(pk, rk)
if (anyNA(idx)) fail("prepared ACS does not cover every packet overlap key")

concordance <- rbindlist(lapply(names(raw_shared), function(nm) {
  a <- packet_overlap[[raw_shared[[nm]]]]
  b <- prep_overlap[[raw_cols[[nm]]]][idx]
  data.table(field = nm, rows = length(a), missing_either = sum(xor(is.na(a), is.na(b))),
             mismatches = sum(num_equal(a, b)))
}))
write.csv(concordance, file.path(outdir, "raw_key_concordance.csv"), row.names = FALSE)
if (any(concordance$mismatches > 0L)) fail("shared raw fields disagree after key join")

# Borrow only the prepared Hispanic field after the raw concordance gate. Keep
# every packet row and every raw field; years outside 2005--2019 stay explicit.
packet[, author_overlap_2005_2019 := FALSE]
packet[match(pk, key3(packet, "YEAR", "SERIAL", "PERNUM")),
       author_overlap_2005_2019 := TRUE]
packet[, prepared_hispanic := NA_real_]
packet[, (paste0("author_", names(author_cols), ".num")) := NA]
prep_cells <- if (use_clean_author_cells) {
  second_birth_author_cells(prep_overlap,
    gender_col = author_cols[["gender"]], edlevel_col = author_cols[["edlevel"]],
    marst_col = author_cols[["marst"]], race_col = author_cols[["race"]],
    statefip_col = author_cols[["statefip"]])
} else {
  second_birth_author_cells_raw(prep_overlap, sex_col = raw_cols[["sex"]],
    educd_col = educd_col, marst_col = raw_cols[["marst"]],
    race_col = raw_cols[["race"]], hispan_col = hisp_col,
    statefip_col = raw_cols[["statefip"]])
}
for (j in seq_len(nrow(packet_overlap))) {
  ii <- match(pk[[j]], key3(packet, "YEAR", "SERIAL", "PERNUM"))
  jj <- idx[[j]]
  packet[ii, prepared_hispanic := if (tolower(hisp_col) == "hispan")
    as.numeric(prep_overlap[[hisp_col]][jj] %in% 1:4) else as.numeric(prep_overlap[[hisp_col]][jj])]
  packet[ii, `:=`(author_gender.num = prep_cells$gender.num[jj],
                  author_edlevel.num = prep_cells$edlevel.num[jj],
                  author_marst.num = prep_cells$marst.num[jj],
                  author_race.num = prep_cells$race.num[jj],
                  author_statefip.num = prep_cells$statefip.num[jj])]
}

# Build and match only on the verified overlap. 2020--2023 packet rows remain
# in the saved join output but cannot enter this author-aligned diagnostic.
dt <- packet[author_overlap_2005_2019 == TRUE]
data.table::setnames(dt, c("author_gender.num", "author_edlevel.num",
  "author_marst.num", "author_race.num", "author_statefip.num"),
  c("gender.num", "edlevel.num", "marst.num", "race.num", "statefip.num"))
proxy <- build_second_birth_proxy(dt,
  fertyr_codes = list(yes = 2L, no = 1L, unknown = c(0L, 8L)),
  match_covariates = c("gender.num", "edlevel.num", "marst.num", "race.num", "statefip.num"),
  event_window = -5:10, age_at_event_bounds = c(25, 45),
  full_pre_min_gap = 5L, reference_min_gap = 2L)
matched <- second_birth_match_exact(proxy$donor_targets, proxy$one_child_donors,
  match_covariates = c("gender.num", "edlevel.num", "marst.num", "race.num", "statefip.num"))

write.csv(proxy$audit, file.path(outdir, "builder_audit.csv"), row.names = FALSE)
write.csv(matched$target_support, file.path(outdir, "target_support.csv"), row.names = FALSE)
write.csv(matched$anchor_support, file.path(outdir, "anchor_support.csv"), row.names = FALSE)
write.csv(matched$donor_weights, file.path(outdir, "donor_weights_by_event_cell.csv"), row.names = FALSE)
saveRDS(packet, file.path(outdir, "packet_with_author_cells.rds"), compress = FALSE)
jsonlite::write_json(list(
  status = "TRANSFORMED_EXACT_MATCH_DIAGNOSTIC_ONLY",
  source_packet = packet_file, prepared_source = prepared_file,
  author_overlap_years = overlap_years,
  outside_author_overlap_rows = sum(!packet$author_overlap_2005_2019),
  raw_concordance = concordance,
  match_covariates = c("gender.num", "edlevel.num", "marst.num", "race.num", "statefip.num"),
  counts = matched$audit,
  housing_estimation_status = "not_run"),
  file.path(outdir, "diagnostic_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
cat("TRANSFORMED_SUPPORT_DIAGNOSTIC_COMPLETE\n")
