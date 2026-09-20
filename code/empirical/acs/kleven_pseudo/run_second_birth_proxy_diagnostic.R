suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("KLEVEN_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
audit_dir <- file.path(root, "forensics/source_audit_extract27_20260919")
packet_file <- file.path(root, "output/kleven_acs_pilot/source_audit_extract27_20260919",
                         "ne_extract27_housing_key_packet.rds")
outdir <- file.path(root, "output/kleven_acs_pilot/second_birth_proxy_diagnostic_20260920")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
progress_file <- file.path(outdir, "progress.log")
checkpoint <- function(stage, detail = "") {
  line <- sprintf("%s\t%s\t%s", format(Sys.time(), "%FT%T%z"), stage, detail)
  write(line, file = progress_file, append = TRUE)
  cat(line, "\n")
}
options(error = function() {
  msg <- geterrmessage()
  writeLines(msg, file.path(outdir, "failure.txt"))
  q(save = "no", status = 1, runLast = FALSE)
})

checkpoint("startup")
stopifnot(file.exists(packet_file), file.exists(file.path(audit_dir, "build_second_birth_proxy.R")))
source(file.path(audit_dir, "build_second_birth_proxy.R"), local = TRUE)
checkpoint("builder_loaded")

dt <- as.data.table(readRDS(packet_file))
setnames(dt, toupper(names(dt)))
required <- c("YEAR", "SAMPLE", "SERIAL", "PERNUM", "MOMLOC", "AGE", "SEX",
              "NCHILD", "FERTYR", "PERWT", "EDUC", "MARST", "RACE", "STATEFIP")
stopifnot(!length(setdiff(required, names(dt))))
numeric_fields <- intersect(required, names(dt))
dt[, (numeric_fields) := lapply(.SD, as.numeric), .SDcols = numeric_fields]
checkpoint("packet_loaded", sprintf("rows=%s cols=%s", nrow(dt), ncol(dt)))

sample_label <- function(x) {
  out <- rep("unknown", length(x))
  out[x %% 100L == 1L & x != 200004L] <- "ACS 1-year"
  out[x %% 100L == 2L] <- "PRCS 1-year"
  out[x %% 100L == 3L] <- "ACS 5-year"
  out[x %% 100L == 4L] <- "PRCS 5-year"
  out[x == 200004L] <- "ACS 2000"
  out
}

proxy <- build_second_birth_proxy(
  dt,
  fertyr_codes = list(yes = 2L, no = 1L, unknown = c(0L, 8L)),
  match_covariates = c("SEX", "EDUC", "MARST", "RACE", "STATEFIP"),
  event_window = -5:10,
  age_at_event_bounds = c(25, 45),
  full_pre_min_gap = 5L,
  reference_min_gap = 2L
)
checkpoint("builder_complete", sprintf("anchors=%s donors=%s targets=%s",
                                        nrow(proxy$anchors), nrow(proxy$one_child_donors),
                                        nrow(proxy$donor_targets)))

write.csv(proxy$audit, file.path(outdir, "builder_audit.csv"), row.names = FALSE)

mother_year <- proxy$mother_rows[, .(
  female_valid_age_mothers = .N,
  strict_eligible_mothers = sum(strict_eligible, na.rm = TRUE),
  event0_anchors = sum(strict_eligible & event_time == 0, na.rm = TRUE),
  full_pre_anchors = sum(gap_full_pre & event_time == 0, na.rm = TRUE),
  reference_anchors = sum(gap_reference & event_time == 0, na.rm = TRUE),
  event0_fertyr_yes = sum(event_time == 0 & FERTYR_status == "yes", na.rm = TRUE),
  event0_fertyr_no = sum(event_time == 0 & FERTYR_status == "no", na.rm = TRUE),
  event0_fertyr_unknown = sum(event_time == 0 & FERTYR_status == "unknown", na.rm = TRUE),
  event0_fertyr_missing = sum(event_time == 0 & is.na(FERTYR_status), na.rm = TRUE)
), by = .(YEAR, SAMPLE)]
mother_year[, sample_label := sample_label(SAMPLE)]
setcolorder(mother_year, c("YEAR", "SAMPLE", "sample_label"))
write.csv(mother_year, file.path(outdir, "support_by_year_sample.csv"), row.names = FALSE)

fertyr_support <- proxy$mother_rows[event_time == 0,
  .(rows = .N, strict_eligible = sum(strict_eligible, na.rm = TRUE)),
  by = .(FERTYR_status)]
write.csv(fertyr_support, file.path(outdir, "fertyr_event0_support.csv"), row.names = FALSE)

gap_support <- proxy$anchors[, .(
  anchors = .N,
  full_pre = sum(gap_full_pre, na.rm = TRUE),
  reference = sum(gap_reference, na.rm = TRUE)
), by = .(birth_gap)][order(birth_gap)]
write.csv(gap_support, file.path(outdir, "anchor_gap_support.csv"), row.names = FALSE)

event_support <- rbindlist(list(
  proxy$post_rows[, .(rows = .N), by = .(event_time)][, source := "post_rows"],
  proxy$donor_targets[, .(rows = .N), by = .(event_time = target_event_time)][, source := "donor_targets"]
), fill = TRUE)
setorder(event_support, event_time, source)
write.csv(event_support, file.path(outdir, "event_time_support.csv"), row.names = FALSE)

link_support <- proxy$links[, .(rows = .N), by = .(link_valid, link_invalid_reason)]
write.csv(link_support, file.path(outdir, "link_support.csv"), row.names = FALSE)

donor_year <- proxy$one_child_donors[, .(one_child_donors = .N), by = .(YEAR, SAMPLE)]
donor_year[, sample_label := sample_label(SAMPLE)]
write.csv(donor_year, file.path(outdir, "one_child_donors_by_year_sample.csv"), row.names = FALSE)

saveRDS(proxy$anchors, file.path(outdir, "anchors.rds"), compress = FALSE)
saveRDS(proxy$post_rows, file.path(outdir, "post_rows.rds"), compress = FALSE)
saveRDS(proxy$one_child_donors, file.path(outdir, "one_child_donors.rds"), compress = FALSE)
saveRDS(proxy$donor_targets, file.path(outdir, "donor_targets.rds"), compress = FALSE)

manifest <- list(
  packet = packet_file,
  rows = nrow(dt),
  match_covariates = c("SEX", "EDUC", "MARST", "RACE", "STATEFIP"),
  fertyr_codes = list(yes = 2L, no = 1L, unknown = c(0L, 8L)),
  event_window = -5:10,
  age_at_event_bounds = c(25L, 45L),
  full_pre_min_gap = 5L,
  reference_min_gap = 2L,
  matching_status = "not_run",
  housing_estimation_status = "not_run",
  counts = list(
    female_valid_age_mothers = nrow(proxy$mother_rows),
    anchors = nrow(proxy$anchors),
    full_pre_anchors = sum(proxy$anchors$gap_full_pre, na.rm = TRUE),
    reference_anchors = sum(proxy$anchors$gap_reference, na.rm = TRUE),
    post_rows = nrow(proxy$post_rows),
    one_child_donors = nrow(proxy$one_child_donors),
    donor_targets = nrow(proxy$donor_targets),
    valid_links = sum(proxy$links$link_valid, na.rm = TRUE)
  ),
  generated = as.character(Sys.time())
)
jsonlite::write_json(manifest, file.path(outdir, "diagnostic_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)
checkpoint("complete")
