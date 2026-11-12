#!/usr/bin/env Rscript
# Household-position key for the ACS first-birth headship designs (2026-11-12).
# Reads only key and position fields from extract27 and writes one key file per
# state plus tabulations of GQ, RELATE and RELATED by year. Reads no matched
# panel and writes nothing outside OUTDIR.
options(stringsAsFactors = FALSE, scipen = 999)
stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
root <- Sys.getenv("PROJECT_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
raw <- Sys.getenv("HOUSING_RAW", file.path(root, "inputs/ACS/local_extract27_20260919/extract27.dta"))
outdir <- Sys.getenv("OUTDIR")
stopifnot(nzchar(outdir), file.exists(raw))
if (dir.exists(outdir) && length(list.files(outdir))) stop("OUTDIR non-empty; refusing overwrite")
dir.create(file.path(outdir, "key"), recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({library(haven); library(data.table)})
message(stamp(), " read start ", raw)
h <- as.data.table(haven::read_dta(raw, col_select = c("year", "sample", "serial", "pernum", "statefip",
                                                         "gq", "relate", "related", "momloc", "poploc", "sploc")))
setnames(h, toupper(names(h)))
for (v in names(h)) h[[v]] <- as.numeric(zap_labels(h[[v]]))
message(stamp(), " read done rows=", nrow(h))
h[, KEY := paste(SAMPLE, YEAR, SERIAL, PERNUM, sep = ":")]
stopifnot(!anyDuplicated(h$KEY))
tab_gq <- h[, .N, by = .(YEAR, GQ)][order(YEAR, GQ)]
tab_rel <- h[, .N, by = .(YEAR, RELATE)][order(YEAR, RELATE)]
tab_related <- h[RELATE %in% c(1, 2, 11), .N, by = .(YEAR, RELATE, RELATED)][order(YEAR, RELATE, RELATED)]
fwrite(tab_gq, file.path(outdir, "tab_gq_by_year.csv"))
fwrite(tab_rel, file.path(outdir, "tab_relate_by_year.csv"))
fwrite(tab_related, file.path(outdir, "tab_related_heads_partners_by_year.csv"))
for (st in sort(unique(h$STATEFIP))) {
  saveRDS(h[STATEFIP == st, .(KEY, YEAR, GQ, RELATE, RELATED, MOMLOC, POPLOC, SPLOC)],
          file.path(outdir, "key", sprintf("statefip_%02d.rds", st)), compress = FALSE)
}
jsonlite::write_json(list(status = "HEADSHIP_KEY_COMPLETE", source = raw,
                          source_bytes = as.numeric(file.info(raw)$size), rows = nrow(h),
                          states = sort(unique(h$STATEFIP)), key_unique = TRUE, generated = stamp()),
                     file.path(outdir, "key_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
message(stamp(), " HEADSHIP_KEY_COMPLETE")
