#!/usr/bin/env Rscript
# Builds a source-partition directory for a TRUE restricted rematch of one state
# (design B: householder, spouse or unmarried partner, households only), used
# only to validate filter-and-reweight. Raw ACS rows covered by the position key
# that fail the rule are dropped; uncovered rows and the CPS file are unchanged.
# Other states' housing packets are symlinked (read-only) so the unchanged
# matcher's packet-reuse gate passes. Writes only under PARTITION_OUTDIR.
options(stringsAsFactors = FALSE, scipen = 999)
src <- Sys.getenv("SOURCE_PARTITIONS", "/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_source_stage_20260921")
dst <- Sys.getenv("PARTITION_OUTDIR"); keydir <- Sys.getenv("KEY_DIR")
st <- as.integer(Sys.getenv("STATEFIP", "50")); hh_gq_file <- Sys.getenv("HH_GQ_FILE", "/scratch/td2248/projects/kleven_acs_pilot_20260917/code_headship_20261112/hh_gq.txt")
hh_gq_txt <- if (file.exists(hh_gq_file)) trimws(readLines(hh_gq_file, warn = FALSE)[1]) else Sys.getenv("HH_GQ", "1,2")
hh_gq <- as.numeric(strsplit(hh_gq_txt, ",")[[1]])
stopifnot(length(hh_gq) >= 2L, !anyNA(hh_gq))
stopifnot(nzchar(dst), nzchar(keydir), !dir.exists(dst) || !length(list.files(dst)))
dir.create(file.path(dst, "partitions"), recursive = TRUE)
for (f in c("partition_receipt.json", "source_schema_receipt.json", "housing_packet_manifest.json"))
  file.copy(file.path(src, f), file.path(dst, f))
for (d in list.dirs(file.path(src, "partitions"), recursive = FALSE)) {
  nd <- file.path(dst, "partitions", basename(d)); dir.create(nd)
  file.symlink(file.path(d, "housing_narrow.rds"), file.path(nd, "housing_narrow.rds"))
}
sd <- file.path(dst, "partitions", sprintf("statefip_%02d", st))
file.symlink(file.path(src, "partitions", sprintf("statefip_%02d", st), "cps_raw.RData"), file.path(sd, "cps_raw.RData"))
e <- new.env(); load(file.path(src, "partitions", sprintf("statefip_%02d", st), "acs_raw.RData"), envir = e); data <- e$data
col <- function(n) names(data)[tolower(names(data)) == n][1]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
k <- paste(num(data[[col("sample")]]), num(data[[col("year")]]), num(data[[col("serial")]]), num(data[[col("pernum")]]), sep = ":")
key <- readRDS(file.path(keydir, sprintf("statefip_%02d.rds", st)))
m <- match(k, key$KEY); covered <- !is.na(m)
elig <- !covered | (key$GQ[m] %in% hh_gq & (key$RELATE[m] %in% c(1, 2) | key$RELATED[m] %in% 1114))
rec <- list(statefip = st, raw_rows = nrow(data), covered = sum(covered), dropped = sum(!elig), hh_gq = hh_gq)
data <- data[elig, , drop = FALSE]
save(data, file = file.path(sd, "acs_raw.RData"), compress = FALSE)
jsonlite::write_json(rec, file.path(dst, "restricted_raw_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
print(rec)
