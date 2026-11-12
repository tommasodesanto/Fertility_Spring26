#!/usr/bin/env Rscript
# Pools one headship design across states and runs the UNCHANGED national
# first-birth housing estimator (snapshot copy) on it. Writes only under OUTDIR.
options(stringsAsFactors = FALSE, scipen = 999); suppressPackageStartupMessages(library(data.table))
stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
design <- Sys.getenv("DESIGN"); base <- Sys.getenv("HEADSHIP_OUT"); estimator <- Sys.getenv("ESTIMATOR_FILE")
outdir <- file.path(base, "estimates", paste0(design, Sys.getenv("OUT_SUFFIX", "")))
reuse <- Sys.getenv("REUSE_DIR", ""); if (!nzchar(reuse)) reuse <- NULL
stopifnot(nzchar(design), dir.exists(base), file.exists(estimator))
if (dir.exists(outdir) && length(list.files(outdir))) stop("design OUTDIR non-empty: ", outdir)
dir.create(outdir, recursive = TRUE)
files <- sort(Sys.glob(file.path(base, "states", "statefip_*", paste0("design_", design, ".rds"))))
stopifnot(length(files) == 51L)
panel <- rbindlist(lapply(files, readRDS), use.names = TRUE)
message(stamp(), " pooled ", design, " rows=", nrow(panel))
summ <- panel[, .(rows = .N, persons = uniqueN(hs_src), households = uniqueN(source_hh_cluster),
                  anchors = uniqueN(hs_anchor[hs_role == "mother" & t_es_lw == 0])), by = .(event_time = t_es_lw, hs_role)]
fwrite(summ, file.path(outdir, "pool_counts_by_event.csv"))
fwrite(data.table(design = design, rows = nrow(panel), persons = uniqueN(panel$hs_src),
                  households = uniqueN(panel$source_hh_cluster),
                  anchors = uniqueN(panel$hs_anchor[panel$hs_role == "mother" & panel$t_es_lw == 0])),
       file.path(outdir, "pool_counts.csv"))
panel <- as.data.frame(panel); invisible(gc())
env <- new.env(parent = globalenv()); sys.source(estimator, envir = env)
cp <- function(x) message(stamp(), " ", jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"))
fit <- env$estimate_national_first_birth_housing(panel, output_dir = outdir, checkpoint = cp, reuse_dir = reuse,
                                                 geography_label = paste("National ACS, design", design))
jsonlite::write_json(list(status = fit$status, design = design, fits = length(fit$fits),
                          estimator = estimator, estimator_sha256 = digest::digest(file = estimator, algo = "sha256"),
                          generated = stamp()), file.path(outdir, "estimation_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
message(stamp(), " DESIGN_ESTIMATION_COMPLETE ", design)
