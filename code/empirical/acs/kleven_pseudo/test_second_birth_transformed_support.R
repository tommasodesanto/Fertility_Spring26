#!/usr/bin/env Rscript
# Tiny end-to-end source-join and transformed-cell diagnostic test.
suppressPackageStartupMessages(library(data.table))

here <- normalizePath(".", mustWork = TRUE)
driver <- file.path(here, "run_second_birth_transformed_support.R")
tmp <- tempfile("second_birth_support_"); dir.create(tmp)

packet <- data.table(
  YEAR = c(rep(2010L, 3), rep(2005L, 2), 2020L),
  SAMPLE = c(rep(201001L, 3), rep(200501L, 2), 202001L),
  SERIAL = c(rep(1L, 3), rep(2L, 2), 9L),
  PERNUM = c(1L, 2L, 3L, 1L, 2L, 1L),
  MOMLOC = c(0L, 1L, 1L, 0L, 1L, 0L),
  AGE = c(30L, 5L, 0L, 25L, 0L, 30L), SEX = c(2L, 1L, 1L, 2L, 1L, 2L),
  NCHILD = c(2L, 0L, 0L, 1L, 0L, 0L), FERTYR = c(2L, NA, NA, NA, NA, NA),
  PERWT = c(100, 50, 50, 40, 50, 100), MARST = 1L, RACE = 1L,
  EDUC = 73L, STATEFIP = 50L
)
prepared <- data.frame(
  year = packet$YEAR[1:5], sample = packet$SAMPLE[1:5], serial = packet$SERIAL[1:5],
  pernum = packet$PERNUM[1:5], sex = packet$SEX[1:5], age = packet$AGE[1:5],
  race = packet$RACE[1:5], educ = packet$EDUC[1:5], educd = 62L,
  marst = packet$MARST[1:5], statefip = packet$STATEFIP[1:5], hispan = 0L
)
packet_file <- file.path(tmp, "packet.rds"); saveRDS(packet, packet_file)
data <- prepared; prepared_file <- file.path(tmp, "prepared.RData"); save(data, file = prepared_file)
outdir <- file.path(tmp, "out")

old <- Sys.getenv(c("SECOND_BIRTH_PACKET", "SECOND_BIRTH_PREPARED", "SECOND_BIRTH_OUTDIR"),
                 unset = NA_character_)
on.exit({
  for (nm in names(old)) {
    if (is.na(old[[nm]])) Sys.unsetenv(nm)
    else do.call(Sys.setenv, setNames(list(old[[nm]]), nm))
  }
}, add = TRUE)
Sys.setenv(SECOND_BIRTH_PACKET = packet_file, SECOND_BIRTH_PREPARED = prepared_file,
           SECOND_BIRTH_OUTDIR = outdir)
status <- system2(Sys.which("Rscript"), driver, stdout = TRUE, stderr = TRUE)
exit_status <- attr(status, "status"); if (is.null(exit_status)) exit_status <- 0L
if (exit_status != 0L) stop(paste(status, collapse = "\n"), call. = FALSE)
manifest <- jsonlite::fromJSON(file.path(outdir, "diagnostic_manifest.json"))
stopifnot(identical(manifest$status, "TRANSFORMED_EXACT_MATCH_DIAGNOSTIC_ONLY"),
          manifest$outside_author_overlap_rows == 1L,
          file.exists(file.path(outdir, "proxy_checkpoint.rds")),
          file.exists(file.path(outdir, "matched_checkpoint.rds")))
conc <- read.csv(file.path(outdir, "raw_key_concordance.csv"), stringsAsFactors = FALSE)
stopifnot(all(conc$mismatches == 0L), all(conc$missing_either == 0L))
support <- read.csv(file.path(outdir, "target_support.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(support) == 5L, sum(support$candidate_donors == 1L) == 1L,
          support$candidate_donors[support$target_event_time == -5L] == 1L)
cat("PASS: transformed-cell source join and tiny end-to-end support diagnostic\n")
