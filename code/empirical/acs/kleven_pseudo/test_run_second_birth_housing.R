#!/usr/bin/env Rscript
# Tiny end-to-end entry test for the second-birth housing driver.  It writes
# temporary proxy, matching, and source-packet checkpoints, then exercises all
# three declared specifications and verifies that complete fit receipts exist.
suppressPackageStartupMessages(library(data.table))

here <- normalizePath(".", mustWork = TRUE)
driver <- file.path(here, "run_second_birth_housing.R")
tmp <- tempfile("second_birth_housing_driver_"); dir.create(tmp)
outdir <- file.path(tmp, "out")

source_rows <- list(); anchors <- list(); post_rows <- list(); links <- list(); target_support <- list()
states <- c(23L, 25L, 33L, 44L, 50L, 9L)
n_anchor <- 24L
add_source <- function(year, sample, serial, pernum, age, state, weight, rooms, bedrooms, own) {
  data.table(YEAR = year, SAMPLE = sample, SERIAL = serial, PERNUM = pernum,
             person_key = paste(year, sample, serial, pernum, sep = "/"),
             AGE_norm = age, STATEFIP = state, PERWT = weight, SEX = 2L,
             ROOMS_RAW = rooms, BEDROOMS_RAW = bedrooms, OWNERSHP_RAW = own)
}
for (a in seq_len(n_anchor)) {
  state <- states[(a - 1L) %% length(states) + 1L]
  age <- 25 + (a %% 5)
  gap <- 2L + (a %% 3)
  fert <- if (a %% 5L == 0L) "no" else "yes"
  post_keys <- character()
  for (e in 0:3) {
    year <- 2005L + ((a + e) %% 8L)
    sample <- 200501L
    serial <- 10L * a + e + 1L
    pk <- paste(year, sample, serial, 1L, sep = "/")
    post_keys <- c(post_keys, pk)
    source_rows[[length(source_rows) + 1L]] <- add_source(year, sample, serial, 1L,
      age + (e %% 2L), state, 10 + a, 3 + (a %% 4L) + e, 2 + (e %% 2L), if (a %% 2L) 1 else 2)
    post_rows[[length(post_rows) + 1L]] <- data.table(person_key = pk, event_time = e,
      birth_gap = gap, PERWT = 10 + a, FERTYR_status = if (e == 0L) fert else NA_character_)
  }
  anchors[[length(anchors) + 1L]] <- data.table(person_key = post_keys[1L], birth_gap = gap,
    FERTYR_status = fert)
  for (et in c(-2L, -1L)) {
    year <- 2005L + ((a + et + 8L) %% 8L)
    sample <- 200501L
    serial <- 1000L + 10L * a + abs(et)
    dpk <- paste(year, sample, serial, 1L, sep = "/")
    source_rows[[length(source_rows) + 1L]] <- add_source(year, sample, serial, 1L,
      age, state, 5 + a, 2 + (a %% 3L), 1 + (a %% 2L), 2L)
    anchor_key <- post_keys[1L]
    links[[length(links) + 1L]] <- data.table(anchor_person_key = anchor_key,
      anchor_household_key = sub("/1$", "", anchor_key), target_event_time = et,
      target_year = year, target_mother_age = age, target_child_age = age - 5,
      donor_person_key = dpk, donor_household_key = sub("/1$", "", dpk),
      donor_PERWT = 5 + a, wgt_match = 1)
    target_support[[length(target_support) + 1L]] <- data.table(
      anchor_person_key = anchor_key, target_event_time = et, has_donor = TRUE)
  }
}
source_rows <- rbindlist(source_rows)
packet <- copy(source_rows)
proxy <- list(input = copy(source_rows), anchors = rbindlist(anchors),
              post_rows = rbindlist(post_rows), config = list())
matched <- list(links = rbindlist(links), target_support = rbindlist(target_support))
proxy_file <- file.path(tmp, "proxy_checkpoint.rds")
matched_file <- file.path(tmp, "matched_checkpoint.rds")
packet_file <- file.path(tmp, "packet_with_author_cells.rds")
saveRDS(proxy, proxy_file, compress = FALSE)
saveRDS(matched, matched_file, compress = FALSE)
saveRDS(packet, packet_file, compress = FALSE)

old <- Sys.getenv(c("SECOND_BIRTH_PROXY_FILE", "SECOND_BIRTH_MATCHED_FILE",
                    "SECOND_BIRTH_PACKET_FILE", "SECOND_BIRTH_HOUSING_OUTDIR"), unset = NA_character_)
on.exit({
  for (nm in names(old)) {
    if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old[[nm]]), nm))
  }
}, add = TRUE)
Sys.setenv(SECOND_BIRTH_PROXY_FILE = proxy_file, SECOND_BIRTH_MATCHED_FILE = matched_file,
           SECOND_BIRTH_PACKET_FILE = packet_file, SECOND_BIRTH_HOUSING_OUTDIR = outdir)
status <- system2(Sys.which("Rscript"), driver, stdout = TRUE, stderr = TRUE)
exit_status <- attr(status, "status"); if (is.null(exit_status)) exit_status <- 0L
if (exit_status != 0L) stop(paste(status, collapse = "\n"), call. = FALSE)
receipt <- jsonlite::fromJSON(file.path(outdir, "housing_run_receipt.json"))
stopifnot(receipt$status == "COMPLETE")
for (spec in c("primary", "joint_negative", "fertyr_event0_yes")) {
  st <- fread(file.path(outdir, spec, "fit_status.csv"))
  stopifnot(nrow(st) == 3L, all(st$status == "FIT_COMPLETE"),
            file.exists(file.path(outdir, spec, "curves.csv")),
            file.exists(file.path(outdir, spec, "contrasts.csv")))
}
cat("PASS: second-birth housing driver dependency closure, three specifications, and nine tiny fits\n")
