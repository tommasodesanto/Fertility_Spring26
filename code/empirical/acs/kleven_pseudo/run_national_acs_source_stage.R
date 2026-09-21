#!/usr/bin/env Rscript
# National ACS source stage for the first-birth housing replication.
#
# This stage inventories the author's national ACS/CPS objects and writes
# state-partitioned checkpoints. It deliberately stops before the expensive
# national matcher unless the inventory and key gates pass. The vendor cleaner
# and matcher are staged separately and are never edited by this script.

options(stringsAsFactors = FALSE, scipen = 999)
started <- Sys.time()
stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
fail <- function(msg, stage = phase) {
  rec <- list(status = "FAILED", stage = stage, error = msg,
              started = stamp0, ended = stamp(),
              elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
  if (dir.exists(outdir) && requireNamespace("jsonlite", quietly = TRUE))
    jsonlite::write_json(rec, file.path(outdir, "failure_receipt.json"),
                         auto_unbox = TRUE, pretty = TRUE)
  message("FAIL: ", msg)
  stop(msg, call. = FALSE)
}
req <- function(ok, msg) if (!isTRUE(ok)) fail(msg)
stamp0 <- stamp()
phase <- Sys.getenv("PHASE", "inventory")
root <- Sys.getenv("PROJECT_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
outdir <- Sys.getenv("OUTDIR", file.path(root, "output", "national_acs_source_stage"))
valid_state_fips <- c(1L, 2L, 4L, 5L, 6L, 8L, 9L, 10L, 11L, 12L, 13L,
                      15L, 16L, 17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L,
                      25L, 26L, 27L, 28L, 29L, 30L, 31L, 32L, 33L, 34L,
                      35L, 36L, 37L, 38L, 39L, 40L, 41L, 42L, 44L, 45L,
                      46L, 47L, 48L, 49L, 50L, 51L, 53L, 54L, 55L, 56L)
state_text <- Sys.getenv("STATEFIP", paste(valid_state_fips, collapse = ","))
states <- suppressWarnings(as.integer(strsplit(state_text, ",", fixed = TRUE)[[1]]))
req(length(states) > 0L && all(is.finite(states) & states %in% valid_state_fips),
    "STATEFIP must be an explicit comma-separated list of valid state FIPS including DC=11")
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)))
  fail("OUTDIR exists and is non-empty; refusing overwrite", "startup")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("jsonlite", "digest")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
req(!length(missing), paste("missing required R packages:", paste(missing, collapse = ", ")))

acs_file <- Sys.getenv("ACS_RAW", file.path(root, "inputs", "ACS", "raw_acs.RData"))
cps_file <- Sys.getenv("CPS_RAW", file.path(root, "inputs", "CPS", "raw_cps.RData"))
vendor_dir <- Sys.getenv("VENDOR_DIR", file.path(root, "vendor"))
housing_packet <- Sys.getenv("HOUSING_PACKET", "")
housing_raw <- Sys.getenv("HOUSING_RAW", file.path(root, "inputs", "ACS", "local_extract27_20260919", "extract27.dta"))
req(file.exists(acs_file), paste("ACS raw input absent:", acs_file))
req(file.exists(cps_file), paste("CPS raw input absent:", cps_file))
req(dir.exists(vendor_dir), paste("vendor directory absent:", vendor_dir))

contract_file <- Sys.getenv("SOURCE_CONTRACT", file.path(root, "code", "empirical", "acs", "kleven_pseudo", "source_contract.json"))
req(file.exists(contract_file), paste("canonical source contract absent:", contract_file), "startup")
contract <- jsonlite::fromJSON(contract_file, simplifyVector = FALSE)
vendor_expected <- unlist(contract$vendor_sha256, use.names = TRUE)
req(identical(sort(names(vendor_expected)), sort(c("clean_acs.R", "clean_cps.R", "matching.R", "setup.R", "functions.R"))),
    "canonical source contract vendor set is incomplete", "startup")
vendor_hash <- vapply(names(vendor_expected), function(nm) {
  p <- file.path(vendor_dir, nm)
  req(file.exists(p), paste("vendor file absent:", p))
  digest::digest(file = p, algo = "sha256")
}, character(1))
req(all(vendor_hash == vendor_expected),
    paste("vendor hash mismatch:", paste(names(vendor_hash)[vendor_hash != vendor_expected], collapse = ", ")))

resolve_one <- function(nms, want, label, required = TRUE) {
  hit <- nms[tolower(nms) == tolower(want)]
  if (!length(hit) && !required) return(NA_character_)
  req(length(hit) == 1L, sprintf("%s field '%s' resolves to %d columns", label, want, length(hit)))
  hit[[1L]]
}
load_object <- function(path, label) {
  env <- new.env(parent = emptyenv())
  objs <- load(path, envir = env)
  req("data" %in% objs, paste(label, "must contain object named data"))
  d <- env$data
  req(is.data.frame(d), paste(label, "data object is not a data.frame"))
  d
}
source_inventory <- function(d, label, need_housing = FALSE) {
  nms <- names(d)
  cols <- list(
    year = resolve_one(nms, "year", label), sample = resolve_one(nms, "sample", label, FALSE),
    serial = resolve_one(nms, "serial", label), pernum = resolve_one(nms, "pernum", label),
    statefip = resolve_one(nms, "statefip", label), sex = resolve_one(nms, "sex", label),
    age = resolve_one(nms, "age", label), race = resolve_one(nms, "race", label, FALSE),
    hispan = resolve_one(nms, "hispan", label, FALSE), birthqtr = resolve_one(nms, "birthqtr", label, FALSE),
    educd = resolve_one(nms, "educd", label, FALSE), marst = resolve_one(nms, "marst", label, FALSE),
    perwt = resolve_one(nms, "perwt", label, FALSE),
    rooms = resolve_one(nms, "rooms", label, FALSE), bedrooms = resolve_one(nms, "bedrooms", label, FALSE),
    ownershp = resolve_one(nms, "ownershp", label, FALSE)
  )
  if (need_housing)
    req(all(!is.na(unlist(cols[c("rooms", "bedrooms", "ownershp")], use.names = FALSE)),
        paste(label, "missing a required housing field")))
  state <- suppressWarnings(as.integer(d[[cols$statefip]]))
  yr <- suppressWarnings(as.integer(d[[cols$year]]))
  list(label = label, rows = nrow(d), columns = ncol(d), fields = nms,
       resolved = cols, years = range(yr, na.rm = TRUE),
       states = sort(unique(state[is.finite(state)])),
       rows_by_state = as.list(table(state[is.finite(state)])),
       hispan_status = ifelse(is.na(cols$hispan), "ABSENT", "PRESENT"),
       birthqtr_status = ifelse(is.na(cols$birthqtr), "ABSENT", "PRESENT"))
}

message("PHASE ", phase, " start ", stamp())
partition_dir <- file.path(outdir, "partitions")
if (phase == "partition_smoke") dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)
progress_file <- file.path(outdir, "progress.log")
checkpoint <- function(...) {
  line <- paste(stamp(), paste(..., collapse = " "))
  cat(line, "\n", file = progress_file, append = TRUE)
  message(line)
}
partition_one <- function(path, label) {
  d <- load_object(path, label)
  meta <- source_inventory(d, label, need_housing = FALSE)
  partition_rows <- setNames(integer(length(states)), as.character(states))
  if (phase == "partition_smoke") {
    s <- suppressWarnings(as.integer(d[[meta$resolved$statefip]]))
    for (st in states) {
      sd <- d[!is.na(s) & s == st, , drop = FALSE]
      req(nrow(sd) > 0L, sprintf("%s state %s partition is empty", label, st))
      partition_rows[as.character(st)] <- nrow(sd)
      st_dir <- file.path(partition_dir, sprintf("statefip_%02d", st))
      dir.create(st_dir, recursive = TRUE, showWarnings = FALSE)
      data <- sd
      save(data, file = file.path(st_dir, paste0(tolower(label), "_raw.RData")), compress = FALSE)
      rm(data, sd)
      checkpoint("state_complete", label, st, partition_rows[as.character(st)])
    }
  }
  rm(d); invisible(gc(verbose = FALSE))
  meta$partition_rows <- partition_rows
  meta
}
acs_meta <- partition_one(acs_file, "ACS")
cps_meta <- partition_one(cps_file, "CPS")
req(!is.na(acs_meta$resolved$sample), "ACS raw source lacks SAMPLE; national bridge cannot be keyed safely")
jsonlite::write_json(list(status = "SCHEMA_PASS", phase = phase,
                          author_inputs = list(ACS = acs_meta, CPS = cps_meta),
                          vendor_sha256 = as.list(vendor_hash),
                          hispan_used_by_author_cleaner = !is.na(acs_meta$resolved$hispan),
                          birthqtr_status = acs_meta$birthqtr_status,
                          states_requested = states, generated = stamp()),
                     file.path(outdir, "source_schema_receipt.json"),
                     auto_unbox = TRUE, pretty = TRUE)

if (phase == "inventory") {
  jsonlite::write_json(list(status = "INVENTORY_COMPLETE", next_phase = "partition_smoke",
                            source_schema = file.path(outdir, "source_schema_receipt.json"),
                            no_matching_run = TRUE, generated = stamp()),
                       file.path(outdir, "stage_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
  message("INVENTORY_COMPLETE ", outdir)
  quit(save = "no", status = 0L)
}

req(identical(phase, "partition_smoke"), "PHASE must be inventory or partition_smoke")

# Housing is attached only after matching, but the smoke also materializes the
# narrow source fields needed by the later bridge when a packet or raw extract
# is available. It never carries the full extract into the matched panel.
housing_partition_receipt <- list(status = "NOT_RUN")
if (nzchar(housing_packet)) {
  req(file.exists(housing_packet), paste("HOUSING_PACKET absent:", housing_packet))
  req(requireNamespace("data.table", quietly = TRUE), "data.table required for HOUSING_PACKET")
  hp <- data.table::as.data.table(readRDS(housing_packet))
  need <- c("YEAR", "SAMPLE", "SERIAL", "PERNUM", "STATEFIP", "ROOMS_RAW", "BEDROOMS_RAW", "OWNERSHP_RAW")
  req(all(need %in% names(hp)), paste("housing packet missing:", paste(setdiff(need, names(hp)), collapse = ", ")))
  for (st in states) {
    hs <- hp[STATEFIP == st, ..need]
    req(nrow(hs) > 0L, sprintf("housing packet state %s partition is empty", st))
    st_dir <- file.path(partition_dir, sprintf("statefip_%02d", st))
    saveRDS(hs, file.path(st_dir, "housing_narrow.rds"), compress = FALSE)
  }
  housing_partition_receipt <- list(status = "PACKET_PARTITION_COMPLETE", rows = nrow(hp),
                                    states = states, fields = need)
  rm(hp); invisible(gc(verbose = FALSE))
} else {
  housing_partition_receipt <- list(status = "PENDING_RAW_EXTRACT_STAGE",
                                    raw_path = housing_raw,
                                    raw_present = file.exists(housing_raw),
                                    required_next = "read only YEAR,SAMPLE,SERIAL,PERNUM,STATEFIP,ROOMS,BEDROOMS,OWNERSHP and partition once")
}
partition_meta <- list(status = "PARTITION_COMPLETE", states = states,
                       ACS_rows_by_state = acs_meta$partition_rows,
                       CPS_rows_by_state = cps_meta$partition_rows,
                       ACS_hispan = acs_meta$hispan_status,
                       ACS_birthqtr = acs_meta$birthqtr_status,
                       housing = housing_partition_receipt,
                       next_gate = "run unchanged vendor cleaner/matcher on partition; then exact four-part housing bridge",
                       generated = stamp())
jsonlite::write_json(partition_meta, file.path(outdir, "partition_receipt.json"),
                     auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(list(status = "PARTITION_SMOKE_READY", no_full_matching = TRUE,
                          output = outdir, generated = stamp()),
                     file.path(outdir, "stage_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
message("PARTITION_SMOKE_READY ", outdir)
