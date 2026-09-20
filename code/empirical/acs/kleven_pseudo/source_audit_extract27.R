suppressPackageStartupMessages({
  library(haven)
  library(data.table)
  library(jsonlite)
})

root <- Sys.getenv("KLEVEN_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
infile <- file.path(root, "inputs/ACS/local_extract27_20260919/extract27.dta")
outdir <- file.path(root, "output/kleven_acs_pilot/source_audit_extract27_20260919")
receipt_file <- file.path(root, "forensics/source_audit_extract27_20260919/source_digest_receipt.json")
prepared_file <- file.path(root, "overnight_ne_benchmark/ne_acs.RData")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

progress_file <- file.path(outdir, "progress.log")
checkpoint <- function(stage, detail = "") {
  line <- sprintf("%s\t%s\t%s", format(Sys.time(), "%FT%T%z"), stage, detail)
  write(line, file = progress_file, append = TRUE)
  cat(line, "\n")
}
failure_file <- file.path(outdir, "source_audit_failure_receipt.json")
options(error = function() {
  msg <- geterrmessage()
  try(write_json(list(status = "FAILED", error = msg, generated = as.character(Sys.time())),
                 failure_file, auto_unbox = TRUE, pretty = TRUE), silent = TRUE)
  q(save = "no", status = 1, runLast = FALSE)
})
checkpoint("startup")

expected_sha256 <- "edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e"
expected_bytes <- 9919999546
ne_states <- c(9L, 23L, 25L, 33L, 44L, 50L)
keycols <- c("year", "sample", "serial", "pernum")
housing_cols <- c("rooms", "bedrooms", "ownershp")

stopifnot(file.exists(infile), file.exists(receipt_file))
stopifnot(file.exists(prepared_file))
checkpoint("prepared_source_path_verified", prepared_file)
receipt <- fromJSON(receipt_file)
stopifnot(identical(receipt$sha256, expected_sha256),
          as.numeric(receipt$bytes) == expected_bytes,
          as.numeric(file.info(infile)$size) == expected_bytes)
checkpoint("source_receipt_verified", sprintf("bytes=%s sha256_len=%s", expected_bytes, nchar(expected_sha256)))

# These names follow extractor27.do. Selection is resolved against the actual
# header below because the Stata file uses lowercase names.
fields <- c(
  "YEAR", "SAMPLE", "SERIAL", "CBSERIAL", "HHWT", "CLUSTER", "STATEFIP",
  "PUMA", "STRATA", "GQ", "OWNERSHP", "OWNERSHPD", "ROOMS", "BEDROOMS",
  "PERNUM", "PERWT", "MOMLOC", "POPLOC", "NCHILD", "NCHLT5", "ELDCH",
  "YNGCH", "RELATE", "SEX", "AGE", "MARST", "FERTYR", "RACE", "EDUC"
)
checkpoint("header_smoke_start")
header <- read_dta(infile, n_max = 1, .name_repair = "minimal")
actual_header_names <- names(header)
header_lookup <- setNames(actual_header_names, toupper(actual_header_names))
selected_fields <- unname(header_lookup[toupper(fields)])
missing_header_fields <- fields[is.na(selected_fields)]
stopifnot(!length(missing_header_fields))
write.csv(data.table(requested_field = fields, actual_field = selected_fields, present = !is.na(selected_fields)),
          file.path(outdir, "header_selection_check.csv"), row.names = FALSE)
rm(header)
checkpoint("header_smoke_pass", sprintf("selected_fields=%s", length(fields)))
d <- as.data.table(read_dta(
  infile,
  col_select = tidyselect::all_of(selected_fields),
  .name_repair = "minimal"
))
setnames(d, tolower(names(d)))
required <- c(keycols, "statefip", "hhwt", "perwt", "gq", "relate", "momloc",
              "sex", "age", "nchild", "eldch", "yngch", "fertyr", housing_cols)
stopifnot(!length(setdiff(required, names(d))))
checkpoint("full_extract_read", sprintf("rows=%s cols=%s", nrow(d), ncol(d)))
numeric_cols <- intersect(names(d), c(keycols, "statefip", "hhwt", "perwt", "gq",
                                       "relate", "momloc", "sex", "age", "nchild",
                                       "eldch", "yngch", "fertyr", housing_cols))
d[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]

rooms_expected <- c(0:9, 10:27, 30)
bedrooms_expected <- 0:22
ownershp_expected <- 0:2
expected_by_variable <- list(
  rooms = rooms_expected,
  bedrooms = bedrooms_expected,
  ownershp = ownershp_expected
)

# Preserve literal source codes and make unknown outcome values missing for
# future outcome work. Rows, keys, weights, and the raw values remain intact.
for (v in housing_cols) {
  raw_name <- paste0(v, "_raw")
  valid_name <- paste0(v, "_code_valid")
  outcome_name <- paste0(v, "_outcome")
  d[, (raw_name) := get(v)]
  d[, (valid_name) := fifelse(is.na(get(raw_name)), NA, get(raw_name) %in% expected_by_variable[[v]])]
  d[, (outcome_name) := fifelse(!is.na(get(valid_name)) & !get(valid_name), NA_real_, get(raw_name))]
}

packet_cols <- intersect(c(keycols, "cbserial", "statefip", "puma", "gq", "ownershp",
                           "ownershpd", "rooms", "bedrooms", "hhwt", "perwt", "cluster",
                           "strata", "relate", "sex", "age", "marst", "fertyr", "race",
                           "educ", "momloc", "poploc", "nchild", "nchlt5", "eldch", "yngch",
                           unlist(lapply(housing_cols, function(v) c(
                             paste0(v, "_raw"), paste0(v, "_code_valid"), paste0(v, "_outcome"))))),
                         names(d))
d_ne <- d[statefip %in% ne_states]
housing_packet <- d_ne[, ..packet_cols]
saveRDS(housing_packet, file.path(outdir, "ne_extract27_housing_key_packet.rds"), compress = FALSE)
checkpoint("ne_packet_saved_immediate", sprintf("rows=%s cols=%s", nrow(housing_packet), ncol(housing_packet)))

sample_product <- function(x) {
  fifelse(x %% 100L == 1L, "ACS 1-year",
  fifelse(x %% 100L == 3L, "ACS 5-year",
  fifelse(x %% 100L == 2L, "PRCS 1-year",
  fifelse(x %% 100L == 4L, "PRCS 5-year", "unknown"))))
}

# Compact national fingerprint: no national key merge or national housing join.
national_ys <- d[, .N, by = .(year, sample)][order(year, sample)]
national_ys[, sample_product := sample_product(sample)]
write.csv(national_ys, file.path(outdir, "national_year_sample_counts.csv"), row.names = FALSE)
write.csv(national_ys[, .N, by = .(sample, sample_product)][order(sample)],
          file.path(outdir, "national_sample_product_counts.csv"), row.names = FALSE)

code_table <- rbindlist(list(
  data.table(variable = "rooms", code = rooms_expected),
  data.table(variable = "bedrooms", code = bedrooms_expected),
  data.table(variable = "ownershp", code = ownershp_expected)
))
write.csv(code_table, file.path(outdir, "authoritative_housing_codes.csv"), row.names = FALSE)
code_counts <- function(x, variable) {
  out <- data.table(code = x)[, .(N = .N), by = code]
  out[, variable := variable]
  out[, .(variable, code, N)]
}
fixture <- data.table(rooms = c(0, 1, 1, NA_real_),
                      bedrooms = c(0, 1, 1, NA_real_),
                      ownershp = c(0, 1, 1, NA_real_))
for (v in housing_cols) {
  tiny <- code_counts(fixture[[v]], v)
  stopifnot(nrow(tiny) == uniqueN(fixture[[v]]), sum(tiny$N) == nrow(fixture))
}
observed_codes <- rbindlist(lapply(housing_cols, function(v) code_counts(d[[v]], v)))
stopifnot(all(vapply(housing_cols, function(v) {
  z <- observed_codes[variable == v]
  nrow(z) == uniqueN(d[[v]]) && sum(z$N) == nrow(d)
}, logical(1))))
observed_codes[, code_authoritative := mapply(
  function(v, z) z %in% code_table[variable == v, code], variable, code
)]
write.csv(observed_codes, file.path(outdir, "housing_code_validation.csv"), row.names = FALSE)
unknown_counts <- rbindlist(lapply(housing_cols, function(v) {
  expected <- expected_by_variable[[v]]
  d[!is.na(get(v)) & !get(v) %in% expected,
    .(N = .N), by = .(year, sample, code = get(v))][, variable := v]
}), fill = TRUE)
if (nrow(unknown_counts)) {
  unknown_counts[, sample_product := sample_product(sample)]
  setcolorder(unknown_counts, c("variable", "year", "sample", "sample_product", "code", "N"))
} else {
  unknown_counts <- data.table(variable = character(), year = integer(), sample = integer(),
                               sample_product = character(), code = numeric(), N = integer())
}
write.csv(unknown_counts, file.path(outdir, "unknown_housing_codes_by_year_sample.csv"), row.names = FALSE)
unknown_ne <- rbindlist(lapply(housing_cols, function(v) {
  expected <- expected_by_variable[[v]]
  d_ne[!is.na(get(v)) & !get(v) %in% expected,
       .(N = .N), by = .(year, sample, code = get(v))][, variable := v]
}), fill = TRUE)
if (nrow(unknown_ne)) {
  unknown_ne[, sample_product := sample_product(sample)]
  setcolorder(unknown_ne, c("variable", "year", "sample", "sample_product", "code", "N"))
} else {
  unknown_ne <- data.table(variable = character(), year = integer(), sample = integer(),
                           sample_product = character(), code = numeric(), N = integer())
}
write.csv(unknown_ne, file.path(outdir, "unknown_housing_codes_ne_by_year_sample.csv"), row.names = FALSE)
write.csv(data.table(
  variable = "rooms",
  source = "IPUMS USA variable page: ROOMS",
  url = "https://usa.ipums.org/usa-action/variables/ROOMS",
  category_definition = "N/A, 1-27, and 30; no category 28 or 29"
), file.path(outdir, "housing_codebook_evidence.csv"), row.names = FALSE)
unknown_total <- sum(unknown_counts$N)
unknown_ne_total <- sum(unknown_ne$N)
checkpoint("national_fingerprint_and_codes",
           sprintf("year_sample_rows=%s unknown_codes=%s unknown_ne=%s", nrow(national_ys), unknown_total, unknown_ne_total))

# Restrict all joins and concordance to the approved Northeast source universe.
rm(housing_packet)
rm(d)
gc()
checkpoint("ne_filter", sprintf("rows=%s states=%s", nrow(d_ne), paste(ne_states, collapse = ",")))

author_env <- new.env(parent = emptyenv())
author_objects <- load(prepared_file, envir = author_env)
data_objects <- author_objects[vapply(author_objects, function(x) {
  inherits(author_env[[x]], "data.frame")
}, logical(1))]
stopifnot(length(data_objects) == 1L)
author <- as.data.table(author_env[[data_objects[[1L]]]])
setnames(author, tolower(names(author)))
stopifnot(!length(setdiff(c(keycols, "statefip", "sex", "age", "ownershp"), names(author))))
author_numeric <- intersect(names(author), c(keycols, "statefip", "sex", "age", "ownershp"))
author[, (author_numeric) := lapply(.SD, as.numeric), .SDcols = author_numeric]
author_ne <- author[statefip %in% ne_states]
rm(author, author_env)
gc()
checkpoint("prepared_ne_loaded", sprintf("rows=%s", nrow(author_ne)))

local_ys <- d_ne[, .N, by = .(year, sample)][order(year, sample)]
local_ys[, sample_product := sample_product(sample)]
author_ys <- author_ne[, .N, by = .(year, sample)][order(year, sample)]
author_ys[, sample_product := sample_product(sample)]
write.csv(merge(local_ys, author_ys, by = c("year", "sample", "sample_product"),
                all = TRUE, suffixes = c("_extract27", "_prepared_ne")),
          file.path(outdir, "ne_year_sample_comparison.csv"), row.names = FALSE)

local_counts <- d_ne[, .N, by = keycols]
author_counts <- author_ne[, .N, by = keycols]
local_unique_keys <- local_counts[N == 1L, ..keycols]
author_unique_keys <- author_counts[N == 1L, ..keycols]
key_overlap <- merge(local_unique_keys, author_unique_keys, by = keycols)
write.csv(data.table(
  extract27_rows = nrow(d_ne), prepared_ne_rows = nrow(author_ne),
  extract27_key_rows = nrow(local_counts), prepared_ne_key_rows = nrow(author_counts),
  extract27_duplicate_key_groups = sum(local_counts$N > 1L),
  prepared_ne_duplicate_key_groups = sum(author_counts$N > 1L),
  unique_key_intersection = nrow(key_overlap),
  extract27_unique_keys = nrow(local_unique_keys),
  prepared_ne_unique_keys = nrow(author_unique_keys)
), file.path(outdir, "ne_key_coverage_summary.csv"), row.names = FALSE)

local_cmp <- d_ne[local_unique_keys, on = keycols, nomatch = 0L,
                  .(year, sample, serial, pernum, sex, age, ownershp)]
author_cmp <- author_ne[author_unique_keys, on = keycols, nomatch = 0L,
                        .(year, sample, serial, pernum, sex, age, ownershp)]
concordance <- merge(local_cmp, author_cmp, by = keycols, all = FALSE,
                     suffixes = c("_extract27", "_prepared_ne"))
concordance_summary <- rbindlist(lapply(c("sex", "age", "ownershp"), function(v) {
  a <- concordance[[paste0(v, "_extract27")]]
  b <- concordance[[paste0(v, "_prepared_ne")]]
  equal <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
  data.table(field = v, compared = length(a), equal = sum(equal), mismatch = sum(!equal),
             extract27_missing = sum(is.na(a)), prepared_ne_missing = sum(is.na(b)))
}))
write.csv(concordance_summary, file.path(outdir, "ne_sex_age_ownership_concordance.csv"),
          row.names = FALSE)

write.csv(d_ne[, .N, by = .(year, sample, rooms)][order(year, sample, rooms)],
          file.path(outdir, "ne_rooms_codes_by_year_sample.csv"), row.names = FALSE)
write.csv(d_ne[, .N, by = .(year, sample, bedrooms)][order(year, sample, bedrooms)],
          file.path(outdir, "ne_bedrooms_codes_by_year_sample.csv"), row.names = FALSE)

people <- unique(d_ne[, .(year, sample, serial, pernum)])
links <- unique(d_ne[momloc > 0 & !is.na(momloc), .(year, sample, serial, momloc)])
setnames(links, "momloc", "pernum")
setkeyv(people, c("year", "sample", "serial", "pernum"))
setkeyv(links, c("year", "sample", "serial", "pernum"))
links[, linked := FALSE]
links[people, on = keycols, linked := TRUE]
mothers <- unique(links[linked == TRUE, .(year, sample, serial, pernum)])
parents <- d_ne[mothers, on = keycols, nomatch = 0L]
parents[, proxy_roster := sex == 2 & age >= 18 & age <= 45 & nchild == 2 & yngch == 0 & eldch >= 1]
parents[, proxy_fertyr_yes := proxy_roster & fertyr == 2]
parents[, proxy_fertyr_known := proxy_roster & fertyr %in% c(1, 2)]
write.csv(parents[, .(
  linked_mother_rows = .N,
  roster_proxy_rows = sum(proxy_roster, na.rm = TRUE),
  roster_proxy_fertyr_yes = sum(proxy_fertyr_yes, na.rm = TRUE),
  roster_proxy_fertyr_known = sum(proxy_fertyr_known, na.rm = TRUE),
  roster_proxy_fertyr_missing = sum(proxy_roster & !fertyr %in% c(1, 2), na.rm = TRUE),
  rooms_positive_proxy = sum(proxy_roster & rooms_outcome > 0, na.rm = TRUE),
  bedrooms_authoritative_code_proxy = sum(proxy_roster & !is.na(bedrooms_outcome), na.rm = TRUE),
  bedrooms_zero_code_proxy = sum(proxy_roster & bedrooms_outcome == 0, na.rm = TRUE)
)], file.path(outdir, "ne_secondbirth_proxy_summary.csv"), row.names = FALSE)
checkpoint("proxy_summary")

manifest <- list(
  input = infile, input_size = expected_bytes, input_sha256 = expected_sha256,
  source_receipt = receipt_file, source_receipt_action = receipt$action,
  prepared_ne_source = prepared_file, states = ne_states,
  key = keycols, rows_ne = nrow(d_ne), prepared_ne_rows = nrow(author_ne),
  duplicate_key_groups_extract27 = sum(local_counts$N > 1L),
  duplicate_key_groups_prepared_ne = sum(author_counts$N > 1L),
  unique_key_intersection = nrow(key_overlap),
  unknown_housing_codes_national = unknown_total,
  unknown_housing_codes_ne = unknown_ne_total,
  housing_code_status = if (unknown_total > 0) "unknown_codes_preserved_as_missing_outcomes" else "all_observed_codes_authoritative",
  generated = as.character(Sys.time())
)
write_json(manifest, file.path(outdir, "source_audit_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
checkpoint("complete")
