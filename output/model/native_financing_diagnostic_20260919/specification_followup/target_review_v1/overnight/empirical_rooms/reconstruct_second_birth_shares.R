#!/usr/bin/env Rscript
# Rebuild only the eligible PSID rooms sample and count later biological births.
# This script fits no regression and writes only to its own output directory.
suppressPackageStartupMessages({library(haven); library(data.table); library(tidyselect)})

out <- normalizePath(dirname(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grep("^--file=", commandArgs(trailingOnly=FALSE))][1])), mustWork=TRUE)
root <- normalizePath(file.path(out, paste(rep("..", 7), collapse="/")), mustWork=TRUE)
source_file <- "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
builder <- file.path(root, "code/data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do")
receipt_file <- file.path(root, "code/data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1/target_receipt.csv")
stopifnot(file.exists(source_file), file.exists(builder), file.exists(receipt_file))

id_vars <- c("ID", "FID", "HHID", "year", "AGEREP", "EDUYEAR", "SEX", "DEATHYEAR", "REL", "CURRENT", "RELCHIREP", "ACTUALROOMS_", "IW")
child_vars <- unlist(lapply(1:20, function(k) c(paste0("RELCHI", k, "TYPE"), paste0("RELCHI", k, "BYEAR"), paste0("RELCHI", k, "ID"))))
vars <- c(id_vars, child_vars)
cat("Reading selected PSID shelf fields; no regression is run.\n", file=stderr())
x <- as.data.table(read_dta(source_file, col_select=all_of(vars)))
for (v in vars) set(x, j=v, value=as.numeric(x[[v]]))
nraw <- nrow(x)

# Reconstruct Stata's person-level first biological birth and reported-history
# flags using the same TYPE==1 and RELCHIREP definitions as the builder.
first_row <- rep(NA_real_, nraw)
for (k in 1:20) {
  typ <- x[[paste0("RELCHI", k, "TYPE")]]
  byr <- x[[paste0("RELCHI", k, "BYEAR")]]
  bio <- !is.na(typ) & typ == 1
  valid <- bio & !is.na(byr)
  y <- byr
  y[!valid] <- NA_real_
  new_first <- !is.na(y) & (is.na(first_row) | y < first_row)
  first_row[new_first] <- y[new_first]
}
x[, first_birth_row := first_row]
first_full <- x[, .(first_birth_year=if (all(is.na(first_birth_row))) NA_real_ else min(first_birth_row, na.rm=TRUE)), by=ID]
reported <- x[, .(max_children_reported=if (all(is.na(RELCHIREP))) NA_real_ else max(RELCHIREP, na.rm=TRUE)), by=ID]

# Current physical dwellings with multiple active family units are excluded,
# counting all current members before restricting to eligible women.
cur_fid <- unique(x[CURRENT == 1 & !is.na(HHID) & HHID > 0 & !is.na(FID) & FID > 0,
                    .(HHID, year, FID)])
n_fid <- cur_fid[, .(n_current_fids=.N), by=.(HHID, year)]
x <- merge(x, n_fid, by=c("HHID", "year"), all.x=TRUE, sort=FALSE)

# Match the builder's row ordering, one-wave rooms shift, year-specific codes,
# and exact eligible-woman / one-FID / one-woman-per-dwelling-year restrictions.
setorderv(x, c("ID", "year"))
x[, rooms := shift(ACTUALROOMS_), by=ID]
x[, rooms_source_year := shift(year), by=ID]
x[year <= 1984 & rooms == 9, rooms := NA_real_]
x[year >= 1985 & year <= 1993 & rooms == 99, rooms := NA_real_]
x[year >= 1994 & rooms %in% c(98, 99), rooms := NA_real_]
x <- x[is.na(DEATHYEAR) | year <= DEATHYEAR]
x <- x[CURRENT == 1 & !is.na(AGEREP) & AGEREP >= 18 & !is.na(IW) & IW > 0 &
       SEX == 2 & REL %in% c(1, 2) & !is.na(HHID) & HHID > 0 & !is.na(FID) & FID > 0]
x <- x[!is.na(n_current_fids) & n_current_fids == 1]
x[, rooms_alignment_gap_years := year - rooms_source_year]
x[!is.na(rooms) & !rooms_alignment_gap_years %in% c(1, 2), rooms := NA_real_]
x <- x[!is.na(rooms) & !is.na(EDUYEAR)]
x[, household_priority := as.integer(REL != 1)]
setorderv(x, c("HHID", "year", "household_priority", "ID"))
x <- x[, .SD[1], by=.(HHID, year)]

# Full-history extrema are computed over all source rows for each ID, matching
# egen over the unfiltered shelf, not over only the selected interview rows.
person <- x[, .(year_entry=min(year)), by=ID]
person <- merge(person, first_full, by="ID", all=TRUE)
person <- merge(person, reported, by="ID", all=TRUE)
rm(first_full); gc()

x <- merge(x, person, by="ID", all.x=FALSE, sort=FALSE)
x <- x[is.na(first_birth_year) | first_birth_year >= year_entry]
x[, `:=`(untimed_known_parent=is.na(first_birth_year) & !is.na(max_children_reported) & max_children_reported > 0,
         unknown_history=is.na(first_birth_year) & is.na(max_children_reported))]
x <- x[!untimed_known_parent & !unknown_history]
x[, K := year - first_birth_year]
input_rows <- nrow(x)
input_ids <- uniqueN(x$ID)
treated_ids <- uniqueN(x[!is.na(first_birth_year), ID])
never_ids <- uniqueN(x[is.na(first_birth_year) & max_children_reported == 0, ID])

# Reconstruct the standard two-way-FE singleton removal on the pre-fit rows.
# This uses only group counts, not the outcome or a regression.
pruned <- copy(x)
prune_steps <- list()
iter <- 0L
repeat {
  iter <- iter + 1L
  id_n <- pruned[, .N, by=ID]
  yr_n <- pruned[, .N, by=year]
  drop_ids <- id_n[N == 1, ID]
  drop_years <- yr_n[N == 1, year]
  if (length(drop_ids) == 0L && length(drop_years) == 0L) break
  n_before <- nrow(pruned)
  pruned <- pruned[!ID %in% drop_ids & !year %in% drop_years]
  prune_steps[[iter]] <- list(iteration=iter, singleton_person_ids=length(drop_ids),
                              singleton_years=length(drop_years), rows_removed=n_before-nrow(pruned))
}
pruned_rows <- nrow(pruned)
pruned_ids <- uniqueN(pruned$ID)
pruned_never_ids <- uniqueN(pruned[is.na(first_birth_year) & max_children_reported == 0, ID])
pruned_treated_ids <- uniqueN(pruned[!is.na(first_birth_year), ID])

# Build child-year inventory only for eligible treated women. Verify relation
# slot identity against the shelf's explicit child IDs before using any shares.
treated <- unique(x[!is.na(first_birth_year), .(ID, first_birth_year)])
ids <- treated$ID
bp <- as.data.table(read_dta(source_file, col_select=all_of(c("ID", child_vars))))
bp[, ID := as.numeric(ID)]
bp <- bp[ID %in% ids]
slot_rows <- vector("list", 20L)
identity_rows <- vector("list", 20L)
for (k in 1:20) {
  tv <- as.numeric(bp[[paste0("RELCHI", k, "TYPE")]])
  by <- as.numeric(bp[[paste0("RELCHI", k, "BYEAR")]])
  child_id <- as.numeric(bp[[paste0("RELCHI", k, "ID")]])
  z <- data.table(ID=bp$ID, slot=k, child_id=child_id, bio=(!is.na(tv) & tv == 1), byear=by)
  z <- z[bio == TRUE]
  slot_rows[[k]] <- z[, .(birth_year=if (all(is.na(byear))) NA_real_ else min(byear, na.rm=TRUE),
                         type1_year_missing=all(is.na(byear))), by=.(ID, slot)]
  identity_rows[[k]] <- z[, .(ID, slot, child_id, byear)]
}
slots <- rbindlist(slot_rows, use.names=TRUE)
bio_identity <- rbindlist(identity_rows, use.names=TRUE)
slot_identity <- bio_identity[!is.na(child_id),
  .(child_id_count=uniqueN(child_id)), by=.(ID, slot)]
child_slot_identity <- bio_identity[!is.na(child_id),
  .(slot_count=uniqueN(slot), birth_year_count=uniqueN(byear[!is.na(byear)])), by=.(ID, child_id)]
identity_audit <- list(
  child_id_variables=20L,
  documented_labels="RELCHI#ID fields are labeled `Ind's child #, unique ID` in the local PSID shelf.",
  treated_biological_slot_id_missing_with_year=bio_identity[is.na(child_id) & !is.na(byear), .N],
  treated_id_slot_pairs_with_multiple_child_ids=slot_identity[child_id_count > 1, .N],
  treated_child_ids_seen_in_multiple_slots=child_slot_identity[slot_count > 1, .N],
  treated_child_ids_with_conflicting_birth_years=child_slot_identity[birth_year_count > 1, .N],
  treated_child_id_observations=bio_identity[!is.na(child_id), .N]
)
coverage <- slots[, .(bio_slots_seen=.N, bio_slots_timed=sum(!is.na(birth_year)),
                      bio_slots_untimed=sum(is.na(birth_year))), by=ID]
dated <- slots[!is.na(birth_year)]
hist_child <- dated[, .(first_from_slots=min(birth_year),
                        second_child_year=if (.N >= 2) sort(birth_year)[2] else NA_real_), by=ID]
hist_later <- dated[, .(second_later_year=if (uniqueN(birth_year) >= 2) sort(unique(birth_year))[2] else NA_real_), by=ID]
hist <- merge(hist_child, hist_later, by="ID", all=TRUE)
hist <- merge(hist, coverage, by="ID", all=TRUE)
hist <- merge(hist, reported, by="ID", all=TRUE)
hist[, `:=`(reported_children_at_least2=!is.na(max_children_reported) & max_children_reported >= 2,
            incomplete_bio_history=(bio_slots_untimed > 0))]
treated <- merge(treated, hist, by="ID", all.x=TRUE)
if (any(treated$first_birth_year != treated$first_from_slots, na.rm=TRUE))
  stop("Child-slot inventory does not reproduce builder's first biological birth year")
id_history <- bio_identity[!is.na(child_id), .(
  birth_year=if (all(is.na(byear))) NA_real_ else min(byear, na.rm=TRUE)), by=.(ID, child_id)]
id_history <- id_history[!is.na(birth_year)]
id_hist_summary <- id_history[, .(first_from_ids=min(birth_year),
  second_child_year_id=if (.N >= 2) sort(birth_year)[2] else NA_real_,
  second_later_year_id=if (uniqueN(birth_year) >= 2) sort(unique(birth_year))[2] else NA_real_), by=ID]
treated <- merge(treated, id_hist_summary, by="ID", all.x=TRUE)
identity_audit$firstbirth_slot_vs_childid_mismatches <- treated[
  xor(is.na(first_from_ids), is.na(first_from_slots)) |
    (!is.na(first_from_ids) & !is.na(first_from_slots) & first_from_slots != first_from_ids), .N]
identity_audit$second_child_year_slot_vs_childid_mismatches <- treated[
  xor(is.na(second_child_year_id), is.na(second_child_year)) |
    (!is.na(second_child_year_id) & !is.na(second_child_year) & second_child_year != second_child_year_id), .N]
identity_audit$later_distinct_year_slot_vs_childid_mismatches <- treated[
  xor(is.na(second_later_year_id), is.na(second_later_year)) |
    (!is.na(second_later_year_id) & !is.na(second_later_year) & second_later_year != second_later_year_id), .N]
identity_audit$share_identity_usable <- (identity_audit$treated_biological_slot_id_missing_with_year == 0 &&
  identity_audit$treated_id_slot_pairs_with_multiple_child_ids == 0 &&
  identity_audit$treated_child_ids_seen_in_multiple_slots == 0 &&
  identity_audit$treated_child_ids_with_conflicting_birth_years == 0 &&
  identity_audit$firstbirth_slot_vs_childid_mismatches == 0 &&
  identity_audit$second_child_year_slot_vs_childid_mismatches == 0 &&
  identity_audit$later_distinct_year_slot_vs_childid_mismatches == 0)

shares <- rbindlist(lapply(c(3,4), function(h) {
  cohort <- treated[, .(ID, first_birth_year, second_child_year,
                        second_later_year, second_child_year_id,
                        second_later_year_id, reported_children_at_least2,
                        incomplete_bio_history)]
  cohort[, `:=`(child_by_cutoff=!is.na(second_child_year) & second_child_year <= first_birth_year+h,
                later_birth_by_cutoff=!is.na(second_later_year) & second_later_year <= first_birth_year+h,
                child_id_by_cutoff=!is.na(second_child_year_id) & second_child_year_id <= first_birth_year+h,
                later_birth_id_by_cutoff=!is.na(second_later_year_id) & second_later_year_id <= first_birth_year+h)]
  at_event <- unique(x[K == h & !is.na(first_birth_year), ID])
  obs <- cohort[ID %in% at_event]
  at_event_pruned <- unique(pruned[K == h & !is.na(first_birth_year), ID])
  obs_pruned <- cohort[ID %in% at_event_pruned]
  cohort_pruned <- cohort[ID %in% unique(pruned[!is.na(first_birth_year), ID])]
  alln <- nrow(cohort); en <- nrow(obs)
  alln_pruned <- nrow(cohort_pruned)
  data.table(
    event_horizon_years=h,
    denominator_all_eligible_treated=alln,
    documented_second_child_by_h_all=sum(cohort$child_by_cutoff, na.rm=TRUE),
    share_second_child_by_h_all=sum(cohort$child_by_cutoff, na.rm=TRUE)/alln,
    documented_later_distinct_birth_by_h_all=sum(cohort$later_birth_by_cutoff, na.rm=TRUE),
    share_later_distinct_birth_by_h_all=sum(cohort$later_birth_by_cutoff, na.rm=TRUE)/alln,
    denominator_singleton_pruned_treated=alln_pruned,
    documented_later_distinct_birth_by_h_all_pruned=sum(cohort_pruned$later_birth_by_cutoff, na.rm=TRUE),
    share_later_distinct_birth_by_h_all_pruned=sum(cohort_pruned$later_birth_by_cutoff, na.rm=TRUE)/alln_pruned,
    event_time_room_observed_women=en,
    documented_second_child_by_h_event=sum(obs$child_by_cutoff, na.rm=TRUE),
    share_second_child_by_h_event=if (en) sum(obs$child_by_cutoff, na.rm=TRUE)/en else NA_real_,
    documented_later_distinct_birth_by_h_event=sum(obs$later_birth_by_cutoff, na.rm=TRUE),
    share_later_distinct_birth_by_h_event=if (en) sum(obs$later_birth_by_cutoff, na.rm=TRUE)/en else NA_real_,
    documented_second_child_by_h_childid_event=sum(obs$child_id_by_cutoff, na.rm=TRUE),
    share_second_child_by_h_childid_event=if (en) sum(obs$child_id_by_cutoff, na.rm=TRUE)/en else NA_real_,
    documented_later_distinct_birth_by_h_childid_event=sum(obs$later_birth_id_by_cutoff, na.rm=TRUE),
    share_later_distinct_birth_by_h_childid_event=if (en) sum(obs$later_birth_id_by_cutoff, na.rm=TRUE)/en else NA_real_,
    event_group_history_incomplete=sum(obs$incomplete_bio_history, na.rm=TRUE),
    event_group_reported_children_ge2=sum(obs$reported_children_at_least2, na.rm=TRUE),
    singleton_pruned_reconstruction_women=length(at_event_pruned),
    documented_second_child_by_h_pruned=sum(obs_pruned$child_by_cutoff, na.rm=TRUE),
    share_second_child_by_h_pruned=if (length(at_event_pruned)) sum(obs_pruned$child_by_cutoff, na.rm=TRUE)/length(at_event_pruned) else NA_real_,
    documented_later_distinct_birth_by_h_pruned=sum(obs_pruned$later_birth_by_cutoff, na.rm=TRUE),
    share_later_distinct_birth_by_h_pruned=if (length(at_event_pruned)) sum(obs_pruned$later_birth_by_cutoff, na.rm=TRUE)/length(at_event_pruned) else NA_real_
  )
}))

target <- fread(receipt_file)
metadata <- jsonlite::fromJSON(file.path(dirname(receipt_file), "metadata.json"))
checks <- list(
  source_path=source_file,
  source_sha256_recorded_in_target_metadata=metadata$source$sha256,
  builder_path=builder,
  source_sha256_recomputed=FALSE,
  raw_rows=nraw,
  reconstructed_prefit_input_rows=input_rows,
  target_receipt_prefit_input_rows=as.integer(target$input_observations),
  reconstructed_prefit_input_ids=input_ids,
  target_receipt_prefit_input_ids=as.integer(target$input_individuals),
  reconstructed_prefit_treated_ids=treated_ids,
  target_receipt_prefit_treated_ids=as.integer(target$treated_individuals),
  reconstructed_never_treated_ids=never_ids,
  target_receipt_confirmed_never_ids=as.integer(target$est_confirmed_never_ids),
  target_estimation_rows=as.integer(target$estimation_observations),
  target_estimation_ids=as.integer(target$estimation_individuals),
  singleton_pruning=list(iterations=prune_steps,
                         reconstructed_rows=pruned_rows,
                         reconstructed_ids=pruned_ids,
                         reconstructed_never_treated_ids=pruned_never_ids,
                         reconstructed_treated_ids=pruned_treated_ids,
                         target_rows=as.integer(target$estimation_observations),
                         target_ids=as.integer(target$estimation_individuals),
                         target_confirmed_never_ids=as.integer(target$est_confirmed_never_ids),
                         counts_match=(pruned_rows == as.integer(target$estimation_observations) &&
                           pruned_ids == as.integer(target$estimation_individuals) &&
                           pruned_never_ids == as.integer(target$est_confirmed_never_ids))),
  prefit_count_consistency=list(rows=input_rows == as.integer(target$input_observations),
                               ids=input_ids == as.integer(target$input_individuals),
                               treated_ids=treated_ids == as.integer(target$treated_individuals)),
  event_window_denominator_definition="Unique treated women with valid rooms and covariates at exact K=+3 or K=+4 after reproducing the target builder's iterative person/year-FE singleton pruning. This matches saved e(sample) row/person/control counts, but no individual e(sample) marker is saved, so identity is not asserted from count agreement alone.",
  e_sample_singleton_reconstruction=prune_steps,
  child_identity_audit=identity_audit,
  candidate_quantitative_status=if (isTRUE(identity_audit$share_identity_usable) &&
    pruned_rows == as.integer(target$estimation_observations) &&
    pruned_ids == as.integer(target$estimation_individuals) &&
    pruned_never_ids == as.integer(target$est_confirmed_never_ids))
    "Candidate descriptive evidence: child IDs/slot histories reconcile and singleton-pruned aggregate sample counts match; no individual e(sample) marker is saved."
  else "Provisional only: child identity or singleton-pruned sample counts failed a consistency check.",
  birth_definition="Biological child TYPE==1. Second child counts a second biological child even if same recorded birth year; later distinct birth requires a strictly later distinct recorded year. No birth-order tie is imputed.",
  limitations="RELCHIREP explicitly counts children with or without records. Documented second-birth shares are observed lower bounds when second-child birth years are absent; table reports sample counts with missing biological birth records and reported child count >=2."
)
fwrite(shares, file.path(out, "second_birth_event_shares.csv"))
jsonlite::write_json(checks, file.path(out, "second_birth_sample_receipt.json"), auto_unbox=TRUE, pretty=TRUE, na="null")
cat(jsonlite::toJSON(list(receipt=checks, shares=shares), auto_unbox=TRUE, pretty=TRUE, na="null"))
