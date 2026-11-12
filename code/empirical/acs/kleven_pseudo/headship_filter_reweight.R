#!/usr/bin/env Rscript
# Household-position designs for the ACS first-birth housing pseudo-panel
# (2026-11-12). One state per call (STATEFIP). Reads the saved full matched
# state panel from the national production run and the household-position key,
# filters people by household position, and recomputes donor matching weights
# 1/k' within each (state, matching run, anchor mother, lag) group.
#
# Equivalence argument: Kleven's matcher (matching.R, sha 7433ec92) matches
# each event-0 mother at lag n exactly, with all ties, to every childless donor
# in her cell, at weight 1/k; donor wgt = wgt_original * 1/k; mothers keep
# their survey weight and stay in the panel even when unmatched. Cells are
# person-level, so dropping people before matching equals dropping them here
# and recounting k. k' is counted on the full Weekly panel before any housing
# missingness. CPS rows, the cleaning rules and min(doiy)+5 are untouched.
#
# Writes only under OUTDIR/states/statefip_XX. Never modifies inputs.
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(data.table))
stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
root <- Sys.getenv("PROJECT_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
prod <- Sys.getenv("PRODUCTION_DIR", file.path(root, "output/national_acs_production_20260921a"))
keydir <- Sys.getenv("KEY_DIR")
outroot <- Sys.getenv("OUTDIR")
st <- as.integer(Sys.getenv("STATEFIP", Sys.getenv("SLURM_ARRAY_TASK_ID")))
hh_gq_file <- Sys.getenv("HH_GQ_FILE", "/scratch/td2248/projects/kleven_acs_pilot_20260917/code_headship_20261112/hh_gq.txt")
hh_gq_txt <- if (file.exists(hh_gq_file)) trimws(readLines(hh_gq_file, warn = FALSE)[1]) else Sys.getenv("HH_GQ", "1,2")
hh_gq <- as.numeric(strsplit(hh_gq_txt, ",")[[1]])
stopifnot(length(hh_gq) >= 2L, !anyNA(hh_gq))
panel_file <- Sys.getenv("PANEL_FILE", file.path(prod, sprintf("statefip_%02d", st), "cps_acs_pseudo-panel_housing.rds"))
stopifnot(nzchar(keydir), nzchar(outroot), is.finite(st), file.exists(panel_file))
out <- file.path(outroot, "states", sprintf("statefip_%02d", st))
if (dir.exists(out) && length(list.files(out))) stop("state OUTDIR non-empty; refusing overwrite: ", out)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fail <- function(msg) {
  jsonlite::write_json(list(status = "FAILED", statefip = st, error = msg, generated = stamp()),
                       file.path(out, "failure.json"), auto_unbox = TRUE, pretty = TRUE)
  stop(msg, call. = FALSE)
}
req <- function(ok, msg) if (!isTRUE(ok)) fail(msg)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

message(stamp(), " state ", st, " read ", panel_file)
p <- readRDS(panel_file)
setDT(p)
p <- p[as.character(matching_sample) == "Weekly"]
keep_cols <- c("wgt", "wgt_original", "wgt_match", "id", "match_id", "child", "eldch", "age1b",
               "age_factor", "doiy_factor", "statefip", "gender", "married", "dataset",
               "rooms_raw", "ownershp_raw", "bedrooms_raw", "source_origin", "from_cps",
               "source_year", "source_sample", "source_serial", "source_pernum", "source_hh_cluster",
               "src_key", "matching_sample", "t_es_lw")
req(all(keep_cols %in% names(p)), paste("panel lacks:", paste(setdiff(keep_cols, names(p)), collapse = ",")))
p <- p[, ..keep_cols]; invisible(gc())
p[, run := as.character(dataset)]
req(all(p$run %in% c("ACS", "CPS")), paste("unexpected dataset values:", paste(unique(p$run), collapse = ",")))
p[, t := num(t_es_lw)]
req(!anyNA(p$t), "Weekly rows with missing event time")
p[, child := as.integer(num(child))]
req(all(p$child %in% c(0L, 1L)), "child not 0/1")
p[, true_acs := as.character(source_origin) == "ACS" & num(from_cps) == 0]
p[, row := .I]

# ---- household-position key -------------------------------------------------
key <- readRDS(file.path(keydir, sprintf("statefip_%02d.rds", st))); setDT(key)
p[, KEY := ifelse(true_acs, paste(num(source_sample), num(source_year), num(source_serial), num(source_pernum), sep = ":"), NA_character_)]
m <- match(p$KEY, key$KEY)
p[, covered := true_acs & !is.na(m)]
p[, `:=`(GQ = key$GQ[m], RELATE = key$RELATE[m], RELATED = key$RELATED[m],
         MOMLOC = key$MOMLOC[m], POPLOC = key$POPLOC[m])]
rm(key); invisible(gc())
# housing observed only on covered rows: assert the bridge agrees
req(!any(!p$covered & !is.na(num(p$rooms_raw))), "rows with rooms but no position key")
p[, hh := covered & GQ %in% hh_gq]
p[, headspouse := covered & RELATE %in% c(1, 2)]
p[, partner := covered & RELATED == 1114]
# eligibility per design; rows without a key (CPS, ACS outside 2005-19) are
# left unfiltered: they carry no housing outcome.
p[, elig_B := !covered | (hh & (headspouse | partner))]
p[, elig_C := !covered | (hh & headspouse)]
p[, elig_D := !covered | hh]

# ---- anchors and donor groups --------------------------------------------------
p[, anchor := child == 1L & t == 0]
p[, grp := ifelse(child == 0L, paste(run, num(match_id), t, sep = "|"), NA_character_)]
p[, anchor_key := ifelse(child == 0L, paste(run, num(match_id), sep = "|"),
                         ifelse(anchor, paste(run, num(id), sep = "|"), NA_character_))]
acs_anchor <- p[anchor & run == "ACS"]
req(!anyDuplicated(acs_anchor$anchor_key), "ACS-run anchor id not unique within state")
req(!anyDuplicated(acs_anchor$src_key), "ACS-run anchor src_key not unique within state")
don <- p[child == 0L]
don[, k := .N, by = grp]
# donors must resolve to an anchor mother in the same run (ACS run checked)
acs_don <- don[run == "ACS"]
miss_anchor <- sum(!acs_don$anchor_key %in% acs_anchor$anchor_key)
req(miss_anchor == 0L, paste("ACS-run donors without anchor mother:", miss_anchor))
# preflight: stored matching weight equals 1/k; stored wgt equals wgt_original/k
dev_match <- max(abs(num(acs_don$wgt_match) - 1 / acs_don$k))
dev_wgt <- max(abs(num(acs_don$wgt) - num(acs_don$wgt_original) / acs_don$k) / num(acs_don$wgt))
dev_match_cps <- if (nrow(don[run == "CPS"])) max(abs(num(don[run == "CPS"]$wgt_match) - 1 / don[run == "CPS"]$k)) else NA_real_
preflight <- list(statefip = st, acs_run_donor_rows = nrow(acs_don), acs_run_anchors = nrow(acs_anchor),
                  acs_run_groups = uniqueN(acs_don$grp), max_abs_dev_wgt_match = dev_match,
                  max_rel_dev_wgt = dev_wgt, cps_run_max_abs_dev_wgt_match = dev_match_cps,
                  pass = dev_match < 1e-9 && dev_wgt < 1e-9)
jsonlite::write_json(preflight, file.path(out, "preflight.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
req(isTRUE(preflight$pass), sprintf("preflight failed: dev_match=%g dev_wgt=%g", dev_match, dev_wgt))
rm(acs_don); invisible(gc())

# anchor eligibility carried to donors
anc_elig <- p[anchor == TRUE, .(anchor_key, a_B = elig_B, a_C = elig_C, a_D = elig_D, a_age1b = num(age1b))]
anc_elig <- anc_elig[!duplicated(anchor_key)]

pool_cols <- c("wgt", "age_factor", "doiy_factor", "statefip", "gender", "event_time",
               "rooms_raw", "ownershp_raw", "bedrooms_raw", "source_origin", "from_cps",
               "source_year", "source_sample", "source_serial", "source_pernum", "source_hh_cluster",
               "matching_sample", "t_es_lw", "hs_role", "hs_anchor", "hs_src")

build <- function(elig_m, elig_d, anchor_col, age_rule = FALSE) {
  # mothers: all event times, survey weights unchanged
  mo <- p[child == 1L & (if (is.null(elig_m)) TRUE else get(elig_m))]
  if (age_rule) mo <- mo[num(age1b) >= 30 & num(age1b) <= 45]
  # donors: drop if donor fails or its anchor mother fails (ACS run only), recount k'
  d <- don[if (is.null(elig_d)) TRUE else get(elig_d)]
  if (!is.null(anchor_col)) {
    ae <- anc_elig[match(d$anchor_key, anc_elig$anchor_key)][[anchor_col]]
    d <- d[run == "CPS" | is.na(ae) | ae]
  }
  if (age_rule) d <- d[num(age1b) >= 30 & num(age1b) <= 45]
  d[, kp := .N, by = grp]
  d[, wgt := num(wgt_original) / kp]
  z <- rbind(mo, d[, names(mo), with = FALSE])
  z
}
narrow <- function(z) {
  z <- z[true_acs & as.character(gender) == "Women"]
  z[, event_time := t_es_lw]
  z[, hs_role := ifelse(child == 1L, "mother", "donor")]
  z[, hs_anchor := paste(st, anchor_key, sep = "|")]
  z[, hs_src := src_key]
  as.data.frame(z[, ..pool_cols])
}
fullpre <- function(z) {
  # anchors with at least one surviving donor in housing coverage at every lag 1..5
  sup <- z[child == 0L & run == "ACS" & covered, .(nl = uniqueN(t)), by = anchor_key]
  ok <- sup[nl == 5L, anchor_key]
  z[(child == 1L & t > 0) | (child == 1L & t == 0 & anchor_key %in% ok) | (child == 0L & anchor_key %in% ok)]
}
designs <- list(
  A = function() build(NULL, NULL, NULL),
  B = function() build("elig_B", "elig_B", "a_B"),
  C = function() build("elig_C", "elig_C", "a_C"),
  D = function() build("elig_D", "elig_D", "a_D"),
  E = function() build(NULL, "elig_B", NULL),
  A_age3045 = function() build(NULL, NULL, NULL, age_rule = TRUE),
  B_age3045 = function() build("elig_B", "elig_B", "a_B", age_rule = TRUE)
)
counts <- list(); support <- list()
anchors_A <- p[anchor == TRUE & run == "ACS" & as.character(gender) == "Women", .(anchor_key, wa = num(wgt), a_B = elig_B, a_C = elig_C, a_D = elig_D)]
donA <- don[run == "ACS", .(nA = .N, covA = sum(covered)), by = .(anchor_key, t)]
for (nm in c(names(designs), "A_fullpre", "B_fullpre")) {
  z <- if (nm == "A_fullpre") fullpre(designs$A()) else if (nm == "B_fullpre") fullpre(designs$B()) else designs[[nm]]()
  if (nm %in% c("A", "B", "C", "D", "E")) {
    # support by lag against the original anchor list
    zd <- z[child == 0L & run == "ACS", .(nS = .N, covS = sum(covered)), by = .(anchor_key, t)]
    acol <- switch(nm, B = "a_B", C = "a_C", D = "a_D", NULL)
    for (lag in 1:5) {
      a <- copy(anchors_A)
      a[, inA := anchor_key %in% donA[t == -lag & nA > 0, anchor_key]]
      a[, kept := if (is.null(acol)) TRUE else get(acol)]
      a[, survS := anchor_key %in% zd[t == -lag & nS > 0, anchor_key]]
      a[, covS := anchor_key %in% zd[t == -lag & covS > 0, anchor_key]]
      a[, cat := fifelse(!inA, "absent_in_A", fifelse(!kept, "anchor_excluded",
                 fifelse(!survS, "all_donors_fail", fifelse(!covS, "calendar_no_housing_coverage", "supported"))))]
      support[[length(support) + 1L]] <- a[, .(anchors = .N, anchor_wgt = sum(wa)), by = cat][, `:=`(design = nm, lag = lag, statefip = st)]
    }
  }
  if (nm == "B" && identical(Sys.getenv("SAVE_MATCH_TABLE"), "1")) {
    mt <- z[run == "ACS", .(src_key, child, t, wgt = num(wgt),
                            anchor_src = acs_anchor$src_key[match(anchor_key, acs_anchor$anchor_key)])]
    saveRDS(mt, file.path(out, "match_table_B.rds"))
  }
  nz <- narrow(z)
  saveRDS(nz, file.path(out, paste0("design_", nm, ".rds")), compress = FALSE)
  ze <- z[true_acs & as.character(gender) == "Women"]
  counts[[nm]] <- data.table(design = nm, statefip = st, rows = nrow(nz),
                             mother_rows = sum(ze$child == 1L), donor_rows = sum(ze$child == 0L),
                             anchors = sum(ze$anchor), donor_rows_rooms_missing = sum(ze$child == 0L & ze$covered & is.na(num(ze$rooms_raw))),
                             unique_persons = uniqueN(ze$src_key), unique_households = uniqueN(ze$source_hh_cluster))
  rm(z, nz, ze); invisible(gc())
}
fwrite(rbindlist(counts), file.path(out, "counts.csv"))
fwrite(rbindlist(support, fill = TRUE), file.path(out, "support_by_lag.csv"))

# ---- diagnostics on baseline A: household position by event time --------------------
A <- designs$A()[true_acs & as.character(gender) == "Women" & covered]
A[, pos := fifelse(!(GQ %in% hh_gq), "group_quarters",
           fifelse(RELATE == 1, "householder", fifelse(RELATE == 2, "spouse",
           fifelse(RELATED == 1114, "unmarried_partner", fifelse(RELATE == 3, "child_of_householder",
           fifelse(RELATE %in% 4:10, "other_relative", "other_nonrelative"))))))]
A[, parent_coresident := (num(MOMLOC) > 0) | (num(POPLOC) > 0)]
A[, rooms9 := pmin(fifelse(num(rooms_raw) %in% c(1:27, 30), num(rooms_raw), NA_real_), 9)]
A[, own := fifelse(num(ownershp_raw) %in% c(1, 2), as.numeric(num(ownershp_raw) == 1), NA_real_)]
A[, role := fifelse(child == 1L, "mother", "donor")]
A[, married_flag := as.character(married)]
diag <- A[, .(n = .N, w = sum(wgt), w_parent = sum(wgt * parent_coresident),
              w_rooms = sum(wgt[!is.na(rooms9)]), wy_rooms = sum((wgt * rooms9)[!is.na(rooms9)]),
              w_own = sum(wgt[!is.na(own)]), wy_own = sum((wgt * own)[!is.na(own)])),
          by = .(t, role, married_flag, pos)][, statefip := st]
fwrite(diag, file.path(out, "diag_positions.csv"))
jsonlite::write_json(list(status = "STATE_HEADSHIP_COMPLETE", statefip = st, panel = panel_file,
                          hh_gq = hh_gq, preflight = preflight, generated = stamp()),
                     file.path(out, "state_receipt.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(stamp(), " STATE_HEADSHIP_COMPLETE ", st)
