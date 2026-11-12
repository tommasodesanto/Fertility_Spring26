#!/usr/bin/env Rscript
# Compares the design-B Weekly ACS-run panel from a true restricted rematch with
# the filter-and-reweight design-B table for the same state. Read-only on inputs.
options(stringsAsFactors = FALSE, scipen = 999); suppressPackageStartupMessages(library(data.table))
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
rp <- readRDS(Sys.getenv("REMATCH_PANEL")); setDT(rp)
rp <- rp[as.character(matching_sample) == "Weekly" & as.character(dataset) == "ACS"]
rp[, t := num(t_es_lw)]; rp[, child := as.integer(num(child))]
anc <- rp[child == 1L & t == 0, .(id = num(id), src_key)]
stopifnot(!anyDuplicated(anc$id))
rt <- rp[, .(src_key, child, t, wgt = num(wgt),
             anchor_src = ifelse(child == 0L, anc$src_key[match(num(match_id), anc$id)], NA_character_))]
ft <- readRDS(Sys.getenv("FILTER_TABLE")); setDT(ft)
key <- c("src_key", "child", "t", "anchor_src")
rt[, dup := .N, by = key]; ft[, dup := .N, by = key]
mg <- merge(rt, ft, by = key, all = TRUE, suffixes = c("_rematch", "_filter"))
out <- list(rows_rematch = nrow(rt), rows_filter = nrow(ft),
            only_rematch = sum(is.na(mg$wgt_filter)), only_filter = sum(is.na(mg$wgt_rematch)),
            dup_keys_rematch = sum(rt$dup > 1), dup_keys_filter = sum(ft$dup > 1),
            max_abs_wgt_diff = max(abs(mg$wgt_rematch - mg$wgt_filter), na.rm = TRUE),
            max_rel_wgt_diff = max(abs(mg$wgt_rematch - mg$wgt_filter) / mg$wgt_rematch, na.rm = TRUE),
            donor_rows_rematch = sum(rt$child == 0L), donor_rows_filter = sum(ft$child == 0L))
print(out)
jsonlite::write_json(out, Sys.getenv("COMPARE_OUT"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
if (out$only_rematch + out$only_filter > 0) fwrite(mg[is.na(wgt_filter) | is.na(wgt_rematch)][1:min(.N, 200)], sub("json$", "mismatch_head.csv", Sys.getenv("COMPARE_OUT")))
