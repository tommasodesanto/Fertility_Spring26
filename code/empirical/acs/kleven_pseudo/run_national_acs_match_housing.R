#!/usr/bin/env Rscript
# State-level national ACS/CPS matcher and housing bridge.
# Consumes state partitions from run_national_acs_source_stage.R, evaluates the
# unchanged Kleven cleaner/matcher in an isolated work directory, and carries
# source identity through the vendor select statements. No job submission,
# downloads, or outcome fitting occurs here.
options(stringsAsFactors = FALSE, scipen = 999)
started <- Sys.time(); stamp0 <- format(started, "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
root <- Sys.getenv("PROJECT_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
source_out <- Sys.getenv("PARTITION_OUTDIR", file.path(root, "output", "national_acs_source_stage"))
outdir <- Sys.getenv("OUTDIR", file.path(root, "output", "national_acs_match_housing"))
phase <- Sys.getenv("PHASE", "preflight")
state_text <- Sys.getenv("STATEFIP", if (phase == "smoke") "50" else "1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56")
states <- suppressWarnings(as.integer(strsplit(state_text, ",", fixed=TRUE)[[1]]))
valid_states <- c(1L,2L,4L,5L,6L,8L,9L,10L,11L,12L,13L,15L,16L,17L,18L,19L,20L,21L,22L,23L,24L,25L,26L,27L,28L,29L,30L,31L,32L,33L,34L,35L,36L,37L,38L,39L,40L,41L,42L,44L,45L,46L,47L,48L,49L,50L,51L,53L,54L,55L,56L)
fail <- function(msg, stage="startup") {
  rec <- list(status="FAILED", stage=stage, error=msg, started=stamp0, ended=stamp())
  if (dir.exists(outdir) && requireNamespace("jsonlite", quietly=TRUE)) try(jsonlite::write_json(rec, file.path(outdir,"failure_receipt.json"), auto_unbox=TRUE, pretty=TRUE), silent=TRUE)
  stop(msg, call.=FALSE)
}
req <- function(ok, msg, stage="startup") if (!isTRUE(ok)) fail(msg, stage)
req(length(states)>0L && all(is.finite(states) & states %in% valid_states), "STATEFIP contains an invalid FIPS")
if (identical(phase, "production")) req(identical(Sys.getenv("RUN_PRODUCTION"), "YES"), "production requires RUN_PRODUCTION=YES")
if (dir.exists(outdir) && length(list.files(outdir, all.files=TRUE, no..=TRUE))) fail("OUTDIR exists and is non-empty; refusing overwrite")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
vendor_dir <- Sys.getenv("VENDOR_DIR", file.path(root,"vendor"))
req(dir.exists(vendor_dir), paste("vendor directory absent:", vendor_dir))
contract_file <- Sys.getenv("SOURCE_CONTRACT", file.path(root,"code","empirical","acs","kleven_pseudo","source_contract.json"))
req(file.exists(contract_file), paste("canonical source contract absent:", contract_file))
contract <- jsonlite::fromJSON(contract_file, simplifyVector=FALSE)
vendor_expected <- unlist(contract$vendor_sha256, use.names=TRUE)
req(identical(sort(names(vendor_expected)), sort(c("clean_acs.R","clean_cps.R","matching.R","setup.R","functions.R"))), "canonical source contract vendor set is incomplete")
req(requireNamespace("digest", quietly=TRUE) && requireNamespace("jsonlite", quietly=TRUE), "digest/jsonlite required")
vh <- vapply(names(vendor_expected), function(nm) { p <- file.path(vendor_dir,nm); req(file.exists(p), paste("missing vendor file", p)); digest::digest(file=p, algo="sha256") }, character(1))
req(all(vh == vendor_expected), paste("vendor hash mismatch:", paste(names(vh)[vh != vendor_expected], collapse=", ")))
packages <- c("Matching","plyr","data.table","dplyr","tidyr","forcats","purrr","labelled","ipumsr","haven","stringr","here")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
req(!length(missing), paste("missing runtime packages; no install is attempted:", paste(missing, collapse=", ")))
suppressPackageStartupMessages({library(Matching); library(plyr); library(data.table); library(dplyr); library(tidyr); library(forcats); library(purrr); library(labelled); library(ipumsr); library(haven); library(stringr); library(here)})
partition_receipt <- file.path(source_out, "partition_receipt.json"); schema_receipt <- file.path(source_out, "source_schema_receipt.json")
req(file.exists(partition_receipt) && file.exists(schema_receipt), "source partition/schema receipts absent", "source_gate")
pr <- jsonlite::fromJSON(partition_receipt, simplifyVector=FALSE); sr <- jsonlite::fromJSON(schema_receipt, simplifyVector=FALSE)
req(identical(pr$status, "PARTITION_COMPLETE"), paste("source stage status is", pr$status), "source_gate")
req(identical(sr$status, "SCHEMA_PASS"), paste("source schema status is", sr$status), "source_gate")
for (st in states) { sd <- file.path(source_out,"partitions",sprintf("statefip_%02d",st)); req(file.exists(file.path(sd,"acs_raw.RData")) && file.exists(file.path(sd,"cps_raw.RData")), paste("partition missing for FIPS",st), "source_gate") }
writeLines(paste(stamp(), "PHASE", phase, "start", paste(states, collapse=",")), file.path(outdir,"progress.log"))
logp <- function(...) { z <- paste(stamp(), ..., collapse=" "); cat(z,"\n",file=file.path(outdir,"progress.log"),append=TRUE); message(z) }
resolve <- function(nms, want, required=TRUE) { h <- nms[tolower(nms) == tolower(want)]; if (!length(h) && !required) return(NA_character_); req(length(h)==1L, paste("field",want,"resolved",length(h),"times"), "source_adapter"); h[[1L]] }
load_data <- function(path) { e <- new.env(parent=emptyenv()); load(path,envir=e); req("data" %in% ls(e), paste("data object absent",path),"source_adapter"); e$data }
add_lineage <- function(d, origin) {
  n <- names(d); cy <- resolve(n,"year"); cs <- resolve(n,"serial"); cp <- resolve(n,"pernum")
  if (origin == "ACS") { ca <- resolve(n,"sample"); d$src_key <- paste(d[[ca]],d[[cy]],d[[cs]],d[[cp]],sep=":"); d$source_doiy <- suppressWarnings(as.integer(d[[cy]])); d$source_origin <- "ACS"; d$from_cps <- 0L } else { d$src_key <- paste("CPS",d[[cy]],d[[cs]],d[[cp]],sep=":"); d$source_doiy <- suppressWarnings(as.integer(d[[cy]])); d$source_origin <- "CPS"; d$from_cps <- 0L }
  ow <- resolve(n,"ownershp",FALSE); rmv <- resolve(n,"rooms",FALSE); bed <- resolve(n,"bedrooms",FALSE)
  d$ownershp_raw <- if (!is.na(ow)) d[[ow]] else NA_integer_; d$rooms_raw <- if (!is.na(rmv)) d[[rmv]] else NA_integer_; d$bedrooms_raw <- if (!is.na(bed)) d[[bed]] else NA_integer_
  req(!anyDuplicated(d$src_key), paste("duplicate source key in",origin), "source_adapter"); d
}
patch_once <- function(txt, pat, repl, tag) { hit <- gregexpr(pat,txt,fixed=TRUE)[[1L]]; nh <- if (length(hit)==1L && hit[1L]==-1L) 0L else length(hit); req(nh==1L, paste(tag,"matched",nh,"times"), "vendor_adapter"); sub(pat,repl,txt,fixed=TRUE) }
fn_env <- new.env(parent=globalenv()); fx <- parse(file=file.path(vendor_dir,"functions.R")); for (e in fx) if (is.call(e) && identical(e[[1L]],as.name("<-")) && is.call(e[[3L]]) && identical(e[[3L]][[1L]],as.name("function"))) eval(e,envir=fn_env)
req(exists("apply_labels",envir=fn_env,inherits=FALSE), "vendor apply_labels missing", "vendor_adapter"); apply_labels <- get("apply_labels",envir=fn_env)
min.age <- 25; max.age <- 45; t.min <- -5; t.max <- 10; ref <- "-2"; age.cutoff <- 44; seed.val <- 9746290
run_state <- function(st) {
  st_dir <- file.path(source_out,"partitions",sprintf("statefip_%02d",st)); work <- file.path(outdir,sprintf("statefip_%02d",st)); dir.create(work,recursive=TRUE,showWarnings=FALSE)
  rawdir <- file.path(work,"Data","Raw_Data"); cleandir <- file.path(work,"Data","Cleaned_Data"); dir.create(file.path(rawdir,"ACS"),recursive=TRUE); dir.create(file.path(rawdir,"CPS"),recursive=TRUE); dir.create(cleandir,recursive=TRUE)
  a <- add_lineage(load_data(file.path(st_dir,"acs_raw.RData")),"ACS"); c <- add_lineage(load_data(file.path(st_dir,"cps_raw.RData")),"CPS")
  data <- a; save(data,file=file.path(rawdir,"ACS","raw_acs.RData"),compress=FALSE); data <- c; save(data,file=file.path(rawdir,"CPS","raw_cps.RData"),compress=FALSE); rm(a,c,data); invisible(gc())
  clean_one <- function(src, kind) {
    txt <- paste(readLines(file.path(vendor_dir,src),warn=FALSE),collapse="\n")
    carry <- if (kind == "ACS") "subset(select=c(id, serial, pernum, doiy, src_key, ownershp_raw, rooms_raw, bedrooms_raw, source_origin, source_doiy, from_cps, wgt," else "subset(select=c(id, serial, pernum, doiy, src_key, source_origin, source_doiy, from_cps, wgt,"
    txt <- patch_once(txt,"subset(select=c(id, serial, pernum, doiy, wgt,", carry, paste(kind,"clean-select"))
    ee <- new.env(parent=globalenv()); ee$rawdir <- rawdir; ee$cleandir <- cleandir; ee$wrkdir <- work; ee$apply_labels <- apply_labels; list2env(list(min.age=min.age,max.age=max.age,t.min=t.min,t.max=t.max,ref=ref,age.cutoff=age.cutoff,seed.val=seed.val),ee)
    ex <- parse(text=txt); req(length(ex)>=3L,"vendor cleaner unexpectedly short","vendor_adapter")
    h1 <- paste(deparse(ex[[1L]]), collapse=" "); h2 <- paste(deparse(ex[[2L]]), collapse=" ")
    req(grepl("^rm\\(list = setdiff", h1), paste(kind,"cleaner head expression 1 changed"), "vendor_adapter")
    req(grepl("setup\\.R", h2), paste(kind,"cleaner head expression 2 changed"), "vendor_adapter")
    for (i in 3:length(ex)) tryCatch(eval(ex[[i]],envir=ee), error=function(e) fail(paste(kind,"cleaner expression",i,conditionMessage(e)),"clean"))
    ep <- file.path(cleandir, if (kind=="ACS") "acs_clean.RData" else "cps_clean.RData"); req(file.exists(ep), paste(kind,"cleaner did not write",ep), "clean")
  }
  clean_one("clean_acs.R","ACS"); clean_one("clean_cps.R","CPS")
  mt <- paste(readLines(file.path(vendor_dir,"matching.R"),warn=FALSE),collapse="\n"); mt <- patch_once(mt,"subset(select=c(age1b,match_bin,", "subset(select=c(age1b,match_bin,src_key,ownershp_raw,rooms_raw,bedrooms_raw,source_origin,source_doiy,from_cps,", "match-source-select"); mt <- patch_once(mt,"subset(select=-c(from_cps))","identity()","match-preserve-from-cps")
  mx <- parse(text=mt); pick <- function(nm) { z <- NULL; for (e in mx) if (is.call(e)&&identical(e[[1L]],as.name("<-"))&&identical(e[[2L]],as.name(nm))&&is.call(e[[3L]])&&identical(e[[3L]][[1L]],as.name("function"))) z <- e; z }; defs <- lapply(c("run_match","fn_match","fn_pseudo_panel"),pick); req(!any(vapply(defs,is.null,logical(1))),"matching API definitions missing","match")
  me <- new.env(parent=globalenv()); me$cleandir <- normalizePath(cleandir); me$wrkdir <- normalizePath(work); list2env(list(min.age=min.age,max.age=max.age,t.min=t.min,t.max=t.max,ref=ref,age.cutoff=age.cutoff,seed.val=seed.val),me); for (fn in ls(fn_env)) assign(fn,get(fn,envir=fn_env),envir=me); for (e in defs) eval(e,envir=me)
  set.seed(seed.val); panel <- tryCatch(get("fn_pseudo_panel",envir=me)(), error=function(e) fail(paste("matcher:",conditionMessage(e)),"match")); req(is.data.frame(panel) && nrow(panel)>0L,"matcher returned no rows","match")
  source_cols <- intersect(c("src_key","src_key.x","src_key.y","source_origin","source_origin.x","source_origin.y","from_cps","from_cps.x","from_cps.y","ownershp_raw"),names(panel)); req(any(grepl("src_key",source_cols)),"matched panel lost source key","lineage")
  saveRDS(panel,file.path(work,"cps_acs_pseudo-panel.rds"),compress=FALSE)
  hp <- file.path(st_dir,"housing_narrow.rds"); housing_status <- "PENDING_RAW_EXTRACT_STAGE"
  if (file.exists(hp)) { h <- data.table::as.data.table(readRDS(hp)); hk <- paste(h$SAMPLE,h$YEAR,h$SERIAL,h$PERNUM,sep=":"); pk <- if ("src_key.x" %in% names(panel)) panel$src_key.x else if ("src_key" %in% names(panel)) panel$src_key else panel$src_key.y; m <- match(pk,hk); panel$rooms_raw <- h$ROOMS_RAW[m]; panel$bedrooms_raw <- h$BEDROOMS_RAW[m]; panel$ownershp_raw <- h$OWNERSHP_RAW[m]; saveRDS(panel,file.path(work,"cps_acs_pseudo-panel_housing.rds"),compress=FALSE); housing_status <- "HOUSING_BRIDGE_COMPLETE" }
  rec <- list(status="STATE_MATCH_COMPLETE", statefip=st, panel_rows=nrow(panel), panel_columns=names(panel), source_columns=source_cols, housing_status=housing_status, generated=stamp()); jsonlite::write_json(rec,file.path(work,"state_receipt.json"),auto_unbox=TRUE,pretty=TRUE); rec
}
prepare_housing <- function() {
  hp <- Sys.getenv("HOUSING_RAW", file.path(root,"inputs","ACS","local_extract27_20260919","extract27.dta"))
  if (!file.exists(hp)) { logp("housing_raw_absent", hp); return(invisible(FALSE)) }
  req(requireNamespace("haven", quietly=TRUE), "haven required for housing raw stage", "housing")
  logp("housing_raw_read_start", hp)
  h <- haven::read_dta(hp, col_select=c(year,sample,serial,pernum,statefip,rooms,bedrooms,ownershp))
  names(h) <- toupper(names(h)); need <- c("YEAR","SAMPLE","SERIAL","PERNUM","STATEFIP","ROOMS","BEDROOMS","OWNERSHP")
  req(all(need %in% names(h)), paste("housing raw missing", paste(setdiff(need,names(h)),collapse=",")), "housing")
  h$ROOMS_RAW <- h$ROOMS; h$BEDROOMS_RAW <- h$BEDROOMS; h$OWNERSHP_RAW <- h$OWNERSHP
  hk <- paste(h$SAMPLE,h$YEAR,h$SERIAL,h$PERNUM,sep=":"); req(!anyDuplicated(hk), "national housing source key is non-unique", "housing")
  for (st in states) { z <- h[h$STATEFIP == st, c(need,"ROOMS_RAW","BEDROOMS_RAW","OWNERSHP_RAW"), drop=FALSE]; req(nrow(z)>0L, paste("housing raw has no rows for FIPS",st), "housing"); dir.create(file.path(source_out,"partitions",sprintf("statefip_%02d",st)),recursive=TRUE,showWarnings=FALSE); saveRDS(z,file.path(source_out,"partitions",sprintf("statefip_%02d",st),"housing_narrow.rds"),compress=FALSE) }
  jsonlite::write_json(list(status="HOUSING_RAW_PARTITION_COMPLETE", rows=nrow(h), key_unique=TRUE, fields=names(h), generated=stamp()), file.path(outdir,"housing_partition_receipt.json"), auto_unbox=TRUE, pretty=TRUE); rm(h); invisible(gc()); TRUE
}
logp("gates passed; matcher adapter staged")
if (phase %in% c("smoke","production")) prepare_housing()
if (phase == "preflight") { jsonlite::write_json(list(status="MATCHER_PREFLIGHT_PASS", states=states, source_out=source_out, vendor_sha256=as.list(vh), no_job_submitted=TRUE, generated=stamp()), file.path(outdir,"stage_receipt.json"), auto_unbox=TRUE, pretty=TRUE); message("MATCHER_PREFLIGHT_PASS ",outdir); quit(save="no",status=0L) }
results <- lapply(states, function(st) { logp("state_start",st); z <- run_state(st); logp("state_complete",st,z$panel_rows); z })
jsonlite::write_json(list(status="MATCH_HOUSING_STAGE_COMPLETE", phase=phase, states=states, results=results, national_fit="NOT_RUN", generated=stamp()), file.path(outdir,"stage_receipt.json"), auto_unbox=TRUE, pretty=TRUE)
message("MATCH_HOUSING_STAGE_COMPLETE ",outdir)
