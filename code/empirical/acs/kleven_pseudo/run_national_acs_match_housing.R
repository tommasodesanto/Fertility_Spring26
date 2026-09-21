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
packages <- c("Matching","plyr","data.table","dplyr","tidyr","forcats","purrr","labelled","ipumsr","haven","stringr","here","fixest")
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
  cm <- resolve(n,"month",FALSE); ca <- resolve(n,"sample",FALSE)
  csex <- resolve(n,"sex",FALSE); cage <- resolve(n,"age",FALSE)
  d$source_year <- suppressWarnings(as.integer(d[[cy]])); d$source_serial <- d[[cs]]; d$source_pernum <- d[[cp]]
  d$source_month <- if (!is.na(cm)) d[[cm]] else NA_integer_
  d$source_sample <- if (!is.na(ca)) d[[ca]] else NA_integer_
  d$source_sex <- if (!is.na(csex)) d[[csex]] else NA_integer_
  d$source_age <- if (!is.na(cage)) d[[cage]] else NA_integer_
  if (origin == "ACS") {
    req(!is.na(ca), "ACS source lacks SAMPLE", "source_adapter")
    d$src_key <- paste(d$source_sample,d$source_year,d$source_serial,d$source_pernum,sep=":")
    d$source_hhcluster <- paste(d$source_sample,d$source_year,d$source_serial,sep=":")
    d$source_doiy <- d$source_year; d$source_origin <- "ACS"; d$from_cps <- 0L
  } else {
    req(!is.na(cm), "CPS source lacks MONTH; annual key is not sufficient", "source_adapter")
    d$src_key <- paste("CPS",d$source_year,d$source_month,d$source_serial,d$source_pernum,sep=":")
    d$source_hhcluster <- paste("CPS",d$source_year,d$source_month,d$source_serial,sep=":")
    d$source_doiy <- d$source_year; d$source_origin <- "CPS"; d$from_cps <- 0L
  }
  ow <- resolve(n,"ownershp",FALSE); rmv <- resolve(n,"rooms",FALSE); bed <- resolve(n,"bedrooms",FALSE)
  d$source_ownershp_raw <- if (!is.na(ow)) d[[ow]] else NA_integer_; d$ownershp_raw <- NA_integer_
  d$rooms_raw <- if (!is.na(rmv)) d[[rmv]] else NA_integer_; d$bedrooms_raw <- if (!is.na(bed)) d[[bed]] else NA_integer_
  req(!anyDuplicated(d$src_key), paste("duplicate source key in",origin), "source_adapter"); d
}
patch_once <- function(txt, pat, repl, tag) { hit <- gregexpr(pat,txt,fixed=TRUE)[[1L]]; nh <- if (length(hit)==1L && hit[1L]==-1L) 0L else length(hit); req(nh==1L, paste(tag,"matched",nh,"times"), "vendor_adapter"); sub(pat,repl,txt,fixed=TRUE) }
coalesce_field <- function(d, stem) {
  candidates <- intersect(c(paste0(stem,".x"), stem, paste0(stem,".y")), names(d))
  req(length(candidates) > 0L, paste("lineage field absent:", stem), "lineage")
  out <- d[[candidates[[1L]]]]
  if (length(candidates) > 1L) for (cc in candidates[-1L]) out <- dplyr::coalesce(out, d[[cc]])
  out
}
normalize_lineage <- function(d) {
  for (stem in c("src_key","source_origin","source_doiy","source_year","source_month",
                 "source_sample","source_serial","source_pernum","source_sex","source_age",
                 "source_hhcluster","source_ownershp_raw","ownershp_raw","rooms_raw","bedrooms_raw"))
    d[[stem]] <- coalesce_field(d, stem)
  raw_from_cps <- if (any(c("from_cps.x","from_cps","from_cps.y") %in% names(d))) coalesce_field(d, "from_cps") else rep(0L, nrow(d))
  cps_flag <- as.character(d$source_origin) == "CPS" |
    suppressWarnings(as.integer(raw_from_cps)) == 1L |
    grepl("^CPS:", as.character(d$src_key))
  d$from_cps <- as.integer(cps_flag)
  d$source_origin <- ifelse(cps_flag, "CPS", as.character(d$source_origin))
  req(all(d$source_origin %in% c("ACS","CPS")), "unresolved source origin", "lineage")
  d$source_hhcluster <- ifelse(d$source_origin == "ACS",
                               paste(d$source_sample,d$source_year,d$source_serial,sep=":"),
                               paste("CPS",d$source_year,d$source_month,d$source_serial,sep=":"))
  d$source_hh_cluster <- d$source_hhcluster
  d
}
fn_env <- new.env(parent=globalenv()); fx <- parse(file=file.path(vendor_dir,"functions.R")); for (e in fx) if (is.call(e) && identical(e[[1L]],as.name("<-")) && is.call(e[[3L]]) && identical(e[[3L]][[1L]],as.name("function"))) eval(e,envir=fn_env)
req(exists("apply_labels",envir=fn_env,inherits=FALSE), "vendor apply_labels missing", "vendor_adapter"); apply_labels <- get("apply_labels",envir=fn_env)
estimator_file <- file.path(root,"code","empirical","acs","kleven_pseudo","estimate_national_first_birth_housing.R")
req(file.exists(estimator_file), paste("national estimator absent:", estimator_file), "estimator")
estimator_env <- new.env(parent=globalenv()); sys.source(estimator_file, envir=estimator_env)
req(exists("estimate_national_first_birth_housing", envir=estimator_env, inherits=FALSE), "national estimator function missing", "estimator")
min.age <- 25; max.age <- 45; t.min <- -5; t.max <- 10; ref <- "-2"; age.cutoff <- 44; seed.val <- 9746290
run_national_estimation <- function(panel, dest, label) {
  dir.create(dest, recursive=TRUE, showWarnings=FALSE)
  cp <- function(x) logp("estimator", label, jsonlite::toJSON(x, auto_unbox=TRUE, null="null"))
  fit <- tryCatch(estimator_env$estimate_national_first_birth_housing(panel, output_dir=dest, checkpoint=cp), error=function(e) fail(paste(label,"estimator:",conditionMessage(e)),"estimator"))
  req(identical(fit$status,"ESTIMATION_COMPLETE_DIAGNOSTIC"), paste(label,"estimator status invalid"), "estimator")
  req(length(fit$fits) > 0L && all(vapply(fit$fits, function(z) is.matrix(z$full_vcov) && nrow(z$full_vcov) > 0L, logical(1))), paste(label,"full covariance output missing"), "estimator")
  expected_outputs <- c("national_event_curves.csv","national_contrasts.csv","national_raw_baselines.csv","national_counts_event_ess.csv","national_fit_status.csv","national_housing_event_curves.png")
  req(all(file.exists(file.path(dest, expected_outputs))), paste(label,"estimator output incomplete"), "estimator")
  req(length(list.files(dest, pattern="^checkpoint_.*\\.rds$")) == length(fit$fits), paste(label,"fit checkpoints incomplete"), "estimator")
  receipt <- list(status=fit$status, label=label, output_dir=dest, fit_count=length(fit$fits), outcomes=unique(fit$curves$outcome), specifications=unique(fit$curves$specification), full_covariance=TRUE, contrast="+3_minus_-1", figure=file.path(dest,"national_housing_event_curves.png"), generated=stamp())
  jsonlite::write_json(receipt, file.path(dest,"estimation_receipt.json"), auto_unbox=TRUE, pretty=TRUE)
  rm(fit); invisible(receipt)
}
run_state <- function(st) {
  st_dir <- file.path(source_out,"partitions",sprintf("statefip_%02d",st)); work <- file.path(outdir,sprintf("statefip_%02d",st)); dir.create(work,recursive=TRUE,showWarnings=FALSE)
  rawdir <- file.path(work,"Data","Raw_Data"); cleandir <- file.path(work,"Data","Cleaned_Data"); dir.create(file.path(rawdir,"ACS"),recursive=TRUE); dir.create(file.path(rawdir,"CPS"),recursive=TRUE); dir.create(cleandir,recursive=TRUE)
  a <- add_lineage(load_data(file.path(st_dir,"acs_raw.RData")),"ACS"); c <- add_lineage(load_data(file.path(st_dir,"cps_raw.RData")),"CPS")
  data <- a; save(data,file=file.path(rawdir,"ACS","raw_acs.RData"),compress=FALSE); data <- c; save(data,file=file.path(rawdir,"CPS","raw_cps.RData"),compress=FALSE); rm(a,c,data); invisible(gc())
  clean_one <- function(src, kind) {
    txt <- paste(readLines(file.path(vendor_dir,src),warn=FALSE),collapse="\n")
    lineage_select <- "src_key, source_origin, source_doiy, source_year, source_month, source_sample, source_serial, source_pernum, source_sex, source_age, source_hhcluster, source_ownershp_raw, from_cps, ownershp_raw, rooms_raw, bedrooms_raw"
    if (kind == "ACS") {
      txt <- patch_once(txt,"subset(select=c(id, serial, pernum, doiy, wgt,",
                         paste0("subset(select=c(id, serial, pernum, doiy, ", lineage_select, ", wgt,"),
                         paste(kind,"clean-select"))
    } else {
      txt <- patch_once(txt,"subset(select=c(doiy, month, serial, pernum, wgt, id,",
                         paste0("subset(select=c(doiy, month, serial, pernum, wgt, id, ", lineage_select, ","),
                         paste(kind,"clean-select"))
    }
    ee <- new.env(parent=globalenv()); ee$rawdir <- rawdir; ee$cleandir <- cleandir; ee$wrkdir <- work; ee$apply_labels <- apply_labels; list2env(list(min.age=min.age,max.age=max.age,t.min=t.min,t.max=t.max,ref=ref,age.cutoff=age.cutoff,seed.val=seed.val),ee)
    ex <- parse(text=txt); req(length(ex)>=3L,"vendor cleaner unexpectedly short","vendor_adapter")
    h1 <- paste(deparse(ex[[1L]]), collapse=" "); h2 <- paste(deparse(ex[[2L]]), collapse=" ")
    req(grepl("^rm\\(list = setdiff", h1), paste(kind,"cleaner head expression 1 changed"), "vendor_adapter")
    req(grepl("setup\\.R", h2), paste(kind,"cleaner head expression 2 changed"), "vendor_adapter")
    for (i in 3:length(ex)) tryCatch(eval(ex[[i]],envir=ee), error=function(e) fail(paste(kind,"cleaner expression",i,conditionMessage(e)),"clean"))
    ep <- file.path(cleandir, if (kind=="ACS") "acs_clean.RData" else "cps_clean.RData"); req(file.exists(ep), paste(kind,"cleaner did not write",ep), "clean")
  }
  clean_one("clean_acs.R","ACS"); clean_one("clean_cps.R","CPS")
  vanilla_check <- function(kind) {
    vroot <- file.path(work, "vanilla"); vraw <- file.path(vroot, "Data", "Raw_Data"); vclean <- file.path(vroot, "Data", "Cleaned_Data")
    dir.create(file.path(vraw, kind), recursive=TRUE, showWarnings=FALSE); dir.create(vclean, recursive=TRUE, showWarnings=FALSE)
    src <- file.path(rawdir, kind, if (kind == "ACS") "raw_acs.RData" else "raw_cps.RData"); file.copy(src, file.path(vraw, kind, basename(src)), overwrite=TRUE)
    txt <- paste(readLines(file.path(vendor_dir, if (kind == "ACS") "clean_acs.R" else "clean_cps.R"), warn=FALSE), collapse="\n"); ex <- parse(text=txt); h1 <- paste(deparse(ex[[1L]]),collapse=" "); h2 <- paste(deparse(ex[[2L]]),collapse=" "); req(grepl("^rm\\(list = setdiff",h1) && grepl("setup\\.R",h2), paste(kind,"vanilla head changed"), "vanilla")
    ee <- new.env(parent=globalenv()); ee$rawdir <- file.path(vroot,"Data","Raw_Data"); ee$cleandir <- vclean; ee$wrkdir <- vroot; ee$apply_labels <- apply_labels; list2env(list(min.age=min.age,max.age=max.age,t.min=t.min,t.max=t.max,ref=ref,age.cutoff=age.cutoff,seed.val=seed.val),ee)
    for (i in 3:length(ex)) tryCatch(eval(ex[[i]],envir=ee), error=function(e) fail(paste(kind,"vanilla cleaner",conditionMessage(e)),"vanilla"))
    ve <- new.env(); load(file.path(vclean, if (kind=="ACS") "acs_clean.RData" else "cps_clean.RData"), envir=ve); pe <- new.env(); load(file.path(cleandir, if (kind=="ACS") "acs_clean.RData" else "cps_clean.RData"), envir=pe); va <- if(kind=="ACS") ve$acs else ve$cps; pa <- if(kind=="ACS") pe$acs else pe$cps
    metadata <- c("src_key","source_origin","source_doiy","source_year","source_month","source_sample","source_serial","source_pernum","source_sex","source_age","source_hhcluster","source_hh_cluster","source_ownershp_raw","from_cps","ownershp_raw","rooms_raw","bedrooms_raw")
    va_core <- setdiff(names(va), metadata); pa_core <- setdiff(names(pa), metadata)
    req(nrow(va)==nrow(pa) && identical(va_core, pa_core) && length(va_core)>0L, paste(kind,"vanilla/adapter clean schema differs"),"vanilla")
    for (cc in va_core) {
      req(identical(attributes(va[[cc]]), attributes(pa[[cc]])), paste(kind,"vanilla/adapter clean attributes differ:",cc), "vanilla")
      req(isTRUE(all.equal(va[[cc]],pa[[cc]],check.attributes=TRUE)), paste(kind,"vanilla/adapter clean column differs:",cc),"vanilla")
    }
    lineage <- c("src_key","source_origin","source_doiy","source_year","source_month","source_sample","source_serial","source_pernum","source_sex","source_age","source_hhcluster","source_ownershp_raw","from_cps")
    req(all(lineage %in% names(pa)), paste(kind,"adapter lost lineage columns:",paste(setdiff(lineage,names(pa)),collapse=",")),"lineage")
    raw_env <- new.env(parent=emptyenv()); load(src, envir=raw_env); raw <- raw_env$data
    ryear <- resolve(names(raw),"year"); rserial <- resolve(names(raw),"serial"); rpernum <- resolve(names(raw),"pernum")
    rmonth <- resolve(names(raw),"month",FALSE); rsample <- resolve(names(raw),"sample",FALSE)
    raw_key <- if (kind == "ACS") paste(raw[[rsample]],raw[[ryear]],raw[[rserial]],raw[[rpernum]],sep=":") else paste("CPS",raw[[ryear]],raw[[rmonth]],raw[[rserial]],raw[[rpernum]],sep=":")
    clean_key <- if (kind == "ACS") paste(pa$source_sample,pa$source_year,pa$source_serial,pa$source_pernum,sep=":") else paste("CPS",pa$source_year,pa$source_month,pa$source_serial,pa$source_pernum,sep=":")
    cid <- match(clean_key, raw_key)
    req(!anyNA(cid), paste(kind,"cleaned source key not found in raw source"), "lineage")
    for (z in c("year","serial","pernum","sex","age")) {
      rr <- resolve(names(raw),z); actual <- pa[[paste0("source_",ifelse(z=="year","year",z))]]
      req(isTRUE(all.equal(as.character(actual),as.character(raw[[rr]][cid]))), paste(kind,"source",z,"changed after clean"), "lineage")
    }
    if (kind == "CPS") {
      rr <- resolve(names(raw),"month"); req(isTRUE(all.equal(as.character(pa$source_month),as.character(raw[[rr]][cid]))), "CPS source month changed after clean", "lineage")
    }
    if (kind == "ACS") {
      rr <- resolve(names(raw),"ownershp",FALSE); if (!is.na(rr)) req(isTRUE(all.equal(as.character(pa$source_ownershp_raw),as.character(raw[[rr]][cid]))), "ACS source OWNERSHP changed after clean", "lineage")
    }
    # clean_cps.R:591 defines wgt=ifelse(is.na(asecwt),wtfinl,asecwt).  The
    # all-original-column vanilla equality above is the authoritative weight
    # guard; comparing CPS wgt only to raw ASECWT incorrectly rejects the
    # documented WTFinL fallback.
    list(status="PASS", kind=kind, rows=nrow(pa), protected_columns=va_core, metadata_columns=metadata, lineage_columns=lineage)
  }
  vanilla_receipts <- list(ACS=vanilla_check("ACS"), CPS=vanilla_check("CPS"))
  mt_vanilla <- paste(readLines(file.path(vendor_dir,"matching.R"),warn=FALSE),collapse="\n")
  mt <- mt_vanilla
  mt <- patch_once(mt,"subset(select=c(age1b,match_bin,", "subset(select=c(age1b,match_bin,src_key,source_hhcluster,source_year,source_month,source_sample,source_serial,source_pernum,source_sex,source_age,source_ownershp_raw,ownershp_raw,rooms_raw,bedrooms_raw,source_origin,source_doiy,from_cps,", "match-source-select")
  mt <- patch_once(mt,"subset(select=-c(from_cps))","identity()","match-preserve-from-cps")
  run_panel <- function(match_text, clean_root, tag) {
    mx <- parse(text=match_text); pick <- function(nm) { z <- NULL; for (e in mx) if (is.call(e)&&identical(e[[1L]],as.name("<-"))&&identical(e[[2L]],as.name(nm))&&is.call(e[[3L]])&&identical(e[[3L]][[1L]],as.name("function"))) z <- e; z }; defs <- lapply(c("run_match","fn_match","fn_pseudo_panel"),pick); req(!any(vapply(defs,is.null,logical(1))),paste(tag,"matching API definitions missing"),"match")
    me <- new.env(parent=globalenv()); me$cleandir <- normalizePath(clean_root); me$wrkdir <- normalizePath(work); list2env(list(min.age=min.age,max.age=max.age,t.min=t.min,t.max=t.max,ref=ref,age.cutoff=age.cutoff,seed.val=seed.val),me); for (fn in ls(fn_env)) assign(fn,get(fn,envir=fn_env),envir=me); for (e in defs) eval(e,envir=me)
    set.seed(seed.val); z <- tryCatch(get("fn_pseudo_panel",envir=me)(), error=function(e) fail(paste(tag,"matcher:",conditionMessage(e)),"match")); req(is.data.frame(z) && nrow(z)>0L,paste(tag,"matcher returned no rows"),"match"); z
  }
  panel <- run_panel(mt, cleandir, "adapter")
  panel <- normalize_lineage(panel)
  match_receipt <- list(status="NOT_RUN", reason="production does not duplicate the author matcher")
  if (phase == "smoke") {
    vanilla_panel <- run_panel(mt_vanilla, file.path(work,"vanilla","Data","Cleaned_Data"), "vanilla")
    match_metadata <- grep("^(src_key|source_.*|from_cps|ownershp_raw|rooms_raw|bedrooms_raw)(\\.|$)", names(panel), value=TRUE)
    vanilla_core <- setdiff(names(vanilla_panel), match_metadata); adapter_core <- setdiff(names(panel), match_metadata)
    req(nrow(vanilla_panel)==nrow(panel) && identical(vanilla_core, adapter_core) && length(vanilla_core)>0L, "vanilla/adapter matching schema differs", "vanilla")
    for (cc in vanilla_core) {
      req(identical(attributes(vanilla_panel[[cc]]), attributes(panel[[cc]])), paste("vanilla/adapter matching attributes differ:",cc), "vanilla")
      req(isTRUE(all.equal(vanilla_panel[[cc]],panel[[cc]],check.attributes=TRUE)), paste("vanilla/adapter matching column differs:",cc), "vanilla")
    }
    match_receipt <- list(status="PASS", rows=nrow(panel), protected_columns=vanilla_core, metadata_columns=match_metadata)
  }
  source_cols <- intersect(c("src_key","source_origin","source_doiy","source_year","source_month","source_sample","source_serial","source_pernum","source_sex","source_age","source_hhcluster","source_hh_cluster","source_ownershp_raw","from_cps","ownershp_raw"),names(panel)); req(any(grepl("src_key",source_cols)),"matched panel lost source key","lineage")
  saveRDS(panel,file.path(work,"cps_acs_pseudo-panel.rds"),compress=FALSE)
  hp <- file.path(st_dir,"housing_narrow.rds"); housing_status <- "PENDING_RAW_EXTRACT_STAGE"
  req(file.exists(hp), paste("housing packet absent for FIPS", st), "housing")
  bridge_receipt <- list(status="NOT_RUN")
  if (file.exists(hp)) {
    h <- data.table::as.data.table(readRDS(hp)); hk <- paste(h$SAMPLE,h$YEAR,h$SERIAL,h$PERNUM,sep=":"); req(!anyDuplicated(hk), "housing partition key duplicated", "housing")
    true_acs <- panel$source_origin == "ACS" & panel$from_cps == 0L
    req(any(true_acs), paste("no true ACS rows for FIPS",st), "housing")
    one_year_sample <- is.finite(suppressWarnings(as.numeric(panel$source_sample))) &
      (suppressWarnings(as.numeric(panel$source_sample)) %% 100 == 1)
    verified_overlap <- true_acs & panel$source_year %in% 2005:2019 & one_year_sample
    out_of_overlap <- true_acs & !verified_overlap
    key_observed <- is.finite(suppressWarnings(as.numeric(panel$source_year))) &
      is.finite(suppressWarnings(as.numeric(panel$source_sample))) &
      is.finite(suppressWarnings(as.numeric(panel$source_serial))) &
      is.finite(suppressWarnings(as.numeric(panel$source_pernum)))
    overlap_missing_key <- sum(verified_overlap & !key_observed)
    req(overlap_missing_key == 0L, paste("verified-overlap rows with missing source key:", overlap_missing_key), "housing")
    pk <- paste(panel$source_sample,panel$source_year,panel$source_serial,panel$source_pernum,sep=":"); m <- match(pk,hk)
    verified_unmatched <- sum(verified_overlap & is.na(m)); req(verified_unmatched == 0L, paste("verified-overlap housing keys unmatched:", verified_unmatched), "housing")
    observed_match <- verified_overlap & !is.na(m)
    concordance <- list()
    for (pair in list(c("source_sex","SEX"),c("source_age","AGE"),c("ownershp_raw","OWNERSHP_RAW"))) {
      lhs <- suppressWarnings(as.numeric(panel[[if (pair[[1L]] == "ownershp_raw") "source_ownershp_raw" else pair[[1L]]]][observed_match])); rhs <- suppressWarnings(as.numeric(h[[pair[[2L]]]][m[observed_match]])); lhs_missing <- is.na(lhs); rhs_missing <- is.na(rhs); equal <- (lhs_missing & rhs_missing) | (!lhs_missing & !rhs_missing & lhs == rhs); mismatch <- sum(!equal); req(mismatch == 0L, paste("housing",pair[[2L]],"concordance failed"), "housing"); concordance[[pair[[2L]]]] <- list(compared=length(lhs), observed_shared=sum(!lhs_missing & !rhs_missing), lhs_missing=sum(lhs_missing), rhs_missing=sum(rhs_missing), both_missing=sum(lhs_missing & rhs_missing), lhs_only_missing=sum(lhs_missing & !rhs_missing), rhs_only_missing=sum(!lhs_missing & rhs_missing), mismatches=mismatch)
    }
    panel$rooms_raw <- NA_real_; panel$bedrooms_raw <- NA_real_; panel$ownershp_raw <- NA_real_; panel$rooms_raw[observed_match] <- h$ROOMS_RAW[m[observed_match]]; panel$bedrooms_raw[observed_match] <- h$BEDROOMS_RAW[m[observed_match]]; panel$ownershp_raw[observed_match] <- h$OWNERSHP_RAW[m[observed_match]]
    saveRDS(panel,file.path(work,"cps_acs_pseudo-panel_housing.rds"),compress=FALSE); housing_status <- "HOUSING_BRIDGE_COMPLETE"; bridge_receipt <- list(status=housing_status, true_acs_rows=sum(true_acs), verified_overlap_rows=sum(verified_overlap), verified_overlap_eligible=sum(observed_match), verified_overlap_unmatched=verified_unmatched, overlap_missing_key=overlap_missing_key, out_of_verified_overlap=sum(out_of_overlap), non_true_acs_rows=sum(!true_acs), key_unique=TRUE, concordance=concordance, source_hhcluster=TRUE)
    if (phase == "smoke") run_national_estimation(panel, file.path(work,"national_first_birth_housing"), paste0("statefip_",st))
  }
  rec <- list(status="STATE_MATCH_COMPLETE", statefip=st, panel_rows=nrow(panel), panel_columns=names(panel), source_columns=source_cols, vanilla_adapter=vanilla_receipts, matching=match_receipt, housing_status=housing_status, housing_bridge=bridge_receipt, generated=stamp()); jsonlite::write_json(rec,file.path(work,"state_receipt.json"),auto_unbox=TRUE,pretty=TRUE); rec
}
prepare_housing <- function() {
  hp <- Sys.getenv("HOUSING_RAW", file.path(root,"inputs","ACS","local_extract27_20260919","extract27.dta"))
  if (!file.exists(hp)) { logp("housing_raw_absent", hp); return(invisible(FALSE)) }
  roster_select <- c("year","sample","serial","cbserial","hhwt","cluster","statefip","puma","strata","gq",
                     "ownershp","ownershpd","rooms","bedrooms","pernum","perwt","momloc","poploc","sploc",
                     "nchild","nchlt5","eldch","yngch","relate","sex","age","marst","fertyr","race","educ")
  expected_packet_fields <- c(toupper(roster_select),"ROOMS_RAW","BEDROOMS_RAW","OWNERSHP_RAW")
  packet_states <- valid_states
  packet_paths <- file.path(source_out,"partitions",sprintf("statefip_%02d",packet_states),"housing_narrow.rds")
  manifest_file <- file.path(source_out,"housing_packet_manifest.json")
  if (all(file.exists(packet_paths))) {
    req(file.exists(manifest_file), "housing packet manifest absent; refusing reuse", "housing")
    man <- jsonlite::fromJSON(manifest_file, simplifyVector=FALSE)
    req(identical(man$status,"HOUSING_PACKET_MANIFEST_COMPLETE") && identical(as.integer(unlist(man$states)), packet_states), "housing packet manifest does not cover all states", "housing")
    req(identical(sort(unlist(man$fields)), sort(expected_packet_fields)), "housing packet manifest fields differ", "housing")
    logp("housing_packet_reuse", length(packet_paths)); return(invisible(TRUE))
  }
  req(!any(file.exists(packet_paths)), "partial housing packet exists; refusing mixed reuse", "housing")
  req(requireNamespace("haven", quietly=TRUE), "haven required for housing raw stage", "housing")
  logp("housing_raw_read_start", hp)
  h <- haven::read_dta(hp, col_select=roster_select)
  names(h) <- toupper(names(h)); need <- toupper(roster_select)
  req(all(need %in% names(h)), paste("housing raw missing", paste(setdiff(need,names(h)),collapse=",")), "housing")
  h$ROOMS_RAW <- h$ROOMS; h$BEDROOMS_RAW <- h$BEDROOMS; h$OWNERSHP_RAW <- h$OWNERSHP
  hk <- paste(h$SAMPLE,h$YEAR,h$SERIAL,h$PERNUM,sep=":"); req(!anyDuplicated(hk), "national housing source key is non-unique", "housing")
  packet_rows <- numeric(length(packet_states))
  for (i in seq_along(packet_states)) { st <- packet_states[[i]]; z <- h[h$STATEFIP == st, c(need,"ROOMS_RAW","BEDROOMS_RAW","OWNERSHP_RAW"), drop=FALSE]; req(nrow(z)>0L, paste("housing raw has no rows for FIPS",st), "housing"); packet_rows[[i]] <- nrow(z); dir.create(file.path(source_out,"partitions",sprintf("statefip_%02d",st)),recursive=TRUE,showWarnings=FALSE); saveRDS(z,packet_paths[[i]],compress=FALSE) }
  manifest_rows <- lapply(seq_along(packet_states), function(i) { p <- packet_paths[[i]]; list(statefip=packet_states[[i]], path=p, rows=as.numeric(packet_rows[[i]]), bytes=as.numeric(file.info(p)$size)) })
  jsonlite::write_json(list(status="HOUSING_PACKET_MANIFEST_COMPLETE", states=packet_states, fields=names(h)[names(h) %in% expected_packet_fields], rows=as.numeric(nrow(h)), key_unique=TRUE, packets=manifest_rows, generated=stamp()), manifest_file, auto_unbox=TRUE, pretty=TRUE)
  jsonlite::write_json(list(status="HOUSING_RAW_PARTITION_COMPLETE", rows=nrow(h), key_unique=TRUE, fields=names(h), roster_fields=need, manifest=manifest_file, generated=stamp()), file.path(outdir,"housing_partition_receipt.json"), auto_unbox=TRUE, pretty=TRUE); rm(h); invisible(gc()); TRUE
}
logp("gates passed; matcher adapter staged")
if (phase %in% c("smoke","production")) prepare_housing()
if (phase == "preflight") { jsonlite::write_json(list(status="MATCHER_PREFLIGHT_PASS", states=states, source_out=source_out, vendor_sha256=as.list(vh), no_job_submitted=TRUE, generated=stamp()), file.path(outdir,"stage_receipt.json"), auto_unbox=TRUE, pretty=TRUE); message("MATCHER_PREFLIGHT_PASS ",outdir); quit(save="no",status=0L) }
results <- lapply(states, function(st) { logp("state_start",st); z <- run_state(st); logp("state_complete",st,z$panel_rows); z })
national_fit_status <- "NOT_RUN"
if (phase == "production") {
  panel_paths <- file.path(outdir, sprintf("statefip_%02d", states), "cps_acs_pseudo-panel_housing.rds")
  req(all(file.exists(panel_paths)), "production state panels incomplete", "pool")
  logp("pool_start", length(panel_paths))
  panels <- lapply(seq_along(panel_paths), function(i) { z <- readRDS(panel_paths[[i]]); logp("pool_state_loaded", states[[i]], nrow(z)); z })
  pooled <- dplyr::bind_rows(panels)
  pooled_file <- file.path(outdir, "national_cps_acs_pseudo-panel_housing.rds")
  saveRDS(pooled, pooled_file, compress=FALSE)
  logp("pool_complete", nrow(pooled), pooled_file)
  run_national_estimation(pooled, file.path(outdir,"national_first_birth_housing"), "national_pooled_women")
  national_fit_status <- "NATIONAL_POOLED_WOMEN_COMPLETE"
  rm(panels, pooled); invisible(gc())
}
if (phase == "smoke") national_fit_status <- "VT_POOLED_WOMEN_COMPLETE"
jsonlite::write_json(list(status="MATCH_HOUSING_STAGE_COMPLETE", phase=phase, states=states, results=results, national_fit=national_fit_status, generated=stamp()), file.path(outdir,"stage_receipt.json"), auto_unbox=TRUE, pretty=TRUE)
message("MATCH_HOUSING_STAGE_COMPLETE ",outdir)
