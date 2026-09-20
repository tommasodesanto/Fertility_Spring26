#!/usr/bin/env Rscript
# Driver for the approved first-birth housing estimator. It consumes only the
# saved v5 panel and verified NE source packet; it never reloads national ACS.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
root <- Sys.getenv("KLEVEN_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
panel_file <- Sys.getenv("PANEL_FILE", file.path(root, "overnight_ne_benchmark/ne_housing_v5/cps_acs_pseudo-panel.rds"))
packet_file <- Sys.getenv("PACKET_FILE", file.path(root, "output/kleven_acs_pilot/source_audit_extract27_20260919/ne_extract27_housing_key_packet.rds"))
audit_file <- Sys.getenv("AUDIT_MANIFEST_FILE", file.path(root, "output/kleven_acs_pilot/source_audit_extract27_20260919/source_audit_manifest.json"))
estimator_file <- Sys.getenv("ESTIMATOR_FILE", file.path(dirname(normalizePath(commandArgs()[1])), "estimate_first_birth_housing.R"))
outdir <- Sys.getenv("OUTDIR", file.path(root, "output/kleven_acs_pilot/first_birth_housing_20260920"))
run_estimation <- identical(Sys.getenv("RUN_ESTIMATION", "1"), "1")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
progress_file <- file.path(outdir, "progress.log")
checkpoint <- function(stage, detail="") { z <- sprintf("%s\t%s\t%s", format(Sys.time(), "%FT%T%z"), stage, detail); write(z, progress_file, append=TRUE); cat(z,"\n") }
options(error=function() { try(write_json(list(status="FAILED",error=geterrmessage(),generated=as.character(Sys.time())), file.path(outdir,"failure_receipt.json"), auto_unbox=TRUE, pretty=TRUE),silent=TRUE); q(save="no",status=1,runLast=FALSE) })
require_cols <- function(x, cols, label) { m <- setdiff(cols,names(x)); if(length(m)) stop(label," missing: ",paste(m,collapse=","),call.=FALSE) }
canonical_panel <- function(x) {
  x <- as.data.table(x); w <- c("src_key","source_origin","from_cps","emp_lw","wgt","t_es_lw","cohort","census","statefip","statename","gender","age_factor","doiy_factor")
  i <- match(toupper(w),toupper(names(x))); if(anyNA(i)) stop("v5 panel lacks: ",paste(w[is.na(i)],collapse=","),call.=FALSE)
  setnames(x,names(x)[i],w); x
}
canonical_source <- function(x) {
  x <- as.data.table(x); names(x) <- toupper(names(x)); k <- c("YEAR","SAMPLE","SERIAL","PERNUM")
  require_cols(x,c(k,"ROOMS_RAW","BEDROOMS_RAW","OWNERSHP_RAW"),"source packet")
  x <- x[,c(k,"ROOMS_RAW","BEDROOMS_RAW","OWNERSHP_RAW"),with=FALSE]
  setnames(x,c("ROOMS_RAW","BEDROOMS_RAW","OWNERSHP_RAW"),c("ROOMS","BEDROOMS","OWNERSHP")); x
}
coding <- list(
  rooms_valid=function(x,year)!is.na(x)&x%in%c(1:27,30), rooms_transform=function(x,year)x,
  rooms_cap=9, rooms_missing_codes=0, rooms_unknown_codes=28,
  bedrooms_valid=function(x,year)!is.na(x)&x%in%c(1:6,22), bedrooms_transform=function(x,year)ifelse(x==22,21,x-1),
  bedrooms_cap=5, bedrooms_missing_codes=0, bedrooms_unknown_codes=7,
  ownership_valid=function(x,year)!is.na(x)&x%in%c(1,2), ownership_missing_codes=0, ownership_unknown_codes=c(3,9),
  allow_uncapped_sensitivity=FALSE)
make_manifest <- function(x) {
  if(!identical(as.character(x$input_sha256),"edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e")) stop("verified source SHA mismatch",call.=FALSE)
  if(as.numeric(x$input_size)!=9919999546) stop("verified source byte mismatch",call.=FALSE)
  list(status="PASS",verified=TRUE,key_columns=c("YEAR","SAMPLE","SERIAL","PERNUM"),source_key_unique=TRUE,overlap_verified=TRUE,source_packet=packet_file,source_sha256=x$input_sha256,source_bytes=as.numeric(x$input_size),unique_key_intersection=as.numeric(x$unique_key_intersection),source_scope="NE states 9,23,25,33,44,50; raw housing codes retained")
}
checkpoint("startup",paste0("run_estimation=",run_estimation))
for(f in c(panel_file,packet_file,audit_file,estimator_file)) if(!file.exists(f)) stop("required input missing: ",f,call.=FALSE)
source(estimator_file,local=TRUE); if(!all(vapply(c("join_first_birth_housing","code_first_birth_housing","estimate_first_birth_housing"),exists,logical(1)))) stop("incomplete estimator interface",call.=FALSE)
checkpoint("estimator_loaded")
manifest <- make_manifest(fromJSON(audit_file)); source_housing <- canonical_source(readRDS(packet_file))
key <- do.call(paste,c(source_housing[,.(YEAR,SAMPLE,SERIAL,PERNUM)],sep=":")); if(anyDuplicated(key)) stop("source packet key duplicated",call.=FALSE)
checkpoint("source_packet_loaded",sprintf("rows=%s cols=%s",nrow(source_housing),ncol(source_housing)))
panel <- canonical_panel(readRDS(panel_file)); require_cols(panel,c("YEAR","SAMPLE","SERIAL","PERNUM","src_key","source_origin","from_cps","emp_lw","wgt","t_es_lw","cohort","census","statefip","statename","gender","age_factor","doiy_factor"),"v5 panel")
if(!is.numeric(panel$wgt)&&!is.integer(panel$wgt)) stop("v5 weights are not numeric",call.=FALSE)
checkpoint("v5_panel_loaded",sprintf("rows=%s cols=%s",nrow(panel),ncol(panel)))
# Real-data interface smoke: exercise the exact join and coding path on up to 2,000 rows.
smoke <- panel[seq_len(min(2000L,nrow(panel)))]; sk <- do.call(paste,c(smoke[,.(YEAR,SAMPLE,SERIAL,PERNUM)],sep=":")); ss <- source_housing[match(sk,key),]; ss <- ss[!is.na(YEAR)]
if(nrow(ss)) { j <- join_first_birth_housing(smoke,ss,manifest); invisible(code_first_birth_housing(j,coding)); if(!identical(j$first_birth_row_id,seq_len(nrow(smoke)))) stop("smoke changed row order",call.=FALSE); checkpoint("interface_smoke_pass",sprintf("panel_rows=%s source_rows=%s",nrow(smoke),nrow(ss))) } else checkpoint("interface_smoke_skipped","no source keys in first panel rows")
ready <- list(status="READY_FOR_LEAD_REVIEW",run_estimation=run_estimation,panel_file=panel_file,packet_file=packet_file,estimator_file=estimator_file,audit_manifest=manifest,panel_rows=nrow(panel),source_packet_rows=nrow(source_housing),coding_config=list(rooms_cap=9,bedrooms_cap=5,rooms_unknown_codes=28,bedrooms_unknown_codes=7,allow_uncapped_sensitivity=FALSE),analysis_scope="NE level housing outcomes with source-household clustering; diagnostic interpretation only",generated=as.character(Sys.time()))
write_json(ready,file.path(outdir,"readiness_receipt.json"),auto_unbox=TRUE,pretty=TRUE)
write.csv(data.table(variable=c("ROOMS","BEDROOMS","OWNERSHP"),cap=c(9,5,NA),unknown_codes=c("28","7","3,9")),file.path(outdir,"coding_config_receipt.csv"),row.names=FALSE)
checkpoint("readiness_receipt_written")
if(run_estimation) {
  checkpoint("estimation_start")
  result <- estimate_first_birth_housing(panel,source_housing,manifest,coding,checkpoint=function(x) { checkpoint("fit_complete",x$name); write_json(list(stage="fit_complete",detail=x$name,generated=as.character(Sys.time())),file.path(outdir,"latest_fit_checkpoint.json"),auto_unbox=TRUE,pretty=TRUE) })
  write.csv(result$join_audit,file.path(outdir,"join_audit.csv"),row.names=FALSE); write.csv(result$code_audit,file.path(outdir,"housing_code_audit.csv"),row.names=FALSE); write.csv(result$support,file.path(outdir,"support.csv"),row.names=FALSE); write.csv(result$curves,file.path(outdir,"curves.csv"),row.names=FALSE); write.csv(result$summary,file.path(outdir,"contrasts.csv"),row.names=FALSE)
  write_json(list(status=result$status,metadata=result$metadata,join_audit=result$join_audit,code_audit=result$code_audit,generated=as.character(Sys.time())),file.path(outdir,"result_receipt.json"),auto_unbox=TRUE,pretty=TRUE); checkpoint("estimation_complete",result$status)
} else checkpoint("complete","readiness_only")
