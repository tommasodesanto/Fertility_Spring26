#!/usr/bin/env Rscript
# Exact-cell donor availability counts only. No matching, assignment, weights,
# or coarsening are performed here.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("KLEVEN_ROOT", "/scratch/td2248/projects/kleven_acs_pilot_20260917")
diag <- Sys.getenv("SECOND_BIRTH_DIAG_DIR", file.path(root,"output/kleven_acs_pilot/second_birth_proxy_diagnostic_20260920"))
outdir <- Sys.getenv("OUTDIR", file.path(root,"output/kleven_acs_pilot/second_birth_candidate_support_20260920")); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
progress <- file.path(outdir,"progress.log"); ck <- function(s,d=""){z<-sprintf("%s\t%s\t%s",format(Sys.time(),"%FT%T%z"),s,d);write(z,progress,append=TRUE);cat(z,"\n")}
options(error=function(){writeLines(geterrmessage(),file.path(outdir,"failure.txt"));q(save="no",status=1,runLast=FALSE)})
tf <- file.path(diag,"donor_targets.rds"); df <- file.path(diag,"one_child_donors.rds"); if(!file.exists(tf)||!file.exists(df))stop("existing builder outputs missing",call.=FALSE)
t <- as.data.table(readRDS(tf)); d <- as.data.table(readRDS(df)); covars <- c("SEX","EDUC","MARST","RACE","STATEFIP")
need_t <- c("target_year","target_mother_age","target_child_age",paste0("target_",covars),"target_event_time","birth_gap","gap_full_pre","gap_reference"); need_d <- c("YEAR","AGE","sole_child_age",covars)
if(length(setdiff(need_t,names(t)))||length(setdiff(need_d,names(d))))stop("builder interface lacks exact support fields",call.=FALSE)
ck("inputs_loaded",sprintf("targets=%s donors=%s",nrow(t),nrow(d)))
tkeys <- c("target_year","target_mother_age","target_child_age",paste0("target_",covars)); dkeys <- c("YEAR","AGE","sole_child_age",covars)
t[,cell:=do.call(paste,c(.SD,sep=":")),.SDcols=tkeys]; d[,cell:=do.call(paste,c(.SD,sep=":")),.SDcols=dkeys]
dc <- d[,.(donor_rows=.N,eligible_donor_rows=if("donor_match_eligible"%in%names(d))sum(donor_match_eligible%in%TRUE)else .N),by=cell]
tc <- t[,.(target_rows=.N,full_pre_target_rows=sum(gap_full_pre%in%TRUE),reference_target_rows=sum(gap_reference%in%TRUE)),by=cell]
cc <- merge(tc,dc,by="cell",all.x=TRUE); cc[is.na(donor_rows),`:=`(donor_rows=0L,eligible_donor_rows=0L)]; cc[,`:=`(candidate_cell=donor_rows>0,eligible_candidate_cell=eligible_donor_rows>0)]
t <- merge(t,dc,by="cell",all.x=TRUE,sort=FALSE); t[is.na(donor_rows),`:=`(donor_rows=0L,eligible_donor_rows=0L)]; t[,`:=`(candidate_cell=donor_rows>0,eligible_candidate_cell=eligible_donor_rows>0)]
overall <- data.table(metric=c("target_rows","target_exact_cells","donor_rows","donor_exact_cells","target_rows_with_any_donor","target_rows_with_any_eligible_donor","target_cells_with_any_donor","target_cells_with_any_eligible_donor","full_pre_target_rows","full_pre_target_rows_with_any_donor","reference_target_rows","reference_target_rows_with_any_donor"),value=c(nrow(t),uniqueN(t$cell),nrow(d),uniqueN(d$cell),sum(t$candidate_cell),sum(t$eligible_candidate_cell),sum(cc$candidate_cell),sum(cc$eligible_candidate_cell),sum(t$gap_full_pre%in%TRUE),sum(t$gap_full_pre%in%TRUE&t$candidate_cell),sum(t$gap_reference%in%TRUE),sum(t$gap_reference%in%TRUE&t$candidate_cell)))
write.csv(overall,file.path(outdir,"candidate_support_overall.csv"),row.names=FALSE); write.csv(cc,file.path(outdir,"candidate_support_cells.csv"),row.names=FALSE)
summary_by <- function(x,g) x[,.(target_rows=.N,candidate_rows=sum(candidate_cell),eligible_candidate_rows=sum(eligible_candidate_cell),exact_cells=uniqueN(cell),candidate_cells=uniqueN(cell[candidate_cell]),fixed_full_pre_rows=sum(gap_full_pre%in%TRUE),fixed_full_pre_candidate_rows=sum(gap_full_pre%in%TRUE&candidate_cell),reference_rows=sum(gap_reference%in%TRUE),reference_candidate_rows=sum(gap_reference%in%TRUE&candidate_cell)),by=g][order(get(g))]
write.csv(summary_by(t,"target_event_time"),file.path(outdir,"candidate_support_by_event_time.csv"),row.names=FALSE); write.csv(summary_by(t,"target_year"),file.path(outdir,"candidate_support_by_target_year.csv"),row.names=FALSE); write.csv(summary_by(t,"target_STATEFIP"),file.path(outdir,"candidate_support_by_statefip.csv"),row.names=FALSE)
jsonlite::write_json(list(status="EXACT_CELL_AVAILABILITY_ONLY",matching_status="not_run",assignment_status="not_run",coarsening_status="not_run",target_file=tf,donor_file=df,exact_cell=c("target_year=donor YEAR","target_mother_age=donor AGE","target_child_age=donor sole_child_age",paste0("target_",covars,"=donor ",covars)),fixed_full_pre_pool="gap_full_pre TRUE; reported separately",reference_pool="gap_reference TRUE; reported separately",generated=as.character(Sys.time())),file.path(outdir,"candidate_support_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
ck("complete",sprintf("target_rows=%s candidate_rows=%s",nrow(t),sum(t$candidate_cell)))
