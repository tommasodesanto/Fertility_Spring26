#!/usr/bin/env Rscript
# Readout-only diagnostics for the completed second-birth housing run.
suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(fixest) })
args <- commandArgs(trailingOnly=TRUE)
script_path <- if(length(args)) args[[1L]] else "run_second_birth_housing_readout.R"
script_dir <- dirname(normalizePath(script_path,mustWork=FALSE))
root <- Sys.getenv("KLEVEN_ROOT","/scratch/td2248/projects/kleven_acs_pilot_20260917")
matched_dir <- Sys.getenv("SECOND_BIRTH_MATCHED_DIR",file.path(root,"output/kleven_acs_pilot/second_birth_transformed_support_18080547_corrected"))
housing_dir <- Sys.getenv("SECOND_BIRTH_HOUSING_OUTDIR",file.path(root,"output/kleven_acs_pilot/second_birth_housing_20260920_job_a218afc6"))
outdir <- Sys.getenv("SECOND_BIRTH_READOUT_OUTDIR",file.path(root,"output/kleven_acs_pilot/second_birth_housing_readout_20260920"))
estimator_file <- Sys.getenv("SECOND_BIRTH_HOUSING_FILE",file.path(script_dir,"second_birth_housing.R"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
progress <- file.path(outdir,"progress.log")
checkpoint <- function(stage,detail=""){z<-sprintf("%s\t%s\t%s",format(Sys.time(),"%FT%T%z"),stage,detail);write(z,progress,append=TRUE);cat(z,"\n")}
fail <- function(...) stop(paste0(...),call.=FALSE)
req <- function(x,cols,label){m<-setdiff(cols,names(x));if(length(m))fail(label," missing: ",paste(m,collapse=","))}
key4 <- function(d) paste(d$YEAR,d$SAMPLE,d$SERIAL,d$PERNUM,sep="\034")
raw_col <- function(d,candidates,label){for(z in candidates){h<-names(d)[toupper(names(d))==toupper(z)];if(length(h)==1L)return(h)};fail(label," missing")}
files <- c(estimator_file,file.path(matched_dir,c("proxy_checkpoint.rds","matched_checkpoint.rds","packet_with_author_cells.rds")))
if(any(!file.exists(files)))fail("missing readout input: ",paste(files[!file.exists(files)],collapse=","))
if(!file.exists(file.path(housing_dir,"primary","support.csv")))fail("completed housing support receipt missing")
source(estimator_file,local=TRUE); checkpoint("startup",paste0("housing_dir=",housing_dir))
proxy <- readRDS(file.path(matched_dir,"proxy_checkpoint.rds")); matched <- readRDS(file.path(matched_dir,"matched_checkpoint.rds")); packet <- as.data.table(readRDS(file.path(matched_dir,"packet_with_author_cells.rds")))
source_rows <- as.data.table(copy(proxy$input)); req(source_rows,c("person_key","YEAR","SAMPLE","SERIAL","PERNUM"),"proxy input"); req(packet,c("YEAR","SAMPLE","SERIAL","PERNUM"),"packet")
pk <- key4(packet); sk <- key4(source_rows); if(anyDuplicated(pk)||anyDuplicated(source_rows$person_key))fail("source identity is not unique"); ix <- match(sk,pk); if(anyNA(ix))fail("proxy input source key absent from packet")
raw_map <- c(ROOMS=raw_col(packet,c("ROOMS_RAW","ROOMS"),"ROOMS"),BEDROOMS=raw_col(packet,c("BEDROOMS_RAW","BEDROOMS"),"BEDROOMS"),OWNERSHP=raw_col(packet,c("OWNERSHP_RAW","OWNERSHP","OWNERSHIP"),"OWNERSHP"))
for(nm in names(raw_map))source_rows[[paste0("raw_",tolower(nm))]]<-packet[[raw_map[[nm]]]][ix]
rr<-suppressWarnings(as.numeric(source_rows$raw_rooms));bb<-suppressWarnings(as.numeric(source_rows$raw_bedrooms));oo<-suppressWarnings(as.numeric(source_rows$raw_ownershp))
source_rows[,rooms9:=ifelse(!is.na(rr)&rr%in%c(1:27,30),pmin(rr,9),NA_real_)];source_rows[,bedrooms5:=ifelse(!is.na(bb)&bb%in%1:22,pmin(ifelse(bb==22,21,bb-1),5),NA_real_)];source_rows[,ownership_lw:=ifelse(!is.na(oo)&oo%in%c(1,2),ifelse(oo==1,1,0),NA_real_)]
proxy$input<-source_rows; checkpoint("source_loaded",paste0("rows=",nrow(source_rows)," years=",min(source_rows$YEAR),"-",max(source_rows$YEAR)))

specs<-list(primary=list(support_spec="all_anchors",fertyr_spec="all"),joint_negative=list(support_spec="joint_negative",fertyr_spec="all"));outcomes<-c("rooms9","bedrooms5","ownership_lw");prepared<-list()
for(nm in names(specs)){prepared[[nm]]<-prepare_second_birth_housing(proxy,matched,source_rows=source_rows,outcomes=outcomes,support_spec=specs[[nm]]$support_spec,fertyr_spec=specs[[nm]]$fertyr_spec);write_json(list(specification=nm,support=prepared[[nm]]$support,contract=prepared[[nm]]$contract),file.path(outdir,paste0(nm,"_prepared_receipt.json")),auto_unbox=TRUE,pretty=TRUE)}

# Replace the buggy curve support metadata by exact outcome/event rows.
curve_rows<-list();for(nm in c("primary","joint_negative","fertyr_event0_yes")){cc<-fread(file.path(housing_dir,nm,"curves.csv"));ss<-fread(file.path(housing_dir,nm,"support.csv"));cc[,event_time:=as.integer(event_time)];ss[,event_time:=as.integer(event_time)];stopifnot(!anyDuplicated(ss[,.(outcome,event_time)]));curve_rows[[nm]]<-merge(cc,ss[,.(outcome,event_time,support_n_rows=n_rows,support_n_observed=n_observed,support_weight_ess=weight_ess,support_clusters=source_household_clusters)],by=c("outcome","event_time"),all.x=TRUE,sort=FALSE)}
fwrite(rbindlist(curve_rows,idcol="specification"),file.path(outdir,"corrected_saved_curves.csv"))

# Weighted age/year/state/marriage/education distributions by event and role.
category_rows<-list();add_categories<-function(d,spec,role){if(!nrow(d))return(NULL);x<-copy(d);src<-source_rows[match(x$person_key,source_rows$person_key)];for(v in c("MARST","EDUC")){h<-names(source_rows)[toupper(names(source_rows))==v];x[[v]]<-if(length(h))src[[h]] else NA};vars<-list(age=as.character(x$AGE_norm),year=as.character(x$YEAR),state=as.character(x$STATEFIP),marriage=as.character(x$MARST),education=as.character(x$EDUC));rbindlist(lapply(names(vars),function(v){x[,category:=vars[[v]]];x[!is.na(category),.(weighted_rows=sum(weight)),by=.(event_time,category)][,`:=`(variable=v,specification=spec,pseudo_role=role,share=weighted_rows/sum(weighted_rows)),by=event_time]}),fill=TRUE)}
for(nm in names(prepared)){p<-prepared[[nm]];category_rows[[paste0(nm,"_post")]]<-add_categories(p$post,nm,"post");category_rows[[paste0(nm,"_donor")]]<-add_categories(p$negative_donor,nm,"negative_donor")}
fwrite(rbindlist(category_rows,fill=TRUE),file.path(outdir,"weighted_composition_distributions.csv"))

# Donor-link gap distribution retains the anchor gap before donor aggregation.
anchors<-as.data.table(copy(proxy$anchors));req(anchors,c("person_key","birth_gap","FERTYR_status"),"proxy anchors");primary_anchor_ids<-anchors[birth_gap>=2&is.finite(birth_gap),person_key];links<-as.data.table(copy(matched$links))[target_event_time%in%c(-2L,-1L)&anchor_person_key%in%primary_anchor_ids];links<-merge(links,anchors[,.(anchor_person_key=person_key,birth_gap)],by="anchor_person_key",all.x=TRUE);links[,weight:=as.numeric(donor_PERWT)*as.numeric(wgt_match)];gap<-links[is.finite(weight)&weight>0,.(weighted_rows=sum(weight)),by=.(event_time=target_event_time,birth_gap)][,`:=`(variable="birth_gap",pseudo_role="negative_donor",specification="primary_gap_ge_2",share=weighted_rows/sum(weighted_rows)),by=event_time];fwrite(gap,file.path(outdir,"weighted_gap_distribution_primary.csv"))
fert<-rbindlist(lapply(names(prepared),function(nm){ids<-prepared[[nm]]$support$selected_ids;z<-anchors[person_key%in%ids,.(anchors=.N),by=FERTYR_status];z[,`:=`(specification=nm,selected_anchor_count=length(ids),all_gap_ge_2_anchor_count=prepared[[nm]]$support$n_all_anchors,eligible_definition="strict gap>=2 anchors; NCHILD=2 and exactly two linked children")]}),fill=TRUE);fwrite(fert,file.path(outdir,"fertyr_eligible_counts.csv"))

# Same-sample sequential specification readout; no main fit receipt is changed.
seq_rows<-list();for(y in outcomes){d<-copy(prepared$primary$panel);d<-d[is.finite(as.numeric(get(y)))&is.finite(weight)&weight>0];base_nobs<-nrow(d);d[,event_factor:=factor(event_time,levels=c(-2L,-1L,0L,1L,2L,3L))];for(m in c("event_only","state_year_fe","age_fe","state_age_year_fe")){f<-switch(m,event_only=as.formula(sprintf("%s ~ i(event_factor, ref = '-2')",y)),state_year_fe=as.formula(sprintf("%s ~ i(event_factor, ref = '-2') | STATEFIP + YEAR",y)),age_fe=as.formula(sprintf("%s ~ i(event_factor, ref = '-2') | AGE_norm",y)),state_age_year_fe=as.formula(sprintf("%s ~ i(event_factor, ref = '-2') | STATEFIP + AGE_norm + YEAR",y)));z<-tryCatch(feols(f,data=d,weights=~weight,cluster=~source_household_cluster),error=function(e)e);if(inherits(z,"error")){seq_rows[[length(seq_rows)+1L]]<-data.table(outcome=y,model=m,status="FIT_FAILED",base_nobs=base_nobs,error=conditionMessage(z));next};b<-coef(z);V<-as.matrix(vcov(z));terms<-paste0("event_factor::",c("3","-1"));if(any(!terms%in%names(b))){seq_rows[[length(seq_rows)+1L]]<-data.table(outcome=y,model=m,status="MISSING_EVENT_TERM",base_nobs=base_nobs,error="event coefficient omitted");next};q<-c(1,-1);est<-sum(q*b[terms]);vv<-as.numeric(t(q)%*%V[terms,terms]%*%q);seq_rows[[length(seq_rows)+1L]]<-data.table(outcome=y,model=m,status="FIT_COMPLETE",base_nobs=base_nobs,nobs=stats::nobs(z),nobs_diff_from_base=stats::nobs(z)!=base_nobs,estimate=est,std_error=sqrt(max(0,vv)),raw_event_neg1=sum(d[event_time==-1,weight]*d[event_time==-1,get(y)])/sum(d[event_time==-1,weight]),raw_event_3=sum(d[event_time==3,weight]*d[event_time==3,get(y)])/sum(d[event_time==3,weight]))}}
fwrite(rbindlist(seq_rows,fill=TRUE),file.path(outdir,"same_sample_specification_readout.csv"))

# Gate the common implied event-year sensitivity before its fits.  The gate
# requires every requested event cell (-2,-1,0,+1,+2,+3), rather than merely
# requiring at least one supported cell for each outcome.
common_cohort_gate <- function(panel, outcomes, cohort_min=2007L,
                               cohort_max=2016L,
                               requested_events=c(-2L,-1L,0L,1L,2L,3L)) {
  d0 <- copy(panel)
  d0[,implied_event_year:=as.integer(YEAR-event_time)]
  d0 <- d0[implied_event_year>=cohort_min & implied_event_year<=cohort_max]
  rows <- rbindlist(lapply(outcomes,function(y){
    d <- copy(d0)
    d[,observed:=is.finite(as.numeric(get(y))) & is.finite(weight) & weight>0]
    d[,.(n_rows=.N,n_observed=sum(observed)),by=.(implied_event_year,event_time)][,outcome:=y]
  }),fill=TRUE)
  rows <- merge(CJ(outcome=outcomes,
                   implied_event_year=seq.int(cohort_min,cohort_max),
                   event_time=requested_events),rows,
                by=c("outcome","implied_event_year","event_time"),
                all.x=TRUE,sort=FALSE)
  rows[is.na(n_rows),n_rows:=0L]
  rows[is.na(n_observed),n_observed:=0L]
  rows[,supported:=n_observed>0L]
  summary <- rows[,.(requested_event_cells=(cohort_max-cohort_min+1L)*length(requested_events),
                     observed_event_cells=sum(n_rows>0L),
                     supported_event_cells=sum(supported),
                     all_six_supported=.N==(cohort_max-cohort_min+1L)*length(requested_events) &&
                       all(supported)),by=outcome]
  summary[,`:=`(cohort_min=cohort_min,cohort_max=cohort_max)]
  list(panel=d0,rows=rows,summary=summary)
}
common<-copy(prepared$primary);gate_obj<-common_cohort_gate(common$panel,outcomes);common$panel<-gate_obj$panel
gate<-gate_obj$rows;gate_receipt<-gate_obj$summary
fwrite(gate,file.path(outdir,"common_cohort_support_gate.csv"));fwrite(gate_receipt,file.path(outdir,"common_cohort_support_gate_summary.csv"))
if(!all(gate_receipt$all_six_supported)){write_json(list(status="SUPPORT_INCOMPLETE",gate=gate_receipt),file.path(outdir,"common_cohort_result_receipt.json"),auto_unbox=TRUE,pretty=TRUE)}else{cr<-list();for(y in outcomes){z<-tryCatch(fit_second_birth_housing(common,outcomes=y),error=function(e)e);if(inherits(z,"error"))cr[[y]]<-data.table(outcome=y,status="FIT_FAILED",error=conditionMessage(z))else{q<-z$contrasts[[y]];cr[[y]]<-data.table(outcome=y,status="FIT_COMPLETE",estimate=q$estimate,std_error=q$se,lower=q$lower,upper=q$upper,nobs=z$fits[[y]]$nobs)}};fwrite(rbindlist(cr,fill=TRUE),file.path(outdir,"common_cohort_contrasts.csv"));write_json(list(status="COMPLETE",cohort=2007:2016,gate=gate_receipt),file.path(outdir,"common_cohort_result_receipt.json"),auto_unbox=TRUE,pretty=TRUE)}
write_json(list(status="COMPLETE",housing_run_dir=housing_dir,fit_rerun_main=FALSE,composition="weighted distributions from saved prepared rows",fertyr_counts="fertyr_eligible_counts.csv",common_cohort="gate written before sensitivity fit"),file.path(outdir,"readout_receipt.json"),auto_unbox=TRUE,pretty=TRUE);checkpoint("complete")
