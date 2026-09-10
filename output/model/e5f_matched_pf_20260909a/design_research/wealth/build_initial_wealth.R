#!/usr/bin/env Rscript
# Diagnostic only: selected-column read, exact existing sample definitions,
# initial-period windows and person-cluster uncertainty. Never writes targets.
suppressPackageStartupMessages({library(haven); library(data.table)})
setDTthreads(1L)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L, dir.exists(args[1L]))
out <- normalizePath(args[1L])
raw_path <- '/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta'
message('START ', format(Sys.time(),tz='UTC'), ' one selected-column raw pass')
cols <- c('ID','year','AGEREP','DEATHYEAR','RELTOHEAD_','RELCHINUM',
          'INCFAMR','NETWORTHR','IW','EARNINDRRC')
raw <- as.data.table(read_dta(raw_path,col_select=all_of(cols)))
for(n in names(raw)) set(raw,j=n,value=as.numeric(raw[[n]]))
setnames(raw,cols,c('id','year','age','death_year','relation_to_head',
                    'n_children_num','income','total_nw','weight','gross_earnings'))
message('READ_DONE ', nrow(raw),' rows ',format(Sys.time(),tz='UTC'))
# Exactly the audit_intergen_bequest_family_size_targets.R base filter.
oldbase <- raw[year>=1984 & year<=2019 &
   (is.na(death_year)|year<=death_year) & age>=65 & age<=84 &
   is.finite(weight) & weight>0 & relation_to_head==10]
# Exactly audit_aggregate_wealth_earnings_ratio.R, with date window later.
aggbase <- raw[age>=18 & age<=85 & relation_to_head==10 &
   is.finite(weight) & weight>0 & is.finite(total_nw) &
   (age>65 | (is.finite(gross_earnings) & gross_earnings>=0))]
raw_n <- nrow(raw); rm(raw); invisible(gc())
wq <- function(x,w,p){
  ok <- is.finite(x)&is.finite(w)&w>0
  x<-x[ok];w<-w[ok]; ord<-order(x);x<-x[ord];w<-w[ord]
  if(!length(x)) return(NA_real_)
  x[which(cumsum(w)/sum(w)>=p)[1L]]
}
oldsample <- function(x){
  x<-copy(x[age>=76 & age<=84 & is.finite(n_children_num)])
  x[,ratio:=total_nw/income]
  x[is.finite(ratio)&income>1000]
}
oldstat <- function(x){
  x<-oldsample(x);p50<-wq(x$ratio,x$weight,.5);p90<-wq(x$ratio,x$weight,.9)
  c(old_p90_p50=p90/p50, old_p50=p50, old_p90=p90)
}
unc <- function(x)c(se=sd(x,na.rm=TRUE),p025=unname(quantile(x,.025,na.rm=TRUE)),
                   p975=unname(quantile(x,.975,na.rm=TRUE)))
old_windows <- list(pooled_1984_2019=c(1984,2019),initial_2005_2007=c(2005,2007),
                    initial_2003_2007=c(2003,2007),initial_2007=c(2007,2007),
                    preannouncement_2003_2005=c(2003,2005),
                    preannouncement_2005=c(2005,2005))
old_results<-list(); old_draws<-list()
for(label in names(old_windows)){
  bounds<-old_windows[[label]]; x<-copy(oldbase[year>=bounds[1]&year<=bounds[2]])
  s<-oldsample(x);point<-oldstat(x)
  # Same sampling population as the authoritative builder: all living weighted
  # reference persons ages 65--84, before the outcome/children-observed filters.
  ids<-unique(x$id); reps<-499L; set.seed(20260715L)
  draws<-matrix(NA_real_,nrow=reps,ncol=3,dimnames=list(NULL,names(point)))
  for(b in seq_len(reps)){
    freq<-tabulate(sample.int(length(ids),length(ids),replace=TRUE),nbins=length(ids))
    boot<-copy(x);boot[,weight:=weight*freq[match(id,ids)]]
    boot<-boot[weight>0];draws[b,]<-oldstat(boot)
  }
  old_results[[label]]<-rbindlist(lapply(seq_along(point),function(j){
    z<-unc(draws[,j]);data.table(window=label,year_min=bounds[1],year_max=bounds[2],
     moment=names(point)[j],estimate=point[j],bootstrap_se=z['se'],
     bootstrap_p025=z['p025'],bootstrap_p975=z['p975'],bootstrap_reps=reps,
     seed=20260715L,family_years=nrow(s),persons=uniqueN(s$id),
     bootstrap_population_persons=length(ids))
  }))
  old_draws[[label]]<-data.table(window=label,rep=seq_len(reps),draws)
  message('OLD_DONE ',label,' point ',point[1],' se ',sd(draws[,1]))
}
fwrite(rbindlist(old_results),file.path(out,'old_wealth_results.csv'))
fwrite(rbindlist(old_draws),file.path(out,'old_wealth_bootstrap_draws.csv'))
agg_windows<-list(pooled_2005_2019=c(2005,2019),initial_2005_2007=c(2005,2007),
                 initial_2003_2007=c(2003,2007),initial_2007=c(2007,2007),
                 preannouncement_2003_2005=c(2003,2005),
                 preannouncement_2005=c(2005,2005))
agg_results<-list();agg_draws<-list();yearly_results<-list()
for(label in names(agg_windows)){
  bounds<-agg_windows[[label]];x<-copy(aggbase[year>=bounds[1]&year<=bounds[2]])
  by_id<-x[,.(num=sum(weight*total_nw),den=sum(weight[age<=65]*gross_earnings[age<=65])),by=id]
  point<-sum(by_id$num)/sum(by_id$den);reps<-999L;set.seed(20260723L)
  draws<-numeric(reps);nids<-nrow(by_id)
  for(b in seq_len(reps)){
    freq<-tabulate(sample.int(nids,nids,replace=TRUE),nbins=nids)
    draws[b]<-sum(freq*by_id$num)/sum(freq*by_id$den)
  }
  z<-unc(draws)
  agg_results[[label]]<-data.table(window=label,year_min=bounds[1],year_max=bounds[2],
   moment='aggregate_wealth_gross_labor_earnings',estimate=point,
   bootstrap_se=z['se'],bootstrap_p025=z['p025'],bootstrap_p975=z['p975'],
   bootstrap_reps=reps,seed=20260723L,family_years=nrow(x),persons=uniqueN(x$id),
   earnings_family_years=sum(x$age<=65),earnings_persons=uniqueN(x[age<=65,id]))
  agg_draws[[label]]<-data.table(window=label,rep=seq_len(reps),draw=draws)
  yearly_results[[label]]<-x[,.(window=label,wealth_households=.N,
   earnings_households=sum(age<=65),aggregate_wealth=sum(weight*total_nw),
   aggregate_gross_labor_earnings=sum(weight[age<=65]*gross_earnings[age<=65])),by=year]
  message('AGG_DONE ',label,' point ',point,' se ',sd(draws))
}
fwrite(rbindlist(agg_results),file.path(out,'aggregate_wealth_results.csv'))
fwrite(rbindlist(agg_draws),file.path(out,'aggregate_wealth_bootstrap_draws.csv'))
fwrite(rbindlist(yearly_results),file.path(out,'aggregate_yearly_components.csv'))
meta<-data.table(raw_path=raw_path,raw_size=file.info(raw_path)$size,
 raw_mtime_utc=format(file.info(raw_path)$mtime,tz='UTC'),raw_rows=raw_n,
 selected_columns=paste(cols,collapse=';'),R_version=R.version.string,
 haven_version=as.character(packageVersion('haven')),
 data_table_version=as.character(packageVersion('data.table')),
 completed_utc=format(Sys.time(),tz='UTC'))
fwrite(meta,file.path(out,'run_metadata.csv'))
message('COMPLETE ',format(Sys.time(),tz='UTC'))
