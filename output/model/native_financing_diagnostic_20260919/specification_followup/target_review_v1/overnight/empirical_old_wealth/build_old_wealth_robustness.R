#!/usr/bin/env Rscript
# Task 3 only: reproduce the 2003/05 old-wealth estimand and two prespecified
# robustness variants. No model, target, or calibration code is changed.
suppressPackageStartupMessages({library(haven); library(data.table)})
setDTthreads(1L)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
out <- normalizePath(args[1L], mustWork=FALSE)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
raw_path <- '/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta'
builder <- '/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/design_research/wealth/build_initial_wealth.R'
stopifnot(file.exists(raw_path), file.exists(builder))
cols <- c('ID','year','AGEREP','DEATHYEAR','RELTOHEAD_','RELCHINUM',
          'INCFAMR','NETWORTHR','IW','EARNINDRRC')
message('READ_START ', format(Sys.time(),tz='UTC'))
raw <- as.data.table(read_dta(raw_path, col_select=all_of(cols)))
for(n in names(raw)) set(raw,j=n,value=as.numeric(raw[[n]]))
setnames(raw,cols,c('id','year','age','death_year','relation_to_head',
                    'n_children','income','wealth','weight','gross_earnings'))
message('READ_DONE rows=',nrow(raw),' ',format(Sys.time(),tz='UTC'))
# Match build_initial_wealth.R oldbase and the preannouncement_2003_2005 window.
base <- raw[year>=1984 & year<=2019 & (is.na(death_year)|year<=death_year) &
            age>=65 & age<=84 & is.finite(weight) & weight>0 & relation_to_head==10 &
            year %in% c(2003,2005)]
rm(raw); invisible(gc())
base_n <- nrow(base)
sample_for <- function(x, children_filter=TRUE) {
  if(children_filter) x <- x[is.finite(n_children)]
  x <- copy(x[age>=76 & age<=84])
  x[, ratio:=wealth/income]
  x[is.finite(ratio) & is.finite(wealth) & is.finite(income) & income>1000]
}
wq <- function(x,w,p) {
  ok <- is.finite(x)&is.finite(w)&w>0
  x<-x[ok]; w<-w[ok]; o<-order(x); x<-x[o]; w<-w[o]
  if(!length(x)) return(NA_real_)
  x[which(cumsum(w)/sum(w)>=p)[1L]]
}
stats <- function(x, children_filter=TRUE) {
  x <- sample_for(x,children_filter)
  c(ratio_p90_p50=wq(x$ratio,x$weight,.9)/wq(x$ratio,x$weight,.5),
    ratio_p50=wq(x$ratio,x$weight,.5), ratio_p90=wq(x$ratio,x$weight,.9),
    wealth_p90_p50=wq(x$wealth,x$weight,.9)/wq(x$wealth,x$weight,.5),
    wealth_p50=wq(x$wealth,x$weight,.5), wealth_p90=wq(x$wealth,x$weight,.9))
}
variants <- list(exact_builder_sample=function(x)stats(x,TRUE),
                 no_children_history_filter=function(x)stats(x,FALSE))
ids <- unique(base$id); reps <- 499L; set.seed(20260715L)
draws <- matrix(NA_real_,nrow=reps,ncol=length(unlist(lapply(variants,function(f)f(base)))),
                dimnames=list(NULL,unlist(lapply(names(variants),function(n)paste(n,names(variants[[n]](base)),sep='__')))))
points <- unlist(lapply(names(variants),function(n)setNames(variants[[n]](base),paste(n,names(variants[[n]](base)),sep='__'))))
for(b in seq_len(reps)) {
  freq <- tabulate(sample.int(length(ids),length(ids),replace=TRUE),nbins=length(ids))
  boot <- copy(base); boot[,weight:=weight*freq[match(id,ids)]]; boot<-boot[weight>0]
  draws[b,] <- unlist(lapply(names(variants),function(n)variants[[n]](boot)))
}
unc <- function(x)c(se=sd(x,na.rm=TRUE),p025=unname(quantile(x,.025,na.rm=TRUE)),p975=unname(quantile(x,.975,na.rm=TRUE)))
sample_counts <- rbindlist(lapply(c(TRUE,FALSE),function(cf){
  x <- base[age>=76 & age<=84]
  after_child <- if(cf) x[is.finite(n_children)] else x
  after_income <- after_child[is.finite(income)&income>1000]
  after_ratio <- after_income[is.finite(wealth/income)]
  data.table(variant=if(cf)'exact_builder_sample' else 'no_children_history_filter',
    base_ages76_84=nrow(x),child_history_observed=nrow(x[is.finite(n_children)]),
    after_child_filter=nrow(after_child),income_missing_or_nonfinite=sum(!is.finite(after_child$income)),
    income_zero=sum(after_child$income==0,na.rm=TRUE),income_negative=sum(after_child$income<0,na.rm=TRUE),
    income_le_1000=sum(is.finite(after_child$income)&after_child$income<=1000),
    after_income_cut=nrow(after_income),wealth_missing_or_nonfinite=sum(!is.finite(after_income$wealth)),
    final_observations=nrow(after_ratio),final_unique_people=uniqueN(after_ratio$id),
    final_weight_sum=sum(after_ratio$weight),person_cluster_population=uniqueN(base$id))
}))
fwrite(data.table(variant=names(points),moment=names(points),estimate=as.numeric(points)),file.path(out,'estimates_wide.csv'))
result <- rbindlist(lapply(seq_along(points),function(j){
  z<-unc(draws[,j]); key<-names(points)[j]; parts<-strsplit(key,'__',fixed=TRUE)[[1]]
  s<-sample_counts[variant==parts[1]]
  data.table(variant=parts[1],moment=parts[2],estimate=unname(points[j]),bootstrap_se=z['se'],
    bootstrap_p025=z['p025'],bootstrap_p975=z['p975'],bootstrap_reps=reps,seed=20260715L,
    family_years=s$final_observations,persons=s$final_unique_people,
    bootstrap_population_persons=s$person_cluster_population)
}))
fwrite(result,file.path(out,'wealth_results.csv'))
fwrite(data.table(rep=seq_len(reps),draws),file.path(out,'bootstrap_draws.csv'))
fwrite(sample_counts,file.path(out,'sample_counts.csv'))
# EARNINDRRC is an earnings variable, not a complete family non-asset-income
# decomposition. Do not use INCFAMR - EARNINDRRC as that denominator.
income_note <- data.table(status='not_computed_component_identity_unavailable',
  reason='The builder carries INCFAMR and EARNINDRRC only. Household gross earnings do not identify total non-asset family income because spouse earnings, pensions, and transfers are not separately observed in this selected-column shelf contract.',
  denominator='INCFAMR')
fwrite(income_note,file.path(out,'nonasset_income_status.csv'))
meta <- data.table(raw_path=raw_path,raw_size=file.info(raw_path)$size,
  raw_mtime_utc=format(file.info(raw_path)$mtime,tz='UTC'),raw_rows_read=base_n,
  raw_columns=paste(cols,collapse=';'),builder_path=builder,
  bootstrap_reps=reps,bootstrap_seed=20260715L,R_version=R.version.string,
  haven_version=as.character(packageVersion('haven')),
  data_table_version=as.character(packageVersion('data.table')),
  completed_utc=format(Sys.time(),tz='UTC'))
fwrite(meta,file.path(out,'run_metadata.csv'))
message('COMPLETE ',format(Sys.time(),tz='UTC'),' exact=',points['exact_builder_sample__ratio_p90_p50'])
