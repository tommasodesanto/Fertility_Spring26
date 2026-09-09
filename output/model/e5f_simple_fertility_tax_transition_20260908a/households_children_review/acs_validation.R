#!/usr/local/bin/Rscript
# Aggregate-only validation from the existing pinned cache; no raw-data reread.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
setDTthreads(1L)
started <- proc.time()[['elapsed']]
Sys.setenv(LC_ALL='C')
root <- normalizePath(getwd())
folder <- file.path(root,'output/model/e5f_simple_fertility_tax_transition_20260908a/households_children_review')
cache_path <- file.path(root,'code/data/moment_standard_errors/cache/acs_analysis_samples.rds')
expected <- '0eae9aaed4e1d3b9655235be967378c1d389f29ec9b64f005a811c2f0ced7df0'
receipt_path <- file.path(folder,'acs_validation_receipt.json')
tryCatch({
  hash_line <- system2('/usr/bin/shasum',c('-a','256',shQuote(cache_path)),stdout=TRUE,stderr=TRUE)
  hash_line <- hash_line[grepl('^[0-9a-f]{64}[[:space:]]',hash_line)]
  stopifnot(length(hash_line)==1L)
  observed <- strsplit(hash_line[1L],'[[:space:]]+')[[1L]][1L]
  stopifnot(identical(observed,expected))
  cache <- readRDS(cache_path) # exactly one cache read
  dt <- as.data.table(cache$targets)
  stopifnot(nrow(dt)==4103889L, uniqueN(dt$met2013)==42L,
            min(dt$year)==2012L,max(dt$year)==2023L,
            all(c('nchild','yngch','eldch','hhwt','age','pernum','ownershp','rooms','relate','due_housing') %in% names(dt)),
            all(dt$pernum==1L),all(is.finite(dt$hhwt)&dt$hhwt>0),
            all(dt$ownershp %in% c(1L,2L)),all(is.finite(dt$rooms)&dt$rooms>0))
  # NCHILD==0 establishes absence even when the associated child-age variable
  # is 99 (not applicable). Positive NCHILD with invalid age is unknown.
  nvalid <- is.finite(dt$nchild)&dt$nchild>=0
  positive <- nvalid & dt$nchild>0
  zero <- nvalid & dt$nchild==0
  dt[, any_coresident_own_child := fifelse(nvalid,positive,NA)]
  dt[, no_coresident_own_child := fifelse(nvalid,zero,NA)]
  yngvalid <- is.finite(dt$yngch)&dt$yngch>=0&dt$yngch<99
  eldvalid <- is.finite(dt$eldch)&dt$eldch>=0&dt$eldch<99
  dt[, any_own_child_under18 := fifelse(zero,FALSE,fifelse(positive&yngvalid,yngch<18,NA))]
  dt[, recent_parent_oldest_child_under4 := fifelse(zero,FALSE,fifelse(positive&eldvalid,eldch<4,NA))]
  metrics <- c('any_coresident_own_child','any_own_child_under18','recent_parent_oldest_child_under4','no_coresident_own_child')
  age_groups <- list(all_18_85=c(18L,85L),young_25_34=c(25L,34L),prime_30_55=c(30L,55L))
  rows <- list(); denominator_audit <- list(); i <- 0L
  for(period in c('pooled_2012_2023','year_2023')) for(label in names(age_groups)) {
    ages <- age_groups[[label]]
    sub <- dt[age>=ages[1L]&age<=ages[2L] & (period=='pooled_2012_2023'|year==2023L)]
    H <- sum(sub$hhwt)
    stopifnot(nrow(sub)>0,H>0)
    denominator_audit[[paste(period,label,sep='/')]] <- list(n_records=nrow(sub),hhwt_sum=H,
      metro_count=uniqueN(sub$met2013),years=sort(unique(sub$year)),
      missing_nchild_records=sum(!is.finite(sub$nchild)),
      positive_nchild_yngch_missing_records=sum(sub$nchild>0 & !is.finite(sub$yngch),na.rm=TRUE),
      positive_nchild_yngch99_records=sum(sub$nchild>0 & sub$yngch==99,na.rm=TRUE),
      positive_nchild_eldch_missing_records=sum(sub$nchild>0 & !is.finite(sub$eldch),na.rm=TRUE),
      positive_nchild_eldch99_records=sum(sub$nchild>0 & sub$eldch==99,na.rm=TRUE))
    for(metric in metrics) {
      x <- sub[[metric]]; valid <- !is.na(x); yes <- valid & x
      wy <- sum(sub$hhwt[yes]); wv <- sum(sub$hhwt[valid]); wu <- sum(sub$hhwt[!valid])
      i<-i+1L
      rows[[i]] <- data.table(period=period,age_group=label,metric=metric,
        n_records=nrow(sub),hhwt_total=H,n_valid=sum(valid),hhwt_valid=wv,n_unknown=sum(!valid),hhwt_unknown=wu,
        unknown_weight_share=wu/H,hhwt_true=wy,share_among_valid=if(wv>0) wy/wv else NA_real_,
        share_full_denominator_lower_bound=wy/H,share_full_denominator_upper_bound=(wy+wu)/H)
    }
  }
  result<-rbindlist(rows)
  fwrite(result,file.path(folder,'acs_validation.csv'))
  due <- dt[age>=30&age<=55 & relate==1L & due_housing==TRUE]
  wm <- function(x,w)sum(x*w)/sum(w)
  own_reproduced <- wm(due[recent_parent_oldest_child_under4==TRUE]$owner,due[recent_parent_oldest_child_under4==TRUE]$hhwt)-wm(due[no_coresident_own_child==TRUE]$owner,due[no_coresident_own_child==TRUE]$hhwt)
  parents <- dt[age>=30&age<=55 & parent_u18==TRUE]
  room_reproduced <- wm(parents[nchild>=3]$rooms,parents[nchild>=3]$hhwt)-wm(parents[nchild>=1&nchild<=2]$rooms,parents[nchild>=1&nchild<=2]$hhwt)
  stopifnot(abs(own_reproduced-0.16766167)<1e-8,abs(room_reproduced-0.36769955881)<1e-10)
  jsonlite::write_json(list(status='complete_awaiting_lead_review',cache_path=cache_path,cache_sha256=observed,
    cache_read_count=1L,threads=getDTthreads(),cache_rows=nrow(dt),metro_count=uniqueN(dt$met2013),
    sample='MMS42-metro household-head cache; owner/renter with positive rooms and HHWT; not a national sample; no sex restriction',
    pooling='Pooled estimates weight all household-year records by HHWT; not an average of annual shares or a unique-household panel count',
    age_windows='Exact empirical integer ages inclusive, not model-node rebinning',
    unknown_rule='NCHILD==0 establishes false for child presence/recent parent; positive NCHILD plus child age NA/99/outside0..98 is unknown. Valid-denominator share and full-denominator bounds both reported.',
    nchild_semantics='Co-resident own children of any age; biological/adopted/step; not children ever born',
    yngch99_semantics='No own child present; compatible with NCHILD0, inconsistent/unknown with NCHILD>0',
    denominator_audit=denominator_audit,
    target_reproduction=list(own_family_gap=list(target=.16766167,reproduced=own_reproduced,absolute_gap=abs(own_reproduced-.16766167)),room_gap=list(target=.36769955881,reproduced=room_reproduced,absolute_gap=abs(room_reproduced-.36769955881))),
    elapsed_seconds=proc.time()[['elapsed']]-started,
    no_se_no_new_target_no_parameter_changes=TRUE),receipt_path,pretty=TRUE,auto_unbox=TRUE,digits=16)
  print(result[,.(period,age_group,metric,share_among_valid,unknown_weight_share)])
},error=function(e){
  jsonlite::write_json(list(status='failed',error=conditionMessage(e),elapsed_seconds=proc.time()[['elapsed']]-started),receipt_path,pretty=TRUE,auto_unbox=TRUE)
  stop(e)
})
