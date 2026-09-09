suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
Sys.setenv(LC_ALL='C')
p <- 'code/data/moment_standard_errors/cache/acs_analysis_samples.rds'
h <- system2('/usr/bin/shasum',c('-a','256',shQuote(p)),stdout=TRUE)
stopifnot(startsWith(h,'0eae9aaed4e1d3b9655235be967378c1d389f29ec9b64f005a811c2f0ced7df0'))
c <- readRDS(p)
d <- as.data.table(c$targets)
stopifnot(nrow(d)==4103889L,uniqueN(d$met2013)==42L)
d <- d[age>=30&age<=55&relate==1L&due_housing==TRUE]
meanown <- function(s)sum(s$hhwt*s$owner)/sum(s$hhwt)
rows <- list()
for(era in c('pooled_2012_2023','year_2023')) {
 a <- if(era=='year_2023') d[year==2023L] else d
 control <- a[nchild==0]
 for(group in c('recent_parent_oldest_under4','any_coresident_own_child','any_own_child_under18')) {
  selected <- switch(group,
    recent_parent_oldest_under4=a[nchild>0 & !is.na(eldch) & eldch!=99 & eldch<4],
    any_coresident_own_child=a[nchild>0],
    any_own_child_under18=a[nchild>0 & !is.na(yngch) & yngch!=99 & yngch<18])
  rows[[length(rows)+1L]] <- data.table(period=era,group=group,comparison_group='NCHILD0, including empty nesters',
    age_window='head age30-55, exact integer ages',sample='MMS42metros DUE head positive rooms owner/renter HHWT',
    group_n=nrow(selected),group_weight=sum(selected$hhwt),control_n=nrow(control),control_weight=sum(control$hhwt),
    group_owner_rate=meanown(selected),control_owner_rate=meanown(control),ownership_gap=meanown(selected)-meanown(control))
 }
}
r<-rbindlist(rows)
stopifnot(abs(r[period=='pooled_2012_2023' & group=='recent_parent_oldest_under4',ownership_gap]-.16766167)<1e-8)
stopifnot(abs(r[period=='pooled_2012_2023' & group=='any_coresident_own_child',ownership_gap]-.152116976514673)<1e-12)
fwrite(r,'output/model/e5f_simple_fertility_tax_transition_20260908a/households_children_review/ownership_target_comparison.csv')
print(r[,.(period,group,group_owner_rate,control_owner_rate,ownership_gap)])
