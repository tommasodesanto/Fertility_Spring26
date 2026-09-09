#!/usr/local/bin/Rscript
# Diagnostic only. Exact masks from build_active_acs_room_target_receipt.R:131-150.
# One cache read, one data.table thread, no raw-data processing or bootstrap.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
started <- proc.time()[['elapsed']]
root <- normalizePath(getwd())
out <- file.path(root, 'output/model/e5f_matched_pf_20260909a/meeting_receipts/acs_date_diagnostic')
cache_path <- file.path(root, 'code/data/moment_standard_errors/cache/acs_analysis_samples.rds')
expected_sha <- '0eae9aaed4e1d3b9655235be967378c1d389f29ec9b64f005a811c2f0ced7df0'
hash_line <- system2('/usr/bin/shasum', c('-a', '256', shQuote(cache_path)), stdout=TRUE)
observed_sha <- strsplit(hash_line[[1L]], '[[:space:]]+')[[1L]][1L]
stopifnot(identical(observed_sha, expected_sha))
cache <- readRDS(cache_path)
stopifnot(is.list(cache), all(c('targets', 'women') %in% names(cache)))
dt <- as.data.table(cache$targets)
required <- c('year','met2013','age','pernum','hhwt','ownershp','rooms',
              'nchild','yngch','parent_u18','owner','renter','mms_location')
stopifnot(all(required %in% names(dt)), nrow(dt)==4103889L,
          min(dt$year)==2012L, max(dt$year)==2023L, uniqueN(dt$met2013)==42L,
          all(dt$pernum==1L), all(is.finite(dt$hhwt)&dt$hhwt>0),
          all(dt$ownershp %in% c(1L,2L)), all(is.finite(dt$rooms)&dt$rooms>0))
# Check, but do not redefine, the cached authoritative parent flag.
expected_parent <- dt$nchild>0 & !is.na(dt$yngch) & dt$yngch!=99 & dt$yngch<18
stopifnot(identical(as.logical(dt$parent_u18), as.logical(expected_parent)))
wm <- function(x) sum(x$rooms*x$hhwt)/sum(x$hhwt)
rows <- list()
for (era in c('pooled_2012_2023', 'year_2023')) {
  sample <- if (era=='year_2023') dt[year==2023L] else dt
  mean_sample <- sample[age>=18L & age<=85L]
  gap_sample <- mean_sample[age>=30L & age<=55L & parent_u18==TRUE & nchild>=1]
  large <- gap_sample[nchild>=3]
  small <- gap_sample[nchild>=1 & nchild<=2]
  stopifnot(nrow(large)>0, nrow(small)>0,
            nrow(large)+nrow(small)==nrow(gap_sample))
  for (key in c('aggregate_mean_occupied_rooms_18_85',
                'prime30_55_parent_3plus_minus_1to2_mean_rooms')) {
    is_gap <- key=='prime30_55_parent_3plus_minus_1to2_mean_rooms'
    group <- if (is_gap) large else mean_sample
    control <- if (is_gap) small else mean_sample[0L]
    group_mean <- wm(group)
    control_mean <- if (is_gap) wm(control) else NA_real_
    estimate <- if (is_gap) group_mean-control_mean else group_mean
    rows[[length(rows)+1L]] <- data.table(moment=key, period=era, estimate=estimate,
      sample_n=if(is_gap)nrow(gap_sample) else nrow(mean_sample),
      sample_hhwt=if(is_gap)sum(gap_sample$hhwt) else sum(mean_sample$hhwt),
      group_n=nrow(group), group_hhwt=sum(group$hhwt), group_mean_rooms=group_mean,
      control_n=nrow(control), control_hhwt=sum(control$hhwt), control_mean_rooms=control_mean,
      metro_count=uniqueN(if(is_gap)gap_sample$met2013 else mean_sample$met2013),
      cache_sha256=observed_sha, status='diagnostic_only_no_new_target_or_uncertainty')
  }
}
result <- rbindlist(rows)
receipt <- fread(file.path(root, 'code/data/moment_standard_errors/output_active_acs_room_target_receipt_20260817/target_receipt.csv'))
for (moment_key in result$moment[result$period=='pooled_2012_2023']) {
  actual <- result[moment==moment_key & period=='pooled_2012_2023']
  reference <- receipt[get('key')==moment_key]
  stopifnot(nrow(reference)==1L, abs(actual$estimate-reference$reproduced_point)<=1e-10,
            abs(actual$estimate-reference$target)<=1e-10,
            actual$sample_n==reference$n_unweighted,
            actual$sample_hhwt==reference$hhwt_sum)
  result[moment==moment_key, pooled_reproduction_abs_gap:=abs(actual$estimate-reference$reproduced_point)]
}
# Publish annual results only after both pooled samples and points pass.
fwrite(result, file.path(out, 'rooms_date_extraction.csv'), na='')
writeLines(c('status=complete_diagnostic_only', 'cache_read_count=1',
             paste0('cache_sha256=',observed_sha), paste0('threads=',getDTthreads()),
             'raw_extract_read=false', 'bootstrap=false', 'pooled_tolerance=1e-10',
             paste0('elapsed_seconds=',proc.time()[['elapsed']]-started)),
           file.path(out,'rooms_date_extraction_receipt.txt'))
print(result[,.(moment,period,estimate,sample_n,sample_hhwt,pooled_reproduction_abs_gap)])
