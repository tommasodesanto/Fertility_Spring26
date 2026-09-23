#!/usr/bin/env Rscript

# Read-only replication and denominator comparison for the July 18-24 PSID
# childless-renter entry sample. Run from this directory or pass PSID_PATH.
suppressPackageStartupMessages({library(data.table); library(haven)})

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
out <- file.path(repo, "output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry")
psid <- if (length(args)) args[[1]] else file.path(dirname(repo), "PSID/PSIDSHELF_MOBILITY.dta")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
weighted_q <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) >= p * sum(w))[1]]
}
rank_bin <- function(x, w, n) {
  # Stable weighted-rank bins: each row is assigned by midpoint cumulative mass.
  ok <- is.finite(x) & is.finite(w) & w > 0
  out <- rep(NA_integer_, length(x))
  ix <- which(ok); o <- ix[order(x[ix], ix)]
  mid <- (cumsum(w[o]) - 0.5 * w[o]) / sum(w[o])
  out[o] <- pmin(n, floor(mid * n) + 1L)
  out
}
weighted_spearman <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  midrank <- function(a, ww) {
    mass <- data.table(a=a, ww=ww)[,.(mass=sum(ww)),by=a][order(a)]
    mass[, r := (cumsum(mass)-0.5*mass)/sum(mass)]
    mass$r[match(a,mass$a)]
  }
  rx <- midrank(x[ok], w[ok])
  ry <- midrank(y[ok], w[ok])
  weighted_cov <- function(a,b) sum(w[ok] * (a-weighted_mean(a,w[ok])) * (b-weighted_mean(b,w[ok]))) / sum(w[ok])
  weighted_cov(rx,ry) / sqrt(weighted_cov(rx,rx) * weighted_cov(ry,ry))
}
quintile_nodes <- function(x, w) {
  b <- rank_bin(x, w, 5L)
  rbindlist(lapply(1:5, function(k) {
    ix <- b == k
    data.table(quintile=k, weight_share=sum(w[ix])/sum(w),
               node=weighted_mean(x[ix],w[ix]), rows=sum(ix))
  }))
}
ratio_summary <- function(label, numerator, denominator, w, sample_label, ids=dt$id) {
  ok <- is.finite(numerator) & is.finite(denominator) & denominator > 0 & is.finite(w) & w > 0
  r <- numerator[ok] / denominator[ok]; ww <- w[ok]
  data.table(sample=sample_label, measure=label, family_years=sum(ok),
    persons=uniqueN(ids[ok]), weighted_mean_row_ratio=weighted_mean(r,ww),
    weighted_median_row_ratio=weighted_q(r,ww,.5),
    ratio_of_weighted_sums=sum(ww*numerator[ok])/sum(ww*denominator[ok]))
}

message("Reading PSID shelf columns for the pre-specified July sample...")
raw <- as.data.table(read_dta(psid, col_select=c("ID","year","AGEREP","DEATHYEAR","RELTOHEAD_","RELCHIREP","INCFAMR","NETWORTH2R","HOMEOWN","IW","EARNINDRRC")))
label_rows <- rbindlist(lapply(c("year","INCFAMR","EARNINDRRC"),function(v) {
  data.table(variable=v,label=as.character(attr(raw[[v]],"label")))
}))
fwrite(label_rows,file.path(out,"source_variable_labels.csv"))
dt <- raw[, .(id=as.numeric(ID), year=as.integer(year), age=as.numeric(AGEREP),
  death_year=as.numeric(DEATHYEAR), rp=as.numeric(RELTOHEAD_)==10,
  children=as.numeric(RELCHIREP), family_income=as.numeric(INCFAMR),
  nonhousing_nw=as.numeric(NETWORTH2R), homeown=as.numeric(HOMEOWN),
  weight=as.numeric(IW), gross_labor=as.numeric(EARNINDRRC))]
rm(raw); invisible(gc())
dt <- dt[year>=1984 & year<=2019 & (is.na(death_year) | year<=death_year) & is.finite(weight) & weight>0]
dt <- dt[rp & is.finite(age) & age>=18 & age<=24 & children==0 & homeown==2 &
  is.finite(family_income) & family_income>1000 & is.finite(nonhousing_nw)]
dt[, legacy_ratio := nonhousing_nw/family_income]
dt <- dt[is.finite(legacy_ratio)]
if (nrow(dt)!=1835L || uniqueN(dt$id)!=1346L) stop(sprintf("July sample mismatch: %d family-years, %d people", nrow(dt), uniqueN(dt$id)))

sample_counts <- data.table(
  item=c("July sample family-years","July sample persons","unique waves","age min","age max","total IW weight",
         "gross labor missing","gross labor zero","gross labor positive <=1000","gross labor >1000","family income missing/<=1000 (should be zero)",
         "gross-denominator excluded IW weight","gross-denominator excluded IW weight share","gross-denominator eligible IW weight share",
         "nonhousing net worth negative","nonhousing net worth zero"),
  value=c(nrow(dt),uniqueN(dt$id),uniqueN(dt$year),min(dt$age),max(dt$age),sum(dt$weight),
    sum(!is.finite(dt$gross_labor)),sum(is.finite(dt$gross_labor)&dt$gross_labor==0),
    sum(is.finite(dt$gross_labor)&dt$gross_labor>0&dt$gross_labor<=1000),
    sum(is.finite(dt$gross_labor)&dt$gross_labor>1000),sum(!is.finite(dt$family_income)|dt$family_income<=1000),
    sum(dt$weight[!is.finite(dt$gross_labor)|dt$gross_labor<=1000]),
    sum(dt$weight[!is.finite(dt$gross_labor)|dt$gross_labor<=1000])/sum(dt$weight),
    sum(dt$weight[is.finite(dt$gross_labor)&dt$gross_labor>1000])/sum(dt$weight),
    sum(dt$nonhousing_nw<0),sum(dt$nonhousing_nw==0)))
fwrite(sample_counts,file.path(out,"sample_counts.csv"))

# Reproduce the original cumulative-weight-cut rule exactly for its primary ratio.
orig_bin <- function(x,w) {
  o <- order(x); cw <- cumsum(w[o])/sum(w[o])
  b <- cut(cw,breaks=c(0,.2,.4,.6,.8,1+1e-12),labels=FALSE,include.lowest=TRUE)
  ans <- rep(NA_integer_,length(x)); ans[o] <- b; ans
}
dt[, legacy_bin := orig_bin(legacy_ratio,weight)]
legacy_nodes <- dt[,.(node=sum(weight*legacy_ratio)/sum(weight),weight_share=sum(weight)/sum(dt$weight),family_years=.N,persons=uniqueN(id)),by=legacy_bin][order(legacy_bin)]
setnames(legacy_nodes,"legacy_bin","quintile")
legacy_nodes[,`:=`(sample="July legacy wealth/family-income ratio",measure="weighted quintile-bin mean row ratio")]
fwrite(legacy_nodes,file.path(out,"legacy_entry_nodes.csv"))

# Common-support candidates retain the exact July sample; EARNINDRRC missing or
# <=$1,000 observations are excluded only from gross-labor ratio calculations.
dt[, gross_ratio := fifelse(is.finite(gross_labor)&gross_labor>1000,nonhousing_nw/gross_labor,NA_real_)]
family_only <- ratio_summary("family income",dt$nonhousing_nw,dt$family_income,dt$weight,"all July rows")
common <- is.finite(dt$gross_labor)&dt$gross_labor>1000
common_family <- ratio_summary("family income",dt$nonhousing_nw[common],dt$family_income[common],dt$weight[common],"common gross-denominator support",dt$id[common])
common_gross <- ratio_summary("RP+spouse gross labor earnings",dt$nonhousing_nw[common],dt$gross_labor[common],dt$weight[common],"common gross-denominator support",dt$id[common])
ratio_rows <- rbind(family_only,common_family,common_gross)
fwrite(ratio_rows,file.path(out,"ratio_summary.csv"))

# Weighted wealth/income levels and levels normalized by income at a common group
# scale. Both income-ranked and wealth-ranked groups are shown to make the grouping explicit.
dt[,wealth_q:=rank_bin(nonhousing_nw,weight,5L)]
dt[,income_q:=rank_bin(family_income,weight,5L)]
level_tables <- rbindlist(lapply(c("wealth_q","income_q"),function(g) {
  rbindlist(lapply(1:5,function(k) {
    ix <- dt[[g]] == k
    data.table(group=g,group_id=k,family_years=sum(ix),persons=uniqueN(dt$id[ix]),
      weight_share=sum(dt$weight[ix])/sum(dt$weight),mean_nonhousing_nw=weighted_mean(dt$nonhousing_nw[ix],dt$weight[ix]),
      mean_family_income=weighted_mean(dt$family_income[ix],dt$weight[ix]),
      ratio_of_group_weighted_means=weighted_mean(dt$nonhousing_nw[ix],dt$weight[ix])/weighted_mean(dt$family_income[ix],dt$weight[ix]),
      ratio_to_full_sample_mean_income=weighted_mean(dt$nonhousing_nw[ix],dt$weight[ix])/weighted_mean(dt$family_income,dt$weight))
  }))
}))
fwrite(level_tables,file.path(out,"group_level_normalization.csv"))

rank_stats <- data.table(
  income_measure=c("INCFAMR family income","EARNINDRRC gross labor earnings; all finite values including zeros","EARNINDRRC gross labor earnings > $1,000"),
  spearman=c(weighted_spearman(dt$nonhousing_nw,dt$family_income,dt$weight),
    weighted_spearman(dt$nonhousing_nw,dt$gross_labor,dt$weight),
    weighted_spearman(dt$nonhousing_nw[common],dt$gross_labor[common],dt$weight[common])),
  family_years=c(nrow(dt),sum(is.finite(dt$gross_labor)),sum(common)),
  persons=c(uniqueN(dt$id),uniqueN(dt$id[is.finite(dt$gross_labor)]),uniqueN(dt$id[common]))
)
fwrite(rank_stats,file.path(out,"weighted_rank_correlations.csv"))

dt[common, `:=`(gross_wealth_q=rank_bin(nonhousing_nw,weight,5L),gross_earnings_tercile=rank_bin(gross_labor,weight,3L))]
joint <- dt[common & !is.na(gross_earnings_tercile)&!is.na(gross_wealth_q),.(weighted_mass=sum(weight),family_years=.N,persons=uniqueN(id)),by=.(earnings_tercile=gross_earnings_tercile,wealth_q=gross_wealth_q)]
joint[,share_of_common_sample:=weighted_mass/sum(weighted_mass)]
setorder(joint,earnings_tercile,wealth_q)
fwrite(joint,file.path(out,"joint_earnings_tercile_wealth_quintile.csv"))

# Exact 499-draw person-cluster bootstrap for mean/median row ratios, ratio of
# weighted sums, gross/family candidate on common support, and node means.
ids <- unique(dt$id); pid <- match(dt$id,ids)
bootfun <- function(freq) {
  w <- dt$weight*freq[pid]
  calc <- function(mask,num,den) {
    ok <- mask&is.finite(num)&is.finite(den)&den>0&w>0
    r <- num[ok]/den[ok]; ww <- w[ok]
    c(mean=weighted_mean(r,ww),median=weighted_q(r,ww,.5),ratio_sums=sum(ww*num[ok])/sum(ww*den[ok]))
  }
  c(legacy=calc(rep(TRUE,nrow(dt)),dt$nonhousing_nw,dt$family_income),
    common_family=calc(common,dt$nonhousing_nw,dt$family_income),
    gross=calc(common,dt$nonhousing_nw,dt$gross_labor),
    setNames(vapply(1:5,function(k) {ok<-dt$legacy_bin==k&w>0;weighted_mean(dt$legacy_ratio[ok],w[ok])},numeric(1)),paste0("fixed_original_bin_node",1:5)))
}
set.seed(20260715L)
boot <- matrix(NA_real_,nrow=499,ncol=length(bootfun(rep(1,length(ids)))))
colnames(boot)<-names(bootfun(rep(1,length(ids))))
for (b in 1:499) {
  freq <- tabulate(sample.int(length(ids),length(ids),replace=TRUE),nbins=length(ids))
  boot[b,] <- bootfun(freq)
}
boot_summary <- data.table(statistic=colnames(boot),estimate=as.numeric(bootfun(rep(1,length(ids)))),
  bootstrap_se=apply(boot,2,sd,na.rm=TRUE),bootstrap_reps=499L,seed=20260715L)
fwrite(boot_summary,file.path(out,"bootstrap_summary.csv"))
fwrite(as.data.table(boot),file.path(out,"bootstrap_draws.csv"))

message(sprintf("Done: %d family-years, %d persons; outputs: %s",nrow(dt),uniqueN(dt$id),out))
