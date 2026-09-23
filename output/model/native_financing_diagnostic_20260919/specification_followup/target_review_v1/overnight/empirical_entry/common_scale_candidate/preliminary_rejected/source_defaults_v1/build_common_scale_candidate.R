#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(haven); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
repo <- normalizePath('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
out <- file.path(repo,'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry/common_scale_candidate')
psid <- if(length(args)) args[[1]] else file.path(dirname(repo),'PSID/PSIDSHELF_MOBILITY.dta')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
wm <- function(x,w){ok<-is.finite(x)&is.finite(w)&w>0;if(!any(ok))return(NA_real_);sum(x[ok]*w[ok])/sum(w[ok])}
weighted_midrank_bin <- function(x,w,n){
  # All observations tied at exactly x share one bin; midrank determines bin.
  ok<-is.finite(x)&is.finite(w)&w>0; ans<-rep(NA_integer_,length(x)); ix<-which(ok)
  groups<-data.table(ix=ix,x=x[ix],w=w[ix])[,.(mass=sum(w)),by=x][order(x)]
  groups[,mid:=(cumsum(mass)-.5*mass)/sum(mass)]
  groups[,bin:=pmin(n,floor(mid*n)+1L)]
  ans[ix]<-groups$bin[match(x[ix],groups$x)]; ans
}
message('Read exact July sample and wave-normalization fields...')
d<-as.data.table(read_dta(psid,col_select=c('ID','year','AGEREP','DEATHYEAR','RELTOHEAD_','RELCHIREP','INCFAMR','NETWORTH2R','NETWORTHR','HOMEOWN','IW','EARNINDRRC')))
dt<-d[,.(id=as.numeric(ID),year=as.integer(year),age=as.numeric(AGEREP),death_year=as.numeric(DEATHYEAR),rp=as.numeric(RELTOHEAD_)==10,children=as.numeric(RELCHIREP),family_income=as.numeric(INCFAMR),nw2=as.numeric(NETWORTH2R),nw_builder=as.numeric(NETWORTHR),homeown=as.numeric(HOMEOWN),w=as.numeric(IW),earn=as.numeric(EARNINDRRC))]
rm(d); invisible(gc())
base<-dt[year>=1984 & year<=2019 & (is.na(death_year)|year<=death_year) & is.finite(w)&w>0]
july<-base[rp & is.finite(age)&age>=18&age<=24&children==0&homeown==2&is.finite(family_income)&family_income>1000&is.finite(nw2)]
if(nrow(july)!=1835L||uniqueN(july$id)!=1346L)stop(sprintf('Exact July sample mismatch: %d rows, %d persons',nrow(july),uniqueN(july$id)))
# July wave mean applies aggregate builder's RP, age, alive, positive-IW, finite NETWORTHR universe;
# worker denominator further requires age<=65 and finite nonnegative gross earnings. Date range extends
# only as necessary to normalize all old July entrant waves (builder's measured aggregate range starts 2005).
wp<-base[rp & is.finite(age)&age>=18&age<=85&is.finite(nw_builder)]
workers<-wp[age<=65 & is.finite(earn)&earn>=0]
wave<-workers[,.(wave_working_rows=.N,wave_working_persons=uniqueN(id),wave_IW=sum(w),gross_earnings_IW=sum(w*earn),Ebar_wave=wm(earn,w),zero_earnings_rows=sum(earn==0),positive_earnings_rows=sum(earn>0)),by=year][order(year)]
# The authoritative aggregate builder does not apply DEATHYEAR; calculate its exact
# 2005-2019 worker universe as an explicit source-support comparison. Candidate uses alive rows.
builder_pop<-dt[year>=2005&year<=2019&rp&is.finite(age)&age>=18&age<=85&is.finite(w)&w>0&is.finite(nw_builder)]
builder_worker<-builder_pop[age<=65&is.finite(earn)&earn>=0]
builder_wave<-builder_worker[,.(builder_working_rows=.N,builder_IW=sum(w),builder_gross_earnings_IW=sum(w*earn),builder_Ebar_no_alive=wm(earn,w)),by=year]
wave<-merge(wave,builder_wave,by='year',all.x=TRUE,sort=FALSE)
wave[,Ebar_alive_minus_builder_no_alive:=Ebar_wave-builder_Ebar_no_alive]
wave[,date_support:=fifelse(year>=2005,'within aggregate builder 2005-2019','extended earlier to cover July entrant waves')]
july<-merge(july,wave[,.(year,Ebar_wave,wave_working_rows,wave_working_persons,wave_IW,date_support)],by='year',all.x=TRUE,sort=FALSE)
if(any(!is.finite(july$Ebar_wave)|july$Ebar_wave<=0))stop('A July entrant wave lacks positive normalization Ebar.')
july[,omega:=nw2/Ebar_wave]
if(any(!is.finite(july$earn))) stop('Missing current gross earnings; cannot define all-row earnings terciles. Report missing separately before revising scope.')
july[,earn_tercile:=weighted_midrank_bin(earn,w,3L)]
july[,wealth_quintile:=NA_integer_]
for(t in 1:3){ix<-which(july$earn_tercile==t); july$wealth_quintile[ix]<-weighted_midrank_bin(july$omega[ix],july$w[ix],5L)}
if(anyNA(july$earn_tercile)|anyNA(july$wealth_quintile))stop('Unexpected unbinned entrant row.')
cell<-july[,.(family_years=.N,persons=uniqueN(id),weighted_mass=sum(w),probability=sum(w)/sum(july$w),conditional_mean_omega=wm(omega,w),mean_nw2=wm(nw2,w),mean_gross_earnings=wm(earn,w),min_gross_earnings=min(earn),max_gross_earnings=max(earn)),by=.(earnings_tercile=earn_tercile,wealth_quintile)][order(earnings_tercile,wealth_quintile)]
cell[,wealth_bin_measure:='omega = NETWORTH2R / same-wave Ebar_wave']
fullgrid<-CJ(earnings_tercile=1:3,wealth_quintile=1:5); cell<-merge(fullgrid,cell,by=c('earnings_tercile','wealth_quintile'),all.x=TRUE,sort=TRUE);cell[is.na(weighted_mass),`:=`(family_years=0L,persons=0L,weighted_mass=0,probability=0,conditional_mean_omega=NA_real_,mean_nw2=NA_real_,mean_gross_earnings=NA_real_,min_gross_earnings=NA_real_,max_gross_earnings=NA_real_)]
cell[,conditional_prob_within_tercile:=probability/sum(probability),by=earnings_tercile]
# Rank-overlap audit: compare the model's product-within-tercile coupling with a
# finer empirical rank coupling that retains the observed wealth-bin mix inside each tied earnings value.
b15_path<-file.path(repo,'output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_entry/model_entry_mapping.json')
b15<-fromJSON(b15_path); zprob<-as.numeric(b15$B15$z_weights); zedge<-c(0,cumsum(zprob)); total_w<-sum(july$w)
tie_groups<-july[,.(mass=sum(w),earn_tercile=earn_tercile[1]),by=earn][order(earn)]
tie_groups[,upper:=cumsum(mass)/total_w];tie_groups[,lower:=shift(upper,fill=0)]
emp_parts<-list(); kk<-0L
for(g in seq_len(nrow(tie_groups))){
 gg<-tie_groups[g]; qmass<-july[earn==gg$earn,.(qmass=sum(w)),by=wealth_quintile]
 for(z in seq_along(zprob)){
  ov<-max(0,min(gg$upper,zedge[z+1])-max(gg$lower,zedge[z]))
  if(ov>0) for(r in seq_len(nrow(qmass))){kk<-kk+1L;emp_parts[[kk]]<-data.table(earnings_tercile=gg$earn_tercile,wealth_quintile=qmass$wealth_quintile[r],z_state=z,empirical_rank_overlap_probability=ov*qmass$qmass[r]/gg$mass)}
 }
}
emp_support<-rbindlist(emp_parts)[,.(empirical_rank_overlap_probability=sum(empirical_rank_overlap_probability)),by=.(earnings_tercile,wealth_quintile,z_state)]
terc_mass<-july[,.(terc_probability=sum(w)/total_w),by=earn_tercile]
model_parts<-list();kk<-0L
for(t in 1:3) for(q in 1:5) for(z in seq_along(zprob)){
 ov<-max(0,min(sum(terc_mass$terc_probability[terc_mass$earn_tercile<=t]),zedge[z+1])-max(sum(terc_mass$terc_probability[terc_mass$earn_tercile<t]),zedge[z]))
 kk<-kk+1L;model_parts[[kk]]<-data.table(earnings_tercile=t,wealth_quintile=q,z_state=z,candidate_product_probability=cell[earnings_tercile==t&wealth_quintile==q,conditional_prob_within_tercile]*ov)
}
model_support<-rbindlist(model_parts);support_cmp<-merge(model_support,emp_support,by=c('earnings_tercile','wealth_quintile','z_state'),all=TRUE);support_cmp[is.na(candidate_product_probability),candidate_product_probability:=0];support_cmp[is.na(empirical_rank_overlap_probability),empirical_rank_overlap_probability:=0];support_cmp[,candidate_combination_absent_in_ranked_microdata:=candidate_product_probability>1e-14&empirical_rank_overlap_probability<=1e-14]
fwrite(support_cmp,file.path(out,'discretization_support.csv'))
# Before discretization, apply the candidate common-scale, tax and return formula directly to each July row.
R_period<-1.04^4; tau_pay<-0.179
july[,direct_resource_common_units:=R_period*omega+4*(1-tau_pay)*earn/Ebar_wave]
direct_negative_mass<-sum(july$w[july$direct_resource_common_units<=0])/total_w
terc<-july[,.(family_years=.N,persons=uniqueN(id),probability=sum(w)/sum(july$w),mean_earnings=wm(earn,w),mean_omega=wm(omega,w),wealth_sd_omega=sqrt(wm((omega-wm(omega,w))^2,w))),by=earn_tercile][order(earn_tercile)]
# weighted Spearman with exact-value pooled midpoint ranks.
wspear<-function(x,y,w){ok<-is.finite(x)&is.finite(y)&is.finite(w)&w>0; rr<-function(a,ww){g<-data.table(a=a,ww=ww)[,. (m=sum(ww)),by=a][order(a)];g[,r:=(cumsum(m)-.5*m)/sum(m)];g$r[match(a,g$a)]};rx<-rr(x[ok],w[ok]);ry<-rr(y[ok],w[ok]); ww<-w[ok]; mx<-wm(rx,ww);my<-wm(ry,ww);sum(ww*(rx-mx)*(ry-my))/sqrt(sum(ww*(rx-mx)^2)*sum(ww*(ry-my)^2))}
# Binning-loss accounting: replace every July omega by its 3x5 conditional mean.
july<-merge(july,cell[,.(earn_tercile=earnings_tercile,wealth_quintile,node=conditional_mean_omega)],by=c('earn_tercile','wealth_quintile'),all.x=TRUE,sort=FALSE)
mu<-wm(july$omega,july$w); variance<-wm((july$omega-mu)^2,july$w); mse<-wm((july$omega-july$node)^2,july$w);mae<-wm(abs(july$omega-july$node),july$w)
metrics<-data.table(metric=c('July family-years','July persons','July IW mass','finite gross-earnings rows','zero gross-earnings rows','positive gross earnings <=1000 rows','waves represented','minimum Ebar_wave','maximum Ebar_wave','mean omega (weighted)','variance omega (weighted)','3x5 node approximation MSE','3x5 node approximation RMSE','3x5 node approximation MAE','variance share left within 3x5 cells (MSE/variance)','variance share explained by node means (1-MSE/variance)','raw weighted Spearman omega/gross earnings','node-imputed weighted Spearman omega/gross earnings','weighted mean absolute omega scale','pooled ratio of weighted sums NETWORTH2R/EARNINDRRC'),value=c(nrow(july),uniqueN(july$id),sum(july$w),sum(is.finite(july$earn)),sum(july$earn==0),sum(july$earn>0&july$earn<=1000),uniqueN(july$year),min(july$Ebar_wave),max(july$Ebar_wave),mu,variance,mse,sqrt(mse),mae,mse/variance,1-mse/variance,wspear(july$omega,july$earn,july$w),wspear(july$node,july$earn,july$w),wm(abs(july$omega),july$w),sum(july$w*july$nw2)/sum(july$w*july$earn)))
metrics<-rbind(metrics,data.table(metric=c('cash-only direct July R*omega + 4*(1-tau)*earn/Ebar_wave <= 0 (descriptive; withdrawn as feasibility)','candidate product mass on (tercile,omega-quintile,z) combinations absent in rank-coupled microdata'),value=c(direct_negative_mass,sum(support_cmp$candidate_product_probability[support_cmp$candidate_combination_absent_in_ranked_microdata]))))
node_mean<-sum(cell$probability*ifelse(is.finite(cell$conditional_mean_omega),cell$conditional_mean_omega,0))
if(abs(sum(cell$probability)-1)>1e-12||abs(node_mean-mu)>1e-12)stop(sprintf('3x5 node mean invariant failed: node %.16g vs entrant E[omega] %.16g',node_mean,mu))
metrics<-rbind(metrics,data.table(metric='weighted 3x5 candidate-node mean omega',value=node_mean))
fwrite(wave,file.path(out,'wave_normalizers.csv')); fwrite(terc,file.path(out,'earnings_terciles.csv')); fwrite(cell,file.path(out,'entry_nodes_3x5.csv')); fwrite(metrics,file.path(out,'candidate_metrics.csv'))
source_counts<-data.table(item=c('July wave rows in aggregate builder 2005-2019','July wave rows in extended 1984-2004 support','unique July waves 1984-2019','alive working denominator rows ages18-65','alive denominator rows within builder dates 2005-2019','alive denominator rows in extended earlier waves','builder exact 2005-2019 working rows without death-year filter','builder wealth-universe rows 2005-2019','candidate sample rows excluded for own earnings missing'),value=c(sum(july$year>=2005),sum(july$year<2005),uniqueN(july$year),nrow(workers),sum(workers$year>=2005),sum(workers$year<2005),nrow(builder_worker),sum(wp$year>=2005),sum(!is.finite(july$earn))))
fwrite(source_counts,file.path(out,'support_counts.csv'))
message(sprintf('Built %d entrants across %d waves; omega mean %.9f; 3x5 MSE %.9f',nrow(july),uniqueN(july$year),mu,mse))
