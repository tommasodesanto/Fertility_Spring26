#!/usr/bin/env python3
"""First-pass population lifecycle dashboard; saved model only, no solves/targets."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
os.environ['MPLBACKEND']='Agg'
import argparse, csv, gzip, hashlib, json, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from inspect_e5f_saved_households import load
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/model/native_financing_diagnostic_20260919/specification_followup/housing_profiles_v1/full'
FERT=ROOT/'output/model/e5f_matched_pf_20260909a/parameter_target_audit/fertility/fertility_availability.json'
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def cps_profiles(out):
    info=json.loads(FERT.read_text()); rows=[]; checks={}
    assert SHA(info['schema_source'])==info['schema_sha256']
    source=Path(info['cps_source']); compressed=Path(str(source)+'.gz')
    stream=source.open('rb') if source.exists() else gzip.open(compressed,'rb')
    with stream:
        for year in (2004,2006):
            q=info['partitions'][str(year)]; stream.seek(q['byte_start'])
            b=stream.read(q['byte_end_exclusive']-q['byte_start'])
            assert hashlib.sha256(b).hexdigest()==q['partition_sha256']
            checks[str(year)]={'sha256':q['partition_sha256'],'bytes':len(b)}
            for pos in range(0,len(b),261):
                r=b[pos:pos+261]
                assert int(r[:4])==year and int(r[9:11])==6
                if int(r[148:149])!=2: continue
                age,n,w=int(r[146:148]),int(r[238:241]),int(r[250:260])/10000
                if 18<=age<=49 and 0<=n<=20 and w>0: rows.append((year,age,min(n,3),w))
    df=pd.DataFrame(rows,columns=['year','age','children','weight'])
    exact=df[df.age==25]; observed=np.average(exact.children,weights=exact.weight)
    assert abs(observed-.8095276384290021)<1e-12
    df['cell']=18+4*((df.age-18)//4)
    # Ages 46--49 would be incomplete if supplement ends at 44: retain full bins only.
    result=[]
    for a,d in df.groupby('cell'):
        if set(d.age.unique())!=set(range(a,a+4)): continue
        result.append(dict(age_lower=int(a),age_upper=int(a+3),value=np.average(d.children,weights=d.weight),n=len(d),weight=d.weight.sum()))
    return result,dict(source=str(source if source.exists() else compressed),partitions=checks,exact_age25_replay=observed,sample='June 2004/2006 women age18-49, FREVER0-20, FRSUPPWT>0; only complete four-year bins18-41 retained (supplement upper age44 leaves42-45 incomplete)',weight='FRSUPPWT/10000 pooled; children min(FREVER,3)',receipt=str(FERT))

def psid_profiles(out):
    dest=out/'psid_four_year_profiles.csv'
    r=r'''
    suppressPackageStartupMessages({library(haven);library(data.table)})
    setDTthreads(1L)
    a<-commandArgs(TRUE); d<-as.data.table(read_dta(a[1],col_select=c('year','RELTOHEAD_','AGEREP','NETWORTHR','EARNINDRRC','IW')))
    for(k in names(d))set(d,j=k,value=as.numeric(d[[k]]))
    d<-d[year %in% c(2005,2007)&RELTOHEAD_==10&AGEREP>=18&AGEREP<=85&is.finite(IW)&IW>0&is.finite(NETWORTHR)&(AGEREP>65|(is.finite(EARNINDRRC)&EARNINDRRC>=0))]
    denom<-d[AGEREP<=65,sum(IW*EARNINDRRC)/sum(IW)]
    ratio<-d[,sum(IW*NETWORTHR)]/d[AGEREP<=65,sum(IW*EARNINDRRC)]
    stopifnot(nrow(d)==11324,abs(ratio-6.92658379107299)<1e-10)
    d[,cell:=18+4*floor((AGEREP-18)/4)]
    z<-d[,.(n=.N,weight=sum(IW),mean_net_worth=sum(IW*NETWORTHR)/sum(IW),mean_worker_earnings=denom,aggregate_ratio=ratio),by=.(age_lower=cell)]
    z[,age_upper:=age_lower+3];z[,value:=mean_net_worth/mean_worker_earnings];setorder(z,age_lower);fwrite(z,a[2])
    '''
    script=out/'psid_selected_columns.R';script.write_text(r)
    source=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta')
    subprocess.run(['Rscript',str(script),str(source),str(dest)],check=True,timeout=240)
    d=pd.read_csv(dest)
    assert d.n.sum()==11324 and np.max(abs(d.aggregate_ratio-6.92658379107299))<1e-10
    return d.to_dict('records'),dict(source=str(source),source_bytes=source.stat().st_size,source_mtime_ns=source.stat().st_mtime_ns,source_rehashed=False,sample='Pooled2005/2007 PSID reference persons RELTOHEAD_=10 ages18-85; IW positive finite, NETWORTHR finite; working ages18-65 additionally finite nonnegative EARNINDRRC. Matches accepted aggregate sample11324.',weight='IW pooled family-years',scaling='Each age-cell mean NETWORTHR divided by weighted mean EARNINDRRC among selected working-age18-65 households; model scaled analogously.',aggregate_ratio_replay=float(d.aggregate_ratio.iloc[0]),mean_working_earnings=float(d.mean_worker_earnings.iloc[0]))

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    ref=json.loads((ROOT/'output/model/daytime_calibration_20260927/households/source_receipt.json').read_text())
    os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']=ref['contract_sha256'];os.environ['E5F_LOCAL_EXECUTION_AUTHORIZATION']='tommaso_authorized_20260927_local_primary_continuation_v1'
    packet,rt,r=load(Path(ref['contract']),Path(ref['case']),out)
    P=packet['parameters'];e=packet['evaluation'];pol=e.policy;g=e.g_current;bg=packet['b_grid'];model=rt['model'];zv,_,_=model.income_transition_values(P)
    modelrows=[];earnnum=0.;earnmass=0.
    for j in range(P.J):
        a=P.age_start+P.da*j;gj=g[:,:,:,j];m=gj.sum();bm=gj.sum(axis=(1,2,3,4,5));tm=gj.sum(axis=(0,2,3,4,5));zm=gj.sum(axis=(0,1,2,4,5))
        asset=e.g_post_fertility[:,:,:,j];ab=asset.sum(axis=(1,2,3,4,5));at=asset.sum(axis=(0,2,3,4,5))
        housewealth=np.dot(at[1:],np.asarray(P.H_own)*pol.price[0]);nw=(np.dot(ab,bg)+housewealth)/asset.sum()
        room=(gj[:,0]*np.minimum(pol.hR_pol[:,0,:,j],9)).sum()+np.dot(tm[1:],np.minimum(P.H_own,9))
        nc=lambda gg:np.dot(np.arange(gg.shape[4]),gg.sum(axis=(0,1,2,3,5)))/gg.sum()
        children=.5*(nc(e.g_pre[:,:,:,j])+nc(gj))
        modelrows.append(dict(age_lower=a,age_upper=a+3,ownership=tm[1:].sum()/m,rooms=room/m,wealth=nw,children=children,children_post=nc(gj),children_pre=nc(e.g_pre[:,:,:,j]),mass=m))
        if a<=65:
            earnnum+=np.dot(zm,[model.annual_gross_income_at_state(P,0,j,z) for z in zv]);earnmass+=m
    meanearn=earnnum/earnmass
    for r in modelrows:r['wealth']/=meanearn
    early=next(r for r in modelrows if r['age_lower']==22)
    assert abs(.125*early['children_pre']+.875*early['children_post']-.5275288287814585)<1e-11
    aggregate_ratio=sum(r['mass']*r['wealth'] for r in modelrows)/sum(r['mass'] for r in modelrows if r['age_lower']<=65)
    assert abs(aggregate_ratio-6.024670022227826)<1e-10
    pd.DataFrame(modelrows).to_csv(out/'model_profiles.csv',index=False)
    h=pd.read_csv(BASE/'housing_profile_by_age.csv');h=h[(h.geography_scope=='national')&(h.age_kind=='four_year')&(h['sample']=='all_structures')]
    h=h.groupby(['age_lower','age_upper'],as_index=False)[['n_records','hhwt','rooms_capped9_sum','owner_hhwt']].sum()
    cps,cp=cps_profiles(out);wealth,wp=psid_profiles(out)
    rows=[]
    for q in modelrows:
        for metric in ('ownership','rooms','wealth','children'):rows.append(dict(series='Model',metric=metric,age_lower=q['age_lower'],age_upper=q['age_upper'],value=q[metric],n='',weight=q['mass']))
    for _,q in h.iterrows():
        for metric,col in [('ownership','owner_hhwt'),('rooms','rooms_capped9_sum')]:rows.append(dict(series='ACS2005-06',metric=metric,age_lower=q.age_lower,age_upper=q.age_upper,value=q[col]/q.hhwt,n=q.n_records,weight=q.hhwt))
    for series,metric,data in [('PSID2005/07','wealth',wealth),('CPS2004/06','children',cps)]:
        rows += [dict(series=series,metric=metric,**{k:q[k] for k in ('age_lower','age_upper','value','n','weight')}) for q in data]
    data=pd.DataFrame(rows);data.to_csv(out/'lifecycle_comparison.csv',index=False)
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12,8))
    titles={'ownership':'Homeownership','rooms':'Housing size','wealth':'Net worth','children':'Children ever born'}
    labels={'ownership':'Share owning','rooms':'Mean rooms, capped at 9','wealth':'Mean / mean working-age annual earnings','children':'Mean children, capped at 3'}
    for ax,metric in zip(axs.flat,titles):
        d=data[data.metric==metric]
        for series,ss in d.groupby('series',sort=False):ax.plot(ss.age_lower+1.5,ss.value,'o-' if series=='Model' else 's--',markersize=4,color='#28658b' if series=='Model' else '#c26a32',label=series)
        ax.set(title=titles[metric],ylabel=labels[metric],xlabel='Age (four-year cell midpoint)',xlim=(18,85));ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=9)
        if metric=='ownership':ax.set_ylim(0,1)
        if metric=='children':ax.set_xlim(18,49)
    fig.suptitle('Lifecycle profiles: saved calibration versus survey data',fontsize=15)
    fig.text(.5,.017,'Diagnostic only; cross-sections, not three simulated lives. All-structure ACS housing; PSID household wealth; CPS women.\nChildren use uniform within-period birth timing. Lines connect four-year cells; no annual interpolated observations.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.068,1,.95]);fig.savefig(out/'lifecycle_dashboard.png',dpi=160);plt.close(fig)
    prov=dict(status='first_pass_diagnostic_not_new_calibration_targets',model=ref,model_mean_working_age_earnings=float(meanearn),model_measurement='Exact g_current population, no simulation draws; physical hR_pol/H_own capped9; beginning-period net worth b+qH weighted by g_post_fertility before housing transactions, matching the calibrated wealth observer; ownership actual tenure. Children arithmetic midpoint of pre-birth and post-birth conditional means, uniform birth-time approximation over each four-year interval; literal cap3, not3.602 top-bin adjustment.',housing=dict(source=str(BASE/'housing_profile_by_age.csv'),sha256=SHA(BASE/'housing_profile_by_age.csv'),source_provenance=str(BASE/'provenance.json'),sample='National ACS2005/2006 household heads ages18-85, all structures, valid tenure/rooms, positiveHHWT; exact precomputed four-year sufficient-statistic bins.',weight='HHWT pooled',caveat='All structures differs from DUE ownership calibration sample; capped9 ACS rooms differs from uncapped AHS2007 calibration target. These are separate diagnostic comparisons.'),wealth=wp,children=cp,caveats=['Cross-sectional age gradients are not cohort trajectories. No standard errors yet; small old-age PSID cells can make mean wealth noisy.','Stationary household reproductive-member ages are compared with CPS women, not an exact sex-composition reconstruction.','Housing is synchronized post-choice; wealth is beginning-period before housing transactions, matching the target observer. Survey interview timing differs. Wealth valuation and earnings-year metadata inherited from PSID shelf.','Model wealth normalization uses its own working-age mean earnings; empirical uses its own. No dollar-level normalization asserted.','No target, weighting, parameter or model change; no native solve.'])
    (out/'provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),rows=len(data),native_solves=0,model_mean_earnings=meanearn)))
if __name__=='__main__':main()
