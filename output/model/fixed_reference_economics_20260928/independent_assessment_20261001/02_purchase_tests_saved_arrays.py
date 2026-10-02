"""Purchase-test tabulations on SAVED winner31 baseline arrays (chain 7 / case 0173, loss 31.284, q=0.719168; phi=.8/.8).
Read-only arithmetic on saved policy and distribution arrays; no model solves. Run from the repository root."""


# ======================================================================
# Part A: origination LTV and end-of-period position of renter-buyers
# ======================================================================
"""Zero-solve tabulation on SAVED winner31 baseline arrays (phi=.8/.8, q=0.719168): origination LTV of renter-buyers.
Group: never-parent renters (inherited tenure 0, n=0, no child at home). Exact timing for this group:
PRE mass -> first-birth attempt (p1) x fecundity -> branch (wait: n=0,cs=0 ; child: n=1,cs=1) -> tenure/size choice."""
import numpy as np, json
np.set_printoptions(linewidth=200, suppress=True, precision=4)
W='output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/'
B=W+'purchase_ltv_v1/local_run/retry5/results/'
a=np.load(B+'baseline_80_80/solution_arrays.npz',allow_pickle=False)
pre=np.load(B+'q0_reference_inherited_states.npz')['g_pre']
s=json.load(open(W+'diagnosis/child_wait_value_v1/summary.json'))
fec=np.zeros(17); fb=s['fecundity_by_age']
print('fecundity_by_age raw:',fb)
vals=list(fb.values()) if isinstance(fb,dict) else list(fb)
fec[:len(vals)]=vals if len(vals)<=17 else vals[:17]
bg=a['b_grid']; q=float(a['p_eq'][0]); H=np.array([2,4,6,8,10.]); Q=q*H; phi=.8; R=1.02**4
w=pre[:,0,0,:,:,0,0]                     # (b, age, z) never-parent renters, pre-fertility
p1=a['fert_probs'][:,0,0,:,:,1]            # attempt prob
tp=a['tenure_probs'][:,0,0,:,:,:,:,:]      # (b,age,z,n,cs,dest)
pi=fec[None,:,None]
m_child=w*pi*p1; m_wait=w-m_child
ages=18+4*np.arange(17)
def tab(mask_age,label):
    out={}
    for br,(m,n,cs) in {'wait':(m_wait,0,0),'child':(m_child,1,1)}.items():
        mm=m*mask_age[None,:,None]
        for k in range(5):
            buy=mm*tp[:,:,:,n,cs,k+1]                 # purchase mass by (b,age,z)
            ltv=1-bg/Q[k]                              # origination LTV at closing for inherited b
            out[(br,k)]=(buy.sum(),buy[ltv>.8+1e-12].sum(),buy[ltv>.9+1e-12].sum(),buy[ltv>=1-1e-12].sum(),(buy.sum((1,2))*np.minimum(np.maximum(ltv,0),2)).sum())
    tot=sum(v[0] for v in out.values()); print(f'--- {label}: group PRE mass {(w*mask_age[None,:,None]).sum():.5f}, purchase mass {tot:.5f} (rate {tot/(w*mask_age[None,:,None]).sum():.3%})')
    print(' branch rooms  purch.mass  share>80%  share>90%  share=100%+  mean orig LTV')
    for (br,k),v in out.items():
        if v[0]>1e-9: print(f' {br:5s} {int(H[k]):3d}   {v[0]:.6f}   {v[1]/v[0]:.3f}     {v[2]/v[0]:.3f}     {v[3]/v[0]:.3f}       {v[4]/v[0]:.3f}')
    for br in ('wait','child'):
        t=np.array([out[(br,k)] for k in range(5)]).sum(0)
        print(f' {br:5s} all   {t[0]:.6f}   {t[1]/t[0]:.3f}     {t[2]/t[0]:.3f}     {t[3]/t[0]:.3f}       {t[4]/t[0]:.3f}')
    t=np.array(list(out.values())).sum(0); print(f' ALL         {t[0]:.6f}   {t[1]/t[0]:.3f}     {t[2]/t[0]:.3f}     {t[3]/t[0]:.3f}       {t[4]/t[0]:.3f}')
tab((ages<=42).astype(float),'never-parent renters, ages 18-42')
tab(((ages<=30)).astype(float),'never-parent renters, ages 18-30')
tab(np.ones(17),'never-parent renters, all ages')
# wealth of never-parent renters relative to the 20% thresholds
for lab,mk in (('18-30',ages<=30),('18-42',ages<=42)):
    ww=(w*mk[None,:,None]).sum((1,2)); tot=ww.sum()
    print(f'never-parent renters {lab}: share with b<20% of Q for 2/4/6/8 rooms:',[round(float(ww[bg<.2*Qk-1e-12].sum()/tot),3) for Qk in Q[:4]],' share with b=0:',round(float(ww[np.isclose(bg,0)].sum()/tot),3),' mean b',round(float((ww*bg).sum()/tot),3),' median b',float(bg[np.searchsorted(np.cumsum(ww)/tot,.5)]))
# end-of-period position of buyers: b' policy of owner size k evaluated at x=b-Q_k
bp=a['bp_pol']
print('min saved bp_pol by owner size (all states) vs -phi*Q:',[ (round(float(np.nanmin(bp[:,k+1][np.isfinite(bp[:,k+1])])),4), round(-phi*Q[k],4)) for k in range(5)])
def endpos(mask_age,label):
    res=[]
    for br,(m,n,cs) in {'wait':(m_wait,0,0),'child':(m_child,1,1)}.items():
        for k in range(5):
            for j in np.where(mask_age>0)[0]:
                for z in range(9):
                    buy=m[:,j,z]*tp[:,j,z,n,cs,k+1]
                    if buy.sum()<=0: continue
                    x=bg-Q[k]; pol=bp[:,k+1,0,j,z,n,cs]
                    bnext=np.interp(x,bg,pol)
                    for i in np.where(buy>0)[0]: res.append((buy[i],-bnext[i]/Q[k],1-bg[i]/Q[k],k,br=='child'))
    r=np.array(res); wt=r[:,0]/r[:,0].sum()
    e=r[:,1]; o=r[:,2]
    print(f'--- {label}: end-of-purchase-period LTV (-b\'/Q): share at floor (>=0.795) {wt[e>=.795].sum():.3f}; >0.7 {wt[e>.7].sum():.3f}; mean {np.sum(wt*e):.3f}; share with positive b\' {wt[e<0].sum():.3f}')
    hi=o>.8+1e-12
    print(f'    among purchases with origination LTV>80% (share {wt[hi].sum():.3f}): share ending at floor {wt[hi&(e>=.795)].sum()/wt[hi].sum():.3f}; mean end LTV {np.sum(wt[hi]*e[hi])/wt[hi].sum():.3f}; mean orig LTV {np.sum(wt[hi]*o[hi])/wt[hi].sum():.3f}')
endpos((ages<=42).astype(float),'never-parent renter buyers 18-42')


# ======================================================================
# Part B: budget-identity check of the array reading; recovered income by age and type
# ======================================================================
import numpy as np, json
np.set_printoptions(linewidth=220, suppress=True, precision=4)
W='output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/'
B=W+'purchase_ltv_v1/local_run/retry5/results/'
a=np.load(B+'baseline_80_80/solution_arrays.npz',allow_pickle=False)
bg=a['b_grid']; q=float(a['p_eq'][0]); H=np.array([2,4,6,8,10.]); Q=q*H; R=1.02**4
kq=0.05545379079326218+0.042393443095490375
bp=a['bp_pol']; cp=a['c_pol']; tp=a['tenure_probs']
# implied income from owner-block budget: y = c + K + b' - R*x   (x = wealth index of the owner block)
print('Implied y = c + K + b\' - R*x in OWNER block (size 6 rooms, n=0,cs=0), by age j and income z; spread across wealth nodes should be ~0 if my reading is right')
k=2
for j in (0,1,2,3):
    row=[]
    for z in range(9):
        x=bg; c=cp[:,k+1,0,j,z,0,0]; b1=bp[:,k+1,0,j,z,0,0]
        ok=np.isfinite(c)&np.isfinite(b1)&(c>0)&(x>-Q[k]-1e-9)&(x<3)
        y=c[ok]+kq*Q[k]+b1[ok]-R*x[ok]
        row.append((round(float(np.median(y)),3),round(float(y.max()-y.min()),3),int(ok.sum())))
    print(' age',18+4*j,row)
# same from renter block: y = c + rent*hR + b' - R*b ; rent unknown -> solve using two unknown? print c + b' - R b by hR to see
print('Renter block: c + b\' - R*b (should equal y - rent*hR):')
hr=a['hR_pol']
for j in (0,1):
    for z in (3,4,5):
        ok=(bg>=0)&(bg<2)&np.isfinite(cp[:,0,0,j,z,0,0])&(cp[:,0,0,j,z,0,0]>0)
        v=cp[ok,0,0,j,z,0,0]+bp[ok,0,0,j,z,0,0]-R*bg[ok]
        print(' age',18+4*j,'z',z,'c+b\'-Rb:',np.round(v[:6],3),'hR:',np.round(hr[ok,0,0,j,z,0,0][:6],2),'c:',np.round(cp[ok,0,0,j,z,0,0][:6],3),"b':",np.round(bp[ok,0,0,j,z,0,0][:6],3))
# concrete buyer example: age 22 (j=1), z=5, wait branch, 6 rooms
j,z=1,5
print('Example j=1 (age 22), z=5, wait branch (n=0,cs=0): b, P(buy 6r), x=b-Q, b\' as buyer, end LTV, c as buyer')
for i in np.where((bg>=0)&(bg<1.8))[0]:
    x=bg[i]-Q[k]; b1=np.interp(x,bg,bp[:,k+1,0,j,z,0,0]); c1=np.interp(x,bg,cp[:,k+1,0,j,z,0,0])
    print(f'  b={bg[i]:.3f} Pbuy6={tp[i,0,0,j,z,0,0,k+1]:.4f} Pbuy4={tp[i,0,0,j,z,0,0,2]:.4f} Prent={tp[i,0,0,j,z,0,0,0]:.4f} x={x:.3f} b\'={b1:.3f} endLTV={-b1/Q[k]:.3f} c={c1:.3f}')


# ======================================================================
# Part C: pass rates of alternative purchase tests
# ======================================================================
"""Zero-solve check on saved winner31 baseline arrays: is the income-inclusive eligibility screen ever the binding restriction?
Implied net income y(age,z) recovered from the owner-block budget identity for z>=4 (exact, zero spread) and scaled by type value for z<4."""
import numpy as np, json
np.set_printoptions(linewidth=220, suppress=True, precision=4)
W='output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/'
B=W+'purchase_ltv_v1/local_run/retry5/results/'
a=np.load(B+'baseline_80_80/solution_arrays.npz',allow_pickle=False)
pre=np.load(B+'q0_reference_inherited_states.npz')['g_pre']
s=json.load(open(W+'diagnosis/child_wait_value_v1/summary.json')); fec=np.array(s['fecundity_by_age'])
bg=a['b_grid']; q=float(a['p_eq'][0]); H=np.array([2,4,6,8,10.]); Q=q*H; R=1.02**4; phi=.8
kq=0.05545379079326218+0.042393443095490375
bp=a['bp_pol']; cp=a['c_pol']; tp=a['tenure_probs'][:,0,0]; tv=a['type_values']
# implied income by (age,z): use z=4..8 exact, infer profile factor f_j = y/tv
k=2; y=np.zeros((17,9))
for j in range(17):
    f=[]
    for z in range(4,9):
        c=cp[:,k+1,0,j,z,0,0]; b1=bp[:,k+1,0,j,z,0,0]; ok=np.isfinite(c)&(c>0)&(bg>-Q[k])&(bg<3)
        if ok.sum()>5: f.append(np.median(c[ok]+kq*Q[k]+b1[ok]-R*bg[ok])/tv[z])
    y[j]=np.median(f)*tv if f else np.nan
print('net four-year income per unit of type value, by age cell:',np.round(y[:,4]/tv[4],4))
print('y(age 18-25) by type:',np.round(y[0],3)); print('y(age 26-33) by type:',np.round(y[2],3))
w=pre[:,0,0,:,:,0,0]; p1=a['fert_probs'][:,0,0,:,:,1]; pi=fec[None,:,None]
m_child=w*pi*p1; m_wait=w-m_child; ages=18+4*np.arange(17); young=(ages<=42)
tot=0; scr_fail=0; lam25_fail=0; due_fail=0; endfeas_fail=0; rows=[]
for br,(m,n,cs) in {'wait':(m_wait,0,0),'child':(m_child,1,1)}.items():
    for kk in range(5):
        for j in np.where(young)[0]:
            for z in range(9):
                buy=m[:,j,z]*tp[:,j,z,n,cs,kk+1]
                if buy.sum()<=0: continue
                b=bg; yy=y[j,z]
                screen=b+yy/R>=(1-phi)*Q[kk]-1e-12
                lam25=b+.25*yy/R>=(1-phi)*Q[kk]-1e-12
                due=b>=(1-phi)*Q[kk]-1e-12
                cmax=R*(b-Q[kk])+yy-kq*Q[kk]+phi*Q[kk]      # max consumption given end floor (ignoring transfers)
                tot+=buy.sum(); scr_fail+=buy[~screen].sum(); lam25_fail+=buy[~lam25].sum(); due_fail+=buy[~due].sum(); endfeas_fail+=buy[cmax<=0].sum()
print(f'never-parent renter buyers 18-42: purchase mass {tot:.5f}')
print(f'  share of purchase mass failing implemented screen b+y/R>=0.2Q: {scr_fail/tot:.4f}')
print(f'  share with max consumption<=0 at the end floor (no-transfer income): {endfeas_fail/tot:.4f}')
print(f'  share that would fail a one-year-income screen b+0.25*y/R>=0.2Q: {lam25_fail/tot:.4f}')
print(f'  share that would fail the stock-only (DUE) screen b>=0.2Q: {due_fail/tot:.4f}')
# incidence among ALL never-parent renters 18-42 (not just buyers): for each house size, share of mass that passes each test
ww=w*young[None,:,None]; T=ww.sum()
print('Share of never-parent renters 18-42 (PRE mass) passing each test, by house size  [implemented screen | end-floor+positive c (c>0.1*4yr-income? no: c>0) | one-year-income screen | stock-only]')
for kk in range(5):
    sc=np.zeros_like(ww,bool); ef=sc.copy(); l25=sc.copy(); du=sc.copy(); ef30=sc.copy()
    for j in np.where(young)[0]:
        for z in range(9):
            yy=y[j,z]; b=bg
            sc[:,j,z]=b+yy/R>=(1-phi)*Q[kk]; l25[:,j,z]=b+.25*yy/R>=(1-phi)*Q[kk]; du[:,j,z]=b>=(1-phi)*Q[kk]
            cmax=R*(b-Q[kk])+yy-kq*Q[kk]+phi*Q[kk]; ef[:,j,z]=cmax>0; ef30[:,j,z]=cmax>=0.5*yy
    print(f'  {int(H[kk]):2d} rooms (Q={Q[kk]:.2f}): screen {ww[sc].sum()/T:.3f} | end-floor c>0 {ww[ef].sum()/T:.3f} | end-floor with c>=half of income {ww[ef30].sum()/T:.3f} | one-year {ww[l25].sum()/T:.3f} | stock-only {ww[du].sum()/T:.3f}')
