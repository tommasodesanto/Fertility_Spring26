"""Zero-solve tabulations on saved winner31 baseline arrays: (1) arithmetic bound for the age-25 children-ever-born observer,
(2) who has first births (income type, wealth, tenure, rooms before birth)."""
import numpy as np, json
np.set_printoptions(linewidth=220, suppress=True, precision=4)
W='output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/'
B=W+'purchase_ltv_v1/local_run/retry5/results/'
a=np.load(B+'baseline_80_80/solution_arrays.npz',allow_pickle=False)
pre=np.load(B+'q0_reference_inherited_states.npz')['g_pre']     # (b,ten,loc,age,z,n,cs) pre-fertility
s=json.load(open(W+'diagnosis/child_wait_value_v1/summary.json')); fec=np.array(s['fecundity_by_age'])
tv=a['type_values']; bg=a['b_grid']; H=np.array([0,2,4,6,8,10.])
# ---------- (1) age-25 observer bound -------------
cell=1
pre_par=pre[:,:,:,cell].sum(axis=(0,1,2,3,5))   # by n
tot=pre_par.sum(); f0=pre_par[1]/tot
p1=a['fert_probs'][...,1]                       # (b,ten,loc,age,z) attempt prob for childless
child0=pre[:,:,:,:,:,0,0]                       # childless, settled state
fb=(child0*p1*fec[None,None,None,:,None]).sum(axis=(0,1,2,4))   # first-birth flow by age cell
at_risk=child0.sum(axis=(0,1,2,4))
h1=fb[cell]/pre[:,:,:,cell,:,0,:].sum()
print('cell [22,26): pre-fertility shares by children ever born (0,1,2,3+)',np.round(pre_par/tot,4),' first-birth hazard among childless in this cell h1=%.4f'%h1)
print('first-birth flow by age cell / total women per cell:',np.round(fb/pre.sum(axis=(0,1,2,4,5,6)),4)[:7])
par2_model=0.0888
for sgm in (None,1.0):
    # implied second-birth share: use observer identity  share2 = 0.875*f0*s
    sv=par2_model/(0.875*f0) if sgm is None else sgm
    share0=(1-f0-pre_par[2:].sum()/tot)*(1-0.875*h1)
    ceb=0.125*(f0)+0.875*(f0*(1-sv)+ (1-f0)*h1 + 2*f0*sv)
    print(f'  second-birth hazard in cell 22 for cell-18 mothers s={sv:.3f}: children ever born at 25.5 = {ceb:.4f}; mothers share = {f0+0.875*(1-f0)*h1:.4f}')
print('  => needed f0 for 0.8095 with s=1 and data mothers share 0.4573: f0 =',round((0.8095-0.4573)/0.875,4),'(model f0 = %.4f)'%f0)
# ---------- (2) who has first births -------------
ages=18+4*np.arange(17)
flow=child0*p1*fec[None,None,None,:,None]        # (b,ten,loc,age,z)
for jj,lab in ((slice(0,2),'ages 18-25'),(slice(0,7),'ages 18-45')):
    fz=flow[:,:,:,jj].sum(axis=(0,1,2,3)); rz=child0[:,:,:,jj].sum(axis=(0,1,2,3))
    print(f'{lab}: first-birth probability per period by income type (type value: prob):',' '.join(f'{t:.2f}:{p:.3f}' for t,p in zip(tv,fz/np.maximum(rz,1e-300))))
    print(f'   share of first births by type: {np.round(fz/fz.sum(),3)}   share of at-risk childless by type: {np.round(rz/rz.sum(),3)}')
    print(f'   mean type value: first-birth households {np.sum(fz*tv)/fz.sum():.3f} vs all at-risk childless {np.sum(rz*tv)/rz.sum():.3f}')
    fb_b=flow[:,:,:,jj].sum(axis=(1,2,3,4)); rb=child0[:,:,:,jj].sum(axis=(1,2,3,4))
    print(f'   mean liquid wealth b: first-birth households {np.sum(fb_b*bg)/fb_b.sum():.3f} vs all at-risk {np.sum(rb*bg)/rb.sum():.3f};  inherited owner share: first-birth {flow[:,1:,:,jj].sum()/flow[:,:,:,jj].sum():.3f} vs at-risk {child0[:,1:,:,jj].sum()/child0[:,:,:,jj].sum():.3f}')
# rooms chosen in the wait branch by first-birth-weighted households (same period), vs child branch
tp=a['tenure_probs']; hr=a['hR_pol']
def rooms(n,cs,wgt):
    # expected rooms after tenure choice given post-fertility family state (n,cs); renters: hR_pol at (b,0,..) ; owners moving: H; stayers keep H
    tot=0; m=0
    for ten in range(6):
        for dest in range(6):
            pr=tp[:,ten,0,:,:,n,cs,dest]
            if dest==0:
                # renter housing policy evaluated at post-transaction wealth: for inherited renter x=b ; for seller approximate with own-node policy
                r=hr[:,0,0,:,:,n,cs]
            else: r=H[dest]
            tot+=(wgt[:,ten,0]*pr*r).sum(); m+=(wgt[:,ten,0]*pr).sum()
    return tot/m
for jj,lab in ((slice(0,2),'ages 18-25'),(slice(0,7),'ages 18-45')):
    w=np.zeros_like(flow); w[:,:,:,jj]=flow[:,:,:,jj]
    print(f'{lab}: first-birth-weighted households: rooms this period if birth (n=1) {rooms(1,1,w):.3f} vs if no birth (n=0) {rooms(0,0,w):.3f}  -> same-period difference {rooms(1,1,w)-rooms(0,0,w):.3f}')
    w2=np.zeros_like(flow); w2[:,:,:,jj]=child0[:,:,:,jj]
    print(f'   all at-risk childless: rooms if no birth {rooms(0,0,w2):.3f}')

# ---------- (3) rental-cap truncation among new parents (first-birth-flow weights, all inherited tenures) -------------
cap=6.0
fl=flow[:,:,:,0:7]                                   # first-birth flow, fertile cells
num_rent=0.; num_cap=0.; tot=fl.sum(); own=0.; rr=0.
for ten in range(6):
    pr_rent=tp[:,ten,0,0:7,:,1,1,0]                  # P(rent | birth) for inherited tenure `ten`
    hrc=hr[:,0,0,0:7,:,1,1]                          # renter rooms policy in the birth state (evaluated at own wealth node)
    m=fl[:,ten,0]*pr_rent
    num_rent+=m.sum(); num_cap+=m[hrc>=cap-1e-6].sum(); rr+=(m*hrc).sum()
    own+=(fl[:,ten,0]*tp[:,ten,0,0:7,:,1,1,1:].sum(-1)).sum()
print(f'new parents (first-birth flow): share renting after the birth {num_rent/tot:.3f}, owning {own/tot:.3f}; among those renting: share at the 6-room cap {num_cap/num_rent:.3f}, mean rented rooms {rr/num_rent:.3f}')
by=[(fl[:,:,0]*tp[:,:,0,0:7,:,1,1,k]).sum()/tot for k in range(1,6)]
print('   owner size shares among new parents (2,4,6,8,10 rooms):',np.round(by,4))
# same households if no birth
nr=0.; nc=0.; r0=0.
for ten in range(6):
    m=fl[:,ten,0]*tp[:,ten,0,0:7,:,0,0,0]; hrc=hr[:,0,0,0:7,:,0,0]
    nr+=m.sum(); nc+=m[hrc>=cap-1e-6].sum(); r0+=(m*hrc).sum()
print(f'same households if no birth: share renting {nr/tot:.3f}; among renters share at cap {nc/nr:.3f}, mean rented rooms {r0/nr:.3f}')
