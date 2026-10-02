"""Local Jacobian and attainable-set analysis from the 1,277 passed candidates of normalized_calibration_v2.
Read-only: regressions on an existing CSV; no model solves. Run from the repository root with code/model/.venv/bin/python."""


# ======================================================================
# Part A: regression Jacobian, bounds, bounded linear least squares
# ======================================================================
import numpy as np, pandas as pd
np.set_printoptions(linewidth=220, suppress=True, precision=4)
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
f='output/model/fixed_reference_economics_20260928/normalized_calibration_v2/deployment/calibration_plateau_diagnosis_v1/all_passed_candidates_10_parameters_10_moments.csv'
d=pd.read_csv(f)
P=[c for c in d.columns if c.startswith('param__')]
M=[c[:-7] for c in d.columns if c.endswith('__model')]
print(len(d),'candidates; params',[p[7:] for p in P]); print('moments',M)
print(d[['base_loss','price','H0_derived']+P].describe().T[['min','50%','max','std']])
tg={m:d[m+'__target'].iloc[0] for m in M}; w={m:d[m+'__weight'].iloc[0] for m in M}
print(pd.DataFrame({'target':tg,'weight':w,'min':{m:d[m+'__model'].min() for m in M},'max':{m:d[m+'__model'].max() for m in M},'std':{m:d[m+'__model'].std() for m in M}}))
# correlation of the three big moments
big=['mean_rooms','first_birth_rooms','early_fertility','recent_parent_ownership','ownership_30_55']
print(d[[m+'__model' for m in big]+['price']].corr().round(3))
best=d.loc[d.base_loss.idxmin()]
X=d[P].values; x0=best[P].values.astype(float)
# local regression Jacobian using all points (weights: none), centered at best
sc=X.std(0)
Z=(X-x0)/sc
A=np.column_stack([np.ones(len(d)),Z])
J=np.zeros((len(M),len(P))); R2={}; res_sd={}
for i,m in enumerate(M):
    y=d[m+'__model'].values
    co,res,rk,sv=np.linalg.lstsq(A,y,rcond=None)
    J[i]=co[1:]/sc
    fit=A@co; R2[m]=1-((y-fit)**2).sum()/((y-y.mean())**2).sum(); res_sd[m]=(y-fit).std()
print('design rank',np.linalg.matrix_rank(A),'cond of Z',np.linalg.cond(Z))
print('R2',{k:round(v,3) for k,v in R2.items()})
print('resid sd',{k:float('%.3g'%v) for k,v in res_sd.items()})
Jd=pd.DataFrame(J,index=M,columns=[p[7:] for p in P])
print('Jacobian (d moment / d param), physical units'); print(Jd.round(3))
# elasticity-like: effect of moving each param by 1 cloud-sd
print('effect of +1 cloud sd of param, in weighted-gap units sqrt(w)*dm'); 
W=np.array([w[m] for m in M]); g=np.array([best[m+'__gap'] for m in M],float)
print((Jd.mul(sc,axis=1).mul(np.sqrt(W),axis=0)).round(3))
print('weighted gaps sqrt(w)*g:',dict(zip(M,(np.sqrt(W)*g).round(3))),'loss',(W*g*g).sum())
# SVD of scaled weighted Jacobian
Jw=np.sqrt(W)[:,None]*J*sc[None,:]
U,s,Vt=np.linalg.svd(Jw)
print('singular values (weighted, per cloud-sd):',s.round(4),'cond',s[0]/s[-1])
print('weakest right-singular vector (param mix, cloud-sd units):',dict(zip([p[7:] for p in P],Vt[-1].round(3))))
# unconstrained GN step
step=-np.linalg.solve(Jw,np.sqrt(W)*g)   # in cloud-sd units
print('unconstrained GN step in cloud-sd units:',dict(zip([p[7:] for p in P],step.round(1))))
print('implied physical new params:',dict(zip([p[7:] for p in P],(x0+step*sc).round(4))))
# bounded linear LS
from scipy.optimize import lsq_linear
lb=np.array([0.94,0.1,0.0,0.0,0.1,0.02,0.02,0.01,0.001,0.0]); ub=np.array([0.99,5,0.8,8,2.3,50,50,0.5,0.1,8])
names=[p[7:] for p in P]
order=['beta_annual','chi','child_benefit_curvature','first_birth_fixed_cost','h_P','kappa_fert','kappa_fert_continuation','psi_child','tenure_choice_kappa','theta0']
assert names==order, names
r=lsq_linear(Jw,-np.sqrt(W)*g,bounds=((lb-x0)/sc,(ub-x0)/sc))
print('bounded linear LS predicted loss:',2*r.cost,'step(sd):',r.x.round(1)); print('new params',dict(zip(names,(x0+r.x*sc).round(4))))
pred=g+J@(r.x*sc); print('predicted gaps at bounded optimum',dict(zip(M,pred.round(4))))
# trust-region-limited versions: restrict step to within k cloud-sd
for k in (3,10,30):
    r=lsq_linear(Jw,-np.sqrt(W)*g,bounds=(np.maximum((lb-x0)/sc,-k),np.minimum((ub-x0)/sc,k)))
    pred=g+J@(r.x*sc)
    print(f'|step|<={k} sd: predicted loss {2*r.cost:.3f}; gaps rooms {pred[M.index("mean_rooms")]:.3f}, fbrooms {pred[M.index("first_birth_rooms")]:.3f}, early {pred[M.index("early_fertility")]:.3f}; active h_P at ub: {abs(x0[4]+r.x[4]*sc[4]-2.3)<1e-9}')
# same but with h_P upper bound removed
for k in (3,10,30):
    ub2=ub.copy(); ub2[4]=10
    r=lsq_linear(Jw,-np.sqrt(W)*g,bounds=(np.maximum((lb-x0)/sc,-k),np.minimum((ub2-x0)/sc,k)))
    pred=g+J@(r.x*sc)
    print(f'[h_P free] |step|<={k} sd: predicted loss {2*r.cost:.3f}; h_P -> {x0[4]+r.x[4]*sc[4]:.3f}; gaps rooms {pred[M.index("mean_rooms")]:.3f}, fbrooms {pred[M.index("first_birth_rooms")]:.3f}, early {pred[M.index("early_fertility")]:.3f}')


# ======================================================================
# Part B: linearity checks, Levenberg-Marquardt path, singular values per 1% parameter change
# ======================================================================
import numpy as np, pandas as pd
np.set_printoptions(linewidth=220, suppress=True, precision=4)
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
f='output/model/fixed_reference_economics_20260928/normalized_calibration_v2/deployment/calibration_plateau_diagnosis_v1/all_passed_candidates_10_parameters_10_moments.csv'
d=pd.read_csv(f)
P=[c for c in d.columns if c.startswith('param__')]; names=[p[7:] for p in P]
M=[c[:-7] for c in d.columns if c.endswith('__model')]
W=np.array([d[m+'__weight'].iloc[0] for m in M]); T=np.array([d[m+'__target'].iloc[0] for m in M])
best=d.loc[d.base_loss.idxmin()]; x0=best[P].values.astype(float); X=d[P].values; sc=X.std(0)
Y=d[[m+'__model' for m in M]].values
# check: does a linear-in-parameters model of moments reproduce the loss?
A=np.column_stack([np.ones(len(d)),(X-x0)/sc]); C=np.linalg.lstsq(A,Y,rcond=None)[0]
Yhat=A@C; loss_hat=((Yhat-T)**2*W).sum(1)
print('corr(actual loss, linear-moment loss)=%.4f; max abs err=%.3f; err at best=%.3f'%(np.corrcoef(d.base_loss,loss_hat)[0,1],np.abs(d.base_loss-loss_hat).max(),loss_hat[d.base_loss.idxmin()]-best.base_loss))
# quadratic check for first_birth_rooms and mean_rooms: add squares
A2=np.column_stack([A,((X-x0)/sc)**2]); C2=np.linalg.lstsq(A2,Y,rcond=None)[0]; Y2=A2@C2
for m in ('first_birth_rooms','mean_rooms','recent_parent_ownership','early_fertility'):
    i=M.index(m); print(m,'R2 lin %.4f quad-diag %.4f'%(1-((Y[:,i]-Yhat[:,i])**2).sum()/((Y[:,i]-Y[:,i].mean())**2).sum(),1-((Y[:,i]-Y2[:,i])**2).sum()/((Y[:,i]-Y[:,i].mean())**2).sum()))
J=C[1:].T   # moments x params, per cloud-sd
g=C[0]-T    # predicted gap at best point (intercept)
print('intercept gaps vs actual gaps at best:'); print(pd.DataFrame({'lin':g,'actual':[best[m+'__gap'] for m in M]},index=M).round(5))
Jw=np.sqrt(W)[:,None]*J; gw=np.sqrt(W)*g
grad=2*Jw.T@gw
print('loss gradient per cloud-sd:',dict(zip(names,grad.round(3))),' |grad|=%.3f'%np.linalg.norm(grad))
print('Levenberg-Marquardt path (Euclidean step length in cloud-sd, predicted loss, gaps for rooms/fbrooms/early/rpo):')
for lam in (100,10,3,1,0.3,0.1,0.03,0.01,0.003,0.001):
    dl=-np.linalg.solve(Jw.T@Jw+lam*np.eye(10),Jw.T@gw); pred=g+J@dl
    print('lam=%7.3f |step|=%7.2f sd  loss=%7.3f  rooms %+.3f fb %+.3f early %+.3f rpo %+.4f own %+.4f  h_P=%.3f psi=%.4f tk=%.4f cbc=%.4f'%(lam,np.linalg.norm(dl),(W*pred**2).sum(),pred[M.index('mean_rooms')],pred[M.index('first_birth_rooms')],pred[M.index('early_fertility')],pred[M.index('recent_parent_ownership')],pred[M.index('ownership_30_55')],x0[4]+dl[4]*sc[4],x0[7]+dl[7]*sc[7],x0[8]+dl[8]*sc[8],x0[2]+dl[2]*sc[2]))
# how far are sampled points from best (in cloud-sd Euclid)? 
r=np.linalg.norm((X-x0)/sc,axis=1); print('sample distance from best (cloud-sd): quantiles',np.quantile(r,[.05,.25,.5,.75,.95]).round(2))
# best attainable for each row alone (linear, |step|<=5sd Euclid): direction of steepest change
print('max change in each moment from a 5-sd Euclidean step vs gap:')
for i,m in enumerate(M):
    print('  %-26s reach=%.4f  gap=%+.4f  ratio gap/reach=%.1f'%(m,5*np.linalg.norm(J[i]),g[i],abs(g[i])/(5*np.linalg.norm(J[i]))))
# percent-scaled Jacobian: effect of +1% in each parameter on weighted gaps
Jp=pd.DataFrame(np.sqrt(W)[:,None]*(J/sc)*(0.01*x0)[None,:],index=M,columns=names)
print('weighted-gap change from +1% of each parameter:'); print(Jp.round(3))
U,s,Vt=np.linalg.svd(Jp.values); print('singular values (per 1% param change):',s.round(4))
for k in (-1,-2,-3):
    print(' weak dir',k,dict(zip(names,Vt[k].round(2))))
print(' strongest moment loadings U[:,0..2]'); print(pd.DataFrame(U[:,:3],index=M).round(2))
print(' moments least spanned (U last cols)'); print(pd.DataFrame(U[:,-3:],index=M).round(2))


# ======================================================================
# Part C: standard errors and subsample stability
# ======================================================================
import numpy as np, pandas as pd
np.set_printoptions(linewidth=220, suppress=True, precision=4)
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
f='output/model/fixed_reference_economics_20260928/normalized_calibration_v2/deployment/calibration_plateau_diagnosis_v1/all_passed_candidates_10_parameters_10_moments.csv'
d=pd.read_csv(f)
P=[c for c in d.columns if c.startswith('param__')]; names=[p[7:] for p in P]
M=[c[:-7] for c in d.columns if c.endswith('__model')]
W=np.array([d[m+'__weight'].iloc[0] for m in M]); T=np.array([d[m+'__target'].iloc[0] for m in M])
def jac(dd):
    X=dd[P].values; x0=X.mean(0); A=np.column_stack([np.ones(len(dd)),X-x0]); Y=dd[[m+'__model' for m in M]].values
    C=np.linalg.lstsq(A,Y,rcond=None)[0]; res=Y-A@C; s2=(res**2).sum(0)/(len(dd)-11); XtXi=np.linalg.inv(A.T@A)
    se=np.sqrt(np.outer(np.diag(XtXi)[1:],s2)).T
    return C[1:].T,se
J,se=jac(d)
rows=['mean_rooms','first_birth_rooms','recent_parent_ownership','ownership_30_55','early_fertility','cps_childlessness']
cols=['h_P','psi_child','tenure_choice_kappa','kappa_fert_continuation','chi','first_birth_fixed_cost','child_benefit_curvature']
print('Jacobian entries with OLS s.e. (all 1277 points)')
for r in rows:
    i=M.index(r); print(' ',r,'; '.join('%s %.3f (%.3f)'%(c,J[i,names.index(c)],se[i,names.index(c)]) for c in cols))
# robustness: by chain-center groups (split by chain mod) and by loss tercile
for lab,sub in (('loss<=median',d[d.base_loss<=d.base_loss.median()]),('loss>median',d[d.base_loss>d.base_loss.median()]),('chains 0-11',d[d.chain<12]),('chains 12-23',d[d.chain>=12])):
    Js,_=jac(sub); print(lab,len(sub),'| d rooms/d h_P %.2f, d fb/d h_P %.3f, d rpo/d h_P %.3f, d rooms/d psi %.1f, d fb/d psi %.2f, d rpo/d psi %.2f, d early/d kcont %.3f'%(Js[M.index('mean_rooms'),4],Js[M.index('first_birth_rooms'),4],Js[M.index('recent_parent_ownership'),4],Js[M.index('mean_rooms'),7],Js[M.index('first_birth_rooms'),7],Js[M.index('recent_parent_ownership'),7],Js[M.index('early_fertility'),6]))
# price relation: regress moments on price + params? simple: elasticity of rooms wrt price in cloud
import numpy.linalg as la
print('share of points with h_P==2.3:',(d.param__h_P>=2.3-1e-12).mean(),' h_P quantiles',d.param__h_P.quantile([.01,.1,.5,.9]).values)
print('corr(h_P, loss)=%.3f'%np.corrcoef(d.param__h_P,d.base_loss)[0,1])
# top 20 by loss: parameters
print(d.nsmallest(8,'base_loss')[['chain','base_loss','price']+P].round(5).to_string())
# contributions at best and sums by block
best=d.loc[d.base_loss.idxmin()]
print({m:round(float(best[m+'__loss_contribution']),3) for m in M})


# ======================================================================
# Part D: bounded descent direction, out-of-sample check, weight illustration
# ======================================================================
import numpy as np, pandas as pd
from scipy.optimize import lsq_linear
np.set_printoptions(linewidth=220, suppress=True, precision=4)
f='output/model/fixed_reference_economics_20260928/normalized_calibration_v2/deployment/calibration_plateau_diagnosis_v1/all_passed_candidates_10_parameters_10_moments.csv'
d=pd.read_csv(f)
P=[c for c in d.columns if c.startswith('param__')]; names=[p[7:] for p in P]
M=[c[:-7] for c in d.columns if c.endswith('__model')]
W=np.array([d[m+'__weight'].iloc[0] for m in M]); T=np.array([d[m+'__target'].iloc[0] for m in M])
lb=np.array([0.94,0.1,0.0,0.0,0.1,0.02,0.02,0.01,0.001,0.0]); ub=np.array([0.99,5,0.8,8,2.3,50,50,0.5,0.1,8])
def fit(dd,x0,sc):
    A=np.column_stack([np.ones(len(dd)),(dd[P].values-x0)/sc]); C=np.linalg.lstsq(A,dd[[m+'__model' for m in M]].values,rcond=None)[0]
    return C[0]-T, C[1:].T
best=d.loc[d.base_loss.idxmin()]; x0=best[P].values.astype(float); sc=d[P].values.std(0)
g,J=fit(d,x0,sc); Jw=np.sqrt(W)[:,None]*J; gw=np.sqrt(W)*g
for k in (2,3,5):
    r=lsq_linear(Jw,-gw,bounds=(np.maximum((lb-x0)/sc,-k),np.minimum((ub-x0)/sc,k)))
    pred=g+J@r.x
    print(f'box +-{k} sd: predicted loss {2*r.cost:.2f} (from {np.sum(gw**2):.2f}); step(sd)',dict(zip(names,r.x.round(1))))
    print('    new params',dict(zip(names,(x0+r.x*sc).round(5))))
    print('    contributions',{m:round(float(W[i]*pred[i]**2),2) for i,m in enumerate(M)})
# out-of-sample check of the linear loss model: fit on chains 0-11, predict chains 12-23 and vice versa
for a,b in ((d.chain<12,d.chain>=12),(d.chain>=12,d.chain<12)):
    tr,te=d[a],d[b]; x0t=tr[P].values.mean(0); sct=tr[P].values.std(0); gt,Jt=fit(tr,x0t,sct)
    pred=gt[None,:]+((te[P].values-x0t)/sct)@Jt.T; lh=(W*pred**2).sum(1)
    print('out-of-sample: corr(pred loss, actual)=%.4f, mean abs err=%.3f, n=%d'%(np.corrcoef(lh,te.base_loss)[0,1],np.abs(lh-te.base_loss).mean(),len(te)))
# how concentrated are the six centers? distance between chain start groups
cent=d.groupby('chain')[P].first(); print('range of chain starting points in sample-sd units per parameter:',((cent.max()-cent.min())/sc).round(1).to_dict())
# inverse-variance reweighting illustration with documented SEs (mean_rooms AHS SE 0.0089; first_birth_rooms A2h SE 0.0503; early fertility bootstrap SE 0.028)
alt=W.copy(); alt[M.index('mean_rooms')]=1/0.0089**2; alt[M.index('first_birth_rooms')]=1/0.0502703781254777**2; alt[M.index('early_fertility')]=1/0.02803450449153517**2
gb=np.array([best[m+'__gap'] for m in M],float)
print('loss at best point under legacy weights: %.2f ; with the three documented sampling SEs substituted: %.1f'%((W*gb**2).sum(),(alt*gb**2).sum()))
print({m:(round(float(W[i]*gb[i]**2),2),round(float(alt[i]*gb[i]**2),1),round(float(abs(gb[i])*np.sqrt(alt[i])),1)) for i,m in enumerate(M) if m in ('mean_rooms','first_birth_rooms','early_fertility')})
