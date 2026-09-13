"""Read saved state only; reduce dependent-child exposure to parent-household exit."""
import gzip, pickle, sys, json, hashlib
from types import SimpleNamespace
import numpy as np
class SnapshotUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module.startswith(('run_', 'e5f_', 'intergen_', 'demographic_transition')):
            return type(name, (SimpleNamespace,), {})
        return super().find_class(module, name)
path=sys.argv[1]
with gzip.open(path,'rb') as f:
    q=SnapshotUnpickler(f).load()
P=q['parameters']; e=q['evaluation']; g=np.asarray(e.g_post_fertility)
assert P.child_state_mode=='independent_count'
assert g.ndim==7 and g.shape[3]==P.J
am=g.sum(axis=(0,1,2,4,5)); m=np.arange(am.shape[1]); C=am@m
s=np.ones(P.J)
if P.use_age_survival:
    s[:-1]=np.asarray(P.survival_probs,dtype=float)[:P.J-1]
s[-1]=0.0
D=(1-s)*C
parents=am[:,1:].sum(axis=1)
# Validate counts are unchanged by current economic choices.
C_current=np.asarray(e.g_current).sum(axis=(0,1,2,4,5))@m
assert np.allclose(C,C_current,atol=1e-8,rtol=1e-8)
# Literal dependents, no 3+ top-bin expansion: same units as household m.
by_anm=g.sum(axis=(0,1,2,4))
mu_losses=np.zeros((P.n_parity,P.n_child_states))
for n in range(P.n_parity):
    for cm in range(P.n_child_states):
        mu_losses[n,cm]=cm-np.dot(P.Pi_child[cm,:,n],m)
M=np.einsum('anm,nm,a->a',by_anm,mu_losses,s)
remaining=s*C-M
assert np.allclose(C,D+M+remaining,atol=1e-10)
rows=[dict(age=float(P.age_start+a*P.da),dependents=float(C[a]),exit_probability=float(1-s[a]),dependent_loss=float(D[a]),parent_households=float(parents[a]),affected_households=float((1-s[a])*parents[a]),maturation=float(M[a])) for a in range(P.J)]
result=dict(first_fertile_age=float(P.age_start+(P.A_f_start-1)*P.da),last_fertile_age=float(P.age_start+(P.A_f_end-1)*P.da),child_maturation_probability=float(P.Pi_child[1,0,1]),checkpoint=path,checkpoint_sha256=hashlib.sha256(open(path,'rb').read()).hexdigest(),period_years=float(P.period_years),children_units='Literal dependent count; 3+ state counts as 3',household_mass=float(g.sum()),dependents=float(C.sum()),births=float(e.births),dependents_losing_parent=float(D.sum()),loss_percent_dependents=float(100*D.sum()/C.sum()),loss_percent_births=float(100*D.sum()/e.births),terminal_age_loss=float(D[-1]),terminal_share_of_loss=float(D[-1]/D.sum()),nonterminal_loss=float(D[:-1].sum()),maturation=float(M.sum()),loss_percent_maturation=float(100*D.sum()/M.sum()),affected_households=float(np.sum((1-s)*parents)),age_rows=rows,max_current_child_count_gap=float(np.max(np.abs(C-C_current))),stock_flow_identity_max_error=float(np.max(np.abs(C-D-M-remaining))))
print(json.dumps(result,indent=2))
