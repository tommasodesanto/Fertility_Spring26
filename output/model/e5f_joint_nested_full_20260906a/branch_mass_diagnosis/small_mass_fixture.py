"""Historical identity-transport diagnosis on explicitly supplied original source.

Use --model-root /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907g/code/model
on Torch. This intentionally demonstrates the original Python scatter loss too.
The complete repaired fixture is test_e5f_small_mass_transport.py.
"""
from pathlib import Path
from types import SimpleNamespace as NS
import argparse,ast,hashlib,inspect,json,sys,textwrap
import numpy as np
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--model-root',type=Path,required=True);args=parser.parse_args()
original=args.model_root/'intergen_eqscale_seq_optimized/solver.py'
assert hashlib.sha256(original.read_bytes()).hexdigest()=='51c921221a37f2c0da1f19b2ec0c1e016d61699e6a47b82d458eeff2e012dee8'
sys.path.insert(0,str(args.model_root))
from intergen_eqscale_seq_optimized import solver
b=np.array([0.,1.]);Nb=2;nt=2;I=1;J=2;Nz=1;npar=2;ncs=2
P=NS(n_house=1,I=I,n_parity=npar,n_child_states=ncs,n_child_stages=0,use_numba_scatter=False)
SD=NS(nc=npar*ncs)
shape=(Nb,nt,I,J,Nz,npar,ncs)
loc=np.ones((Nb,nt,I,I,J,Nz,npar,ncs));ten=np.zeros(shape+(nt,))
ten[...,0]=1-1e-7;ten[...,1]=1e-7
bp=np.broadcast_to(b.reshape(Nb,1,1,1,1,1,1),shape).copy()
lidx=np.zeros((I,nt,Nb),dtype=np.int64);lwt=np.broadcast_to(b,(I,nt,Nb)).copy()
tidx=np.zeros((I,nt,nt,npar,ncs,Nb),dtype=np.int64);twt=np.broadcast_to(b,tidx.shape).copy()
g=np.zeros((Nb,nt,I,Nz,npar,ncs));g[0,0,0,0,1,1]=1
args=(0,loc,np.zeros(shape,dtype=np.int64),ten,bp,P,b,SD,lidx,lwt,tidx,twt,False,None,np.ones((1,1)))
tree=ast.parse(textwrap.dedent(inspect.getsource(solver.advance_cohort_one_period_markov_income)))
changes=0
for node in ast.walk(tree):
 if isinstance(node,ast.If) and isinstance(node.test,ast.Compare):
  t=node.test
  if len(t.ops)==1 and isinstance(t.ops[0],ast.Lt) and len(t.comparators)==1 and isinstance(t.comparators[0],ast.Constant) and t.comparators[0].value==1e-15:
   t.ops=[ast.LtE()];t.comparators=[ast.Constant(0.)];changes+=1
assert changes==3
scope=dict(solver.__dict__);exec(compile(ast.fix_missing_locations(tree),'<zero-only-transport>','exec'),scope)
strict=scope['advance_cohort_one_period_markov_income']
rows=[]
for mass in [1.,3.55e-9,1e-13]:
 orig=solver.advance_cohort_one_period_markov_income(mass*g,*args)
 zero=strict(mass*g,*args)
 expected=mass*g.copy();expected[0,0,0,0,1,1]=mass*(1-1e-7);expected[0,1,0,0,1,1]=mass*1e-7
 rows.append(dict(mass=mass,original_relative_mass_gap=abs(float(orig.sum())-mass)/mass,zero_only_relative_mass_gap=abs(float(zero.sum())-mass)/mass,zero_only_relative_l1_error=float(np.abs(zero-expected).sum())/mass,original_owner_mass=float(orig[:,1].sum()),expected_owner_mass=mass*1e-7,zero_only_owner_mass=float(zero[:,1].sum())))
print(json.dumps(dict(status='closed_form_fixture',transition_comparisons_changed=changes,compiled_scatter=False,rows=rows),indent=2))
