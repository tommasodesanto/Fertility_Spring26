from pathlib import Path
import importlib.util,sys,unittest
import numpy as np
HERE=Path(__file__).resolve().parents[1];ROOT=HERE.parents[3]
sys.path.insert(0,str(HERE))
from refactor_lab.engine import shared,kernels

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
class ScaleBinding(unittest.TestCase):
    def test_source_bundle_current_children_direct_benefit(self):
        inputs=module(ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_round2_v1/inputs.py','input_scale_test')
        P,bg=inputs.proposal();P.ces_enabled=True;P.ces_eta=.487;P.lambda_housing=.2;P.alpha_cons=.733
        P.delta_alpha=P.delta_alpha_jump=0.;P.compensated_child_housing_shares=False;P.child_room_floor=False;P.hbar_first_child_jump=P.hbar_child_rooms=0.
        P.child_benefit_curvature=.15;P.psi_child=.139
        sd=shared.precompute_shared(P,bg)
        self.assertTrue(np.all(sd.h_bar==0));self.assertTrue(np.all(sd.c_bar==0));self.assertTrue(np.all(sd.escale_flat==1));self.assertTrue(np.all(sd.alpha_flat==.733))
        ec=sd.ces_ec_flat.reshape((P.n_parity,P.n_child_states),order='F');eh=sd.ces_eh_flat.reshape(ec.shape,order='F')
        for n in range(P.n_parity):
            for m in range(P.n_child_states):
                valid=m if m<=n else 0
                self.assertEqual(ec[n,m],((2+.7*valid)/2)**.7)
                self.assertEqual(eh[n,m],ec[n,m]*(1+.2*(valid>0)))
                if m<=n:self.assertEqual(sd.psi_v[n,m],.139*m**.85 if m else 0.)
    def test_original_cd_exhaustive_exact_regression(self):
        old=module(ROOT/'code/model/refactor_lab/engine/kernels.py','old_cd_regression')
        bg=np.array([0.,.1,.2,.5,1.,2.]);V=np.array([-2.,-1.8,-1.5,-1.1,-.8,-.5])
        for owner in (False,True):
            args=(0.,1.5,2.,V,bg,.3,0.,0.,.15,8.,.733,-1.,.97,1.,.1,4.**(-.267),owner)
            self.assertEqual(kernels.exhaustive_saving_scalar(*args),old.exhaustive_saving_scalar(*args))
if __name__=='__main__':unittest.main()
