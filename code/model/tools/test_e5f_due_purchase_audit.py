"""Small direct purchase/stayer accounting fixtures; no model solves."""
import unittest
from types import SimpleNamespace as NS
import numpy as np
import e5f_due_purchase_audit as audit


def fixture(*, wealth=-100., saving=-100., age=0, psi=0., R=1.1):
    bg=np.array([-100.,-60.,0.,100.]);shape=(4,2,1,2,1,1,1)
    g=np.zeros(shape);stay=np.zeros(shape);post=np.zeros(shape)
    b=int(np.flatnonzero(bg==wealth)[0]);g[b,1,0,age,0,0,0]=.5;stay[b,1,0,age,0,0,0]=.5
    post[b,1,0,age,0,0,0]=.5
    # Buyer enters same final wealth -100 from initial renter wealth zero.
    g[0,1,0,age,0,0,0]+=.5;post[2,0,0,age,0,0,0]=.5
    bp=np.zeros(shape);bp[:,1]=-80.;bs=bp.copy();bs[b,1,0,age,0,0,0]=saving
    probs=np.zeros((4,2,1,2,1,1,1,2));probs[...,1]=1.
    maps=np.zeros((1,2,2,1,1,4),int);weights=np.zeros_like(maps,float)
    for old in range(2):
      for new in range(2):
        x=bg if old==new else bg+(1-psi)*100*old-100*new
        idx=np.clip(np.searchsorted(bg,x,side='right')-1,0,2)
        maps[0,old,new,0,0]=idx;weights[0,old,new,0,0]=np.clip((x-bg[idx])/(bg[idx+1]-bg[idx]),0,1)
    policy=NS(price=np.array([1.]),bp_pol=bp,bp_pol_stay=bs,c_pol=np.ones(shape),c_pol_stay=np.ones(shape),
              tenure_probs=probs,loc_probs=np.ones((4,2,1,1,2,1)),maps=NS(tmx_idx=maps,tmx_wt=weights))
    e=NS(g_current=g,g_stay_distribution=stay,g_post_fertility=post,policy=policy)
    p=NS(native_due_stayer_credit=True,I=1,J=2,z_grid=[1.],H_own=[100.],psi=psi,R_gross=R,
         n_parity=1,n_child_states=1,use_age_survival=True,survival_probs=np.array([1.]))
    shared=NS(phi_choice=np.full((1,2,1,1),.8))
    return e,p,shared,bg,NS(income_at_state=lambda *args:200.)


class DuePurchaseTests(unittest.TestCase):
    def test_mixed_origins_same_node(self):
        receipt=audit.audit_purchase_accounting(*fixture())
        self.assertEqual(receipt['stayer_mass'],.5)
        self.assertEqual(receipt['end_mortgage_floor_violation_mass'],0.)
        self.assertEqual(receipt['stayer_principal_floor_violation_mass'],0.)

    def test_buyer_not_grandfathered(self):
        args=fixture();args[0].policy.bp_pol[0,1]=-100.
        with self.assertRaises(RuntimeError):audit.audit_purchase_accounting(*args)

    def test_within_ltv_cashout(self):
        audit.audit_purchase_accounting(*fixture(wealth=-60.,saving=-80.))
        with self.assertRaises(RuntimeError):audit.audit_purchase_accounting(*fixture(wealth=-60.,saving=-81.))

    def test_above_ltv_no_additional_principal(self):
        with self.assertRaises(RuntimeError):audit.audit_purchase_accounting(*fixture(saving=-101.))

    def test_interest_not_capitalized(self):
        audit.audit_purchase_accounting(*fixture(saving=-100.,R=1.2))
        with self.assertRaises(RuntimeError):audit.audit_purchase_accounting(*fixture(saving=-120.,R=1.2))

    def test_death_estate_floor_only_at_risk(self):
        audit.audit_purchase_accounting(*fixture(saving=-100.,psi=.1))
        with self.assertRaises(RuntimeError):audit.audit_purchase_accounting(*fixture(saving=-100.,age=1,psi=.1))
        audit.audit_purchase_accounting(*fixture(saving=-90.,age=1,psi=.1))
        args=fixture(saving=-100.,psi=.1);args[1].survival_probs[0]=.99
        with self.assertRaises(RuntimeError):audit.audit_purchase_accounting(*args)

    def test_missing_and_invalid_split(self):
        args=fixture();args[0].g_stay_distribution=None
        with self.assertRaises(ValueError):audit.audit_purchase_accounting(*args)
        args=fixture();args[0].g_stay_distribution*=3
        with self.assertRaises(ValueError):audit.audit_purchase_accounting(*args)

    def test_default_off_exact_receipt(self):
        args=fixture(saving=-80.);args[1].native_due_stayer_credit=False
        self.assertEqual(audit.audit_purchase_accounting(*args),audit._inherited_audit()(*args))

    def test_transaction_grid_check_retained(self):
        args=fixture();args[0].policy.maps.tmx_wt[0,0,1,0,0,2]=.2
        with self.assertRaises(RuntimeError):audit.audit_purchase_accounting(*args)


if __name__=='__main__':unittest.main()
