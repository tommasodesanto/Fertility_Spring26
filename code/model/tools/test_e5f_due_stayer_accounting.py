"""Tiny origin-specific accounting fixtures; no Bellman/stationary solves."""
import ast
from pathlib import Path
from types import SimpleNamespace as NS
import unittest
import numpy as np
from e5f_overnight_estate_audit import policy_mass_branches, branch_estate_accounts

ROOT = Path(__file__).parent


def fixture():
    g = np.zeros((1, 2, 1, 1, 1, 1, 1)); g[:, 1] = 1
    stay = .5*g
    bp = np.zeros_like(g); bp[:, 1] = -9
    bs = np.zeros_like(g); bs[:, 1] = -7
    p = NS(bp_pol=bp, bp_pol_stay=bs, c_pol=np.ones_like(g),
           c_pol_stay=np.ones_like(g), hR_pol=np.zeros_like(g), price=np.array([1.]))
    return NS(g_current=g, g_stay_distribution=stay, policy=p), NS(native_due_stayer_credit=True)


def budget_function():
    # Execute only the pure accounting function, never the model-importing driver.
    tree = ast.parse((ROOT/'run_e5f_matched_pf_smoke.py').read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'dated_budget')
    namespace = dict(np=np, model=NS(income_at_state=lambda *a: 0.))
    exec(compile(ast.Module(body=[node], type_ignores=[]), '<budget fixture>', 'exec'), namespace)
    return namespace['dated_budget']


class DueAccountingTests(unittest.TestCase):
    def test_signed_estates_are_integrated_before_netting(self):
        e,p = fixture()
        result = branch_estate_accounts(e,p,[1.],[0.,10.],.2)['totals']
        # Buyers leave net -1, stayers net +1: neither flow may be averaged away.
        self.assertEqual(result['net_positive'],.5)
        self.assertEqual(result['net_negative'],.5)
        self.assertEqual(result['negative_estate_death_mass'],.5)
        self.assertEqual(result['death_mass'],1.)
        self.assertEqual(result['signed_cost_identity_residual'],0.)

    def test_default_off_exact_legacy(self):
        from audit_e5f_estate_resource_account import signed_accounts
        e,p=fixture(); p.native_due_stayer_credit=False
        e.g_stay_distribution=None
        self.assertEqual(branch_estate_accounts(e,p,[1.],[0.,10.],.2),
                         signed_accounts(e.g_current,e.policy.bp_pol,[1.],[0.,10.],.2))

    def test_missing_split_fails(self):
        for field in ('g_stay_distribution','bp_pol_stay','c_pol_stay'):
            e,p=fixture()
            setattr(e if field=='g_stay_distribution' else e.policy,field,None)
            with self.assertRaises(ValueError): policy_mass_branches(e,p)

    def test_invalid_split_fails(self):
        for bad in ('excess','renter','nan','shape'):
            e,p=fixture()
            if bad=='excess':e.g_stay_distribution*=3
            if bad=='renter':e.g_stay_distribution[:,0]=.1
            if bad=='nan':e.g_stay_distribution[:,1]=np.nan
            if bad=='shape':e.g_stay_distribution=e.g_stay_distribution[:,1:]
            with self.assertRaises(ValueError):policy_mass_branches(e,p)

    def test_distinct_consumption_budget(self):
        e,p=fixture()
        p.J=p.n_parity=p.n_child_states=1;p.z_grid=[1.];p.R_gross=1.
        p.H_own=np.array([10.]);p.delta=p.tau_H=0.
        e.policy.c_pol[:,1]=9.;e.policy.c_pol_stay[:,1]=7.
        result=budget_function()(e,p,NS(gb_flat=np.array([0.])),np.array([0.]),1.)
        self.assertEqual(result['budget_excess_mass'],0.)
        e.policy.c_pol_stay[:,1]=8.
        with self.assertRaises(RuntimeError):
            budget_function()(e,p,NS(gb_flat=np.array([0.])),np.array([0.]),1.)

    def test_dated_split_uses_post_fertility_and_owned_policy(self):
        tree=ast.parse((ROOT/'run_dynamic_population_transition.py').read_text())
        node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate_period')
        e,p=fixture(); p._bp_pol_stay=np.full_like(e.g_current,999.)
        p._c_pol_stay=np.full_like(e.g_current,999.)
        policy=e.policy
        policy.loc_probs=policy.tenure_choice=policy.tenure_probs=policy.fert_probs=None
        policy.maps=NS(lmm_idx=None,lmm_wt=None,tmx_idx=None,tmx_wt=None)
        seen=[]
        def stay(g,*args):
            seen.append(g.copy()); return .5*g
        namespace=dict(np=np, joint_nested_enabled=lambda p:False,
            gate_pre_fertility_distribution=lambda g,*args:(g.copy(),0.),
            policy_continuation_birth_probs=lambda *args:None,
            apply_fertility=lambda g,*args:(2*g,0.,np.zeros(1)),
            model=NS(realize_current_cross_section=lambda g,*a,**k:g.copy(),
                     realize_stayer_cross_section=stay),
            housing_demand_by_location=lambda *args:np.ones(1),
            PeriodEvaluation=lambda *args,**kw:NS(**kw))
        exec('from __future__ import annotations\n'+ast.unparse(node),namespace)
        result=namespace['evaluate_period'](np.array([1.]),e.g_current,p,None,None,None,
            supply_rule=NS(quantity=lambda price:np.ones(1)),supplied_policy=policy)
        np.testing.assert_array_equal(seen[0],2*e.g_current)
        np.testing.assert_array_equal(result.g_stay_distribution,e.g_current)
        self.assertIs(p._bp_pol_stay,policy.bp_pol_stay)
        self.assertIs(p._c_pol_stay,policy.c_pol_stay)

    def test_old_pinned_ledger_rejects_due(self):
        tree=ast.parse((ROOT/'e5f_current_transition_runtime.py').read_text())
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='EstateAuditContract')
        namespace={}; exec(compile(ast.Module(body=[node],type_ignores=[]),'<audit guard>','exec'),namespace)
        old=NS(EstateFundingShortfall=RuntimeError,audit=lambda *a,**k:'legacy')
        guarded=namespace['EstateAuditContract'](old)
        e,p=fixture()
        with self.assertRaises(RuntimeError):guarded.audit(e,p,None)
        p.native_due_stayer_credit=False
        self.assertEqual(guarded.audit(e,p,None),'legacy')

    def test_no_input_mutation(self):
        e,p=fixture();old=e.g_current.copy();stay=e.g_stay_distribution.copy()
        policy_mass_branches(e,p)
        np.testing.assert_array_equal(e.g_current,old)
        np.testing.assert_array_equal(e.g_stay_distribution,stay)


if __name__=='__main__':unittest.main()
