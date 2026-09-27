"""Zero-solve tests for permanent-price-fall bounds and audit gates."""
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import run_e5f_due_price_fall_check as r

class PriceFallPlan(unittest.TestCase):
    def plan(self):
        return dict(arms=['baseline','due'],price_multiplier=.9,expectations='permanent_constant_price',
                    fiscal_closure='fixed_tax_fixed_pension_partial_equilibrium',seconds_per_arm=300,total_seconds=900,
                    initial_state='exact_original_stationary_g_pre',target_reporting='none_dated_distribution')

    def test_exact_scope_and_bounds(self):
        p=self.plan();r.validate_spec(p)
        for key,value in [('price_multiplier',.8),('expectations','one_date_surprise'),('seconds_per_arm',301),
                          ('total_seconds',901),('initial_state','reweighted'),('target_reporting','stationary_fit'),
                          ('fiscal_closure','balanced')]:
            q=dict(p);q[key]=value
            with self.assertRaises(ValueError):r.validate_spec(q)

    def test_requires_live_parent_and_global_end(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'parent.json'
            receipt=dict(status='supervising',maximum_active_children=1,seconds_per_arm=300,
                         start_epoch=100.,absolute_end_epoch=1000.,pid=333)
            path.write_text(json.dumps(receipt))
            p=self.plan();p.update(execution_authorized=True,absolute_end_epoch=1000.,parent_supervisor_receipt=str(path),
                                  parent_supervisor_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            with patch.object(r.os,'kill') as alive:
                r.verify_supervision(p,now=500.);alive.assert_called_once_with(333,0)
            for now in (99.,1000.,1001.):
                with self.assertRaises(ValueError):r.verify_supervision(p,now=now)
            p['execution_authorized']=False
            with self.assertRaises(ValueError):r.verify_supervision(p,now=500.)

    def test_signed_estates_and_stayer_violation_cannot_pass(self):
        budget={'budget_excess_mass':0.}
        purchase={'maximum_occupied_transaction_wealth_error':0.,'stayer_death_solvency_violation_mass':0.,'end_mortgage_floor_violation_mass':0.}
        estate={'status':'funded','audit_id':'estate_funded_dated_entry_provisional_net_v1','estate':{'totals':{'net_negative':0.}}}
        self.assertTrue(all(r.audit_gates(budget,purchase,estate,0.).values()))
        bad=copy.deepcopy(estate);bad['estate']['totals']['net_negative']=1e-8
        self.assertFalse(r.audit_gates(budget,purchase,bad,0.)['negative_estates'])
        badp=dict(purchase,stayer_death_solvency_violation_mass=1e-8)
        self.assertFalse(r.audit_gates(budget,badp,estate,0.)['stayer_death_solvency_violation_mass'])
        bad=copy.deepcopy(estate);bad['audit_id']='stationary_entry'
        self.assertFalse(r.audit_gates(budget,purchase,bad,0.)['dated_entry'])

    def test_source_has_no_stationary_solve_or_target_scoring(self):
        source=Path(r.__file__).read_text()
        self.assertNotIn('solve_markov_income_at_prices(',source)
        self.assertNotIn('score_targets(',source)
        self.assertIn('next_entrant_cohort=next_cohort',source)
        self.assertIn('cal._require_exact_inherited_distribution(g0,policy,P,grid)',source)

if __name__=='__main__':unittest.main()
