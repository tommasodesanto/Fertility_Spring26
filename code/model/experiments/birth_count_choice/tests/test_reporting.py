"""Zero-solve tests of count-menu observer weights and input preservation."""
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
from model.engine.birth_count import birth_count_transition, identity_realized_probabilities
from model.observer_adapters import apply_count_fertility, first_birth_accounting_by_age, count_flows_and_risks


class ReportingTests(unittest.TestCase):
    def test_cap_two_observers_use_configured_support(self):
        pre = np.zeros((1, 1, 1, 1, 1, 4, 4))
        pre[0, 0, 0, 0, 0, 0, 0] = 1.
        rp = identity_realized_probabilities(pre.shape)
        rp[0, 0, 0, 0, 0, 0, 0] = [0., .25, .75, 0.]
        P = SimpleNamespace(J=1, birth_count_choice_cap=2)
        post, births, _ = apply_count_fertility(pre, None, P, birth_count_realized_probs=rp)
        self.assertAlmostEqual(births, 1.75)
        evaluation = SimpleNamespace(g_pre=pre, g_post_fertility=post,
            policy=SimpleNamespace(birth_count_realized_probs=rp))
        flows, risks = count_flows_and_risks(evaluation, P)
        np.testing.assert_array_equal(flows[0], [1., .75, 0.])
        np.testing.assert_array_equal(risks[0], [1., 1., 0.])
        self.assertAlmostEqual(first_birth_accounting_by_age(evaluation, P)['flow'][0], 1.)
        rp[0, 0, 0, 0, 0, 0, 0] = [0., .25, .5, .25]
        with self.assertRaises(ValueError):
            apply_count_fertility(pre, None, P, birth_count_realized_probs=rp)

    def test_crossed_orders_count_children_and_first_events_once(self):
        pre = np.zeros((2, 2, 1, 17, 1, 4, 4))
        pre[0, 0, 0, 3, 0, 0, 0] = 1
        rp = identity_realized_probabilities(pre.shape)
        rp[0, 0, 0, 3, 0, 0, 0] = [.1, .2, .3, .4]
        post = birth_count_transition(pre, rp)
        np.testing.assert_allclose(post['births_by_order'], [.9, .7, .4])
        self.assertAlmostEqual(post['first_birth_tagged_post'].sum(), .9)
        self.assertAlmostEqual(post['expected_births'], 2.)
        self.assertAlmostEqual(post['post'].sum(), 1.)
        evaluation = SimpleNamespace(g_pre=pre, g_post_fertility=post['post'],
                                     policy=SimpleNamespace(birth_count_realized_probs=rp))
        P = SimpleNamespace(J=17)
        first = first_birth_accounting_by_age(evaluation, P)
        self.assertAlmostEqual(first['flow'][3], .9)
        self.assertAlmostEqual(first['hazard'][3], .9)
        flows, risks = count_flows_and_risks(evaluation, P)
        np.testing.assert_allclose(flows[3], [.9, .7, .4])
        np.testing.assert_allclose(risks[3], [1, 1, 1])

    def test_recent_parent_weights_use_households_and_origin_status(self):
        from model.inputs import load_inputs
        from model.reporting import build_context
        P, grid = load_inputs(); P.birth_count_choice_enabled = True
        P.birth_count_choice_cap = 3
        with tempfile.TemporaryDirectory() as directory:
            functions=[]; codes=[]
            for index in range(3):
                out=Path(directory)/str(index); out.mkdir()
                context = build_context(P, grid, out, price_start=.776,
                                        deadline=time.time()+120, max_lifecycle=32, closure='fixed_h0')
                wrapped=context['prepared'].rt['observe_recent_parent_flow']
                module=next(cell.cell_contents for cell in wrapped.__closure__
                            if hasattr(cell.cell_contents,'observe_recent_parent_flow'))
                functions.append(module.observe_recent_parent_flow)
                codes.append(module.observe_recent_parent_flow.__code__)
                self.assertIs(module.observe_recent_parent_flow.__globals__['_production_model_facade'],
                              context['prepared'].rt['model'])
                self.assertIn('selected_households',module.observe_recent_parent_flow.__code__.co_varnames)
            self.assertTrue(all(fn is functions[0] for fn in functions))
            self.assertTrue(all(code is codes[0] for code in codes))
            rt = context['prepared'].rt; facade = rt['model']
            nb, nt, loc, ages, nz = 2, 6, 1, 17, 1
            shape = (nb, nt, loc, ages, nz, 4, 4)
            pre = np.zeros(shape)
            # One childless renter reaches n=1/2/3 with equal probability;
            # one former parent owner receives one child; two controls wait.
            pre[0, 0, 0, 3, 0, 0, 0] = 1
            pre[0, 1, 0, 3, 0, 1, 0] = 1
            pre[1, 0, 0, 3, 0, 0, 0] = 2
            rp = identity_realized_probabilities(shape)
            rp[0, 0, 0, 3, 0, 0, 0] = [0, 1/3, 1/3, 1/3]
            rp[0, 1, 0, 3, 0, 1, 0] = [0, 1, 0, 0]
            result = birth_count_transition(pre, rp)
            policy = SimpleNamespace(fert_probs=np.zeros(shape[:-2]+(4,)),
                fert2_probs=np.zeros(shape[:-2]+(2,2,4)), birth_count_realized_probs=rp,
                V=np.zeros(shape), joint_choice=None, loc_probs=np.ones((nb,nt,loc,loc,ages,nz,4,4)),
                tenure_choice=np.zeros(shape,dtype=int), tenure_probs=None,
                price=np.array([.776]), maps=SimpleNamespace(
                    lmm_idx=np.zeros((loc,nt,nb),dtype=int),lmm_wt=np.zeros((loc,nt,nb)),
                    tmx_idx=np.zeros((loc,nt,nt,4,4,nb),dtype=int),tmx_wt=np.zeros((loc,nt,nt,4,4,nb))))
            evaluation = SimpleNamespace(policy=policy,g_pre=pre,g_post_fertility=result['post'],
                                         g_current=result['post'],births=result['expected_births'],
                                         feasibility_projection_mass=0.)
            original = facade.realize_current_cross_section
            facade.realize_current_cross_section = lambda mass,*args,**kwargs: mass.copy()
            try:
                recent = rt['observe_recent_parent_flow'](evaluation,P,diagnostic_enabled=True,
                    snapshot=rt['SNAPSHOT'],age_projection=rt['AGE_PROJECTION'],
                    diagnostic_allow_residence_proxy=True)
            finally:
                facade.realize_current_cross_section = original
            self.assertAlmostEqual(recent['groups']['selected_birth']['denominator'], 2.)
            self.assertAlmostEqual(recent['groups']['first_birth']['denominator'], 1.)
            self.assertAlmostEqual(recent['groups']['continuation_birth']['denominator'], 1.)
            self.assertAlmostEqual(recent['model_value'], .5)
            self.assertAlmostEqual(recent['accounting']['all_births'], 3.)

    def test_driver_preserves_every_canonical_input(self):
        import json
        import subprocess
        receipt = json.loads(subprocess.check_output([sys.executable, str(BASE / 'run.py'), '--preflight'], text=True))
        self.assertEqual(receipt['lifecycle_solves'], 0)
        self.assertEqual(receipt['canonical_input_fields'], 245)
        self.assertEqual(receipt['changed_input_fields'], ['birth_count_choice_enabled'])
        self.assertEqual((receipt['target_rows'], receipt['parameter_rows']), (14,31))


if __name__ == '__main__':
    unittest.main()
