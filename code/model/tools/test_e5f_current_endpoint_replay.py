import unittest
from types import SimpleNamespace as NS
import numpy as np
import run_e5f_current_endpoint_replay as r

class EndpointReplay(unittest.TestCase):
    def test_pinned_reference_is_repeated_and_closed(self):
        value=r.authenticate_endpoint()
        self.assertTrue(value['repeat_verified'])
        self.assertLess(abs(value['endpoint']['renewal_residual']),2.5e-5)

    def test_parameter_changes_not_hidden_by_flags(self):
        old=NS(psi_child=.1,H0=np.array([6.]),eq_iter=2)
        new=NS(psi_child=.1,H0=np.array([6.]),eq_iter=10,native_solvency_credit=True)
        r.unchanged_parameters(old,new)
        new.psi_child=.2
        with self.assertRaises(RuntimeError):r.unchanged_parameters(old,new)

    def test_closed_accounting_uses_supply_demand_scale(self):
        P=NS(H0=np.array([r.POPULATION]),user_cost_rate=1.,r_bar=np.array([r.PRICE]),xi_supply=np.array([.63]))
        sol=NS(entry_rate=.1,adult_entry_adjusted_birth_children=.21,housing_demand=np.array([1.]))
        a=r.closed_accounting(sol,P)
        self.assertAlmostEqual(a['renewal_residual'],0.)
        self.assertEqual(a['population_scale'],r.POPULATION)
        sol.adult_entry_adjusted_birth_children=.3
        with self.assertRaises(RuntimeError):r.closed_accounting(sol,P)

if __name__=='__main__':unittest.main()
