"""Contract checks: additive shocks cannot replace the maintained coordinates."""
import sys,unittest,copy
from pathlib import Path
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
import run_e5f_transition_calibration as cal

class ContractTests(unittest.TestCase):
    def args(self):return SimpleNamespace(two_shock_choice=True,fixed_tenure_choice_kappa=.005,
        housing_supply_elasticity=.63,model_profile=cal.REPAIRED_MODEL_PROFILE,
        estimate_first_child_room_jump=True,target_profile=cal.E5_BASELINE_TARGET_PROFILE)
    def test_original_parameters_retained(self):
        theta={'kappa_fert':2.1,'kappa_fert_continuation':1.7}
        domain,overrides,profile=cal.activate_model_profile(cal.REPAIRED_MODEL_PROFILE,theta)
        domain,overrides,profile=cal.configure_first_child_room_jump(model_profile_name=cal.REPAIRED_MODEL_PROFILE,
            theta=theta,active_domain=domain,profile_overrides=overrides,model_profile=profile,fixed_jump=None,estimate_jump=True)
        before=copy.deepcopy((domain,theta));cal.configure_two_shock_calibration(self.args(),overrides,profile)
        self.assertEqual((domain,theta),before)
        self.assertEqual(len(domain),11)
        names={r[0] for r in domain}
        self.assertTrue({'kappa_fert','kappa_fert_continuation'}<=names)
        self.assertTrue({'tenure_choice_kappa','joint_nest_lambda'}.isdisjoint(names))
        self.assertTrue(overrides['two_shock_choice'] and overrides['joint_nested_choice'])
    def test_rejects_changed_external_restrictions_and_target_system(self):
        for key,value in [('fixed_tenure_choice_kappa',.01),('housing_supply_elasticity',1.),
                          ('estimate_first_child_room_jump',False),('target_profile','young_ownership')]:
            args=self.args();setattr(args,key,value)
            with self.assertRaises(ValueError):cal.configure_two_shock_calibration(args,{}, {})

if __name__=='__main__':unittest.main()
