"""Contract-level scoring and estate tests; run inside Torch only."""
import copy
import os
import sys
import unittest

if sys.platform!='linux' or not os.environ.get('SLURM_JOB_ID'):
    raise RuntimeError('Run calibration runtime tests only in a Torch allocation')

from e5f_calibration_runtime import score_targets, validate_normalization, NonpositiveNormalizedBenefit
from e5f_overnight_estate_audit import run_self_tests


def inputs():
    fertility={'uniform_birth_time':{'moments':dict(
        childless_rate_40_44=.2,exactly_one_among_mothers_40_44=.21,
        period_mean_age_first_birth=26.,period_share_first_births_age30plus=.25,
        mean_children_ever_born_capped3_age25=.8)}}
    housing={'moments':dict(aggregate_wealth_to_annual_gross_labor_earnings=6.9,
        annual_bequest_flow_to_aggregate_wealth=.007,old_total_wealth_to_annual_income_p90_p50_7684=3.5,
        aggregate_mean_occupied_rooms_ahs_uncapped_18_85=5.8,
        aggregate_mean_occupied_rooms_capped9_18_85=5.6,
        own_rate_30_55=.67,housing_increment_0to1=1.4,
        prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9=.38)}
    names=['initial_normalization','cps_childlessness','cps_exactly_one','nchs_mean_age','nchs_share30',
           'early_fertility','wealth_earnings','bequest_wealth','old_dispersion','mean_rooms',
           'ownership_30_55','first_birth_rooms','family_rooms','recent_parent_ownership']
    target=dict(cps_projection='uniform_birth_time',target_rows=[dict(restriction_id=name,
        target=2.1 if name=='initial_normalization' else 1.,
        actual_weight=None if name=='initial_normalization' else 100.) for name in names])
    return target,fertility,housing,.12,2.1


class RegistryTests(unittest.TestCase):
    def test_all_fourteen_rows_and_thirteen_scores_preserve_definitions(self):
        args=inputs(); before=copy.deepcopy(args); rows=score_targets(*args)
        self.assertEqual(args,before)
        self.assertEqual(len(rows),14)
        self.assertEqual(sum(row['weight']!='' for row in rows),13)
        by_name={row['moment']:row for row in rows}
        self.assertEqual(by_name['mean_rooms']['model'],5.8)
        self.assertEqual(by_name['family_rooms']['model'],.38)
        self.assertEqual(by_name['early_fertility']['model'],.8)
        self.assertAlmostEqual(by_name['early_fertility']['loss_contribution'],4.)
        self.assertEqual(by_name['initial_normalization']['loss_contribution'],'')

    def test_missing_duplicate_unknown_and_unavailable_rows_fail(self):
        for mode in ('missing','duplicate','unknown','unavailable','nan'):
            args=list(inputs())
            if mode=='missing': args[0]['target_rows'].pop()
            if mode=='duplicate': args[0]['target_rows'].append(args[0]['target_rows'][-1])
            if mode=='unknown': args[0]['target_rows'][0]['restriction_id']='invented'
            if mode in ('unavailable','nan'):
                args[1]['uniform_birth_time']['moments']['mean_children_ever_born_capped3_age25']=None if mode=='unavailable' else float('nan')
            with self.assertRaises((ValueError,RuntimeError),msg=mode): score_targets(*args)

    def test_no_target_is_silently_unweighted_and_normalization_not_double_scored(self):
        for index,weight in [(1,None),(0,1.),(1,0.),(1,-1.),(1,float('nan'))]:
            args=list(inputs()); args[0]['target_rows'][index]['actual_weight']=weight
            with self.assertRaises(ValueError): score_targets(*args)

    def test_signed_estate_funding_account(self):
        self.assertEqual(run_self_tests()['status'],'passed')

    def test_only_nonpositive_final_benefit_is_candidate_rejection(self):
        audit=dict(stationary_solves=3,psi_child=-.01)
        with self.assertRaises(NonpositiveNormalizedBenefit) as raised:
            validate_normalization(audit,-.01,3)
        self.assertEqual(raised.exception.audit['normalization'],audit)
        validate_normalization(dict(stationary_solves=3,psi_child=.01),.01,3)
        with self.assertRaises(RuntimeError) as bad_count:
            validate_normalization(audit,-.01,2)
        self.assertIsNot(type(bad_count.exception),NonpositiveNormalizedBenefit)


if __name__=='__main__': unittest.main()
