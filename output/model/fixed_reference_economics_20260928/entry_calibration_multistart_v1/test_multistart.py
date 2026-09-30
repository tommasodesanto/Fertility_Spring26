"""Extra-start identity and provenance checks; no native lifecycle solves."""
import hashlib,json,unittest
from pathlib import Path
import inputs,runner


class MultistartChecks(unittest.TestCase):
    def test_dispatch_exactly_nine_distinct_extra_starts(self):
        dispatch=runner.PLAN['dispatch_lanes']
        self.assertEqual(len(dispatch),9)
        self.assertEqual(len(set(dispatch)),9)
        self.assertEqual(set(dispatch),set(inputs.LANES)-set(runner.PLAN['existing_lanes_not_dispatched']))
        self.assertEqual(runner.PLAN['common_deadline_iso'],'2026-10-01T02:38:17Z')
        self.assertEqual(runner.PLAN['common_deadline_epoch'],1790822297.0)
        for base in runner.PLAN['existing_lanes_not_dispatched']:
            seeds=[tuple(inputs.LANES[base+'_'+s]['seed'][k] for k in inputs.PARAMETERS) for s in ['s1','s2','s3']]
            self.assertEqual(len(set(seeds)),3)
            self.assertTrue(all(seed!=tuple(inputs.LANES[base]['seed'][k] for k in inputs.PARAMETERS) for seed in seeds))

    def test_offsets_bounds_economic_contract_and_parent_proof(self):
        for lane in runner.PLAN['dispatch_lanes']:
            cfg=inputs.LANES[lane];parent=inputs.LANES[cfg['parent_lane']]
            seed,bounds,_=inputs.seed_and_bounds(lane);inputs.check_point(seed,bounds)
            proof=cfg['starting_point_construction']
            self.assertFalse(proof['native_verified']);self.assertEqual(proof['starting_point_full_ge_repeats'],0)
            self.assertIn('parent pilot',cfg['seed_provenance']['refers_to'])
            self.assertEqual(proof['parent_seed'],parent['seed'])
            self.assertEqual(proof['step_sizes'],runner.step_sizes(parent['seed'],bounds))
            self.assertEqual(proof['parent_plan_sha256'],runner.PLAN['parent_plan_sha256'])
            self.assertEqual(cfg['parent_seed_parameter_table'],parent['seed_parameter_table'])
            for key in inputs.PARAMETERS:
                requested=3*proof['step_sizes'][key]*proof['direction'][list(inputs.PARAMETERS).index(key)]
                self.assertEqual(proof['requested_offsets'][key],requested)
                self.assertEqual(seed[key],min(bounds[key][1],max(bounds[key][0],parent['seed'][key]+requested)))
                self.assertEqual(proof['effective_offsets'][key],seed[key]-parent['seed'][key])
            for key in ['arm','dimensions','size']:
                self.assertEqual(cfg[key],parent[key])
            self.assertEqual(runner.PLAN['arms'][lane],runner.PLAN['arms'][cfg['parent_lane']])

    def test_same_nonnegative_starts_across_grids(self):
        for suffix in ['s1','s2','s3']:
            coarse=inputs.LANES['nonnegative_mean_120x9_'+suffix]
            fine=inputs.LANES['nonnegative_mean_160x15_'+suffix]
            self.assertEqual(coarse['seed'],fine['seed'])
            self.assertEqual(coarse['starting_point_construction'],fine['starting_point_construction'])

    def test_unchanged_runners_and_sources(self):
        here=Path(__file__).resolve().parent;parent=here.with_name('entry_calibration_round1_v1')
        for name in ['runner.py','inputs.py','phase_b_pilot.py','test_round1.py']:
            self.assertEqual((here/name).read_bytes(),(parent/name).read_bytes())
        source=json.loads((parent/'plan.json').read_text())
        for key in ['target_contract','fixed','annual_interest_rate','bundle_sources','maximum_full_ge','maximum_rounds','per_lifecycle_seconds','maximum_lifecycle_per_ge','selected_full_ge_repeats','final_reserve_seconds','minimum_ge_seconds','method']:
            self.assertEqual(runner.PLAN[key],source[key])
        runner.verify_sources()


if __name__=='__main__':unittest.main()
