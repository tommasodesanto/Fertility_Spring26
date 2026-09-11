"""Pure contract/source/parent receipt tests; no model evaluations."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import run_e5f_balanced_history_probe as driver


class HistoryProbeContractTest(unittest.TestCase):
    def setUp(self):
        pin=lambda p,h='1'*64:dict(path=p,sha256=h)
        self.c=dict(schema=driver.SCHEMA,source_sha256={},originating_source_commit='c6dd3508',count=6,
            initial_checkpoint=pin('/initial/new_balanced_smoke/repetition_02/initial_state.pkl.gz',driver.terminal_driver.INITIAL_SHA256),
            initial_summary=pin('/initial/new_balanced_smoke/repetition_02/summary.json'),
            initial_contract=pin('/initial/input_contract.json'),terminal_checkpoint=pin('/terminal/terminal_state.pkl.gz'),
            terminal_summary=pin('/terminal/summary.json'),terminal_contract=pin('/terminal_input.json'),
            terminal_root_receipt=pin('/terminal/root_receipt.json'),outside_origin_entry_share=.169,
            outside_entry_status='diagnostic_outstanding_not_estimated',psi_change_from_initial=-.25,
            preference_rule='announced_linear_2007_2023_then_constant_diagnostic',
            demographic_seed_mode='serialized_frozen_2023_primitives_without_realignment',terminal_demographic_year=2100,
            seconds=1800,root_seconds=1500,identity_tolerance=2e-9,initial_fertility_tolerance=5e-4,
            audit_controls=dict(reconstruction_tolerance=5e-9,feasibility_projection_tolerance=1e-6,
                probability_tolerance=1e-12,occupied_mass_tolerance=1e-12,value_drop_tolerance=1e-7),
            root_controls=dict(price_bounds=[.1,2.],pension_bounds=[.1,4.],market_tolerance=2e-4,
                fiscal_tolerance=1e-6,market_slope=1.63,fiscal_slope=1.,max_log_step=.1,damping=1.,
                max_evaluations=2,max_condition_number=1e10,worsening_factor=1.5,
                final_reproduction_tolerance=2e-10,initial_jacobian=None),
            initial_prices=[.6]*6,initial_pensions=[2.]*6,standard_graph_count=17)
        for name in ('e5f_parenthood_utility.py','e5f_social_security.py','e5f_balanced_terminal.py',
            'e5f_matched_pf_endpoint.py','e5f_social_security_root.py','e5f_matched_pf_path_root.py',
            'run_dynamic_population_transition.py','run_e5f_open_population_transition.py',
            'run_e5f_perfect_foresight_transition.py','run_e5f_perfect_foresight_person_demography.py'):
            self.c['source_sha256']['code/model/tools/'+name]='a'*64

    def test_explicit_six_date_replay_contract_passes_without_mutation(self):
        before=copy.deepcopy(self.c)
        driver.validate_contract(self.c)
        self.assertEqual(before,self.c)

    def test_wrong_timing_budget_entry_or_preference_rejected(self):
        for key,value in [('count',5),('seconds',1801),('root_seconds',1700),
                          ('outside_origin_entry_share',0.),('psi_change_from_initial',-.3),
                          ('outside_entry_status','estimated'),('terminal_demographic_year',2023),
                          ('standard_graph_count',16)]:
            c=copy.deepcopy(self.c);c[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):driver.validate_contract(c)

    def test_no_implicit_guesses_or_optimizer_step(self):
        for key,value in [('initial_prices',[.6]*5),('initial_pensions',[float('nan')]*6),
                          ('initial_prices',[10.]*6)]:
            c=copy.deepcopy(self.c);c[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):driver.validate_contract(c)
        for key,value in [('max_evaluations',3),('initial_jacobian',[[1.]])]:
            c=copy.deepcopy(self.c);c['root_controls'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):driver.validate_contract(c)

    def test_pending_hash_wrong_initial_or_cross_packet_receipts_rejected(self):
        for name,key,value in [('terminal_checkpoint','sha256','PENDING'),
            ('initial_checkpoint','sha256','0'*64),('terminal_summary','path','/another/summary.json'),
            ('terminal_root_receipt','path','relative/root_receipt.json')]:
            c=copy.deepcopy(self.c);c[name][key]=value
            with self.subTest(name=name),self.assertRaises(ValueError):driver.validate_contract(c)

    def test_source_manifest_requires_complete_unchanged_source(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d).resolve();source=root/'code/model/tools/probe.py';source.parent.mkdir(parents=True)
            source.write_text('x=1\n');relative=str(source.relative_to(root))
            c={'source_sha256':{relative:driver.digest(source)}}
            with mock.patch.object(driver.terminal_driver,'ROOT',root):
                self.assertEqual(driver.terminal_driver.verify_sources(c),1)
                source.write_text('x=2\n')
                with self.assertRaisesRegex(ValueError,'SHA256'):driver.terminal_driver.verify_sources(c)
                with self.assertRaisesRegex(ValueError,'every'):driver.terminal_driver.verify_sources({'source_sha256':{}})

    def test_hash_is_checked_before_pinned_json_is_read(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'receipt.json';p.write_text('{}')
            self.assertEqual(driver.pinned_json(dict(path=str(p),sha256=driver.digest(p))),{})
            with self.assertRaisesRegex(ValueError,'SHA256'):
                driver.pinned_json(dict(path=str(p),sha256='0'*64))

    def parents(self):
        c=self.c
        initial={'contract_sha256':c['initial_contract']['sha256']}
        initial_summary=dict(status='passed_initial_diagnostic',checkpoint_sha256=c['initial_checkpoint']['sha256'])
        initial_contract=dict(schema='e5f_parenthood_initial_probe_v1',normalize=True)
        terminal=dict(schema='e5f_balanced_terminal_probe_checkpoint_v1',contract_sha256=c['terminal_contract']['sha256'])
        terminal_summary=dict(status='passed_terminal_root_diagnostic',endpoint_numerically_verified=True,
            final_checks={'all':True},checkpoint_sha256=c['terminal_checkpoint']['sha256'])
        terminal_contract={k:copy.deepcopy(c[k]) for k in ('initial_checkpoint','initial_summary','initial_contract',
            'psi_change_from_initial','terminal_demographic_year','demographic_seed_mode')}
        terminal_contract['schema']=driver.terminal_driver.SCHEMA
        terminal_contract['source_sha256']=copy.deepcopy(c['source_sha256'])
        root=dict(schema='e5f_balanced_terminal_v1',converged=True,endpoint_production_eligible=True,
            fresh_endpoint_matches_final=True,gates={'all':True})
        return [initial,initial_summary,initial_contract,terminal,terminal_summary,terminal_contract,root]

    def test_original_initial_and_terminal_receipt_chain(self):
        compared=driver.validate_parents(self.c,*self.parents())
        self.assertEqual(compared,[dict(path=k,sha256=v) for k,v in sorted(self.c['source_sha256'].items())])
        for index,key,value in [(0,'contract_sha256','bad'),(2,'normalize',False),
            (3,'contract_sha256','bad'),(4,'status','failed_terminal_root_diagnostic'),
            (5,'psi_change_from_initial',-.1),(6,'converged',False),(6,'fresh_endpoint_matches_final',False)]:
            parents=self.parents();parents[index][key]=value
            with self.subTest(index=index,key=key),self.assertRaises(ValueError):
                driver.validate_parents(self.c,*parents)

    def test_terminal_cannot_have_another_initial_parent(self):
        parents=self.parents();parents[5]['initial_contract']['sha256']='bad'
        with self.assertRaisesRegex(ValueError,'another normalized initial'):
            driver.validate_parents(self.c,*parents)

    def test_changed_terminal_economic_source_is_rejected(self):
        parents=self.parents()
        parents[5]['source_sha256']['code/model/tools/e5f_social_security.py']='b'*64
        with self.assertRaisesRegex(ValueError,'economic source'):
            driver.validate_parents(self.c,*parents)

    def test_omitted_tools_source_cannot_escape_terminal_chain(self):
        parents=self.parents()
        parents[5]['source_sha256']['code/model/tools/some_other_economic_helper.py']='c'*64
        with self.assertRaisesRegex(ValueError,'economic source'):
            driver.validate_parents(self.c,*parents)


class ActualRootContractTest(unittest.TestCase):
    setUp=HistoryProbeContractTest.setUp
    def root_contract(self,count=6):
        c=copy.deepcopy(self.c);c.update(schema=driver.ROOT_SCHEMA,mode='actual_root',count=count,
            initial_prices=[.6]*count,initial_pensions=[2.]*count)
        c['root_controls']['max_evaluations']=8
        c['root_seconds']=1620 if count==6 else 6600
        c['seconds']=1800 if count==6 else 7200
        return c

    def test_explicit_actual_root_six_and_long_contracts(self):
        for count in (6,28):
            c=self.root_contract(count);before=copy.deepcopy(c)
            driver.validate_contract(c);self.assertEqual(c,before)

    def test_root_mode_count_budget_and_incomplete_vectors_rejected(self):
        for key,value in [('mode','replay'),('count',5),('count',100),('seconds',1801),
                          ('root_seconds',1621),('initial_pensions',[2.]*5)]:
            c=self.root_contract();c[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):driver.validate_contract(c)
        c=self.root_contract();del c['mode']
        with self.assertRaises(ValueError):driver.validate_contract(c)
        c=self.root_contract(28);c['seconds']=7201
        with self.assertRaises(ValueError):driver.validate_contract(c)
        for budget in (2,9):
            c=self.root_contract();c['root_controls']['max_evaluations']=budget
            with self.assertRaises(ValueError):driver.validate_contract(c)

    def replay_fixture(self):
        c=self.root_contract();count=c['count']
        # Trial2 is chosen; trial3 worsens; trial4 freshly repeats2.
        prices=[.6,.7,.8,.7];pensions=[2.,2.2,2.4,2.2]
        signatures=[];records=[]
        for i,(p,b) in enumerate(zip(prices,pensions),1):
            signatures.append([dict(year=y,policy={'V':str(p)},distributions={'g_pre':str(b)},
                fiscal={'pension':b}) for y in range(2007,2031,4)])
            records.append(dict(evaluation=i,phase='final' if i==4 else ('initial' if i==1 else 'iterate'),
                prices=[p]*count,fiscal_values=[b]*count,market_residual=[p-.72]*count,
                fiscal_residual=[b-2.22]*count,bellman_solves=2*count))
        def point(i):
            x=copy.deepcopy(records[i-1]);x.update(mapping_valid=True,payload=dict(trial=i,bellman_solves=2*count,
                mapping_gates={'mass':{'passed':True}},
                dated_household_audits=[{'gates':{'household':True}} for _ in range(count)]))
            return x
        receipt=dict(evaluations=4,best=point(2),final=point(4),history=copy.deepcopy(records),
            fresh_path_matches_final=True,gates=dict(mapping=True,market_replay=True,fiscal_replay=True,
                housing=False,social_security=False))
        return c,receipt,signatures,records

    def test_final_replay_uses_selected_trial_not_initial_or_last_iterate(self):
        c,r,s,records=self.replay_fixture()
        got=driver.validate_replay(c,r,s,records)
        self.assertEqual(got,dict(verified=True,selected_trial=2,fresh_final_trial=4,
            completed_mappings=4,bellman_solves=48))
        # Residual gates remain false: a valid replay must not claim equilibrium.
        self.assertFalse(r['gates']['housing']);self.assertFalse(r['gates']['social_security'])

    def test_initial_instead_of_selected_signature_rejected(self):
        c,r,s,records=self.replay_fixture();s[-1]=copy.deepcopy(s[0])
        with self.assertRaisesRegex(RuntimeError,'whole-path'):driver.validate_replay(c,r,s,records)

    def test_wrong_selected_trial_and_missing_final_rejected(self):
        for change in ('trial','final','association','phase'):
            c,r,s,records=self.replay_fixture()
            if change=='trial':r['best']['payload']['trial']=1
            elif change=='final':r['final']=None
            elif change=='association':r['fresh_path_matches_final']=False
            else:records[-1]['phase']='iterate'
            with self.subTest(change=change),self.assertRaises(RuntimeError):driver.validate_replay(c,r,s,records)

    def test_missing_or_duplicate_mapping_and_dates_rejected(self):
        for change in ('missing','duplicate','date','bellman'):
            c,r,s,records=self.replay_fixture()
            if change=='missing':records.pop()
            elif change=='duplicate':records[-1]['evaluation']=3
            elif change=='date':s[1].pop()
            else:records[1]['bellman_solves']=10
            with self.subTest(change=change),self.assertRaises(RuntimeError):driver.validate_replay(c,r,s,records)

    def test_coordinate_residual_household_or_fiscal_drift_rejected(self):
        for change in ('coordinates','residual','audit','fiscal','nonfinite'):
            c,r,s,records=self.replay_fixture()
            if change=='coordinates':r['final']['fiscal_values'][0]=3.
            elif change=='residual':r['final']['fiscal_residual'][0]=.1
            elif change=='audit':r['best']['payload']['dated_household_audits'][0]['gates']['household']=False
            elif change=='fiscal':s[-1][0]['fiscal']['pension']=3.
            else:
                for point in (r['best'],r['final'],records[1],records[-1],r['history'][1],r['history'][-1]):
                    point['market_residual'][0]=float('inf')
            with self.subTest(change=change),self.assertRaises(RuntimeError):driver.validate_replay(c,r,s,records)

    def test_original_two_mapping_receipt_still_passes(self):
        c,r,s,records=self.replay_fixture();c=copy.deepcopy(self.c)
        records=[records[1],records[-1]]
        for i,row in enumerate(records,1):row['evaluation']=i;row['phase']='initial' if i==1 else 'final'
        r['evaluations']=2;r['history']=copy.deepcopy(records)
        r['best']['payload']['trial']=1;r['final']['payload']['trial']=2
        got=driver.validate_replay(c,r,[s[1],s[-1]],records)
        self.assertEqual(got['bellman_solves'],24)


class ContinuationContractTest(unittest.TestCase):
    def setUp(self):
        fixture=ActualRootContractTest();fixture.setUp()
        previous,root,signatures,records=fixture.replay_fixture()
        self.directory=tempfile.TemporaryDirectory();self.addCleanup(self.directory.cleanup)
        folder=Path(self.directory.name)
        root.update(schema='e5f_balanced_history_v1',converged=False,
            finite_horizon_market_fiscal_converged=False,years=list(range(2007,2028,4)),
            closure='fixed_tax',payroll_tax=.179,property_tax_period=.04,
            property_tax_rebate=0.,horizon_verified=False,final_damping=.5,
            final_jacobian=[[float(i==j) for j in range(12)] for i in range(12)])
        checkpoint=folder/'dated_2023.pkl.gz';checkpoint.write_bytes(b'opaque checkpoint evidence only')
        summary=dict(status='incomplete_finite_history_root_diagnostic',
            finite_horizon_market_fiscal_converged=False,mapping_replay_verified=True,
            checkpoint_reload_verified=True,standard_graph_count=17,
            checkpoint_sha256=driver.digest(checkpoint),root_evaluations=4)
        objects=dict(contract=dict(previous,contract_sha256='a'*64),root_receipt=root,
            summary=summary,dated_reproduction=dict(signatures=signatures),root_evaluations=records)
        for name,obj in objects.items():
            (folder/driver.CONTINUATION_FILES[name]).write_text(json.dumps(obj))
        self.c=copy.deepcopy(previous)
        self.c.update(schema=driver.CONTINUATION_SCHEMA,
            initial_prices=root['best']['prices'],initial_pensions=root['best']['fiscal_values'],
            continuation={name:dict(path=str(folder/file),sha256=driver.digest(folder/file))
                          for name,file in driver.CONTINUATION_FILES.items()})
        self.c['root_controls'].update(initial_jacobian=root['final_jacobian'],damping=.5)
        self.root=root

    def test_verified_restart_keeps_original_initial_state_and_physical_jacobian(self):
        before=copy.deepcopy(self.c);driver.validate_contract(self.c)
        result=driver.load_continuation(self.c)
        self.assertEqual(self.c,before)
        self.assertEqual(result['prior_final'],self.root['final'])
        self.assertIn('original2007 initial input remains unchanged',result['checkpoint_role'])
        self.assertNotEqual(self.c['initial_checkpoint']['path'],self.c['continuation']['checkpoint']['path'])

    def test_invalid_or_unpinned_matrix_rejected(self):
        for matrix in ([[1.]],[[float('nan')]*12 for _ in range(12)],
                       [[True]*12 for _ in range(12)],[["1"]*12 for _ in range(12)]):
            c=copy.deepcopy(self.c);c['root_controls']['initial_jacobian']=matrix
            with self.subTest(matrix=str(matrix)[:20]),self.assertRaises(ValueError):driver.validate_contract(c)
        c=copy.deepcopy(self.c);c['schema']=driver.ROOT_SCHEMA;del c['continuation']
        with self.assertRaises(ValueError):driver.validate_contract(c)

    def test_changed_science_point_jacobian_damping_or_gate_rejected(self):
        for change in ('source','entry','price','pension','jacobian','damping','tolerance'):
            c=copy.deepcopy(self.c)
            if change=='source':c['source_sha256']['code/model/tools/e5f_social_security.py']='b'*64
            elif change=='entry':c['outside_origin_entry_share']=.17
            elif change=='price':c['initial_prices'][0]+=.01
            elif change=='pension':c['initial_pensions'][0]+=.01
            elif change=='jacobian':c['root_controls']['initial_jacobian'][0][0]=2.
            elif change=='damping':c['root_controls']['damping']=.6
            else:c['root_controls']['market_tolerance']=.001
            with self.subTest(change=change),self.assertRaises(ValueError):driver.load_continuation(c)

    def test_prior_hash_drift_rejected_before_loading(self):
        c=copy.deepcopy(self.c);c['continuation']['summary']['sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'SHA256'):driver.load_continuation(c)

    def test_first_mapping_replays_before_a_second_update_can_run(self):
        continuation=driver.load_continuation(self.c)
        record=copy.deepcopy(self.root['final']);record.update(evaluation=1,phase='initial')
        signatures=continuation['prior_signatures']
        passed=driver.verify_continuation_initial(continuation,record,signatures,2e-10)
        self.assertTrue(passed['verified']);self.assertTrue(passed['initial_mapping_counts_in_budget'])
        for change in ('coordinates','residual','signature','phase','invalid'):
            r=copy.deepcopy(record);s=copy.deepcopy(signatures)
            if change=='coordinates':r['prices'][0]+=.01
            elif change=='residual':r['market_residual'][0]+=.001
            elif change=='signature':s[0]['policy']['V']='changed'
            elif change=='phase':r['evaluation']=2
            else:r['mapping_valid']=False
            calls=[]
            def next_update():
                driver.verify_continuation_initial(continuation,r,s,2e-10)
                calls.append('unsafe next update')
            with self.subTest(change=change),self.assertRaises(RuntimeError):next_update()
            self.assertEqual(calls,[])


if __name__=='__main__':unittest.main()
