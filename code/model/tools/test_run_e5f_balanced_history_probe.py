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


if __name__=='__main__':unittest.main()
