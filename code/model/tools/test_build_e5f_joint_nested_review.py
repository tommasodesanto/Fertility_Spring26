"""Read-only report validation against a complete collected scientific fixture."""
import argparse
import copy
import csv
import json
import tempfile
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).parent))
import build_e5f_joint_nested_review as report

class ReportValidationTests(unittest.TestCase):
    def setUp(self):
        self.selected = FIXTURE/'smoke/smoke_anchor/task_001'
        self.policy = FIXTURE/'parallel_policy_smoke'
        self.summary = report.read_json(self.selected/'summary.json')
        self.selected_sha = report.sha(self.selected/'summary.json')

    def test_complete_selected_and_policy_fixture(self):
        fit, pars, summary, receipt, checked = report.validate_selected(self.selected)
        self.assertEqual((len(fit),len(pars),checked),(12,14,21))
        _, checks, status = report.validate_policy(self.policy,self.selected_sha,self.summary)
        self.assertEqual(len(checks),8)
        self.assertIn('two-date smoke',status)

    def test_selected_rejects_false_receipts_and_missing_gate_evidence(self):
        original = report.read_json
        for mutation in ('status','missing_hash','nonfinite','missing_gate','probability'):
            def altered(path):
                obj=original(path)
                if path == self.selected/'case_receipt.json':
                    obj=copy.deepcopy(obj)
                    if mutation=='status':obj['status']='failed'
                    elif mutation=='missing_hash':del obj['artifact_sha256']['target_fit_long.csv']
                    elif mutation=='nonfinite':obj['gates']['market']['value']=float('nan')
                    elif mutation=='missing_gate':del obj['gates']['mass']
                    else:
                        bounds=next(iter(obj['policy_array_diagnostic']['probabilities'].values()))
                        bounds['maximum']=1.01
                return obj
            with self.subTest(mutation=mutation),patch.object(report,'read_json',side_effect=altered):
                with self.assertRaises(RuntimeError):report.validate_selected(self.selected)

    def test_policy_rejects_mixed_partial_and_invalid_evidence(self):
        original=report.read_json
        for mutation in ('partial','bundle','target','handoff','horizon','probability','nonfinite'):
            def altered(path):
                obj=copy.deepcopy(original(path))
                if path == self.policy/'equilibrium_receipt.json':
                    if mutation=='partial':obj['status']='partial_policy_failures'
                    elif mutation=='bundle':obj['scientific_bundle']='wrong'
                    elif mutation=='target':obj['target_fingerprint']='wrong'
                    elif mutation=='handoff':obj['inherited_state_verification_sha256']='wrong'
                    elif mutation=='horizon':obj['smoke']=False
                if path.name=='policy_array_summary.json' and mutation=='probability':
                    next(iter(obj['probabilities'].values()))['maximum']=1.01
                if path.name=='budget_summary.json' and mutation=='nonfinite':obj['budget_excess_mass']=float('nan')
                return obj
            with self.subTest(mutation=mutation),patch.object(report,'read_json',side_effect=altered):
                with self.assertRaises(RuntimeError):report.validate_policy(self.policy,self.selected_sha,self.summary)

    def test_full_writer_horizon_uses_endpoint_effects_and_all_dated_paths(self):
        # Synthetic interface fixture only: repeated data tests the writer's
        # 11-path-row / 2-effect-row schema, not any economic solution.
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            overall=copy.deepcopy(report.read_json(self.policy/'equilibrium_receipt.json'))
            overall['smoke']=False
            (root/'inherited_state_verification.json').write_bytes((self.policy/'inherited_state_verification.json').read_bytes())
            def write_rows(path,rows):
                with path.open('w',newline='') as stream:
                    writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
            for name in report.EXPECTED_POLICIES:
                folder=root/name;folder.mkdir()
                overall['cases'][name]['dates']=11
                (folder/'receipt.json').write_text(json.dumps(overall['cases'][name]))
                original=report.read_csv(self.policy/name/'policy_path.csv')
                rows=[]
                for year in range(2023,2064,4):
                    source_year=2023 if year==2023 else 2027
                    row=dict(original[0 if year==2023 else 1],calendar_year=str(year));rows.append(row)
                    (folder/f'date_{year}').symlink_to(self.policy/name/f'date_{source_year}',target_is_directory=True)
                write_rows(folder/'policy_path.csv',rows)
            (root/'equilibrium_receipt.json').write_text(json.dumps(overall))
            effects=report.read_csv(self.policy/'policy_effects.csv')
            for row in effects:
                if int(row['year'])==2027:row['year']='2063'
            write_rows(root/'policy_effects.csv',effects)
            _,checks,status=report.validate_policy(root,self.selected_sha,self.summary)
            self.assertEqual(len(checks),44)
            self.assertIn('44-date',status)
            path=root/'baseline/policy_path.csv';rows=report.read_csv(path);rows[-1]['calendar_year']='2064';write_rows(path,rows)
            with self.assertRaisesRegex(RuntimeError,'inconsistent dates'):
                report.validate_policy(root,self.selected_sha,self.summary)

    def test_benchmark_requires_complete_comparable_fit_and_parameters(self):
        root=FIXTURE.parents[1]/'e5f_joint_nested_experiment_20260906a'
        fit=report.read_csv(self.selected/'target_fit_long.csv')
        ref=report.validate_reference(root/'reference_target_fits.csv',root/'reference_parameters.csv',fit)
        self.assertAlmostEqual(ref['loss'],30.482966707698903,places=10)
        with self.assertRaisesRegex(RuntimeError,'complete fit and parameter'):
            report.validate_reference(root/'reference_target_fits.csv',None,fit)
        fit[0]['weight']=str(2*float(fit[0]['weight']))
        with self.assertRaisesRegex(RuntimeError,'target or weight differs'):
            report.validate_reference(root/'reference_target_fits.csv',root/'reference_parameters.csv',fit)

    def test_policy_rejects_duplicate_effect_rows(self):
        original=report.read_csv
        def altered(path):
            rows=original(path)
            return rows+[rows[0]] if path.name=='policy_effects.csv' else rows
        with patch.object(report,'read_csv',side_effect=altered):
            with self.assertRaisesRegex(RuntimeError,'incomplete or duplicated'):
                report.validate_policy(self.policy,self.selected_sha,self.summary)

    def _partial_policy_fixture(self, root):
        """Copied synthetic evidence; the bundled scientific fixture stays read-only."""
        overall=copy.deepcopy(report.read_json(self.policy/'equilibrium_receipt.json'))
        overall['smoke']=False
        overall['status']='partial_policy_failures'
        overall['failures']={'property-tax-2pct-no-rebate': {'error': 'Housing market did not clear: residual=2.791e-04.'}}
        (root/'inherited_state_verification.json').write_bytes((self.policy/'inherited_state_verification.json').read_bytes())
        def write_rows(path,rows):
            with path.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        for name in report.EXPECTED_POLICIES:
            folder=root/name;folder.mkdir()
            original=report.read_csv(self.policy/name/'policy_path.csv')
            if name == 'property-tax-2pct-no-rebate':
                rows=[]
                for year in range(2023,2048,4):
                    row=dict(original[0 if year==2023 else 1],calendar_year=str(year),policy_case=name,
                             relative_market_residual='0.00001',mass_accounting_residual='0',nonfinite_distribution_count='0')
                    rows.append(row)
                    (folder/f'date_{year}').symlink_to(self.policy/name/f'date_{2023 if year==2023 else 2027}',target_is_directory=True)
                write_rows(folder/'policy_path_progress.csv',rows)
                (folder/'failure.json').write_text(json.dumps(overall['failures'][name]))
                continue
            overall['cases'][name]['dates']=11
            (folder/'receipt.json').write_text(json.dumps(overall['cases'][name]))
            rows=[]
            for year in range(2023,2064,4):
                source_year=2023 if year==2023 else 2027
                rows.append(dict(original[0 if year==2023 else 1],calendar_year=str(year)))
                (folder/f'date_{year}').symlink_to(self.policy/name/f'date_{source_year}',target_is_directory=True)
            write_rows(folder/'policy_path.csv',rows)
        del overall['cases']['property-tax-2pct-no-rebate']
        (root/'equilibrium_receipt.json').write_text(json.dumps(overall))
        effects=[row for row in report.read_csv(self.policy/'policy_effects.csv') if row['policy'] != 'property-tax-2pct-no-rebate']
        for row in effects:
            if int(row['year'])==2027:row['year']='2063'
        write_rows(root/'policy_effects.csv',effects)

    def test_partial_policy_requires_opt_in_and_validates_failure_prefix_and_gates(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);self._partial_policy_fixture(root)
            with self.assertRaisesRegex(RuntimeError,'allow-partial-policies'):
                report.validate_policy(root,self.selected_sha,self.summary)
            _,checks,status=report.validate_policy(root,self.selected_sha,self.summary,allow_partial_policies=True)
            self.assertEqual(len(checks),40);self.assertIn('33 certified full-branch dates plus 7 valid property-tax prefix dates, 40 total, 4 unavailable',status)
            for mutation, expected in (('missing_failure','failure.json'),('incorrect_failure','does not match'),('hash','different selected'),('dated_gate','dated gate'),('prefix','strict consecutive')):
                with self.subTest(mutation=mutation):
                    with tempfile.TemporaryDirectory() as altered_temp:
                        altered=Path(altered_temp);self._partial_policy_fixture(altered)
                        if mutation=='missing_failure': (altered/'property-tax-2pct-no-rebate/failure.json').unlink()
                        elif mutation=='incorrect_failure': (altered/'property-tax-2pct-no-rebate/failure.json').write_text('{}')
                        elif mutation=='hash':
                            receipt=report.read_json(altered/'equilibrium_receipt.json');receipt['selected_summary_sha256']='wrong';(altered/'equilibrium_receipt.json').write_text(json.dumps(receipt))
                        else:
                            path=altered/'property-tax-2pct-no-rebate/policy_path_progress.csv';rows=report.read_csv(path)
                            if mutation=='dated_gate':rows[0]['relative_market_residual']='0.0003'
                            else:rows[-1]['calendar_year']='2051'
                            with path.open('w',newline='') as stream:
                                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
                        with self.assertRaisesRegex(RuntimeError,expected):
                            report.validate_policy(altered,self.selected_sha,self.summary,allow_partial_policies=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--fixture-root',type=Path,required=True)
    args=parser.parse_args();FIXTURE=args.fixture_root.resolve()
    unittest.main(argv=[sys.argv[0]])
