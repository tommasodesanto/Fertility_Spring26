"""Complete scores for previously verified saved cases; no model solves or promotion."""
import csv
import hashlib
import json
import math
from pathlib import Path
from score_initial import fingerprint, score_initial

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
CONTRACT_HASH = 'c0e266d3a0d430343c469d780d1aedb45fa87f8763c9c938889e0c37daa31de2'
CONTRACT_BYTES = 'e0bd8316a19bb197ab0fe9adaf25cb3173ff3b4072c14515cdd4e34b516ff43c'

def read(p): return json.loads(p.read_text())
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p, x): p.write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')
def rows(p):
    with p.open() as f: return list(csv.DictReader(f))

def main():
    path = HERE / 'working_contract.json'
    assert digest(path) == CONTRACT_BYTES
    contract = read(path)
    assert fingerprint(contract) == CONTRACT_HASH
    for name in ('scorer', 'target_fit', 'working_weights', 'lead_decision',
                 'parameter_table', 'fertility_provenance', 'housing_wealth_provenance'):
        item = contract['input_artifacts'][name]
        assert digest(Path(item['path'])) == item['sha256'] == contract['source_fingerprints'][name + '_file_sha256']
    panel_dir = BASE / 'initial_sensitivity_panel'
    panel = read(panel_dir / 'raw_case_summary.json')
    assert digest(panel_dir / 'raw_case_summary.json') == read(panel_dir / 'panel_receipt.json')['raw_case_summary_sha256']
    baseline = next(c for c in panel['cases'] if c['case_id'] == 'baseline')
    joint_dir = BASE / 'initial_joint_round_01/smoke_readout'
    joint_raw = read(joint_dir / 'raw_case_summary.json')
    assert digest(joint_dir / 'raw_case_summary.json') == read(joint_dir / 'array_collection_receipt.json')['raw_case_summary_sha256']
    joint = next(c for c in joint_raw['cases'] if c['smoke_reuse'])
    recent_dir = BASE / 'initial_fit_readout/recent_parent_probe/completed_17362130'
    recent_panel = read(recent_dir / 'recent_parent_panel.json')
    recent_baseline = next(c for c in recent_panel['rows'] if c['case_id'] == 'baseline')
    recent_joint = read(recent_dir / 'recent_parent_joint.json')['row']
    cases = [
        ('baseline', baseline, recent_baseline, panel_dir / 'collected/baseline/repetition_01', 'numerical_gates_verified'),
        ('joint_smoke', joint, recent_joint, joint_dir / 'array_collected' / joint['case_id'] / 'repetition_02', 'numeric_gates_verified'),
    ]
    out = HERE / 'saved_case_scores'
    out.mkdir(exist_ok=True)
    results = []
    for name, raw, recent_record, directory, gate_key in cases:
        assert raw[gate_key] is True
        assert raw['checkpoint']['sha256'] == recent_record['checkpoint_sha256']
        recent_path = Path(recent_record['observation_path'])
        assert digest(recent_path) == recent_record['observation_sha256']
        recent = read(recent_path)
        early = read(directory / 'early_measurement.json')
        assert early == raw['early_measurement']
        summary = read(directory / 'summary.json')
        assert summary['checkpoint_sha256'] == raw['checkpoint']['sha256']
        parameters = rows(directory / 'parameters.csv')
        assert len(parameters) == 17
        for parameter in parameters:
            assert float(parameter['estimate']) == raw['parameters'][parameter['parameter']]
        inputs = dict(early_measurement=early, recent_parent_observation=recent,
                      normalization=summary['normalization'], parameters=parameters)
        receipt = dict(schema='e5f_initial_score_receipt_v1', status='verified',
            source_fingerprints=contract['source_fingerprints'],
            checkpoint_sha256=raw['checkpoint']['sha256'], numerical_gates_verified=True,
            recent_parent_certified=True,
            recent_parent_approximation_id=contract['recent_parent_approximation']['approximation_id'],
            input_sha256={k: fingerprint(v) for k, v in inputs.items()},
            verification_basis='Previously verified cluster collection and checked saved observation hashes; remote bytes not re-read during current authentication outage.',
            initial_collection=raw, recent_observation_collection=recent_record)
        result = score_initial(contract, expected_contract_sha256=CONTRACT_HASH,
                               evaluation_receipt=receipt, **inputs)
        # Independent direct calculation from the saved fitted rows.
        independently_summed = math.fsum((r['gap']/r['working_scale'])**2
                                        for r in result['target_fit'] if r['scored'])
        assert math.isclose(result['loss'], independently_summed, rel_tol=1e-14)
        save(out / (name + '_receipt.json'), receipt)
        save(out / (name + '_score.json'), result)
        for suffix, table in [('fit', result['target_fit']), ('parameters', result['parameters'])]:
            with (out / (name + '_' + suffix + '.csv')).open('w', newline='') as f:
                fields = list(dict.fromkeys(k for row in table for k in row))
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader(); writer.writerows(table)
        results.append(dict(case=name, loss=result['loss'],
                            recent_parent_model=recent['model_value'],
                            recent_parent_loss=next(r['loss_contribution'] for r in result['target_fit']
                                                   if r['restriction_id']=='recent_parent_ownership')))
    saved_panel = read(HERE / 'complete_panel_analysis.json')
    panel_baseline = next(c for c in saved_panel['cases'] if c['case_id']=='baseline')
    assert results[0]['loss'] == panel_baseline['complete_12_moment_loss']
    save(out / 'comparison.json', dict(contract_sha256=CONTRACT_HASH, cases=results,
        source_packets_unchanged=True, new_model_solves=0, production_benchmark_selected=False))
    print(json.dumps(results))

if __name__ == '__main__': main()
