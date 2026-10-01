"""Render saved native dated packets without rebuilding or solving the model."""
from __future__ import annotations
import argparse
import contextlib
import gzip
import hashlib
import json
import pickle
from pathlib import Path

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def pinned(path, digest):
    path = Path(path)
    require(sha(path) == digest, 'Input hash differs: ' + str(path))
    return path

def load_record(path, digest):
    record = json.loads(pinned(path, digest).read_text())
    rows = record['rows']
    require(len(rows) > 0, 'Empty native dated mapping')
    require([r['period'] for r in record['diagnostic_packets']] == sorted({0, len(rows)//2, len(rows)-1}),
            'Original first/middle/last diagnostic dates required')
    for item in record['diagnostic_packets']:
        require(set(item) == {'path', 'sha256', 'period', 'dated_rent'}, 'Unexpected native packet receipt fields')
        pinned(item['path'], item['sha256'])
    require(len(record['market_residual']) == len(rows) == len(record['fiscal_residual']),
            'Native residual dates differ')
    require(isinstance(record['gates'], dict) and bool(record['gates']), 'Native gates missing')
    return record

@contextlib.contextmanager
def forbid_solves(runtime):
    def denied(*args, **kwargs):
        raise RuntimeError('A native solve is forbidden during saved-packet reporting')
    originals = {name: getattr(runtime.model, name) for name in dir(runtime.model)
                 if name.startswith('solve_') and callable(getattr(runtime.model, name))}
    require(bool(originals), 'Native solver facade required for zero-call guard')
    try:
        for name in originals:
            setattr(runtime.model, name, denied)
        yield
    finally:
        for name, value in originals.items():
            setattr(runtime.model, name, value)

def verify_outputs(output, rendered, names, packets, rows):
    require(len(names) == 17 and len(set(names)) == 17, 'Exactly 17 original plot names required')
    require([r['period'] for r in rendered['sampled_dates']] == [p['period'] for p in packets],
            'Rendered native dates differ')
    retained = []
    for item in rendered['sampled_dates']:
        period = item['period']
        folder = Path(output)/f'date_{period:03d}'/'standard_diagnostics'
        files = sorted(folder.glob('*.png'))
        require(len(files) == 17 and {p.name for p in files} == set(names), 'Rendered standard 17-plot set changed')
        require(set(item['plots']) == set(names), 'Renderer receipt plot set changed')
        require(all(sha(p) == item['plots'][p.name] for p in files), 'Rendered PNG hash differs')
        retained.append(dict(period=period, calendar_year=rows[period]['calendar_year'],
                             plots={p.name: dict(path=str(p), sha256=sha(p)) for p in files}))
    require({p.name for p in Path(output).glob('date_*')} == {f'date_{p["period"]:03d}' for p in packets},
            'Unexpected diagnostic output dates')
    return retained

def render(*, handoff, handoff_sha256, mapping_record, mapping_record_sha256, output,
           run_identity=None, run_identity_sha256=None, runtime_class=None):
    output = Path(output)
    require(not output.exists(), 'Fresh report output directory required')
    # Authenticate every serialized input before any native object is opened.
    pinned(handoff, handoff_sha256)
    record = load_record(mapping_record, mapping_record_sha256)
    require((run_identity is None) == (run_identity_sha256 is None), 'Run identity path and hash must be paired')
    plan = None
    if run_identity is not None:
        plan = json.loads(pinned(run_identity, run_identity_sha256).read_text())
        for pin in plan.get('source_files', {}).values():
            pinned(pin['path'], pin['sha256'])
    if runtime_class is None:
        from floor_runtime import FloorRuntime
        runtime_class = FloorRuntime
    runtime = runtime_class.from_handoff(dict(path=str(handoff), sha256=handoff_sha256), output/'authentication')
    require(runtime.total_native_calls == 0, 'Authentication must not solve the model')
    identity = runtime.identity()
    if plan is not None:
        require(identity == plan.get('identity', plan), 'Run numerical identity differs')
        if 'handoff' in plan:
            require(plan['handoff']['sha256'] == handoff_sha256, 'Run handoff identity differs')
    with runtime.native_bindings(), forbid_solves(runtime):
        # Preserve the actual dated parameters, including shocked psi, in each pickle.
        for item in record['diagnostic_packets']:
            with gzip.open(pinned(item['path'], item['sha256']), 'rb') as stream:
                packet = pickle.load(stream)
            require(packet['period'] == item['period'] and packet['dated_rent'] == item['dated_rent'],
                    'Actual native packet date/rent differs')
            require(all(k in packet for k in ('parameters', 'b_grid', 'evaluation', 'shared')),
                    'Actual native diagnostic objects missing')
            del packet
        names = runtime.ctx['manifest']['standard_diagnostic_names']
        rendered = runtime.scaffold.render_diagnostics(record['diagnostic_packets'], output/'diagnostics',
                                                       runtime.rt['audit'], names)
    require(runtime.total_native_calls == 0, 'Diagnostic reporting must make zero native model calls')
    retained = verify_outputs(output/'diagnostics', rendered, names, record['diagnostic_packets'], record['rows'])
    receipt = dict(schema='current_floor_saved_diagnostic_report_v1', diagnostic_only=True,
                   fit_certified=False, production_certified=False, visual_review_pending=True,
                   native_model_calls=0, native_gates=record['gates'],
                   accounting_valid=record.get('accounting_valid'), identity=identity,
                   reporter_sha256=sha(__file__),
                   source_inventory_canonical_sha256=hashlib.sha256(json.dumps(identity['source_pins'], sort_keys=True,
                       separators=(',', ':')).encode()).hexdigest(),
                   handoff=dict(path=str(handoff), sha256=handoff_sha256),
                   mapping_record=dict(path=str(mapping_record), sha256=mapping_record_sha256),
                   diagnostic_packets=record['diagnostic_packets'], retained_dates=retained,
                   market_residual=record['market_residual'], fiscal_residual=record['fiscal_residual'])
    provenance = runtime.handoff['checkpoint_and_sources']
    receipt['source_pin_manifest_sha256'] = (provenance['source_pin_manifest_local']['sha256']
        if 'source_pin_manifest_local' in provenance else provenance['source_pins_manifest_sha256'])
    if run_identity is not None:
        receipt['run_identity'] = dict(path=str(run_identity), sha256=run_identity_sha256)
    (output/'report.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    return receipt

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('handoff', 'mapping-record'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--run-identity', type=Path)
    parser.add_argument('--run-identity-sha256')
    result = render(**vars(parser.parse_args()))
    print(json.dumps(dict(status='diagnostic_rendered', native_model_calls=0,
                          retained_dates=len(result['retained_dates']), fit_certified=False)))

if __name__ == '__main__':
    main()
