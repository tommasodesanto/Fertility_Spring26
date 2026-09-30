"""Deterministically prepare extra nearby starts; never run a native model."""
import copy
import datetime
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT = HERE.with_name('entry_calibration_round1_v1')
PARAMETERS = ('beta_annual','chi','first_birth_fixed_cost','kappa_fert',
              'kappa_fert_continuation','theta0','delta_alpha_jump',
              'child_benefit_curvature','tenure_choice_kappa')
DIRECTIONS = {
    's1': [1,-1,-1,1,-1,1,-1,1,-1],
    's2': [-1,1,1,-1,1,-1,1,-1,1],
    's3': [-1,-1,1,1,1,-1,1,-1,1],
}
COMMON_END = '2026-10-01T02:38:17Z'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare():
    parent_plan = PARENT/'plan.json'
    plan = json.loads(parent_plan.read_text())
    base_lanes = copy.deepcopy(plan['lanes'])
    bounds = {r['parameter']: (float(r['lower']),float(r['upper']))
              for r in plan['reference_parameter_table'] if r['parameter'] in PARAMETERS}
    dispatch = []
    for base, config in base_lanes.items():
        seed = config['seed']
        h = {k: max(1e-6,min(.02*max(abs(seed[k]),.01),.005*(bounds[k][1]-bounds[k][0])))
             for k in PARAMETERS}
        for suffix, direction in DIRECTIONS.items():
            lane = base+'_'+suffix
            candidate = copy.deepcopy(config)
            raw = {k: 3*h[k]*direction[i] for i,k in enumerate(PARAMETERS)}
            point = {k: min(bounds[k][1],max(bounds[k][0],seed[k]+raw[k])) for k in PARAMETERS}
            candidate['seed'] = point
            candidate['parent_lane'] = base
            candidate['parent_seed_parameter_table'] = candidate.pop('seed_parameter_table')
            candidate['seed_provenance']['refers_to'] = 'parent pilot selected point, not this perturbed starting point'
            candidate['starting_point_construction'] = {
                'method': 'Parent selected seed plus three finite-difference steps times a fixed direction, clipped to existing bounds',
                'parameter_order': list(PARAMETERS), 'direction': direction,
                'step_sizes': h, 'requested_offsets': raw,
                'effective_offsets': {k: point[k]-seed[k] for k in PARAMETERS},
                'parent_seed': seed, 'parent_plan_sha256': sha(parent_plan),
                'native_verified': False, 'starting_point_full_ge_repeats': 0,
                'required_gate': 'Fresh native baseline plus independent baseline repeat before any search; no fallback on failure',
            }
            plan['lanes'][lane] = candidate
            plan['arms'][lane] = plan['arms'][base]
            dispatch.append(lane)
    plan['experiment'] = 'entry_calibration_multistart_v1'
    plan['dispatch_lanes'] = dispatch
    plan['existing_lanes_not_dispatched'] = list(base_lanes)
    plan['common_deadline_iso'] = COMMON_END
    plan['common_deadline_epoch'] = datetime.datetime.fromisoformat(COMMON_END.replace('Z','+00:00')).timestamp()
    plan['deadline_contract'] = 'Same absolute hard end as the already-running three lanes; launcher must pass this epoch, not a fresh four-hour clock'
    plan['parent_plan_sha256'] = sha(parent_plan)
    (HERE/'plan.json').write_text(json.dumps(plan,sort_keys=True,indent=2)+'\n')
    old = json.loads((PARENT/'source_pins.json').read_text())
    pins = {}
    for rel,digest in old.items():
        replaced = rel.replace(str(PARENT.relative_to(ROOT)),str(HERE.relative_to(ROOT)))
        pins[replaced] = sha(ROOT/replaced) if replaced!=rel else digest
        # The independent identity test reads the original five compact files.
        # Pin and stage them too; the frozen cluster checkout predates round 1.
        if replaced!=rel:
            pins[rel] = digest
    for name in ('prepare_multistart.py','test_multistart.py'):
        path = HERE/name
        pins[str(path.relative_to(ROOT))] = sha(path)
    (HERE/'source_pins.json').write_text(json.dumps(pins,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'lanes_in_plan':len(plan['lanes']),'dispatch_lanes':dispatch,
                      'source_pins':len(pins),'common_deadline':COMMON_END}))


if __name__ == '__main__':
    prepare()
