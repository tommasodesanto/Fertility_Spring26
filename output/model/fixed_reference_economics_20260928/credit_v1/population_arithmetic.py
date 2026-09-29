"""Small fixed-age-rate population illustration; no household or equilibrium solve."""
import csv
import hashlib
import json
import math
from pathlib import Path
import signal
import sys

LABEL = '2007 stationary reference — block0506, September 28 verified export'
ROOT = Path('/scratch/td2248/projects/fixed_reference_credit_20260929')
FROZEN = Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def main():
    assert sys.platform == 'linux'
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('30-second budget')))
    signal.alarm(30)
    target = ROOT / 'results/population_arithmetic_v1.json'
    assert not target.exists()
    manifest_path = FROZEN / 'output/model/fertility_identification_20260928/fixed_reference_manifest.json'
    assert sha(manifest_path) == '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
    manifest = read(manifest_path)
    # Read actual saved parameter values; never source-constructor defaults.
    params = manifest['actual_serialized_parameters']
    survival = params['survival_probs']
    top = params['tfr_top_bin_weight']
    assert params['adult_entry_clock'] == 'split_birth_vintage' and params['period_years'] == 4
    run = ROOT / 'results/solve_v1'
    complete = read(run / 'completed.json')
    assert complete['status'] == 'passed'
    schedules, pins = {}, {}
    for name in ('control', 'grid_control', 'credit'):
        case = run / name
        receipt = read(case / 'receipt.json')
        record = next(r for r in complete['completed'] if r['case'] == name)
        assert sha(case / 'receipt.json') == record['receipt_sha256']
        assert sha(case / 'observers.json') == receipt['artifact_hashes']['observers.json']
        rows = read(case / 'observers.json')['fertility']['uniform_birth_time']['accounting']['parity_birth_flows_by_age']
        lifecycle = list(csv.DictReader((case / 'lifecycle_2023.csv').open()))
        ages = [float(r['pre_fertility_mass']) for r in lifecycle]
        births = [r[0] + r[1] + (top-2)*r[2] for r in rows]
        assert abs(sum(births)/ages[0] - receipt['completed_fertility']) < 1e-10
        schedules[name] = dict(rates=[b/n for b,n in zip(births,ages)], ages=ages,
            births=sum(births), fertility=receipt['completed_fertility'])
        pins[name] = dict(receipt_sha256=record['receipt_sha256'], observers_sha256=sha(case/'observers.json'),
            lifecycle_sha256=sha(case/'lifecycle_2023.csv'))
    initial = schedules['control']['ages']
    half_prehistory = schedules['control']['births']/2.1/2
    results = {}
    for name, spec in schedules.items():
        age = initial.copy()
        q16, q20 = [half_prehistory]*3, [half_prehistory]*4
        population = [sum(age)]
        for t in range(100):
            births = sum(n*b for n,b in zip(age,spec['rates']))
            entry = q16.pop(0)+q20.pop(0)
            q16.append(births/2.1/2)
            q20.append(births/2.1/2)
            age = [entry] + [age[j]*survival[j] for j in range(len(survival))]
            population.append(sum(age))
        life = [1.]
        for s in survival: life.append(life[-1]*s)
        weights = [s*b/2.1 for s,b in zip(life,spec['rates'])]
        def characteristic(lam):
            return sum(w*.5*(lam**(-j-4)+lam**(-j-5)) for j,w in enumerate(weights))
        lo, hi = .9, 1.1
        for _ in range(100):
            mid = (lo+hi)/2
            if characteristic(mid)>1:lo=mid
            else:hi=mid
        lam=(lo+hi)/2
        assert abs(characteristic(lam)-1)<1e-12
        assert abs(population[-1]/population[-2]-lam)<1e-5
        results[name] = dict(annual_growth_pct=100*(lam**.25-1), reproduction_ratio=sum(weights),
            mean_native_generation_years=sum(w*(4*j+18) for j,w in enumerate(weights))/sum(weights),
            household_population={str(y):population[y//4] for y in (0,4,16,20,40,60,80)})
    assert max(abs(n-1) for n in results['control']['household_population'].values())<1e-5
    increases={str(y):100*(results['credit']['household_population'][str(y)]/
        results['grid_control']['household_population'][str(y)]-1) for y in (0,4,16,20,40,60,80)}
    output=dict(reference_label=LABEL,status='passed',source_sha256=sha(__file__),input_pins=pins,
        budget=dict(cases=3,periods_per_case=100,seconds=30,model_solves=0),
        assumption='Hold each case age-specific adjusted birth rates fixed at its recomputed cohort averages. Begin with exact reference age masses and reference birth prehistory. Retain survival and half-16/half-20-year entry with births/2.1 conversion.',
        limitation='Household population units, not resident-person headcounts. Approximate demographic projection: within-age wealth/family composition and prices/pensions do not evolve consistently. Not a native policy transition or equilibrium steady state.',
        results=results,credit_vs_matched_baseline_population_pct=increases)
    target.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(output,indent=2))


if __name__ == '__main__':main()
