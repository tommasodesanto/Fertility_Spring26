"""Read actual frozen seed and utility bindings; no household/equilibrium solve."""
from pathlib import Path
import copy, gzip, hashlib, json, pickle, sys, types
import numpy as np

ROOT = Path(__file__).resolve().parents[8]
HERE = Path(__file__).resolve().parent
BUNDLE = ROOT / 'tmp/utility_overnight_20260923_v1'
# Compatibility only for the serialized checkpoint's Python/plot imports.
import pathlib
pathlib.__path__ = []
local = types.ModuleType('pathlib._local')
for name in ('Path', 'PosixPath', 'PurePath'):
    setattr(local, name, getattr(pathlib, name))
sys.modules['pathlib._local'] = local
mat = types.ModuleType('matplotlib'); mat.__path__ = []; mat.use = lambda *a, **k: None
plt = types.ModuleType('matplotlib.pyplot'); mat.pyplot = plt
sys.modules['matplotlib'] = mat; sys.modules['matplotlib.pyplot'] = plt
sys.path[:0] = [str(BUNDLE / p) for p in (
    'source/code/model/tools', 'source/code/model',
    'source/code/model/intergen_eqscale_seq_optimized', 'tools')]
from intergen_eqscale_seq_optimized import solver as model
from intergen_eqscale_seq_optimized.parameters import build_debt_caps
import e5f_parenthood_utility as parent
from e5f_stationary_paygo import stationary_age_income_mass, bind_initial_balanced_pension
import build_period_earnings_process as income
import run_e5f_preference_share_candidate as shares

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

checkpoint = BUNDLE / 'inputs/normalized_checkpoint/normalized_old.pkl.gz'
with gzip.open(checkpoint, 'rb') as f:
    old = pickle.load(f)['old']
arms = {}
for arm in ('B_floor', 'B_shares'):
    plan_path = BUNDLE / f'inputs/utility_templates/{arm}/template_plan.json'
    plan = json.loads(plan_path.read_text())
    if arm == 'B_shares':
        shares.install_share_utility(parent)
    P = parent.initialize_parenthood_utility(copy.deepcopy(old.parameters))
    P = parent.bind_parenthood_utility(P, plan['structural_parameters'])
    overrides, _ = income.build_period_earnings_process(
        **plan['income_specification']['constructor_arguments'])
    for key, value in overrides.items():
        setattr(P, key, copy.deepcopy(value))
    P = build_debt_caps(P)
    # Same analytic fiscal binding as native run; no Bellman or GE call.
    P, _ = bind_initial_balanced_pension(P, payroll_tax=.179)
    mass = stationary_age_income_mass(P).sum(axis=(0, 2))
    z = np.asarray(P.z_grid); weights = np.asarray(P.z_weights)
    gross = np.array([[model.annual_gross_income_at_state(P, 0, j, float(v))
                       for v in z] for j in range(int(P.J_R))])
    annual_mean = float(mass[:P.J_R] @ (gross @ weights) / mass[:P.J_R].sum())
    sample_b = np.array([-2., -1., -.5, 0., 1.])
    floor = model.renter_borrowing_floor(P, sample_b, 0)
    assert np.array_equal(floor, np.minimum(sample_b, 0.))
    assert P.transfer_floor_G0 == P.transfer_floor_Gn == P.lambda_d == 0.
    assert P.c_bar_0 == P.c_bar_n == 0. and P.preference_spec == 'eqscale'
    assert int(P.J_R) == 12 and P.period_years == 4 and P.age_start == 18
    row = dict(template_sha256=sha(plan_path), R_gross=float(P.R_gross),
        annual_return=float(P.R_gross ** (1 / P.period_years) - 1),
        lambda_d=float(P.lambda_d), debt_taper_start_age=float(P.debt_taper_start_age),
        debt_taper_end_age=float(P.debt_taper_end_age),
        next_period_taper_s=float(P.debt_taper_weights[1]),
        next_period_debt_cap_D=float(P.debt_caps[1]),
        transfer_floor_G0=float(P.transfer_floor_G0), transfer_floor_Gn=float(P.transfer_floor_Gn),
        childless_cbar=0., childless_hbar=0., use_age_survival=bool(P.use_age_survival),
        survival_probs=np.asarray(P.survival_probs).tolist(),
        income_age_profile=np.asarray(P.income_age_profile).tolist(),
        income_age_breaks=np.asarray(P.income_age_breaks).tolist(),
        income_age_values=np.asarray(P.income_age_values).tolist(),
        working_age_mass=mass[:P.J_R].tolist(),
        working_age_nodes=(P.age_start + np.arange(P.J_R) * P.period_years).tolist(),
        working_aftertax_period_income=np.asarray(P.income)[0, :P.J_R].tolist(),
        annual_gross_income_by_age_state=gross.tolist(),
        mean_annual_gross_working_income=annual_mean,
        z_grid=z.tolist(), z_weights=weights.tolist(),
        period_years=float(P.period_years), tau_pay=float(P.tau_pay))
    arms[arm] = row
for key in arms['B_floor']:
    if key != 'template_sha256':
        assert arms['B_floor'][key] == arms['B_shares'][key], key
paths = [checkpoint] + [BUNDLE / p for p in (
    'source/code/model/intergen_eqscale_seq_optimized/solver.py',
    'source/code/model/intergen_eqscale_seq/parameters.py',
    'source/code/model/tools/e5f_parenthood_utility.py',
    'source/code/model/tools/e5f_stationary_paygo.py',
    'source/code/model/tools/e5f_social_security.py',
    'tools/run_e5f_earnings_wealth_candidate.py',
    'tools/run_e5f_preference_share_candidate.py')]
receipt = dict(schema='actual_frozen_parameter_extraction_v1',
    status='passed_no_solve_inspection', household_solves=0, equilibrium_solves=0,
    method='Actual serialized seed; frozen utility/income/debt and analytic pension binding. Supply rebase omitted because it changes only H0 and elasticity, neither used in this mapping.',
    source_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths}, arms=arms)
(HERE / 'actual_frozen_parameters.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({k: arms['B_floor'][k] for k in (
    'R_gross', 'mean_annual_gross_working_income', 'use_age_survival',
    'next_period_taper_s', 'next_period_debt_cap_D')}, indent=2))
