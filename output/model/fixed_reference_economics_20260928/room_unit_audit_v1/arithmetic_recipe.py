"""Standalone arithmetic on saved text inputs: no model imports or solves."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
OUT = ROOT / 'output/model/fixed_reference_economics_20260928/room_unit_audit_v1'
OUT.mkdir(exist_ok=True)
REF = ROOT / 'output/model/fertility_identification_20260928'
LIVE = ROOT / 'output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT'
P = json.loads((REF / 'fixed_reference_manifest.json').read_text())['actual_serialized_parameters']
PS = json.loads((REF / 'resume_v1/selected_export/primary/standard_diagnostics/summary.json').read_text())
F = {x['parameter']: float(x['estimate']) for x in csv.DictReader((LIVE / 'parameters.csv').open())}
FC = json.loads((LIVE / 'closure.json').read_text())
FS = json.loads((LIVE / 'standard_diagnostics/summary.json').read_text())
PINS = json.loads((ROOT / 'output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/source_pins.json').read_text())
pin_checks = {}
for f in ['household.py', 'shared.py', 'kernels.py', 'child_preferences.py', 'equilibrium.py', 'distribution.py']:
    p = Path('code/model/refactor_lab/engine') / f
    digest = hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
    pin_checks[str(p)] = {'sha256': digest, 'matches_floor_source_pin': digest == PINS[str(p)]}

checks = []
def check(name, original, converted, expected_factor=1., **detail):
    expected = original * expected_factor
    error = abs(converted - expected) / max(abs(expected), 1e-300)
    checks.append(dict(name=name, original=original, converted=converted,
                       expected_factor=expected_factor, relative_error=error, **detail))

def compensation(a, a0, rent):
    logk = a*math.log(a)+(1-a)*math.log((1-a)/rent)
    logk0 = a0*math.log(a0)+(1-a0)*math.log((1-a0)/rent)
    return math.exp(logk0-logk)

def material(c, H, m, owner, a0, jump, rent, compensated, floor, chi, sigma):
    a = a0-jump if m > 0 else a0
    physical_surplus = H-(floor if m > 0 else 0.)
    if physical_surplus <= 0:
        return None
    services = physical_surplus*(chi if owner else 1.)
    A = compensation(a,a0,rent) if compensated else 1.
    e = ((2+.7*m)/2)**.7 if m > 0 else 1.
    return (A*c**a*services**(1-a)/e)**(1-sigma)/(1-sigma)

def probabilities(values, scale):
    z=max(values)
    e=[math.exp((x-z)/scale) for x in values]
    s=sum(e)
    return [x/s for x in e]

for label, par, price, Q, supply, pop in [
    ('adopted_block0506', P, PS['owner_asset_price'][0], PS['aggregate_housing_demand'], PS['aggregate_housing_supply'], 1.),
    ('experimental_floor_chain7_0173', F, FC['price'], FC['normalized_housing_demand'], FC['physical_housing_supply'], FC['population_scale']),
]:
    a0=par['alpha_cons']; sigma=par['sigma']; xi=.63
    H0=P['H0'][0]; rbar=P['r_bar'][0]; rent=P['user_cost_rate']*price
    jump=P['delta_alpha_jump'] if label.startswith('adopted') else 0.
    compensated=label.startswith('adopted')
    floor=0. if compensated else F['h_P']
    rstar=P['utility_reference_rent']; chi=par['chi']
    psi=par['psi_child']; gamma=par['child_benefit_curvature']
    for lam in [.1, 10.]:
        K=lam**((1-a0)*(1-sigma))
        for H in P['H_own']:
            check('purchase_value',price*H,(price/lam)*(lam*H),specification=label,lambda_h=lam,H=H)
            check('down_payment',.2*price*H,.2*(price/lam)*(lam*H),specification=label,lambda_h=lam,H=H)
            check('rent_payment',rent*H,(rent/lam)*(lam*H),specification=label,lambda_h=lam,H=H)
        for m in range(4):
            for owner in [False,True]:
                for c in [.1,1.,10.]:
                    for H in [.35,2.,2.3,2.4,4.,6.,10.]:
                        old=material(c,H,m,owner,a0,jump,rstar,compensated,floor,chi,sigma)
                        new=material(c,lam*H,m,owner,a0,jump,rstar/lam,compensated,lam*floor,chi,sigma)
                        if old is None:
                            if new is not None:
                                raise ArithmeticError('Physical feasibility failed to transform')
                            continue
                        check('material_utility',old,new,K,specification=label,lambda_h=lam,m=m,owner=owner,c=c,H=H)
                        benefit=psi*m**(1-gamma) if m else 0.
                        check('material_plus_child_benefit',old+benefit,new+K*benefit,K,specification=label,lambda_h=lam,m=m,owner=owner,c=c,H=H)
        check('supply_at_saved_price',H0*(rent/rbar)**xi,(lam*H0)*((rent/lam)/(rbar/lam))**xi,lam,specification=label,lambda_h=lam)
        check('inverse_supply_price',rbar*(supply/H0)**(1/xi),(rbar/lam)*((lam*supply)/(lam*H0))**(1/xi),1/lam,specification=label,lambda_h=lam)
        check('endogenous_population_ratio',supply/Q,(lam*supply)/(lam*Q),specification=label,lambda_h=lam)
        check('inverse_supply_slope',rbar*(supply/H0)**(1/xi)/(xi*supply),(rbar/lam)*((lam*supply)/(lam*H0))**(1/xi)/(xi*lam*supply),1/lam**2,specification=label,lambda_h=lam)
        check('supply_5percent_expansion_price_factor',1.05**(1/xi),((1.05*lam*supply)/(lam*supply))**(1/xi),specification=label,lambda_h=lam)
        for b in [0.,1.,10.]:
            old=par['theta0']*((par['theta1']+b)**(1-sigma)-par['theta1']**(1-sigma))/(1-sigma)
            new=K*par['theta0']*((par['theta1']+b)**(1-sigma)-par['theta1']**(1-sigma))/(1-sigma)
            check('bequest_utility',old,new,K,specification=label,lambda_h=lam,b=b)
        # Illustrative algebraic value contrasts, not stored Bellman values.
        vals=[-2.,-.7-par['first_birth_fixed_cost'],-1.1]
        for key in ['kappa_fert','kappa_fert_continuation','tenure_choice_kappa']:
            old=probabilities(vals,par[key]); new=probabilities([K*x for x in vals],K*par[key])
            for i in range(len(old)):
                check('illustrative_softmax_'+key,old[i],new[i],specification=label,lambda_h=lam,alternative=i)
        gap=.123; w=128.02070205233477
        check('room_moment_loss_units',w*gap*gap,(w/lam**2)*(lam*gap)**2,specification=label,lambda_h=lam)

results = {
    'scope':'Read-only equation audit and standalone arithmetic on saved text inputs. No model imports, optimization, simulation, policy recomputation, test suite or equilibrium solves.',
    'reference':'2007 stationary reference block0506, September 28 verified export',
    'separate_experiment':'October 1 09:41 NY floor chain7 case0173_nm; not adopted; corrected D=0 nonnegative entry; not D=.25 or D=.53',
    'source_pin_checks':pin_checks,
    'arithmetic_count':len(checks),
    'maximum_relative_error':max(x['relative_error'] for x in checks),
    'all_arithmetic_within_1e_minus_12':all(x['relative_error']<1e-12 for x in checks),
    'factor_for_ten_room_units':.1**((1-P['alpha_cons'])*(1-P['sigma'])),
    'price_percent_for_5percent_quantity_expansion':100*(1.05**(1/.63)-1),
    'limitations':['No rescaled equilibrium or policy-array comparison was run.','Fixed feasibility sentinels, regularizers and solver tolerances require a numerical audit before certifying arbitrary-unit execution.','Arithmetic proves coordinate equivalence, not empirical validation of service mapping, price levels, or elasticities.'],
    'checks':checks,
}
(OUT/'arithmetic.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps({k:v for k,v in results.items() if k!='checks'},indent=2))
