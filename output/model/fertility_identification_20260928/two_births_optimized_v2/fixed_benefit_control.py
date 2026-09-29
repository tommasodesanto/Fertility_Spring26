"""One fixed-benefit mechanism control; not a replacement-stationary calibration.

Reuse the verified v2 household model and all native non-renewal gates. Report
rather than enforce the renewal residual because this control fixes reference
child benefit. Target2.1 remains visible. No target or accepted model is changed.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[key]='1'
import argparse,copy,csv,gzip,pickle,signal,subprocess,sys,time,types
from pathlib import Path
import driver

HERE=Path(__file__).resolve().parent
OUT=HERE/'fixed_benefit_v1'


def worker(out):
    out.mkdir(parents=True,exist_ok=False)
    smoke=driver.read(HERE/'smoke_v1/worker/receipt.json')
    verified=driver.read(HERE/'run_v1/worker/receipt.json')
    assert smoke['status']==verified['status']=='passed'
    assert smoke['source_pins']==verified['source_pins']
    for name,digest in smoke['source_pins'].items(): assert driver.sha(HERE/name)==digest,name
    c,obj,evaluator,reference,point,ref_receipt=driver.load_runtime(out)
    driver.install_overlay(out,evaluator)
    from intergen_eqscale_seq_optimized.adult_entry import require_closed_stationary_renewal
    reference_psi=float(reference['parameters'].psi_child)
    old_bind=evaluator.bind
    def bind(self,trial):
        P=old_bind(trial);P.two_births_per_period=True;P.psi_child=reference_psi
        return P
    evaluator.bind=types.MethodType(bind,evaluator)
    calls=[]
    def fixed(self,trial,folder,deadline):
        assert not calls, 'Only one equilibrium authorized'
        calls.append(reference_psi)
        P=self.bind(trial)
        payroll,fiscal_rule=self.adapter.pension_tax_from_demographics(P)
        ledger=[dict(status='started',psi_child=reference_psi,scope='fixed-benefit control',epoch=time.time())]
        driver.write(folder/'stationary_solves.json',ledger)
        started=time.monotonic()
        sol,P,price,fiscal=self.rt['solve_balanced_initial_equilibrium'](
            model=self.rt['model'],parameters=P,b_grid=self.selected['b_grid'],
            initial_prices=reference['solution'].p_eq,payroll_tax=payroll,
            marginal_tolerance=1e-9,fiscal_tolerance=1e-6,warm_price_state={})
        self.adapter.verify_pension_ratio(fiscal)
        assert float(P.psi_child)==reference_psi, 'Fixed-benefit control changed child benefit'
        elapsed=time.monotonic()-started
        completed=float(self.rt['chain'].extract_moments(sol,P)['tfr'])
        E=float(sol.entry_rate);births=float(sol.adult_entry_adjusted_birth_children)
        assert E>0 and abs(births/E-completed)<5e-9
        try:
            renewal=require_closed_stationary_renewal(E,births,5e-4)
            renewal.update(status='passes_at_fixed_benefit',enforced_for_this_diagnostic=False)
        except ValueError as error:
            if not str(error).startswith('Closed birth-entry renewal failed:'): raise
            renewal=dict(status='fails_replacement_reported_control_only',enforced_for_this_diagnostic=False,
                entry_E=E,potential_B=births/2.1,births_per_entry=births/E,
                fertility_gap=births/E-2.1,entry_residual=E-births/2.1,native_gate_error=str(error))
        norm=dict(status='fixed_reference_benefit_diagnostic_not_normalized',initial_psi=reference_psi,
            psi_child=reference_psi,completed_fertility=completed,target=2.1,
            absolute_gap=abs(completed-2.1),stationary_solves=1,stationary_solve_seconds=elapsed)
        ledger[0].update(status='completed',seconds=elapsed,completed_fertility=completed,
            price=float(price[0]),market_residual=float(sol.timings['best_eq_error']))
        driver.write(folder/'stationary_solves.json',ledger)
        return (sol,P,price,elapsed,norm),fiscal_rule,renewal
    evaluator.normalize=types.MethodType(fixed,evaluator)
    evaluator.c=copy.deepcopy(evaluator.c)
    evaluator.c['normalization'].update(initial_psi=reference_psi,maximum_stationary_solves=1,
        diagnostic_fixed_benefit=True,demographic_renewal_enforced=False)
    evaluator.c['economic_changes']=[
        'Experimental two-birth household choice tree identical to verified v2 normalized point.',
        'Mechanism control: hold reference child benefit fixed; report fertility-target and renewal misses rather than adjusting benefit.',
        'All ten reference calibration coordinates and remaining primitives fixed; no target, weight, bound, household or market tolerance change.'
    ]
    case=out/'case'
    receipt=evaluator.evaluate(point,case,deadline_epoch=time.time()+840,graphs=True)
    assert len(calls)==1
    with gzip.open(case/'initial_state.pkl.gz','rb') as stream: packet=pickle.load(stream)
    cache=driver.audit_extra_cache(packet,evaluator.rt['model'])
    driver.write(case/'extra_cache_audit.json',cache)
    receipt.update(status='audited_fixed_benefit_control_not_demographic_steady_state',
        reference_label=driver.LABEL,reference_checkpoint_sha256=ref_receipt['case_checkpoint_sha256'],
        effective_source_manifest_sha256=driver.sha(out/'effective_source_manifest.json'),
        extra_cache_audit=cache,free_count=0,estimated_coordinates_this_experiment=0,
        normalized_count=0,normalized_coordinates=0,held_reference_coordinate_count=10,
        held_reference_child_benefit=True,scientific_promotion=False,
        control_source_sha256=driver.sha(__file__),source_pins=smoke['source_pins'])
    driver.write(case/'receipt.json',receipt)
    parameters=list(csv.DictReader((case/'parameters.csv').open()))
    for row in parameters:
        if row['parameter'] in point: row['status']='held at reference estimate during fixed-benefit control'
        if row['parameter']=='psi_child': row['status']='held at reference child benefit; diagnostic only'
        if row['parameter']=='child_benefit_CRRA_coefficient': row['status']='derived from held reference child benefit'
    with (case/'parameters.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(parameters[0]),lineterminator='\n');writer.writeheader();writer.writerows(parameters)
    driver.write(case/'artifact_hashes.json',{str(path.relative_to(case)):driver.sha(path)
        for path in sorted(case.rglob('*')) if path.is_file() and path.name!='artifact_hashes.json'})
    driver.write(out/'receipt.json',dict(status='passed_fixed_benefit_diagnostic_checks',case=str(case),
        calibration_candidate=False,demographic_renewal_enforced=False,model_solves=1,
        target_fertility=2.1,actual_fertility=receipt['normalization']['completed_fertility'],
        source_pins=smoke['source_pins'],source_sha256=driver.sha(__file__),scientific_promotion=False))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--worker',action='store_true');args=parser.parse_args()
    if args.worker: worker(OUT/'worker');return
    OUT.mkdir(exist_ok=False)
    for name in ('latest_completed.json','best_so_far.json'):
        driver.write(OUT/name,dict(status='none_completed',scope='one fixed-benefit diagnostic'))
    started=time.time()
    with (OUT/'worker.log').open('w') as stream:
        process=subprocess.Popen([sys.executable,__file__,'--worker'],stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        while process.poll() is None:
            elapsed=time.time()-started
            driver.write(OUT/'heartbeat.json',dict(status='running',elapsed_seconds=elapsed,budget_seconds=900,pid=process.pid))
            if elapsed>=900:
                os.killpg(process.pid,signal.SIGTERM)
                try: process.wait(timeout=10)
                except subprocess.TimeoutExpired: os.killpg(process.pid,signal.SIGKILL);process.wait()
                driver.write(OUT/'completion.json',dict(status='timeout',elapsed_seconds=time.time()-started));raise SystemExit(2)
            time.sleep(15)
    if process.returncode==0:
        for name in ('latest_completed.json','best_so_far.json'): driver.write(OUT/name,driver.read(OUT/'worker/receipt.json'))
    driver.write(OUT/'completion.json',dict(status='completed' if process.returncode==0 else 'failed',exit_code=process.returncode,elapsed_seconds=time.time()-started))
    raise SystemExit(process.returncode)

if __name__=='__main__':main()
