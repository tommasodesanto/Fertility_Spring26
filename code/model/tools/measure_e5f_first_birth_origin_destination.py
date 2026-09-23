#!/usr/bin/env python3
"""Read the active stationary first-birth branch from a pinned checkpoint.

This is a read-only behavioral diagnostic.  D1 is the active scorer's branch;
D0 is an origin-date realized-choice diagnostic, not an annual event study.
"""
from __future__ import annotations
import argparse, dataclasses, hashlib, json, math, time
from pathlib import Path
from typing import Any
import numpy as np
import audit_e5f_earnings_entry_checkpoint as checkpoint_loader
from audit_e5f_earnings_entry_checkpoint import _load_checkpoint, _get

def validate_source(root, manifest, expected_hash):
    if sha256(manifest) != expected_hash:
        raise ValueError('source manifest SHA-256 mismatch')
    rows = [r for r in json.loads(manifest.read_text())['files'] if r['label'].startswith('source:')]
    if not rows:
        raise ValueError('empty frozen source manifest')
    for row in rows:
        relative = Path(row['label'].split(':', 1)[1])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('invalid relative source path')
        if sha256(root / relative) != row['sha256']:
            raise ValueError(f'frozen source mismatch: {relative}')
    return {'source_files_verified': len(rows), 'manifest_sha256': expected_hash}

def verified_gap(observed, expected, tolerance):
    if not all(math.isfinite(x) for x in (observed, expected, tolerance)) or tolerance <= 0:
        raise ValueError('finite D1 values and positive finite tolerance required')
    gap = abs(observed - expected)
    if gap > tolerance:
        raise RuntimeError(f'saved D1 mismatch: observed={observed:.16g}, expected={expected:.16g}, gap={gap:.3e}')
    return gap

def sha256(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def require(v: Any, name: str) -> Any:
    if v is None: raise ValueError(f'checkpoint lacks {name}')
    return v

def assert_profile(P: Any) -> None:
    if (int(getattr(P,'housing_event_horizon',-1))!=1 or not bool(getattr(P,'sequential_births',False))
        or any(bool(getattr(P,n,False)) for n in ('joint_nested_choice','fertility_nest_choice','two_shock_choice'))
        or str(getattr(P,'child_state_mode',''))!='independent_count'):
        raise ValueError('requires active one-period sequential independent-child observer without a joint nest')
    if float(getattr(P, 'period_years', -1)) != 4 or float(getattr(P, 'da', -1)) != 4:
        raise ValueError('requires frozen four-year age and decision periods')

def origin_branches(e, P, native, smoke: bool):
    """The exact non-joint birth/control loop in active ``begin``."""
    g=np.asarray(require(_get(e,'g_pre'),'evaluation.g_pre'),dtype=float); pol=require(_get(e,'policy'),'evaluation.policy')
    if g.ndim!=7 or g.shape[3]!=int(P.J): raise ValueError('expected seven-axis pre-fertility g_pre')
    tr,co=np.zeros_like(g),np.zeros_like(g); fec=native.get_fecundity_by_age(P); settled=native.readiness_settled_state(P)
    found=False; positive=0
    for j in range(int(P.J)-1):
        if not int(P.A_f_start)<=j+1<=int(P.A_f_end): continue
        for z in range(g.shape[4]):
            x=float(fec[j])*g[:,:,:,j,z,0,settled]*pol.fert_probs[:,:,:,j,z,1]; mass=float(x.sum())
            if mass>1e-15: positive+=1
            tr[:,:,:,j,z,1,1]=x; co[:,:,:,j,z,0,settled]=x; found=found or mass>1e-15
    if not np.isfinite(tr).all() or np.any(tr < 0): raise ValueError('invalid first-birth mass')
    if smoke and not found: raise RuntimeError('no positive age-income first-birth block')
    counts={'positive_age_income_blocks':positive,'selected_blocks':1 if smoke else positive}
    if smoke:
        # Select by birth mass alone, before viewing housing; tiny isolated
        # cohorts would needlessly stress absolute pruning in this loop smoke.
        masses=tr.sum(axis=(0,1,2,5,6))
        j,z=np.unravel_index(np.argmax(masses),masses.shape)
        tblock=tr[:,:,:,j,z].copy(); cblock=co[:,:,:,j,z].copy()
        tr.fill(0); co.fill(0); tr[:,:,:,j,z]=tblock; co[:,:,:,j,z]=cblock
        counts.update(smoke_selection='largest origin birth mass before housing measurement',smoke_age_index=int(j),smoke_income_index=int(z))
    return tr,co,settled,counts

def origin_means(tr,co,e,P,measurement,native):
    pol=e.policy; mass=float(tr.sum())
    if mass<=1e-14 or not math.isclose(mass,float(co.sum()),rel_tol=0,abs_tol=2e-13): raise RuntimeError('origin treated/control mass gate failed')
    def current(g,label):
        x=native.realize_current_cross_section(g,pol.loc_probs,pol.tenure_choice,pol.tenure_probs,pol.maps.lmm_idx,pol.maps.lmm_wt,pol.maps.tmx_idx,pol.maps.tmx_wt,use_compiled_scatter=bool(getattr(P,'use_numba_scatter',False)))
        return measurement.normalize_branch_transport_mass(x,expected_mass=mass,stage=f'{label}_origin_current_choice')
    tc,tg=current(tr,'treated'); cc,cg=current(co,'control'); tm,cm=float(tc.sum()),float(cc.sum())
    if not math.isclose(tm,cm,rel_tol=0,abs_tol=2e-11): raise RuntimeError('origin realized masses differ')
    th=float(np.sum(measurement.calendar.housing_demand_by_location(tc,pol.hR_pol,P))/tm); ch=float(np.sum(measurement.calendar.housing_demand_by_location(cc,pol.hR_pol,P))/cm)
    if not all(math.isfinite(x) for x in (th,ch)): raise RuntimeError('nonfinite origin housing')
    return {'origin_mass':mass,'treated_origin_mean_housing':th,'control_origin_mean_housing':ch,'d0_origin_postbirth_minus_control':th-ch,'origin_mass_gates':[tg,cg]}

def observe(packet: Any, saved_d1: float, tolerance: float, smoke: bool, source_root: Path):
    from intergen_eqscale_seq_optimized import solver as native
    import run_e5f_transition_calibration as measurement
    for module in (native, measurement, measurement.calendar, measurement.transition):
        if not Path(module.__file__).resolve().is_relative_to(source_root.resolve()):
            raise ValueError(f'module loaded outside frozen source: {module.__file__}')
    P=require(_get(packet,'parameters'),'parameters'); e=require(_get(packet,'evaluation'),'evaluation'); bg=np.asarray(require(_get(packet,'b_grid'),'b_grid'),dtype=float); shared=require(_get(packet,'shared'),'shared')
    assert_profile(P); measurement.calendar.model=native
    tr,co,settled,counts=origin_branches(e,P,native,smoke); origin=origin_means(tr,co,e,P,measurement,native)
    e1=e
    if smoke:
        try:
            gmask=np.zeros_like(e.g_pre); chosen=tr[...,1,1]>0
            gmask[...,0,settled]=np.where(chosen,e.g_pre[...,0,settled],0.0)
            e1=dataclasses.replace(e,g_pre=gmask)
        except TypeError as exc: raise TypeError('evaluation must be a dataclass for --smoke-one-cohort') from exc
    print('origin cohorts and choice-mass gates checked; advancing saved-policy branches', flush=True)
    branch=measurement.begin_dated_first_birth_housing_branch(e1,P,bg,shared,origin_period=0)
    if not math.isclose(float(branch['origin_mass']),origin['origin_mass'],rel_tol=0,abs_tol=2e-11): raise RuntimeError('helper origin mass differs from active begin branch')
    dst=measurement.finish_dated_first_birth_housing_branch(branch,e1,P,bg,shared,destination_period=1)
    for k in ('housing_response','treated_mean_housing','control_mean_housing','treated_continuation_births','destination_mass'):
        if not math.isfinite(float(dst[k])): raise RuntimeError(f'nonfinite active D1 field {k}')
    if smoke:
        return {'status':'smoke_partial','origin_means_finite':True,'active_branch_origin_mass':float(branch['origin_mass']),'destination_mass':float(dst['destination_mass']),'branch_counts':counts,'d1_aggregate_gate':'NOT_RUN','d0_full_result':'NOT_EXPOSED'}
    gap=verified_gap(float(dst['housing_response']), saved_d1, tolerance)
    return {'status':'complete','origin':origin,'destination':dst,'branch_counts':counts,'d1_destination_birth_minus_control':float(dst['housing_response']),'saved_d1':saved_d1,'d1_absolute_gap':gap,'d1_gate':'PASS','d0_origin_postbirth_minus_control':origin['d0_origin_postbirth_minus_control']}

def main():
    a=argparse.ArgumentParser(description=__doc__); a.add_argument('--checkpoint',type=Path,required=True); a.add_argument('--checkpoint-sha256',required=True); a.add_argument('--source-root',type=Path,required=True); a.add_argument('--saved-d1',type=float,required=True); a.add_argument('--d1-tolerance',type=float,default=2e-10); a.add_argument('--smoke-one-cohort',action='store_true'); a.add_argument('--output',type=Path,required=True)
    a.add_argument('--source-manifest', type=Path, required=True)
    a.add_argument('--source-manifest-sha256', required=True)
    a.add_argument('--loader-sha256', required=True)
    x=a.parse_args(); start=time.monotonic()
    verified_gap(x.saved_d1, x.saved_d1, x.d1_tolerance)
    loader_path=Path(checkpoint_loader.__file__).resolve()
    if sha256(loader_path) != x.loader_sha256:
        raise ValueError('external checkpoint loader SHA-256 mismatch')
    source_receipt = validate_source(x.source_root, x.source_manifest, x.source_manifest_sha256)
    if sha256(x.checkpoint)!=x.checkpoint_sha256: raise ValueError('checkpoint SHA-256 mismatch')
    if x.output.exists(): raise FileExistsError(x.output)
    print('frozen source and checkpoint hashes verified; loading saved policies', flush=True)
    r=observe(_load_checkpoint(x.checkpoint,x.source_root),x.saved_d1,x.d1_tolerance,x.smoke_one_cohort,x.source_root)
    r.update(schema='e5f_first_birth_origin_destination_v2',checkpoint=str(x.checkpoint),checkpoint_sha256=x.checkpoint_sha256,source_root=str(x.source_root),source_validation=source_receipt,helper_sha256=sha256(Path(__file__)),elapsed_seconds=time.monotonic()-start,diagnostic_only=True,limitation='Active stationary behavioral branch diagnostic only; no annual interpolation, regression, target replacement, household solve, or GE solve.')
    r['external_checkpoint_loader']={'path':str(loader_path),'sha256':x.loader_sha256}
    x.output.write_text(json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+'\n')
if __name__=='__main__': main()
