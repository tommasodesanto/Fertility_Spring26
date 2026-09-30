"""Synthetic shape/feasibility check; zero model calls."""
import json, shutil, tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import run_comparison as runner

def main():
    with tempfile.TemporaryDirectory(dir=runner.HERE) as scratch:
        out=Path(scratch)
        for arm,nb,nz in [('control_160x15',160,15),('proposal_120x9',120,9)]:
            shape=(nb,2,1,2,nz,2,2)
            value=np.ones(shape)
            value[:,:,:,1]+=2 # isolate age, not location
            value[:,:,:,:,:,0,1]=-1e10 # inadmissible m>n
            if nb==120:value[0,0,0,0,0,1,0]=-1e10 # one feasibility mismatch
            tp=np.zeros(shape+(2,));tp[...,0]=.3;tp[...,1]=.7
            sol=SimpleNamespace(V=value,bp_pol=value.copy(),c_pol=value.copy(),hR_pol=value.copy(),tenure_probs=tp)
            live=dict(sol=sol,b_grid=np.arange(nb,dtype=float),P=SimpleNamespace(z_grid=np.linspace(0,1,nz)))
            snapshot=runner.policy_snapshot(live)
            assert snapshot['owner_choice_probability'].shape==shape
            assert np.all(snapshot['owner_choice_probability']==.7)
            folder=out/arm/'phase_b_ge/selected_root';folder.mkdir(parents=True)
            np.savez(folder/'common_support_policies.npz',**snapshot)
            for name in ['target_fit.csv','parameters.csv']:shutil.copy(runner.INPUT_PACKET/('reference_'+name),folder/name)
            (folder/'closure.json').write_text(json.dumps(dict(price=1,population_scale=1)))
        runner.compare(out)
        result=json.loads((out/'comparison_common_support_policies.json').read_text())
        assert result['feasibility']['old_only_state_count']==1
        assert result['V']['max_abs']==0 and result['V']['youngest_max_abs']==0 and result['V']['oldest_max_abs']==0
        assert result['owner_choice_probability']['max_abs']<1e-15
        print(json.dumps(dict(status='synthetic_policy_shape_feasibility_passed',lifecycle_solves=0,tenure_input_dimensions=8,owner_output_dimensions=7)))
if __name__=='__main__':main()
