"""Bounded component checks; no lifecycle or equilibrium solve."""
from pathlib import Path
from types import SimpleNamespace
import copy
import hashlib
import json
import numpy as np
import pytest
from experiments.stationary_single_market import inputs
from experiments.stationary_single_market.contract import validate_contract
from experiments.stationary_single_market.engine import parameters, household, shared, kernels
from refactor_lab.engine import parameters as old_parameters
from refactor_lab.engine import household as old_household
from refactor_lab.engine import shared as old_shared
from refactor_lab.engine import kernels as old_kernels

ROOT = Path(__file__).resolve().parents[5]
PACKAGE = Path(__file__).resolve().parents[1]
BUNDLE_SHA = '427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7'

@pytest.fixture(scope='module')
def reference():
    return inputs.load_inputs(ROOT/'output/model/publication_refactor_20260929/local_export_v1/inputs', ROOT, BUNDLE_SHA)

def test_input_identity_and_contract(reference):
    validate_contract(reference.parameters)
    assert reference.parameters.I == 1
    assert reference.parameters.n_child_states == 4
    assert reference.parameters.Pi_z.shape == (15, 15)

def test_original_sources_preserved():
    receipt = json.loads((PACKAGE/'source_provenance.json').read_text())
    for path, sha in receipt['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha

@pytest.mark.parametrize('n',range(4))
def test_all_valid_child_states(n,reference):
    p=reference.parameters
    for m in range(n+1):
        assert parameters.children_at_home_count(n,m,p) == old_parameters.children_at_home_count(n,m,p)
        assert shared.get_completed_fertility(n,m,p) == old_shared.get_completed_fertility(n,m,p)
        assert household.birth_destination_child_state(p,m) == old_household.birth_destination_child_state(p,m)
        for cutoff in (2,3):
            assert household.current_child_bin_dt(n,m,2,cutoff,'independent_count') == old_household.current_child_bin_dt(n,m,2,cutoff,'independent_count')

@pytest.mark.parametrize('field,value',[('I',2),('child_state_mode','shared_clock'),('sequential_births',False),('joint_nested_choice',True),('permanent_income_levels_enabled',True),('use_loc_kernel',False)])
def test_unsupported_contract_rejected(field,value,reference):
    p=copy.deepcopy(reference.parameters); setattr(p,field,value)
    with pytest.raises(ValueError): validate_contract(p)

def test_income_nodes_not_pinned(reference):
    p=copy.deepcopy(reference.parameters); p.z_grid=np.linspace(.5,1.5,9); p.Pi_z=np.eye(9); p.Nz=9; p.z_weights=np.ones(9)/9
    validate_contract(p)
    bad=copy.deepcopy(p); bad.Nz=15
    with pytest.raises(ValueError): validate_contract(bad)
    bad=copy.deepcopy(p); bad.z_weights=np.ones(15)/15
    with pytest.raises(ValueError): validate_contract(bad)
    p.Pi_z[0,0]=.5
    with pytest.raises(ValueError): validate_contract(p)

def test_precompute_exact(reference):
    # Pure shared construction on authenticated full inputs; does not solve Bellman.
    p=reference.parameters; grid=reference.b_grid
    a=vars(old_shared.precompute_shared(copy.deepcopy(p),grid))
    b=vars(shared.precompute_shared(copy.deepcopy(p),grid))
    assert a.keys()==b.keys()
    for key in a:
        np.testing.assert_array_equal(a[key],b[key],err_msg=key)

@pytest.mark.parametrize('compiled',[False,True])
def test_single_destination_kernel_exact(compiled):
    rng=np.random.default_rng(218); v=rng.normal(size=(7,3,1,4,4))
    v[0]= -1e10; v[1,0,0,0,0] = -1e9; v[2,0,0,0,0] = -1e9+1
    idx=np.zeros((7,1,3),dtype=np.int64); wt=np.zeros_like(idx,dtype=float)
    for shift in (0.,.123):
        args=(v,idx,wt,np.array([[shift]]),.217)
        new=kernels.location_logit_kernel if compiled else kernels.location_logit_kernel.py_func
        old=old_kernels.location_logit_kernel if compiled else old_kernels.location_logit_kernel.py_func
        a=old(*args); b=new(*args)
        for x,y in zip(a,b): np.testing.assert_array_equal(x,y)

@pytest.mark.parametrize('mode,d',[('reference',None),('corrected',0.),('corrected',.14)])
def test_credit_binding_unchanged(mode,d,reference):
    from experiments.stationary_single_market import credit
    from refactor_lab import credit as old_credit
    a=copy.deepcopy(reference.parameters); b=copy.deepcopy(reference.parameters)
    credit.bind_engine_credit(a,mode,d); old_credit.bind_engine_credit(b,mode,d)
    assert inputs.serialized(vars(a)) == inputs.serialized(vars(b))
