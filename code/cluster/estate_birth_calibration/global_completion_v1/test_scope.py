"""Small zero-solve checks for the derivative's numerical classification and point routing."""
import ast
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OUT=ROOT/'output/model/experiments/birth_count_choice/estate_a_global_completion_20261005_v1/deployment'

def main():
    source=(OUT/'explore.py').read_text()
    tree=ast.parse(source)
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='classify_native_runtime')
    scope={};exec(compile(ast.Module(body=[node],type_ignores=[]),str(OUT/'explore.py'),'exec'),scope)
    classify=scope['classify_native_runtime']
    price=classify(RuntimeError('native GE acceptance failed: uncomputed_price_unbracketed'))
    assert price['status']=='inadmissible_numerical' and price['rejection_kind']=='price_unbracketed_diagnostic_caps' and 'loss' not in price
    budget=classify(RuntimeError('native GE acceptance failed: uncomputed_bounded_budget'))
    assert budget['status']=='budget_exhausted' and 'loss' not in budget
    class InheritedDistributionInfeasible(RuntimeError):
        classification='inherited_distribution_infeasible'
        audit={'inherited_infeasible_mass':1.77536618261e-12}
    inherited=classify(InheritedDistributionInfeasible('native inherited gate'))
    assert inherited['status']=='inadmissible_numerical' and inherited['rejection_kind']=='typed_inherited_distribution_gate' and 'loss' not in inherited
    assert classify(RuntimeError('native GE acceptance failed: uncomputed_price_other')) is None
    class FakeInheritedDistributionInfeasible(RuntimeError): pass
    assert classify(FakeInheritedDistributionInfeasible('native inherited gate')) is None
    new=json.loads((OUT/'control/plan.json').read_text())
    old=json.loads((OUT/'control/original_plan.json').read_text())
    assert new['points']==old['points'] and len(new['points'])==64
    selected=sum(new['chunks'],[])
    assert len(selected)==36 and len(set(selected))==36
    assert sorted(selected+new['previous_completed_indices']+new['previous_fatal_attempted_indices'])==list(range(64))
    print('zero_solve_scope_tests_passed: 36 unique original points; recognized errors rejected without loss; unknown errors fatal')

if __name__=='__main__':main()
