"""Zero-solve recovery budget boundary controls and driver diff scope."""
import ast,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
archive=ROOT/'output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/deployment/stage.tar.gz'
with tarfile.open(archive) as tar:
    source=tar.extractfile('source/code/model/experiments/birth_count_choice/cluster_calibrate.py').read().decode()
tree=ast.parse(source)
nodes=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in ('BudgetStop','ensure_native_launch_budget','authenticated_budget_error')]
module=ast.Module(body=nodes,type_ignores=[]);space={};exec(compile(module,'<budget>', 'exec'),space)
BudgetStop=space['BudgetStop'];guard=space['ensure_native_launch_budget'];classify=space['authenticated_budget_error']
guard(100,901)
try:guard(101,901)
except BudgetStop as e:assert str(e)=='native_700s_minimum_plus_100s_guard'
else:raise AssertionError('Prelaunch guard failed')
known=('native GE acceptance failed: uncomputed_bounded_budget','No time reserve for selected reporting and exact repeat','No time for exact repeat')
assert all(classify(m,700,True) for m in known)
assert all(not classify(m,701,True) for m in known)
assert all(not classify(m,699,False) for m in known)
assert not classify('native GE acceptance failed: inadmissible_numerical',0,True)
assert not classify('unexpected solver exception',0,True)
assert source.count('ensure_native_launch_budget(time.time(),deadline-RESERVE)')==1
assert source.count('if authenticated_budget_error(message,remaining,(out/label).is_dir()):')==1
assert source.count("verify=evaluate('selected_postcheck',chosen['parameters'],deadline)")==1
assert source.count("raise BudgetStop('native_evaluation_budget_exhausted') from exc")==1
print('budget_boundary_tests_passed_zero_solves')
