"""Stdlib-only immutable pair preparation; zero model imports/solves."""
import ast, hashlib, json, shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OLD=ROOT/'output/model/fixed_reference_economics_20260928/credit_no_taper_v1/small_credit_v1'
SCALAR=ROOT/'output/model/publication_refactor_20260929/scalar_src_gridfix/refactor_lab/engine/kernels.py'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def segments(p):
    text=p.read_text(); lines=text.splitlines(keepends=True)
    return {n.name: ''.join(lines[min([n.lineno]+[d.lineno for d in n.decorator_list])-1:n.end_lineno]) for n in ast.parse(text).body if isinstance(n,ast.FunctionDef)}
indexed=OLD/'source/small_credit_lab/engine/kernels.py'
assert sha(indexed)=='fe7d43afdac234af71edd09e1260666e21c309dda83de4adc6a3abff35c3d5d4'
assert sha(SCALAR)=='19dceb70a6b4684ad8da04e119198bef0d88d51d5ee0489debcb297d4434a99e'
x,y=segments(indexed),segments(SCALAR)
assert [k for k in y if x.get(k)!=y[k]]==['exhaustive_saving_scalar']
assert set(x)-set(y)=={'_rank','_interp_ranked','_renter_value','_owner_value'}
assert not (set(y)-set(x))
for arm in ['scalar','indexed']:
    dest=HERE/'arms'/arm
    assert not dest.exists(), 'Refusing overwrite'
    shutil.copytree(OLD/'source',dest/'source',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for name in ['driver.py','single_price.py','phase_a.py','phase_b_ge.py']:
        shutil.copy2(OLD/name,dest/name)
    if arm=='scalar': shutil.copy2(SCALAR,dest/'source/small_credit_lab/engine/kernels.py')
    (dest/'run_arm.py').write_text('''import sys\nimport driver\nOriginalBudget=driver.Budget\nclass PairBudget(OriginalBudget):\n    def __init__(self,*args,**kwargs):\n        super().__init__(*args,**kwargs)\n        if self.smoke: (self.out / "phase_b_ge").mkdir(parents=True, exist_ok=True)\n        else: self.max_lifecycle=6\ndriver.Budget=PairBudget\ndriver.main()\n''')
# Every shared file must be byte-identical; only saving kernels differ.
a,b=HERE/'arms/scalar',HERE/'arms/indexed'
paths=[p.relative_to(a) for p in a.rglob('*') if p.is_file()]
assert set(paths)=={p.relative_to(b) for p in b.rglob('*') if p.is_file()}
differences=[str(p) for p in paths if sha(a/p)!=sha(b/p)]
assert differences==['source/small_credit_lab/engine/kernels.py']
receipt=dict(status='prepared_zero_solves',differences=differences,scalar_kernel_sha256=sha(SCALAR),indexed_kernel_sha256=sha(indexed),bundle_sha256='427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7',new_budget_seconds=2400,max_lifecycle_per_arm=6,source_experiment='18869900',comparison_scope='saving optimization only; same corrected D=.14 full renewal/population GE')
(HERE/'preparation.json').write_text(json.dumps(receipt,indent=2)+'\n')
with (HERE/'source.sha256').open('w') as out:
    for p in sorted(HERE.rglob('*')):
        if p.is_file() and p.name!='source.sha256' and '__pycache__' not in p.parts:
            out.write(f'{sha(p)}  {p.relative_to(HERE)}\n')
print(json.dumps(receipt))
