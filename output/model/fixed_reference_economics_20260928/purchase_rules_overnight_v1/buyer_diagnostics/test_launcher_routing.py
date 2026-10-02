"""Zero-solve guard for reviewed source and immutable selection-store routing."""
from pathlib import Path

here = Path(__file__).resolve().parent
launch = (here / 'launch_readout.sh').read_text()
preflight = (here / 'preflight_torch.sh').read_text()
stage = (here / 'stage_torch.sh').read_text()
source = '/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a'
store = '/scratch/td2248/projects/purchase_mechanism_v1'
buyer = '/scratch/td2248/projects/purchase_buyer_diagnostics_v6'
extension = '/scratch/td2248/projects/purchase_mechanism_horizon_extension_v1'
assert f'mechanism={source}' in launch and f'mechanism={extension}' in launch
assert f'mechanism={extension}' in preflight
assert f'selection_store={store}' in launch
assert f'remote={buyer}' in launch and f'remote={buyer}' in preflight and f'remote={buyer}' in stage
assert '"$mechanism/source/$packet:$repo/$packet:ro"' in launch
assert '"$selection_store/selected_postchecks:$repo/$packet/results:ro"' in launch
assert '"$selection_store/selection:$repo/$packet/collection/readout:ro"' in launch
assert '"$mechanism/results:$repo/$packet/mechanism/results:ro"' in launch
assert 'verify_dated_extension.py' in launch
assert 'dated_receipt=$("$python"' in launch
assert 'case_06_quarter_control_h48' in launch and 'case_07_quarter_temporary_h48' in launch
assert '"$date" == date_000' in launch
assert '"$selection_store/selection/manifest.json"' in launch
assert '"$selection_store/selection/selected_${arm}.json"' in launch
assert '"$selection_store" "$arm"' in launch
assert '/purchase_rules_overnight_v1/results:' not in launch
readout = (here / 'run_selected.py').read_text()
assert 'root / "standard_diagnostics"' in readout
assert 'repeat / "standard_diagnostics"' not in readout
print('reviewed source and immutable selection routing passed; zero solves')
