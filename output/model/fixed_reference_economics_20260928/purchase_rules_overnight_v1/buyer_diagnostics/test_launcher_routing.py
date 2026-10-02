"""Zero-solve guard for reviewed source and immutable selection-store routing."""
from pathlib import Path

here = Path(__file__).resolve().parent
launch = (here / 'launch_readout.sh').read_text()
preflight = (here / 'preflight_torch.sh').read_text()
stage = (here / 'stage_torch.sh').read_text()
source = '/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a'
store = '/scratch/td2248/projects/purchase_mechanism_v1'
buyer = '/scratch/td2248/projects/purchase_buyer_diagnostics_v3'
assert f'mechanism={source}' in launch and f'mechanism={source}' in preflight
assert f'selection_store={store}' in launch
assert f'remote={buyer}' in launch and f'remote={buyer}' in preflight and f'remote={buyer}' in stage
assert '"$mechanism/source/$packet:$repo/$packet:ro"' in launch
assert '"$selection_store/selected_postchecks:$repo/$packet/results:ro"' in launch
assert '"$selection_store/selection:$repo/$packet/collection/readout:ro"' in launch
assert '"$mechanism/results:$repo/$packet/mechanism/results:ro"' in launch
assert '"$selection_store/selection/manifest.json"' in launch
assert '"$selection_store/selection/selected_${arm}.json"' in launch
assert '"$selection_store" "$arm"' in launch
assert '/purchase_rules_overnight_v1/results:' not in launch
print('reviewed source and immutable selection routing passed; zero solves')
