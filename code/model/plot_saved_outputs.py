"""Regenerate standard PDFs from the saved 2007 and 2023 snapshots; no solves."""
from pathlib import Path
import gzip, pickle, shutil, sys, json
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'code/model'), str(ROOT/'code/model/tools')]
from production.storage import load_case
from model_data_assessment import prepare_assessment
OUT = ROOT/'output/model/production'
result, case = load_case(OUT/'2007/solution')
prepare_assessment(result, case, output=OUT/'2007', save_pages=False)
snapshot_dir=(OUT/'2023/solution.pkl.gz').resolve().parent
if not json.loads((snapshot_dir/'state_receipt.json').read_text()).get('actual_dated_state_exact'):
    raise RuntimeError('2023 dated snapshot needs preparation with tools/plot_transition_model_data.py')
with gzip.open(OUT/'2023/solution.pkl.gz', 'rb') as stream:
    result = pickle.load(stream)
case = (OUT/'2023/solution.pkl.gz').resolve().parent
prepare_assessment(result, case, output=OUT/'2023', data_year=2023,
                   psid_path=case/'data/psid_recent.csv',
                   acs_path=case/'data/housing_profile_by_age.csv', save_pages=False)
shutil.copyfile(OUT/'transition/solution/slide_plot/fertility_2007_2063.pdf',
                OUT/'transition/fertility_path.pdf')
print(OUT)
