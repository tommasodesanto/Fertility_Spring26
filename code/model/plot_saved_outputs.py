"""Regenerate standard PDFs from the saved 2007 and 2023 snapshots; no solves."""
from pathlib import Path
import gzip, pickle, shutil, sys, json, csv
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
# Concise 2023 fit readout from exactly the numbers drawn in the PDF.
shutil.copyfile(OUT/'transition/solution/fertility_fit.csv',OUT/'2023/target_fit.csv')
with (OUT/'2023/plotted_long.csv').open() as stream:
    plotted=list(csv.DictReader(stream))
def points(panel,series):
    return [(float(r['x']),float(r['value'])) for r in plotted if r['panel']==panel and r['series']==series]
def mean_count(panel,series):
    return sum(x*y for x,y in points(panel,series))
def cdf_mean(panel,series):
    last=0.; total=0.
    for x,f in points(panel,series): total+=x*(f-last); last=f
    return total
coverage=json.loads((OUT/'2023/metadata.json').read_text())['coverage']
def mean_children(panel,series):
    return coverage[panel]['mean_model' if series=='Model' else 'mean_data']
checks=[]
for label,panel,operator in [
    ('Children ever born, ages 22–25 (model 3+ adjusted)','children_22_25',mean_children),
    ('Childlessness, ages 22–25','children_22_25',lambda p,s:points(p,s)[0][1]),
    ('Children ever born, ages 40–44 (model 3+ adjusted)','children_40_44',mean_children),
    ('Childlessness, ages 40–44','children_40_44',lambda p,s:points(p,s)[0][1]),
    ('First-birth mean age (band midpoints)','first_birth_age',lambda p,s:sum((20+4*x)*y for x,y in points(p,s))),
    ('First-birth share age 30+','first_birth_age',lambda p,s:sum(y for x,y in points(p,s) if x>=3)),
    ('Mean rooms, all households (PSID)','rooms_all',cdf_mean),
    ('Mean rooms, owners (PSID)','rooms_owner',cdf_mean),
    ('Mean rooms, renters (PSID)','rooms_renter',cdf_mean),
    ('Mean total net wealth / own mean earnings','resource_totalwealth_cdf',cdf_mean),
    ('Negative financial position, all ages','resource_financial_cdf',lambda p,s:max([y for x,y in points(p,s) if x<0],default=0))]:
    checks.append((label,operator(panel,'Data'),operator(panel,'Model')))
text=['# 2023 snapshot: model and data','',
      'These cross-sectional comparisons are validation moments, not fitted targets. CPS is June 2024; first-birth timing is NCHS 2023; resources and room distributions use PSID 2019.','',
      '| Validation moment | Data | Model | Model − data |','|---|---:|---:|---:|']
text += [f'| {label} | {data:.4f} | {model:.4f} | {model-data:+.4f} |' for label,data,model in checks]
text += ['', '## Complete target system for the retained transition', '',
         'Only the 2020–2023 window was fitted. These are the original household-rate fertility statistics, distinct from children-ever-born stocks above.', '',
         '| Birth window | Target | Model | Gap | Weight | Loss contribution |','|---|---:|---:|---:|---:|---:|']
with (OUT/'2023/target_fit.csv').open() as stream:
    for r in csv.DictReader(stream):
        text.append(f"| {r['birth_year_start']}–{r['birth_year_end']} | {float(r['target']):.6f} | {float(r['model']):.6f} | {float(r['gap']):+.6f} | {float(r['weight']):g} | {float(r['loss_contribution']):.8g} |")
text += ['', 'Estimated shock: psi_child=0.1199969464, bounds [0.0017892072, 0.3578414413], away from either bound. [All 31 supplied baseline parameters and reference bounds](../transition/solution/baseline_parameters.csv).', '',
         '## CPS data check','',
         'The weighted five-year age groups reproduce Census Table 1 population totals and children-count shares to its published rounding. Ages 35–45 have about 600–700 respondents per single age. The annual-age line is a cross-section of different cohorts, not a trajectory for the same women. No monotonicity is imposed.', '',
         'Mean graphs use the retained model 3+ weight, 3.602359422009, and CPS public counts (five or more coded five). Distribution bars still group 3+. Applying the completed-fertility tail mean at younger ages and within income groups is a reporting approximation: fourth and later birth dates are not separately modeled. The previous capped-at-three view was comparable on its own terms, but omitted this existing measurement adjustment.', '',
         '[Official Census Table 1](https://www2.census.gov/programs-surveys/demo/tables/fertility/2024/am-women-fertility/t1.xlsx). Census documentation identifies PRTAGE as the masked public-use age variable; masking is a possible contributor to single-age irregularity, not an established explanation of these specific jumps. Mean graph definitions now include the existing top-bin adjustment; model behavior and calibration targets are unchanged.']
(OUT/'2023/model_data_fit.md').write_text('\n'.join(text)+'\n')
print(OUT)
