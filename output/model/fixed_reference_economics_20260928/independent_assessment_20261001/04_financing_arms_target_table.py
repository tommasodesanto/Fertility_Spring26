"""Fourteen target rows across the saved purchase-financing arms at the winner31 reference (fixed price q0 and GE roots).
Read-only: prints existing target_fit.csv files. Run from the repository root."""
import pandas as pd, glob, os
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
L='output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/purchase_ltv_v1/local_run/'
arms={'base 80/80 (q0)':'retry5/results/baseline_80_80','buy90/stay80 (q0)':'retry7/results/purchase_90_stayer_80','both 90/90 (q0)':'retry7/results/both_90_90','buy100/stay80 (q0)':'retry8/results/purchase_100_stayer_80','both 100/100 (q0)':'retry9/results/both_100_100'}
out={}
for k,p in arms.items():
    t=pd.read_csv(L+p+'/target_fit.csv'); out[k]=t.set_index('moment')['model']; tg=t.set_index('moment')['target']
df=pd.DataFrame(out); df.insert(0,'target',tg); print(df.round(5).to_string())
print('\nGE selected repeats (price and population re-solved; completed fertility 2.1 imposed by the renewal closure):')
for g in sorted(glob.glob(L+'ge_retry*/*/selected_repeat_*/target_fit.csv')):
    t=pd.read_csv(g).set_index('moment')['model']
    print(g.replace(L,''),{m:round(float(t[m]),5) for m in ('mean_rooms','first_birth_rooms','recent_parent_ownership','ownership_30_55','early_fertility','cps_childlessness','nchs_mean_age')})
