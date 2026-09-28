"""Supplemental lifecycle figure from verified small tables; Torch, zero solves."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['MPLBACKEND'] = 'Agg'
import sys
import csv
import json
import hashlib
from pathlib import Path
assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit()
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

out = Path(__file__).resolve().parent
source = out/'fertility_lifecycle_matched_windows.csv'
rows = list(csv.DictReader(source.open()))
assert [(r['age_lower'], r['age_upper']) for r in rows] == [
    ('20','24'), ('25','29'), ('30','34'), ('35','39'), ('40','44')]
x = np.arange(len(rows))
labels = [r['age_lower']+'–'+r['age_upper'] for r in rows]
plt.rcParams.update({'font.size': 12, 'axes.spines.top': False,
                     'axes.spines.right': False, 'font.family': 'DejaVu Sans'})
fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.5))
series = []
for ax, field, title, ymax in zip(axes, ['capped3','mother_share','given_mother'],
        ['Children per woman', 'Women who are mothers', 'Children among mothers'], [2., 1., 2.5]):
    data = np.array([float(r['data_'+field]) for r in rows])
    model = np.array([float(r['model_'+field]) for r in rows])
    ld, = ax.plot(x, data, 's--', color='#b65e27', linewidth=2.2, markersize=6, label='Data: CPS 2004/2006')
    lm, = ax.plot(x, model, 'o-', color='#126c8c', linewidth=2.4, markersize=6, label='Model: 2007 stationary reference')
    np.testing.assert_array_equal(ld.get_ydata(), data)
    np.testing.assert_array_equal(lm.get_ydata(), model)
    ax.set(title=title, ylim=(0, ymax), xlabel='Age group', xticks=x, xticklabels=labels)
    ax.grid(axis='y', alpha=.2)
    ax.tick_params(axis='both', labelsize=11)
    if field == 'mother_share':
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    else:
        ax.set_ylabel('Children ever born, capped at 3')
    series.append(dict(metric=field, data=data.tolist(), model=model.tolist(),
                       gaps=(model-data).tolist()))
axes[0].annotate('Largest shortfall', xy=(2, 1.0960299155087398), xytext=(1., .20),
    fontsize=10, color='#126c8c', arrowprops=dict(arrowstyle='->',color='#126c8c'))
axes[0].annotate('Catch-up by 40–44', xy=(4,1.7307694210670574), xytext=(2.1,1.91),
    fontsize=10, color='#303030', arrowprops=dict(arrowstyle='->',color='#303030'))
fig.suptitle('Fertility over the lifecycle: model versus data', fontsize=20, x=.075, ha='left', y=.985)
fig.text(.075, .914, '2007 stationary reference — block0506, September 28 verified export', fontsize=11, color='#555555')
handles, legend_labels = axes[0].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc='lower center', ncol=2, bbox_to_anchor=(.5,.093), frameon=False, fontsize=11)
fig.text(.5, .039,
    'Each point averages the same five-year age group in data and model. Lines join those averages.\n'
    'Different women are observed at different ages; this is not a cohort followed over time. Supplemental figure; no new calibration.',
    ha='center', va='center', fontsize=10, color='#555555')
fig.subplots_adjust(left=.075, right=.985, top=.81, bottom=.25, wspace=.31)
figure = out/'fertility_lifecycle_comparison.png'
fig.savefig(figure, dpi=160, facecolor='white')
plt.close(fig)
qa = dict(status='saved_table_plot_verified', model_solves=0, slurm_job=os.environ['SLURM_JOB_ID'],
    input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    image_sha256=hashlib.sha256(figure.read_bytes()).hexdigest(),
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    plotted_series=series, standard_17_plots_changed=False,
    interpretation='Cross-sectional age profiles with identical windows, not a cohort catch-up estimate')
(out/'lifecycle_plot_qa.json').write_text(json.dumps(qa, indent=2)+'\n')
print(json.dumps(qa))
