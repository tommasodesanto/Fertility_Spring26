"""Illustrative cohort completion with 2023 age-specific birth rates frozen.

This empirical accounting scenario is separate from both the economic model and
the CPS children-ever-born survey series. It does not estimate future rates.
"""
from pathlib import Path
import csv
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence'
SOURCE=BASE/'source/stock_forecast/nchs_asfr_1990_2023.csv'
AGE_COLUMNS=['age10_14','age15_19','age20_24','age25_29','age30_34','age35_39','age40_44','age45plus']


def main():
    raw=list(csv.DictReader(SOURCE.open()))
    rates={int(r['year']):np.array([float(r[k])/1000 for k in AGE_COLUMNS]) for r in raw}
    assert sorted(rates)==list(range(1990,2024))
    tfr_checks=[dict(year=int(r['year']),reconstructed=5*float(rates[int(r['year'])].sum()),
                     published=float(r['published_tfr'])) for r in raw]
    assert max(abs(r['reconstructed']-r['published']) for r in tfr_checks)<=.006
    fixed=rates[2023]
    terminal=5*float(fixed.sum())
    assert abs(terminal-1.621)<1e-12
    records=[];cells=[]
    for cohort in range(1980,2016):
        observed=future=0.
        for age in range(10,50):
            year=cohort+age
            source_year=min(year,2023)
            rate=float(rates[source_year][(age-10)//5])
            is_observed=year<=2023
            if is_observed:observed+=rate
            else:future+=rate
            cells.append(dict(birth_cohort=cohort,age=age,calendar_year=year,source_year=source_year,
                              rate=rate,historical_rate=is_observed))
        total=observed+future
        # Independently integrate one age group at a time.
        direct=sum(sum(float(rates[min(cohort+a,2023)][g]) for a in range(10+5*g,15+5*g)) for g in range(8))
        assert abs(total-direct)<1e-12
        if cohort>=2013:assert abs(total-terminal)<1e-12
        records.append(dict(birth_cohort=cohort,age_in_2023=2023-cohort,
                            year_reaching_50=cohort+50,historical_rate_contribution=observed,
                            projected_remaining_contribution=future,projected_completed_fertility=total))
    # A rate schedule constant over the whole life reproduces its period TFR.
    assert abs(sum(float(fixed[(age-10)//5]) for age in range(10,50))-terminal)<1e-12
    fig,ax=plt.subplots(figsize=(10.5,6.1))
    plt.rcParams.update({'font.size':11})
    years=[r['year_reaching_50'] for r in records]
    values=[r['projected_completed_fertility'] for r in records]
    line,=ax.plot(years,values,color='#245f99',lw=2.4,label='Projected completed cohort fertility')
    assert list(line.get_ydata())==values
    ax.axhline(2.1,color='#a77a49',ls='--',lw=1)
    ax.text(2050,2.115,'Approximate replacement: 2.1',color='#906437',fontsize=10)
    ax.axhline(terminal,color='#bd433c',ls='--',lw=1.3,label='Period fertility fixed at 2023: 1.62')
    selected=[r for r in records if r['birth_cohort'] in (1980,1990,2000,2010)]
    for r in selected:
        x,y=r['year_reaching_50'],r['projected_completed_fertility']
        ax.plot(x,y,'o',color='#245f99',ms=5)
        ax.annotate(f'{y:.2f}',(x,y),xytext=(0,11),textcoords='offset points',ha='center',color='#245f99',fontsize=11)
    ax.set(xlim=(2029,2066),ylim=(1.50,2.30),xlabel='Year the cohort reaches age 50',ylabel='Births per woman',
           xticks=[2030,2035,2040,2045,2050,2055,2060,2065])
    top=ax.secondary_xaxis('top',functions=(lambda x:x-50,lambda x:x+50))
    top.set_xlabel('Mother’s birth cohort',labelpad=9)
    top.set_xticks([1980,1985,1990,1995,2000,2005,2010,2015])
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.16)
    ax.legend(loc='upper right',frameon=False,fontsize=9.5)
    fig.suptitle('Completed fertility if 2023 birth rates persisted',x=.09,y=.98,ha='left',fontsize=17)
    fig.text(.09,.04,
             'Illustrative calculation from NCHS annual age-specific birth rates, 1990–2023 (Driscoll and Hamilton, 2025, Tables 2 and 4).\n'
             'Historical rates through 2023; each age-specific rate held at its 2023 level thereafter. Five-year age groups; approximate completion at 50.\n'
             'This projects cohort fertility from birth rates; it is a different measure from the CPS ages 40–44 survey series in the presentation.',
             fontsize=8.4,color='#555555',linespacing=1.5)
    fig.subplots_adjust(left=.09,right=.97,bottom=.24,top=.78)
    out=BASE/'figures'
    for ext in ('png','pdf'):
        fig.savefig(out/f'completed_fertility_constant_rates.{ext}',dpi=170,facecolor='white')
    plt.close(fig)
    for name,rows in [('completed_fertility_constant_rates',records),('completed_fertility_constant_rates_cells',cells)]:
        with (out/f'{name}.csv').open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    receipt=dict(status='PASS',source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                 source_method='Transcribed from accessible NCBI rendering of NCHS Tables2and4; rates per1000 converted to perwoman. Direct HTTP download was unavailable.',
                 source_urls=['https://www.ncbi.nlm.nih.gov/books/NBK617829/table/nvsr74-3.t2/',
                              'https://www.ncbi.nlm.nih.gov/books/NBK617829/table/nvsr74-3.t4/'],
                 cutoff=2023,fixed_schedule_tfr=terminal,cohort_definition='Nominal birth cohort c; age a in calendar year c+a. Uniform rate within each five-year age group; single-age/Lexis reconstruction is approximate.',
                 completion_definition='Sum ages10–49, allocating the published45+ category to five ages45–49, as in its TFR convention. Births outside these bins are not separately modeled.',
                 scope='Mechanical constant-rate scenario; not an official forecast, not the economic model, not a continuation of CPS survey means. Does not model migration selection or birth-history-specific catch-up.',
                 historical_tfr_checks=tfr_checks,independent_grouped_cohort_sums=True,
                 fixed_schedule_completion_identity=True,artist_values_exact=True,selected=selected)
    (out/'completed_fertility_constant_rates_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(selected,indent=2))


if __name__=='__main__':
    main()
