"""Rebuild the complete overnight readout from collected receipts; no model solves."""
from pathlib import Path
from datetime import datetime,timezone
import csv,hashlib,json,math
from xml.sax.saxutils import escape
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912'
def read(p):return json.loads(p.read_text())
def csvrows(p):return list(csv.DictReader(p.open()))
def fmt(v):return f'{float(v):.6g}' if v not in ('',None) else '-'
def main():
    out=ROOT/'output/pdf/e5f_overnight_review.pdf';out.parent.mkdir(exist_ok=True)
    support=BASE/'report_support';support.mkdir(exist_ok=True)
    fit=csvrows(BASE/'verified_final/selected_target_fit.csv');params=csvrows(BASE/'verified_final/selected_parameters.csv')
    summary=read(BASE/'verified_final/summary.json');fiscal=read(BASE/'verified_final/initial_raw_summary.json')['fiscal']['actual_accounts']
    loss=sum(float(r['loss_contribution']) for r in fit if r['loss_contribution'])
    assert len(fit)==13 and len(params)==17 and math.isclose(loss,summary['best_loss'],rel_tol=1e-12)
    assert summary['selected_exact_repetitions_verified']
    assert float(next(r['estimate'] for r in params if r['parameter']=='pension_period'))==fiscal['pension_period_units']
    assert abs(fiscal['scaled_pension_budget_residual'])<1e-6
    short=read(BASE/'matched_short/summary.json');receipt=read(BASE/'matched_short/root_receipt.json')
    assert receipt['finite_horizon_market_fiscal_converged'];m=max(map(abs,receipt['final']['market_residual']));b=max(map(abs,receipt['final']['fiscal_residual']))
    assert m<=2e-4 and b<=1e-6
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    for folder,label,color in [('matched_short','First window matched; decline 0.01414','#17698c'),('native_smoke','Larger decline 0.045','#a4a8ad')]:
        rows=csvrows(BASE/folder/'expected_transition.csv');fert=read(BASE/folder/'fertility.json');years=np.array([int(r['calendar_year']) for r in rows])
        axes[0,0].plot(years+4,[r['period_tfr_topcode_adjusted'] for r in fert],'o--',color=color,label=label)
        axes[0,1].plot(years,[float(r['asset_price']) for r in rows],'o-',color=color,label=label)
        axes[1,0].plot(years,[100*float(r['owner_rate']) for r in rows],'o-',color=color)
    axes[0,0].plot([2011,2015,2019,2023],[1.974875,1.861,1.755375,1.64575],'ks-',label='Data: four-year fertility average')
    axes[0,0].set(title='Fertility: expected paths, not fitted history',ylabel='Period fertility',xlabel='End of birth window');axes[0,0].legend(fontsize=7,loc='lower left')
    axes[0,1].set(title='Expected house prices',ylabel='Model units',xlabel='Year')
    axes[1,0].set(title='Expected ownership',ylabel='Percent of household heads',xlabel='Year')
    axes[1,1].semilogy(years,np.maximum(np.abs(receipt['final']['fiscal_residual']),1e-16),'o-',color='#17698c')
    axes[1,1].axhline(1e-6,color='black',ls=':',label='PAYGO tolerance');axes[1,1].set(title='Pension balance: matched short forecast',ylabel='Absolute relative residual',xlabel='Year');axes[1,1].legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.2);ax.set_xticks(years+4 if ax is axes[0,0] else years);ax.tick_params(axis='x',labelsize=8)
    fig.tight_layout();chart=support/'forecast_comparison.png';fig.savefig(chart,dpi=180);plt.close(fig)
    styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='BodySmall',fontName='Helvetica',fontSize=9,leading=12,spaceAfter=8,textColor=colors.HexColor('#243343')))
    styles.add(ParagraphStyle(name='CellSmall',fontName='Helvetica',fontSize=7.5,leading=9.5));styles.add(ParagraphStyle(name='DeckTitle',fontName='Helvetica-Bold',fontSize=24,leading=28,spaceAfter=15,textColor=colors.HexColor('#19384c')))
    P=lambda text,style='BodySmall':Paragraph(text,styles[style])
    def table(headers,rows,widths):
        data=[[P(escape(str(x)),'CellSmall') for x in row] for row in [headers,*rows]]
        t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e4edf2')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#809aaa')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8f9')])]))
        return t
    when=datetime.now(timezone.utc).strftime('%d September %Y, %H:%M UTC')
    story=[P('Quantitative model<br/>Overnight review','DeckTitle'),P(when),P('The first fertility window now matches. The complete historical transition and new policies are not yet verified.','Heading2')]
    story+=[P(f'<b>Initial calibration:</b> loss {loss:.6f}, down 0.438% from 159.238986. Three bounded batches completed 252 new trials; the final winner reproduced twice. Annual beta remains estimated at its 0.99 cap. This is a modest refinement, not a resolved fit.'),
        P(f'<b>First historical window:</b> the six-date forecast with a permanent preference decline of 0.01414 produces {short["model"]:.6f} against {short["data"]:.6f}. The absolute error is {abs(short["gap"]):.8f}. Housing and pension residuals pass their unchanged gates. This is a short-horizon result, not a four-window historical fit.'),
        P('<b>Long forecast:</b> the earlier, larger decline of 0.045 also cleared the 28-date housing/PAYGO root. It failed terminal-distance checks: terminal household mass remained 12.56% from the stationary endpoint and the pension benefit 7.78% away. A finite equilibrium is not sufficient for horizon certification.'),
        P('<b>Numerical diagnosis:</b> changing only the terminal root starting guesses allowed smaller preference declines to solve. Earlier failures at the initial guess did not prove equilibrium nonexistence. Separate age-advancement mass failures remain unresolved; their tolerance was not relaxed.'),
        P('<b>Still missing:</b> an accepted sequence of four unexpected shocks through 2023, horizon verification of that sequence, and the new property-tax/rebate policy comparisons. No new policy effect should enter slides as an established result.'),
        P('<b>Current follow-up:</b> independent 28- and 56-date checks of the preference that matches the first short window. They use verified numerical starts and unchanged source/target/fiscal contracts. The initial calibration candidate is not silently swapped beneath a historical run.'),
        P('Reading order: page 2 shows every target; page 3 every parameter/restriction; page 4 the actual forecast shapes; page 5 the interpretation limits. The appendix preserves the complete standard initial-calibration graph set.'),PageBreak()]
    story += [P('Complete initial target fit','Heading1'),P('All 12 scored moments plus the separate 2.1 normalization. Weight means the actual objective weight; the normalization has no loss contribution. Labels and values are read directly from the collected table.')]
    story += [table(['Moment','Target','Model','Gap','Weight','Loss'],[[r['label'],*[fmt(r[k]) for k in ('target','model','gap','actual_weight','loss_contribution')]] for r in fit],[170,61,61,61,61,61]),Spacer(1,10),P(f'Sum of the 12 scored contributions: <b>{loss:.9f}</b>. Complete empirical definitions, uncertainty and provenance are retained in selected_target_fit.csv. The first-birth rooms target remains 0.7202462623815278 with its original weight. No target was dropped or reweighted.'),PageBreak()]
    story += [P('Parameters and restrictions','Heading1'),P('Nine structural coordinates remain free. Beta is estimated over [0.94, 0.99], not fixed. The table preserves the search bounds and near-bound flags; derived or external restrictions have no search bounds.')]
    story += [table(['Parameter','Estimate','Lower','Upper','Near bound','Status'],[[r['parameter'],fmt(r['estimate']),fmt(r['lower']),fmt(r['upper']),r['near_bound'],r['status']] for r in params],[129,67,47,47,60,125]),Spacer(1,10),P(f'The stationary pension benefit {fiscal["pension_period_units"]:.6f} is budget-derived. Its actual payroll revenue is {fiscal["payroll_tax_revenue"]:.9f}, pension outlays {fiscal["pension_outlays"]:.9f}, and relative gap {fiscal["scaled_pension_budget_residual"]:.3g}. The historical bridge changes age weights; its balanced pension need not equal this stationary benefit.'),PageBreak()]
    story += [P('What the short forecasts actually show','Heading1'),P('Every model line is the forecast made after one permanent preference change. Future surprises are absent. Matching the first point does not match the later data. Neither six-date line passes the terminal-horizon checks.'),Image(str(chart),width=510,height=357),Spacer(1,10),P(f'Matched-window native root: maximum housing residual {m:.6g} (tolerance 0.0002); maximum PAYGO residual {b:.6g} (tolerance 0.000001). The full root receipt and dated source rows are saved with the figure.'),P('The grey curve uses the larger decline that was extended to 28 dates. The blue curve uses the smaller decline that matches the first window and is now being checked at longer horizons. These are different preference experiments under the same pinned historical initial condition.'),PageBreak()]
    frontier=[]
    for name,delta in [('arm_0',-.0225),('arm_1',-.0275)]:
        d=read(BASE/'frontier'/name/'summary.json');frontier.append([str(delta),f'{d["model"]:.6f}',f'{d["gap"]:.6f}','Pass','Not certified'])
    frontier.insert(0,['-0.01414',f'{short["model"]:.6f}',f'{short["gap"]:.6f}','Pass','Not certified'])
    story += [P('Interpretation and remaining checks','Heading1'),table(['Preference change','First-window fertility','Gap to data','Short root','Horizon'],frontier,[100,110,90,65,110]),Spacer(1,14),P('<b>Information:</b> at each shock date households learn the current preference and expect it to persist. They do not anticipate later preference surprises. A full history must carry the actual household distribution and birth queues from each implemented first period.'),P('<b>Measurement:</b> the fitted flow uses the retained household-rate analogue of female TFR. The 2.1 stationary normalization is a different object. Published TFR targets are four-year averages for windows ending 2011, 2015, 2019 and 2023. This measurement approximation remains explicit.'),P('<b>Population:</b> historical head-age totals are imposed through the existing bridge and the 2023 person anchor. They are not independent model predictions. The outside-origin entry share 0.169 remains a diagnostic/outstanding closure object, not an estimated production normalization.'),P('<b>Fiscal closure:</b> payroll revenue and pension outlays must balance at every accepted date. Property-tax rebates require a separate budget equation. A verified stationary or short-path fiscal balance is not evidence that all policy paths are solved.'),P('<b>Calibration provenance:</b> the final overnight initial refinement and the historical workers use distinct, explicitly pinned candidates. Their losses and parameter tables must not be silently combined. The historical workers retain the previously repeated capped candidate; the new initial refinement is a separate result.'),P('<b>What deserves immediate attention:</b> whether the first-window match survives a longer horizon; then fitting the remaining unexpected shocks; then running the tax/rebate comparisons on that accepted inherited 2023 state. The unresolved mass-gate failures and measurement/entry conventions remain visible.'),PageBreak()]
    images=sorted((BASE/'verified_final/selected_standard_diagnostics').glob('*.png'));assert len(images)==17
    for i in range(0,len(images),2):
        story += [P('Standard initial-calibration diagnostics','Heading1'),P('Complete unchanged graph set. These are stationary candidate diagnostics, not evidence of a fitted 2007-2023 path.')]
        for path in images[i:i+2]:
            img=Image(str(path));img._restrictSize(510,285)
            story += [P(escape(path.stem),'Heading3'),img,Spacer(1,9)]
        if i+2<len(images):story.append(PageBreak())
    def footer(canvas,doc):
        canvas.saveState();canvas.setFont('Helvetica',8);canvas.setFillColor(colors.HexColor('#60717b'));canvas.drawString(35,22,'Diagnostic research output | September 12, 2026');canvas.drawRightString(A4[0]-35,22,str(doc.page));canvas.restoreState()
    SimpleDocTemplate(str(out),pagesize=A4,rightMargin=35,leftMargin=35,topMargin=35,bottomMargin=38,title='Quantitative model - overnight review',author='Research review').build(story,onFirstPage=footer,onLaterPages=footer)
    sources=[BASE/'verified_final/selected_target_fit.csv',BASE/'verified_final/selected_parameters.csv',BASE/'verified_final/summary.json',BASE/'verified_final/initial_raw_summary.json',BASE/'matched_short/root_receipt.json',BASE/'matched_short/summary.json',*images]
    (support/'verification.json').write_text(json.dumps(dict(created=when,target_rows=13,parameter_rows=17,standard_graphs=17,recomputed_loss=loss,matched_short_housing=m,matched_short_paygo=b,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},pdf_sha256=hashlib.sha256(out.read_bytes()).hexdigest()),indent=2)+'\n')
    print(out)
if __name__=='__main__':main()
