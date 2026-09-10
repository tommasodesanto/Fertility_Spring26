"""Research decision memo, sourced appendix and frozen numerical receipt tables.

Read-only with respect to model/data/production contracts. No model solve.
Run with the bundled Python runtime. Source notes and CSVs remain alongside it.
"""
from pathlib import Path
import csv
import json
import re
from html import unescape
from xml.sax.saxutils import escape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[3]
OUT = ROOT / 'output/pdf/calibration_design_decision_20260910.pdf'
FONTS = Path('/System/Library/Fonts/Supplemental')
for name, filename in [('Arial','Arial.ttf'),('Arial-Bold','Arial Bold.ttf'),('Arial-Italic','Arial Italic.ttf')]:
    pdfmetrics.registerFont(TTFont(name, str(FONTS/filename)))
pdfmetrics.registerFontFamily('Arial', normal='Arial', bold='Arial-Bold', italic='Arial-Italic', boldItalic='Arial-Bold')
styles = {
    'title': ParagraphStyle('title',fontName='Arial-Bold',fontSize=20,leading=24,spaceAfter=15),
    'body': ParagraphStyle('body',fontName='Arial',fontSize=10.5,leading=14.2,spaceAfter=9),
    'heading': ParagraphStyle('heading',fontName='Arial-Bold',fontSize=12,leading=16,spaceBefore=8,spaceAfter=7),
    'small': ParagraphStyle('small',fontName='Arial',fontSize=9,leading=12,spaceAfter=8),
    'cell': ParagraphStyle('cell',fontName='Arial',fontSize=8.6,leading=11),
    'compact': ParagraphStyle('compact',fontName='Arial',fontSize=7.7,leading=9.7),
}
WIDTH = A4[0]-88
story, md = [], []
def plain(s):
    s=re.sub(r'<b>(.*?)</b>',r'**\1**',s)
    s=re.sub(r'<i>(.*?)</i>',r'*\1*',s)
    s=re.sub(r'<link href="([^"]+)">(.*?)</link>',r'[\2](\1)',s)
    s=s.replace('<br/>',' ')
    return unescape(re.sub(r'<[^>]+>','',s))
def rich(s):
    # Arial lacks Unicode subscript digits; use positioned ordinary glyphs.
    for digit,ch in enumerate('₀₁₂₃₄₅₆₇₈₉'):
        s=s.replace(ch,f'<sub>{digit}</sub>')
    return s
def para(s,style='body',markdown=None):
    story.append(Paragraph(rich(s),styles[style])); md.append(markdown if markdown is not None else plain(s))
def page(title):
    if story: story.append(PageBreak())
    para(title,'title', '# '+title)
def heading(s): para(s,'heading','## '+s)
def table(rows, fractions, compact=False):
    st=styles['compact' if compact else 'cell']
    cells=[[Paragraph(rich(escape(str(v))),st) for v in row] for row in rows]
    t=Table(cells,colWidths=[WIDTH*x/sum(fractions) for x in fractions],repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#ededed')),
                          ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),5),
                          ('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),3 if compact else 4),
                          ('BOTTOMPADDING',(0,0),(-1,-1),3 if compact else 4),
                          ('LINEBELOW',(0,0),(-1,0),.5,colors.HexColor('#777777')),
                          ('LINEBELOW',(0,1),(-1,-1),.2,colors.HexColor('#dddddd'))]))
    story.append(t);story.append(Spacer(1,9))
    md.append('\n'.join(['| '+' | '.join(map(str,rows[0]))+' |','| '+' | '.join(['---']*len(rows[0]))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows[1:]]))
def rows(path): return list(csv.DictReader(path.open()))
housing={(r['window'],r['moment']):r for r in rows(BASE/'housing/early_housing_target_candidates.csv')}
wealth={(r['window'],r['moment']):r for fn in ['aggregate_wealth_results.csv','old_wealth_results.csv'] for r in rows(BASE/'wealth'/fn)}
def h(moment): return housing['2005_2006',moment]
def w(moment): return wealth['preannouncement_2003_2005',moment]
def fm(x): return f'{float(x):.6f}'
def point_se(r,point='estimate',se='bootstrap_se'): return fm(r[point])+' ('+fm(r[se])+')'
psid_ratio=w('aggregate_wealth_gross_labor_earnings')
psid_tail=w('old_p90_p50')
psid_level=w('old_p50')

page('A pre-decline benchmark, then a fertility shock')
para('Decision report for Tommaso De Santo | 10 September 2026','small')
para('<b>I recommend calibrating an approximate pre-2007 economy first, then holding its structural parameters fixed while fitting an announced fertility-preference decline. Policy scenarios should begin from the resulting 2023 distribution.</b> This is a recommendation for a new empirical design; no target, parameter, model specification or cluster job was changed during this review.')
para('The purpose is exactly the exercise you described: assume that child preferences declined, condition on the transition this produces, and assess housing policies. The fitted preference change absorbs what this maintained model needs to reproduce the selected fertility decline. It does not identify why preferences changed, or establish that preferences caused the decline.')
heading('The design in four steps')
para('<b>1. Fit the initial benchmark.</b> Use observations from before the announcement: NCHS 2003-2006 fertility, CPS June 2004/2006 parity, ACS 2005/2006 housing, and PSID 2003/2005 wealth. Keep the reviewed birth-housing response as a pooled structural restriction, with its original estimator and an explicit stability assumption. These windows describe an approximate early economy; none is relabeled as a literal 2007 cross-section.')
para('<b>2. Assess that fit before freezing parameters.</b> Show current birth rates, older cohorts\' completed fertility, wealth, ownership and housing together. A stationary approximation is acceptable only if its discrepancies are understood and do not drive the policy mechanism. A low optimizer loss alone cannot establish that.')
para('<b>3. Fit the preference decline with perfect foresight.</b> Households learn the entire path in 2007. Fit one amplitude to the 2020-2023 average period TFR, 1.64575, using the calendar mapping on page 6. Retain middle windows as checks and keep the fitted path fixed across policies. The existing linear 2007-2023 path with a flat continuation is a parsimonious first specification, not an estimated path shape. [10]')
para('<b>4. Solve matched policies from 2023.</b> The policy and baseline must inherit the same historical state and future preference path. The main property-tax comparison still needs its agreed equal-rebate treatment in both regimes, verified market clearing and horizon stability.')
heading('What this investigation changed')
para('<b>The early data are available.</b> We constructed housing and wealth candidates, including uncertainty. Earlier statements that these targets were unavailable were too strong: the existing production cache did not contain them. The remaining issues are measurement, the initial stationary approximation and numerical certification.')
para('The literature provides precedents for both initial calibration followed by a transition and joint estimation along a dated path. Neither forces us to fit a 2023 steady state. My preference for the staged design is driven by your policy question, the newly established early-data feasibility and its lower computational cost. [1-5]','small')

page('The economic choice and its main limitation')
para('Two decisions have been getting conflated. <b>The initial-state assumption</b> determines the distribution of wealth, parity and housing inherited in 2007. <b>The estimation design</b> determines whether later observations can change the common structural parameters. Choosing joint estimation does not automatically change the initial-state assumption.')
table([['Design','What it does','Assessment'],
 ['Initial fit, then shock','Fit common parameters to early moments; fit the shock on the transition.','Recommended first route. Fast initial diagnosis and a clear separation between benchmark fit and shock fitting.'],
 ['Joint dated transition fit','Fit common parameters and the shock together, comparing every observation with its own date/cohort/window.','Coherent alternative if later outcomes must discipline parameters. A larger estimator; the same stationary initial-state restriction remains.'],
 ['Empirical initial state','Construct the initial joint distribution from data or explicit statistical matching; estimate behavior from flows and subsequent evolution.','Addresses missing cohort histories directly. Requires additional state construction and a new identification design.']],[1.1,2.2,2.5])
heading('The fertility discrepancy is real, but it is not a contradiction in the data')
para('Period fertility summarizes births under current age-specific rates. Completed fertility sums births experienced by a particular cohort over its life. In the early data, the period index is <b>2.0605</b>; women aged 40-44 in the 2004/2006 CPS report <b>1.8566</b> births capped at five, or <b>1.8784</b> uncapped. They are different women and histories. [6,9; local fertility receipt]')
para('Under matching populations, exposures and age ranges, an age-specific schedule constant over time makes the full period and cohort integrals equal:',markdown='Under matching populations, exposures and age ranges, a time-invariant fertility schedule makes the full period and cohort integrals equal:')
para('TFR(t) = ∫ f(a,t) da;   CF(b) = ∫ f(a,b+a) da.','body',r'\[TFR(t)=\int f(a,t)\,da,\qquad CF(b)=\int f(a,b+a)\,da.\]')
para('This identity does not literally equate a capped stock at ages 40-44 with full lifetime births. But removing the cap closes only 0.0218 of the gap. At frozen early rates, even the entire remaining age-40-49 fertility integral is only 0.0479 births. We cannot reconcile the gap by capping and incomplete completion alone under that schedule. Nor have we established that postponement explains it.')
para('<b>My proposed approximation is period-oriented.</b> Anchor the initial current fertility level; use early childlessness, one-child share and birth timing as additional discipline; always display old completed fertility and higher-parity outcomes as checks. The parity restrictions may force an excessive higher-parity tail. This is a testable risk, not a promise that the initial fit will work.')
para('If that approximation fails on the states that matter for housing and policy, the substantive fallback is a historically heterogeneous initial state. A joint fit with the same stationary distribution changes the compromise; it cannot manufacture the missing pre-2007 histories. [3,6]','small')

page('The complete proposed initial target set')
para('Each row states the main parameter discipline, not exclusive identification. Parentheses contain newly measured uncertainty where available. These are candidate observations, not active calibration weights. All housing values use the same 42 MET2013 city codes, a new sample definition explained on the next page.','small')
target_rows=[['Parameter / restriction','Empirical moment and window','Candidate value (SE)'],
 ['Initial child preference ψ₀','Period TFR, NCHS 2003-2006; proposed initial anchor','2.060500; weight not set'],
 ['First-birth cost F','Childless, women 40-44, CPS 2004/2006','0.198279; SE pending'],
 ['First-birth dispersion κE','Mean first-birth age in model midpoint bins, period 2003-2006','25.976264; SE/scale pending'],
 ['First-birth dispersion κE','First births at 30+, same period pool','0.249278; SE/scale pending'],
 ['Further-birth dispersion κC','Exactly one child among mothers 40-44, CPS 2004/2006','0.213655; SE pending'],
 ['Patience β','Total net worth / gross labor earnings, PSID 2003/2005',point_se(psid_ratio)],
 ['Bequest strength θ₀','Annual bequests / wealth; inherited external restriction','0.008800; historical external'],
 ['Bequest wealth shift θ₁','Old wealth/income p90/median, ages 76-84, PSID 2003/2005',point_se(psid_tail)],
 ['Housing supply scale H₀','Mean min(rooms,9), heads 18-85, ACS 2005/2006',point_se(h('aggregate_mean_occupied_rooms_capped9_18_85'),'point','metro_bootstrap_se')],
 ['Owner preference χ','Ownership, heads 30-55, standard structures, ACS 2005/2006',point_se(h('own_rate_30_55'),'point','metro_bootstrap_se')],
 ['First-child housing jump hJ','Reviewed first-birth PSID Sun-Abraham contrast, -1 to +3','0.720246 (0.085260) rooms'],
 ['Per-child housing floor hC','3+ minus 1-2 resident-child rooms, ACS 2005/2006',point_se(h('prime30_55_resident_3plus_minus_1to2_rooms_capped9'),'point','metro_bootstrap_se')],
 ['Joint housing/fertility restriction','Recent-parent minus no-resident-child ownership, ACS 2005/2006',point_se(h('recent_parent_minus_no_resident_child_ownership_30_55'),'point','metro_bootstrap_se')]]
table(target_rows,[1.3,2.55,1.6])
para('<b>Count:</b> 13 restrictions for 10 common structural parameters plus the initial preference level. If the period-rate anchor recovers ψ₀ internally, the remaining objective has 12 rows for 10 searched coordinates. This passes only the count requirement; the local weighted Jacobian must establish informative variation. The shock amplitude is fitted subsequently, not in this initial objective.','small')
para('<b>Keep visible as initial checks:</b> CPS completed fertility 1.856608 capped at five; young ownership 25-34 = '+fm(h('own_rate_25_34')['point'])+' (SE '+fm(h('own_rate_25_34')['metro_bootstrap_se'])+'); old wealth/income median '+point_se(psid_level)+'. The first two diagnose fertility history and young housing access. The last is a candidate extra bequest restriction if the existing block is weak; adoption requires a new rank check.','small')

page('What the new data permit, and what must be aligned')
heading('A genuine early housing sample, with a declared geographic change')
para('The new ACS construction uses the same 42 city identifiers in early and later years. It does <b>not</b> preserve the current custom set of admitted geographic areas within those cities. Restoring the old filter reproduces all four current 2023 targets to numerical precision. The old filter excludes 48.9% of the full 42-city head weight in 2012, and 4.54% in 2023. The new construction therefore requires a new target contract and consistent remeasurement of later comparisons.')
para('The IPUMS city assignment itself is approximate because public-use areas change. A fixed city code is not an exact fixed land area. ACS also removed the national nine-room topcode in 2008. For comparable room-level observations, apply min(rooms,9) to both data and model <i>before averaging</i>; this changes the observation rule, not the model\'s housing choices. The reviewed uncapped PSID event-study response retains its own observer. [7,8]')
table([['Same proposed observation rule','ACS 2005/2006','ACS 2023'],
 ['Ownership, heads 30-55','64.8334%','58.7409%'],
 ['Ownership, heads 25-34','43.1158%','35.6027%'],
 ['Mean rooms, capped at nine','5.561097','5.583723']],[2.5,1.2,1.2])
para('These differences are descriptive, not causal estimates. They show why later ownership levels cannot simply be used as initial levels under an unnamed pooling assumption. The new metro-bootstrap errors describe resampling 42 cities; they are not official ACS replicate-weight survey errors. Full covariance and paired draws are saved.','small')
heading('Wealth is measurable before the shock')
para('The PSID 2003/2005 pool supplies wealth/earnings and old wealth dispersion with fresh person-cluster bootstrap uncertainty. Both the original long-pool point estimates and their bootstrap errors reproduce the authoritative builders. Wealth/earnings is 6.146 in that early pool, 6.452 in 2005 alone and 7.364 in 2007: the window is economically consequential. Use the saved alternatives as initial-state sensitivities. These are living-household stocks, not estates; distinguish survey-wave dates from income reference years.')
heading('Identification changes must be explicit')
para('The staged proposal replaces later fertility discipline with early period level, timing and parity shape for ψ₀, F, κE and κC; remeasures housing levels and family groups for the housing block; and replaces long-pool wealth/dispersion with early observations for the saving/bequest block. The reviewed childbirth response and external bequest-flow restriction remain. Parameters are jointly disciplined within these blocks. Later outcomes become stated validation under this <i>new</i> design; that has not happened to the live objective.')
para('The Sun-Abraham coefficient is an empirical housing response, not itself a deep preference parameter. Applying it to the initial economy requires an explicit assumption that this response transfers across periods; its estimator, sample and reviewed value remain unchanged.','small')
para('The supply elasticity 0.63, tenure dispersion 0.005 and child-invariant bequests remain maintained restrictions. The recent CEX preference estimate α₀ = 0.733 requires an explicit time-stability assumption; the source of 0.63 remains unresolved. The bequest wealth shift is weakly disciplined in existing evidence. No sensitivity matrix for the proposed system has yet established identification.','small')
para('National fertility and PSID restrictions combined with metro housing remain a maintained geographic approximation. Choosing this new metro target sample does not by itself certify the population or geographic closure of policy experiments.','small')

page('What the primary literature actually establishes')
para('The methodological lesson is to match an observation to the model object that generated it. There is no universal requirement that every estimate come from one calendar year. A maintained preference elasticity and an equilibrium ownership level have different reasons for being pooled. The following distinctions were checked in primary texts.','small')
table([['Paper','Verified method','What it supports here'],
 ['Sommer, Sullivan & Verbrugge (2013), JME [1]','Section 3.5 calibrates four parameters to stationary moments from several vintages. Section 6 starts from an initial steady state, applies an unexpected permanent change in rates/downpayments, then solves perfect foresight.','A close housing precedent for initial calibration plus transition. Their experiment is stylized; it does not jointly estimate a complete annual historical path.'],
 ['Greenwood, Seshadri & Vandenbroucke (2005), AER [2]','Section III.B jointly fits preference/technology parameters and dated household-technology levels to fertility along 1800-1990; market productivity is supplied.','A genuine deterministic fertility-path estimation precedent. Joint fitting is legitimate, but their causal technology story and simpler OLG structure differ from ours.'],
 ['De Nardi, French & Jones (2010), JPE [3]','Section IV initializes from the observed 1996 state distribution, then matches cohort/age/income-specific wealth profiles using simulated moments.','Observed initial distributions are a coherent alternative to stationary initialization. This is a retiree lifecycle model, not our fertility/housing GE transition.'],
 ['Borella, De Nardi, Pak, Russo & Yang (2023), JEEA [4]','Section 6 and Appendix C.9 estimate a cohort lifecycle model under dated tax regimes with explicit perfect foresight; 19 parameters, 448 moments.','Calendar-specific observations can enter one objective. Taxes are measured inputs; this is not joint estimation of an endogenous aggregate price path.'],
 ['Kaplan, Mitman & Violante (2020), JPE [5]','Section III calibrates a stochastic ergodic economy and aggregate regime/belief processes; the boom-bust is a realization of a Markov process.','Useful for disciplined initial micro moments and aggregate drivers. Its expectations structure is not the deterministic announcement assumed here.'],
 ['Kohler & Ortega (2002), Demographic Research [6]','Sections 1 and 3.1 distinguish period/cohort fertility and age-parity exposure; adjusted period rates do not identify actual completed cohorts without additional assumptions.','The old CPS stock and early period index require different observers. A generic appeal to tempo effects does not resolve our measured discrepancy.']],[1.25,2.25,2.25])
para('<b>My judgment:</b> the staged design is a credible baseline for the conditional policy question. Joint dated estimation is a possible robustness exercise if later observations must help discipline common parameters. Neither is a shortcut around initial-state specification or imperfect empirical measurement. No literature claim here certifies the current implementation.','small')

page('Implementation, numerical status and a bounded next run')
heading('Four prerequisites before activating the proposed objective')
para('<b>Fertility measurement.</b> The current period diagnostic divides births by adult-household mass, not by correctly aligned female exposure. The annual population bridge does not yet provide the maternal age/parity mapping needed for an exact NCHS comparison. First-birth timing must use period birth flows, while CPS parity uses the relevant completed cohorts and age window. The model\'s inherited representative count for the 3+ bin also differs from the early CPS count.')
para('<b>Calendar mapping.</b> The annual bridge allocates a decision dated t to births in t+1,...,t+4. Under that convention, a 2020-2023 rate window belongs to the 2019 decision block, not the 2023 decision. The historical observer must implement and verify that convention; simply renaming its 2023 row would be wrong.')
para('<b>Initial demographic closure.</b> The existing 2.1 normalization also sets births-to-household conversion, initial birth queues and an entry-flow gate. Replacing it with a female period TFR is not a justified one-line substitution. The new fertility anchor and household renewal law must be specified separately. The initial supply construction also inherits elasticity 1.75 while the dated restriction is 0.63; its intended role must be reconciled.')
para('<b>Family groups and source contracts.</b> Recent-parent versus no-resident-child data are currently compared with different model groups; resident own children and model dependents are also different. Correct the observers, finish CPS/timing uncertainty, pin the new target/weight fingerprint, then test the exact candidate loop. No missing uncertainty should inherit an unrelated old standard error.')
heading('There has been real progress on the price solver')
para('The inherited-parameter sequential 100-date path now clears markets: maximum residual <b>0.003673%</b>, below the unchanged <b>0.020000%</b> gate. The root records a zero-difference final replay; collected bytes and fit arithmetic have been checked. The terminal unit-rent gap is still <b>1.088827%</b> against a <b>1%</b> requirement. Horizon stability, re-estimation and matched policies remain outstanding. Complete inherited fits and parameters are on the final page.')
table([['Work unit','Observed / conditional budget','Interpretation'],
 ['Normalized initial candidate','About 5 minutes; four stationary GE solves in the observed case','An initial 21-case sensitivity panel is about 1.7 core-hours, potentially one parallel wave. Not a full calibration.'],
 ['One 100-date policy/value mapping','42-56 minutes','The entire backward/forward path at supplied prices.'],
 ['One warm-started equilibrium candidate','Roughly 2.2-2.9 hours if three mappings suffice','Optimistic: three calls allow only one price update and a replay.'],
 ['23-case joint parameter panel','Roughly 51-68 core-hours under the same optimistic root assumption','Several hours with sufficient independent workers; one derivative round, not finished estimation.']],[1.3,1.7,2.65])
para('<b>Next numerical decision:</b> after the measurement contract is complete, run the small initial stationary panel and fit first; show every target and diagnostic, inspect parity tails and boundary policies, and assess local rank. Only then freeze the common parameters and launch the scalar-shock PF search. Smoke-test each loop, retain per-case checkpoints and stop on failed gates. This research review launched no new numerical jobs.','small')

page('Sources and reproducibility')
sources=[
 ('1','Sommer, K., P. Sullivan and R. Verbrugge (2013). The equilibrium effect of fundamentals on house prices and rents. Journal of Monetary Economics 60, 854-870. Sections 3.5 and 6; printed pp. 860-861 and 867-868.','https://www.kamilasommer.net/RentPriceRatio.pdf'),
 ('2','Greenwood, J., A. Seshadri and G. Vandenbroucke (2005). The Baby Boom and Baby Bust. American Economic Review 95(1), 183-207. Section III.B, printed pp. 189-190. Main transition exercise, distinct from Section IV\'s illustrative steady states.','https://www.jeremygreenwood.net/papers/bb.pdf'),
 ('3','De Nardi, M., E. French and J. B. Jones (2010). Why Do the Elderly Save? The Role of Medical Expenses. Journal of Political Economy 118(1), 39-75. Section IV, printed pp. 46-48.','https://users.nber.org/~denardim/research/De_Nardi_French_Jones_JPE_2010.pdf'),
 ('4','Borella, M., M. De Nardi, M. Pak, N. Russo and F. Yang (2023). FBBVA Lecture 2023. The Importance of Modeling Income Taxes over Time: U.S. Reforms and Outcomes. JEEA 21(6), 2237-2286. Section 6 and Appendix C.9.','https://academic.oup.com/jeea/article/21/6/2237/7275099'),
 ('5','Kaplan, G., K. Mitman and G. L. Violante (2020). The Housing Boom and Bust: Model Meets Evidence. Journal of Political Economy 128(9), 3285-3345. Section III, especially pp. 3302-3310.','https://violante.economics.princeton.edu/sites/g/files/toruqf5621/files/documents/kaplan-et-al-2020-the-housing-boom-and-bust-model-meets-evidence.pdf'),
 ('6','Kohler, H.-P. and J. A. Ortega (2002). Tempo-Adjusted Period Parity Progression Measures, Fertility Postponement and Completed Cohort Fertility. Demographic Research 6(6), 91-144. Sections 1 and 3.1, pp. 92-93 and 102-103.','https://www.demographic-research.org/volumes/vol6/6/6-6.pdf'),
 ('7','IPUMS USA. MET2013 documentation: geographic assignment, mismatch threshold and comparability across PUMA vintages. Accessed 10 September 2026.','https://usa.ipums.org/usa-action/variables/MET2013'),
 ('8','IPUMS USA. ROOMS documentation, Comparability: national topcode of nine removed in 2008. Accessed 10 September 2026.','https://usa.ipums.org/usa-action/variables/ROOMS'),
 ('9','NCHS. Births: Final Data for 2007, Table 4, printed p. 25. Source of 2003-2006 period rates; 2.0605 is our equal-year average.','https://www.cdc.gov/nchs/data/nvsr/nvsr58/nvsr58_24.pdf'),
 ('10','NCHS. Births: Final Data for 2023, Table 2, printed p. 14. The 2020-2023 average is 1.64575, versus 1.621 in 2023 alone; window averaging is our calculation.','https://www.cdc.gov/nchs/data/nvsr/nvsr74/nvsr74-1.pdf'),
]
for n,desc,url in sources: para('['+n+'] '+escape(desc)+' <link href="'+url+'" color="#222222"><u>Primary source</u></link>.','small')
heading('Local evidence')
para('All local evidence is indexed in <b>output/model/e5f_matched_pf_20260909a/design_research/README.md</b>. Housing components, bootstrap draws and source receipts are under housing/; wealth equivalents are under wealth/. The fertility extraction and initial coherence review remain under the adjacent parameter_target_audit/ folder. The active isolated model snapshot is 96a41873. The final numerical receipts are under computation/final_replay/.','small')
para('The lead checked primary papers and source definitions, reviewed empirical masks, independently recalculated saved components and uncertainty, and checked numerical receipt hashes and fit arithmetic. Claude supplied a public-literature review using two research agents in safe mode; its unverified leads were not promoted. No private project material was exported to Claude. Research statements are separated from unimplemented proposals throughout.','small')

page('Appendix: inherited parameters at the converged price path')
numerical=BASE/'computation/final_replay/evaluation_003'
summary=json.loads((numerical/'summary.json').read_text())
fits=rows(numerical/'target_fit.csv');params=rows(numerical/'parameters.csv')
LABELS={'tfr':'Completed fertility','childless_rate':'Childless share','mean_age_first_birth':'Mean first-birth age','share_first_birth_30plus':'First births at 30+','first_birth_housing_response':'First-birth room response'}
# Existing label map is a read-only report helper; import without running its main.
import sys
sys.path.insert(0,str(BASE.parent))
from review_horizon100_root_progress import LABELS as EXISTING_LABELS
LABELS.update(EXISTING_LABELS)
para(f'Objective {summary["loss"]:.8f}. These are the unchanged inherited inputs and current 12-row objective, <b>not estimates under the proposed early target set</b>. The finite price path converged; horizon certification remains open. Shares and gaps use fraction units.','small')
def num(v): return '-' if v=='' else f'{float(v):.6g}'
table([['Moment','Target','Model','Gap','Weight','Loss']]+[[LABELS.get(r['moment'],r['moment'])]+[num(r[k]) for k in ['target','model','gap','weight','loss_contribution']] for r in fits],[2.9,1,1,1,1,1],compact=True)
para('All 11 free coordinates are inherited rather than re-estimated. ψ₀ is normalized to old completed fertility 2.1, which is distinct from the proposed female period-rate anchor. Bounds and near-bound flags reproduce the existing receipt.','small')
short={'beta_annual':'β','kappa_fert':'κE','kappa_fert_continuation':'κC','chi':'χ','H0':'H₀','theta0':'θ₀','theta1':'θ₁','hbar_child_rooms':'hC','first_birth_fixed_cost':'F','hbar_first_child_jump':'hJ','psi_child_change_2023':'Δψ','psi_child_2007':'ψ₀','psi_child_2023':'ψ2023','tenure_choice_kappa':'Tenure dispersion','housing_supply_elasticity':'Supply elasticity'}
prows=[['Parameter','Value','Lower','Upper','Role','Near bound']]
for r in params:
    role='Free input' if r['is_free_parameter']=='True' else 'Normalized' if r['parameter']=='psi_child_2007' else 'Derived' if r['parameter']=='psi_child_2023' else 'Fixed'
    prows.append([short.get(r['parameter'],r['parameter']),num(r['value']),num(r['lower_bound']),num(r['upper_bound']),role,'Yes' if r['near_bound']=='True' else 'No'])
table(prows,[1.7,1.4,1,1,1.5,1.1],compact=True)
para('The full parameter names, unrounded values, every target weight and loss contribution are preserved in the source CSVs. No inference about fit improvement across target systems is made. Existing ACS date/group problems still apply to these inherited rows; see the main report before interpreting them.','small')

def footer(canvas,doc):
    canvas.setFont('Arial',8)
    canvas.setFillColor(colors.HexColor('#555555'))
    canvas.drawRightString(A4[0]-44,23,str(doc.page))
OUT.parent.mkdir(parents=True,exist_ok=True)
SimpleDocTemplate(str(OUT),pagesize=A4,title='A pre-decline benchmark, then a fertility shock',author='Fertility research project',leftMargin=44,rightMargin=44,topMargin=35,bottomMargin=36).build(story,onFirstPage=footer,onLaterPages=footer)
(BASE/'DECISION_REPORT.md').write_text('\n\n'.join(md)+'\n')
print(OUT)
