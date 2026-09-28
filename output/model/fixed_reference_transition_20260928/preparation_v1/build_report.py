"""Render the readiness memo and complete reference tables on Torch; no model imports."""
from pathlib import Path
import json
import os
from xml.sax.saxutils import escape
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether

HERE = Path(__file__).resolve().parent
ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
M = json.loads((ROOT/'output/model/fertility_identification_20260928/fixed_reference_manifest.json').read_text())
OUT = HERE/'report'
OUT.mkdir(exist_ok=True)
navy, teal = colors.HexColor('#172d40'), colors.HexColor('#216a78')
body = ParagraphStyle('body', fontName='Helvetica', fontSize=9.1, leading=12.2, textColor=navy, spaceAfter=7)
small = ParagraphStyle('small', parent=body, fontSize=8, leading=10.2, spaceAfter=5)
title = ParagraphStyle('title', parent=body, fontName='Helvetica-Bold', fontSize=20, leading=24, spaceAfter=12)
heading = ParagraphStyle('heading', parent=body, fontName='Helvetica-Bold', fontSize=11, leading=14, textColor=teal, spaceBefore=8, spaceAfter=6)
cell = ParagraphStyle('cell', parent=small, fontSize=7.8, leading=10, spaceAfter=0)
story = []

def p(text, style=body):
    story.append(Paragraph(text, style))

def h(text):
    p(text, heading)

def t(rows, widths, fs=None):
    style = cell if fs is None else ParagraphStyle('cell'+str(fs), parent=cell, fontSize=fs, leading=fs+2)
    items = [[Paragraph(escape(str(v)), style) for v in row] for row in rows]
    table = Table(items, colWidths=widths, repeatRows=1, hAlign='LEFT')
    table.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'), ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8f1f3')),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f9')]),
        ('LINEBELOW',(0,0),(-1,0),0.7,teal), ('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),
        ('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]))
    story.extend([table, Spacer(1,8)])

def page():
    story.append(PageBreak())

def number(x):
    if x is None or x == '': return '-'
    x = float(x)
    if x != 0 and (abs(x)<.001 or abs(x)>=10000): return f'{x:.2e}'
    return f'{x:.3f}'.rstrip('0').rstrip('.')

summary_path = HERE/'runs/18738593/suite_result.json'
test = json.loads(summary_path.read_text()) if summary_path.exists() else {'status':'PENDING','completed':[]}
completed = test.get('completed',[])
one_path = HERE/'runs/18739319/suite_result.json'
one = json.loads(one_path.read_text()) if one_path.exists() else {'status':'PENDING','completed':[]}

p('Transition readiness', title)
p('<b>2007 stationary reference — block0506, September 28 verified export</b><br/>September 28, 2026 | Preparation and validation; no policy effect claimed', small)
h('Decision')
p('The dated household and population code can be reused, but the previous transition launcher cannot be run unchanged. It selects an older calibration and an approximate credit rule. The new isolated bridge authenticates all 260 saved parameter fields and imports the frozen September 28 source. Preferences and measurements remain fixed.')
p('Your two definitions are recorded: <b>fixed physical housing stock, with prices and rents clearing</b>; and <b>removal of artificial borrowing/down-payment limits while preserving repayment and lifetime solvency</b>. Individual housing and tenure choices remain free.')
h('Numerical preparation')
rows = [['Check','Status','Evidence']]
control = completed[0] if completed else None
rows.append(['Exact reference control', 'PASS' if control else 'Unverified',
             '113 exact arrays; 14 fit rows, 31 parameters, 17 identical plots. Imported control; no new solve.'])
rows.append(['Two-date no shock','Timed out','180-second cap reached before the first dated audit. No operator acceptance.'])
one_completed = one.get('completed',[])
if one_completed:
    r=one_completed[0]
    rows.append(['One-date diagnostic',r['status'],f"{r['bellman_solves']} calls; {r['seconds']:.1f} s; housing {r['market_residual']:.3g}; PAYGO {r['fiscal_residual']:.3g}."])
else:
    rows.append(['One-date diagnostic',one['status'],'Instrumented, two-call maximum; six-minute cap. Preserved receipt gives the stopping point.'])
rows.append(['Six dates / fixed stock','Not run','Longer demographic and fixed-stock operator gates remain open.'])
t(rows,[122,68,350])
p('The source project is read-only. Earlier failures exposed an outdated serializer and a diagnostic path pointing into the reference. The bridge now checks all fields itself and redirects only diagnostic output. No shocked equilibrium or natural-credit implementation is certified.',small)
h('Meaning of the old estimated transition')
p('The short historical fit estimated four successive unexpected preference changes: after each surprise households expected its new level to last forever. A separate announced-path experiment revealed all four changes at once. These are different expectations assumptions.')
t([['Surprise','Fertility window','Old preference','Target','Short fit'],
   ['2007','2008-2011','0.129','1.975','1.973'],['2011','2012-2015','0.117','1.861','1.861'],
   ['2015','2016-2019','0.106','1.755','1.755'],['2019','2020-2023','0.092','1.646','1.646']], [61,124,113,121,121])
p('The subsequent long-horizon refit accepted <b>zero fitted shocks</b>. One 104-date market/fiscal root passed but missed fertility and remained 1.742% away in population and 9.446% in distribution distance from its endpoint. Old absolute preference levels cannot silently transfer to the new utility specification.',small)

page()
p('Closures that must remain explicit',title)
p('Estimated and normalized inputs are frozen at their saved values. An endogenous outcome is solved after a shock; an outstanding closure is not silently replaced by a convenient default.',small)
t([['Object / classification','Maintained treatment'],
 ['Preferences / estimated or calibrated, then fixed','All saved preferences, including child-benefit scale. No fertility renormalization after a credit shock.'],
 ['Earnings / externally estimated','Approved B15 Markov process, 15 income states and saved age profile. Reuse the completed measurement audit.'],
 ['Population / normalized initial state','Exact saved pre-choice distribution, mass one. Carry level-valued population forward; no reset or age reweighting.'],
 ['Adult entry / author-fixed normalization','Half of a birth vintage enters after 16 years, half after 20; adjusted births map to households at 1/2.1. Keep raw and adjusted queues.'],
 ['Entry wealth / empirical normalization','Saved conditional wealth/income distribution; retain negative entry assets and inherited rank coupling. No zero-wealth substitute.'],
 ['Survival and geography / fixed','Saved age survival; one pooled nationally calibrated market. Closed diagnostic path, with no new migration, retention or old quota defaults.'],
 ['PAYGO / empirical baseline, fixed tax','Hold payroll tax at 8.028%; solve dated pension benefits from actual worker and retiree masses. Baseline benefit 0.918 is a starting value.'],
 ['Property tax / externally fixed','Saved annual rate 1.060%, period rate 4.239%; zero household rebate. Do not import old equal rebates.'],
 ['Estate settlement / outstanding','Retain provisional net-estate funding of actual next entrants and residual sink. Shortfalls and negative estates fail. Counterparty and physical settlement remain unresolved.'],
 ['Housing / author-defined experiments','Credit experiment retains the absolute elastic supply curve. Fixed-stock arm holds actual reference supply (5.848), not the intercept (6.294). Replacement of depreciation is implicit.'],
 ['Terminal population / endogenous','Solve renewal, PAYGO and housing jointly. A unit-mass stationary distribution alone is not a demographic equilibrium.']], [155,385])
h('Natural solvency is a separate implementation gate')
p('Construct feasible states backward at each decision node, preserving the existing information and fertility/housing timing. Every reachable future income/family outcome must permit repayment. At each possible death, post-saving liquid wealth plus net housing liquidation value must be nonnegative. The same condition applies at terminal age. With positive death risk and no default or life insurance, it can itself exclude unsecured borrowing; positive future earnings do not guarantee repayment after an early death.',small)
p('The old option infers feasibility from a value cutoff and the first feasible grid node. It conflicts with the saved incumbent-owner debt rule. Replace the artificial purchaser and incumbent limits explicitly, then verify worst-income/death boundaries and an expanded/refined debt grid. Financed share equal to one is insufficient.',small)

page()
p('Architecture, gates and finite budgets',title)
h('Endpoint and path are different problems')
p('At fixed preferences, solve the price and pension for demographic renewal and balanced PAYGO. With B adjusted births, E entrant households and d housing demand per normalized household, renewal requires B/(2.1 E) = 1. Housing then determines the population level N = H<super>S</super>(q)/d; under fixed stock, replace H<super>S</super>(q) by the reference stock. Verify native one-step distribution and both entry queues. Do not restore fertility mechanically if a root is absent.')
t([['Calculation','Interpretation'],
 ['Impact / prescribed continuation','Household response from the exact inherited population. No market-clearing claim.'],
 ['Heuristic short path','A stated guessed price/pension sequence, with population carried forward. Not an equilibrium price path.'],
 ['Cleared finite path','Households, housing and pensions clear at every included date. Endpoint and horizon adequacy still require verification.'],
 ['Genuine transition','Finite clearing plus terminal distribution/queue approach and stable early responses under horizon extension.']], [145,395])
p('A pure 10% supply-intercept increase has an algebraic endpoint candidate: unchanged prices, pension and per-household policies, with population and all aggregate flows 10% higher. Renewal, PAYGO and the provisional estate ledger scale together. Native scaling checks remain required; this does not establish uniqueness, stability or a transition.',small)
h('Bounded sequence')
t([['Stage','Maximum / stop condition'],
 ['Current preparation','45 active minutes, excluding queue time. Serial Torch checks, one CPU/24 GiB. After a 180-second two-date timeout, one instrumented date is capped at six minutes and two calls. No policy run.'],
 ['Natural-solvency verification (proposed)','12 constructed boundary cases plus two fixed-price household solves, 20 minutes. Stop at support/repayment disagreement.'],
 ['One credit endpoint (proposed)','16 stationary evaluations including repeat/one-step checks; 45 minutes. No preference search. Stop without a valid renewal bracket or at budget.'],
 ['Six-date cleared diagnostic (proposed)','First smoke six dates (about 11 min; 15-min cap). Then at most eight root mappings plus one replay: 108 calls, about 94 min; proposed 2-hour cap. Count Jacobian work separately.'],
 ['Historical fixed stock / long horizon','Not launched. Requires an agreed shock contract and a separate finite plan. The other chat owns the +10% supply comparison.']], [156,384])
p('Use measured replay/mapping time before approving the next budget. The independent fixed-price control took about 66 seconds. The six-date mapping crosses both maturation lags; the initial entry mismatch is only 4.889e-8 households, but subsequent drift must be reported rather than repaired by normalization.',small)
h('Acceptance and next decision')
p('Preserve source/checkpoint identities, all preference fields, exact initial mass, occupied household budgets and debt feasibility, probabilities, value monotonicity, estate funding, and zero inherited-state projection. Retain housing tolerance 2e-4, PAYGO 1e-6, mass accounting 2e-8 and backward/forward reproduction 1e-10. Keep the 17 standard plots and add dated price, rent, quantities, population and boundary-distance panels.',small)
p('For the historical fixed-stock comparison, the author has been asked whether to leave historical shocks pending or use old proportional preference declines only as an illustrative diagnostic. No old shocks, new estimates, preferences or fiscal closures have been adopted here.',small)
p('Traceable sources: transition_readiness.md; historical_user_evidence.jsonl; historical_source_identity.json; the frozen reference manifest; September 13 original-queue and September 16 recovery receipts. The exact output paths and commands are indexed in the packet README.',small)

page()
p('Complete reference target fit',title)
p('Reference loss 19.581. These are unchanged reference observations, not counterfactual results. All 14 rows are retained, including the normalization and three zero-weight validation rows.',small)
labels={'initial_normalization':'Completed fertility', 'cps_childlessness':'Childlessness', 'cps_exactly_one':'Exactly one child',
 'nchs_mean_age':'Mean first-birth age','nchs_share30':'First birth at 30+', 'wealth_earnings':'Wealth / earnings',
 'bequest_wealth':'Bequests / wealth','old_dispersion':'Old-age dispersion', 'mean_rooms':'Mean rooms',
 'family_rooms':'Family rooms gap','recent_parent_ownership':'Recent-parent ownership gap','early_fertility':'Children at age 25',
 'own_rate_3055':'Ownership ages 30-55','ownership_30_55':'Ownership ages 30-55','first_birth_rooms':'First-birth rooms response'}
rows=[['Moment','Target','Model','Gap','Weight','Loss','Role']]
for r in M['full_target_table']:
    rows.append([labels.get(r['moment'],r['moment'].replace('_',' ')),*[number(r[k]) for k in ('target','model','gap','weight','loss_contribution')],r['role']])
t(rows,[150,55,55,65,75,65,75],7.8)
p('The full-precision target table is preserved in the selected export. Loss comparisons require this exact target/weight system. Completed fertility was a calibration normalization; credit-policy evaluation will display its gap instead of imposing it.',small)
h('Identity and interpretation')
p('Checkpoint: b15ba92d...7309d. Source manifest: 07d84336...496d. Identification contract: 68323aad...5abf. Full hashes and all 260 serialized fields are in fixed_reference_manifest.json. Two exact repetitions and the standard plot identities were authenticated before this preparation.',small)
p('Existing limitations remain: high-wealth ownership and age-30 housing profiles, buyer-conditional policy plots, held-out earnings covariance miss, and provisional estate/entry interpretation. Reproduction verifies implementation identity; it does not establish global calibration optimality or resolve these economic issues.',small)

page()
p('Estimated and normalized quantities',title)
p('All are frozen for the requested credit experiment. Search bounds describe the inherited calibration; no search takes place in this preparation. A near-bound flag means within 1% of the inherited search interval.',small)
rows=[['Parameter','Estimate','Bounds / restriction','Near bound']]
for r in M['full_parameter_table'][:11]:
    bounds = f"[{number(r['lower'])}, {number(r['upper'])}]" if r['lower'] else 'Replacement normalization; now frozen'
    rows.append([r['parameter'].replace('_',' '),number(r['estimate']),bounds,{'True':'Yes','False':'No'}.get(r['near_bound'],r['near_bound'] or '-')])
t(rows,[180,77,211,72],8.4)
p('The two fertility taste scales are near their lower search bounds. This is an inherited calibration caveat, not permission to drop targets, change weights or renormalize preferences after a shock.',small)
h('Diagnostics')
p('The 17 standard PNGs remain unchanged in the selected-export standard_diagnostics folder. The imported exact control reproduces every PNG hash. This readiness memo does not replace that full diagnostic packet; new solved cases must produce it again, with supplemental transition panels clearly labeled.',small)

page()
p('All remaining reference parameters',title)
p('Continuation of the complete 31-row parameter table. Values, restrictions and meanings come from the saved export; constructor defaults are not used.',small)
rows=[['Parameter','Value','Inherited restriction / meaning']]
for r in M['full_parameter_table'][11:]:
    rows.append([r['parameter'].replace('_',' '),number(r['estimate']),r['status']])
t(rows,[195,72,273],8.2)
p('In a counterfactual, a retained financed-share number is metadata if its operative artificial constraint is replaced. Payroll tax stays fixed but the pension is solved from dated PAYGO. Housing supply changes only in the explicitly specified fixed-stock arm.',small)

def footer(canvas, doc):
    canvas.setStrokeColor(teal)
    canvas.line(36,36,576,36)
    canvas.setFont('Helvetica',7)
    canvas.setFillColor(navy)
    canvas.drawString(36,24,'September 28, 2026 | Frozen block0506 | Preparation; no policy effect certified')
    canvas.drawRightString(576,24,str(doc.page))

pdf = OUT/'transition_readiness.pdf'
SimpleDocTemplate(str(pdf), pagesize=letter, leftMargin=36,rightMargin=36,topMargin=38,bottomMargin=48,
    title='Transition readiness - frozen block0506',author='Fertility research project').build(story,onFirstPage=footer,onLaterPages=footer)
print(json.dumps(dict(pdf=str(pdf),target_rows=len(M['full_target_table']),parameter_rows=len(M['full_parameter_table']),two_date_status=test['status'],one_date_status=one['status'])))
