"""Diagnostic: sale-to-rent screen with current income (R b + S + y >= 0) vs adopted (R b + S >= 0), revolving rule, LTV 80 and 95."""
import sys; from pathlib import Path
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): import analyze as A
import numpy as np
C = {c: A.load(c) for c in ['R80', 'R95', 'R80I', 'R95I']}; O = {c: A.outcomes(S) for c, S in C.items()}
def imp(a, b):
    pre = C[a]['birth_count_pre_distribution']; return 100 * float((pre * (C[b]['f'] - C[a]['f'])).sum()) / O[a]['B']
rows = []
for lab, a, b in [('Adopted screen (no income)', 'R80', 'R95'), ('Screen with income (diagnostic)', 'R80I', 'R95I')]:
    oa, ob = O[a], O[b]
    rows.append([lab, f'{imp(a, b):+.2f}%', f"{oa['ceb']:.3f} -> {ob['ceb']:.3f} ({100*(ob['ceb']/oa['ceb']-1):+.2f}%)",
                 f"{100*oa['childless']:.1f} -> {100*ob['childless']:.1f}%", f"{100*oa['own1829']:.1f} -> {100*ob['own1829']:.1f}% ({100*(ob['own1829']-oa['own1829']):+.1f} pp)",
                 f"{100*oa['own']:.1f} -> {100*ob['own']:.1f}%", f"{oa['rooms']:.3f} -> {ob['rooms']:.3f} ({100*(ob['rooms']/oa['rooms']-1):+.2f}%)"])
lvl = [f"Screen effect at LTV 80 (R80I vs R80): births on impact {imp('R80','R80I'):+.2f}%, completed {100*(O['R80I']['ceb']/O['R80']['ceb']-1):+.2f}%; "
       f"at LTV 95 (R95I vs R95): births on impact {imp('R95','R95I'):+.2f}%, completed {100*(O['R95I']['ceb']/O['R95']['ceb']-1):+.2f}%."]
hdr = ['LTV 80 -> 95', 'births on impact', 'completed fertility 46-49', 'childless 46-49', 'ownership 18-29', 'ownership all', 'mean rooms']
md = '| ' + ' | '.join(hdr) + ' |\n|' + '---|' * len(hdr) + '\n' + ''.join('| ' + ' | '.join(r) + ' |\n' for r in rows)
txt = ('DIAGNOSTIC on the copied engine tmp/origination_only_20261008 (P.sale_screen_includes_income); not a production change. '
       'Fixed price 0.779414, rebate T held, base 14.402, revolving stayer rule.\n\n' + md + '\n' + lvl[0] + '\n')
(HERE / 'sale_screen_income_diagnostic.md').write_text(txt); print(txt)
