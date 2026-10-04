"""Three-page comparison from a validated saved solution and cached survey data."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
PSID = ROOT / 'code/data/psid_followup_mar2026/output/model_assessment/psid_2005_2007.csv'
CPS24 = ROOT / 'code/data/cps_fertility/cache/jun24pub.csv'
ACS = ROOT / 'output/model/native_financing_diagnostic_20260919/specification_followup/housing_profiles_v1/full/housing_profile_by_age.csv'
CPS_AVAIL = ROOT / 'output/model/e5f_matched_pf_20260909a/parameter_target_audit/fertility/fertility_availability.json'
NCHS = ROOT / 'code/data/nchs_natality_timing/first_birth_counts_year_age.csv'
BLUE, ORANGE = '#246aa4', '#d47a2a'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1 << 20), b''):
            h.update(part)
    return h.hexdigest()


def require_sha256(path, expected, label):
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(f'{label} SHA-256 mismatch: expected {expected}; found {actual}')
    return actual


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def weighted_ecdf(values, weights):
    v, w = np.asarray(values, float).ravel(), np.asarray(weights, float).ravel()
    if v.size != w.size:
        raise ValueError('ECDF values and weights differ in length')
    ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if not ok.any():
        raise ValueError('ECDF has no finite positive mass')
    v, w = v[ok], w[ok]
    order = np.argsort(v, kind='stable')
    v, w = v[order], w[order]
    starts = np.r_[0, np.flatnonzero(np.diff(v)) + 1]
    mass = np.add.reduceat(w, starts)
    return v[starts], np.cumsum(mass) / mass.sum()


def thirds(category_mass):
    """Fractionally split tied ordered income categories at exact weighted thirds."""
    mass = np.asarray(category_mass, float)
    total = mass.sum()
    if total <= 0:
        raise ValueError('Income groups have zero mass')
    lower = np.cumsum(mass) - mass
    upper = np.cumsum(mass)
    return np.array([
        np.maximum(0, np.minimum(upper, (k + 1) * total / 3)
                   - np.maximum(lower, k * total / 3)) / np.maximum(mass, 1e-300)
        for k in range(3)
    ])


def within_cell_stock(pre, post, cell_starts, ages, weights, period_years):
    """Interpolate pre/post stocks in the same cell at actual interview ages."""
    ages, weights = np.asarray(ages, float), np.asarray(weights, float)
    starts = np.asarray(cell_starts, float)
    idx = np.floor((ages + .5 - starts[0]) / period_years).astype(int)
    if np.any(idx < 0) or np.any(idx >= len(starts)):
        raise ValueError('Interview age outside model support')
    frac = (ages + .5 - starts[idx]) / period_years
    return np.average((1 - frac) * np.asarray(pre)[idx] + frac * np.asarray(post)[idx],
                      weights=weights, axis=0)


def invert_independent_birth_map(post, first, later, fecundity, fertile_indices):
    post = np.asarray(post, float)
    pre = post.copy()
    rate = np.zeros_like(post)
    if post.ndim != 7 or post.shape[-2:] != (4, 4):
        raise ValueError('Expected seven-axis independent-count model with 4x4 child states')
    for j in fertile_indices:
        rate[:, :, :, j, :, 0, 0] = fecundity[j] * first[:, :, :, j, :, 1]
        for n in (1, 2):
            for m in range(n + 1):
                rate[:, :, :, j, :, n, m] = fecundity[j] * later[:, :, :, j, :, 1, n - 1, m]
        for n in range(4):
            for m in range(n + 1):
                inflow = (0 if n == 0 or m == 0 else
                          pre[:, :, :, j, :, n - 1, m - 1] * rate[:, :, :, j, :, n - 1, m - 1])
                denominator = 1 - rate[:, :, :, j, :, n, m]
                if np.any(denominator <= 0):
                    raise ValueError('Birth map cannot be inverted')
                pre[:, :, :, j, :, n, m] = (post[:, :, :, j, :, n, m] - inflow) / denominator
    rebuilt = pre * (1 - rate)
    for n in range(1, 4):
        for m in range(1, n + 1):
            rebuilt[..., n, m] += pre[..., n - 1, m - 1] * rate[..., n - 1, m - 1]
    if pre.min() < -1e-12 or np.max(np.abs(rebuilt - post)) > 1e-10:
        raise ValueError('Saved birth-map inversion fails mass identity')
    return pre, rate


def financial_positions(net_worth, gross_home_value):
    wealth = np.asarray(net_worth, float)
    return wealth - np.asarray(gross_home_value, float), wealth


def select_grouped_national_acs(source_rows):
    selected = [r for r in source_rows if r['geography_scope'] == 'national'
                and r['sample'] == 'all_structures' and r['age_kind'] == 'annual']
    if not selected:
        raise ValueError('No annual national all-structures ACS rows')
    return selected


def baseline_cps():
    manifest = json.loads(CPS_AVAIL.read_text())
    require_sha256(Path(manifest['schema_source']), manifest['schema_sha256'], 'CPS schema')
    source = Path(manifest['cps_source'])
    open_source = source.open if source.exists() else lambda mode: gzip.open(str(source) + '.gz', mode)
    records, pins = [], {}
    with open_source('rb') as stream:
        for year in (2004, 2006):
            part = manifest['partitions'][str(year)]
            stream.seek(part['byte_start'])
            raw = stream.read(part['byte_end_exclusive'] - part['byte_start'])
            digest = hashlib.sha256(raw).hexdigest()
            if digest != part['partition_sha256']:
                raise RuntimeError(f'CPS {year} partition SHA-256 mismatch')
            pins[str(year)] = digest
            for offset in range(0, len(raw), 261):
                line = raw[offset:offset + 261]
                if len(line) != 261 or line[9:11] != b'06' or line[148:149] != b'2':
                    continue
                age, count, weight = int(line[146:148]), int(line[238:241]), int(line[250:260]) / 10000
                if 18 <= age <= 49 and 0 <= count <= 20 and weight > 0:
                    records.append((age, min(count, 3), weight))
    return np.asarray(records, float), pins


def fecundity_by_age(P, ages):
    f = np.ones(len(ages))
    if P.fecundity_omega1 != 0:
        f = np.clip(1 - P.fecundity_omega1 * np.exp(P.fecundity_omega2 * (ages - P.age_start)), 0, 1)
        f[ages >= P.fecundity_terminal_age] = 0
    if P.fecundity_terminal_decay:
        f *= np.exp(-P.fecundity_terminal_decay * np.maximum(ages - P.fecundity_tail_start_age, 0))
    return f


def model_arrays(result):
    s, P = result.solution, result.P
    if int(P.I) != 1:
        raise ValueError('Assessment gross-earnings grid currently requires exactly one location (P.I == 1)')
    if (P.child_state_mode != 'independent_count' or not P.sequential_births
            or P.joint_nested_choice or P.readiness_gate_enabled):
        raise ValueError('Saved fertility architecture differs from exact inversion contract')
    start = P.age_start + np.arange(P.J) * P.period_years
    post = np.asarray(s.g_beginning_distribution, float)
    pre, rate = invert_independent_birth_map(
        post, np.asarray(s.fert_probs), np.asarray(s.fert2_probs),
        fecundity_by_age(P, start), range(P.A_f_start - 1, P.A_f_end))
    flow = np.sum(pre * rate, axis=(0, 1, 2, 4, 6))[:, 0]
    native = np.asarray(P._first_births_by_age)
    if np.max(np.abs(flow - native)) > 1e-10:
        raise ValueError('Inverted first-birth flow disagrees with saved native flow')
    return post, pre, start, native


def model_stock_at_age(pre, post, starts, age, period):
    j = int(np.floor((age + .5 - starts[0]) / period))
    f = (age + .5 - starts[j]) / period
    return (1 - f) * pre[:, :, :, j] + f * post[:, :, :, j]


def profiles(result, empirical_ages, empirical_weights, pre, post, starts):
    out = []
    for age in sorted(set(int(x) for x in empirical_ages)):
        if age < starts[0] or age >= starts[-1] + result.P.period_years:
            continue
        q = model_stock_at_age(pre, post, starts, age, result.P.period_years)
        n = q.sum(axis=(0, 1, 2, 3, 5))
        out.append((age, n / n.sum()))
    return out


def prepare_assessment(result, run_directory, *, output=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from production.engine.shared import annual_gross_income_at_state

    case = Path(run_directory)
    out = Path(output) if output else case / 'aggregate_plots/model_data_assessment'
    out.mkdir(parents=True, exist_ok=True)
    psid_meta = json.loads(PSID.with_name('metadata.json').read_text())
    require_sha256(PSID, psid_meta['selected_cache_sha256'], 'PSID cache')
    cps_manifest = json.loads((ROOT / 'code/data/cps_fertility/source_manifest.json').read_text())
    require_sha256(CPS24, cps_manifest['sha256'], 'CPS 2024')
    baseline, baseline_pins = baseline_cps()
    post, pre, starts, first_flow = model_arrays(result)
    P, s = result.P, result.solution
    g = np.asarray(s.g, float)
    b = np.asarray(s.b_grid, float)
    ages = starts + P.period_years / 2
    long, figures, coverage = [], [], {}

    def add(panel, series, x, y):
        long.extend({'page': page, 'panel': panel, 'series': series, 'x': float(a), 'value': float(z)}
                    for a, z in zip(np.ravel(x), np.ravel(y)))

    def line(ax, panel, title, model_x, model_y, data_x, data_y, ylabel):
        ax.plot(model_x, model_y, 'o-', color=BLUE, lw=1.7, ms=3, label='Model')
        ax.plot(data_x, data_y, 'o-', color=ORANGE, lw=1.7, ms=3, label='Data')
        ax.set(title=title, xlabel='Age' if 'age' in panel else '', ylabel=ylabel)
        ax.grid(alpha=.2)
        add(panel, 'Model', model_x, model_y)
        add(panel, 'Data', data_x, data_y)

    def bars(ax, panel, title, labels, model_values, data_values, ylabel='Share'):
        x = np.arange(len(labels))
        ax.bar(x - .18, model_values, .36, color=BLUE, label='Model')
        ax.bar(x + .18, data_values, .36, color=ORANGE, label='Data')
        ax.set(title=title, xticks=x, xticklabels=labels, ylabel=ylabel)
        add(panel, 'Model', x, model_values)
        add(panel, 'Data', x, data_values)

    def cdf(ax, panel, title, model_v, model_w, data_v, data_w, xlabel):
        mv, mc = weighted_ecdf(model_v, model_w)
        dv, dc = weighted_ecdf(data_v, data_w)
        ax.step(mv, mc, where='post', color=BLUE, label='Model')
        ax.step(dv, dc, where='post', color=ORANGE, label='Data')
        ax.set(title=title, xlabel=xlabel, ylabel='CDF', ylim=(0, 1.02))
        ax.grid(alpha=.2)
        add(panel, 'Model', mv, mc)
        add(panel, 'Data', dv, dc)
        if 'resource' in panel or panel.startswith('rooms_'):
            lower_quantile = .01 if 'resource' in panel else 0
            upper_quantile = .975 if 'resource' in panel else .99
            low = min(np.interp(lower_quantile, mc, mv), np.interp(lower_quantile, dc, dv))
            high = max(np.interp(upper_quantile, mc, mv), np.interp(upper_quantile, dc, dv))
            if high > low:
                ax.set_xlim(low, high)
                coverage[panel] = {'model_offscreen_mass': float(np.interp(low, mv, mc) + 1 - np.interp(high, mv, mc)),
                                   'data_offscreen_mass': float(np.interp(low, dv, dc) + 1 - np.interp(high, dv, dc)),
                                   'visible_range': [float(low), float(high)],
                                   'display_quantiles': [lower_quantile, upper_quantile]}

    def newpage(title, footnote):
        fig, axs = plt.subplots(3, 2, figsize=(8.5, 11))
        fig.suptitle(title, fontsize=15, y=.987)
        fig.text(.045, .012, footnote, ha='left', va='bottom', fontsize=8, wrap=True)
        return fig, axs.ravel()

    def finish(fig, page):
        for ax in fig.axes:
            if ax.get_legend_handles_labels()[0]:
                ax.legend(frameon=False, fontsize=8)
            ax.tick_params(labelsize=8)
            ax.title.set_fontsize(10)
            ax.xaxis.label.set_size(9)
            ax.yaxis.label.set_size(9)
        fig.subplots_adjust(left=.11, right=.97, top=.94, bottom=.065, hspace=.42, wspace=.32)
        figures.append(fig)
        fig.savefig(out / f'page_{page}.png', dpi=180)

    # Fertility: same-cell pre/post stocks, weighted by the survey's exact age exposure.
    page = 1
    fig, axs = newpage('Fertility: model versus data',
        'CPS June 2004/06 women (orange); saved model reproductive-member proxy (blue). Counts cap at 3+. '
        'First births: NCHS 2003–06 counts. Income: CPS June 2024 family money income versus model current labor earnings ranks; external validation.')
    cx = np.arange(18, 46)
    data_mean, model_mean = [], []
    for age in cx:
        r = baseline[baseline[:, 0] == age]
        data_mean.append(np.average(r[:, 1], weights=r[:, 2]) if len(r) else np.nan)
        q = model_stock_at_age(pre, post, starts, age, P.period_years)
        nm = q.sum(axis=(0, 1, 2, 3, 5))
        model_mean.append(nm @ np.arange(4) / nm.sum())
    line(axs[0], 'children_by_age', 'Mean children ever born by age', cx, model_mean, cx, data_mean, 'Mean, capped at 3')
    for ax, lo, hi, panel in ((axs[1], 22, 25, 'children_22_25'), (axs[2], 40, 44, 'children_40_44')):
        r = baseline[(baseline[:, 0] >= lo) & (baseline[:, 0] <= hi)]
        d = np.bincount(r[:, 1].astype(int), weights=r[:, 2], minlength=4); d /= d.sum()
        m = np.zeros(4)
        for age in range(lo, hi + 1):
            w = r[r[:, 0] == age, 2].sum()
            q = model_stock_at_age(pre, post, starts, age, P.period_years)
            nm = q.sum(axis=(0, 1, 2, 3, 5)); m += w * nm / nm.sum()
        m /= m.sum()
        bars(ax, panel, f'Children ever born, ages {lo}–{hi}', ['0', '1', '2', '3+'], m, d)
    nr = [r for r in rows(NCHS) if 2003 <= int(r['year']) <= 2006]
    n_age = np.array([int(r['age']) for r in nr]); n_w = np.array([float(r['n_first_births']) for r in nr])
    n_bins = np.clip((n_age - 18) // 4, 0, 6)
    n_mass = np.bincount(n_bins, weights=n_w, minlength=7); n_mass /= n_mass.sum()
    m_mass = first_flow[:7] / first_flow[:7].sum()
    bars(axs[3], 'first_birth_age', 'First births by age band', ['20', '24', '28', '32', '36', '40', '44'], m_mass, n_mass)
    axs[3].set_xlabel('Four-year band midpoint')
    coverage['nchs_boundary_collapsed'] = {'below_18_share': float(n_w[n_age < 18].sum() / n_w.sum()),
                                          'above_45_share': float(n_w[n_age > 45].sum() / n_w.sum()),
                                          'operator': 'ages <=21 to 18–21; >=42 to 42–45'}
    # Read only fields needed from the verified 2024 source.
    income_rows = []
    with CPS24.open(newline='') as stream:
        for r in csv.DictReader(stream):
            if r['PESEX'] != '2':
                continue
            age, count, wt, category = map(number, (r['PRTAGE'], r['PTSF1'], r['PWSSWGT'], r['HEFAMINC']))
            if ((24 <= age <= 26 or 40 <= age <= 44) and 0 <= count <= 5
                    and wt > 0 and 1 <= category <= 16):
                income_rows.append((int(age), min(int(count), 3), wt, int(category)))
    income_rows = np.asarray(income_rows, float)
    earnings_grid = np.array([[annual_gross_income_at_state(P, 0, j, float(z)) if j < P.J_R else np.nan
                               for z in P.z_grid] for j in range(P.J)])
    for ax, lo, hi, panel in ((axs[4], 24, 26, 'income_24_26'), (axs[5], 40, 44, 'income_40_44')):
        rr = income_rows[(income_rows[:, 0] >= lo) & (income_rows[:, 0] <= hi)]
        category_mass = np.bincount(rr[:, 3].astype(int), weights=rr[:, 2], minlength=17)[1:]
        alloc = thirds(category_mass)
        dv = [float(np.sum(rr[:, 2] * rr[:, 1] * alloc[k, rr[:, 3].astype(int) - 1]) /
                    np.sum(rr[:, 2] * alloc[k, rr[:, 3].astype(int) - 1])) for k in range(3)]
        state = []
        for age in range(lo, hi + 1):
            w = rr[rr[:, 0] == age, 2].sum()
            q = model_stock_at_age(pre, post, starts, age, P.period_years)
            j = int((age + .5 - starts[0]) // P.period_years)
            yz = q.sum(axis=(0, 1, 2, 4, 5))
            for z in range(P.Nz):
                state.append((earnings_grid[j, z], w * yz[z] / q.sum(),
                              w * np.sum(q[:, :, :, z] * np.arange(4)[None, None, None, :, None]) / q.sum()))
        # Rank unique gross earnings values; split all exact ties fractionally.
        earnings = sorted(set(x[0] for x in state))
        masses = np.array([sum(x[1] for x in state if x[0] == e) for e in earnings])
        numerators = np.array([sum(x[2] for x in state if x[0] == e) for e in earnings])
        a = thirds(masses)
        mv = np.sum(a * numerators, axis=1) / np.sum(a * masses, axis=1)
        bars(ax, panel, f'Children by income third, ages {lo}–{hi}', ['Low', 'Middle', 'High'], mv, dv, 'Mean, cap 3')
    finish(fig, page)

    # Housing: realized tenure and renter policy; annual ACS cells are disjoint.
    page = 2
    fig, axs = newpage('Housing: model versus data',
        'Room CDFs: PSID 2005/07; display to 99th percentile. Means: ACS 2005/06 national all structures. '
        'Family rooms use ages 30–55 and resident children (ACS) versus children at home (model). Ownership target has a distinct sample.')
    pr = rows(PSID)
    ps = np.array([[number(r[k]) for k in ('age', 'weight', 'total_net_wealth', 'annual_gross_labor_earnings',
                                            'owner', 'gross_home_value', 'rooms')] for r in pr])
    pg = g[:, 0]
    renter_policy = np.asarray(s.hR_pol)[:, 0]
    rv, rw = renter_policy.ravel(), pg.ravel()
    ov = np.broadcast_to(np.asarray(P.H_own)[None, :, None, None, None, None, None], g[:, 1:].shape).ravel()
    ow = g[:, 1:].ravel()
    for ax, key, mask_model, mask_data, title in (
        (axs[0], 'rooms_all', None, np.isfinite(ps[:, 6]), 'Rooms: all households'),
        (axs[1], 'rooms_owner', True, (ps[:, 4] == 1) & np.isfinite(ps[:, 6]), 'Rooms: owners'),
        (axs[2], 'rooms_renter', False, (ps[:, 4] == 0) & np.isfinite(ps[:, 6]), 'Rooms: renters')):
        if mask_model is None:
            v, w = np.r_[rv, ov], np.r_[rw, ow]
        elif mask_model:
            v, w = ov, ow
        else:
            v, w = rv, rw
        cdf(ax, key, title, v, w, ps[mask_data, 6], ps[mask_data, 1], 'Rooms')
    acs = select_grouped_national_acs(rows(ACS))
    age_data = []
    for start in starts:
        group = [r for r in acs if start <= int(r['age_lower']) < start + P.period_years]
        mass = sum(float(r['hhwt']) for r in group)
        if mass:
            age_data.append((start + P.period_years / 2,
                             sum(float(r['owner_hhwt']) for r in group) / mass,
                             sum(float(r['rooms_capped9_sum']) for r in group) / mass))
    age_model = []
    for j, age in enumerate(ages):
        q = g[:, :, :, j]
        own = q[:, 1:].sum() / q.sum()
        rent_q = q[:, 0]
        rent_h = renter_policy[:, :, j]
        valid = rent_q > 0
        room_total = np.sum(rent_q[valid] * np.minimum(rent_h[valid], 9))
        for ten, room in enumerate(P.H_own, 1):
            room_total += q[:, ten].sum() * min(room, 9)
        age_model.append((age, own, room_total / q.sum()))
    line(axs[3], 'ownership_by_age', 'Homeownership by age', [x[0] for x in age_model],
         [x[1] for x in age_model], [x[0] for x in age_data], [x[1] for x in age_data], 'Owner share')
    line(axs[4], 'rooms_by_age', 'Mean rooms by age (cap 9)', [x[0] for x in age_model],
         [x[2] for x in age_model], [x[0] for x in age_data], [x[2] for x in age_data], 'Rooms')
    children_data, children_model = [], []
    for m, label in enumerate(('0', '1', '2', '3+')):
        group = [r for r in acs if r['current_children'] == label and 30 <= int(r['age_lower']) <= 55]
        mass = sum(float(r['hhwt']) for r in group)
        children_data.append(sum(float(r['rooms_capped9_sum']) for r in group) / mass)
        age_fraction = np.array([max(0, min(start + P.period_years - 1, 55) - max(start, 30) + 1)
                                 / P.period_years for start in starts])
        q = g[..., m] * age_fraction[None, None, None, :, None, None]
        renter_rooms = np.asarray(s.hR_pol)[..., m][:, 0]
        room_total = np.sum(q[:, 0][q[:, 0] > 0] * np.minimum(renter_rooms[q[:, 0] > 0], 9))
        room_total += sum(q[:, ten].sum() * min(room, 9) for ten, room in enumerate(P.H_own, 1))
        children_model.append(room_total / q.sum())
    bars(axs[5], 'rooms_by_children', 'Mean rooms by children at home, ages 30–55', ['0', '1', '2', '3+'],
         children_model, children_data, 'Rooms, cap 9')
    finish(fig, page)

    # Resources: beginning net financial position includes mortgage debt;
    # total net wealth adds the value of the home held at the same instant.
    page = 3
    fig, axs = newpage('Resources: model versus data',
        'PSID 2005/07 reference persons. Each domain divides earnings and wealth by its own weighted mean gross working-age earnings. '
        'Model wealth: beginning net financial position plus old home value. CDF tails retained; display 1st–97.5th percentile range.')
    age_p, wt_p, total_p, earn_p, owner_p, home_p, room_p = ps.T
    valid_work = (age_p <= 65) & np.isfinite(earn_p) & (earn_p >= 0)
    valid_wealth = np.isfinite(total_p) & ((age_p > 65) | valid_work)
    valid_finance = valid_wealth & np.isfinite(home_p)
    normalizer_sample = valid_work & valid_wealth
    denom_p = np.average(earn_p[normalizer_sample], weights=wt_p[normalizer_sample])
    finance_p, _ = financial_positions(total_p, home_p)
    begin = np.asarray(s.g_beginning_distribution, float)
    working_mass = 0.; working_income = 0.
    for j in range(P.J_R):
        mass_z = begin[:, :, :, j].sum(axis=(0, 1, 2, 4, 5))
        working_mass += mass_z.sum()
        working_income += sum(mass_z[z] * earnings_grid[j, z] for z in range(P.Nz))
    denom_m = working_income / working_mass
    ev = earnings_grid[:P.J_R].ravel() / denom_m
    ew = begin[:, :, :, :P.J_R].sum(axis=(0, 1, 2, 5, 6)).ravel()
    bv = np.broadcast_to(b[:, None, None, None, None, None, None], begin.shape)
    p = np.asarray(s.p_eq)
    old_rooms = np.r_[0, P.H_own]
    home_value = p[None, None, :, None, None, None, None] * old_rooms[None, :, None, None, None, None, None]
    fin_m = bv / denom_m
    total_m = (bv + home_value) / denom_m
    cdf(axs[0], 'resource_earnings_cdf', 'Gross labor earnings / own mean', ev, ew,
        earn_p[valid_work] / denom_p, wt_p[valid_work], 'Multiple of mean earnings')
    cdf(axs[1], 'resource_financial_cdf', 'Net financial position / own mean',
        fin_m.ravel(), begin.ravel(), finance_p[valid_finance] / denom_p, wt_p[valid_finance], 'Multiple of mean earnings')
    cdf(axs[2], 'resource_totalwealth_cdf', 'Total net wealth / own mean',
        total_m.ravel(), begin.ravel(), total_p[valid_wealth] / denom_p, wt_p[valid_wealth], 'Multiple of mean earnings')
    model_earn, model_wealth, model_neg = [], [], []
    data_earn, data_wealth, data_neg = [], [], []
    for j, start in enumerate(starts):
        q = begin[:, :, :, j]
        model_earn.append((age_model[j][0],
                           sum(q[:, :, :, z].sum() * earnings_grid[j, z] for z in range(P.Nz)) / q.sum() / denom_m)
                          if j < P.J_R else None)
        model_wealth.append((age_model[j][0], np.sum(q * (bv[:, :, :, j] + home_value[:, :, :, 0])) / q.sum() / denom_m))
        model_neg.append((age_model[j][0], q[b < 0].sum() / q.sum()))
        keep = (age_p >= start) & (age_p < start + P.period_years)
        k = keep & valid_work
        if np.any(k):
            data_earn.append((age_model[j][0], np.average(earn_p[k], weights=wt_p[k]) / denom_p))
        k = keep & valid_wealth
        if np.any(k):
            data_wealth.append((age_model[j][0], np.average(total_p[k], weights=wt_p[k]) / denom_p))
        k = keep & valid_finance
        if np.any(k):
            data_neg.append((age_model[j][0], np.average(finance_p[k] < 0, weights=wt_p[k])))
    model_earn = [x for x in model_earn if x is not None]
    for ax, panel, title, model_rows, data_rows, ylabel in (
        (axs[3], 'resource_earnings_age', 'Mean gross earnings by age', model_earn, data_earn, 'Multiple of own mean'),
        (axs[4], 'resource_totalwealth_age', 'Mean total wealth by age', model_wealth, data_wealth, 'Multiple of own mean'),
        (axs[5], 'resource_negative_age', 'Negative financial position by age', model_neg, data_neg, 'Share')):
        line(ax, panel, title, [x[0] for x in model_rows], [x[1] for x in model_rows],
             [x[0] for x in data_rows], [x[1] for x in data_rows], ylabel)
    finish(fig, page)

    pdf = out / 'model_data_assessment.pdf'
    with PdfPages(pdf) as book:
        for fig in figures:
            book.savefig(fig)
            plt.close(fig)
    with (out / 'plotted_long.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['page', 'panel', 'series', 'x', 'value'])
        writer.writeheader(); writer.writerows(long)
    receipt = {
        'status': 'complete', 'case': str(case.resolve()),
        'case_files_sha256': {name: sha256(case / name) for name in
                              ('native_result.npz', 'metadata.json', 'target_fit.csv', 'parameters.csv')},
        'sources_sha256': {'psid_cache': sha256(PSID), 'psid_metadata': sha256(PSID.with_name('metadata.json')),
                           'cps_2024': sha256(CPS24), 'cps_2004_2006_partitions': baseline_pins,
                           'acs_profile': sha256(ACS), 'nchs_first_birth_counts': sha256(NCHS)},
        'normalizers': {'psid_annual_gross_labor_age18_65': float(denom_p),
                        'model_annual_gross_labor_before_retirement': float(denom_m)},
        'coverage': {'psid_rows': len(ps), 'psid_earnings_rows': int(valid_work.sum()),
                     'psid_normalizer_rows': int(normalizer_sample.sum()),
                     'psid_wealth_rows': int(valid_wealth.sum()),
                     'psid_financial_rows': int(valid_finance.sum()), **coverage},
        'definitions': {'fertility': 'prebirth/postbirth stock interpolation within same four-year age cell, empirical integer-age exposure',
                        'income': 'fractional weighted thirds of HEFAMINC (CPS) and current gross labor income (model)',
                        'model_financial': 'beginning b including mortgage debt; total b + price times old owned rooms',
                        'psid_financial': 'NETWORTHR minus gross HOMEVALUER',
                        'housing': 'realized g, renter hR_pol, owner H_own; ACS annual disjoint cells; room means cap 9; family rooms ages 30–55, with 54–57 model cell weighted 2/4',
                        'nchs': '2003–06 first-birth counts, boundary-collapsed four-year age bands; counts not hazards'},
        'sources_limitations': ['CPS 2024 income comparison is external validation, not matched-vintage calibration.',
                                'ACS ownership target DUE sample differs from national all-structures profiles.',
                                'AHS 2007 has a room mean but no local room distribution; PSID provides CDF validation.',
                                'PSID reference-person earnings are household earnings; model is reproductive-member proxy.'],
        'files': {'pdf': str(pdf), 'pages': [str(out / f'page_{i}.png') for i in (1, 2, 3)],
                  'plotted_long': str(out / 'plotted_long.csv')}
    }
    (out / 'metadata.json').write_text(json.dumps(receipt, indent=2, allow_nan=False) + '\n')
    return out
