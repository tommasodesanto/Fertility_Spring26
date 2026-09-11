"""Build the return-home assessment from saved receipts, without model solves."""
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEST = HERE/'return_home_20260911'
BASELINE = HERE.parent/'initial_calibration_contract/extended_refinement/collected_17378993'


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def fmt(value):
    if value in ('', None):
        return '—'
    try:
        return f'{float(value):.7g}'
    except (TypeError, ValueError):
        return str(value).replace('|', '/')


def main():
    long_dir = DEST/'paths/delta_m005_long'
    long_summary = json.loads((long_dir/'summary.json').read_text())
    long_root = json.loads((long_dir/'root_receipt.json').read_text())
    assert not long_summary['finite_horizon_market_fiscal_converged']
    long_market = max(map(abs,long_root['final']['market_residual']))
    long_fiscal = max(map(abs,long_root['final']['fiscal_residual']))
    lines = ['# Quantitative assessment after the commute — September 11', '',
        'The utility and pension implementation has produced verified initial and terminal solutions and three cleared short historical paths. Historical preference estimation, continuation-horizon verification, and policies under the revised calibration remain unfinished. No further search was submitted during this assessment.', '',
        '## Main conclusions', '',
        '- Social Security: the actual-budget repair is verified in the initial economy, all three trial terminal equilibria, and all three short historical paths. This does not certify the still-solving longer paths.',
        '- Initial calibration: the unrestricted reference loss is 158.0764. Fixing annual beta at0.98 and reoptimizing gives222.5244, with two exact repetitions. The0.99 profile reached168.2011 before a population-mass check stopped it; it has no final exact repetitions or joint refinement and is not an optimum.',
        '- Lower beta has a real trade-off. At0.98 the first-birth rooms response is0.7575 against0.7202, but wealth/earnings falls to3.8823 against6.1459 and ownership at30–55 to49.60% against64.83%. All rows below must be considered together.',
        '- Preference sensitivity: all three short roots clear and their reported rate profiles reproduce exactly. These are fixed amplitudes, not adaptive estimation. Correct female/maternal-age measurement and the outer preference update remain missing.',
        f'- Longer paths: the central28-date run finished its eight-mapping budget without convergence (maximum housing residual{100*long_market:.4f}%, pension residual{100*long_fiscal:.6f}%). Exact replay/checkpoint checks pass. Both alternatives were still solving at collection. None supplies horizon certification.',
        '- Presentation: the revised specification and empirical evidence can be developed now. New numerical historical-fit and policy claims are not ready. Do not combine old policy results with the revised utility/pension calibration.', '',
        '## Reproduced short-path fertility diagnostics', '',
        '**These columns do not measure exactly the same object.** Data are equal-weight averages of published annual female period TFR. Model columns sum four-year birth flows divided by model-household mass in each age cell. Their numerical proximity is suggestive, not a certified data fit. The paths also use only six dates, so horizon sensitivity remains open.', '',
        '| Birth-data window | Model decision | Female TFR, data | Smaller decline −0.025 | Middle decline −0.05 | Larger decline −0.10 |',
        '|---|---:|---:|---:|---:|---:|']
    proof = {'files': {}, 'short_paths': {}, 'calibrations': {}}
    rates = {}
    for label in ['delta_m0025','delta_m005','delta_m010']:
        d = DEST/'paths'/label
        summary = json.loads((d/'summary.json').read_text())
        root = json.loads((d/'root_receipt.json').read_text())
        recorded = json.loads((d/'dated_household_fertility.json').read_text())['mappings']
        assert all(summary[k] for k in ['finite_horizon_market_fiscal_converged','mapping_replay_verified','checkpoint_reload_verified'])
        chosen, fresh = root['best']['payload']['trial'], root['final']['payload']['trial']
        assert recorded[chosen-1] == recorded[fresh-1]
        assert len(recorded[fresh-1]) == 6
        rates[label] = {r['calendar_year']:r['diagnostics']['period_tfr_topcode_adjusted'] for r in recorded[fresh-1]}
        proof['short_paths'][label] = {'selected_trial':chosen,'fresh_trial':fresh,'observer_exactly_reproduced':True}
        for file in d.iterdir():
            if file.is_file():proof['files'][str(file.relative_to(DEST))] = hashlib.sha256(file.read_bytes()).hexdigest()
    for row in read_csv(HERE/'inputs/empirical_blocks.csv'):
        year = int(row['decision_year'])
        cells = [f"{row['birth_year_start']}–{row['birth_year_end']}",year,row['period_tfr_arithmetic_mean'],*[rates[k][year] for k in ['delta_m0025','delta_m005','delta_m010']]]
        lines.append('| '+' | '.join(fmt(x) for x in cells)+' |')
    lines += ['', 'The2023 decision generates2024–2027 births; it is not compared with the2020–2023 data block. Different2007 outcomes reflect anticipation of different announced future paths despite the same inherited state.', '',
        '## Complete initial-calibration fits and parameters', '',
        'All three profiles use the same12 scored moments/weights and separate2.1 normalization. These are the approved working minimum-distance diagnostics, not a claim of fully certified empirical SMM or global identification. Baseline has nine free structural coordinates; fixed-beta profiles have eight. Near-bound flags reproduce the saved convention and need not indicate a binding constraint.', '']
    baseline_weights = None
    for label, directory, name in [
        ('reference', BASELINE, 'Unrestricted reference: reproduced r5_joint_09'),
        ('beta_098', DEST/'beta_098', 'Annual beta0.98: completed bounded search and two exact repetitions'),
        ('beta_099', DEST/'beta_099', 'Annual beta0.99: stopped first-round profile; unrepeated best candidate')]:
        rows = read_csv(directory/'selected_target_fit.csv')
        scored = [r for r in rows if r['scored'].lower()=='true']
        assert len(rows)==13 and len(scored)==12
        contract = [(r['restriction_id'],r['target'],r['actual_weight']) for r in rows]
        if baseline_weights is None:baseline_weights=contract
        else:assert contract==baseline_weights
        total = 0.0
        for row in scored:
            target, model, gap, weight, contribution = [float(row[k]) for k in ['target','model','gap','actual_weight','loss_contribution']]
            assert math.isclose(model-target,gap,abs_tol=1e-12)
            assert math.isclose(gap*gap*weight,contribution,abs_tol=1e-8)
            total+=contribution
        if label!='reference':
            summary=json.loads((directory/'summary.json').read_text())
            assert math.isclose(total,summary['best_loss'],abs_tol=1e-8)
        lines += [f'### {name}', '', f'Loss: **{total:.9f}**.', '', '| Moment | Target | Model | Gap | Weight | Loss contribution |', '|---|---:|---:|---:|---:|---:|']
        for row in rows:
            lines.append('| '+' | '.join([row['label'],*[fmt(row[k]) for k in ['target','model','gap','actual_weight','loss_contribution']]])+' |')
        lines += ['', '| Parameter | Value | Lower | Upper | Near bound | Restriction/status |', '|---|---:|---:|---:|---|---|']
        parameters = read_csv(directory/'selected_parameters.csv')
        assert len(parameters)==17
        for row in parameters:
            fixed=label!='reference' and row['parameter']=='beta_annual'
            lower,upper=(row['estimate'],row['estimate']) if fixed else (row['lower'],row['upper'])
            status='Externally fixed in this profile' if fixed else row['status']
            lines.append('| '+' | '.join([row['parameter'],fmt(row['estimate']),fmt(lower),fmt(upper),row['near_bound'],status])+' |')
        lines += ['',f'Source tables: `{directory}/selected_target_fit.csv` and `{directory}/selected_parameters.csv`.', '']
        proof['calibrations'][label]={'rows':len(rows),'scored':len(scored),'parameters':len(parameters),'recomputed_loss':total,'same_targets_and_weights':True}
        for name in ['selected_target_fit.csv','selected_parameters.csv']:
            path=directory/name;proof['files'][str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    lines += ['## Next decision', '',
        'Keep the announced-path experiment and finish its measurement objective and adaptive preference fit. Use the completed trial results as initial evaluations. Decide between the high-beta reference and a completed0.99 profile using the full fit, rather than demanding another initial-calibration optimum before historical fitting. Establish the horizon accuracy of the selected history before promoting policy results.', '',
        'The slides task can independently revise the exposition, utility-floor equation and pension description, while scientific claims and numerical replacements are obtained from the technical tasks.']
    (DEST/'READOUT.md').write_text('\n'.join(lines)+'\n')
    (DEST/'verification.json').write_text(json.dumps(proof,indent=2)+'\n')
    print('Verified three short-path rate replays, 39 target rows, 51 parameter rows, equal targets/weights and all loss contributions.')


if __name__=='__main__':
    main()
