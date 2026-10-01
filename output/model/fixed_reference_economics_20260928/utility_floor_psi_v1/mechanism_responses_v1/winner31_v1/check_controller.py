"""Zero-lifecycle tests of the six-cell order and fail-closed clocks."""
import json, tempfile, time
from pathlib import Path
import fixed_price_responses as driver

def fake(_auth, regime, factor, out, _deadline):
    driver.write(out/'mock.json', {'regime': regime, 'factor': factor})
    return {'lifecycle_seconds': 0., 'support_diagnostic': {'status': 'mock_no_certificate'}}

with tempfile.TemporaryDirectory(prefix='winner31_mock_') as folder:
    base=Path(folder)
    done=driver.run(base/'six',time.time()+1200,auth={},evaluate=fake,mock=True)
    assert [(x['regime'],x['price_factor']) for x in done]==[
        ('reference',1.),('reference',.99),('reference',1.01),
        ('lifetime_repayment_only',.99),('lifetime_repayment_only',1.),('lifetime_repayment_only',1.01)]
    launch=json.loads((base/'six/launch.json').read_text())
    assert launch['lifecycle_solve_cap']==6 and launch['case_budget_seconds']==300
    assert launch['global_budget_seconds']==1200 and launch['memory_gib_cap']==24
    assert launch['deadline_epoch']-launch['started_epoch']<=1200.0001
    for row in done:
        started=json.loads((base/'six'/row['label']/'latest.json').read_text())
        assert 0<started['case_deadline_epoch']-started['started_epoch']<=300.0001
    def slow(_auth,_regime,_factor,_out,_deadline):
        time.sleep(.2)
        return fake(_auth,_regime,_factor,_out,_deadline)
    try:driver.run(base/'deadline',time.time()+.03,auth={},evaluate=slow,mock=True)
    except driver.CaseDeadline:pass
    else:raise AssertionError('First case did not honor deadline')
    failed=json.loads((base/'deadline/completed.json').read_text())
    assert failed['status']=='fatal_q0_reference_failure' and failed['case_attempts']==1 and failed['no_auto_retry']
print('PASS exact six-cell zero-LC loop, 300/1200/24 caps, q0 timeout and no retry')
