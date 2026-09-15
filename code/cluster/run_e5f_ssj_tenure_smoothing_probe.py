"""Experiment-only probe: evaluate one fixed 104-date path under a larger tenure-smoothing scale.

Reuses the announced batch's context, endpoint, preference routing and the
frozen native path operator.  The only change is ``tenure_choice_kappa`` on a
deep copy of the initial-state parameters (the logit scale of the own/rent
mixture).  Prices, pensions, rebates, preference path, initial 2007 state,
queue law, gates and residual definitions are the warm-start checkpoint's.
No root is solved; one mapping per job.  The terminal continuation remains
the kappa=0.005 endpoint, which is recorded as a caveat.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from unittest.mock import patch

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=lambda v: v.tolist()) + "\n")
    tmp.replace(path)


def summarize_rows(rows, residual, horizon, jump_threshold=0.004):
    """Ownership-jump and residual summary for one evaluated path (pure numpy)."""
    own = np.asarray([r["owner_rate"] for r in rows], dtype=float)
    births = np.asarray([r["birth_children"] for r in rows], dtype=float)
    mass = np.asarray([r["adult_population"] for r in rows], dtype=float)
    revenue = np.asarray([r["property_tax_revenue"] for r in rows], dtype=float)
    jumps = np.diff(own)
    residual = np.asarray(residual, dtype=float)
    h, p, reb = residual[:horizon], residual[horizon:2 * horizon], residual[2 * horizon:]
    return dict(
        owner_rate=own.tolist(), births=births.tolist(), adult_population=mass.tolist(), property_tax_revenue=revenue.tolist(),
        owner_rate_first=float(own[0]), owner_rate_mean=float(own.mean()), owner_rate_last=float(own[-1]),
        max_abs_owner_jump_pp=float(100 * np.max(np.abs(jumps))), median_abs_owner_jump_pp=float(100 * np.median(np.abs(jumps))),
        owner_jump_dates_above_threshold=[(int(t + 1), float(100 * jumps[t])) for t in np.where(np.abs(jumps) > jump_threshold)[0]],
        births_total=float(births.sum()),
        max_abs_housing=float(np.max(np.abs(h))), max_abs_paygo=float(np.max(np.abs(p))), max_abs_rebate=float(np.max(np.abs(reb))),
        rebate_argmax_date=int(np.argmax(np.abs(reb))),
        rebate_residual_at=dict((str(t), float(reb[t])) for t in (43, 44, 68, 69, 70, 71, 99, 100, 101, 102)),
        rebate_spikes_above_0p15=[(int(t), float(reb[t])) for t in np.where(np.abs(reb) > 0.15)[0]],
        residual=residual.tolist())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()
    m = read(args.manifest)
    for path, digest in m["file_sha256"].items():
        if sha(path) != digest:
            raise ValueError("Changed pinned input: " + path)
    announced = read(m["announced_manifest"])
    checkpoint = read(m["warm_start_checkpoint"])
    kappa = float(m["tenure_choice_kappa"])
    out = Path(m["output"])
    if (out / "probe_contract.json").exists():
        raise ValueError("Refusing to overwrite a started probe")
    sys.path.insert(0, str(Path(m["announced_source_dir"])))
    import run_e5f_announced_original_queue as ann
    spec = read(announced["spec"])
    sys.path.insert(0, str(Path(spec["batch"]) / "source"))
    import run_e5f_original_queue_experiments as runner
    c = runner.load_context(announced["spec"])
    c.spec_path = Path(announced["spec"])
    endpoint, receipt = ann.endpoint_from_manifest(announced)
    horizon = 104
    x = np.asarray(checkpoint["prices"], dtype=float)
    if x.shape != (3 * horizon,):
        raise ValueError("checkpoint has the wrong horizon")
    prices, pensions, transfers = x[:horizon], x[horizon:2 * horizon], x[2 * horizon:]
    baseline_kappa = float(getattr(c.old.parameters, "tenure_choice_kappa", 0.0))
    started = time.time()
    save(out / "probe_contract.json", dict(manifest_sha256=sha(args.manifest), tenure_choice_kappa=kappa,
         baseline_tenure_choice_kappa=baseline_kappa, horizon=horizon, checkpoint_score=checkpoint.get("score"),
         terminal_continuation="kappa-0.005 endpoint retained (caveat: last dates mildly inconsistent at other scales)",
         policy_cache_used=True, root_solved=False, production_eligible=False, started_unix=started))
    probe_old = copy.copy(c.old)
    probe_old.parameters = copy.deepcopy(c.old.parameters)
    probe_old.parameters.tenure_choice_kappa = kappa
    psi_path = np.r_[announced["psi_levels"], np.full(100, ann.FINAL_PSI)]
    evaluate, _first = ann.announced_queue_path(c, psi_path)
    from types import SimpleNamespace as NS
    terminal = NS(parameters=endpoint.parameters, policy=endpoint.policy, asset_price=endpoint.asset_price)
    try:
        with c.queue.original_queue_adapter(), c.cache.policy_cache(c.joined.pf, max_bytes=12 * 1024**3):
            began = time.monotonic()
            result = evaluate(inherited=c.rebated.InheritedState(2007, c.old.initial_state), old_state=probe_old,
                              prices=prices, pensions=pensions, transfers=transfers, psi=float(psi_path[-1]),
                              terminal=terminal)
            seconds = time.monotonic() - began
        rows = result.rows
        if len(rows) != horizon:
            raise ValueError("Path did not return 104 rows")
        blocks = [c.rebated.dated_residual(demand=r["housing_demand"], supply=r["housing_supply"],
                                           payroll_accounts=r, tax_accounts=r) for r in rows]
        residual = c.rebated.stack_dated_residuals(blocks)
        summary = summarize_rows(rows, residual, horizon)
        summary.update(tenure_choice_kappa=kappa, mapping_seconds=seconds,
                       gates=dict(maximum_mass_accounting_error=float(result.maximum_mass_accounting_error),
                                  maximum_policy_reproduction_error=float(result.maximum_policy_reproduction_error),
                                  maximum_feasibility_projection_mass=float(result.maximum_feasibility_projection_mass)),
                       checkpoint_residual_max_abs=float(np.max(np.abs(np.asarray(checkpoint["residual"], dtype=float)))),
                       reproduces_checkpoint_residual=(float(np.max(np.abs(residual - np.asarray(checkpoint["residual"], dtype=float))))
                                                       if kappa == baseline_kappa else None))
        save(out / "rows.json", rows)
        save(out / "summary.json", summary)
        save(out / "probe_complete.json", dict(completed=True, tenure_choice_kappa=kappa, elapsed_seconds=time.time() - started,
             max_abs_owner_jump_pp=summary["max_abs_owner_jump_pp"], max_abs_rebate=summary["max_abs_rebate"]))
    except BaseException as exc:
        save(out / "probe_failure.json", dict(error_type=type(exc).__name__, error=str(exc), elapsed_seconds=time.time() - started))
        raise


if __name__ == "__main__":
    main()
