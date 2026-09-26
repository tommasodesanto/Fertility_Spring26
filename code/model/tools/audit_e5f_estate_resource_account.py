#!/usr/bin/env python3
"""Read-only estate/entry stock-flow audit of the retained three-case experiment.

All imports and execution belong on Torch. This does not solve households,
change a receipt law, or certify an aggregate national resource constraint.
"""
from __future__ import annotations

import argparse
import gc
import gzip
import json
import math
from pathlib import Path
import pickle


def signed_accounts(g, bp, death, housing_values, selling_cost):
    """Integrate signed estate values without hiding negative estates."""
    import numpy as np
    if g.shape != bp.shape or g.ndim != 7 or g.shape[2] != 1:
        raise ValueError("Expected the retained one-location seven-axis arrays")
    if len(death) != g.shape[3] or len(housing_values) != g.shape[1]:
        raise ValueError("Age/tenure supports differ")
    if not 0 <= selling_cost < 1:
        raise ValueError("Invalid selling cost")
    if (not np.isfinite(g).all() or not np.isfinite(bp).all()
            or np.any(g < 0) or not np.isfinite(death).all()
            or np.any(np.asarray(death) < 0) or np.any(np.asarray(death) > 1)
            or not np.isfinite(housing_values).all()
            or np.any(np.asarray(housing_values) < 0)):
        raise ValueError("Invalid distribution, policy, death or house values")
    rows = []
    for j, d in enumerate(death):
        row = dict(age_index=j, death_mass=0., gross_positive=0.,
                   gross_negative=0., net_positive=0., net_negative=0.,
                   negative_estate_death_mass=0., housing_sale_cost=0.)
        for ten, house in enumerate(housing_values):
            mass = g[:, ten, 0, j] * float(d)
            gross = bp[:, ten, 0, j] + float(house)
            net = gross - selling_cost * float(house)
            row["death_mass"] += float(mass.sum())
            for label, values in (("gross", gross), ("net", net)):
                row[label + "_positive"] += float((mass * np.maximum(values, 0)).sum())
                row[label + "_negative"] += float((mass * np.maximum(-values, 0)).sum())
            row["negative_estate_death_mass"] += float(mass[net < 0].sum())
            row["housing_sale_cost"] += float(mass.sum()) * selling_cost * float(house)
        rows.append(row)
    totals = {key: sum(row[key] for row in rows) for key in rows[0] if key != "age_index"}
    totals["gross_signed"] = totals["gross_positive"] - totals["gross_negative"]
    totals["net_signed"] = totals["net_positive"] - totals["net_negative"]
    totals["signed_cost_identity_residual"] = (
        totals["gross_signed"] - totals["net_signed"] - totals["housing_sale_cost"])
    if abs(totals["signed_cost_identity_residual"]) > 1e-10:
        raise AssertionError("Signed liquidation-cost identity failed")
    return dict(totals=totals, by_age=rows)


def fixtures():
    import numpy as np
    g = np.zeros((2, 2, 1, 2, 1, 1, 1))
    bp = np.zeros_like(g)
    # A renter dies owing 2; an owner leaves gross 1, net -1 after costs.
    g[0, 0, 0, 1, 0, 0, 0] = .25
    bp[0, 0, 0, 1, 0, 0, 0] = -2
    g[0, 1, 0, 1, 0, 0, 0] = .5
    bp[0, 1, 0, 1, 0, 0, 0] = -9
    # Another owner leaves gross 15, net 13.
    g[1, 1, 0, 1, 0, 0, 0] = .25
    bp[1, 1, 0, 1, 0, 0, 0] = 5
    # This large negative-wealth state survives and must not enter deaths.
    g[0, 0, 0, 0, 0, 0, 0] = 1
    bp[0, 0, 0, 0, 0, 0, 0] = -100
    a = signed_accounts(g, bp, np.array([0., 1.]), [0., 10.], .2)["totals"]
    expected = dict(death_mass=1., gross_positive=4.25, gross_negative=.5,
                    net_positive=3.25, net_negative=1.,
                    negative_estate_death_mass=.75, housing_sale_cost=1.5,
                    gross_signed=3.75, net_signed=2.25)
    for key, value in expected.items():
        if not math.isclose(a[key], value, abs_tol=1e-14):
            raise AssertionError((key, a[key], value))
    return dict(status="passed", checked_scalars=len(expected),
                cases="renter debt; owner crossing zero after costs; positive owner; surviving debt")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    import run_e5f_estate_receiver_probe as probe
    probe.require_torch()
    destination = args.output_dir / "resource_account.json"
    if destination.exists():
        raise FileExistsError("Refusing to replace retained resource audit")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    test = fixtures()
    # Recreate the hash-pinned import/runtime context required by the pickles.
    context = probe.verified_runtime(args.output_dir / "runtime")
    del context
    gc.collect()
    import numpy as np
    cases = []
    for name in probe.CASES:
        case_dir = args.run_root / name
        with gzip.open(case_dir / "initial_state.pkl.gz", "rb") as stream:
            saved = pickle.load(stream)
        P, sol = saved["parameters"], saved["solution"]
        retained = json.loads((case_dir / "receipt.json").read_text())
        estate = retained["estate_accounts"]
        price = float(np.asarray(sol.p_eq).reshape(-1)[0])
        survival = (np.asarray(P.survival_probs, dtype=float)
                    if P.use_age_survival else np.ones(int(P.J)-1))
        death = np.r_[1-survival, 1.]
        audit = signed_accounts(np.asarray(sol.g), np.asarray(sol.bp_pol), death,
                                np.r_[0., price*np.asarray(P.H_own)], float(P.psi))
        for key, observed in (("net_positive", estate["generated_net_period"]),
                              ("gross_positive", estate["generated_gross_period"])):
            if abs(audit["totals"][key] - observed) > 1e-10:
                raise AssertionError("Retained estate observer differs: " + key)
        pre = np.asarray(saved["stationary_g_pre"])
        entrant = pre[:, :, :, 0]
        if float(entrant[:, 1:].sum()) > 1e-10:
            raise AssertionError("Entrant pre-choice measure unexpectedly contains owners")
        b = np.asarray(saved["b_grid"]).reshape((-1, 1, 1, 1, 1, 1))
        entry_mass = float(entrant.sum())
        if abs(entry_mass-retained["normalized_entrant_flow"]) > 1e-10:
            raise AssertionError("Entrant flow differs from retained receipt")
        entry_positive = float((entrant * np.maximum(b, 0)).sum())
        entry_negative = float((entrant * np.maximum(-b, 0)).sum())
        entry = dict(mass=entry_mass, positive_financial_endowment=entry_positive,
                     negative_financial_position=entry_negative,
                     signed_financial_endowment=entry_positive-entry_negative,
                     mean_signed_wealth=(entry_positive-entry_negative)/entry_mass)
        paid = float(estate["paid_period"])
        cases.append(dict(case=name, estate=audit, entry=entry, paid=paid,
                          positive_net_estates_not_paid=audit["totals"]["net_positive"]-paid,
                          checkpoint_sha256=probe.sha(case_dir / "initial_state.pkl.gz")))
        del saved, P, sol, pre, entrant
        gc.collect()
    result = dict(status="retained_stock_flow_audit_completed", fixtures=test, cases=cases,
                  interpretation="All values are per model four-year period. Negative net estates identify unpaid liabilities at exit, not an adopted creditor-loss rule. Entrant financial positions are a separate external boundary stock. This is not an aggregate national resource-constraint certification.")
    probe.write(destination, result)
    print(json.dumps({"status":result["status"],"cases":[{"case":c["case"],
          "negative_net_estates":c["estate"]["totals"]["net_negative"],
          "mean_entry_wealth":c["entry"]["mean_signed_wealth"]} for c in cases]}))


if __name__ == "__main__":
    main()
