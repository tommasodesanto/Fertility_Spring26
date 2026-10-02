"""Zero-solve native helper smoke under the staged Apptainer source overlays."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from financial_access import matched_first_birth_access

PACKET = Path(__file__).resolve().parent.parent
ROOT = PACKET.parents[3]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=("hard", "quarter"), required=True)
    args = ap.parse_args()
    engine = PACKET / "engines" / args.arm
    sys.path.insert(0, str(engine))
    from refactor_lab import inputs
    from refactor_lab.engine import solver
    if not Path(solver.__file__).resolve().is_relative_to(engine):
        raise RuntimeError("Native smoke imported the wrong isolated engine")
    bundle = ROOT / "output/model/publication_refactor_20260929/local_export_v1/inputs"
    loaded = inputs.load_inputs(bundle, ROOT, inputs.sha256_file(bundle / "bundle.json"))
    P = loaded.parameters
    P.experimental_purchase_saving_fraction = .25 if args.arm == "quarter" else 1.
    grid = loaded.b_grid
    SD = solver.precompute_shared(P, grid)
    shape = (len(grid), 1 + P.n_house, P.I, P.J, len(P.z_grid), P.n_parity, P.n_child_states)
    pre = np.zeros(shape)
    fert = np.zeros(shape[:-2] + (P.n_parity,))
    i = int(np.argmin(np.abs(grid - .1)))
    pre[i, 0, 0, 2, 4, 0, 0] = 1.
    fert[i, 0, 0, 2, 4, 1] = .5
    result = matched_first_birth_access(P, SD, grid, loaded.reference_price,
                                        pre, fert, rule=args.arm)
    if not 0 < result["first_birth_origin_renter_flow_mass"] < .5:
        raise RuntimeError("Fecundity-adjusted first-birth flow is invalid")
    if result["share_ineligible_at_80_eligible_at_100"] is None:
        raise RuntimeError("Matched financed-share access was not evaluated")
    print(json.dumps({"status": "zero_solve_native_access_passed", "arm": args.arm,
                      "engine": str(Path(solver.__file__).resolve()),
                      "first_birth_mass": result["first_birth_origin_renter_flow_mass"],
                      "access_gain_share": result["share_ineligible_at_80_eligible_at_100"]}))


if __name__ == "__main__":
    main()
