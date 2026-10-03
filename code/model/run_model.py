"""Editable local stationary-GE runner.  Importing this file is inert."""
from __future__ import annotations

import argparse

# Change this one name to select another data file in code/model/parameters.
PARAMETER_FILE = "best_params.py"


def main(argv=None):
    # Make the executable safe to inspect: argparse handles --help and rejects
    # unsupported arguments before importing the solver/workflow or reserving a case.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params", default=PARAMETER_FILE, help="Parameter file name or absolute path")
    args = parser.parse_args(argv)
    from production.workflow import run_stationary
    from production.parameter_files import load_parameter_file, claim_output_root
    config = load_parameter_file(args.params)
    output_root = claim_output_root(args.params)
    result, case = run_stationary(
        config["parameters"], config["external_inputs"], config["native_overrides"],
        price_guess=config["price_guess"], budget_seconds=config["budget_seconds"],
        closure=config["closure"], output_root=output_root,
        parameter_file_metadata={key: config[key] for key in
                                 ("config_source", "config_sha256", "config_text", "provenance")})
    print("Completed local stationary GE.")
    status = (result.closure.get("status", "converged")
              if isinstance(result.closure, dict) else "converged")
    print(f"Status: {status}; price={result.price:.12g}")
    if isinstance(result.closure, dict):
        for key in ("renewal_residual", "absolute_housing_residual",
                    "actual_paygo_residual", "population_scale",
                    "fixed_h0_population_scale", "implied_H0_at_population_one"):
            if key in result.closure:
                print(f"{key}: {result.closure[key]}")
        if {"fixed_h0_population_scale", "implied_H0_at_population_one"} <= set(result.closure):
            print("Scale interpretations: fixed-H0 population scale and population-one implied H0.")
    print(f"Completed case: {case}")
    print(f"Latest cached case: {output_root / 'latest'}")
    return case


if __name__ == "__main__":
    main()
