"""Stage the frozen paper-baseline initial scorer for a local, read-only replay.

This performs no model solve.  It copies the 641-file source snapshot and the
archived scored-run template, relocates path-bearing JSON contracts, and checks
source fingerprints. The scored wrapper performs independent preflight at launch.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
ARCH = ROOT / "output/model/paper_baseline_sep14"
NATIVE = ARCH / "native/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913"
TMP = ROOT / "tmp/paper_baseline_sep14"
OUT = ROOT / "output/model/local_mini_calibration_20260919/staged"
SOURCE = OUT / "source"
TEMPLATE = OUT / "template"
RECIPE = OUT / "recipe"
RUNTIME = OUT / "runtime"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"refusing to overwrite existing staging directory: {OUT}")
    OUT.mkdir(parents=True)

    # The tmp baseline has the exact 641-file inventory and executable pins;
    # restore the pinned README from the source commit because that one file is
    # documentation-only and differs in the worktree copy.
    native_initial = NATIVE / "corrected_initial_template_v6/seed_case/initial_contract.json"
    initial = json.loads(native_initial.read_text())
    source_pins = initial["source_sha256"]
    readme_tmp = TMP / "code/model/README.md"
    readme_old = subprocess.check_output(
        ["git", "show", "70abd4a80df979903b05d361f47ac2bf5aab2c23:code/model/README.md"],
        cwd=ROOT,
    )
    readme_expected = source_pins["code/model/README.md"]
    if hashlib.sha256(readme_old).hexdigest() != readme_expected:
        raise RuntimeError("git source README does not match frozen source pin")
    readme_tmp_sha = sha(readme_tmp)

    for rel, pin in source_pins.items():
        src = TMP / rel
        if not src.is_file():
            raise FileNotFoundError(src)
        dst = SOURCE / rel
        copy_file(src, dst)
        if rel == "code/model/README.md":
            dst.write_bytes(readme_old)
        if sha(dst) != pin:
            raise RuntimeError(f"source pin mismatch: {rel}")

    # Copy the archived template inputs and contracts.  Only JSON path fields
    # are rewritten; executable files retain their archived hashes.
    t_rel = Path("corrected_initial_template_v6")
    for rel in [
        "inputs/run_scored_candidate.py", "inputs/score_initial.py", "inputs/panel_validator.py",
        "inputs/working_contract.json", "seed_case/initial_contract.json", "seed_case/run_contract.json",
    ]:
        copy_file(NATIVE / t_rel / rel, TEMPLATE / rel)
    for name in ("plan_capped_beta_099.json", "run_capped_beta.py", "run_profile.py"):
        copy_file(ARCH / "initial_recipe" / name, TEMPLATE / name)
    for name in ("run_e5f_rebated_initial_overnight.py", "run_e5f_joint_rebated_initial_probe.py",
                 "run_e5f_joint_rebated_initial_scored.py", "e5f_stationary_paygo.py", "proposal.json"):
        copy_file(ARCH / "initial_recipe" / name, RECIPE / name)
    seed_score = NATIVE / "corrected_initial_template_v6/seed_case/evaluation/scored_repetition_01/score.json"
    copy_file(seed_score, TEMPLATE / "seed_case/evaluation/scored_repetition_01/score.json")
    plan = json.loads((TEMPLATE / "plan_capped_beta_099.json").read_text())
    plan["resume_score_path"] = str(TEMPLATE / "seed_case/evaluation/scored_repetition_01/score.json")
    plan["source_root"] = str(SOURCE)
    write_json(TEMPLATE / "plan_capped_beta_099.json", plan)

    # Preserve the runtime pickle namespace compatibility used by the verified
    # Torch replay.  This is loaded through PYTHONPATH by the lead's runner.
    (RUNTIME / "sitecustomize.py").parent.mkdir(parents=True, exist_ok=True)
    (RUNTIME / "sitecustomize.py").write_text(
        "import importlib, pathlib, sys, numpy as np\n"
        "sys.modules.setdefault('numpy._core', np.core)\n"
        "sys.modules.setdefault('numpy._core.multiarray', importlib.import_module('numpy.core.multiarray'))\n"
        "sys.modules.setdefault('numpy._core.numeric', importlib.import_module('numpy.core.numeric'))\n"
        "sys.modules.setdefault('pathlib._local', pathlib)\n"
    )

    # Copy the exact normalized seed when the lead has placed it; staging can
    # still be inspected before that dependency arrives.
    seed = OUT.parent / "normalized_old.pkl.gz"
    seed_pin = initial["normalized_checkpoint_sha256"]
    staged_seed = OUT / "normalized_old.pkl.gz"
    if seed.is_file():
        if sha(seed) != seed_pin:
            raise RuntimeError("normalized seed exists but has the wrong SHA-256")
        copy_file(seed, staged_seed)

    initial_path = TEMPLATE / "seed_case/initial_contract.json"
    initial = json.loads(initial_path.read_text())
    initial["normalized_checkpoint"] = str(staged_seed)
    initial["source_sha256"] = source_pins
    write_json(initial_path, initial)

    # The two canonical-JSON manifest inputs are derived directly from the
    # frozen objective's embedded provenance maps.  This preserves their
    # canonical fingerprints while giving the relocated contract real files.
    objective = json.loads((TEMPLATE / "inputs/working_contract.json").read_text())
    write_json(TEMPLATE / "inputs/economic_source.json",
               objective["source_provenance"]["economic_source"]["source_sha256"])
    write_json(TEMPLATE / "inputs/observation_snapshot.json", source_pins)

    run_path = TEMPLATE / "seed_case/run_contract.json"
    run = json.loads(run_path.read_text())
    run["source_root"] = str(SOURCE)
    run["initial_solve_contract"]["path"] = str(initial_path)
    run["initial_solve_contract"]["sha256"] = sha(initial_path)
    for key in ("working_objective", "scorer", "validator"):
        old = Path(run[key]["path"])
        local = TEMPLATE / "inputs" / old.name
        run[key]["path"] = str(local)
        run[key]["sha256"] = sha(local)
    for key, item in run["objective_source_files"].items():
        # Source/objective files are embedded in the source snapshot where
        # possible; otherwise retain the repository-local archived evidence.
        old = Path(item["path"])
        candidates = [
            TEMPLATE / "inputs/economic_source.json" if key == "economic_source_manifest_7e872053" else None,
            TEMPLATE / "inputs/observation_snapshot.json" if key == "observation_snapshot_manifest_70abd4a8" else None,
            ROOT / old.relative_to(ROOT) if old.is_absolute() and str(old).startswith(str(ROOT)) else None,
            ROOT / "output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/fertility_target_contract.json" if key == "fertility_provenance_file_sha256" else None,
            ROOT / "output/model/e5f_matched_pf_20260909a/design_research/observer_contract/observer_contract.json" if key == "housing_wealth_provenance_file_sha256" else None,
            ROOT / "output/model/e5f_matched_pf_20260909a/initial_calibration_contract/lead_methodological_decision.json" if key == "lead_decision_file_sha256" else None,
            ROOT / "output/model/e5f_matched_pf_20260909a/initial_fit_readout/parameters.csv" if key == "parameter_table_file_sha256" else None,
            TEMPLATE / "inputs/score_initial.py" if key == "scorer_file_sha256" else None,
            ROOT / "output/model/e5f_matched_pf_20260909a/initial_fit_readout/target_fit.csv" if key == "target_fit_file_sha256" else None,
            ROOT / "output/model/e5f_matched_pf_20260909a/initial_calibration_contract/working_weights.csv" if key == "working_weights_file_sha256" else None,
        ]
        local = next((p for p in candidates if p is not None and p.is_file()), None)
        if local is None:
            raise FileNotFoundError(f"cannot relocate objective source {key}: {old}")
        item["path"] = str(local)
    write_json(run_path, run)

    # The working objective remains byte-identical and therefore retains its
    # approved canonical hash; only its path in the run contract is relocated.
    objective_src = NATIVE / t_rel / "inputs/working_contract.json"
    copy_file(objective_src, TEMPLATE / "inputs/working_contract.json")
    contract = dict(run)
    contract["case_id"] = "local_staging_preflight"
    contract["seconds"] = 2100
    contract["wrapper_sha256"] = sha(TEMPLATE / "inputs/run_scored_candidate.py")
    contract["working_objective"]["path"] = str(TEMPLATE / "inputs/working_contract.json")
    contract["working_objective"]["sha256"] = sha(TEMPLATE / "inputs/working_contract.json")
    write_json(TEMPLATE / "run_contract.json", contract)

    receipt = {
        "status": "staged_preflight_only",
        "source_root": str(SOURCE),
        "source_file_count": len(source_pins),
        "readme_worktree_sha256": readme_tmp_sha,
        "readme_frozen_sha256": readme_expected,
        "working_objective_canonical_sha256": "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1",
        "normalized_seed_present": staged_seed.is_file(),
        "normalized_seed_sha256": sha(staged_seed) if staged_seed.is_file() else None,
        "model_solve_run": False,
        "run_command": ["code/model/.venv/bin/python", str(RECIPE / "run_e5f_joint_rebated_initial_scored.py"),
                        "--helper", str(RECIPE / "run_e5f_rebated_initial_overnight.py"),
                        "--helper-sha256", sha(RECIPE / "run_e5f_rebated_initial_overnight.py"),
                        "--joint", str(RECIPE / "run_e5f_joint_rebated_initial_probe.py"),
                        "--joint-sha256", sha(RECIPE / "run_e5f_joint_rebated_initial_probe.py"),
                        "--template", str(TEMPLATE), "--proposal", str(RECIPE / "proposal.json"),
                        "--output", str(OUT / "candidate")],
    }
    write_json(OUT / "README.json", receipt)
    (OUT / "README.md").write_text(
        "# Local paper-baseline staging\n\n"
        "Read-only portable staging for the frozen initial scored candidate. "
        "No model solve is run by the staging script. The source inventory is "
        "641 files; executable pins are preserved. The worktree README differed "
        "from the frozen source and was restored from commit `70abd4a8`; see "
        "`README.json` for hashes and seed status.\n"
    )


if __name__ == "__main__":
    main()
