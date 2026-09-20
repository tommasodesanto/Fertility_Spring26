#!/usr/bin/env bash
set -euo pipefail
E="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/rental_wedge_exhaustive_v1"; F="/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source_v2"
check() { local got; got="$(sha256sum "$1" | cut -d' ' -f1)"; [[ "$got" == "$2" ]] || { echo "hash mismatch: $1 got $got expected $2" >&2; exit 4; }; }
check "$E/source_manifest.json" "f9cec6e92ddecebedfa6780449700191fef317f910920ea94eb98c699eb838db"
check "$E/plan.json" "1ae0ba4dfc1a66088e1f3803e80a0ba888f08b46874da3daf0ef4fab6ef76c4a"
check "$E/code/model/tools/run_e5f_isolated_rental_wedge.py" "9b059d7301ab450bf8cab401a73199850dd328b79566491dca08fd86983aad7b"
for pin in build_e5f_native_financing_report.py:beccdda8ac758326c5ff75943831ef402c1e82d1e8a185eebaaac6c079ca65b5 build_persistent_transitory_income_candidate.py:7c9571fb38697e9433b32baabbd630656bd935179a4325f7b1e3bbb70e94d087 run_e5f_financing_factorial.py:0cb6d5b1db5f1ab7b77912099f27ffb8908f8e446574ffcb3d69d724eda364b6 run_e5f_income_candidate_calibration.py:e688e4fa401097993dcd6cbca1b7fc9afe6911a2086e7bc74925c962bbaf92f2 run_e5f_income_candidate_search.py:1286ba54f1d2e5a64798150de3c6786c109d1f28733c0e47a8e769dad339fa18 run_e5f_income_overnight_search.py:111776865158779fa9ab9de05308a0b7d3e98adf8b0280cfaf2ed159be58eb48 run_e5f_native_financing_diagnostic.py:ff65058d8e796069f82fcaf463256e906a500af876bb91d91c646c200d3582f2 run_e5f_native_income_cohort_diagnostic.py:1ca604b4b6e83f5cd43b48163a4889512165caf5a4ffb88f6b8759c60a0e1d50 run_e5f_native_rental_access_diagnostic.py:90d48059a3f76668b5e6c4ea0f789de14c2280e21ded7493dc41f67ca367a924; do check "$E/code/model/tools/${pin%%:*}" "${pin##*:}"; done
[[ "$(ls "$E/code/model/tools/"*.py | wc -l)" == 10 ]] || { echo "unexpected helper count in $E/code/model/tools" >&2; exit 4; }
python3 - "$E" "$F" <<'PY'
import hashlib, json, pathlib, sys
E, F = map(pathlib.Path, sys.argv[1:3])
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
port = json.loads((E / "source_manifest.json").read_text()); launch = json.loads((E / "launch_manifest.json").read_text())
pins = port["source_files"]; relevant = set(port["relevant_files"]); base = port["base_hashes_frozen"]
present = {str(p.relative_to(E / "source")) for p in (E / "source/code/model").rglob("*.py") if "__pycache__" not in p.parts}
if present != set(pins): raise SystemExit(f"staged Python set differs from the port manifest: extra={sorted(present - set(pins))} missing={sorted(set(pins) - present)}")
for rel, pin in pins.items():
    if sha(E / "source" / rel) != pin: raise SystemExit(f"staged source hash mismatch: {rel}")
rows = {}
for rel, pin in pins.items():
    fpath = F / rel
    if rel in relevant and rel not in port.get("owned_tests", []):
        if not fpath.is_file(): raise SystemExit(f"frozen root lacks relevant file {rel}")
        fh = sha(fpath)
        if fh != base[rel]: raise SystemExit(f"frozen root {rel} hash {fh} != port base {base[rel]}")
        rows[rel] = {"frozen": fh, "ported": pin, "changed": fh != pin}
    elif fpath.is_file():
        fh = sha(fpath)
        if fh != pin: raise SystemExit(f"non-relevant snapshot file differs from the frozen root: {rel}")
        rows[rel] = {"frozen": fh, "ported": pin, "changed": False}
    else:
        if not rel.startswith("code/model/intergen_eqscale_seq_optimized/tests/"): raise SystemExit(f"snapshot-only file outside tests/: {rel}")
        rows[rel] = {"frozen": None, "ported": pin, "changed": True, "snapshot_only": True}
changed = sorted(r for r, v in rows.items() if v["changed"] and not v.get("snapshot_only"))
expected = sorted(r for r, v in launch["changed_source_vs_frozen"].items() if v["changed"] and v["frozen"] is not None)
if changed != expected: raise SystemExit(f"changed-file set {changed} differs from the launch manifest {expected}")
for rel in changed:
    if launch["changed_source_vs_frozen"][rel] != {"frozen": rows[rel]["frozen"], "ported": rows[rel]["ported"], "changed": True}: raise SystemExit(f"changed-file hashes differ from the launch manifest: {rel}")
nonpy = {str(p.relative_to(E / "source")): sha(p) for p in sorted((E / "source").rglob("*")) if p.is_file() and p.suffix != ".py"}
(E / "frozen_comparison.json").write_text(json.dumps({"frozen_root": str(F), "files": rows, "changed_vs_frozen": changed, "snapshot_only": sorted(r for r, v in rows.items() if v.get("snapshot_only")), "python_files": len(rows), "nonpython_files": len(nonpy)}, indent=2, sort_keys=True) + "\n")
(E / "source_nonpython_manifest.json").write_text(json.dumps(nonpy, indent=2, sort_keys=True) + "\n")
print(f"verified {len(rows)} Python files ({len(changed)} patched vs frozen) and {len(nonpy)} non-Python frozen files")
PY
[[ "$(find "$E/source" -name __pycache__ | wc -l)" == 0 ]] || { echo "cache directories present in staged source" >&2; exit 4; }
chmod -R a-w "$E/code" "$E/source" "$E/plan.json" "$E/launch_manifest.json" "$E/source_manifest.json" "$E/smoke.sbatch" "$E/production.sbatch" "$E/remote_stage.sh" "$E/remote_finalize.sh" "$E/driver_plan.json" "$E/frozen_comparison.json" "$E/source_nonpython_manifest.json"
echo "finalized $E"
