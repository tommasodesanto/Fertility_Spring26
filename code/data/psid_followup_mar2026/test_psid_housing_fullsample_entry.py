"""Structural and no-data tests for the portable PSID staging entry point."""
import subprocess
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
DRIVER = (HERE / "psid_housing_fullsample_driver.do").read_text()
LAUNCHER = (HERE / "launch_psid_housing_fullsample_torch.sh").read_text()


def test_driver_is_portable_and_prioritizes_exact_first_arms():
    assert 'args source outroot arm variant' in DRIVER
    assert 'sa_rooms_first_birth_household_aligned_v1.do' in LAUNCHER
    assert 'sa_replication_own_only.do' in LAUNCHER
    assert 'eventstudyinteract' in DRIVER
    assert 'source_sha256_known' in LAUNCHER
    assert '[pw=IW]' in DRIVER
    assert 'unweighted corrected extension; F6 included' in DRIVER
    assert 'assert inlist(HOMEOWN, 0, 1) | missing(HOMEOWN)' in DRIVER
    assert 'STATA_COMPLETE' in LAUNCHER
    assert 'rooms_s_c_y_`variant\'_covariance.csv' in LAUNCHER
    assert 'gen double estimate = b[1,`hook_bp\'] - b[1,`hook_bm\']' in LAUNCHER


def test_author_ownership_is_marked_unweighted():
    marker = 'exact author command intentionally has no [pw=IW]'
    assert marker in LAUNCHER
    assert 'first_birth_ownership' in LAUNCHER


def test_event_contrast_and_full_covariance_are_required():
    assert 'e(V_iw)' in DRIVER
    assert 'L3event' in DRIVER and 'F1event' in DRIVER
    assert 'contrast_ci_lo' in DRIVER and 'contrast_ci_hi' in DRIVER


def test_launcher_requires_remote_source_and_stata_without_ssh_side_effects():
    assert '--source' in LAUNCHER and '--out' in LAUNCHER
    assert '[[ -f "$SOURCE" ]]' in LAUNCHER
    assert 'command -v "$STATA_BIN"' in LAUNCHER
    assert 'ssh ' not in LAUNCHER and 'scp ' not in LAUNCHER


def test_author_script_generation_resolves_anchored_paths():
    old_source = "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
    old_output = "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
    new_source = "/scratch/td2248/projects/kleven_acs_pilot_20260917/inputs/PSID/PSIDSHELF_MOBILITY.dta"
    new_output = "/scratch/td2248/projects/kleven_acs_pilot_20260917/output/psid_full"
    staged_project = "/scratch/td2248/projects/Fertility_Spring26"
    for name in (
        "sa_rooms_first_birth_household_aligned_v1.do",
        "sa_replication_own_only.do",
        "sa_rooms_second_birth_with_onechild_controls_v1.do",
    ):
        text = (HERE / name).read_text()
        text = text.replace(old_source, new_source).replace(old_output, new_output)
        text = text.replace('local dta  "`root\'/PSID/PSIDSHELF_MOBILITY.dta"', f'local dta  "{new_source}"')
        text = text.replace('local outroot "`project\'/code/data/psid_followup_mar2026/output"', f'local outroot "{new_output}"')
        text = text.replace('local out_root "' + old_output + '"', f'local out_root "{new_output}"')
        text = text.replace("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26", staged_project)
        text = text.replace('local root "/Users/tommasodesanto/Desktop/Projects/Fertility"', 'local root "/scratch/td2248/projects"')
        text = text.replace("/Users/tommasodesanto/Desktop/Projects/Fertility", "/scratch/td2248/projects")
        assert "/Users/" not in text
        assert new_source in text
        assert new_output in text


def test_actual_launcher_generator_runs_for_all_author_templates():
    """Execute the launcher heredoc itself, rather than reimplementing it."""
    old_source = "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
    old_output = "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
    new_source = "/scratch/td2248/projects/kleven_acs_pilot_20260917/inputs/PSID/PSIDSHELF_MOBILITY.dta"
    new_output = "/scratch/td2248/projects/kleven_acs_pilot_20260917/output/psid_full"
    staged_project = "/scratch/td2248/projects/Fertility_Spring26"
    py_source = LAUNCHER.split("<<'PY'\n", 1)[1].split("\nPY\n", 1)[0]
    for name in (
        "sa_rooms_first_birth_household_aligned_v1.do",
        "sa_replication_own_only.do",
        "sa_rooms_second_birth_with_onechild_controls_v1.do",
    ):
        with tempfile.TemporaryDirectory() as td:
            dst = Path(td) / name
            subprocess.run(
                ["python3", "-", str(HERE / name), str(dst), old_source, new_source,
                 old_output, new_output, staged_project, new_output + "/ado"],
                input=py_source,
                text=True,
                check=True,
            )
            text = dst.read_text()
        assert "/Users/" not in text
        assert new_source in text
        assert new_output in text
        assert text.count("mata: mata mlib index") == 1
        assert "clear all\nsysdir set PLUS" in text


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
    print("PSID full-sample entry tests passed")
