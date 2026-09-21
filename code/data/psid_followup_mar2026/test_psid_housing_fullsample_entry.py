"""No-data structural tests for the portable PSID staging entry point."""
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
    assert 'unweighted author-style extension' in DRIVER


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


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
    print("PSID full-sample entry tests passed")
