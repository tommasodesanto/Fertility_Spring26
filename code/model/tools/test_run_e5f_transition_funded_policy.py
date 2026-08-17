from __future__ import annotations

import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import run_e5f_transition_funded_policy as funded


def test_selected_cases_always_put_status_quo_first() -> None:
    assert funded.selected_cases(["rebated_tax2"]) == [
        "unrebated_tax1_status_quo",
        "rebated_tax2",
    ]


def test_smoke_always_includes_first_reform_date() -> None:
    assert funded.effective_post_2023_periods(40, True) == 1
    assert funded.effective_post_2023_periods(0, True) == 1
    assert funded.effective_post_2023_periods(40, False) == 40


def test_policy_timing_starts_after_matched_2023() -> None:
    assert funded.POLICY_START_PERIOD == funded.calibration.TRANSITION_PERIODS + 1
    assert funded.POLICY_START_YEAR == 2027
    assert funded.POLICY_CASES["unrebated_tax1_status_quo"][2] is False
    assert funded.POLICY_CASES["rebated_tax1_baseline"][2] is True


def _synthetic_paths() -> dict[str, list[dict[str, float | int | bool]]]:
    paths: dict[str, list[dict[str, float | int | bool]]] = {}
    for case, (_, _, is_funded) in funded.POLICY_CASES.items():
        rows = []
        for period in range(funded.POLICY_START_PERIOD + 1):
            rows.append(
                {
                    "period": period,
                    "year": funded.calibration.START_YEAR + 4 * period,
                    "policy_active": bool(is_funded and period >= funded.POLICY_START_PERIOD),
                    "population_index": 1.2 if period == funded.calibration.TRANSITION_PERIODS else 1.0,
                    "asset_price_index": 1.1 if period == funded.calibration.TRANSITION_PERIODS else 1.0,
                    "asset_price": 2.0,
                    "adult_population": 3.0,
                    "birth_children_topcode_adjusted": 0.4,
                    "entry_flow_E": 0.2,
                    "owner_rate": 0.6,
                    "dependent_child_owner_rate": 0.7,
                    "topcode_adjusted_births_per_adult": 0.1,
                }
            )
        paths[case] = rows
    return paths


def test_paths_are_identical_through_2023_and_reforms_start_in_2027() -> None:
    paths = _synthetic_paths()
    assert funded.prepolicy_identity_gap(paths) == 0.0
    gates = funded.matched_2023_and_policy_timing_gates(
        paths,
        {
            "matched_2023_population_index": 1.2,
            "matched_2023_housing_cost_index": 1.1,
        },
    )
    assert gates == {
        "matched_2023_population_gap": 0.0,
        "matched_2023_housing_cost_gap": 0.0,
    }
    assert paths["unrebated_tax1_status_quo"][5]["policy_active"] is False
    assert paths["rebated_tax1_baseline"][5]["policy_active"] is True


def test_code_bundle_hash_is_deterministic_sha256() -> None:
    first = funded.code_bundle_sha256()
    second = funded.code_bundle_sha256()
    assert first == second
    assert re.fullmatch(r"[0-9a-f]{64}", first)
    required = {
        "code/model/tools/audit_closed_reproductive_closure.py",
        "code/model/intergen_eqscale_seq_optimized/e5f_income_entry_profile.py",
        "code/model/intergen_eqscale_seq_optimized/e6b_profile.py",
        "code/model/intergen_eqscale_seq_optimized/e5f_floor_profile.py",
        "code/model/intergen_eqscale_seq_optimized/e5_profile.py",
    }
    assert required.issubset(set(funded.CODE_BUNDLE_FILES))


def test_bundle_contains_all_top_level_package_python_and_launcher_imports_contract() -> None:
    expected_package = tuple(
        path.relative_to(funded.ROOT).as_posix()
        for path in sorted(
            (
                funded.ROOT
                / "code/model/intergen_eqscale_seq_optimized"
            ).glob("*.py"),
            key=lambda item: item.name,
        )
    )
    assert funded.PACKAGE_CODE_BUNDLE_FILES == expected_package
    assert funded.CODE_BUNDLE_FILES == funded.TOOL_CODE_BUNDLE_FILES + expected_package
    launcher = (
        funded.ROOT / "code/cluster/submit_e5f_transition_funded_policy.sh"
    ).read_text(encoding="utf-8")
    assert "funded.code_bundle_sha256(root)" in launcher
    assert "list(funded.CODE_BUNDLE_FILES)" in launcher
    assert "files = (" not in launcher


def test_stationary_endpoint_uses_the_status_quo_transition_contract(
    monkeypatch,
    tmp_path: Path,
) -> None:
    captured = {}

    def fake_endpoint(**kwargs):
        captured.update(kwargs)
        return {
            "status": "complete",
            "asset_price": 0.7,
            "stationary_population_scale": 0.9,
            "renewal_denominator": 0.1,
            "fixed_stock_relative_market_gap": 1.0e-6,
            "psi_child": -0.25,
        }

    monkeypatch.setattr(funded.transition, "solve_stationary_open_endpoint", fake_endpoint)
    inputs = {
        "chain": object(),
        "base_overrides": {"property_tax_lump_sum_transfer": 0.0},
        "baseline_price": 0.8,
        "outside_flow": 0.01,
        "retention": 0.7,
        "conversion": 0.5,
        "maturation_survival_yield": 0.8,
        "raw_maturation_survival_yield": 0.8,
        "supply_rule": SimpleNamespace(mode="static-elastic"),
    }
    endpoint = funded.solve_status_quo_stationary_endpoint(
        inputs,
        {"new_psi_child": -0.25},
        tmp_path,
    )
    assert captured["policy_case"] == "none"
    assert captured["new_psi_child"] == -0.25
    assert captured["outside_flow"] == inputs["outside_flow"]
    assert captured["retention"] == inputs["retention"]
    assert captured["conversion"] == inputs["conversion"]
    assert captured["maturation_survival_yield"] == 0.8
    assert endpoint["maturation_survival_yield"] == 0.8


def test_persistent_entry_types_clear_jump_and_each_branch_renews(tmp_path: Path) -> None:
    search = tmp_path / "stationary_endpoint_search.csv"
    common = {
        "entry_households_per_adult": 0.1,
        "birth_children_per_adult": 0.04,
        "top_bin_entry_birth_flow_per_adult": 0.01,
        "topcode_adjusted_birth_children_per_adult": 0.05,
        "extra_children_per_top_bin_family": 1.0,
        "renewal_child_accounting_method": "synthetic_common_yield",
    }
    funded.calendar.write_csv(
        search,
        [
            {
                **common,
                "asset_price": 0.499999999,
                "fixed_stock_market_gap": 1.0,
                "fixed_stock_relative_market_gap": 0.25,
                "housing_demand_per_adult": 10.0,
                "housing_user_cost": 1.0,
            },
            {
                **common,
                "asset_price": 0.500000001,
                "fixed_stock_market_gap": -1.0,
                "fixed_stock_relative_market_gap": -0.25,
                "housing_demand_per_adult": 5.0,
                "housing_user_cost": 2.0,
            },
        ],
    )
    endpoint = funded.refine_stationary_jump_persistent_types(
        {
            "status": "complete",
            "asset_price": 0.5,
            "fixed_stock_relative_market_gap": -0.25,
            "psi_child": -0.2,
        },
        search,
        {
            "supply_rule": SimpleNamespace(
                mode="static-elastic",
                quantity=lambda prices: np.array([4.0]),
            ),
            "outside_flow": 0.05,
            "retention": 0.5,
            "conversion": 0.5,
            "maturation_survival_yield": 1.0,
            "baseline_price": 1.0,
            "parameters": SimpleNamespace(user_cost_rate=0.3),
        },
    )
    assert abs(endpoint["fixed_stock_relative_market_gap"]) <= 1.0e-10
    weight = endpoint["persistent_type_population_share_positive_branch"]
    assert 0.0 < weight < 1.0
    assert endpoint["market_clearing_method"].startswith("persistent entry-type")
    assert endpoint["positive_branch_renewal_denominator"] > 0.0
    assert endpoint["negative_branch_renewal_denominator"] > 0.0
    assert abs(endpoint["positive_branch_renewal_identity_residual"]) <= 1.0e-12
    assert abs(endpoint["negative_branch_renewal_identity_residual"]) <= 1.0e-12
    assert abs(endpoint["outside_entry_flow_identity_residual"]) <= 1.0e-12
    assert abs(
        endpoint["positive_branch_outside_entry_flow"]
        + endpoint["negative_branch_outside_entry_flow"]
        - 0.05
    ) <= 1.0e-12
    assert abs(
        endpoint["implied_outside_entrant_share_positive_branch"]
        + endpoint["implied_outside_entrant_share_negative_branch"]
        - 1.0
    ) <= 1.0e-12
    assert endpoint["housing_user_cost"] == 0.3 * endpoint["asset_price"]
    raw_B_over_E = (
        endpoint["mature_entrant_households_per_adult"]
        / endpoint["entry_households_per_adult"]
    )
    adjusted_B_over_E = (
        endpoint["topcode_adjusted_mature_entrant_households_per_adult"]
        / endpoint["entry_households_per_adult"]
    )
    assert endpoint["reproduction_ratio_B_over_E"] == raw_B_over_E
    assert endpoint["raw_state_B_over_E"] == raw_B_over_E
    assert endpoint["topcode_adjusted_B_over_E_diagnostic"] == adjusted_B_over_E
    assert endpoint["topcode_consistent_B_over_E"] == adjusted_B_over_E
    assert endpoint["queue_mature_entrant_flow_B"] == endpoint[
        "topcode_adjusted_mature_entrant_households_per_adult"
    ]
    assert endpoint["queue_birth_children_raw_explicit_states"] == endpoint[
        "birth_children_per_adult"
    ]
    assert endpoint["queue_birth_children_topcode_adjusted"] == endpoint[
        "topcode_adjusted_birth_children_per_adult"
    ]
    assert endpoint["topcode_adjusted_birth_children_per_adult"] == (
        endpoint["birth_children_per_adult"]
        + endpoint["extra_children_per_top_bin_family"]
        * endpoint["top_bin_entry_birth_flow_per_adult"]
    )
    assert endpoint["mature_entrant_households_per_adult"] == (
        endpoint["entrant_conversion_factor"]
        * endpoint["maturation_exit_yield"]
        * endpoint["birth_children_per_adult"]
    )
    assert endpoint["topcode_adjusted_mature_entrant_households_per_adult"] == (
        endpoint["entrant_conversion_factor"]
        * endpoint["maturation_exit_yield"]
        * endpoint["topcode_adjusted_birth_children_per_adult"]
    )
    assert endpoint["threeplus_share_of_mature_flow"] == (
        endpoint["top_bin_entry_birth_flow_per_adult"]
        / endpoint["birth_children_per_adult"]
    )
    assert endpoint["scale_mapping_Hs_over_D"] == (
        endpoint["housing_supply"] / endpoint["housing_demand_per_adult"]
    )
    assert abs(
        endpoint["scale_mapping_Hs_over_D"]
        - endpoint["stationary_population_scale"]
    ) <= 1.0e-10
    assert endpoint["stationary_housing_supply"] == endpoint["housing_supply"]
    assert endpoint["stationary_total_housing_demand"] == (
        endpoint["stationary_population_scale"]
        * endpoint["housing_demand_per_adult"]
    )
    assert "solve_seconds" not in endpoint


def test_transition_diagnostic_writes_stable_png_and_pdf(tmp_path: Path) -> None:
    outputs = funded.write_transition_diagnostic(_synthetic_paths(), tmp_path)
    assert Path(outputs["png"]).name == "funded_policy_transition_paths.png"
    assert Path(outputs["pdf"]).name == "funded_policy_transition_paths.pdf"
    assert Path(outputs["png"]).stat().st_size > 0
    assert Path(outputs["pdf"]).stat().st_size > 0


def test_legacy_report_defaults_to_floor_profile() -> None:
    profile = funded.report_model_profile({}, {})
    assert profile["name"] == funded.calibration.BASELINE_MODEL_PROFILE
    assert profile["legacy_default"] is True


def test_report_profile_rejects_conflicting_contracts() -> None:
    try:
        funded.report_model_profile(
            {"model_profile": {"name": funded.calibration.REPAIRED_MODEL_PROFILE}},
            {"model_profile": funded.calibration.BASELINE_MODEL_PROFILE},
        )
    except RuntimeError as error:
        assert "Conflicting model profiles" in str(error)
    else:
        raise AssertionError("conflicting report profiles must fail")


def test_income_entry_contract_preserves_cost_and_measured_income() -> None:
    theta = {"first_birth_fixed_cost": 1.75}
    _, _, declared = funded.calibration.activate_model_profile(
        funded.calibration.REPAIRED_MODEL_PROFILE,
        dict(theta),
    )
    contract = {
        "model_profile": {
            "name": funded.calibration.REPAIRED_MODEL_PROFILE,
            "report_metadata": declared,
        }
    }
    domain, overrides, metadata = funded.activate_contract_profile(contract, theta)
    assert domain[-1][0] == "first_birth_fixed_cost"
    assert theta["first_birth_fixed_cost"] == 1.75
    assert metadata["name"] == funded.calibration.REPAIRED_MODEL_PROFILE
    assert overrides["permanent_income_levels_enabled"] is True
    assert len(np.asarray(overrides["z_grid"])) == 15


def test_income_entry_report_variance_nz_and_cost_semantics_are_exact() -> None:
    theta = {"first_birth_fixed_cost": 1.75}
    _, _, declared = funded.calibration.activate_model_profile(
        funded.calibration.REPAIRED_MODEL_PROFILE,
        dict(theta),
    )
    for field, replacement in (
        ("permanent_income_log_variance", float(declared["permanent_income_log_variance"]) + 0.01),
        ("income_state_count", int(declared["income_state_count"]) - 1),
        ("first_birth_fixed_cost_semantics", "paid every period"),
    ):
        malformed = dict(declared)
        malformed[field] = replacement
        contract = {
            "model_profile": {
                "name": funded.calibration.REPAIRED_MODEL_PROFILE,
                "report_metadata": malformed,
            }
        }
        try:
            funded.activate_contract_profile(contract, dict(theta))
        except RuntimeError:
            pass
        else:
            raise AssertionError(f"mismatched repaired-profile field must fail: {field}")


def test_report_provenance_requires_selected_share_and_calibration_bundle() -> None:
    live = {"bundle_sha256": "a" * 64, "files": {"solver": "b" * 64}}
    payload = {
        "outside_origin_entry_share": 0.169,
        "code_fingerprints": live,
    }
    profile = {
        "name": funded.calibration.REPAIRED_MODEL_PROFILE,
        "legacy_default": False,
    }
    provenance = funded.validate_transition_report_provenance(
        payload,
        {},
        profile,
        cli_outside_origin_share=0.169,
        live_calibration_code_fingerprints=live,
        expected_model_profile=funded.calibration.REPAIRED_MODEL_PROFILE,
    )
    assert provenance["outside_origin_entry_share"] == 0.169
    assert provenance["legacy_provenance_exception"] is False
    try:
        funded.validate_transition_report_provenance(
            payload,
            {},
            profile,
            cli_outside_origin_share=0.17,
            live_calibration_code_fingerprints=live,
            expected_model_profile=funded.calibration.REPAIRED_MODEL_PROFILE,
        )
    except RuntimeError as error:
        assert "Outside-origin" in str(error)
    else:
        raise AssertionError("CLI/report outside-share mismatch must fail")
    mismatched_live = {"bundle_sha256": "c" * 64, "files": {"solver": "d" * 64}}
    try:
        funded.validate_transition_report_provenance(
            payload,
            {},
            profile,
            cli_outside_origin_share=0.169,
            live_calibration_code_fingerprints=mismatched_live,
            expected_model_profile=funded.calibration.REPAIRED_MODEL_PROFILE,
        )
    except RuntimeError as error:
        assert "code-bundle" in str(error)
    else:
        raise AssertionError("report/live calibration-code mismatch must fail")


def test_report_provenance_reader_accepts_canonical_nested_fields() -> None:
    payload = {
        "contract": {
            "outside_origin_entry_share": 0.169,
            "code_fingerprints": {
                "bundle_sha256": "a" * 64,
                "files": {"solver": "b" * 64},
            },
        }
    }
    assert funded.report_outside_origin_share(payload, {}) == 0.169
    assert funded.report_calibration_code_fingerprints(payload, {})["bundle_sha256"] == "a" * 64


def test_missing_provenance_is_allowed_only_for_explicit_legacy_report() -> None:
    live = {"bundle_sha256": "a" * 64, "files": {"solver": "b" * 64}}
    legacy = {
        "name": funded.calibration.BASELINE_MODEL_PROFILE,
        "legacy_default": True,
    }
    provenance = funded.validate_transition_report_provenance(
        {},
        {},
        legacy,
        cli_outside_origin_share=0.169,
        live_calibration_code_fingerprints=live,
        expected_model_profile=funded.calibration.BASELINE_MODEL_PROFILE,
    )
    assert provenance["legacy_provenance_exception"] is True
    assert set(provenance["legacy_missing_fields"]) == {
        "outside_origin_entry_share",
        "code_fingerprints.bundle_sha256",
    }
    try:
        funded.validate_transition_report_provenance(
            {},
            {},
            legacy,
            cli_outside_origin_share=0.169,
            live_calibration_code_fingerprints=live,
            expected_model_profile=None,
        )
    except RuntimeError as error:
        assert "explicitly" in str(error)
    else:
        raise AssertionError("implicit legacy routing must fail")


def test_income_entry_report_cannot_silently_default_missing_cost() -> None:
    contract = {
        "model_profile": {
            "name": funded.calibration.REPAIRED_MODEL_PROFILE,
            "report_metadata": {},
        }
    }
    try:
        funded.activate_contract_profile(contract, {})
    except RuntimeError as error:
        assert "must estimate and save first_birth_fixed_cost" in str(error)
    else:
        raise AssertionError("missing repaired-profile cost must fail")


def test_income_entry_cost_must_respect_declared_domain() -> None:
    _, _, declared = funded.calibration.activate_model_profile(
        funded.calibration.REPAIRED_MODEL_PROFILE,
        {"first_birth_fixed_cost": 0.0},
    )
    contract = {
        "model_profile": {
            "name": funded.calibration.REPAIRED_MODEL_PROFILE,
            "report_metadata": declared,
        }
    }
    for value in (-0.01, 8.01):
        try:
            funded.activate_contract_profile(
                contract,
                {"first_birth_fixed_cost": value},
            )
        except RuntimeError as error:
            assert "[0,8]" in str(error)
        else:
            raise AssertionError("out-of-domain first-birth cost must fail")


def test_fiscal_policy_uses_period_tax_and_six_room_grant_threshold() -> None:
    parameters = SimpleNamespace(
        period_years=4,
        q=0.10,
        delta=0.08,
        H_own=np.array([2.0, 4.0, 6.0, 8.0, 10.0]),
    )
    funded.set_fiscal_policy(
        parameters,
        annual_tax=0.02,
        grant=True,
        transfer=0.3,
    )
    assert parameters.tau_H == 0.08
    assert parameters.user_cost_rate == 0.26
    assert parameters.property_tax_lump_sum_transfer == 0.3
    assert parameters.birth_entry_grant_amount == 0.4
    np.testing.assert_array_equal(
        parameters.birth_entry_grant_owner_rungs,
        np.array([3, 4, 5]),
    )
