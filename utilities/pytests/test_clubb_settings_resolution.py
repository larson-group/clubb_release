"""Focused unit tests for settings implications used by Dash and Tune."""

from utilities.clubb_settings_validation import (
    apply_linked_parameter_values,
    build_settings_schema,
    evaluate_settings,
    flag_relationship_groups,
    linked_parameter_groups,
    resolve_clubb_settings,
    setting_key,
)


def test_predict_wind_constrains_pdf_without_choosing_one():
    resolution = resolve_clubb_settings({"l_predict_upwp_vpwp": True})
    assert resolution.allowed_flag_values == {"iiPDF_type": [1, 7]}
    assert not [issue for issue in resolution.issues if issue.severity == "error"]


def test_incompatible_pdf_and_wind_prediction_remain_visible_not_rewritten():
    resolution = resolve_clubb_settings(
        {"iiPDF_type": 4, "l_predict_upwp_vpwp": True},
        auto_correct=True,
    )
    assert resolution.forced_flags == {}
    assert resolution.effective_flags["l_predict_upwp_vpwp"] is True
    assert {"flags:iiPDF_type", "flag:l_predict_upwp_vpwp"} <= set(resolution.conflicting_setting_keys)


def test_hybrid_pdf_remains_compatible_with_wind_flux_prediction():
    resolution = resolve_clubb_settings(
        {"iiPDF_type": 7, "l_predict_upwp_vpwp": True},
    )

    assert resolution.forced_flags == {}
    assert not [issue for issue in resolution.issues if issue.code == "predict_wind_pdf_type"]


def test_equal_parameter_input_exposes_a_follower():
    resolution = resolve_clubb_settings(parameters={"C6rt": 2.0})
    assert resolution.forced_parameters == {"C6thl": 2.0}
    assert resolution.parameter_states["C6thl"]["state"] == "linked"


def test_linked_parameter_groups_include_the_conditional_em_pair():
    base_groups = linked_parameter_groups()
    em_groups = linked_parameter_groups({"l_damp_wp2_using_em": True})

    assert ("C6rt", "C6thl") in base_groups
    assert ("C1", "C14") not in base_groups
    assert ("C1", "C14") in em_groups


def test_boolean_relationships_remain_declarative_and_surface_constrained_members():
    relationships = flag_relationship_groups()
    assert any(item["label"] == "exclusive" for item in relationships)
    assert any(item["members"] == ["l_damp_wp2_using_em", "l_stability_correct_tau_zm"] for item in relationships)

    clipping = resolve_clubb_settings(
        {"l_min_xp2_from_corr_wx": True, "l_enable_relaxed_clipping": True},
        auto_correct=True,
    )
    assert clipping.forced_flags == {}
    assert {"flag:l_min_xp2_from_corr_wx", "flag:l_enable_relaxed_clipping"} <= set(clipping.conflicting_setting_keys)

    damping = resolve_clubb_settings({"l_damp_wp2_using_em": True}, auto_correct=True)
    assert damping.effective_flags["l_damp_wp2_using_em"] is True
    assert damping.forced_flags == {}


def test_fortran_setting_names_are_case_insensitive_at_the_python_boundary():
    resolution = resolve_clubb_settings({"iipdf_type": 7})
    assert resolution.effective_flags["iiPDF_type"] == 7


def test_pdf_four_and_seven_owners_and_c10_are_explicitly_classified():
    new_pdf = resolve_clubb_settings({"iiPDF_type": 4})
    hybrid_pdf = resolve_clubb_settings({"iiPDF_type": 7})
    assert new_pdf.parameter_states["coef_spread_DG_means_rt"]["state"] == "active"
    assert new_pdf.parameter_states["slope_coef_spread_DG_means_w"]["state"] == "inactive-mode"
    assert hybrid_pdf.parameter_states["slope_coef_spread_DG_means_w"]["state"] == "active"
    assert hybrid_pdf.parameter_states["coef_spread_DG_means_rt"]["state"] == "inactive-mode"
    assert hybrid_pdf.parameter_states["C10"]["state"] == "unused"


def test_evaluator_owns_overrides_link_state_and_invalid_precedence():
    schema = build_settings_schema(
        {"l_predict_upwp_vpwp": False},
        {
            "flags": {"iiPDF_type": "1"},
            "tunable": {"C6rt": "1.0", "C6thl": "1.0"},
            "silhs": {},
        },
        [
            {"file": "flags", "name": "iiPDF_type"},
            {"file": "tunable", "name": "C6rt"},
            {"file": "tunable", "name": "C6thl"},
        ],
    )
    values = apply_linked_parameter_values(
        {
            setting_key("flags", "iiPDF_type"): "4",
            setting_key("tunable", "C6rt"): "2.0",
            setting_key("tunable", "C6thl"): "1.0",
        },
        {"C6rt=C6thl": "2.0"},
    )
    result = evaluate_settings(
        schema,
        flag_values={"l_predict_upwp_vpwp": True},
        parameter_values=values,
    )

    assert result["overrides"]["flags"]["l_predict_upwp_vpwp"] == ".true."
    assert result["overrides"]["flags"]["iiPDF_type"] == "4"
    assert result["overrides"]["tunable"] == {"C6rt": "2.0", "C6thl": "2.0"}
    assert result["control_states"]["flag:l_predict_upwp_vpwp"]["state"] == "invalid"
    assert result["control_states"]["flags:iiPDF_type"]["state"] == "invalid"
    assert result["linked_group_states"]["C6rt=C6thl"]["state"] == "linked-changed"


def test_evaluator_marks_nonboolean_flag_value_errors_on_the_actual_control():
    schema = build_settings_schema(
        {},
        {"flags": {"iiPDF_type": "1"}, "tunable": {}, "silhs": {}},
        [{"file": "flags", "name": "iiPDF_type"}],
    )

    result = evaluate_settings(
        schema,
        parameter_values={setting_key("flags", "iiPDF_type"): "99"},
    )

    assert result["control_states"]["flags:iiPDF_type"]["state"] == "invalid"
