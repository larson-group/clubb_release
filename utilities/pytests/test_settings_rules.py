"""Focused tests for the maintained no-build settings rule list."""

from utilities.clubb_settings_validation import (
    audit_gaps,
    parameter_activity,
    resolve_clubb_settings,
    validate_clubb_settings,
)


def _codes(**kwargs):
    return {issue.code for issue in validate_clubb_settings(**kwargs)}


def test_c10_is_explicitly_unused():
    assert resolve_clubb_settings(parameters={"C10": 1.0}).parameter_states["C10"]["state"] == "unused"
    assert resolve_clubb_settings(parameters={"C13": 1.0}).parameter_states["C13"]["state"] == "unused"


def test_hand_maintained_audit_covers_each_parameter_once():
    assert audit_gaps() == {
        "missing_parameters": (),
        "duplicate_parameters": (),
        "missing_flags": (),
        "duplicate_flags": (),
    }


def test_pdf_parameter_activity():
    assert resolve_clubb_settings({"iiPDF_type": 4}).parameter_states["coef_spread_DG_means_rt"]["state"] == "active"
    assert resolve_clubb_settings({"iiPDF_type": 7}).parameter_states["coef_spread_DG_means_rt"]["state"] == "inactive-mode"


def test_required_equal_parameter_pairs():
    assert "required_equal_parameters" in _codes(parameters={"C6rt": 1.0, "C6thl": 2.0})


def test_em_wp2_constraint():
    assert "em_wp2_damping_inconsistent" in _codes(parameters={"C1": 1.0, "C14": 2.0}, flags={"l_damp_wp2_using_em": True})


def test_xp2_clipping_constraint():
    assert "incompatible_xp2_clipping" in _codes(parameters={}, flags={"l_min_xp2_from_corr_wx": True, "l_enable_relaxed_clipping": True})


def test_wind_prediction_constraint():
    assert "predict_wind_pdf_type" in _codes(parameters={}, flags={"l_predict_upwp_vpwp": True, "iiPDF_type": 8})
    assert "predict_wind_requires_implicit_ta" in _codes(parameters={}, flags={"l_predict_upwp_vpwp": True}, explicit_turbulent_adv_wpxp=True)


def test_pdf_input_field_constraint():
    assert "pdf_requires_input_fields" in _codes(parameters={}, flags={"iiPDF_type": 2})


def test_range_constraints():
    assert "invalid_pdf_call_placement" in _codes(parameters={}, flags={"ipdf_call_placement": 9})
    assert "invalid_saturation_formula" in _codes(parameters={}, flags={"saturation_formula": 9})


def test_advance_order_constraint():
    assert "duplicate_advance_order" in _codes(parameters={}, advance_orders={"order_xm_wpxp": 1, "order_xp2_xpyp": 1})


def test_tau_lscale_warning():
    assert "tau_lscale_parameter_warning" in _codes(parameters={"C1": 2.0}, flags={"l_diag_Lscale_from_tau": True})


def test_richardson_closures_replace_their_tunable_functions():
    c7 = parameter_activity(flags={"l_use_C7_Richardson": True})
    assert all(c7[name][0] == "inactive-mode" for name in ("C7", "C7b", "C7c", "C7_Lscale0"))
    c11 = parameter_activity(flags={"l_use_C11_Richardson": True})
    assert all(c11[name][0] == "inactive-mode" for name in ("C11", "C11b", "C11c"))


def test_lscale_mode_switches_the_correct_parameter_families():
    direct = parameter_activity(flags={"l_diag_Lscale_from_tau": False})
    assert all(direct[name][0] == "inactive-mode" for name in (
        "C_invrs_tau_bkgnd", "C_invrs_tau_N2", "C_invrs_tau_wpxp_N2_thresh",
    ))
    tau = parameter_activity(flags={"l_diag_Lscale_from_tau": True})
    assert all(tau[name][0] == "inactive-mode" for name in (
        "Lscale_mu_coef", "Lscale_pert_coef", "C6rt_Lscale0", "C6thl_Lscale0",
    ))


def test_branch_controls_disable_only_their_owned_parameters():
    precipitation = parameter_activity(flags={"l_use_precip_frac": False})
    assert precipitation["omicron"][0] == "inactive-mode"
    assert precipitation["zeta_vrnce_rat"][0] == "inactive-mode"
    radiation = parameter_activity(flags={"l_calc_thlp2_rad": False})
    assert radiation["thlp2_rad_coef"][0] == "inactive-mode"
    assert radiation["thlp2_rad_cloud_frac_thresh"][0] == "inactive-mode"


def test_lscale_plume_centered_constraint():
    assert "lscale_plume_centered_requires_avg" in _codes(parameters={}, flags={"l_Lscale_plume_centered": True}, l_avg_Lscale=False)


def test_radiation_constraints():
    assert "thlp2_rad_silhs_rad_conflict" in _codes(parameters={}, flags={"l_calc_thlp2_rad": True, "l_silhs_rad": True}, rad_scheme="rrtmg")
    assert "soil_veg_requires_radiation" in _codes(parameters={}, flags={"l_soil_veg": True}, rad_scheme="none")


def test_microphysics_constraints():
    assert "coamps_requires_predict_nc" in _codes(parameters={}, flags={"l_predict_Nc": False}, microphys_scheme="coamps")
    assert "kk_requires_diagnosed_nc" in _codes(parameters={}, flags={"l_predict_Nc": True}, microphys_scheme="kk")


def test_host_constraints():
    assert "host_incompatible_flag" in _codes(parameters={}, flags={"l_rtm_nudge": True}, l_implemented=True)
