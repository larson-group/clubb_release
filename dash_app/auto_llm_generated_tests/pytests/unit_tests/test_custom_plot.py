import numpy as np
import pytest

from dash_app.plot_tab.plot_types import custom_plot


def test_expression_parser_collects_clubb_and_les_variables():
    tree, names = custom_plot.expression_variable_names(
        "np.abs(wpthlp) + 0.61 * thv_ds_zm * les.wprtp + rc_coef_zm * les.wprcp"
    )

    assert names == ["les.wprcp", "les.wprtp", "rc_coef_zm", "thv_ds_zm", "wpthlp"]
    values = {
        "wpthlp": np.array([[-1.0, 2.0]]),
        "thv_ds_zm": np.array([[300.0, 301.0]]),
        "les.wprtp": np.array([[0.001, 0.002]]),
        "rc_coef_zm": np.array([[4.0, 5.0]]),
        "les.wprcp": np.array([[0.01, 0.02]]),
    }
    expected = np.abs(values["wpthlp"]) + 0.61 * values["thv_ds_zm"] * values["les.wprtp"]
    expected += values["rc_coef_zm"] * values["les.wprcp"]

    np.testing.assert_allclose(custom_plot._eval_expression_node(tree, values), expected)


def test_gaussian_wprcp_expression_supports_erf_and_derived_les_chi_fields():
    expression = "les.wpchi * 0.5 * (1 + erf(les.chi / sqrt(2 * les.chip2)))"
    tree, names = custom_plot.expression_variable_names(expression)
    values = {
        "les.wpchi": np.array([[8.0e-5, 4.0e-5]]),
        "les.chi": np.array([[0.0, -2.0e-4]]),
        "les.chip2": np.array([[1.0e-8, 4.0e-8]]),
    }

    assert names == ["les.chi", "les.chip2", "les.wpchi"]
    result = custom_plot._eval_expression_node(tree, values)
    expected_cloud_fraction = 0.5 * (
        1.0
        + custom_plot._erf(
            values["les.chi"] / np.sqrt(2.0 * values["les.chip2"])
        )
    )
    np.testing.assert_allclose(result, values["les.wpchi"] * expected_cloud_fraction)


@pytest.mark.parametrize("expression", ["les", "foo.bar", "les.wprcp.real", "les.__class__"])
def test_expression_parser_rejects_invalid_namespaces(expression):
    with pytest.raises(custom_plot.ExpressionError):
        custom_plot.expression_variable_names(expression)


def test_mixed_expression_interpolates_and_broadcasts_les_profile(monkeypatch):
    plot = custom_plot.CustomPlotType()
    tree, names = custom_plot.expression_variable_names("wp2 + les.wprcp")
    clubb_profiles = {
        "profiles": np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        "z_values": np.array([0.0, 1.0, 2.0]),
        "labels": ["column 1", "column 2"],
        "z_units": "m",
    }
    les_profiles = {
        "wprcp": {
            "profile": np.array([10.0, 30.0]),
            "z_values": np.array([0.0, 2.0]),
            "z_units": "m",
        }
    }

    monkeypatch.setattr(custom_plot.shared, "selected_time_indices", lambda *args: [0])
    monkeypatch.setattr(
        custom_plot.shared,
        "extract_time_avg_profile_for_path",
        lambda *args, **kwargs: clubb_profiles,
    )
    monkeypatch.setattr(plot, "_extract_les_profiles", lambda *args: (les_profiles, "sam"))

    result = plot._extract_expression_profiles(
        "stats.nc",
        "wp2 + les.wprcp",
        tree,
        names,
        {
            "case_data": {},
            "column_mode": "all",
            "selected_column": 0,
        },
    )

    np.testing.assert_allclose(
        result["profiles"],
        np.array([[11.0, 22.0, 33.0], [14.0, 25.0, 36.0]]),
    )
    assert result["labels"] == ["column 1", "column 2"]


def test_les_only_expression_uses_les_grid_and_label(monkeypatch):
    plot = custom_plot.CustomPlotType()
    tree, names = custom_plot.expression_variable_names("2 * les.wprcp")
    les_profiles = {
        "wprcp": {
            "profile": np.array([1.0, 2.0]),
            "z_values": np.array([100.0, 200.0]),
            "z_units": "m",
        }
    }

    monkeypatch.setattr(custom_plot.shared, "selected_time_indices", lambda *args: [0])
    monkeypatch.setattr(plot, "_extract_les_profiles", lambda *args: (les_profiles, "sam"))

    result = plot._extract_expression_profiles(
        "stats.nc",
        "2 * les.wprcp",
        tree,
        names,
        {"case_data": {}},
    )

    np.testing.assert_allclose(result["profiles"], [[2.0, 4.0]])
    np.testing.assert_allclose(result["z_values"], [100.0, 200.0])
    assert result["labels"] == ["LES (SAM)"]
