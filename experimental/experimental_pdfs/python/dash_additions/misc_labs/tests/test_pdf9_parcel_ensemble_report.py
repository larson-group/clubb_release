"""Focused contracts for the PDF-9 parcel-ensemble Misc laboratory."""

from __future__ import annotations

import numpy as np
import pytest

from .reports import pdf9_parcel_ensemble as report


def _cube():
    parcel_w = np.zeros((1, 2, 3))
    parcel_rt = np.zeros_like(parcel_w)
    parcel_thl = np.zeros_like(parcel_w)
    parcel_buoyancy = np.zeros_like(parcel_w)
    status = np.zeros_like(parcel_w, dtype=np.int8)
    parcel_w[0, 1] = (1.0, 2.0, 3.0)
    parcel_rt[0, 1] = (0.010, 0.012, 0.014)
    parcel_thl[0, 1] = (299.0, 300.0, 301.0)
    parcel_buoyancy[0, 1] = (0.1, 0.2, -0.1)
    status[0, 1] = 1
    down_w = -0.5 * parcel_w
    down_rt = parcel_rt - 0.001
    down_thl = parcel_thl + 0.5
    down_buoyancy = -parcel_buoyancy
    down_status = status.copy()
    candidate_valid_up = np.array(((0.0, 1.0),))
    return {
        "parcel_w": parcel_w,
        "parcel_rt": parcel_rt,
        "parcel_thl": parcel_thl,
        "parcel_buoyancy": parcel_buoyancy,
        "parcel_status": status,
        "down_parcel_w": down_w,
        "down_parcel_rt": down_rt,
        "down_parcel_thl": down_thl,
        "down_parcel_buoyancy": down_buoyancy,
        "down_parcel_status": down_status,
        "candidate_valid_up": candidate_valid_up,
        "candidate_launch_from_g3_up": np.array(((0.0, 1.0),)),
        "candidate_launch_from_g3_down": np.array(((0.0, 1.0),)),
        "candidate_crossing_count_up": np.array(((0.0, 3.0),)),
        "candidate_weighted_support_up": np.array(((0.0, 2.4),)),
        "candidate_donor_distance_up": np.array(((0.0, 80.0),)),
        "candidate_w_up": np.array(((0.0, 2.0),)),
        "candidate_rt_up": np.array(((0.011, 0.012),)),
        "candidate_thl_up": np.array(((299.5, 300.0),)),
        "candidate_var_w_up": np.array(((0.0, 2.0 / 3.0),)),
        "candidate_var_rt_up": np.array(((0.0, 8.0e-6 / 3.0),)),
        "candidate_var_thl_up": np.array(((0.0, 2.0 / 3.0),)),
        "candidate_covar_w_rt_up": np.array(((0.0, 4.0e-3 / 3.0),)),
        "candidate_covar_w_thl_up": np.array(((0.0, 2.0 / 3.0),)),
        "candidate_covar_rt_thl_up": np.array(((0.0, 4.0e-3 / 3.0),)),
        "candidate_corr_w_rt_up": np.array(((0.0, 1.0),)),
        "candidate_corr_w_thl_up": np.array(((0.0, 1.0),)),
        "candidate_corr_rt_thl_up": np.array(((0.0, 1.0),)),
        "candidate_valid_down": np.array(((0.0, 1.0),)),
        "candidate_crossing_count_down": np.array(((0.0, 3.0),)),
        "candidate_weighted_support_down": np.array(((0.0, 1.6),)),
        "candidate_donor_distance_down": np.array(((0.0, 120.0),)),
        "candidate_w_down": np.array(((0.0, -1.0),)),
        "candidate_w_down_uncapped": np.array(((0.0, -1.5),)),
        "candidate_down_cap_fraction": np.array(((0.0, 0.75),)),
        "candidate_destination_sigma_w": np.array(((0.4, 0.8),)),
        "candidate_rt_down": np.array(((0.011, 0.011),)),
        "candidate_thl_down": np.array(((299.5, 300.5),)),
        "candidate_var_w_down": np.array(((0.0, 1.0 / 6.0),)),
        "candidate_var_rt_down": np.array(((0.0, 8.0e-6 / 3.0),)),
        "candidate_var_thl_down": np.array(((0.0, 2.0 / 3.0),)),
        "candidate_covar_w_rt_down": np.array(((0.0, -2.0e-3 / 3.0),)),
        "candidate_covar_w_thl_down": np.array(((0.0, -1.0 / 3.0),)),
        "candidate_covar_rt_thl_down": np.array(((0.0, 4.0e-3 / 3.0),)),
        "candidate_corr_w_rt_down": np.array(((0.0, -1.0),)),
        "candidate_corr_w_thl_down": np.array(((0.0, -1.0),)),
        "candidate_corr_rt_thl_down": np.array(((0.0, 1.0),)),
        "candidate_branch_prob_up": np.array(((0.0, 0.5),)),
        "candidate_branch_prob_down": np.array(((0.0, 0.5),)),
        "candidate_valid_combined": candidate_valid_up.copy(),
        "candidate_crossing_count_combined": np.array(((0.0, 6.0),)),
        "candidate_weighted_support_combined": np.array(((0.0, 4.0),)),
        "candidate_donor_distance_combined": np.array(((0.0, 96.0),)),
        "candidate_w_combined": np.array(((0.0, 0.5),)),
        "candidate_rt_combined": np.array(((0.011, 0.0115),)),
        "candidate_thl_combined": np.array(((299.5, 300.25),)),
        "candidate_var_w_combined": np.array(((0.0, 8.0 / 3.0),)),
        "candidate_var_rt_combined": np.array(((0.0, 2.9166666667e-6),)),
        "candidate_var_thl_combined": np.array(((0.0, 0.7291666667),)),
        "candidate_covar_w_rt_combined": np.array(((0.0, 13.0e-3 / 12.0),)),
        "candidate_covar_w_thl_combined": np.array(((0.0, -5.0 / 24.0),)),
        "candidate_covar_rt_thl_combined": np.array(((0.0, 0.0012083333),)),
        "candidate_corr_w_rt_combined": np.array(((0.0, 0.44),)),
        "candidate_corr_w_thl_combined": np.array(((0.0, -0.26),)),
        "candidate_corr_rt_thl_combined": np.array(((0.0, 0.83),)),
        "g12_mixt_frac": np.array(((0.5, 0.25),)),
        "g1_w": np.array(((0.0, 1.0),)),
        "g2_w": np.array(((0.0, -0.5),)),
        "g1_rt": np.array(((0.011, 0.013),)),
        "g2_rt": np.array(((0.011, 0.0105),)),
        "g1_var_w": np.array(((0.1, 0.4),)),
        "g2_var_w": np.array(((0.1, 0.2),)),
        "g1_var_rt": np.array(((1.0e-6, 2.0e-6),)),
        "g2_var_rt": np.array(((1.0e-6, 1.5e-6),)),
        "g1_corr_w_rt": np.array(((0.0, 0.6),)),
        "g2_corr_w_rt": np.array(((0.0, -0.2),)),
        "g1_chi": np.array(((0.0, 0.0015),)),
        "g2_chi": np.array(((0.0, 0.0005),)),
        "g1_stdev_chi": np.array(((1.0e-5, 3.0e-4),)),
        "g2_stdev_chi": np.array(((1.0e-5, 2.0e-4),)),
        "g1_corr_w_chi": np.array(((0.0, 0.4),)),
        "g2_corr_w_chi": np.array(((0.0, -0.3),)),
        "g1_rc": np.array(((0.0, 0.0012),)),
        "g2_rc": np.array(((0.0, 0.0003),)),
        "g3_weight": np.array(((0.0, 0.2),)),
        "g3_w": np.array(((0.0, 0.8),)),
        "g3_rt": np.array(((0.011, 0.012),)),
        "g3_thl": np.array(((299.5, 300.1),)),
        "g3_var_w": np.array(((0.1, 0.3),)),
        "g3_var_rt": np.array(((1.0e-6, 1.5e-6),)),
        "g3_var_thl": np.array(((0.1, 0.3),)),
        "g3_corr_w_rt": np.array(((0.0, 0.5),)),
        "g3_corr_w_thl": np.array(((0.0, 0.1),)),
        "g3_corr_rt_thl": np.array(((0.0, 0.2),)),
        "g3_chi": np.array(((0.0, 0.0010),)),
        "g3_stdev_chi": np.array(((1.0e-5, 2.5e-4),)),
        "g3_corr_w_chi": np.array(((0.0, 0.6),)),
        "g3_rc": np.array(((0.0, 0.0008),)),
        "mean_rt": np.array(((0.011, 0.011),)),
        "mean_thl": np.array(((299.5, 300.5),)),
        "entrain_g3_weight": np.array(((0.0, 0.2),)),
        "entrain_rt_up": np.array(((0.011, 0.0118),)),
        "entrain_rt_down": np.array(((0.011, 0.0106),)),
        "entrain_thl_up": np.array(((299.5, 300.2),)),
        "entrain_thl_down": np.array(((299.5, 300.8),)),
        "local_var_w": np.array(((0.2, 0.3),)),
        "local_var_rt": np.array(((1.0e-6, 2.0e-6),)),
        "local_cov_w_rt": np.array(((1.0e-4, 2.0e-4),)),
        "height_m": np.array((100.0, 200.0)),
        "time_seconds": np.array((60.0,)),
        "launch_height_m": np.array((40.0, 80.0, 120.0)),
        "sam_reference_mean": np.zeros((1, 2, len(report.REFERENCE_NAMES), 3)),
        "sam_reference_covariance": np.tile(
            np.eye(3), (1, 2, len(report.REFERENCE_NAMES), 1, 1)
        ),
        "sam_reference_available": np.ones(
            (1, 2, len(report.REFERENCE_NAMES)), dtype=np.int8
        ),
    }


def test_speed_weighted_ensemble_reports_effective_support():
    result = report._ensemble(_cube(), 0, 1, "speed")
    assert result["available"]
    np.testing.assert_allclose(result["mean"], (14.0 / 6.0, 0.0126666667, 300.3333333))
    assert result["count"] == 3
    assert result["neff"] == 36.0 / 14.0
    assert np.min(np.linalg.eigvalsh(result["covariance"])) >= -1.0e-12


def test_occupancy_reduction_uses_one_consistent_vote_per_crossing():
    cube = _cube()
    result = report._ensemble(cube, 0, 1, "uniform")
    assert result["available"]
    np.testing.assert_array_equal(result["weight"], (1.0, 1.0, 1.0))
    np.testing.assert_allclose(result["mean"], (2.0, 0.012, 300.0))
    values = np.column_stack(
        (cube["parcel_w"][0, 1], cube["parcel_rt"][0, 1], cube["parcel_thl"][0, 1])
    )
    np.testing.assert_allclose(result["covariance"], np.cov(values, rowvar=False, bias=True))
    assert result["neff"] == 3.0


def test_occupancy_is_explicitly_named_in_weighting_controls():
    assert report.WEIGHTINGS["uniform"] == "Occupancy (equal per crossing)"


def test_masked_cache_scalars_have_safe_boundary_defaults():
    assert report._finite_int(np.nan) == 0
    assert report._finite_int(np.ma.masked) == 0
    assert report._finite_float(np.nan) == 0.0
    assert report._finite_float(np.ma.masked) == 0.0
    assert report._distance_label(report._finite_float(np.nan, np.nan)) == "unavailable"


def test_residual_component_weights_include_g3_outer_weight():
    cube = _cube()
    cube["g3_weight"] = np.array(((0.0, 0.2),))
    g1 = report._active_residual_component(cube, 0, 1, 1)
    g2 = report._active_residual_component(cube, 0, 1, 2)
    assert g1["available"] and g2["available"]
    assert np.isclose(g1["weight"], 0.2)
    assert np.isclose(g2["weight"], 0.6)
    assert g1["covariance"][0, 1] > 0.0
    assert g2["covariance"][0, 1] < 0.0


def test_component_weight_profile_partitions_the_pdf_mass_by_height():
    figure = report.component_weights_figure(_cube(), 0)
    assert [trace.name for trace in figure.data] == ["G1", "G2", "G3"]
    np.testing.assert_allclose(figure.data[0].x[1], 0.2)
    np.testing.assert_allclose(figure.data[1].x[1], 0.6)
    np.testing.assert_allclose(figure.data[2].x[1], 0.2)
    np.testing.assert_allclose(figure.data[0].y, (100.0, 200.0))


def test_scalar_component_geometry_uses_chi_covariance_and_rc_center():
    cube = _cube()
    chi = report._active_scalar_component(cube, 0, 1, 1, "chi")
    rc = report._active_scalar_component(cube, 0, 1, 3, "rc")
    assert chi["available"] and chi["covariance"][0, 1] > 0.0
    assert rc["available"] and rc["covariance"] is None
    np.testing.assert_allclose(rc["mean"], (0.8, 0.0008))


def test_whole_case_occupancy_comparison_has_both_center_metrics(monkeypatch):
    monkeypatch.setattr(report, "_manifest", lambda: {"cases": [{"name": "arm"}]})
    monkeypatch.setattr(
        report,
        "_whole_case_center_metrics",
        lambda _case, scheme: (1.0, 2.0, 3) if scheme == "speed" else (0.5, 2.2, 3),
    )
    figure = report.occupancy_comparison_figure()
    assert len(figure.data) == 4
    assert figure.data[0].name == "Crossing speed w"
    assert figure.data[2].name == "Occupancy (equal crossing)"
    assert figure.layout.yaxis.title.text == "RMSE [m/s]"
    assert figure.layout.yaxis2.title.text == "RMSE [g/kg]"


def test_moist_flux_reduction_excludes_dry_anomaly():
    result = report._ensemble(_cube(), 0, 1, "moist_flux")
    assert result["available"]
    assert result["count"] == 2
    assert result["mean"][1] > 0.012


def test_downward_ensemble_preserves_negative_w_and_positive_support():
    result = report._downward_ensemble(_cube(), 0, 1)
    assert result["available"]
    assert result["mean"][0] < 0.0
    assert result["count"] == 3
    assert result["neff"] == 36.0 / 14.0


def test_ellipse_uses_rt_on_x_and_w_on_y():
    mean = np.array((2.0, 0.012, 300.0))
    covariance = np.diag((0.25, 1.0e-6, 0.1))
    x, y = report._ellipse(mean, covariance)
    assert len(x) == len(y) == 121
    assert np.isclose(np.mean(x), mean[1], atol=2.0e-5)
    assert np.isclose(np.mean(y), mean[0], atol=2.0e-2)


def test_report_is_registered():
    assert report.REPORT.slug == "pdf9-parcel-ensemble"
    assert report.REPORT.register_callbacks is report.register_callbacks


def test_incomplete_cache_is_reported_before_rendering(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "CACHE_DIR", tmp_path)
    with pytest.raises(FileNotFoundError, match="Generate / refresh PDF-9 data"):
        report._validate_cache_files(
            {"cases": [{"name": "arm", "data_file": "arm.npz"}]}
        )


def test_tracks_include_every_launch_and_two_sigma_band():
    figure = report.tracks_figure(_cube(), 0, 1)
    names = {trace.name for trace in figure.data}
    assert "All 3 upward launch parcels" in names
    assert "All 3 downward launch parcels" in names
    assert "Downward crossings at selected altitude" in names
    assert "CLUBB grid mean r_t +/- 2 sigma" in names
    assert "Crossings at selected altitude" in names
    assert "SAM cloud-water-weighted center" in names
    assert "SAM clear moist updraft" in names
    assert figure.layout.xaxis2.title.text == "θₗ [K]"


def test_individual_crossings_are_drawn_as_distinct_footprints():
    cube = _cube()
    candidate = report._ensemble(cube, 0, 1, "speed")
    figure = report.go.Figure()
    report._add_parcel_footprints(figure, cube, 0, 1, candidate)
    names = {trace.name for trace in figure.data}
    assert "Individual parcel footprints (illustrative local width)" in names
    assert "Individual parcel crossing centers" in names


def test_individual_downward_crossings_are_drawn_as_distinct_footprints():
    cube = _cube()
    result = report._downward_ensemble(cube, 0, 1)
    figure = report.go.Figure()
    report._add_downward_parcel_footprints(figure, cube, 0, 1, result)
    names = {trace.name for trace in figure.data}
    assert "Downward parcel footprints (illustrative local width)" in names
    assert "Downward parcel crossing centers" in names


def test_downward_tracks_connect_launch_to_lower_destinations_in_order():
    cube = _cube()
    cube["down_parcel_status"][:] = 0
    cube["down_parcel_status"][0, 0, 2] = 1
    cube["down_parcel_w"][0, 0, 2] = -0.5
    cube["down_parcel_rt"][0, 0, 2] = 0.010
    cube["down_parcel_thl"][0, 0, 2] = 300.0
    figure = report.tracks_figure(cube, 0, 0)
    trace = next(
        item
        for item in figure.data
        if item.name == "All 1 downward launch parcels"
    )
    np.testing.assert_allclose(np.asarray(trace.y[:2], float), (120.0, 100.0))


def test_occupancy_center_marker_tracks_equal_crossing_mean():
    cube = _cube()
    figure = report.go.Figure()
    occupancy = report._add_occupancy_center(figure, cube, 0, 1)
    assert occupancy["available"]
    assert len(figure.data) == 1
    marker = figure.data[0]
    assert marker.name == "Occupancy-weighted candidate center"
    assert marker.marker.symbol == "diamond-open"
    np.testing.assert_allclose(marker.x, (12.0,))
    np.testing.assert_allclose(marker.y, (2.0,))


def test_recursive_candidate_reconstructs_fortran_mean_and_covariance():
    result = report._recursive_candidate(_cube(), 0, 1)
    assert result["available"]
    np.testing.assert_allclose(result["mean"], (0.5, 0.0115, 300.25))
    np.testing.assert_allclose(
        result["covariance"],
        np.array(
            (
                (8.0 / 3.0, 13.0e-3 / 12.0, -5.0 / 24.0),
                (13.0e-3 / 12.0, 2.9166666667e-6, 0.0012083333),
                (-5.0 / 24.0, 0.0012083333, 0.7291666667),
            )
        ),
    )


def test_downward_candidate_reconstructs_separate_fortran_geometry():
    result = report._stored_candidate(_cube(), 0, 1, "down")
    assert result["available"]
    np.testing.assert_allclose(result["mean"], (-1.0, 0.011, 300.5))
    assert result["covariance"][0, 1] < 0.0


def test_tracks_and_recursive_audit_show_zero_weight_g3_diagnosis():
    cube = _cube()
    tracks = report.tracks_figure(cube, 0, 1)
    assert "Fortran pooled G3 candidate" in {trace.name for trace in tracks.data}
    audit = report.recursive_candidate_figure(cube)
    assert len(audit.data) == 9
    assert "Distance-weighted bidirectional G3" in audit.layout.title.text
    assert "Branch/distance-weighted support" in {
        annotation.text for annotation in audit.layout.annotations
    }
    assert "Prior G3 used for downward launch" in {
        annotation.text for annotation in audit.layout.annotations
    }


def test_entrainment_environment_shows_directional_scalar_targets():
    figure = report.entrainment_environment_figure(_cube(), 0, 1)
    names = {trace.name for trace in figure.data}
    assert "Grid-mean environment" in names
    assert "Upward parcel target" in names
    assert "Downward parcel target" in names
    assert "PDF G3 weight" in names
    assert "Weight used by entrainment" in names
    np.testing.assert_allclose(figure.data[1].x[1], 11.8)
    np.testing.assert_allclose(figure.data[2].x[1], 10.6)


def test_downward_arrival_cap_audit_preserves_raw_and_capped_fields():
    cube = _cube()
    figure = report.downward_arrival_cap_figure(cube)
    assert len(figure.data) == 6
    np.testing.assert_allclose(figure.data[0].z[1, 0], -1.5)
    np.testing.assert_allclose(figure.data[1].z[1, 0], -1.0)
    np.testing.assert_allclose(figure.data[3].z[1, 0], 0.5)
    np.testing.assert_allclose(figure.data[4].z[1, 0], 0.75)
    np.testing.assert_allclose(figure.data[5].z[1, 0], -0.8)


def test_closure_profile_overlays_raw_and_reconstructed_values(monkeypatch):
    cube = _cube()
    cube.update(
        {
            "sam_time_seconds": np.array((60.0,)),
            "sam_height_m": np.array((100.0, 200.0)),
            "sam_mean": np.array((((0.0, 0.011, 299.5), (0.1, 0.012, 300.0)),)),
            "sam_covariance": np.tile(
                np.diag((0.25, 1.0e-6, 0.5)), (1, 2, 1, 1)
            ),
            "sam_targets": np.array((((1.0e-5, 1.0e-5, 1.0e-8, 1.0e-5), (2.0e-5, 2.0e-5, 2.0e-8, 2.0e-5)),)),
        }
    )
    monkeypatch.setattr(
        report,
        "load_snapshot_pressure_profile",
        lambda _run_dir, _seconds: (60, np.array((100.0, 200.0)), np.array((100000.0, 99000.0))),
    )

    figure = report.closure_profile_figure(
        cube, {"sam_run_dir": "/tmp/sam"}, 0, report.CLOSURE_PACKET_FRACTION
    )

    assert len(figure.data) == 8
    np.testing.assert_allclose(figure.data[0].y, cube["height_m"])
    np.testing.assert_allclose(figure.data[1].y, cube["height_m"])


def test_layout_remains_reloadable_when_cache_is_missing(monkeypatch):
    report._clear_data_caches()

    def unavailable():
        raise FileNotFoundError("test cache absent")

    monkeypatch.setattr(report, "_manifest", unavailable)
    layout = report.build_layout()
    rendered = repr(layout)
    assert report.REFRESH_ID in rendered
    assert report.CACHE_STATUS_ID in rendered
    assert report.OVERLAY_VARIABLE_ID in rendered
    assert report.COMPONENT_WEIGHTS_ID in rendered
    assert report.CAP_ID in rendered
    assert report.ERROR_ID not in rendered
    assert report.ANOMALY_ID not in rendered
    assert report.BAKEOFF_ID not in rendered
    assert report.OCCUPANCY_COMPARISON_ID not in rendered
    assert report.CONTINUITY_ID not in rendered
    assert report.ENTRAINMENT_ID not in rendered
