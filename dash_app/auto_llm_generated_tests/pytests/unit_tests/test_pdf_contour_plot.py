import numpy as np
from netCDF4 import Dataset

from dash_app.plot_tab.plot_types import pdf_contour_plot


def test_raw_3d_runs_are_discovered_by_case_name() -> None:
    assert (pdf_contour_plot.SHARED_3D_SAM_ROOT / "arm_3d").name == "arm_3d"
    assert (pdf_contour_plot.SHARED_3D_SAM_ROOT / "bomex_3d").name == "bomex_3d"
    assert (pdf_contour_plot.SHARED_3D_SAM_ROOT / "lba_3d").name == "lba_3d"


def _write_pdf_stats(path, *, three_gaussian=False):
    time_count, height_count = 2, 3
    with Dataset(path, "w") as dataset:
        dataset.createDimension("time", time_count)
        dataset.createDimension("zt", height_count)
        dataset.createDimension("col", 1)
        time = dataset.createVariable("time", "f8", ("time",))
        time.units = "seconds"
        time[:] = [60.0, 120.0]
        zt = dataset.createVariable("zt", "f8", ("zt",))
        zt.units = "m"
        zt[:] = [100.0, 500.0, 1000.0]
        col = dataset.createVariable("col", "f8", ("col",))
        col[:] = [0.0]

        fields = {
            "mixt_frac": 0.4,
            "w_1": 1.0,
            "w_2": -0.5,
            "varnce_w_1": 0.20,
            "varnce_w_2": 0.30,
            "chi_1": 6.0e-4,
            "chi_2": -3.0e-4,
            "stdev_chi_1": 2.0e-4,
            "stdev_chi_2": 3.0e-4,
            "corr_w_chi_1": 0.0,
            "corr_w_chi_2": 0.0,
            "rt_1": 0.012,
            "rt_2": 0.009,
            "varnce_rt_1": 2.0e-6,
            "varnce_rt_2": 3.0e-6,
            "corr_w_rt_1": 0.0,
            "corr_w_rt_2": 0.0,
            "thl_1": 300.0,
            "thl_2": 298.5,
            "varnce_thl_1": 0.20,
            "varnce_thl_2": 0.30,
            "corr_w_thl_1": -0.25,
            "corr_w_thl_2": 0.10,
            "corr_rt_thl_1": -0.35,
            "corr_rt_thl_2": -0.20,
            "crt_1": 1.0,
            "crt_2": 1.0,
            "cthl_1": 0.0,
            "cthl_2": 0.0,
            "wm_zt": 0.1,
            "rcm": 2.0e-4,
            "rc_1": 5.0e-4,
            "rc_2": 1.0e-5,
            "cloud_frac_1": 0.8,
            "cloud_frac_2": 0.1,
        }
        if three_gaussian:
            fields.update(
                {
                    "mixt_frac_3": 0.2,
                    "w_3": 2.0,
                    "varnce_w_3": 0.1,
                    "rt_3": 0.016,
                    "varnce_rt_3": 1.0e-6,
                    "varnce_thl_3": 0.1,
                    "corr_rt_thl_3": 0.0,
                    "corr_w_rt_3": 0.7,
                    "chi_3": 1.0e-3,
                    "stdev_chi_3": 1.5e-4,
                    "corr_w_chi_3": 0.6,
                    "crt_3": 1.0,
                    "cthl_3": 0.0,
                    "rc_3": 1.0e-3,
                    "cloud_frac_3": 1.0,
                }
            )
        for name, value in fields.items():
            variable = dataset.createVariable(name, "f8", ("time", "zt", "col"))
            variable[:] = np.full((time_count, height_count, 1), value)


def _case_data(path):
    names = set(pdf_contour_plot.COMMON_PDF_FIELDS)
    for projection in pdf_contour_plot.PROJECTIONS:
        names.update(pdf_contour_plot._projection_fields(projection))
    names.update(("wm_zt", "rcm", "rc_1", "rc_2", "cloud_frac_1", "cloud_frac_2"))
    return {
        "name": "synthetic",
        "files": [str(path)],
        "compare_mode": False,
        "profile_vars": [{"label": name, "value": name} for name in sorted(names)],
        "time_len": 2,
        "time_controls_physical": False,
        "default_time_start": 1,
        "default_time_duration": 2,
        "height_slider_min": 100.0,
        "height_slider_max": 1000.0,
        "height_step": 10.0,
        "default_height_range": [100.0, 1000.0],
        "benchmarks": {
            "available_sources": ["sam", "coamps"],
            "sam": {"default": "sam.nc"},
            "coamps": {"default": "coamps.nc"},
        },
    }


def _benchmark_profiles(field_names):
    values = {
        "wm": 0.15,
        "rtm": 0.011,
        "thlm": 298.0,
        "rcm": 2.5e-4,
        "wp2": 0.5,
        "rtp2": 2.5e-6,
        "thlp2": 0.25,
        "rcp2": 4.0e-8,
        "rtpthlp": -1.0e-4,
        "wprtp": 2.0e-4,
        "wpthlp": -3.0e-2,
        "wprcp": 2.0e-5,
        "p_in_Pa": 9.0e4,
    }
    return {
        name: {
            "z_values": np.array([100.0, 1000.0]),
            "profile": np.array([values[name], values[name]]),
            "z_units": "m",
            "units": "",
            "long_name": name,
        }
        for name in field_names
    }


def test_covariance_ellipse_uses_requested_rms_radius():
    x_values, y_values = pdf_contour_plot.covariance_ellipse_points(
        [3.0, -2.0],
        [[4.0, 0.0], [0.0, 9.0]],
        radius=2.0,
    )

    assert np.isclose(np.min(x_values), -1.0)
    assert np.isclose(np.max(x_values), 7.0)
    assert np.isclose(np.min(y_values), -8.0)
    assert np.isclose(np.max(y_values), 4.0)


def test_bivariate_density_integrates_to_one():
    x_values = np.linspace(-8.0, 8.0, 501)
    y_values = np.linspace(-8.0, 8.0, 501)
    x_grid, y_grid = np.meshgrid(x_values, y_values)
    density = pdf_contour_plot._bivariate_normal_density(
        x_grid,
        y_grid,
        mean_x=0.5,
        mean_y=-0.25,
        var_x=2.0,
        var_y=1.0,
        correlation=0.4,
    )

    integral = np.trapezoid(np.trapezoid(density, x_values, axis=1), y_values)
    assert np.isclose(integral, 1.0, atol=2.0e-5)


def test_enclosed_mass_thresholds_identify_highest_density_regions():
    values = np.linspace(-5.0, 5.0, 401)
    x_grid, y_grid = np.meshgrid(values, values)
    density = pdf_contour_plot._bivariate_normal_density(
        x_grid, y_grid, 0.0, 0.0, 1.0, 1.0, 0.0
    )

    thresholds = pdf_contour_plot.enclosed_mass_thresholds(density, values, values)

    assert len(thresholds) == len(pdf_contour_plot.ENCLOSED_MASSES) == 3
    assert thresholds[0] > thresholds[1] > thresholds[2]
    total = np.sum(density)
    for mass, threshold in zip(pdf_contour_plot.ENCLOSED_MASSES, thresholds):
        enclosed = np.sum(density[density >= threshold]) / total
        assert enclosed >= mass
        assert enclosed < mass + 0.02


def test_smoothed_raw_enclosure_keeps_raw_weight_targets():
    raw_weight = np.zeros((31, 31))
    raw_weight[5, 5] = 2.0
    raw_weight[9, 21] = 3.0
    raw_weight[22, 8] = 5.0
    raw_weight[24, 25] = 1.0

    smooth, thresholds = pdf_contour_plot.smoothed_raw_enclosure(raw_weight, 1.25)

    assert len(thresholds) == len(pdf_contour_plot.ENCLOSED_MASSES) == 3
    assert not np.allclose(smooth, raw_weight)
    total = np.sum(raw_weight)
    for mass, threshold in zip(pdf_contour_plot.ENCLOSED_MASSES, thresholds):
        enclosed = np.sum(raw_weight[smooth >= threshold]) / total
        assert enclosed >= mass
        assert enclosed < mass + np.max(raw_weight) / total + 1.0e-12


def test_component_ellipse_includes_between_time_mean_variance():
    component = pdf_contour_plot.aggregate_gaussian_component(
        weights=[0.5, 0.5],
        mean_x=[-1.0, 1.0],
        mean_y=[0.0, 0.0],
        var_x=[1.0, 1.0],
        var_y=[0.5, 0.5],
        correlation=[0.0, 0.0],
    )

    assert np.allclose(component["mean"], [0.0, 0.0])
    assert np.isclose(component["covariance"][0, 0], 2.0)
    assert np.isclose(component["covariance"][1, 1], 0.5)


def test_three_gaussian_plot_detects_component_and_uses_absolute_weights(tmp_path):
    stats_path = tmp_path / "three_gaussian_stats.nc"
    _write_pdf_stats(stats_path, three_gaussian=True)
    case_data = _case_data(stats_path)

    plot = pdf_contour_plot.PdfContourPlotType()
    bundle = plot._density_bundle(
        str(stats_path),
        "w_rt",
        500.0,
        case_data,
        {
            "time_range": [1, 2],
            "time_point": 1,
            "time_mode": "range",
            "selected_column": 0,
        },
    )

    assert len(bundle["components"]) == 3
    assert np.isclose(bundle["components"][0]["weight"], 0.32)
    assert np.isclose(bundle["components"][1]["weight"], 0.48)
    assert np.isclose(bundle["components"][2]["weight"], 0.20)
    assert "component_3" in bundle

    figure = plot.build_figure(
        {"var": "w_rt", "height": 500.0, "size": "large"},
        {
            "case_data": case_data,
            "time_range": [1, 2],
            "time_point": 1,
            "selected_column": 0,
            "enabled_benchmarks": [],
            "size": "large",
        },
    )
    component_3 = next(
        trace
        for trace in figure.data
        if trace.name == "CLUBB Gaussian 3 shape (2× RMS)"
    )
    assert component_3.line.color == "#22c55e"


def test_pdf_plot_removes_les_proxy_and_uses_one_signed_clubb_background(
    tmp_path, monkeypatch
):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)

    def _extract(_case_data, _source_name, field_names, *_time_args):
        return _benchmark_profiles(field_names)

    monkeypatch.setattr(
        pdf_contour_plot.benchmark_overlay,
        "extract_normalized_benchmark_profiles",
        _extract,
    )
    plot = pdf_contour_plot.PdfContourPlotType()
    figure = plot.build_figure(
        {"var": "w_chi", "height": 500.0, "size": "large"},
        {
            "case_data": case_data,
            "time_range": [1, 2],
            "time_point": 1,
            "selected_column": 0,
            "enabled_benchmarks": ["sam", "coamps"],
            "size": "large",
        },
    )

    heatmaps = [trace for trace in figure.data if trace.type == "heatmap"]
    assert [trace.name for trace in heatmaps] == [
        "CLUBB probability + signed w′r꜀′ cloud-water covariance",
    ]
    assert heatmaps[0].customdata.shape[-1] == 4
    assert "w′r꜀′ cloud-water covariance" in heatmaps[0].hovertemplate
    assert "normalized signed" not in heatmaps[0].hovertemplate
    assert not any(annotation.x == 0.99 for annotation in figure.layout.annotations)
    trace_names = {trace.name for trace in figure.data}
    expected_contour_label = pdf_contour_plot.CLUBB_ENCLOSED_WEIGHTED_LABEL.format(
        pdf_contour_plot.COLOR_SIGNALS["wprcp"]["label"]
    )
    assert expected_contour_label in trace_names
    assert "CLUBB Gaussian 1 shape (2× RMS)" in trace_names
    component_2 = next(
        trace for trace in figure.data if trace.name == "CLUBB Gaussian 2 shape (2× RMS)"
    )
    assert component_2.line.color == "#06b6d4"
    assert component_2.line.width == 3.0
    mass_contours = [
        trace
        for trace in figure.data
        if trace.type == "contour"
        and trace.name == expected_contour_label
    ]
    assert len(mass_contours) == len(pdf_contour_plot.ENCLOSED_MASSES)
    for trace in mass_contours:
        assert trace.colorscale[0][1] == pdf_contour_plot.MASS_CONTOUR_COLOR
        assert trace.colorscale[-1][1] == pdf_contour_plot.MASS_CONTOUR_COLOR
        assert trace.line.dash == "solid"
    assert [trace.opacity for trace in mass_contours] == [1.0, 0.72, 0.48]
    legend_entries = [
        trace
        for trace in figure.data
        if trace.name == expected_contour_label
        and bool(trace.showlegend)
    ]
    assert len(legend_entries) == 1


def test_matching_raw_sam_plane_replaces_clubb_heatmap(tmp_path, monkeypatch):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)
    sample_count = 200
    rt = np.linspace(0.008, 0.014, sample_count)
    w = np.linspace(-1.0, 2.0, sample_count)
    rc = np.maximum(rt - 0.011, 0.0)

    class Snapshot:
        samples = np.column_stack((w, rt, np.full(sample_count, 298.0)))
        chi_samples = rt - 0.011
        rc_samples = rc
        mean = np.mean(samples, axis=0)
        rcm = float(np.mean(rc))
        elapsed_minutes = 510.0
        height_m = 500.0

    monkeypatch.setattr(
        pdf_contour_plot,
        "_matching_raw_snapshot",
        lambda *_args, **_kwargs: Snapshot(),
    )
    figure = pdf_contour_plot.PdfContourPlotType().build_figure(
        {"var": "w_rt", "height": 500.0, "size": "large"},
        {
            "case_data": case_data,
            "time_range": 10.0,
            "time_point": 500.0 * 60.0,
            "selected_column": 0,
            "enabled_benchmarks": ["sam"],
            "size": "large",
        },
    )

    heatmaps = [trace for trace in figure.data if trace.type == "heatmap"]
    assert [trace.name for trace in heatmaps] == [
        "Raw SAM 3-D probability + signed w′r꜀′ cloud-water covariance"
    ]
    assert "w′r꜀′ cloud-water covariance" in heatmaps[0].hovertemplate
    assert not any(annotation.x == 0.99 for annotation in figure.layout.annotations)
    assert np.isclose(np.sum(np.asarray(heatmaps[0].customdata)[..., 0]), 1.0)
    raw_contours = [
        trace
        for trace in figure.data
        if trace.type == "contour"
        and trace.name.startswith("Raw SAM enclosed |w′r꜀′ cloud-water covariance")
    ]
    assert len(raw_contours) == len(pdf_contour_plot.ENCLOSED_MASSES)
    assert all(trace.line.color == pdf_contour_plot.RAW_MASS_CONTOUR_COLOR for trace in raw_contours)


def test_raw_color_signals_preserve_the_selected_plane_moment():
    sample_count = 240
    rt = np.linspace(0.008, 0.014, sample_count)
    w = np.linspace(-1.2, 2.4, sample_count)
    thl = 298.0 + 0.4 * np.sin(np.linspace(0.0, 3.0, sample_count))

    class Snapshot:
        source_file = "synthetic-raw-signal.nc"
        elapsed_seconds = 1800
        height_m = 500.0
        samples = np.column_stack((w, rt, thl))
        chi_samples = rt - 0.011
        rc_samples = np.maximum(rt - 0.011, 0.0)
        mean = np.mean(samples, axis=0)
        rcm = float(np.mean(rc_samples))
        elapsed_minutes = 30.0

    pdf_contour_plot._RAW_HISTOGRAM_CACHE.clear()
    x = np.linspace(7.5, 14.5, 41)
    y = np.linspace(-1.5, 2.7, 41)
    bundle = pdf_contour_plot._raw_histogram_bundle(
        Snapshot(), "w_rt", x, y, "wprtp2"
    )
    expected = np.mean(
        (w - np.mean(w)) * np.square(rt * 1000.0 - np.mean(rt * 1000.0))
    )

    assert bundle["color_signal"] == "wprtp2"
    assert "third-moment" in bundle["color_signal_label"]
    assert np.isclose(np.sum(bundle["total"]), 1.0)
    assert np.isclose(np.sum(bundle["signed_transport"]), expected)
    assert not np.allclose(bundle["enclosed_mass_field"], np.abs(bundle["signed_transport"]))
    assert len(bundle["enclosed_mass_thresholds"]) == len(pdf_contour_plot.ENCLOSED_MASSES)
    assert "|w′rₜ′² total-water third-moment contribution|" in bundle["enclosed_mass_label"]

    probability_only = pdf_contour_plot._raw_histogram_bundle(
        Snapshot(), "w_rt", x, y, "probability"
    )
    assert probability_only["probability_only"]
    assert not probability_only["color_available"]
    assert np.allclose(probability_only["signed_transport"], 0.0)
    assert not np.allclose(probability_only["enclosed_mass_field"], probability_only["total"])
    assert probability_only["enclosed_mass_label"] == pdf_contour_plot.RAW_ENCLOSED_PROBABILITY_LABEL


def test_clubb_enclosures_follow_the_selected_local_signal(tmp_path):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)
    context = {
        "time_range": [1, 2],
        "time_point": 1,
        "time_mode": "range",
        "selected_column": 0,
    }
    plot = pdf_contour_plot.PdfContourPlotType()

    probability = plot._density_bundle(
        str(stats_path), "w_rt", 500.0, case_data, context,
        color_signal="probability",
    )
    weighted = plot._density_bundle(
        str(stats_path), "w_rt", 500.0, case_data, context,
        color_signal="wprtp2",
    )

    assert probability["enclosed_mass_label"] == pdf_contour_plot.ENCLOSED_PROBABILITY_LABEL
    assert np.allclose(probability["enclosed_mass_field"], probability["total"])
    assert weighted["enclosed_mass_label"] == pdf_contour_plot.CLUBB_ENCLOSED_WEIGHTED_LABEL.format(
        pdf_contour_plot.COLOR_SIGNALS["wprtp2"]["label"]
    )
    assert not np.allclose(weighted["enclosed_mass_field"], weighted["total"])
    assert len(weighted["mass_thresholds"]) == len(pdf_contour_plot.ENCLOSED_MASSES)

    integrated = plot._density_bundle(
        str(stats_path), "rt_thl", 500.0, case_data, context,
        color_signal="wpthlp",
    )
    assert integrated["enclosed_mass_label"] == pdf_contour_plot.CLUBB_ENCLOSED_WEIGHTED_LABEL.format(
        pdf_contour_plot.COLOR_SIGNALS["wpthlp"]["label"]
    )
    assert not np.allclose(integrated["enclosed_mass_field"], integrated["total"])


def test_non_cloud_signal_is_probability_only_without_raw_sam(tmp_path):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)
    figure = pdf_contour_plot.PdfContourPlotType().build_figure(
        {"var": "w_rt", "color_signal": "wpthlp", "height": 500.0, "size": "large"},
        {
            "case_data": case_data,
            "time_range": [1, 2],
            "time_point": 1,
            "selected_column": 0,
            "enabled_benchmarks": [],
            "size": "large",
        },
    )
    heatmap = next(trace for trace in figure.data if trace.type == "heatmap")
    assert heatmap.name == "CLUBB probability"
    assert np.allclose(np.asarray(heatmap.customdata)[..., 1], 0.0)
    assert "probability density" in heatmap.hovertemplate
    assert "w′θₗ′ liquid-water-potential-temperature covariance contribution=unavailable" in heatmap.hovertemplate


def test_w_rc_plot_uses_finite_density_and_cloud_water_range(tmp_path, monkeypatch):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)
    monkeypatch.setattr(
        pdf_contour_plot.benchmark_overlay,
        "extract_normalized_benchmark_profiles",
        lambda _case, _source, fields, *_args: _benchmark_profiles(fields),
    )

    figure = pdf_contour_plot.PdfContourPlotType().build_figure(
        {"var": "w_rc", "height": 500.0, "size": "large"},
        {
            "case_data": case_data,
            "time_range": [1, 2],
            "time_point": 1,
            "selected_column": 0,
            "enabled_benchmarks": ["sam"],
            "size": "large",
        },
    )

    x_values = np.asarray(figure.data[0].x, dtype=float)
    density = np.asarray(figure.data[0].customdata, dtype=float)[..., 0]
    assert np.all(np.isfinite(density))
    assert x_values[0] == 0.0
    assert not any(annotation.x == 0.99 for annotation in figure.layout.annotations)
    assert not any("LES moment proxy" in (trace.name or "") for trace in figure.data)


def test_all_projections_map_nonzero_clubb_probability_and_signed_transport(
    tmp_path, monkeypatch
):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)
    monkeypatch.setattr(
        pdf_contour_plot.benchmark_overlay,
        "extract_normalized_benchmark_profiles",
        lambda _case, _source, fields, *_args: _benchmark_profiles(fields),
    )

    for projection in ("w_chi", "w_rt", "w_rc"):
        figure = pdf_contour_plot.PdfContourPlotType().build_figure(
            {"var": projection, "height": 500.0, "size": "large"},
            {
                "case_data": case_data,
                "time_range": [1, 2],
                "time_point": 1,
                "selected_column": 0,
                "enabled_benchmarks": ["sam"],
                "size": "large",
            },
        )
        heatmaps = [trace for trace in figure.data if trace.type == "heatmap"]
        assert len(heatmaps) == 1
        for trace in heatmaps:
            customdata = np.asarray(trace.customdata, dtype=float)
            assert customdata.shape[-1] == 4
            assert np.nanmax(customdata[..., 0]) > 0.0
            assert np.any(np.isfinite(customdata[..., 1]))


def test_transport_fields_separate_upward_downward_and_signed_covariance():
    density = np.ones((1, 4))
    x_grid = np.array([[2.0, 0.0, 2.0, 0.0]])
    w_grid = np.array([[1.0, 1.0, -1.0, -1.0]])

    upward, downward, signed = pdf_contour_plot._cloud_transport_fields(
        density,
        x_grid,
        w_grid,
        "w_rc",
        mean_w=0.0,
        mean_rc=1.0,
    )

    assert np.allclose(upward, [[2.0, 0.0, 0.0, 0.0]])
    assert np.allclose(downward, [[0.0, 0.0, 2.0, 0.0]])
    assert np.allclose(signed, [[1.0, -1.0, -1.0, 1.0]])


def test_pdf_transport_background_is_always_combined_probability_and_signed(
    tmp_path, monkeypatch
):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)
    monkeypatch.setattr(
        pdf_contour_plot.benchmark_overlay,
        "extract_normalized_benchmark_profiles",
        lambda _case, _source, fields, *_args: _benchmark_profiles(fields),
    )
    plot = pdf_contour_plot.PdfContourPlotType()
    context = {
        "case_data": case_data,
        "time_range": [1, 2],
        "time_point": 1,
        "selected_column": 0,
        "enabled_benchmarks": ["sam"],
        "size": "large",
    }

    figure = plot.build_figure(
        {
            "var": "w_chi",
            "height": 500.0,
            "transport_view": "downward",
            "size": "large",
        },
        context,
    )
    heatmaps = [trace for trace in figure.data if trace.type == "heatmap"]
    assert [trace.name for trace in heatmaps] == [
        "CLUBB probability + signed w′r꜀′ cloud-water covariance"
    ]
    assert not any(bool(trace.showscale) for trace in heatmaps)
    assert "w′r꜀′ cloud-water covariance" in heatmaps[0].hovertemplate
    assert not any(annotation.x == 0.99 for annotation in figure.layout.annotations)


def test_w_rc_view_uses_hover_only_color_details(tmp_path, monkeypatch):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)
    monkeypatch.setattr(
        pdf_contour_plot.benchmark_overlay,
        "extract_normalized_benchmark_profiles",
        lambda _case, _source, fields, *_args: _benchmark_profiles(fields),
    )

    figure = pdf_contour_plot.PdfContourPlotType().build_figure(
        {
            "var": "w_rc",
            "height": 500.0,
            "transport_view": "signed",
            "size": "large",
        },
        {
            "case_data": case_data,
            "time_range": [1, 2],
            "time_point": 1,
            "selected_column": 0,
            "enabled_benchmarks": ["sam"],
            "size": "large",
        },
    )

    assert not any(annotation.x == 0.99 for annotation in figure.layout.annotations)


def test_pdf_projection_options_require_parameters_in_single_or_compare_mode(tmp_path):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    plot = pdf_contour_plot.PdfContourPlotType()
    case_data = _case_data(stats_path)

    assert [option["value"] for option in plot.case_data_options(case_data)] == [
        "w_chi",
        "w_rt",
        "w_rc",
        "w_thl",
        "rt_thl",
        "chi_thl",
    ]
    case_data["compare_mode"] = True
    assert [option["value"] for option in plot.case_data_options(case_data)] == [
        "w_chi",
        "w_rt",
        "w_rc",
        "w_thl",
        "rt_thl",
        "chi_thl",
    ]


def test_thermodynamic_marginals_integrate_w_transport(tmp_path):
    stats_path = tmp_path / "synthetic_stats.nc"
    _write_pdf_stats(stats_path)
    case_data = _case_data(stats_path)
    context = {
        "case_data": case_data,
        "time_range": [1, 2],
        "time_point": 1,
        "selected_column": 0,
        "enabled_benchmarks": [],
        "size": "large",
    }

    for projection in ("w_thl", "rt_thl", "chi_thl"):
        figure = pdf_contour_plot.PdfContourPlotType().build_figure(
            {"var": projection, "height": 500.0, "size": "large"}, context
        )
        heatmap = next(trace for trace in figure.data if trace.type == "heatmap")
        assert heatmap.name == "CLUBB probability + signed w′r꜀′ cloud-water covariance"
        assert np.all(np.isfinite(np.asarray(heatmap.customdata)))
        assert pdf_contour_plot.PROJECTIONS[projection]["y_title"] in figure.layout.yaxis.title.text


def test_pdf_compare_mode_builds_shared_axis_directory_panels(tmp_path, monkeypatch):
    first_path = tmp_path / "first_stats.nc"
    second_path = tmp_path / "second_stats.nc"
    _write_pdf_stats(first_path)
    _write_pdf_stats(second_path)
    case_data = _case_data(first_path)
    case_data.update(
        {
            "files": [str(first_path), str(second_path)],
            "compare_mode": True,
            "source_labels": ["ADG1", "New PDF"],
        }
    )
    monkeypatch.setattr(
        pdf_contour_plot.benchmark_overlay,
        "extract_normalized_benchmark_profiles",
        lambda _case, _source, fields, *_args: _benchmark_profiles(fields),
    )

    figure = pdf_contour_plot.PdfContourPlotType().build_figure(
        {"var": "w_chi", "height": 500.0, "size": "large"},
        {
            "case_data": case_data,
            "time_range": [1, 2],
            "time_point": 1,
            "selected_column": 0,
            "enabled_benchmarks": ["sam"],
            "size": "large",
        },
    )

    density_traces = [trace for trace in figure.data if trace.type == "heatmap"]
    assert [trace.name for trace in density_traces] == [
        "ADG1 probability + signed w′r꜀′ cloud-water covariance",
        "New PDF probability + signed w′r꜀′ cloud-water covariance",
    ]
    assert list(figure.layout.xaxis.range) == list(figure.layout.xaxis2.range)
    assert list(figure.layout.yaxis.range) == list(figure.layout.yaxis2.range)
    assert not any(bool(trace.showscale) for trace in density_traces)


def test_pdf_contour_defers_while_live_stats_output_is_being_replaced(tmp_path):
    """An in-place rerun must not turn a transient invalid NetCDF into a 500."""
    stats_path = tmp_path / "synthetic_stats.nc"
    stats_path.write_bytes(b"")
    plot = pdf_contour_plot.PdfContourPlotType()
    case_data = _case_data(stats_path)

    paths, unreadable_output = plot._paths_for_projection([str(stats_path)], "w_chi")
    figure = plot.build_figure(
        {"var": "w_chi", "height": 500.0, "size": "large"},
        {"case_data": case_data, "enabled_benchmarks": [], "size": "large"},
    )

    assert paths == []
    assert unreadable_output
    assert "being replaced by a run" in figure.layout.annotations[0].text


def test_pdf_height_constraints_follow_the_global_height_window():
    minimums, maximums, heights = pdf_contour_plot.constrained_pdf_heights(
        [600.0, 1200.0], [500.0, 800.0, 1500.0]
    )

    assert minimums == [600.0, 600.0, 600.0]
    assert maximums == [1200.0, 1200.0, 1200.0]
    assert heights == [600.0, 800.0, 1200.0]
