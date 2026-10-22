"""Check the saved-column affine/coherence calculation and its Dash interface.

Run through ``bash tests/run_pytests.sh -dash -k affine_mixing_length``.
Analytic/synthetic arrays check entrainment, endpoint barriers, inverse-square
layer integrals, aggregation weights, controls, and figure values. Synthetic
NetCDF files check migration of saved parameter names. The optional local ARM
checks use output/dash_default/arm_stats.nc at record 612 and skip when that
file is absent; they do not generate a new model run or prove Fortran parity.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from dash import Dash

from dash_app.misc_tab.affine_mixing_length_trajectories import SUBTAB
from dash_app.misc_tab.affine_mixing_length_trajectories import callbacks
from dash_app.misc_tab.affine_mixing_length_trajectories.analysis import (
    COHERENCE_AVERAGED_OVERWRITE,
    LEGACY_MAX_OVERWRITE,
    NO_OVERWRITE,
    _coherence_weighted_averages,
    affine_parcel_matrix,
    compute_inverse_square_mixing_length,
    compute_record_fields,
    compute_required_energy_barriers,
    inverse_square_coherence_fraction,
    load_dataset_record,
)
from dash_app.misc_tab.affine_mixing_length_trajectories.callbacks import (
    register_callbacks,
)
from dash_app.misc_tab.affine_mixing_length_trajectories.figures import (
    coherence_length_profile_figure,
    coherence_matrix_figure,
    coherence_sum_figure,
    coherence_temporal_figure,
    parcel_matrix_figure,
)
from dash_app.misc_tab.affine_mixing_length_trajectories.layout import (
    _case_options,
    _directory_options,
    build_layout,
)
from dash_app.misc_tab.mixing_length_trajectories.analysis import (
    _entrain_one_layer,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATION_FILE = REPO_ROOT / "output" / "dash_default" / "arm_stats.nc"


def _recursive_parcel_matrix(environment, z, mu):
    """Return a serial entrainment reference for every launch/target pair on ascending z."""
    reference = np.empty((z.size, z.size))
    for launch in range(z.size):
        reference[launch, launch] = environment[launch]
        parcel = environment[launch]
        for target in range(launch + 1, z.size):
            parcel = _entrain_one_layer(
                parcel,
                environment[target - 1],
                environment[target],
                z[target] - z[target - 1],
                mu,
            )
            reference[launch, target] = parcel
        parcel = environment[launch]
        for target in range(launch - 1, -1, -1):
            parcel = _entrain_one_layer(
                parcel,
                environment[target + 1],
                environment[target],
                z[target + 1] - z[target],
                mu,
            )
            reference[launch, target] = parcel
    return reference


def _walk(component):
    """Yield a Dash component and each nested component that supports Plotly serialization."""
    yield component
    children = getattr(component, "children", None)
    if children is None:
        return
    if not isinstance(children, (list, tuple)):
        children = [children]
    for child in children:
        if hasattr(child, "to_plotly_json"):
            yield from _walk(child)


def _registered_callback(app, name):
    """Return the undecorated registered callback with the requested function name."""
    return next(
        entry["callback"].__wrapped__
        for entry in app.callback_map.values()
        if entry["callback"].__name__ == name
    )


def test_affine_parcel_matrix_matches_explicit_layer_walk_both_directions():
    """Compare affine and serial entrainment on a nonuniform grid and require unchanged
    environmental diagonals.
    """
    z = np.array([12.0, 27.0, 61.0, 109.0, 188.0, 303.0])
    environment = np.array([299.2, 300.7, 298.9, 302.4, 301.1, 306.8])
    mu = 8.7e-4

    affine = affine_parcel_matrix(environment, z, mu)
    recursive = _recursive_parcel_matrix(environment, z, mu)

    np.testing.assert_allclose(affine, recursive, rtol=2.0e-15, atol=2.0e-12)
    np.testing.assert_array_equal(np.diag(affine), environment)


def test_record_fields_launch_without_velocity_moment_inputs():
    """Use directional PDF temperature/water departures and verify launches without
    velocity-moment inputs.
    """
    z = np.array([20.0, 60.0, 100.0, 140.0])
    thlm = np.array([299.0, 299.2, 299.4, 299.6])
    rtm = np.array([0.0100, 0.0098, 0.0096, 0.0094])
    record = SimpleNamespace(
        z=z,
        zm=np.array([0.0, 40.0, 80.0, 120.0, 160.0]),
        thvm=np.full(4, 301.0),
        thlm=thlm,
        rtm=rtm,
        em=np.full(5, 1.0),
        pressure=np.array([100000.0, 99500.0, 99000.0, 98500.0]),
        exner=np.full(4, 1.0),
        thv_ds=np.full(4, 300.0),
        mu=1.0e-3,
        pdf_launch_frac=0.5,
        mixt_frac=np.full(4, 0.01),
        thl_1=thlm + 1.0,
        thl_2=thlm - 2.0,
        rt_1=rtm + 0.001,
        rt_2=rtm - 0.002,
    )

    fields = compute_record_fields(record)

    expected_thl = affine_parcel_matrix(
        thlm,
        z,
        record.mu,
        thlm + 0.5,
        thlm - 1.0,
    )
    expected_rt = affine_parcel_matrix(
        rtm,
        z,
        record.mu,
        rtm + 0.0005,
        rtm - 0.001,
    )
    np.testing.assert_allclose(fields.parcel_thl_matrix, expected_thl)
    np.testing.assert_allclose(fields.parcel_rt_matrix, expected_rt)


def test_inverse_square_coherence_and_exact_first_layer_integral():
    """Check known survival fractions and a three-level analytic layer integral with the
    default aggregation.
    """
    np.testing.assert_allclose(
        inverse_square_coherence_fraction([0.0, 1.0, 9.0], 1.0),
        [1.0, 0.25, 0.01],
    )
    fields = SimpleNamespace(
        z=np.array([0.0, 40.0, 80.0]),
        tke=np.ones(3),
    )
    barriers = np.zeros((3, 3))
    barriers[0, 1:] = 9.0
    result = compute_inverse_square_mixing_length(
        fields, 20.0, barriers=barriers
    )

    assert result.raw_up[0] == pytest.approx(4.0)
    assert result.tke_multiplier == pytest.approx(1.0)
    assert result.overwrite_mode == COHERENCE_AVERAGED_OVERWRITE
    assert not hasattr(result, "pooled_up")
    assert not hasattr(result, "enveloped_up")


@pytest.mark.parametrize(
    "overwrite_mode",
    [
        NO_OVERWRITE,
        LEGACY_MAX_OVERWRITE,
        COHERENCE_AVERAGED_OVERWRITE,
    ],
)
def test_each_overwrite_choice_reuses_the_same_parcel_solution(overwrite_mode):
    """Reuse supplied three-level barriers for each aggregation and require finite
    lengths and unchanged barriers.
    """
    fields = SimpleNamespace(
        z=np.array([0.0, 40.0, 80.0]),
        tke=np.ones(3),
    )
    barriers = np.zeros((3, 3))
    barriers[0, 1:] = [1.0, 9.0]

    result = compute_inverse_square_mixing_length(
        fields,
        20.0,
        barriers=barriers,
        overwrite_mode=overwrite_mode,
    )

    assert result.overwrite_mode == overwrite_mode
    np.testing.assert_array_equal(result.barrier_matrix, barriers)
    assert np.all(np.isfinite(result.lscale))


def test_required_energy_barrier_uses_only_layer_endpoints():
    """Use work that recovers in a later layer and require the earlier endpoint barrier
    to persist.
    """
    fields = SimpleNamespace(
        z=np.array([0.0, 1.0, 2.0]),
        buoyancy_matrix=np.array(
            [
                [0.0, -2.0, 4.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ]
        ),
        launch_buoyancy_up=np.zeros(3),
        launch_buoyancy_down=np.zeros(3),
    )

    barriers = compute_required_energy_barriers(fields)

    np.testing.assert_allclose(barriers[0], [0.0, 1.0, 1.0])


def test_output_choices_group_cases_under_their_directory():
    """Build synthetic output choices and require cases to remain grouped and sorted
    within their directory.
    """
    records = [
        {
            "directory": "/tmp/output/a",
            "case": "arm",
            "path": "/tmp/output/a/arm_stats.nc",
        },
        {
            "directory": "/tmp/output/a",
            "case": "bomex",
            "path": "/tmp/output/a/bomex_stats.nc",
        },
        {
            "directory": "/tmp/output/b",
            "case": "rico",
            "path": "/tmp/output/b/rico_stats.nc",
        },
    ]

    directories = _directory_options(records)
    cases = _case_options(records, "/tmp/output/a")

    assert [option["value"] for option in directories] == [
        "/tmp/output/a",
        "/tmp/output/b",
    ]
    assert cases == [
        {"label": "arm", "value": "/tmp/output/a/arm_stats.nc"},
        {"label": "bomex", "value": "/tmp/output/a/bomex_stats.nc"},
    ]


def test_retained_average_weights_by_arriving_coherence():
    """Check a hand-computed pooled reach and remove remote arrivals to recover the local length."""
    z = np.array([0.0, 10.0, 20.0, 30.0])
    coherence = np.array([
        [1.0, 0.25, 0.0625, 0.0],
        [1.0, 1.0, 1.0, 0.0],
        [0.0625, 0.25, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])
    raw = np.ones(4)
    upward, downward = _coherence_weighted_averages(z, coherence, raw, raw)

    # At level 1, the remote population has coherence 0.25 and conditional
    # remaining reach 5 m; the local population has coherence 1 and reach 10 m.
    expected = (0.25 * 5.0 + 1.0 * 10.0) / 1.25
    assert upward[1] == pytest.approx(expected)
    assert downward[1] == pytest.approx(expected)

    coherence[0, 1:] = 0.0
    coherence[2, :2] = 0.0
    upward, downward = _coherence_weighted_averages(z, coherence, raw, raw)
    assert upward[1] == pytest.approx(10.0)
    assert downward[1] == pytest.approx(10.0)


def test_file_selection_uses_saved_parameters(monkeypatch):
    """Stub a loaded record and require the callback to initialize all four parcel
    controls from its saved values.
    """
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_callbacks(app)
    update = _registered_callback(app, "update_affine_mu_control")
    monkeypatch.setattr(
        callbacks,
        "load_dataset_record",
        lambda *_args: SimpleNamespace(
            mu=8.0e-4,
            tke_multiplier=0.9,
            pdf_launch_frac=1.0,
            up_weight=0.6,
        ),
    )

    values = update("/tmp/output/arm_stats.nc", 0)

    assert values[2:] == (8.0e-4, 0.9, 1.0, 0.6)


@pytest.mark.skipif(
    not VALIDATION_FILE.exists(), reason="Local ARM validation output is not present."
)
def test_dense_record_matrix_is_cached_once(monkeypatch):
    """Prepare the same optional ARM record twice and require one field/barrier build
    and the same cached result.
    """
    record = load_dataset_record(VALIDATION_FILE, 612, 0)
    build_counts = {"fields": 0, "barriers": 0}
    original_fields = callbacks.compute_record_fields
    original_barriers = callbacks.compute_required_energy_barriers

    def counted_fields(dataset_record, pdf_launch_frac):
        """Count parcel-field builds while delegating to the original calculation."""
        build_counts["fields"] += 1
        return original_fields(dataset_record, pdf_launch_frac)

    def counted_barriers(fields):
        """Count barrier builds while delegating to the original calculation."""
        build_counts["barriers"] += 1
        return original_barriers(fields)

    monkeypatch.setattr(callbacks, "compute_record_fields", counted_fields)
    monkeypatch.setattr(
        callbacks, "compute_required_energy_barriers", counted_barriers
    )
    callbacks._clear_prepared_record_cache()
    try:
        first = callbacks._prepare_record(
            VALIDATION_FILE, 612, 0, record.mu, record.pdf_launch_frac, 2500.0
        )
        second = callbacks._prepare_record(
            VALIDATION_FILE, 612, 0, record.mu, record.pdf_launch_frac, 2500.0
        )
    finally:
        callbacks._clear_prepared_record_cache()

    assert first is second
    assert first.fields.buoyancy_matrix.shape == (
        first.fields.z.size,
        first.fields.z.size,
    )
    assert build_counts == {"fields": 1, "barriers": 1}


@pytest.mark.skipif(
    not VALIDATION_FILE.exists(), reason="Local ARM validation output is not present."
)
def test_record_fields_feed_the_single_retained_reducer():
    """Use the optional ARM record to check matrix dimensions, energy diagonals, stop
    directions, and finite reduced lengths.
    """
    record = load_dataset_record(VALIDATION_FILE, 612, 0)
    fields = compute_record_fields(record)
    barriers = compute_required_energy_barriers(fields)
    result = compute_inverse_square_mixing_length(
        fields, record.lmin, barriers=barriers
    )

    nzt = fields.z.size
    assert fields.parcel_thl_matrix.shape == (nzt, nzt)
    assert fields.parcel_rt_matrix.shape == (nzt, nzt)
    assert fields.buoyancy_matrix.shape == (nzt, nzt)
    assert fields.energy_matrix.shape == (nzt, nzt)
    np.testing.assert_allclose(np.diag(fields.energy_matrix), fields.tke)
    assert np.all(fields.legacy_upward_stops[:-2] >= fields.z[:-2])
    assert np.all(fields.legacy_downward_stops[1:] <= fields.z[1:])
    np.testing.assert_array_equal(result.barrier_matrix, barriers)
    assert np.all(np.isfinite(result.lscale))


def test_retained_figures_show_barrier_coherence_profiles_and_time():
    """Use synthetic fields to verify physical hover data, directional sums, and
    calculated/saved figure traces.
    """
    z = np.array([20.0, 60.0, 100.0])
    fields = SimpleNamespace(
        z=z,
        tke=np.ones(3),
        energy_matrix=np.array(
            [[1.0, 0.5, -2.0], [0.8, 1.0, 0.4], [-1.0, 0.2, 1.0]]
        ),
        buoyancy_matrix=np.array(
            [[0.0, -0.01, 0.02], [0.01, 0.0, -0.02], [0.02, -0.01, 0.0]]
        ),
        legacy_upward_stops=np.array([60.0, np.nan, np.nan]),
        legacy_downward_stops=np.array([np.nan, 20.0, 60.0]),
    )
    barriers = np.zeros((3, 3))
    barriers[0, 1:] = [1.0, 4.0]
    result = compute_inverse_square_mixing_length(
        fields, 20.0, barriers=barriers
    )
    record = SimpleNamespace(
        fortran_up=result.lscale_up.copy(),
        fortran_down=result.lscale_down.copy(),
        fortran_lscale=result.lscale.copy(),
    )
    series = {
        "record_indices": np.array([0, 1, 2]),
        "time_minutes": np.array([0.0, 1.0, 2.0]),
        "selected_time_minutes": 1.0,
        "sensitive_z": 60.0,
        "inverse_square_at_level": np.array([10.0, 12.0, 11.0]),
        "fortran_at_level": np.array([9.0, 11.0, 10.0]),
        "inverse_square_roughness": np.array([1.0, 0.5]),
        "fortran_roughness": np.array([0.8, 0.4]),
    }

    parcel_matrix = parcel_matrix_figure(fields)
    matrix = coherence_matrix_figure(result)
    coherence_sums = coherence_sum_figure(result)
    profiles = coherence_length_profile_figure(record, result)
    temporal = coherence_temporal_figure(series)

    assert [trace.type for trace in parcel_matrix.data] == [
        "heatmap",
        "heatmap",
        "scattergl",
        "scattergl",
        "scattergl",
        "scattergl",
    ]
    np.testing.assert_allclose(
        np.asarray(parcel_matrix.data[0].customdata)[:, :, 1].astype(float),
        fields.energy_matrix.T,
    )
    np.testing.assert_allclose(
        np.asarray(parcel_matrix.data[1].customdata)[:, :, 1].astype(float),
        fields.buoyancy_matrix.T,
    )
    assert parcel_matrix.data[2].marker.color == "#000000"
    assert parcel_matrix.data[4].marker.color == "#064e3b"
    assert parcel_matrix.layout.legend.y < 0.0
    assert [trace.type for trace in matrix.data] == [
        "heatmap",
        "heatmap",
        "contour",
    ]
    np.testing.assert_allclose(
        np.asarray(matrix.data[0].customdata)[:, :, 1].astype(float),
        barriers.T,
    )
    expected_upward_sum = np.array(
        [
            np.sum(result.coherence_matrix[: target + 1, target])
            for target in range(z.size)
        ]
    )
    np.testing.assert_allclose(coherence_sums.data[0].x, expected_upward_sum)
    np.testing.assert_allclose(coherence_sums.data[0].y, z)
    expected_downward_sum = np.array(
        [
            np.sum(result.coherence_matrix[target:, target])
            for target in range(z.size)
        ]
    )
    np.testing.assert_allclose(coherence_sums.data[1].x, expected_downward_sum)
    np.testing.assert_allclose(coherence_sums.data[1].y, z)
    np.testing.assert_allclose(
        coherence_sums.data[2].x,
        np.sum(result.coherence_matrix, axis=0),
    )
    np.testing.assert_allclose(coherence_sums.data[2].y, z)
    assert len(coherence_sums.data) == 3
    assert {trace.legendgroup for trace in profiles.data} == {
        "inverse-square",
        "fortran",
    }
    assert len(profiles.data) == 6
    assert len(temporal.data) == 4


@pytest.mark.skipif(
    not VALIDATION_FILE.exists(), reason="Local ARM validation output is not present."
)
def test_diagnostic_callback_populates_the_reduced_panel_set():
    """Use the optional ARM record and require five figures, four metrics, and the
    selected-record label.
    """
    record = load_dataset_record(VALIDATION_FILE, 612, 0)
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_callbacks(app)
    update = _registered_callback(app, "update_affine_diagnostics")

    result = update(
        str(VALIDATION_FILE),
        612,
        0,
        record.mu,
        2500.0,
        record.tke_multiplier,
        record.pdf_launch_frac,
        record.up_weight,
        COHERENCE_AVERAGED_OVERWRITE,
    )

    assert len(result) == 7
    assert all(hasattr(figure, "data") for figure in result[:5])
    assert len(result[5]) == 4
    assert "record 612" in result[6]


def test_affine_tab_builds_only_the_retained_controls_and_panels():
    """Build the page and verify its control/graph IDs, callback registration, defaults,
    and unbounded nonnegative multiplier.
    """
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_callbacks(app)
    components = list(_walk(build_layout()))
    component_ids = {getattr(component, "id", None) for component in components}
    components_by_id = {
        component.id: component
        for component in components
        if getattr(component, "id", None) is not None
    }

    assert SUBTAB.title == "Affine Coherence Mixing Length"
    assert len(app.callback_map) == 6
    assert {
        "amlt-parcel-matrix-figure",
        "amlt-matrix-figure",
        "amlt-profile-figure",
        "amlt-temporal-figure",
        "amlt-coherence-sum-figure",
        "amlt-tke-multiplier",
        "amlt-pdf-launch-frac",
        "amlt-overwrite-mode",
        "amlt-up-weight",
        "amlt-altitude-max",
        "amlt-output-directory",
        "amlt-case",
        "amlt-file-catalog",
    } <= component_ids
    assert {
        "amlt-upward-figure",
        "amlt-downward-figure",
        "amlt-survival-matrix-figure",
    }.isdisjoint(component_ids)
    assert components_by_id["amlt-tke-multiplier"].value == 1.0
    assert components_by_id["amlt-tke-multiplier"].min == 0.0
    assert getattr(components_by_id["amlt-tke-multiplier"], "max", None) is None
    assert components_by_id["amlt-pdf-launch-frac"].value == 0.0
    assert "amlt-energy-spread" not in component_ids
    assert (
        components_by_id["amlt-overwrite-mode"].value
        == COHERENCE_AVERAGED_OVERWRITE
    )


@pytest.mark.parametrize("multiplier", [0.0, 1000.0])
def test_tke_multiplier_has_no_upper_bound_and_allows_zero(multiplier):
    """Reduce zero and large launch-energy multipliers and check finite lengths, zero
    travel, and numerical floors.
    """
    fields = SimpleNamespace(z=np.array([0.0, 40.0, 80.0]), tke=np.ones(3))
    result = compute_inverse_square_mixing_length(
        fields, 20.0, multiplier, barriers=np.zeros((3, 3))
    )
    assert np.all(np.isfinite(result.lscale))
    if multiplier == 0.0:
        np.testing.assert_array_equal(result.layer_integrals, 0.0)
        np.testing.assert_array_equal(result.raw_up, 0.1)
        np.testing.assert_array_equal(result.raw_down, 0.1)
    else:
        assert result.raw_up[0] == pytest.approx(40.0)


@pytest.mark.parametrize("multiplier", [-1.0, np.nan, np.inf])
def test_tke_multiplier_rejects_invalid_values(multiplier):
    """Require negative and nonfinite multipliers to fail before attempting parcel reduction."""
    with pytest.raises(ValueError, match="finite and nonnegative"):
        compute_inverse_square_mixing_length(None, 1.0, multiplier)


@pytest.mark.parametrize(
    "parameters, expected",
    [
        ({}, (1.0, 0.0, 0.5)),
        ({"Lscale_survival_fraction": 25.0 / 49.0,
          "Lscale_pdf_launch_frac": 0.7, "Lscale_up_strength": 0.3}, (2.5, 0.7, 0.3)),
        ({"Lscale_tke_multiplier": 1000.0, "Lscale_pdf_launch_frac": 0.0,
          "Lscale_up_weight": 0.8}, (1000.0, 0.0, 0.8)),
        ({"Lscale_tke_multiplier": 0.0, "Lscale_survival_fraction": 0.25,
          "Lscale_up_weight": 0.0, "Lscale_up_strength": 0.6}, (0.0, 0.0, 0.0)),
    ],
)
def test_saved_parameter_names_migrate_without_changing_energy(tmp_path, parameters, expected):
    """Write a tiny NetCDF parameter table and check old-name conversion, new-name
    precedence, and defaults.
    """
    import netCDF4 as nc
    from dash_app.misc_tab.mixing_length_trajectories.analysis import REQUIRED_VARIABLES

    path = tmp_path / "saved_stats.nc"
    with nc.Dataset(path, "w") as dataset:
        for name, size in (("zt", 3), ("zm", 4), ("time", 1), ("col", 1)):
            dataset.createDimension(name, size)
        for name in REQUIRED_VARIABLES:
            if name in {"zt", "zm", "time"}:
                dataset.createVariable(name, "f8", (name,))[:] = np.arange(
                    len(dataset.dimensions[name])
                )
            else:
                vertical = "zm" if name == "em" else "zt"
                dataset.createVariable(name, "f8", ("time", "col", vertical))[:] = 1.0
        if parameters:
            dataset.createDimension("param", len(parameters))
            dataset.createDimension("name_len", 40)
            dataset.createVariable("param_name", "S1", ("param", "name_len"))[:] = (
                nc.stringtochar(np.asarray(list(parameters), dtype="S40"))
            )
            dataset.createVariable("clubb_params", "f8", ("param", "col"))[:, 0] = (
                list(parameters.values())
            )
    record = load_dataset_record(path, 0)
    assert (record.tke_multiplier, record.pdf_launch_frac, record.up_weight) == pytest.approx(expected)
