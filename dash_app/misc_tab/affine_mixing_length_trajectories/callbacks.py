"""Connect saved-column controls to the coherent mixing-length figures.

Usage:
  ``register_callbacks(app)`` installs the six callbacks for the Affine
  Coherence Mixing Length page. The page layout uses the ``amlt-`` IDs.

Understanding:
  File selection sets record/column bounds and loads saved parameters.
  Record preparation truncates height, builds directional parcel fields,
  and caches their endpoint barriers. Changing only energy multiplier,
  aggregation, or directional weight reuses those matrices. File path,
  modification time/size, zero-based record/column, entrainment rate, launch
  fraction, and altitude ceiling identify a cached entry; at most 64 remain.
  A lock serializes access to this shared cache.

  Diagnostics reduce the selected record and a nearby time window into five
  figures and four metric cards. Expected input/file errors leave existing
  figures visible and return an error message. These callbacks read saved
  NetCDF data; they neither advance the model nor modify its files.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from threading import RLock

import numpy as np
from dash import Input, Output, State, ctx, html, no_update

from .analysis import (
    compute_inverse_square_mixing_length,
    compute_record_fields,
    compute_required_energy_barriers,
    discover_netcdf_file_records,
    inspect_dataset,
    load_dataset_record,
    profile_metrics,
)
from .layout import _case_options, _directory_options
from .figures import (
    coherence_length_profile_figure,
    coherence_matrix_figure,
    coherence_sum_figure,
    coherence_temporal_figure,
    parcel_matrix_figure,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
_MATRIX_CACHE_LOCK = RLock()


@dataclass(frozen=True)
class _PreparedRecord:
    """Keep a loaded/truncated record, parcel fields, barriers, and its original file
    entrainment rate.
    """

    record: object
    fields: object
    barriers: np.ndarray
    file_mu: float


def _truncate_record(record, altitude_max):
    """Return a record truncated through altitude_max [m].

    Keep at least three thermodynamic levels and one more momentum level,
    slicing state, PDF, and saved-length fields consistently without editing
    record. The retained top becomes the recomputed calculation's boundary.
    """

    level_count = int(np.searchsorted(record.z, float(altitude_max), side="right"))
    level_count = int(np.clip(level_count, 3, record.z.size))
    thermodynamic_slice = slice(0, level_count)
    momentum_slice = slice(0, level_count + 1)
    return replace(
        record,
        z=record.z[thermodynamic_slice],
        zm=record.zm[momentum_slice],
        thlm=record.thlm[thermodynamic_slice],
        rtm=record.rtm[thermodynamic_slice],
        mean_thlm=record.mean_thlm[thermodynamic_slice],
        mean_rtm=record.mean_rtm[thermodynamic_slice],
        mixt_frac=record.mixt_frac[thermodynamic_slice],
        thl_1=record.thl_1[thermodynamic_slice],
        thl_2=record.thl_2[thermodynamic_slice],
        rt_1=record.rt_1[thermodynamic_slice],
        rt_2=record.rt_2[thermodynamic_slice],
        thvm=record.thvm[thermodynamic_slice],
        pressure=record.pressure[thermodynamic_slice],
        exner=record.exner[thermodynamic_slice],
        thv_ds=record.thv_ds[thermodynamic_slice],
        em=record.em[momentum_slice],
        fortran_lscale=record.fortran_lscale[thermodynamic_slice],
        fortran_up=record.fortran_up[thermodynamic_slice],
        fortran_down=record.fortran_down[thermodynamic_slice],
    )


@lru_cache(maxsize=64)
def _cached_prepared_record(
    path,
    record_index,
    column,
    mu,
    pdf_launch_frac,
    altitude_max,
    file_mtime_ns,
    file_size,
):
    """Load and prepare one cache entry for zero-based record/column indices.

    Override mu [1/m], launch fraction, and altitude ceiling [m], then return
    _PreparedRecord. File modification time/size participate only in the cache
    key so an updated file gets a new entry.
    """

    del file_mtime_ns, file_size
    record = load_dataset_record(path, record_index, column)
    file_mu = float(record.mu)
    record = _truncate_record(replace(record, mu=mu), altitude_max)
    fields = compute_record_fields(record, pdf_launch_frac)
    return _PreparedRecord(
        record,
        fields,
        compute_required_energy_barriers(fields),
        file_mu,
    )


def _prepare_record(
    path,
    record_index,
    column,
    mu,
    pdf_launch_frac,
    altitude_max,
):
    """Return a cached prepared record, resolving the path and current file metadata
    under the cache lock.
    """
    resolved_path = Path(path).expanduser().resolve()
    file_status = resolved_path.stat()
    with _MATRIX_CACHE_LOCK:
        return _cached_prepared_record(
            str(resolved_path),
            int(record_index),
            int(column),
            float(mu),
            float(pdf_launch_frac),
            float(altitude_max),
            file_status.st_mtime_ns,
            file_status.st_size,
        )


def _clear_prepared_record_cache():
    """Discard all prepared records under the cache lock; return no value."""
    with _MATRIX_CACHE_LOCK:
        _cached_prepared_record.cache_clear()


def _custom_file_record(path):
    """Return the file-catalog entry for a resolved user-supplied NetCDF path."""
    path = Path(path).expanduser().resolve()
    return {
        "path": str(path),
        "key": str(path),
        "directory": str(path.parent),
        "root": str(path.parent),
        "relative_path": path.name,
        "label": str(path),
        "case": path.name.removesuffix("_stats.nc"),
    }


def _slider_marks(times):
    """Return record-index slider labels with times converted from seconds to minutes."""
    count = len(times)
    if count <= 1:
        return {0: "0"}
    indices = sorted(
        {0, count - 1, *(round(i * (count - 1) / 5) for i in range(1, 5))}
    )
    return {
        index: f"{index}\n({times[index] / 60.0:.0f}m)" for index in indices
    }


def _altitude_marks(altitude_top):
    """Return altitude slider labels from 500 m through altitude_top [m]."""
    values = {500.0, min(3500.0, altitude_top), altitude_top}
    for fraction in (0.25, 0.5, 0.75):
        values.add(500.0 * round(fraction * altitude_top / 500.0))
    return {
        value: f"{value / 1000.0:g} km" if value >= 1000.0 else f"{value:g} m"
        for value in sorted(values)
        if 500.0 <= value <= altitude_top
    }


def _metric_card(label, value, detail, color):
    """Return a Dash metric card with a label, value, explanatory detail, and value color."""
    return html.Div(
        [
            html.Div(label, className="mlt-metric-label"),
            html.Div(value, className="mlt-metric-value", style={"color": color}),
            html.Div(detail, className="mlt-metric-detail"),
        ],
        className="mlt-metric-card",
    )


def _temporal_series(
    path,
    selected_record,
    column,
    mu,
    pdf_launch_frac,
    tke_multiplier,
    up_weight,
    overwrite_mode,
    altitude_max,
):
    """Return recomputed/saved lengths and timestep changes near a selected record.

    Read one column and apply the supplied parcel controls and altitude ceiling
    [m] to each record. Use up to ten neighboring records on either side,
    except the retained 600--620 window when selected_record lies there and
    those records exist. Return times [min], the level with the largest saved
    profile jump [m], its length series [m], and profile RMS step changes [m].
    """

    metadata = inspect_dataset(path)
    record_count = metadata["record_count"]
    if record_count >= 621 and 600 <= selected_record <= 620:
        first_record, final_record = 600, 620
    else:
        first_record = max(0, selected_record - 10)
        final_record = min(record_count - 1, selected_record + 10)
    record_indices = np.arange(first_record, final_record + 1, dtype=int)

    times = []
    inverse_square_profiles = []
    fortran_profiles = []
    z = None
    for record_index in record_indices:
        prepared = _prepare_record(
            path,
            record_index,
            column,
            mu,
            pdf_launch_frac,
            altitude_max,
        )
        inverse_square = compute_inverse_square_mixing_length(
            prepared.fields,
            prepared.record.lmin,
            tke_multiplier,
            barriers=prepared.barriers,
            up_weight=up_weight,
            overwrite_mode=overwrite_mode,
        )
        times.append(prepared.record.time_seconds / 60.0)
        inverse_square_profiles.append(inverse_square.lscale)
        fortran_profiles.append(prepared.record.fortran_lscale)
        z = inverse_square.z

    time_minutes = np.asarray(times)
    inverse_square_profiles = np.asarray(inverse_square_profiles)
    fortran_profiles = np.asarray(fortran_profiles)
    if record_indices.size > 1:
        fortran_steps = np.diff(fortran_profiles, axis=0)
        sensitive_index = int(
            np.unravel_index(
                np.argmax(np.abs(fortran_steps)), fortran_steps.shape
            )[1]
        )
        inverse_square_roughness = np.sqrt(
            np.mean(np.diff(inverse_square_profiles, axis=0) ** 2, axis=1)
        )
        fortran_roughness = np.sqrt(np.mean(fortran_steps**2, axis=1))
    else:
        sensitive_index = int(
            np.argmax(np.abs(inverse_square_profiles[0] - fortran_profiles[0]))
        )
        inverse_square_roughness = np.empty(0)
        fortran_roughness = np.empty(0)
    selected_offset = int(
        np.clip(selected_record - first_record, 0, record_indices.size - 1)
    )
    return {
        "record_indices": record_indices,
        "time_minutes": time_minutes,
        "selected_time_minutes": float(time_minutes[selected_offset]),
        "sensitive_z": float(z[sensitive_index]),
        "inverse_square_at_level": inverse_square_profiles[:, sensitive_index],
        "fortran_at_level": fortran_profiles[:, sensitive_index],
        "inverse_square_roughness": inverse_square_roughness,
        "fortran_roughness": fortran_roughness,
    }


def register_callbacks(app):
    """Register file selection, parameter, label, and diagnostic callbacks on app.

    Return no value. The nested callbacks read/update the page's amlt-* controls
    and figures; parcel matrices are shared through the prepared-record cache.
    """
    @app.callback(
        Output("amlt-file-catalog", "data"),
        Output("amlt-output-directory", "options"),
        Output("amlt-output-directory", "value"),
        Input("amlt-load-path", "n_clicks"),
        Input("amlt-refresh-files", "n_clicks"),
        State("amlt-custom-path", "value"),
        State("amlt-output-directory", "value"),
        State("amlt-case", "value"),
        prevent_initial_call=True,
    )
    def update_affine_output_choices(
        _load_clicks,
        _refresh_clicks,
        custom_path,
        current_directory,
        current_path,
    ):
        """Return a refreshed file catalog and directory choices, preferring a loaded
        custom path or newest available file.
        """
        records = discover_netcdf_file_records(REPO_ROOT)
        preferred_path = current_path
        selected_directory = current_directory
        if ctx.triggered_id == "amlt-load-path" and custom_path:
            candidate = Path(custom_path).expanduser().resolve()
            inspect_dataset(candidate)
            if str(candidate) not in {str(record["path"]) for record in records}:
                records.append(_custom_file_record(candidate))
            preferred_path = str(candidate)
            selected_directory = str(candidate.parent)
        elif current_path not in {str(record["path"]) for record in records}:
            preferred_path = (
                str(
                    max(
                        records,
                        key=lambda record: float(record.get("modified", 0.0)),
                    )["path"]
                )
                if records
                else None
            )
        directories = {str(record["directory"]) for record in records}
        if selected_directory not in directories:
            selected_directory = (
                str(Path(preferred_path).parent)
                if preferred_path is not None
                else None
            )
        return (
            {"records": records, "preferred_path": preferred_path},
            _directory_options(records),
            selected_directory,
        )

    @app.callback(
        Output("amlt-case", "options"),
        Output("amlt-case", "value"),
        Input("amlt-file-catalog", "data"),
        Input("amlt-output-directory", "value"),
    )
    def update_affine_case_choices(catalog, directory):
        """Return case options in the chosen directory and the preferred or first available path."""
        catalog = catalog or {}
        records = catalog.get("records", [])
        options = _case_options(records, directory)
        available = {str(option["value"]) for option in options}
        preferred = catalog.get("preferred_path")
        selected = preferred if preferred in available else None
        if selected is None and options:
            selected = str(options[0]["value"])
        return options, selected

    @app.callback(
        Output("amlt-record", "max"),
        Output("amlt-record", "marks"),
        Output("amlt-record", "value"),
        Output("amlt-column", "options"),
        Output("amlt-column", "value"),
        Output("amlt-altitude-max", "max"),
        Output("amlt-altitude-max", "marks"),
        Output("amlt-altitude-max", "value"),
        Output("amlt-file-status", "children"),
        Input("amlt-case", "value"),
    )
    def update_affine_file_metadata(path):
        """Return record/column/altitude controls and a status message for path.

        Use zero-based indices, default to column zero and record 612 or the last
        available record, and cap the initial altitude at 3500 m. Missing/invalid
        files return placeholder controls and a readable error.
        """
        if not path:
            return (
                1,
                {0: "0", 1: "1"},
                0,
                [],
                None,
                20000.0,
                {3500.0: "3.5 km", 20000.0: "20 km"},
                3500.0,
                "Choose a file.",
            )
        try:
            metadata = inspect_dataset(path)
            first_record = load_dataset_record(path, 0, 0)
        except Exception as error:
            return (
                1,
                {0: "0", 1: "1"},
                0,
                [],
                None,
                20000.0,
                {3500.0: "3.5 km", 20000.0: "20 km"},
                3500.0,
                html.Span(str(error), className="mlt-error"),
            )
        columns = [
            {"label": str(index), "value": index}
            for index in range(metadata["column_count"])
        ]
        altitude_top = float(np.max(first_record.z))
        return (
            metadata["record_count"] - 1,
            _slider_marks(metadata["times"]),
            min(612, metadata["record_count"] - 1),
            columns,
            0,
            altitude_top,
            _altitude_marks(altitude_top),
            min(3500.0, altitude_top),
            f"{metadata['path']} · {metadata['record_count']} records · "
            f"{metadata['column_count']} column(s)",
        )

    @app.callback(
        Output("amlt-mu", "max"),
        Output("amlt-mu", "marks"),
        Output("amlt-mu", "value"),
        Output("amlt-tke-multiplier", "value"),
        Output("amlt-pdf-launch-frac", "value"),
        Output("amlt-up-weight", "value"),
        Input("amlt-case", "value"),
        Input("amlt-column", "value"),
    )
    def update_affine_mu_control(path, column):
        """Return saved parcel parameter values and entrainment slider bounds for the
        selected file/column.
        """
        if path is None or column is None:
            return 3.0e-3, {}, 1.0e-3, 1.0, 0.0, 0.5
        try:
            record = load_dataset_record(path, 0, column)
        except (IndexError, OSError, RuntimeError, TypeError, ValueError):
            return 3.0e-3, {}, 1.0e-3, 1.0, 0.0, 0.5
        slider_max = max(3.0e-3, 2.0 * record.mu)
        return (
            slider_max,
            {
                1.0e-5: "0.00001",
                record.mu: f"{record.mu:.4g} file",
                slider_max: f"{slider_max:.4g}",
            },
            record.mu,
            record.tke_multiplier,
            record.pdf_launch_frac,
            record.up_weight,
        )

    @app.callback(
        Output("amlt-mu-label", "children"),
        Output("amlt-altitude-max-label", "children"),
        Output("amlt-tke-multiplier-label", "children"),
        Output("amlt-pdf-launch-frac-label", "children"),
        Output("amlt-up-weight-label", "children"),
        Input("amlt-mu", "value"),
        Input("amlt-altitude-max", "value"),
        Input("amlt-tke-multiplier", "value"),
        Input("amlt-pdf-launch-frac", "value"),
        Input("amlt-up-weight", "value"),
    )
    def update_affine_slider_labels(
        mu,
        altitude_max,
        tke_multiplier,
        pdf_launch_frac,
        up_weight,
    ):
        """Return five readable parameter/altitude labels, using empty text for unset controls."""
        return (
            "" if mu is None else f"{float(mu):.6g} m⁻¹",
            "" if altitude_max is None else f"{float(altitude_max):.0f} m",
            ""
            if tke_multiplier is None
            else f"Ē/TKE={float(tke_multiplier):.3g}",
            ""
            if pdf_launch_frac is None
            else f"{float(pdf_launch_frac):.2f}",
            "" if up_weight is None else f"{float(up_weight):.2f}",
        )

    @app.callback(
        Output("amlt-parcel-matrix-figure", "figure"),
        Output("amlt-matrix-figure", "figure"),
        Output("amlt-profile-figure", "figure"),
        Output("amlt-temporal-figure", "figure"),
        Output("amlt-coherence-sum-figure", "figure"),
        Output("amlt-metrics", "children"),
        Output("amlt-time-label", "children"),
        Input("amlt-case", "value"),
        Input("amlt-record", "value"),
        Input("amlt-column", "value"),
        Input("amlt-mu", "value"),
        Input("amlt-altitude-max", "value"),
        Input("amlt-tke-multiplier", "value"),
        Input("amlt-pdf-launch-frac", "value"),
        Input("amlt-up-weight", "value"),
        Input("amlt-overwrite-mode", "value"),
    )
    def update_affine_diagnostics(
        path,
        record_index,
        column,
        mu,
        altitude_max,
        tke_multiplier,
        pdf_launch_frac,
        up_weight,
        overwrite_mode,
    ):
        """Return five figures, four metric cards, and the selected record label.

        Prepare/cache the chosen column, reduce its coherent travel, then compute
        a nearby time series. Missing controls preserve existing figures. Expected
        file/value errors return a message while preserving the previous figures.
        """
        if any(
            value is None
            for value in (
                path,
                record_index,
                column,
                mu,
                altitude_max,
                tke_multiplier,
                pdf_launch_frac,
                up_weight,
                overwrite_mode,
            )
        ):
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                [],
                "",
            )
        try:
            prepared = _prepare_record(
                path,
                record_index,
                column,
                mu,
                pdf_launch_frac,
                altitude_max,
            )
            result = compute_inverse_square_mixing_length(
                prepared.fields,
                prepared.record.lmin,
                tke_multiplier,
                barriers=prepared.barriers,
                up_weight=up_weight,
                overwrite_mode=overwrite_mode,
            )
            temporal = _temporal_series(
                path,
                int(record_index),
                int(column),
                float(mu),
                float(pdf_launch_frac),
                tke_multiplier,
                float(up_weight),
                overwrite_mode,
                float(altitude_max),
            )
        except (IndexError, OSError, RuntimeError, TypeError, ValueError) as error:
            message = f"Could not compute the affine coherence closure: {error}"
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                html.Span(message, className="mlt-error"),
                message,
            )

        record = prepared.record
        difference = profile_metrics(result.lscale, record.fortran_lscale)
        calculated_jump = (
            float(np.max(np.abs(np.diff(temporal["inverse_square_at_level"]))))
            if temporal["inverse_square_at_level"].size > 1
            else 0.0
        )
        fortran_jump = (
            float(np.max(np.abs(np.diff(temporal["fortran_at_level"]))))
            if temporal["fortran_at_level"].size > 1
            else 0.0
        )
        fractional_cells = int(
            np.count_nonzero(
                (result.coherence_matrix > 0.05)
                & (result.coherence_matrix < 0.95)
            )
        )
        metrics = [
            _metric_card(
                "Profile difference from saved Fortran",
                f"{difference['rmse']:.2f} m RMS",
                f"max |Δ| {difference['max_abs']:.2f} m · "
                f"r {difference['correlation']:.5f}",
                "#38bdf8",
            ),
            _metric_card(
                "Maximum timestep roughness",
                "inverse-square "
                f"{np.max(temporal['inverse_square_roughness'], initial=0.0):.1f} m",
                f"saved Fortran "
                f"{np.max(temporal['fortran_roughness'], initial=0.0):.1f} m",
                "#059669",
            ),
            _metric_card(
                f"Jump at z={temporal['sensitive_z']:.0f} m",
                f"inverse-square {calculated_jump:.1f} m",
                f"saved Fortran {fortran_jump:.1f} m",
                "#f59e0b",
            ),
            _metric_card(
                "Fractional coherence cells",
                f"{fractional_cells:,}",
                f"0.05 < q < 0.95 of {result.z.size**2:,} cells · "
                f"Nz={result.z.size}",
                "#8b5cf6",
            ),
        ]
        time_label = (
            f"record {record.record} · model time "
            f"{record.time_seconds / 60.0:.2f} min · column {record.column}"
        )
        return (
            parcel_matrix_figure(prepared.fields),
            coherence_matrix_figure(result),
            coherence_length_profile_figure(record, result),
            coherence_temporal_figure(temporal),
            coherence_sum_figure(result),
            metrics,
            time_label,
        )
