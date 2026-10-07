"""In-memory CLUBB-vs-benchmark loss evaluation for tuning workflows.

This module owns the reusable loss path used by the tuner and standalone loss
executable. It initializes CLUBB once, prepares benchmark truth profiles from a
converted LES stats file, and returns per-window, per-variable, per-column
profile-comparison metrics. Loss settings come from ``&tuner_loss_nl``.

Benchmark profiles are interpolated onto the native CLUBB zt or zm source grid.
Benchmark records use the open/closed convention ``(window_start, window_end]``.
Successive windows continue from the previous model state. Nonfinite profiles
receive finite bad metrics so tuner ranking and output do not inherit NaN/Inf.

Adaptations: derived types are dataclasses, output arguments are return values,
indices are zero based, and driver globals live in the JAX state dictionary.
Host I/O and lifecycle surround JAX profile arithmetic and shared jitted model
kernels. Short candidate batches are padded to the initialized runtime shape.
See ``clubb_jax/TUNER_PORT.md`` for these boundaries and supported combinations.
"""
from dataclasses import dataclass, field
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
from netCDF4 import Dataset

from clubb_jax.src.Input_fields.namelist import _read_namelist_groups, read_namelist
from clubb_jax.src.clubb_case_initalization import (
    init_clubb_case,
    set_case_initial_conditions,
    clean_up_clubb,
)
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end
from clubb_jax.src.CLUBB_core.constants_clubb import g_per_kg, sec_per_day
from clubb_jax.src.CLUBB_core.interpolation import lin_interpolate_two_points
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.parameters_tunable import (
    PARAM_NAMES,
    PARAMETER_HARD_BOUNDS,
    PNAME_IDX,
)

loss_name_len = 64                 # Maximum stored variable-name length.
max_loss_variables = 32           # Maximum simultaneously scored variables.
invalid_loss_penalty = 1.0e30     # Finite penalty for nonfinite profiles.


@dataclass
class loss_field_type:
    """One CLUBB-vs-benchmark profile comparison."""

    clubb_var_name: str                           # Requested stats variable.
    benchmark_var_name: str                       # Paired benchmark variable.
    # Adaptation: grid-bank/slot pair replaces the Fortran stats registry id.
    stats_var_id: tuple = ()
    k_min: int = 0                               # Lower compared model level.
    k_max: int = 0                               # Upper compared model level.
    truth_profile: object = None                 # Benchmark by window/level.


@dataclass
class loss_request_type:
    """Full configuration for one loss evaluation."""

    num_variables: int = 0                        # Number of scored variables.
    time_average_range: tuple = (0, 0)            # Absolute averaging window [s].
    num_time_windows: int = 1                     # Number of equal subwindows.
    time_window_ranges: tuple = ()                # Absolute (start, end) [s].
    total_param_sets: int = 0                     # Configured parameter columns.
    dt_main_seconds: float = 0.0                  # Runtime main timestep [s].
    time_initial_seconds: float = 0.0             # Runtime initial time [s].
    altitude_comparison_range: tuple = (0.0, 0.0)  # Compared height window [m].
    les_stats_file: str = ''                      # Benchmark stats file.
    l_initialized: bool = False                   # Ready for repeated scoring.
    fields: list = field(default_factory=list)    # Metadata and prepared truth.
    # Adaptation: the Fortran driver owns these globals in clubb_driver.
    state: dict = field(default_factory=dict)


# Module-owned prepared request reused across manual loss calls.
active_request = loss_request_type()


def is_finite_core_value(value):
    """Return true only for finite CLUBB-core real values.

    Input: value, model or diagnostic real value.
    Return: finite predicate, avoiding NaN/Inf in tuner ranking and output.
    """
    return jnp.isfinite(value)


def set_invalid_field_metric_outputs():
    """Fill one field/column result with finite values that rank as very bad.

    Returns: scaled_rmse_value, correlation_value, std_ratio_value,
    centered_rmse_norm_value, bias_norm_value, in source output order.
    """
    return invalid_loss_penalty, 0.0, 0.0, invalid_loss_penalty, invalid_loss_penalty


@jax.jit
def calculate_taylor_metrics(
    model_profile, benchmark_profile,                         # Intent(in)
):
    """Compute Taylor diagnostics for a model and paired benchmark profile.

    Inputs: model_profile and benchmark_profile over compared levels.
    Returns: correlation, model/benchmark standard-deviation ratio,
    centered_rmse_norm and bias_norm, normalized by benchmark standard deviation.
    Adaptation: the last axis contains levels; leading axes vectorize columns.
    """
    if model_profile.shape[-1] != benchmark_profile.shape[-1]:
        stop_with_error('Taylor metric profiles must have matching sizes')
    num_levels = model_profile.shape[-1]
    if num_levels <= 0:
        stop_with_error('Taylor metric profiles must contain at least one level')

    model_mean = jnp.sum(model_profile, axis=-1) / num_levels
    benchmark_mean = jnp.sum(benchmark_profile, axis=-1) / num_levels

    # Adaptation: source level accumulations become whole-profile reductions.
    model_centered = model_profile - model_mean[..., None]
    benchmark_centered = benchmark_profile - benchmark_mean[..., None]
    model_centered_sumsq = jnp.sum(model_centered**2, axis=-1)
    benchmark_centered_sumsq = jnp.sum(benchmark_centered**2, axis=-1)
    covariance_sum = jnp.sum(model_centered * benchmark_centered, axis=-1)
    centered_diff_sumsq = jnp.sum((model_centered - benchmark_centered)**2, axis=-1)

    model_stddev = jnp.sqrt(model_centered_sumsq / num_levels)
    benchmark_stddev = jnp.sqrt(benchmark_centered_sumsq / num_levels)
    centered_rmse = jnp.sqrt(centered_diff_sumsq / num_levels)
    bias = model_mean - benchmark_mean

    # A flat benchmark has no vertical shape or variability amplitude to
    # compare. Treat those Taylor components as neutral and let centered
    # RMSE and bias carry any actual mismatch.
    # Adaptation: safe denominators keep inactive JAX branches finite.
    norm = jnp.where(benchmark_stddev > 0.0, benchmark_stddev, 1.0)
    correlation = jnp.where(
        benchmark_stddev <= 0.0,
        1.0,
        jnp.where(
            model_stddev > 0.0,
            jnp.clip(
                covariance_sum / jnp.where(
                    model_centered_sumsq * benchmark_centered_sumsq > 0.0,
                    jnp.sqrt(model_centered_sumsq * benchmark_centered_sumsq),
                    1.0,
                ),
                -1.0, 1.0,
            ),
            0.0,
        ),
    )
    std_ratio = jnp.where(benchmark_stddev <= 0.0, 1.0, model_stddev / norm)
    return correlation, std_ratio, centered_rmse / norm, bias / norm


def stop_with_error(message):
    """Report a fatal configuration error at the host boundary.

    Input: message. Adaptation: a Python exception replaces the source stop.
    """
    raise ValueError(message)


def init_loss_request(
    runfile, total_param_sets,                                # Intent(in)
    request,                                                 # Intent(out)
):
    """Read tuner_loss_nl and prepare reusable benchmark truth profiles.

    Inputs: aggregate runfile and total requested parameter-set columns.
    Output: request, populated in place; return requested_clubb_var_names in
    printed row order. Adaptation: request.state supplies driver-owned globals.
    """
    # Find the tuner-specific namelist inside the aggregate runfile.
    groups = _read_namelist_groups(runfile)
    if 'tuner_loss_nl' not in groups:
        stop_with_error(f'Missing &tuner_loss_nl in {runfile}')
    nml = groups['tuner_loss_nl']
    les_stats_file = str(nml.get('les_stats_file', '')).strip()
    altitude_comparison_range = tuple(
        nml.get('altitude_comparison_range', (0.0, 0.0))
    )
    time_average_range = tuple(nml.get('time_average_range', (0, 0)))
    num_time_windows = nml.get('num_time_windows', 1)
    if int(num_time_windows) != num_time_windows:
        stop_with_error('num_time_windows must be an integer')
    num_time_windows = int(num_time_windows)

    # Request validation
    # Reject malformed high-level loss input before touching CLUBB state.
    if not les_stats_file:
        stop_with_error('tuner_loss_nl requires les_stats_file')
    if (
        len(altitude_comparison_range) != 2
        or not np.isfinite(altitude_comparison_range).all()
        or not 0 <= altitude_comparison_range[0] <= altitude_comparison_range[1]
    ):
        stop_with_error(
            'tuner_loss_nl requires ordered non-negative altitude_comparison_range'
        )
    if (
        len(time_average_range) != 2
        or not np.isfinite(time_average_range).all()
        or not 0 <= time_average_range[0] < time_average_range[1]
    ):
        stop_with_error(
            'tuner_loss_nl requires an ordered non-negative time_average_range'
        )
    if any(int(value) != value for value in time_average_range):
        stop_with_error('time_average_range must use integer seconds')
    if num_time_windows < 1:
        stop_with_error('tuner_loss_nl requires num_time_windows >= 1')
    if (time_average_range[1] - time_average_range[0]) % num_time_windows:
        stop_with_error('time_average_range must divide evenly by num_time_windows')

    # Collapse the fixed-size namelist arrays down to the active requests.
    # Adaptation: the minimal namelist parser retains indexed keys; f90nml
    # returns lists. Normalize the fixed-size source arrays at this I/O boundary.
    names = {}
    for key in ('clubb_var_names', 'benchmark_var_name'):
        raw = nml.get(key, [])
        if isinstance(raw, str):
            raw = [raw]
        first = (getattr(nml, 'start_index', {}).get(key, [1])[0] or 1) - 1
        if first < 0:
            stop_with_error('Loss variable array indices must be one based')
        entries = [''] * first + list(raw)
        for nml_key, value in nml.items():
            if nml_key.startswith(key + '('):
                index = int(nml_key[len(key) + 1:-1]) - 1
                values = value if isinstance(value, list) else [value]
                if not 0 <= index or index + len(values) > max_loss_variables:
                    stop_with_error('Too many loss variables')
                entries.extend([''] * max(0, index + len(values) - len(entries)))
                entries[index:index + len(values)] = values
        if len(entries) > max_loss_variables:
            stop_with_error('Too many loss variables')
        names[key] = [str(value or '').strip()[:loss_name_len] for value in entries]
    clubb_var_names = names['clubb_var_names']
    benchmark_var_name = names['benchmark_var_name']
    nvars = 0
    l_seen_blank = False
    for i in range(max(len(clubb_var_names), len(benchmark_var_name))):
        name = clubb_var_names[i] if i < len(clubb_var_names) else ''
        benchmark_name = benchmark_var_name[i] if i < len(benchmark_var_name) else ''
        if not name:
            l_seen_blank = True
            if benchmark_name:
                stop_with_error('benchmark_var_name contains an unpaired entry')
            continue
        if l_seen_blank:
            stop_with_error('clubb_var_names contains entries after the first blank entry')
        if not benchmark_name:
            stop_with_error('Each CLUBB variable requires a paired benchmark_var_name')
        if name in clubb_var_names[:i]:
            stop_with_error(f'Duplicate CLUBB variable requested: {name}')
        nvars += 1

    if nvars == 0:
        stop_with_error('clubb_var_names must contain at least one CLUBB variable')

    # Store only the active request entries.
    request.les_stats_file = str(Path(les_stats_file).resolve())
    request.altitude_comparison_range = altitude_comparison_range
    request.time_average_range = time_average_range
    request.num_time_windows = num_time_windows
    request.num_variables = nvars
    request.total_param_sets = total_param_sets
    request.l_initialized = False
    window_width = (time_average_range[1] - time_average_range[0]) // num_time_windows
    request.time_window_ranges = tuple(
        (
            time_average_range[0] + i * window_width,
            time_average_range[0] + (i + 1) * window_width,
        )
        for i in range(num_time_windows)
    )

    # Store each requested comparison in one field descriptor.
    request.fields = [
        loss_field_type(clubb_var_names[i], benchmark_var_name[i])
        for i in range(nvars)
    ]

    # CLUBB-dependent preparation
    # Snapshot the configured stats registry used for name lookup and scoring.
    request.dt_main_seconds = request.state['dt_main']
    request.time_initial_seconds = request.state['time_initial']
    writer = request.state['stats_writer']
    if writer is None:
        stop_with_error('Loss evaluation requires enabled statistics')
    # Adaptation: JAX exposes the active accumulation bank, rather than Fortran's
    # averaged stats snapshot. Require one output window per scored subwindow.
    if (
        writer.stats_tstart != time_average_range[0]
        or writer.stats_tend != time_average_range[1]
        or writer.stats_nout * writer.dt_main != window_width
    ):
        stop_with_error(
            'Loss windows must match stats_tstart, stats_tend and stats_tout'
        )
    if time_average_range[1] > request.state['time_final']:
        stop_with_error('Loss window ends after the model run')

    # Finish the CLUBB-dependent request preparation before batching starts.
    prepare_loss_request_for_scoring(request, writer)

    # Return row labels in request order.
    requested_clubb_var_names = [item.clubb_var_name for item in request.fields]
    request.l_initialized = True
    return requested_clubb_var_names


def init_clubb_loss(
    runfile,                                                 # Intent(in)
    return_default_params=False,                             # Optional output selector
):
    """Initialize CLUBB once and prepare the reusable loss request.

    Input: aggregated runfile. Return: clubb_var_names and, when requested,
    clubb_params_all. Adaptation: the shared Python backend API selects the
    source's optional parameter output with return_default_params; runtime
    errors live in state['err_info'] and fatal setup errors raise exceptions.
    """
    global active_request

    # Phase 1: initialize CLUBB and prepare the loss request.
    if active_request.l_initialized:
        stop_with_error('init_clubb_loss called while a prepared request is still active')
    # TODO: define post-restart loss windows; the source loss loop starts at 1.
    # Reject this combination before standalone initialization restores a file
    # or creates output. Normal standalone restarts are supported independently.
    if read_namelist(runfile).get('l_restart', False):
        stop_with_error('JAX loss evaluation does not support l_restart')
    state = init_clubb_case(runfile)
    active_request = loss_request_type(state=state)
    try:
        names = init_loss_request(runfile, state['total_param_sets'], active_request)
    except Exception:
        clean_up_clubb(state)
        active_request = loss_request_type()
        raise
    if return_default_params:
        return names, state['clubb_params_all']
    return names


def get_loss_time_window_count():
    """Return the time-window count for the active prepared loss request."""
    return active_request.num_time_windows if active_request.l_initialized else 1


def clubb_get_loss_for_params(clubb_params_all):
    """Reuse initialized CLUBB/loss state to score parameter columns.

    Input: clubb_params_all, the full (parameter column, parameter) matrix.
    Returns: scaled_rmse, correlation, std_ratio, centered_rmse_norm and
    bias_norm, shaped (time window, requested variable, parameter column).
    Runtime errors remain in the driver state and penalize affected columns.

    Adaptation: in-memory reruns may vary the candidate count. Pad the final
    batch and discard spare columns, preserving one compiled runtime shape.
    NetCDF output retains the source's configured total-column count.
    """
    if not active_request.l_initialized:
        stop_with_error('clubb_get_loss_for_params requires init_clubb_loss first')
    request = active_request
    state = request.state
    clubb_params_all = jnp.asarray(clubb_params_all, dtype=state['clubb_params'].dtype)
    if (
        clubb_params_all.ndim != 2
        or clubb_params_all.shape[1] != len(PARAM_NAMES)
        or not clubb_params_all.shape[0]
    ):
        stop_with_error(
            f'clubb_params_all must have shape (candidates, {len(PARAM_NAMES)})'
        )

    # Adaptation: variable candidate counts fit the reusable in-memory bank,
    # but cannot change the initialized NetCDF column dimension.
    if (
        state['stats_output_path']
        and clubb_params_all.shape[0] != request.total_param_sets
    ):
        stop_with_error(
            'Variable candidate counts require in-memory statistics; '
            'set stats_output_filename=""'
        )

    runtime_batch_size = state['ngrdcol']
    batch_metrics = []

    # Run each parameter batch through the full case.
    for batch_start in range(0, clubb_params_all.shape[0], runtime_batch_size):
        params = clubb_params_all[batch_start:batch_start + runtime_batch_size]
        num_columns = params.shape[0]
        params = jnp.concatenate((
            params,
            jnp.repeat(params[-1:], runtime_batch_size - num_columns, axis=0),
        ))

        # Adaptation: validate candidates before the source's host lmin fatal
        # check. Substitute defaults for invalid columns during advancement,
        # then return their finite penalties without aborting healthy neighbors.
        valid = jnp.all(
            jnp.isfinite(params)
            & (params >= PARAMETER_HARD_BOUNDS[0])
            & (params <= PARAMETER_HARD_BOUNDS[1]),
            axis=1,
        )
        for left, right in (
            ('C6rt', 'C6thl'),
            ('C6rtb', 'C6thlb'),
            ('C6rtc', 'C6thlc'),
            ('C6rt_Lscale0', 'C6thl_Lscale0'),
        ):
            a, b = params[:, PNAME_IDX[left]], params[:, PNAME_IDX[right]]
            valid &= (
                jnp.abs(a - b)
                <= jnp.abs(a + b) * 0.5 * np.finfo(np.dtype(params.dtype)).eps
            )
        valid &= params[:, PNAME_IDX['lmin_coef']] * 40.0 >= 1.0
        params = jnp.where(valid[:, None], params, state['_initial_state']['clubb_params'])

        # Reset the state and activate the current parameter batch.
        if clubb_params_all.shape[0] == request.total_param_sets:
            set_case_initial_conditions(
                state,                                                        # InOut
                params, batch_start // runtime_batch_size + 1,                 # Optional in
            )
        else:
            set_case_initial_conditions(state, params)

        itime_start = 1
        window_metrics = []
        for window_idx, (_, end) in enumerate(request.time_window_ranges):
            itime_end = (end - request.time_initial_seconds) / request.dt_main_seconds
            # Fortran NINT rounds half-integers away from zero.
            itime_end = int(
                np.floor(itime_end + 0.5) if itime_end >= 0.0
                else np.ceil(itime_end - 0.5)
            )

            # Advance through this subwindow for the current batch.
            advance_clubb_to_end(
                state, False,                                                # InOut / in
                itime_start=itime_start, itime_end=itime_end,                  # Optional in
            )

            # Snapshot and score this completed stats subwindow.
            window_metrics.append(calculate_field_loss(
                request, state['_jax_stats'], window_idx,                     # In
            ))
            itime_start = itime_end + 1

        metrics = tuple(jnp.stack([values[i] for values in window_metrics]) for i in range(5))

        # Runtime errors remain column-local at debug level -1. Reject only
        # affected candidates after every healthy neighbor has finished once.
        valid &= ~state['err_info'].fatal_mask()
        metrics = tuple(
            jnp.where(valid[None, None, :], values, penalty)[..., :num_columns]
            for values, penalty in zip(metrics, set_invalid_field_metric_outputs())
        )
        batch_metrics.append(metrics)

    return tuple(
        jnp.concatenate([values[i] for values in batch_metrics], axis=2)
        for i in range(5)
    )


def finalize_clubb_loss():
    """Release CLUBB state and clear the prepared module-local loss request."""
    global active_request
    if active_request.state:
        clean_up_clubb(active_request.state)
    active_request = loss_request_type()


def clubb_get_loss(runfile):
    """Initialize, run all default parameter batches, and clean up.

    Input: aggregate runfile. Return: clubb_var_names and the five metric arrays,
    in source output order. Fatal runtime setup errors raise exceptions.
    """
    # Phase 1: initialize CLUBB and prepare the reusable loss request.
    names, params = init_clubb_loss(runfile, True)
    try:
        # Phase 2: run the default parameter matrix through the reusable loss path.
        return (names, *clubb_get_loss_for_params(params))
    finally:
        # Phase 3: release CLUBB state after the one-shot evaluation.
        finalize_clubb_loss()


def prepare_loss_request_for_scoring(
    request,                                                 # Intent(inout)
    stats_snapshot,                                          # Intent(in)
):
    """Resolve stats bindings, model levels and benchmark truth profiles.

    Input/output: request, enriched with bindings, height bounds and truth.
    Input: stats_snapshot, the configured native-grid stats registry.
    Adaptation: NetCDF metadata/indexing stay on the host; profile arithmetic
    uses JAX arrays and the public interpolation kernel.
    """
    # Use the native CLUBB grids for height-range resolution and interpolation.
    runtime_zt, runtime_zm = stats_snapshot.get_source_grid()
    banks = JaxStats.from_layout(
        stats_snapshot.get_jax_layout(),
        ncol=request.state['ngrdcol'],
    )

    # Benchmark time setup
    # Open the benchmark dataset once for the whole request.
    with Dataset(request.les_stats_file) as benchmark_file:
        # nf90_get_var reads raw values; Python must not mask valid zero fill
        # values or apply scale_factor/add_offset attributes implicitly.
        benchmark_file.set_auto_mask(False)
        benchmark_file.set_auto_scale(False)

        # input_netcdf.F90 resolves coordinate names from the first field.
        first_var = benchmark_file[request.fields[0].benchmark_var_name]
        time_name = next(
            (name for name in first_var.dimensions if name in ('T', 't', 'time')),
            None,
        )
        z_name = next(
            (name for name in first_var.dimensions if name in (
                'Z', 'z', 'zt', 'zm', 'altitude', 'height', 'lev',
                'lh_zt', 'rad_zt', 'rad_zm',
            )),
            None,
        )
        if time_name is None or z_name is None:
            stop_with_error('Benchmark variable requires time and altitude coordinates')

        # Convert the benchmark time coordinate to seconds.
        time_units = benchmark_file[time_name].units.split()[0]
        multiplier = {'hours': 3600.0, 'minutes': 60.0, 'seconds': 1.0}.get(time_units)
        if multiplier is None:
            stop_with_error(f'Benchmark time units are unsupported: {time_units}')

        raw_time_values = np.asarray(benchmark_file[time_name][:]).reshape(-1)
        benchmark_time_values_seconds = raw_time_values * multiplier
        if (
            np.any(np.diff(benchmark_time_values_seconds) <= 0)
            or not np.isfinite(benchmark_time_values_seconds).all()
        ):
            stop_with_error('Benchmark times must be finite and strictly ascending')
        z = np.asarray(benchmark_file[z_name][:]).reshape(-1)
        if np.any(np.diff(z) <= 0) or not np.isfinite(z).all():
            stop_with_error('Benchmark altitudes must be finite and strictly ascending')

        z_min, z_max = request.altitude_comparison_range

        # Per-field setup
        # Prepare each requested field once before batching starts.
        for field_idx, field in enumerate(request.fields):
            # Resolve the requested CLUBB name to one stats field and native grid.
            if field.clubb_var_name not in stats_snapshot.registry:
                stop_with_error(
                    f'Requested stats variable was not found: {field.clubb_var_name}'
                )
            field.stats_var_id = banks.name_to_slot[field.clubb_var_name]
            grid_name = stats_snapshot.registry[field.clubb_var_name][0]

            # Pick the native CLUBB grid used by this field.
            if grid_name not in ('zt', 'zm'):
                stop_with_error(
                    f'Requested CLUBB variable must live on zt or zm: {field.clubb_var_name}'
                )
            active_grid = (runtime_zt if grid_name == 'zt' else runtime_zm)[0]

            # Convert the requested height range into model levels.
            levels = np.flatnonzero((active_grid >= z_min) & (active_grid <= z_max))
            if not levels.size:
                stop_with_error('Requested height window does not include any model levels')
            field.k_min, field.k_max = int(levels[0]), int(levels[-1])

            if active_grid[field.k_min] < z[0] or active_grid[field.k_max] > z[-1]:
                stop_with_error(
                    'Requested height window is outside benchmark domain: '
                    + field.benchmark_var_name
                )

            # Build the vertical interpolation stencil once for this field.
            # Adaptation: only compared levels are stored in the truth profile.
            compared_grid = active_grid[field.k_min:field.k_max + 1]
            upper_idx = np.searchsorted(z, compared_grid).clip(0, len(z) - 1)
            lower_idx = np.maximum(upper_idx - 1, 0)
            exact_upper = np.abs(z[upper_idx] - compared_grid) < 1.e-6
            exact_lower = np.abs(z[lower_idx] - compared_grid) < 1.e-6
            upper_idx = np.where(exact_lower, lower_idx, upper_idx)
            lower_idx = np.where(exact_upper, upper_idx, lower_idx)

            # Adaptation: inactive exact-level interpolation uses a safe upper
            # height. The source copies those levels without dividing by zero.
            height_high = np.where(
                upper_idx == lower_idx, z[lower_idx] + 1.0, z[upper_idx]
            )

            var = benchmark_file[field.benchmark_var_name]
            if (
                var.ndim not in (3, 4)
                or np.dtype(var.dtype).kind != 'f'
                or np.dtype(var.dtype).itemsize not in (4, 8)
            ):
                stop_with_error(
                    "input_netcdf.get_var: The netCDF data doesn't conform to "
                    'expected precision, shape, or dimensions'
                )
            if time_name not in var.dimensions or z_name not in var.dimensions:
                stop_with_error(
                    'Benchmark variable requires time/altitude dimensions: '
                    + field.benchmark_var_name
                )
            if any(name not in (
                time_name, z_name, 'col', 'column',
                'X', 'x', 'longitude', 'lon', 'Y', 'y', 'latitude', 'lat',
            ) for name in var.dimensions):
                stop_with_error('Benchmark variable has an unsupported dimension')

            # input_netcdf selects the first benchmark column. Host buffer
            # normalization and MKS conversions are the NetCDF I/O boundary.
            profile = np.moveaxis(
                np.asarray(var[:]),
                (var.dimensions.index(time_name), var.dimensions.index(z_name)),
                (0, 1),
            )
            if (
                not any(name in var.dimensions for name in ('col', 'column'))
                and (
                    np.prod(profile.shape[2:]) != 1
                    or not any(name in var.dimensions for name in (
                        'X', 'x', 'longitude', 'lon',
                    ))
                    or not any(name in var.dimensions for name in (
                        'Y', 'y', 'latitude', 'lat',
                    ))
                )
            ):
                stop_with_error(
                    'Benchmark file must be a normalized single-column profile dataset'
                )
            profile = jnp.asarray(
                profile.reshape(len(benchmark_time_values_seconds), len(z), -1)[:, :, 0]
            )
            units = var.units.strip()
            if units == 'g/kg':
                profile = profile / g_per_kg
            elif units == 'K/day':
                profile = profile / sec_per_day
            elif units == 'W/m2':
                stop_with_error(
                    'get_netcdf_var: Unable to convert variables of this type to MKS units'
                )

            truth_profiles = []
            for window_idx, (start, end) in enumerate(request.time_window_ranges):
                records = np.flatnonzero(
                    (benchmark_time_values_seconds > start + 1.e-6)
                    & (benchmark_time_values_seconds <= end + 1.e-6)
                )
                if not records.size:
                    stop_with_error(
                        'Requested time subwindow does not align with benchmark time coordinate'
                    )

                if field_idx == 0:
                    print(
                        'benchmark timing:', 'window', window_idx + 1,
                        'records', records[0] + 1, records[-1] + 1,
                        'absolute_window', start, end,
                        'benchmark_window', start, end,
                        'matched_window', benchmark_time_values_seconds[records[0]],
                        benchmark_time_values_seconds[records[-1]],
                    )

                # Accumulate this benchmark time window directly on the CLUBB grid.
                # Adaptation: records and compared levels use fixed array slices.
                interpolated = jnp.where(
                    lower_idx == upper_idx,
                    profile[records[:, None], lower_idx],
                    lin_interpolate_two_points(
                        compared_grid, height_high, z[lower_idx],              # In
                        profile[records[:, None], upper_idx],                  # In
                        profile[records[:, None], lower_idx],                  # In
                    ),
                )
                truth_profiles.append(jnp.sum(interpolated, axis=0) / len(records))
            field.truth_profile = jnp.stack(truth_profiles)

    # Close the shared benchmark file (owned by the context manager).


def calculate_field_loss(
    request, stats_snapshot, window_idx,                      # Intent(in)
):
    """Compute loss and Taylor diagnostics for a completed stats window.

    Inputs: prepared request, completed stats_snapshot and zero-based window_idx.
    Returns: scaled_rmse, correlation, std_ratio, centered_rmse_norm and
    bias_norm, each shaped (requested variable, active batch column).
    Adaptation: divide JAX accumulation banks by their sample counts to obtain
    the source's averaged buffer; columns and levels are vectorized.
    """
    field_metrics = []

    # Score each prepared field independently.
    for field in request.fields:
        benchmark_profile = field.truth_profile[window_idx]
        norm = jnp.max(benchmark_profile) - jnp.min(benchmark_profile)
        norm = jnp.where(is_finite_core_value(norm) & (norm > 0.0), norm, 1.0)

        # Sum the profile mismatch for each parameter column.
        # Collapse the requested height window for this one column.
        bank, slot = field.stats_var_id
        model_profile = stats_snapshot.buffers[bank][
            slot, :, field.k_min:field.k_max + 1
        ]
        counts = stats_snapshot.nsamples[bank][
            slot, :, field.k_min:field.k_max + 1
        ]
        model_profile = jnp.where(
            counts > 0, model_profile / jnp.maximum(counts, 1), 0.0
        )
        scaled_rmse = jnp.sum(((model_profile - benchmark_profile) / norm)**2, axis=-1)

        correlation, std_ratio, centered_rmse_norm, bias_norm = calculate_taylor_metrics(
            model_profile, benchmark_profile,                                # In
        )
        metrics = scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm
        valid = (
            jnp.all(is_finite_core_value(model_profile), axis=-1)
            & jnp.all(is_finite_core_value(benchmark_profile))
        )
        valid &= jnp.all(
            jnp.stack([is_finite_core_value(values) for values in metrics]), axis=0
        )
        field_metrics.append(tuple(
            jnp.where(valid, values, penalty)
            for values, penalty in zip(metrics, set_invalid_field_metric_outputs())
        ))

    return tuple(jnp.stack([values[i] for values in field_metrics]) for i in range(5))
