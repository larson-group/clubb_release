"""Manage time-dependent case input on the CLUBB physics grid.

TODO(port-scope): The Fortran module can ingest forcing on a host model's
dycore grid and remap it during the run. The standalone JAX driver does not own
that grid or its remapping state, so this port intentionally omits the dycore
routine and its arguments instead of exposing a permanently unsupported API.
Restore that source path only with an owning driver and end-to-end tests.

References:
    None.
"""

from __future__ import annotations

from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.constants_clubb import grav, pascal_per_mb, sec_per_hr
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor, zlinterp_fnc
from clubb_jax.src.Input_fields.input_names import (
    CO2_umol_name,
    T_sfc_name,
    latent_ht_name,
    omega_mb_hr_name,
    omega_name,
    pressure_name,
    rt_f_name,
    rt_name,
    sens_ht_name,
    sp_humidity_f_name,
    temperature_f_name,
    thetal_name,
    thetal_f_name,
    theta_f_name,
    time_name,
    ug_name,
    um_f_name,
    um_ref_name,
    upwp_sfc_name,
    vg_name,
    vm_f_name,
    vm_ref_name,
    vpwp_sfc_name,
    wm_name,
    wpqtp_sfc_name,
    wpthlp_sfc_name,
    z_name,
)
from clubb_jax.src.Input_fields.input_reader import (
    count_columns,
    deallocate_one_dim_vars,
    fill_blanks_one_dim_vars,
    fill_blanks_two_dim_vars,
    get_target_index,
    one_dim_read_var,
    read_one_dim_file,
    read_two_dim_file,
    read_x_profile,
    read_x_table,
    two_dim_read_var,
)


nforcings = 10

# Module variables used to describe the surface over time.
time_sfc_given = None
latent_ht_given = None
sens_ht_given = None
thlm_sfc_given = None
rtm_sfc_given = None
CO2_sfc_given = None
upwp_sfc_given = None
vpwp_sfc_given = None
T_sfc_given = None
wpthlp_sfc_given = None
wpqtp_sfc_given = None

t_dependent_forcing_data: list[two_dim_read_var] = []
t_dependent_forcing_data_f_grid: list[two_dim_read_var] = []
dimension_var: one_dim_read_var | None = None

l_t_dependent = False
l_input_xpwp_sfc = False
l_ignore_forcings = False
l_sfc_already_initialized = False

# Resolve the source's ../input/case_setups/ relative path independently of the
# process working directory at this host file-I/O boundary.
input_path = Path(__file__).parents[3] / "input" / "case_setups"
forcings_path = "_forcings.in"
sfc_path = "_sfc.in"


# ================================================================================================
def initialize_t_dependent_input(
    iunit,
    runtype,
    nz,
    grid,
    p_in_Pa,
    grid_adapt_in_time_method,
):
    """This subroutine reads in time dependent information about a
    case that is stored inside the module.

    References:
        None.
    """
    global l_sfc_already_initialized

    # ----------------- Begin Code --------------------
    if not l_ignore_forcings:
        input_file = input_path / f"{runtype}{forcings_path}"
        initialize_t_dependent_forcings(
            iunit,
            input_file,
            nz,
            grid,
            p_in_Pa,
            grid_adapt_in_time_method,
        )

    if not l_sfc_already_initialized:
        input_file = input_path / f"{runtype}{sfc_path}"
        initialize_t_dependent_sfc(iunit, input_file)


# ================================================================================================
def finalize_t_dependent_input():
    """This subroutine frees memory stored after initilizing the
    time dependent data of this module.

    References:
        None.
    """
    # ----------------- Begin Code --------------------
    if not l_ignore_forcings:
        finalize_t_dependent_forcings()
    finalize_t_dependent_sfc()


# ================================================================================================
def initialize_t_dependent_sfc(iunit, input_file):
    """This subroutine reads in a file that details time dependent input values
    that vary in one dimension.
    """
    global time_sfc_given, latent_ht_given, sens_ht_given
    global thlm_sfc_given, rtm_sfc_given, CO2_sfc_given
    global upwp_sfc_given, vpwp_sfc_given, T_sfc_given
    global wpthlp_sfc_given, wpqtp_sfc_given, l_sfc_already_initialized

    # ----------------- Begin Code --------------------
    n_sfc_forcings = count_columns(iunit, input_file)

    # Read the surface.in file and store the necessary input information in retVars.
    retVars = read_one_dim_file(iunit, n_sfc_forcings, input_file)

    # Fill blank values stored as -999.9 using linear interpolation.
    fill_blanks_one_dim_vars(n_sfc_forcings, retVars)

    # dim_size is the number of values input for a particular variable.
    dim_size = retVars[0].values.size

    # Store the data read from the file in each [variable]_sfc_given
    # Target indexes are zero-based, so a found variable has index >= 0 rather
    # than the source's one-based index > 0.
    if get_target_index(n_sfc_forcings, time_name, retVars) >= 0:
        time_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, time_name, retVars, input_file)
        )
    if get_target_index(n_sfc_forcings, latent_ht_name, retVars) >= 0:
        latent_ht_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, latent_ht_name, retVars, input_file)
        )
    if get_target_index(n_sfc_forcings, sens_ht_name, retVars) >= 0:
        sens_ht_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, sens_ht_name, retVars, input_file)
        )
    if get_target_index(n_sfc_forcings, thetal_name, retVars) >= 0:
        thlm_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, thetal_name, retVars, input_file)
        )
    if get_target_index(n_sfc_forcings, rt_name, retVars) >= 0:
        rtm_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, rt_name, retVars, input_file)
        )

    # As of July 2010, this is only in cobra.
    if get_target_index(n_sfc_forcings, CO2_umol_name, retVars) >= 0:
        CO2_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, CO2_umol_name, retVars, input_file)
        )

    # As of July 2010, this is only in gabls3_night
    if get_target_index(n_sfc_forcings, upwp_sfc_name, retVars) >= 0:
        upwp_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, upwp_sfc_name, retVars, input_file)
        )
    if get_target_index(n_sfc_forcings, vpwp_sfc_name, retVars) >= 0:
        vpwp_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, vpwp_sfc_name, retVars, input_file)
        )

    # As of July 2010, this is only in astex_a209
    if get_target_index(n_sfc_forcings, T_sfc_name, retVars) >= 0:
        T_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, T_sfc_name, retVars, input_file)
        )
    if get_target_index(n_sfc_forcings, wpthlp_sfc_name, retVars) >= 0:
        wpthlp_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, wpthlp_sfc_name, retVars, input_file)
        )
    if get_target_index(n_sfc_forcings, wpqtp_sfc_name, retVars) >= 0:
        wpqtp_sfc_given = jnp.asarray(
            read_x_profile(n_sfc_forcings, dim_size, wpqtp_sfc_name, retVars, input_file)
        )

    # Deallocate memory.
    deallocate_one_dim_vars(n_sfc_forcings, retVars)
    l_sfc_already_initialized = True


# ================================================================================================
def initialize_t_dependent_forcings(
    iunit,
    input_file,
    nz,
    grid,
    p_in_Pa,
    grid_adapt_in_time_method,
):
    """This subroutine reads in a file that details time dependent input values
    that vary in two dimensions.
    """
    global t_dependent_forcing_data, t_dependent_forcing_data_f_grid, dimension_var

    # ----------------- Begin Code --------------------
    if not t_dependent_forcing_data_f_grid:
        t_dependent_forcing_data_f_grid, dimension_var = read_two_dim_file(
            iunit, nforcings, input_file
        )

    n_f_grid_z = t_dependent_forcing_data_f_grid[0].values.shape[0]
    n_f_grid_t = dimension_var.values.size

    # Fill in blanks with linear interpolation. Whole profiles of -999.9 will
    # remain that way, thus marking them blank.
    fill_blanks_two_dim_vars(
        nforcings,
        dimension_var,
        t_dependent_forcing_data_f_grid,
    )

    if t_dependent_forcing_data_f_grid[0].name == z_name:
        model_grid = np.asarray(grid)
    elif t_dependent_forcing_data_f_grid[0].name == pressure_name:
        model_grid = -np.asarray(p_in_Pa)
        t_dependent_forcing_data_f_grid[0].values *= -1.0
    else:
        raise ValueError("Incompatible grid type in first element of t_dependent_forcings")

    t_dependent_forcing_data = [
        two_dim_read_var(
            name=t_dependent_forcing_data_f_grid[0].name,
            dim1_name=time_name,
            dim2_name=t_dependent_forcing_data_f_grid[0].name,
            values=jnp.broadcast_to(jnp.asarray(model_grid)[:, None], (nz, n_f_grid_t)),
        )
    ]

    # Host-only setup: forcing names are file metadata, so this loop builds the
    # heterogeneous named records once before any compiled timestep runs.
    for i in range(1, nforcings):
        target_name = t_dependent_forcing_data_f_grid[i].name
        values = read_to_grid(
            nforcings,
            n_f_grid_z,
            n_f_grid_t,
            nz,
            model_grid,
            t_dependent_forcing_data_f_grid,
            target_name,
        )
        t_dependent_forcing_data.append(
            two_dim_read_var(
                name=target_name,
                dim1_name=time_name,
                dim2_name=t_dependent_forcing_data_f_grid[0].name,
                values=jnp.asarray(values),
            )
        )

    dimension_var.values = jnp.asarray(dimension_var.values)

    # Keep the source-grid data only when grid adaptation may require the input
    # to be interpolated again during the run.
    if grid_adapt_in_time_method == 0:
        t_dependent_forcing_data_f_grid = []


# ================================================================================================
def finalize_t_dependent_forcings():
    """Clear memory initialized in initialize_t_dependent_forcings.
    This should be called at the end of the model.
    """
    global t_dependent_forcing_data, t_dependent_forcing_data_f_grid, dimension_var

    # ----------------- Begin Code --------------------
    t_dependent_forcing_data_f_grid = []
    t_dependent_forcing_data = []
    dimension_var = None


# ================================================================================================
def finalize_t_dependent_sfc():
    """Clear memory initialized in initialize_t_dependent_surface.
    This should be called at the end of the model.
    """
    global time_sfc_given, latent_ht_given, sens_ht_given
    global thlm_sfc_given, rtm_sfc_given, CO2_sfc_given
    global upwp_sfc_given, vpwp_sfc_given, T_sfc_given
    global wpthlp_sfc_given, wpqtp_sfc_given, l_sfc_already_initialized

    # ----------------- Begin Code --------------------
    time_sfc_given = None
    latent_ht_given = None
    sens_ht_given = None
    thlm_sfc_given = None
    rtm_sfc_given = None
    CO2_sfc_given = None
    upwp_sfc_given = None
    vpwp_sfc_given = None
    T_sfc_given = None
    wpthlp_sfc_given = None
    wpqtp_sfc_given = None
    # The Python process may initialize more than one SCM run; unlike a
    # standalone Fortran process, it must reset this module lifecycle flag.
    l_sfc_already_initialized = False


# ================================================================================================
def read_to_grid(
    ntwo_dim_vars,
    dim_size,
    other_dim_size,
    nz,
    grid,
    two_dim_vars,
    target_name,
):
    """This is a helper function for doing the translation from the forcing
    grid to the model grid.
    """
    # ----------------- Begin Code --------------------
    temp_var = read_x_table(
        ntwo_dim_vars,
        dim_size,
        other_dim_size,
        target_name,
        two_dim_vars,
    )
    # Each time slice is interpolated independently onto the same model grid.
    # Vectorizing here also keeps the host-side input path consistent with the
    # column-vectorized device routines.
    interpolated = np.asarray(
        jax.vmap(
            lambda source_grid, source_values: zlinterp_fnc(
                grid, source_grid, source_values
            ),
            in_axes=(1, 1),
            out_axes=1,
        )(jnp.asarray(two_dim_vars[0].values), jnp.asarray(temp_var))
    )
    return interpolated.reshape(nz, other_dim_size)


# ================================================================================================
def apply_time_dependent_forcings_from_array(
    ngrdcol,
    nzm,
    nzt,
    sclr_dim,
    edsclr_dim,
    sclr_idx,
    gr,
    rtm,
    rho,
    exner,
    forcings_array,
    thlm_f,
    rtm_f,
    um_ref,
    vm_ref,
    um_f,
    vm_f,
    wm_zt,
    wm_zm,
    ug,
    vg,
    sclrm_forcing,
    edsclrm_forcing,
):
    """This helper subroutine takes the forcings in a 2D matrix and initializes
    the forcing variables.
    """

    # --------------------- Begin Code ---------------------
    # This fixed-size Python loop dispatches heterogeneous forcing metadata
    # (temperature, moisture, winds, and omega have different updates). It is
    # resolved while JAX traces the ten-entry forcing schema; it is not a loop
    # over columns, vertical levels, or another runtime-sized device dimension.
    for n in range(1, nforcings):
        name = t_dependent_forcing_data[n].name
        temp_array = forcings_array[n, :]
        eps = jnp.finfo(temp_array.dtype).eps

        # Check to see if temp_array is an actual profile or a dummy profile
        # If it is a dummy profile we dont want it to apply itself as it may
        # overwrite legitimate information from another source.
        l_dummy = jnp.any(
            jnp.abs(temp_array - (-999.9))
            < jnp.abs(temp_array + (-999.9)) * 0.5 * eps
        )

        if name in (temperature_f_name, theta_f_name, thetal_f_name):
            if name == temperature_f_name:
                candidate = temp_array[jnp.newaxis, :] / exner
            else:
                candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            thlm_f = jnp.where(l_dummy, thlm_f, candidate)

            if sclr_idx.iisclr_thl > 0:
                old = sclrm_forcing[:, :, sclr_idx.iisclr_thl - 1]
                sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_thl - 1].set(
                    jnp.where(l_dummy, old, thlm_f)
                )
            if sclr_idx.iiedsclr_thl > 0:
                old = edsclrm_forcing[:, :, sclr_idx.iiedsclr_thl - 1]
                edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_thl - 1].set(
                    jnp.where(l_dummy, old, thlm_f)
                )

        elif name in (rt_f_name, sp_humidity_f_name):
            if name == sp_humidity_f_name:
                candidate = temp_array[jnp.newaxis, :] * (1.0 + rtm) ** 2
            else:
                candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            rtm_f = jnp.where(l_dummy, rtm_f, candidate)

            if sclr_idx.iisclr_rt > 0:
                old = sclrm_forcing[:, :, sclr_idx.iisclr_rt - 1]
                sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_rt - 1].set(
                    jnp.where(l_dummy, old, rtm_f)
                )
            if sclr_idx.iiedsclr_rt > 0:
                old = edsclrm_forcing[:, :, sclr_idx.iiedsclr_rt - 1]
                edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_rt - 1].set(
                    jnp.where(l_dummy, old, rtm_f)
                )

        elif name == um_ref_name:
            candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            um_ref = jnp.where(l_dummy, um_ref, candidate)
        elif name == vm_ref_name:
            candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            vm_ref = jnp.where(l_dummy, vm_ref, candidate)
        elif name == um_f_name:
            candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            um_f = jnp.where(l_dummy, um_f, candidate)
        elif name == vm_f_name:
            candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            vm_f = jnp.where(l_dummy, vm_f, candidate)
        elif name in (wm_name, omega_name, omega_mb_hr_name):
            if name == wm_name:
                candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            else:
                if name == omega_mb_hr_name:
                    temp_array = temp_array * pascal_per_mb / sec_per_hr
                candidate = -temp_array[jnp.newaxis, :] / (grav * rho)
            wm_zt = jnp.where(l_dummy, wm_zt, candidate)
            wm_zm_candidate = zt2zm(nzm, nzt, ngrdcol, gr, wm_zt)
            wm_zm = jnp.where(l_dummy, wm_zm, wm_zm_candidate)
        elif name == ug_name:
            candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            ug = jnp.where(l_dummy, ug, candidate)
        elif name == vg_name:
            candidate = jnp.broadcast_to(temp_array, (ngrdcol, temp_array.size))
            vg = jnp.where(l_dummy, vg, candidate)
        else:
            raise ValueError(f"Incompatible forcing type: {name}")

    return (
        thlm_f,
        rtm_f,
        um_ref,
        vm_ref,
        um_f,
        vm_f,
        wm_zt,
        wm_zm,
        ug,
        vg,
        sclrm_forcing,
        edsclrm_forcing,
    )


# ================================================================================================
def apply_time_dependent_forcings(
    ngrdcol,
    nzm,
    nzt,
    sclr_dim,
    edsclr_dim,
    sclr_idx,
    gr,
    time,
    rtm,
    rho,
    exner,
    thlm_f,
    rtm_f,
    um_ref,
    vm_ref,
    um_f,
    vm_f,
    wm_zt,
    wm_zm,
    ug,
    vg,
    sclrm_forcing,
    edsclrm_forcing,
):
    """This subroutine converts the time dependent information stored in
    memory (time_dependent_forcing_data) into the format used by CLUBB.
    """

    # --------------------- Begin Code ---------------------
    time_frac = -1.0
    before_time, after_time, time_frac = time_select(
        time,
        dimension_var.values.size,
        dimension_var.values,
    )

    forcing_values = jnp.stack(
        [forcing.values for forcing in t_dependent_forcing_data[1:nforcings]]
    )
    interpolated = linear_interp_factor(
        time_frac,
        forcing_values[:, :, after_time],
        forcing_values[:, :, before_time],
    )
    # The source leaves row 1 unused; JAX arrays initialize it explicitly.
    forcings_array = jnp.concatenate(
        (jnp.zeros((1, nzt), dtype=thlm_f.dtype), interpolated), axis=0
    )

    return apply_time_dependent_forcings_from_array(
        ngrdcol,
        nzm,
        nzt,
        sclr_dim,
        edsclr_dim,
        sclr_idx,
        gr,
        rtm,
        rho,
        exner,
        forcings_array,
        thlm_f,
        rtm_f,
        um_ref,
        vm_ref,
        um_f,
        vm_f,
        wm_zt,
        wm_zm,
        ug,
        vg,
        sclrm_forcing,
        edsclrm_forcing,
    )


# ================================================================================================
def time_select(time, nvar, time_array):
    """Determine which indexes of time_array should be used when interpolating
    a value at the specified time and the location of time between them.

    References:
        None.
    """
    times = jnp.asarray(time_array)

    # Default initialization
    before_time = jnp.asarray(-999, dtype=jnp.int32)
    after_time = jnp.asarray(-999, dtype=jnp.int32)

    # convert time to a real so it has the same precision as the values
    # in time_array
    # TODO(port-mirror): Compiled JAX calls cannot raise the source ERROR STOP.
    # Eager calls retain it; remove this split when the driver carries a
    # functional diagnostic status out of the forcing kernel.
    if not isinstance(time, jax.core.Tracer):
        times_host = np.asarray(times)
        if time < times_host[0]:
            # If time is less than the lowest value in time_array, an invalid
            # time has been provided. Stop execution.
            raise ValueError("Selected time is before the first time at which data are available")
        if time > times_host[-1]:
            # If time is greater than the highest value in time_array, an invalid
            # time has been provided. Stop execution.
            raise ValueError("Selected time is after the last time at which data are available")
        if np.any(times_host[: nvar - 1] >= times_host[1:nvar]):
            # If times are not increasing then they aren't sorted properly.
            raise ValueError("times are not sorted; check the case surface and forcing files")

    # searchsorted(side="left") selects the source's first inclusive bracket
    # at exact interior time nodes.
    after_time = jnp.searchsorted(times[:nvar], time, side="left")
    after_time = jnp.clip(after_time, 1, nvar - 1)
    before_time = after_time - 1

    # Compute the position of time between before_time and after_time
    # as a fraction.
    time_frac = (
        (time - times[before_time])
        / (times[after_time] - times[before_time])
    )
    return before_time, after_time, time_frac


__all__ = [
    "finalize_t_dependent_input",
    "time_select",
    "apply_time_dependent_forcings",
    "initialize_t_dependent_input",
]
