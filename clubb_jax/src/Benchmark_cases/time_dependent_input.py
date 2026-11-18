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

from functools import partial

import jax
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.constants_clubb import grav, pascal_per_mb, sec_per_hr
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor, zlinterp
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
    read_one_dim_file,
    read_two_dim_file,
    read_x_profile,
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

# Dense forcing values follow the source's (column, level, variable, time) axes.
t_dependent_forcing_data = None
t_dependent_forcing_data_f_grid = None
time_forcing_given = None
forcing_names = ()

l_t_dependent = False
l_input_xpwp_sfc = False
l_ignore_forcings = False
l_sfc_already_initialized = False


# ================================================================================================
def initialize_t_dependent_input(
    ngrdcol,
    iunit,
    forcing_file, sfc_file, iop_file,
    nz,
    grid,
    p_in_Pa,
    grid_adapt_in_time_method,
):
    """This subroutine reads in time dependent information about a
    case that is stored inside the module.

    Fortran declarations (initialize_t_dependent_input in src/Benchmark_cases/time_dependent_input.F90):

        ! Input Variable(s)
        integer, intent(in) :: ngrdcol ! Number of model columns [-]
        integer, intent(in) :: iunit ! File I/O

        character(len=*), intent(in) :: &
          forcing_file, & ! Time-dependent forcing profiles
          sfc_file, & ! Time-dependent surface data
          iop_file ! Optional IOP NetCDF dataset

        integer, intent(in) :: nz ! Size of the model grid

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: grid ! Model grid

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: p_in_Pa ! Pressure[Pa]

        logical, intent(in) :: &
          l_add_dycore_grid, &     ! Flag to set remapping from dycore on or off
          l_readiopdata            ! Flag whether or not we are using a netcdf file for input

        real(kind=core_rknd), dimension(:,:), intent(in) :: rho

        integer, intent(in) :: &
          grid_adapt_in_time_method   ! flag that says whether the grid should be adapted over
    """
    global l_sfc_already_initialized

    # ----------------- Begin Code --------------------
    if not l_ignore_forcings:
        initialize_t_dependent_forcings(
            ngrdcol,
            iunit,
            forcing_file,
            nz,
            grid,
            p_in_Pa,
            grid_adapt_in_time_method,
        )

    if not l_sfc_already_initialized:
        initialize_t_dependent_sfc(iunit, sfc_file)


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
    """This subroutine reads in a file that details time dependent
      input values that vary in one dimension.

    Fortran declarations (initialize_t_dependent_sfc in src/Benchmark_cases/time_dependent_input.F90):

        ! Input Variable(s)
        integer, intent(in) :: iunit ! File I/O unit

        character(len=*), intent(in) :: input_file ! Path to surface input file
        character(len=*), intent(in) :: iop_file ! Optional IOP NetCDF dataset

        logical, intent(in) :: l_readiopdata
        real(kind=core_rknd), dimension(:,:), intent(in) :: rho
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
    ngrdcol,
    iunit,
    input_file,
    nz,
    grid,
    p_in_Pa,
    grid_adapt_in_time_method,
):
    """This subroutine reads in a file that details time dependent
      input values that vary in two dimensions.

    JAX batches fields and model columns with vmap around the compiled
    read_to_grid routine. The result stays on the device with axes
    (column, level, variable, time). The dycore-only l_add_dycore_grid argument
    is omitted as described in the module support note.

    Fortran declarations (initialize_t_dependent_forcings in src/Benchmark_cases/time_dependent_input.F90):

        ! Input Variable(s)
        integer, intent(in) :: ngrdcol ! Number of model columns [-]
        integer, intent(in) :: iunit ! File I/O

        character(len=*), intent(in) :: input_file ! Path to the input file

        integer, intent(in) :: nz  ! Size of Model Grid [-]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: grid ! Altitudes of Grid [m]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: p_in_Pa ! Pressure [Pa]

        logical, intent(in) :: &
          l_add_dycore_grid     ! Flag to set remapping from dycore on or off

        integer, intent(in) :: &
          grid_adapt_in_time_method   ! flag that says whether the grid should be adapted over
    """
    global t_dependent_forcing_data, t_dependent_forcing_data_f_grid
    global time_forcing_given, forcing_names

    if t_dependent_forcing_data_f_grid is None:
        forcing_input, forcing_times = read_two_dim_file(
            iunit, nforcings, input_file,
        )
        # Fill source-file gaps before packing the temporary records into arrays.
        fill_blanks_two_dim_vars(nforcings, forcing_times, forcing_input)
        forcing_names = tuple(field.name for field in forcing_input)
        time_forcing_given = forcing_times.values
        t_dependent_forcing_data_f_grid = np.stack(
            [field.values for field in forcing_input], axis=1,
        )

    n_f_grid_z = t_dependent_forcing_data_f_grid.shape[0]
    n_f_grid_t = time_forcing_given.size
    if forcing_names[0] == z_name:
        model_grid = jnp.asarray(grid)
    elif forcing_names[0] == pressure_name:
        model_grid = -jnp.asarray(p_in_Pa)
        t_dependent_forcing_data_f_grid[:, 0, :] = -np.abs(t_dependent_forcing_data_f_grid[:, 0, :])
    else:
        raise ValueError("Incompatible grid type in first element of t_dependent_forcings")

    source_indices = np.array([forcing_names.index(name) for name in forcing_names[1:]])
    source_values = jnp.asarray(t_dependent_forcing_data_f_grid[:, source_indices, :])
    source_grid = jnp.asarray(t_dependent_forcing_data_f_grid[:, 0, :])

    # Batch independent fields and columns around the compiled interpolation.
    # No column result is converted back to NumPy during table construction.
    interpolate_columns = jax.vmap(read_to_grid, in_axes=(None, None, None, 0, None, None))
    interpolate_variables = jax.vmap(
        interpolate_columns, in_axes=(None, None, None, None, None, 1),
    )
    profiles = interpolate_variables(
        n_f_grid_z, n_f_grid_t, nz, model_grid, source_grid, source_values,
    )
    coordinates = jnp.broadcast_to(model_grid[:, :, None, None], (ngrdcol, nz, 1, n_f_grid_t))
    t_dependent_forcing_data = jnp.concatenate(
        (coordinates, jnp.transpose(profiles, (1, 2, 0, 3))), axis=2,
    )
    time_forcing_given = jnp.asarray(time_forcing_given)

    if grid_adapt_in_time_method == 0:
        t_dependent_forcing_data_f_grid = None


# ================================================================================================
def finalize_t_dependent_forcings():
    """Clear memory initialized in initialize_t_dependent_forcings.
    This should be called at the end of the model.
    """
    global t_dependent_forcing_data, t_dependent_forcing_data_f_grid, time_forcing_given

    # ----------------- Begin Code --------------------
    t_dependent_forcing_data_f_grid = None
    t_dependent_forcing_data = None
    time_forcing_given = None


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
@partial(jax.jit, static_argnums=(0, 1, 2))
def read_to_grid(dim_size, other_dim_size, nz, grid, grid_src, var_src):
    """Interpolate a source field onto one column's model grid at every input time.
    Independent time samples form the batch passed to vertical interpolation.

    Returns:
        var

    Fortran declarations (read_to_grid in src/Benchmark_cases/time_dependent_input.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: dim_size ! Source level count [-]
        integer, intent(in) :: other_dim_size ! Time sample count [-]
        integer, intent(in) :: nz ! Model level count [-]
        real( kind = core_rknd ), intent(in) :: grid(nz) ! Model heights [m] or negative pressures [Pa]
        real( kind = core_rknd ), intent(in) :: &
          grid_src(dim_size,other_dim_size), & ! Source heights [m] or negative pressures [Pa]
          var_src(dim_size,other_dim_size) ! Source field [units vary]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: var(nz,other_dim_size) ! Interpolated field [units vary]
    """
    var = zlinterp(
        other_dim_size, nz, dim_size,
        jnp.broadcast_to(grid, (other_dim_size, nz)), grid_src.T, var_src.T,
    ).T
    return var


# ================================================================================================
def apply_time_dependent_forcings_from_array(
    ngrdcol,
    nzm,
    nzt,
    sclr_dim,
    edsclr_dim,
    sclr_idx,
    l_do_expldiff_rtm_thlm,
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
    """This subroutine is a helper subroutine, that takes the forcings
                   in a 3D array (one set of forcing profiles per column) and
                   initializes the forcing variables.

    JAX uses masked array updates for the Fortran per-column validity
    checks. OpenACC data regions and host/device updates are unnecessary.

    Returns:
        (thlm_f, rtm_f, um_ref, vm_ref, um_f, vm_f, wm_zt, wm_zm, ug, vg,
          sclrm_forcing, edsclrm_forcing)

    Fortran declarations (apply_time_dependent_forcings_from_array in src/Benchmark_cases/time_dependent_input.F90):

        !--------------------- Input Variables ---------------------
        integer, intent(in) :: &
          ngrdcol, &
          nzm, &
          nzt, &
          sclr_dim, &
          edsclr_dim

        type (sclr_idx_type), intent(in) :: &
          sclr_idx

        logical, intent(in) :: &
          l_do_expldiff_rtm_thlm ! thlm/rtm eddy scalars are copies diffused in
                                 ! advance_windm_edsclrm

        type (grid), intent(in) :: &
          gr

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          exner,   & ! Exner Function                             [-]
          rho,     & ! Air Density                                [kg/m^3]
          rtm        ! Total Water Mixing Ratio                   [kg/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,nforcings), intent(in) :: &
          forcings_array ! Forcing profiles (column, level, variable) [units vary]

        !--------------------- Output Variables ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
          thlm_f, & ! Potential Temperature forcing     [K/s]
          rtm_f,  & ! Total Water Mixing Ration forcing [kg/kg/s]
          um_ref, & ! um reference                      [m/s]
          vm_ref, & ! vm reference                      [m/s]
          um_f,   & ! um tendency                       [m/s/s]
          vm_f,   & ! vm tendency                       [m/s/s]
          wm_zt,  & ! subsidence on zt grid             [m/s]
          ug,     & ! u geostrophic wind                [m/s]
          vg        ! v geostrophic wind                [m/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
          wm_zm     ! subsidense on zm grid             [m/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt, sclr_dim), intent(inout) :: &
          sclrm_forcing ! Scalar forcing [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzt, edsclr_dim), intent(inout) :: &
          edsclrm_forcing ! Edscalar forcing [-]
    """

    # --------------------- Begin Code ---------------------
    # With l_do_expldiff_rtm_thlm, advance_windm_edsclrm puts copies of thlm and
    # rtm in the iiedsclr_thl and iiedsclr_rt eddy scalars. thlm and rtm are
    # already forced with thlm_f and rtm_f in advance_xm_wpxp, so the copies
    # must not be forced again.
    l_expldiff_copies = l_do_expldiff_rtm_thlm and edsclr_dim > 1

    # This fixed-size Python loop dispatches heterogeneous forcing metadata
    # (temperature, moisture, winds, and omega have different updates). It is
    # resolved while JAX traces the ten-entry forcing schema; it is not a loop
    # over columns, vertical levels, or another runtime-sized device dimension.
    for n in range(1, nforcings):
        name = forcing_names[n]
        temp_array = forcings_array[:, :, n]
        eps = jnp.finfo(temp_array.dtype).eps

        # Check to see if temp_array is an actual profile or a dummy profile
        # If it is a dummy profile we dont want it to apply itself as it may
        # overwrite legitimate information from another source.
        l_dummy = jnp.any(
            jnp.abs(temp_array - (-999.9))
            < jnp.abs(temp_array + (-999.9)) * 0.5 * eps,
            axis=1, keepdims=True,
        )

        if name in (temperature_f_name, theta_f_name, thetal_f_name):
            if name == temperature_f_name:
                candidate = temp_array / exner
            else:
                candidate = temp_array
            thlm_f = jnp.where(l_dummy, thlm_f, candidate)

            if sclr_idx.iisclr_thl > 0:
                old = sclrm_forcing[:, :, sclr_idx.iisclr_thl - 1]
                sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_thl - 1].set(
                    jnp.where(l_dummy, old, thlm_f)
                )
            if sclr_idx.iiedsclr_thl > 0:
                old = edsclrm_forcing[:, :, sclr_idx.iiedsclr_thl - 1]
                edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_thl - 1].set(
                    jnp.where(l_dummy, old, 0.0 if l_expldiff_copies else thlm_f)
                )

        elif name in (rt_f_name, sp_humidity_f_name):
            if name == sp_humidity_f_name:
                candidate = temp_array * (1.0 + rtm) ** 2
            else:
                candidate = temp_array
            rtm_f = jnp.where(l_dummy, rtm_f, candidate)

            if sclr_idx.iisclr_rt > 0:
                old = sclrm_forcing[:, :, sclr_idx.iisclr_rt - 1]
                sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_rt - 1].set(
                    jnp.where(l_dummy, old, rtm_f)
                )
            if sclr_idx.iiedsclr_rt > 0:
                old = edsclrm_forcing[:, :, sclr_idx.iiedsclr_rt - 1]
                edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_rt - 1].set(
                    jnp.where(l_dummy, old, 0.0 if l_expldiff_copies else rtm_f)
                )

        elif name == um_ref_name:
            candidate = temp_array
            um_ref = jnp.where(l_dummy, um_ref, candidate)
        elif name == vm_ref_name:
            candidate = temp_array
            vm_ref = jnp.where(l_dummy, vm_ref, candidate)
        elif name == um_f_name:
            candidate = temp_array
            um_f = jnp.where(l_dummy, um_f, candidate)
        elif name == vm_f_name:
            candidate = temp_array
            vm_f = jnp.where(l_dummy, vm_f, candidate)
        elif name in (wm_name, omega_name, omega_mb_hr_name):
            if name == wm_name:
                candidate = temp_array
            else:
                if name == omega_mb_hr_name:
                    temp_array = temp_array * pascal_per_mb / sec_per_hr
                candidate = -temp_array / (grav * rho)
            wm_zt = jnp.where(l_dummy, wm_zt, candidate)
            wm_zm_candidate = zt2zm(nzm, nzt, ngrdcol, gr, wm_zt)
            wm_zm = jnp.where(l_dummy, wm_zm, wm_zm_candidate)
        elif name == ug_name:
            candidate = temp_array
            ug = jnp.where(l_dummy, ug, candidate)
        elif name == vg_name:
            candidate = temp_array
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
    l_do_expldiff_rtm_thlm,
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

    Fortran declarations (apply_time_dependent_forcings in src/Benchmark_cases/time_dependent_input.F90):

        !--------------------- Input Variables ---------------------
        integer, intent(in) :: &
          ngrdcol, &
          nzm, &
          nzt, &
          sclr_dim, &
          edsclr_dim

        type (sclr_idx_type), intent(in) :: &
          sclr_idx

        logical, intent(in) :: &
          l_do_expldiff_rtm_thlm ! thlm/rtm eddy scalars are copies diffused in
                                 ! advance_windm_edsclrm

        type (grid), intent(in) :: &
          gr

        real(kind=time_precision), intent(in) :: &
          time ! Model Time [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          exner,   & ! Exner Function                             [-]
          rho,     & ! Air Density                                [kg/m^3]
          rtm        ! Total Water Mixing Ratio                   [kg/kg]

        !--------------------- Output Variables ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
          thlm_f, & ! Potential Temperature forcing     [K/s]
          rtm_f,  & ! Total Water Mixing Ration forcing [kg/kg/s]
          um_ref, & ! um reference                      [m/s]
          vm_ref, & ! vm reference                      [m/s]
          um_f,   & ! um tendency                       [m/s/s]
          vm_f,   & ! vm tendency                       [m/s/s]
          wm_zt,  & ! subsidence on zt grid             [m/s]
          ug,     & ! u geostrophic wind                [m/s]
          vg        ! v geostrophic wind                [m/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
          wm_zm     ! subsidense on zm grid             [m/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt, sclr_dim), intent(inout) :: &
          sclrm_forcing ! Scalar forcing [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzt, edsclr_dim), intent(inout) :: &
          edsclrm_forcing ! Edscalar forcing [-]
    """

    # --------------------- Begin Code ---------------------
    time_frac = -1.0
    before_time, after_time, time_frac = time_select(
        time, time_forcing_given.size, time_forcing_given,
    )
    interpolated = linear_interp_factor(
        time_frac,
        t_dependent_forcing_data[:, :, 1:, after_time],
        t_dependent_forcing_data[:, :, 1:, before_time],
    )
    # The source leaves the coordinate field unused during forcing application.
    forcings_array = jnp.concatenate(
        (jnp.zeros((ngrdcol, nzt, 1), dtype=thlm_f.dtype), interpolated), axis=2,
    )

    return apply_time_dependent_forcings_from_array(
        ngrdcol,
        nzm,
        nzt,
        sclr_dim,
        edsclr_dim,
        sclr_idx,
        l_do_expldiff_rtm_thlm,
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
    """This subroutine determines which indexes of the given
    time_array should be used when interpolating a value
    at the specified time and the location of time between
    these indexes.

    Returns:
        (before_time, after_time, time_frac)

    Fortran declarations (time_select in src/Benchmark_cases/time_dependent_input.F90):

        ! Input Variable(s)

        integer, intent(in) :: nvar                     ! Number of array elements [-]

        real(kind=time_precision), intent(in) :: time   ! Target time              [s]

        real( kind = core_rknd ), dimension(nvar), intent(in) :: time_array ! Array of times [s]

        ! Output Variable(s)

        integer, intent(out) :: &
          after_time, &  ! Index of a time later than the target time [-]
          before_time      ! Index of time before the target time       [-]

        real( kind = core_rknd ), intent(out) :: &
          time_frac      ! The fraction representing the point where time
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
