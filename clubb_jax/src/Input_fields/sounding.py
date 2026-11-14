"""Read a case's sounding files and initialize the model columns from them.

JAX port of sounding.F90. A case describes its initial atmosphere in text files
under input/case_data. This module reads those files and interpolates them onto
the model grid.

Usage:
    read_sounding: read the sounding and any scalar soundings, and return the
        initial model profiles and surface values. Call it once when a case is
        initialized, after setting up the grid.
    read_sounding_file: read the configured main sounding onto its own levels.
    read_sclr_sounding_file, read_edsclr_sounding_file: read the scalar
        soundings, which share the levels of the main sounding.
    convert_snd2extended_atm, load_extended_std_atm: build the extended
        atmosphere that BUGSrad uses above the model top.

File format:
    Each file has a header naming its columns, followed by one row per level.
    Comment lines at the top contain "!", and -999.9 marks a missing value. The
    main sounding has 8 columns, in any order. For example:

        z[m]  thm[K]  'rt[kg/kg]'  'u[m/s]'  'v[m/s]'  'w[m/s]'  'ug[m/s]'  'vg[m/s]'

    Its first column may be altitude (z[m]) or pressure (Press[Pa]), its
    temperature may be thm[K], thlm[K] or T[K], and its vertical velocity may be
    w[m/s] or omega[Pa/s]. input_interpret.py interprets these columns.

Porting deviations:
    * Fortran returns outputs through arguments and stores the extended
      atmosphere in module variables. These routines return them instead.
    * Fortran sizes its sounding arrays with nmaxsnd. Here they have nlevels
      values, except sclr_snd and edsclr_snd, which keep Fortran's zero padding
      for a scalar file with fewer rows than the main sounding.

References:
    McClatchey, et al., (1972) _Environmental Research Papers_, No. 411, p.94
    Wyant, et al. (2007), for the DYCOMS II RF02 profiles
"""

from pathlib import Path

import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.constants_clubb import (
    g_per_kg,
    kappa,
    p0,
    pascal_per_mb,
)
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.interpolation import (
    lin_interpolate_two_points,
    mono_cubic_interp,
)
from clubb_jax.src.Input_fields.input_interpret import (
    read_z_profile,
    read_theta_profile,
    read_subs_profile,
)
from clubb_jax.src.Input_fields.input_names import (
    z_name,
    theta_name,
    thetal_name,
    temperature_name,
    wm_name,
    rt_name,
    um_name,
    vm_name,
    ug_name,
    vg_name,
    CO2_name,
    ozone_name,
    press_mb_name,
    sp_humidity_name,
)
from clubb_jax.src.Input_fields.input_reader import (
    deallocate_one_dim_vars,
    fill_blanks_one_dim_vars,
    read_one_dim_file,
    read_x_profile,
)

n_snd_var = 8  # Constant parameter

sclr_max = 1000  # Maximum number of scalars


def read_sounding(ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
                  gr, iunit, runtype, sounding_file, sclr_sounding_file,
                  edsclr_sounding_file, ozone_sounding_file, iop_file, zm_init,
                  saturation_formula,
                  l_readiopdata, l_modify_ic_with_cubic_int,
                  p_sfc):
    """Read the configured sounding files and interpolate initial fields onto the model grid.

    Returns:
        (thlm, temperature_type, rtm, um, vm, ugm, vgm, altitude_type, press,
          subs_type, wm, rtm_sfc, thlm_sfc, sclrm, edsclrm)

    Fortran declarations (read_sounding in src/Input_fields/sounding.F90):

        !--------------------- Input Variables ---------------------
        integer, intent(in) :: &
          ngrdcol, &
          sclr_dim, &
          edsclr_dim

        type (sclr_idx_type), intent(in) :: &
          sclr_idx

        type(grid), intent(in) :: &
          gr

        integer, intent(in) :: iunit ! File unit to use for namelist

        character(len=*), intent(in) ::  &
          runtype ! Used to determine if this in a DYCOMS II RF02 simulation

        character(len=*), intent(in) :: &
          sounding_file, & ! Main initial profile
          sclr_sounding_file, & ! Higher-order passive scalar profiles
          edsclr_sounding_file, & ! Eddy-diffusivity passive scalar profiles
          ozone_sounding_file, & ! Ozone profile for the extended atmosphere
          iop_file ! IOP NetCDF dataset

        real( kind = core_rknd ), intent(in) :: &
          zm_init ! Height at zm(1)         [m]

        integer, intent(in) :: &
          saturation_formula ! Integer that stores the saturation formula to be used

        ! Flag for interpolating the sounding profile with Steffen's monotone cubic
        ! method to obtain smoother initial condition profile, which is found to be
        ! beneficial to achive a better numerical solution convergence. If this flag
        ! is turned off, the initial conditions will be generated with linear interpolation.
        ! This is done on a case-by-case basis, since using the monotone cubic method
        ! requires a special sounding.in file with many additional sounding levels.
        logical, intent(in) :: &
          l_modify_ic_with_cubic_int, &
          l_readiopdata

        real( kind = core_rknd ), intent(inout) :: &
          p_sfc ! Pressure at the surface [Pa]

        !--------------------- Output Variables ---------------------
        real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) ::  &
          thlm,  & ! Liquid water potential temperature    [K]
          rtm,   & ! Total water mixing ratio              [kg/kg]
          um,    & ! u wind component                      [m/s]
          vm,    & ! v wind component                      [m/s]
          ugm,   & ! u geostrophic wind component          [m/s]
          vgm,   & ! v geostrophic wind component          [m/s]
          press, & ! Pressure                              [Pa]
          wm       ! Subsidence                            [m/s or Pa/s]

        real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  &
          rtm_sfc,  & ! Initial surface rtm                [kg/kg]
          thlm_sfc    ! Initial surface thlm               [K]

        character(len=*), intent(out) :: &
          temperature_type, &     ! Type of temperature sounding
          altitude_type, &       ! Type of independent coordinate
          subs_type         ! Type of subsidence

        real( kind = core_rknd ), intent(out), dimension(ngrdcol, gr%nzt, sclr_dim) ::  &
          sclrm   ! Passive scalar output      [units vary]

        real( kind = core_rknd ), intent(out), dimension(ngrdcol, gr%nzt, edsclr_dim) ::  &
          edsclrm ! Eddy Passive scalar output [units vary]
    """
    # Constant parameter
    nmaxsnd = 10000

    # Sounding profile retVars stay empty for a scalar file that is not read.
    sclr_sounding_retVars = []

    # ---- Begin Code ----

    temperature_type = theta_name  # Default value
    altitude_type = z_name  # Default value
    subs_type = wm_name  # Default value

    # Determine which files exist ahead of time to allow for a graceful exit if
    # one is missing.
    l_sounding_exists = bool(sounding_file) and Path(sounding_file).is_file()

    l_sclr_sounding_exists = bool(sclr_sounding_file) and Path(sclr_sounding_file).is_file()

    l_edsclr_sounding_exists = bool(edsclr_sounding_file) and Path(edsclr_sounding_file).is_file()

    # ---------------------------------------------------------------------------------------------
    # Status Message
    if clubb_at_least_debug_level(1):
        print("Path to sounding: ", sounding_file)
        print("File exists? ", l_sounding_exists)
        print("Path to sclr_sounding: ", sclr_sounding_file)
        print("File exists? ", l_sclr_sounding_exists)
        print("Path to edsclr_sounding: ", edsclr_sounding_file)
        print("File exists? ", l_edsclr_sounding_exists)
    # ---------------------------------------------------------------------------------------------

    if l_readiopdata:
        # TODO(port-mirror): readiopdata_module.F90 is not ported, as listed in
        # test_file_name_parity.py. Port it to support this branch.
        raise ValueError("l_readiopdata = true is not supported.")
    elif l_sounding_exists:
        # Read the configured main sounding.
        (nlevels, z, theta, rt, u, v,
         ug, vg, p_in_Pa, subs,
         temperature_type, altitude_type, subs_type,
         sounding_retVars) = read_sounding_file(iunit, sounding_file, nmaxsnd, p_sfc,
                                                saturation_formula, zm_init)
    else:
        raise ValueError(f"Cannot open {sounding_file}")
        # sounding namelist is no longer used.
        # Joshua Fasching April 2009

    # Read in a passive scalar sounding, if enabled
    if sclr_dim > 0 or edsclr_dim > 0:

        # Allocate large arrays, using only the columns that hold scalars. Initialize to zero
        sclr_snd = np.zeros((nmaxsnd, sclr_dim))
        edsclr_snd = np.zeros((nmaxsnd, edsclr_dim))

        # Read the configured passive and eddy-diffusivity scalar soundings.
        if sclr_dim > 0:
            if l_sclr_sounding_exists:
                sclr_snd, sclr_sounding_retVars = read_sclr_sounding_file(
                    sclr_dim, sclr_idx, iunit, sclr_sounding_file, nmaxsnd,
                    sclr_snd,
                )
            else:
                raise ValueError(f"Cannot open {sclr_sounding_file}")

        if edsclr_dim > 0:
            if l_edsclr_sounding_exists:
                edsclr_snd = read_edsclr_sounding_file(
                    edsclr_dim, sclr_idx, iunit, edsclr_sounding_file, nmaxsnd,
                    edsclr_snd,
                )
            else:
                raise ValueError(f"Cannot open {edsclr_sounding_file}")

    if nlevels > nmaxsnd:
        raise ValueError(
            f"Error in sounding: nlevels > nmaxsnd; nlevels = {nlevels}, nmaxsnd = {nmaxsnd}"
        )

    zt = np.asarray(gr.zt)
    zm = np.asarray(gr.zm)

    # Error check: if lowest above-model-surface themodynamic grid height
    # (gr%zt(1,1)) is lower than the lowest value from the input sounding,
    # then the linear interpolation scheme will fail.

    if zt[0, 0] < z[0]:
        raise ValueError(
            "Lowest level of input sounding, z(1), must be"
            " below the first above-model-surface thermodynamic level, gr%zt(1,1)."
            f" First sounding level z(1) = {z[0]}."
            f" First thermodynamic level gr%zt(1,1) = {zt[0, 0]}."
        )

    # First sounding level should be near ground value

    # dschanen 1 May 2007
    # We have changed this for Nov. 11 and June 25, both of which
    # begin above the ground.
    nzt = zt.shape[1]
    um = np.full((ngrdcol, nzt), u[0])
    vm = np.full((ngrdcol, nzt), v[0])
    ugm = np.full((ngrdcol, nzt), ug[0])
    vgm = np.full((ngrdcol, nzt), vg[0])
    thlm = np.full((ngrdcol, nzt), theta[0])
    rtm = np.full((ngrdcol, nzt), rt[0])
    press = np.full((ngrdcol, nzt), p_in_Pa[0])
    wm = np.full((ngrdcol, nzt), subs[0])

    sclrm = np.zeros((ngrdcol, nzt, sclr_dim))
    edsclrm = np.zeros((ngrdcol, nzt, edsclr_dim))

    if sclr_dim > 0:
        sclrm[:, 0, :] = sclr_snd[0, :]

    if edsclr_dim > 0:
        edsclrm[:, 0, :] = edsclr_snd[0, :]

    if clubb_at_least_debug_level(1):

        print("Reading in sounding information")
        # ------------Printing Model Inputs-------------------------------
        print("z = ", z)
        print("u = ", u)
        print("v = ", v)
        print("ug = ", ug)
        print("vg = ", vg)
        print("theta = ", theta)
        print("rt = ", rt)
        print("p_in_Pa = ", p_in_Pa)
        print("subs = ", subs)

        for sclr in range(sclr_dim):
            print(f"sclr_snd({sclr + 1}) = ", sclr_snd[:nlevels, sclr])

        for edsclr in range(edsclr_dim):
            print(f"edsclr_snd({edsclr + 1}) = ", edsclr_snd[:nlevels, edsclr])

    # ----------------------------------------------------------------------
    # Use linear interpolation from two nearest prescribed grid points
    # (one above and one below) to initialize mean quantities in the model
    # Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )

    # Fortran walks k_above up from the bottom and interpolates at every step, so
    # only the first level with z(k_above) >= gr%zt(i,k) sets the result. Here
    # k_above is that level, found once for all columns and heights.
    l_z_above = z[None, None, :] >= zt[:, :, None]

    if not np.all(np.any(l_z_above, axis=2)):
        raise ValueError(
            "STOP Not enough sounding data to initialize grid:"
            f" highest sounding level {z[nlevels - 1]}"
            f" should be higher than highest thermodynamic point {zt.max()}"
        )

    k_above = np.argmax(l_z_above, axis=2)

    # A level at or below z(1) never enters the Fortran loop, so it keeps the
    # values set above.
    l_interp = k_above > 0
    k_above = np.maximum(k_above, 1)

    # situation w/ cubic int. (achieve better numerical solution convergence)
    if l_modify_ic_with_cubic_int:
        # use Steffen's monotone cubic interpolation method to obtain
        # smoothing initial condition profile for convergence test
        # note: vertical levels in sounding file need to be not too coarse
        #
        # k_above is zero-based here, so 1 is the Fortran k_above == 2 branch.
        km1 = np.where(k_above == 1, 0, k_above - 2)
        kp1 = k_above
        kp2 = np.where(k_above == 1, 2, k_above + 1)
        k00 = np.where(k_above == 1, 0, k_above - 1)

        # if z(k_above) reaches at the top level in sounding profile,
        # then use the nearest levels for interpolation
        l_top = (k_above > 1) & (z[k_above] >= z[nlevels - 1])
        km1 = np.where(l_top, nlevels - 3, km1)
        kp1 = np.where(l_top, nlevels - 1, kp1)
        kp2 = np.where(l_top, nlevels - 1, kp2)
        k00 = np.where(l_top, nlevels - 2, k00)

        # mono_cubic_interp takes one height at a time.
        mono_cubic_interp_vec = jnp.vectorize(mono_cubic_interp)

        um = np.where(l_interp, mono_cubic_interp_vec(zt, km1, k00, kp1, kp2, z[km1], z[k00], z[kp1],
                                                      z[kp2], u[km1], u[k00], u[kp1], u[kp2]), um)
        vm = np.where(l_interp, mono_cubic_interp_vec(zt, km1, k00, kp1, kp2, z[km1], z[k00], z[kp1],
                                                      z[kp2], v[km1], v[k00], v[kp1], v[kp2]), vm)
        ugm = np.where(l_interp, mono_cubic_interp_vec(zt, km1, k00, kp1, kp2, z[km1], z[k00], z[kp1],
                                                       z[kp2], ug[km1], ug[k00], ug[kp1], ug[kp2]), ugm)
        vgm = np.where(l_interp, mono_cubic_interp_vec(zt, km1, k00, kp1, kp2, z[km1], z[k00], z[kp1],
                                                       z[kp2], vg[km1], vg[k00], vg[kp1], vg[kp2]), vgm)
        thlm = np.where(l_interp, mono_cubic_interp_vec(zt, km1, k00, kp1, kp2, z[km1], z[k00], z[kp1],
                                                        z[kp2], theta[km1], theta[k00], theta[kp1], theta[kp2]), thlm)
        rtm = np.where(l_interp, mono_cubic_interp_vec(zt, km1, k00, kp1, kp2, z[km1], z[k00], z[kp1],
                                                       z[kp2], rt[km1], rt[k00], rt[kp1], rt[kp2]), rtm)
        press = np.where(l_interp, mono_cubic_interp_vec(zt, km1, k00, kp1, kp2, z[km1], z[k00], z[kp1],
                                                         z[kp2], p_in_Pa[km1], p_in_Pa[k00], p_in_Pa[kp1],
                                                         p_in_Pa[kp2]), press)
        wm = np.where(l_interp, mono_cubic_interp_vec(zt, km1, k00, kp1, kp2, z[km1], z[k00], z[kp1],
                                                      z[kp2], subs[km1], subs[k00], subs[kp1], subs[kp2]), wm)

        if runtype != "dycoms2_rf02":
            # initial condition for tracers
            if sclr_dim > 0:
                for sclr in range(sclr_dim):
                    sclrm[:, :, sclr] = np.where(
                        l_interp,
                        mono_cubic_interp_vec(zt, km1, k00, kp1, kp2,
                                              z[km1], z[k00], z[kp1], z[kp2],
                                              sclr_snd[km1, sclr], sclr_snd[k00, sclr],
                                              sclr_snd[kp1, sclr], sclr_snd[kp2, sclr]),
                        sclrm[:, :, sclr],
                    )
            if edsclr_dim > 0:
                for edsclr in range(edsclr_dim):
                    edsclrm[:, :, edsclr] = np.where(
                        l_interp,
                        mono_cubic_interp_vec(zt, km1, k00, kp1, kp2,
                                              z[km1], z[k00], z[kp1], z[kp2],
                                              edsclr_snd[km1, edsclr], edsclr_snd[k00, edsclr],
                                              edsclr_snd[kp1, edsclr], edsclr_snd[kp2, edsclr]),
                        edsclrm[:, :, edsclr],
                    )
        else:
            # DYCOMS2_RF02 case (use the same treatment as in regular situation w/ linear int.)
            ugm = np.where(l_interp, um, ugm)
            vgm = np.where(l_interp, vm, vgm)
            # Passive Scalars
            # Change this if they are not equal to theta_l and rt in RF02
            if sclr_idx.iisclr_thl > 0:
                sclrm[:, :, sclr_idx.iisclr_thl - 1] = np.where(l_interp, thlm, sclrm[:, :, sclr_idx.iisclr_thl - 1])
                edsclrm[:, :, sclr_idx.iisclr_thl - 1] = np.where(l_interp, thlm, edsclrm[:, :, sclr_idx.iisclr_thl - 1])
            if sclr_idx.iisclr_rt > 0:
                sclrm[:, :, sclr_idx.iisclr_rt - 1] = np.where(l_interp, rtm, sclrm[:, :, sclr_idx.iisclr_rt - 1])
                edsclrm[:, :, sclr_idx.iisclr_rt - 1] = np.where(l_interp, rtm, edsclrm[:, :, sclr_idx.iisclr_rt - 1])

    else:  # default model setup

        # Regular situation w/ linear int.
        if runtype != "dycoms2_rf02":

            um = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], u[k_above], u[k_above - 1]), um)
            vm = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], v[k_above], v[k_above - 1]), vm)
            ugm = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], ug[k_above], ug[k_above - 1]), ugm)
            vgm = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], vg[k_above], vg[k_above - 1]), vgm)
            thlm = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], theta[k_above], theta[k_above - 1]), thlm)
            rtm = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], rt[k_above], rt[k_above - 1]), rtm)
            press = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], p_in_Pa[k_above],
                                                                  p_in_Pa[k_above - 1]), press)
            wm = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], subs[k_above], subs[k_above - 1]), wm)

            if sclr_dim > 0:
                for sclr in range(sclr_dim):
                    sclrm[:, :, sclr] = np.where(
                        l_interp,
                        lin_interpolate_two_points(zt, z[k_above], z[k_above - 1],
                                                   sclr_snd[k_above, sclr], sclr_snd[k_above - 1, sclr]),
                        sclrm[:, :, sclr],
                    )
            if edsclr_dim > 0:
                for edsclr in range(edsclr_dim):
                    edsclrm[:, :, edsclr] = np.where(
                        l_interp,
                        lin_interpolate_two_points(zt, z[k_above], z[k_above - 1],
                                                   edsclr_snd[k_above, edsclr], edsclr_snd[k_above - 1, edsclr]),
                        edsclrm[:, :, edsclr],
                    )

        else:  # DYCOMS II RF02 case

            l_below_795 = zt < 795.0

            # (Wyant, et al. 2007, eq 1--4)
            um = np.where(l_interp, 3.0 + (4.3 * zt) / 1000.0, um)  # Known magic number
            vm = np.where(l_interp, -9.0 + (5.6 * zt) / 1000.0, vm)  # Known magic number
            ugm = np.where(l_interp, um, ugm)
            vgm = np.where(l_interp, vm, vgm)
            # The clip keeps the cube root defined where l_below_795 selects 288.3.
            thlm = np.where(
                l_interp,
                np.where(l_below_795, 288.3, 295.0 + np.clip(zt - 795.0, 0.0, None)**(1.0 / 3.0)),  # Known magic number
                thlm,
            )
            rtm = np.where(
                l_interp,
                np.where(
                    l_below_795,
                    9.45 / g_per_kg,  # Known magic number
                    (5.0 - 3.0 * (1.0 - np.exp((795.0 - zt) / 500.0))) / g_per_kg,  # Known magic number
                ),
                rtm,
            )
            # Passive Scalars
            # Change this if they are not equal to theta_l and rt in RF02
            if sclr_idx.iisclr_thl > 0:
                sclrm[:, :, sclr_idx.iisclr_thl - 1] = np.where(l_interp, thlm, sclrm[:, :, sclr_idx.iisclr_thl - 1])
                edsclrm[:, :, sclr_idx.iisclr_thl - 1] = np.where(l_interp, thlm, edsclrm[:, :, sclr_idx.iisclr_thl - 1])
            if sclr_idx.iisclr_rt > 0:
                sclrm[:, :, sclr_idx.iisclr_rt - 1] = np.where(l_interp, rtm, sclrm[:, :, sclr_idx.iisclr_rt - 1])
                edsclrm[:, :, sclr_idx.iisclr_rt - 1] = np.where(l_interp, rtm, edsclrm[:, :, sclr_idx.iisclr_rt - 1])
            press = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1],
                                                                  p_in_Pa[k_above], p_in_Pa[k_above - 1]), press)
            wm = np.where(l_interp, lin_interpolate_two_points(zt, z[k_above], z[k_above - 1], subs[k_above], subs[k_above - 1]), wm)

    # The sounding will be initialized to thermodynamic grid levels successfully
    # as long as the thermodynamic level 1 (at altitude gr%zt(1,1)) is at or above
    # the lowest sounding level (z(1)).  However, it is advantageous to know the
    # initial surface values of a few variables, as long as the sounding extends
    # to the surface, which is found at momentum level 1 (at altitude gr%zm(1,1)).
    zm_sfc = zm[:, 0]

    # Find the two sounding levels that the surface (gr%zm(1,1)) is found
    # between.  idx is the index of the closest value greater than or equal to
    # gr%zm(1,1), and at least the second level, as in Fortran's binary_search.
    idx = np.maximum(np.argmax(z[None, :] >= zm_sfc[:, None], axis=1), 1)

    # The surface (or model lower boundary) can be below the lowest sounding
    # level.  Initialize the values of rtm_sfc and thlm_sfc to negative
    # values that will be overwritten later.  Otherwise the surface is found
    # between sounding levels idx and idx-1, so use linear interpolation to find
    # the values of rtm_sfc and thlm_sfc.
    rtm_sfc = np.where(zm_sfc < z[0], -999.0,
                       lin_interpolate_two_points(zm_sfc, z[idx], z[idx - 1],
                                                  rt[idx], rt[idx - 1]))

    thlm_sfc = np.where(zm_sfc < z[0], -999.0,
                        lin_interpolate_two_points(zm_sfc, z[idx], z[idx - 1],
                                                   theta[idx], theta[idx - 1]))

    # TODO(port-mirror): Fortran prepares the extended sounding for radiation here
    # when rad_scheme == "bugsrad", through load_extended_std_atm or
    # convert_snd2extended_atm. clubb_case_initialization.py rejects bugsrad, so
    # this branch cannot run. Add the calls when JAX supports bugsrad.

    # Deallocate sounding and scalar sounding profiles.
    deallocate_one_dim_vars(n_snd_var, sounding_retVars)
    deallocate_one_dim_vars(len(sclr_sounding_retVars), sclr_sounding_retVars)

    return (thlm, temperature_type, rtm, um, vm, ugm, vgm,
            altitude_type, press, subs_type, wm,
            rtm_sfc, thlm_sfc, sclrm, edsclrm)


# -----------------------------------------------------------------------------
def read_sounding_file(iunit, sounding_file, nmaxsnd, p_sfc,
                       saturation_formula, zm_init):
    """This subroutine reads in a configured sounding.in file and
      returns the values contained in that file.

    Returns:
        (nlevels, z, theta, rt, u, v, ug, vg, p_in_Pa, subs, temperature_type,
          altitude_type, subs_type, retVars)

    Fortran declarations (read_sounding_file in src/Input_fields/sounding.F90):

        ! Input Variable(s)
        integer, intent(in) :: iunit ! I/O unit

        character(len=*), intent(in) :: &
          sounding_file ! Path to the configured input file

        real( kind = core_rknd ), intent(in) :: &
          p_sfc, & ! Pressure at the surface [Pa]
          zm_init ! Height at zm(1)         [m]

        integer, intent(in) :: &
          nmaxsnd   ! Maximum number of levels allowed in the sounding.in file

        integer, intent(in) :: &
          saturation_formula ! Integer that stores the saturation formula to be used

        ! Output Variable(s)
        integer, intent(out) :: nlevels ! Number of levels from the sounding.in file

        real( kind = core_rknd ), intent(out), dimension(nmaxsnd) :: &
          z,      & ! Altitude                               [m]
          theta,  & ! Liquid potential temperature sounding  [K]
          rt,     & ! Total water mixing ratio sounding      [kg/kg]
          u,      & ! u wind sounding                        [m/s]
          v,      & ! v wind sounding                        [m/s]
          ug,     & ! u geostrophic wind sounding            [m/s]
          vg,     & ! v geostrophic wind sounding            [m/s]
          p_in_Pa,& ! Pressure sounding                      [Pa]
          subs      ! Subsidence sounding                    [m/s or Pa/s]

        type(one_dim_read_var), intent(out), dimension(n_snd_var) :: &
          retVars ! Structure containing sounding profile

        character(len=*), intent(out) :: &
          temperature_type, &     ! Type of temperature sounding
          altitude_type, &       ! Type of independent coordinate
          subs_type         ! Type of subsidence
    """
    del nmaxsnd

    # ---- Begin Code ----

    retVars = read_one_dim_file(iunit, n_snd_var,
                                sounding_file)

    fill_blanks_one_dim_vars(n_snd_var, retVars)

    # The arrays below have one value per sounding level rather than nmaxsnd values.
    nlevels = retVars[0].values.size

    z, p_in_Pa, altitude_type = read_z_profile(n_snd_var, nlevels, retVars, p_sfc, zm_init,
                                               saturation_formula)

    temperature_type, theta = read_theta_profile(n_snd_var, nlevels, retVars)

    rt = read_x_profile(n_snd_var, nlevels, rt_name, retVars)

    u = read_x_profile(n_snd_var, nlevels, um_name, retVars)

    v = read_x_profile(n_snd_var, nlevels, vm_name, retVars)

    ug = read_x_profile(n_snd_var, nlevels, ug_name, retVars)

    vg = read_x_profile(n_snd_var, nlevels, vg_name, retVars)

    subs_type, subs = read_subs_profile(n_snd_var, nlevels, retVars)

    return (nlevels, z, theta, rt, u, v,
            ug, vg, p_in_Pa, subs,
            temperature_type, altitude_type, subs_type,
            retVars)


# -------------------------------------------------------------------------------------------------
def read_sclr_sounding_file(sclr_dim, sclr_idx, iunit, sclr_sounding_file, nmaxsnd,
                            sclr_snd):
    """This subroutine reads in a configured sclr_sounding.in file and
       returns the values contained in that file.

    Returns:
        (sclr_snd, retVars)

    Fortran declarations (read_sclr_sounding_file in src/Input_fields/sounding.F90):

        !--------------------- Input Variables ---------------------
        integer, intent(in) :: &
          sclr_dim

        type (sclr_idx_type), intent(in) :: &
          sclr_idx

        integer, intent(in) :: iunit ! I/O unit

        character(len=*), intent(in) :: sclr_sounding_file ! String identifying the model case;
        !                                         e.g. bomex

        integer, intent(in) :: &
          nmaxsnd   ! Maximum number of levels allowed in the sounding.in file

        !--------------------- InOut Variables ---------------------
        real( kind = core_rknd ), intent(inout), dimension(nmaxsnd,sclr_max) :: &
          sclr_snd        ! Scalar sounding [?]

        !--------------------- Output Variables --------------------
        type(one_dim_read_var), dimension(sclr_dim), intent(out) :: &
          retVars ! Structure containing scalar sounding
    """
    del nmaxsnd

    retVars = read_one_dim_file(iunit, sclr_dim,
                                sclr_sounding_file)

    for sclr in range(sclr_dim):
        name = retVars[sclr].name
        if name == CO2_name:
            if sclr + 1 != sclr_idx.iisclr_CO2 and sclr_idx.iisclr_CO2 > 0:
                raise ValueError("iisclr_CO2 index does not match column.")
        elif name == rt_name:
            if sclr + 1 != sclr_idx.iisclr_rt and sclr_idx.iisclr_rt > 0:
                raise ValueError("iisclr_rt index does not match column.")
        elif name in (theta_name, thetal_name, temperature_name):
            if sclr + 1 != sclr_idx.iisclr_thl and sclr_idx.iisclr_thl > 0:
                raise ValueError("iisclr_thl index does not match column.")
        sclr_snd[:retVars[sclr].values.size, sclr] = retVars[sclr].values

    return sclr_snd, retVars


# -------------------------------------------------------------------------------------------------
def read_edsclr_sounding_file(edsclr_dim, sclr_idx, iunit, edsclr_sounding_file, nmaxsnd,
                              edsclr_snd):
    """This subroutine reads in a configured edsclr_sounding.in file and
      returns the values contained in that file.

    Returns:
        edsclr_snd

    Fortran declarations (read_edsclr_sounding_file in src/Input_fields/sounding.F90):

        !--------------------- Input Variables ---------------------
        integer, intent(in) :: &
          edsclr_dim

        type (sclr_idx_type), intent(in) :: &
          sclr_idx

        integer, intent(in) :: iunit ! I/O unit

        character(len=*), intent(in) :: &
          edsclr_sounding_file ! Path to the configured input file

        integer, intent(in) :: &
          nmaxsnd   ! Maximum number of levels allowed in the sounding.in file

        !--------------------- Output Variables ---------------------
        real( kind = core_rknd ), intent(inout), dimension(nmaxsnd,sclr_max) :: &
        edsclr_snd ! Eddy Scalars [?]
    """
    del nmaxsnd

    # ---- Begin Code ----

    retVars = read_one_dim_file(iunit, edsclr_dim,
                                edsclr_sounding_file)

    for edsclr in range(edsclr_dim):

        name = retVars[edsclr].name
        if name == CO2_name:
            if edsclr + 1 != sclr_idx.iiedsclr_CO2 and sclr_idx.iiedsclr_CO2 > 0:
                raise ValueError("iisclr_CO2 index does not match column.")
        elif name == rt_name:
            if edsclr + 1 != sclr_idx.iiedsclr_rt and sclr_idx.iiedsclr_rt > 0:
                raise ValueError("iisclr_rt index does not match column.")
        elif name in (theta_name, thetal_name, temperature_name):
            if edsclr + 1 != sclr_idx.iiedsclr_thl and sclr_idx.iiedsclr_thl > 0:
                raise ValueError("iisclr_thl index does not match column.")
        edsclr_snd[:retVars[edsclr].values.size, edsclr] = retVars[edsclr].values

    deallocate_one_dim_vars(edsclr_dim, retVars)

    return edsclr_snd


# read_profile is commented out in sounding.F90, which has no callers, and is
# not ported.


# -------------------------------------------------------------------------------------------------
def convert_snd2extended_atm(iunit, ozone_sounding_file, n_snd_var, p_sfc, zm_init,
                             sounding_profiles, saturation_formula):
    """This subroutine converts information retrieved from the
        sounding files of a case into a format usable for an extended atmosphere.
        The extended atmosphere profile is stored in module variables.

    Fortran declarations (convert_snd2extended_atm in src/Input_fields/sounding.F90):

        ! Input Variable(s)

        integer, intent(in) :: iunit  ! Fortran file unit

        character(len=*), intent(in) ::  ozone_sounding_file ! Ozone profile path

        integer, intent(in) :: n_snd_var   ! Number of variables from sounding [-]

        real( kind = core_rknd ), intent(in) :: p_sfc  ! Pressure at the surface [Pa]

        real( kind = core_rknd ), intent(in) :: zm_init ! Height at zm(1) [m]

        type(one_dim_read_var), dimension(n_snd_var), intent(in) :: &
          sounding_profiles ! Sounding profile

        integer, intent(in) :: &
          saturation_formula ! Integer that stores the saturation formula to be used
    """
    # Constant Parameters
    n_rad_scalars = 1

    # -- Begin Code --

    # Determine the size of the extended atmosphere buffer
    extended_atmos_dim = sounding_profiles[0].values.size

    # Either convert to pressure or from pressure

    extended_alt, extended_p_in_Pa, altitude_type = read_z_profile(
        n_snd_var, extended_atmos_dim, sounding_profiles, p_sfc, zm_init,
        saturation_formula,
    )

    if altitude_type == z_name:
        raise ValueError("Fatal error in convert_snd2extended_atm. Feature not implemented")

    extended_p_in_mb = extended_p_in_Pa / pascal_per_mb

    # Convert to temperature from thlm or theta

    temperature_type, extended_T_in_K = read_theta_profile(n_snd_var, extended_atmos_dim,
                                                           sounding_profiles)

    if temperature_type != temperature_name:
        extended_exner = (extended_p_in_Pa / p0)**kappa
        extended_exner[0] = (p_sfc / p0)**kappa
        extended_T_in_K = extended_T_in_K * extended_exner

    # Convert rtm to specific humidity

    extended_sp_hmdty = read_x_profile(n_snd_var, extended_atmos_dim, rt_name,
                                       sounding_profiles)
    extended_sp_hmdty = extended_sp_hmdty / (extended_sp_hmdty + 1.0)

    # Read in radiation scalars sounding (currently it only holds O3)
    extended_rad_scalars = read_one_dim_file(iunit, n_rad_scalars,
                                             ozone_sounding_file)

    # Set the array holding the values of o3l (ozone)
    extended_o3l = read_x_profile(1, extended_atmos_dim, ozone_name,
                                  extended_rad_scalars)
    # We would add the setting of new radiation scalar arrays like O3 right
    # after this. New variable names would have to be added to input_names.

    return dict(alt=extended_alt, T_in_K=extended_T_in_K, sp_hmdty=extended_sp_hmdty,
                p_in_mb=extended_p_in_mb, o3l=extended_o3l)


# -------------------------------------------------------------------------------------------------
def load_extended_std_atm(iunit):
    """Loads in the U.S. Standard atmosphere data from a file.

    References:
    McClatchey, et al., (1972) _Environmental Research Papers_,
    No. 411, p.94

    Fortran declarations (load_extended_std_atm in src/Input_fields/sounding.F90):

        ! Input Variable(s)
        integer, intent(in) :: iunit ! File I/O unit [-]
    """
    # Constant Parameters
    nCol = 5  # Number of columns in the text file
    std_atmos_dim = 50  # Dimension of the standard atmosphere table

    atm_input_file = Path(__file__).resolve().parents[3] / "src" / "Radiation" / "BUGSrad" / "atmosphere.in"

    # -- Begin Code --

    extended_atmos_dim = std_atmos_dim

    retVars = read_one_dim_file(iunit, nCol, atm_input_file)

    extended_alt = read_x_profile(nCol, extended_atmos_dim, z_name, retVars,
                                  atm_input_file)

    extended_T_in_K = read_x_profile(nCol, extended_atmos_dim, temperature_name, retVars,
                                     atm_input_file)

    extended_sp_hmdty = read_x_profile(nCol, extended_atmos_dim, sp_humidity_name, retVars,
                                       atm_input_file)

    extended_p_in_mb = read_x_profile(nCol, extended_atmos_dim, press_mb_name, retVars,
                                      atm_input_file)

    extended_o3l = read_x_profile(nCol, extended_atmos_dim, ozone_name, retVars,
                                  atm_input_file)

    return dict(alt=extended_alt, T_in_K=extended_T_in_K, sp_hmdty=extended_sp_hmdty,
                p_in_mb=extended_p_in_mb, o3l=extended_o3l)
