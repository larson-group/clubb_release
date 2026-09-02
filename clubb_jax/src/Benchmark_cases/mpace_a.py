"""Contains subroutines for the M-PACE A intercomparison.

References:
    http://science.arm.gov/wg/cpm/scm/scmic5/index.html
"""

from __future__ import annotations

from pathlib import Path

import jax
import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.time_dependent_input import time_select
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv, Rd, g_per_kg, sec_per_hr
from clubb_jax.src.CLUBB_core.file_functions import file_read_1d, file_read_2d
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor, zlinterp_fnc


# These variables were moved here so that they could be
# accessible to all subroutines in mpace_a
# Joshua Fasching December 2007
file_ntimes = 139
file_nlevels = 38
per_line = 5

file_pressure = None
file_heights = None
file_times = None
dTdt_forcing = None
dqdt_forcing = None
vertT_forcing = None
vertq_forcing = None
um_obs = None
vm_obs = None
file_latent_ht = None
file_sens_ht = None


# ----------------------------------------------------------------------
def mpace_a_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, time, p_in_Pa):
    """References:
        Liou, Wallace and Hobbs, Shettle and Weinman
        http://science.arm.gov/wg/cpm/scm/scmic5/index.html
    """

    p_sfc = 101000.0

    # --------------------- Begin Code ---------------------
    # Use time_select to get the indexes before and after the specified time
    # and to get the ratio necessary for interpolation.
    before_time, after_time, ratio = time_select(time, file_ntimes, file_times)

    dTdt_column = linear_interp_factor(
        ratio,
        dTdt_forcing[:, after_time],
        dTdt_forcing[:, before_time],
    )
    dqdt_column = linear_interp_factor(
        ratio,
        dqdt_forcing[:, after_time],
        dqdt_forcing[:, before_time],
    )
    vertT_column = linear_interp_factor(
        ratio,
        vertT_forcing[:, after_time],
        vertT_forcing[:, before_time],
    )
    vertq_column = linear_interp_factor(
        ratio,
        vertq_forcing[:, after_time],
        vertq_forcing[:, before_time],
    )
    um_column = linear_interp_factor(
        ratio,
        um_obs[:, after_time],
        um_obs[:, before_time],
    )
    vm_column = linear_interp_factor(
        ratio,
        vm_obs[:, after_time],
        vm_obs[:, before_time],
    )

    # Do linear interpolation in space using zlinterp_fnc
    # TODO(port-mirror): vmap expresses the source loop over CLUBB columns.
    dTdt_hoc_grid = jax.vmap(zlinterp_fnc, in_axes=(0, None, None))(
        gr.zt, file_heights, dTdt_column
    )
    dqdt_hoc_grid = jax.vmap(zlinterp_fnc, in_axes=(0, None, None))(
        gr.zt, file_heights, dqdt_column
    )
    vertT_hoc_grid = jax.vmap(zlinterp_fnc, in_axes=(0, None, None))(
        gr.zt, file_heights, vertT_column
    )
    vertq_hoc_grid = jax.vmap(zlinterp_fnc, in_axes=(0, None, None))(
        gr.zt, file_heights, vertq_column
    )
    um_hoc_grid = jax.vmap(zlinterp_fnc, in_axes=(0, None, None))(
        gr.zt, file_heights, um_column
    )
    vm_hoc_grid = jax.vmap(zlinterp_fnc, in_axes=(0, None, None))(
        gr.zt, file_heights, vm_column
    )

    # Compute vertical motion
    # Michael Falk removed omega from this case.

    # Interpolation
    # no need to interpolate since wm_zt is set to 0 above

    # Boundary condition
    wm_zm = jnp.zeros((ngrdcol, gr.nzm), dtype=p_in_Pa.dtype)
    wm_zt = jnp.zeros((ngrdcol, gr.nzt), dtype=p_in_Pa.dtype)

    # Compute large-scale tendencies
    thlm_forcing = (
        (dTdt_hoc_grid + vertT_hoc_grid)
        * (p_sfc / p_in_Pa) ** (Rd / Cp)
        / sec_per_hr
    )
    rtm_forcing = (dqdt_hoc_grid + vertq_hoc_grid) / g_per_kg / sec_per_hr

    sclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, sclr_dim), dtype=p_in_Pa.dtype)
    edsclrm_forcing = jnp.zeros((ngrdcol, gr.nzt, edsclr_dim), dtype=p_in_Pa.dtype)

    # Test scalars with thetal and rt if desired
    if sclr_dim > 0:
        if sclr_idx.iisclr_thl > 0:
            sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_thl - 1].set(thlm_forcing)
        if sclr_idx.iisclr_rt > 0:
            sclrm_forcing = sclrm_forcing.at[:, :, sclr_idx.iisclr_rt - 1].set(rtm_forcing)
    if edsclr_dim > 0:
        if sclr_idx.iiedsclr_thl > 0:
            edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_thl - 1].set(
                thlm_forcing
            )
        if sclr_idx.iiedsclr_rt > 0:
            edsclrm_forcing = edsclrm_forcing.at[:, :, sclr_idx.iiedsclr_rt - 1].set(
                rtm_forcing
            )
    return (
        wm_zt,
        wm_zm,
        thlm_forcing,
        rtm_forcing,
        um_hoc_grid,
        vm_hoc_grid,
        sclrm_forcing,
        edsclrm_forcing,
    )


# ----------------------------------------------------------------------
def mpace_a_sfclyr(ngrdcol, time, rho_sfc):
    """Surface forcing subroutine for mpace_a case. Written
    October 2007 by Michael Falk.

    References:
        http://science.arm.gov/wg/cpm/scm/scmic5/index.html
    """
    # choose which times to use
    before_time, after_time, ratio = time_select(time, file_ntimes, file_times)

    latent_heat_flx = linear_interp_factor(
        ratio,
        file_latent_ht[after_time],
        file_latent_ht[before_time],
    )
    sensible_heat_flx = linear_interp_factor(
        ratio,
        file_sens_ht[after_time],
        file_sens_ht[before_time],
    )

    # Compute heat and moisture fluxes
    wpthlp_sfc = sensible_heat_flx / (rho_sfc * Cp)
    wprtp_sfc = latent_heat_flx / (rho_sfc * Lv)

    # Declare the value of ustar.
    ustar = jnp.full(ngrdcol, 0.25)
    return wpthlp_sfc, wprtp_sfc, ustar


# ----------------------------------------------------------------
def mpace_a_init(iunit, file_path):
    """This subroutine initializes the module by reading in forcing
    data used in the tndcy and sfclyr subroutines.
    """
    global file_pressure, file_heights, file_times
    global dTdt_forcing, dqdt_forcing, vertT_forcing, vertq_forcing
    global um_obs, vm_obs, file_latent_ht, file_sens_ht

    # ---- Begin Code ----
    # Path joining replaces the source's character concatenation at this host
    # file-I/O boundary.
    file_path = Path(file_path)
    file_pressure = jnp.asarray(
        file_read_1d(iunit, file_path / "mpace_a_press.dat", file_nlevels, per_line)
    )
    file_heights = jnp.asarray(
        file_read_1d(iunit, file_path / "mpace_a_heights.dat", file_nlevels, per_line)
    )
    file_times = jnp.asarray(
        file_read_1d(iunit, file_path / "mpace_a_times.dat", file_ntimes, per_line)
    )

    dTdt_forcing = jnp.asarray(
        file_read_2d(
            iunit,
            file_path / "mpace_a_dTdt.dat",
            file_nlevels,
            file_ntimes,
            per_line,
        )
    )
    dqdt_forcing = jnp.asarray(
        file_read_2d(
            iunit,
            file_path / "mpace_a_dqdt_horiz.dat",
            file_nlevels,
            file_ntimes,
            per_line,
        )
    )
    vertT_forcing = jnp.asarray(
        file_read_2d(
            iunit,
            file_path / "mpace_a_verts.dat",
            file_nlevels,
            file_ntimes,
            per_line,
        )
    )
    vertq_forcing = jnp.asarray(
        file_read_2d(
            iunit,
            file_path / "mpace_a_vertq.dat",
            file_nlevels,
            file_ntimes,
            per_line,
        )
    )
    um_obs = jnp.asarray(
        file_read_2d(
            iunit,
            file_path / "mpace_a_um_obs.dat",
            file_nlevels,
            file_ntimes,
            per_line,
        )
    )
    vm_obs = jnp.asarray(
        file_read_2d(
            iunit,
            file_path / "mpace_a_vm_obs.dat",
            file_nlevels,
            file_ntimes,
            per_line,
        )
    )
    file_latent_ht = jnp.asarray(
        file_read_1d(iunit, file_path / "mpace_a_lh.dat", file_ntimes, per_line)
    )
    file_sens_ht = jnp.asarray(
        file_read_1d(iunit, file_path / "mpace_a_sh.dat", file_ntimes, per_line)
    )


__all__ = ["mpace_a_tndcy", "mpace_a_sfclyr", "mpace_a_init"]
