"""CLUBB ↔ BUGSrad interface — port of bugsrad_driver.F90.

  * `compute_bugsrad_radiation`  ← bugsrad_driver.F90:compute_bugsrad_radiation (maps CLUBB state onto the
                                   top-down radiation grid — vertical flip + buffer + std-atm extension — calls
                                   bugs_rad per column, flips the heating rates back to the CLUBB grid).
  * `build_rad_grid_setup` / `_assemble_zt` — the grid-setup helpers this driver owns.

The std/extended-atmosphere readers it consumes live in their Fortran-home modules: `load_extended_std_atm` /
`convert_snd2extended_atm` / `read_ozone_sounding` in Input_fields/sounding.py (sounding.F90), and
`determine_extended_atmos_bounds` in Radiation/extended_atmosphere_module.py (extended_atmosphere_module.F90).

The radiation grid runs TOA→surface (index 1 = top), opposite to CLUBB's surface→top, hence the `[::-1]`
flips. `complete_momentum`/`complete_alt` from the Fortran are only used for optional diagnostic output and
are intentionally omitted. With the CLUBB build (-Dradoffline -Dnooverlap) bugs_rad runs on this grid directly.
"""
import numpy as np
import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

import functools
from clubb_jax.src.Radiation.BUGSrad.bugs_rad import bugs_rad
# The Fortran driver (bugsrad_driver.F90:357) passes constants_clubb grav/Cp to bugs_rad for the heating
# rate (heat_fac = grav·0.01/cp), NOT BUGSrad's own physconst (9.80665/1004.0) — a uniform 3.26e-4 factor.
from clubb_jax.src.CLUBB_core.constants_clubb import (
    grav as _GRAV_CLUBB, Cp as _CP_CLUBB,
    pascal_per_mb as PASCAL_PER_MB, g_per_kg as G_PER_KG, cloud_frac_min as CLOUD_FRAC_MIN,
)   # bugsrad_driver.F90:38-39 `use constants_clubb`


@functools.lru_cache(maxsize=None)
def _jitted_bugs_rad(s0, grav, cp, umco2, umch4, umn2o):
    """Cached jit of bugs_rad with the scalar gas/solar/gravity args bound (static). MUST be jitted: the
    eager 18-band×k-interval dispatch leaks ~700 MB/call (un-jitted device buffers accumulate → OOM after
    ~6 timesteps); jit compiles once per grid shape → bounded memory + ~100× faster, bit-exact to ~1e-13."""
    return jax.jit(functools.partial(bugs_rad, s0=s0, grav=grav, cp=cp,
                                     umco2=umco2, umch4=umch4, umn2o=umn2o))


# ── load_extended_std_atm / read_ozone_sounding / convert_snd2extended_atm now live in their Fortran-home
#    module Input_fields/sounding.py (sounding.F90:load_extended_std_atm / convert_snd2extended_atm, mirror-refactor
#    iter 215); radiation_module / clubb_driver / tests import them from there directly. ──


# determine_extended_atmos_bounds now lives in its Fortran-home module
# Radiation/extended_atmosphere_module.py (mirror-refactor iter 113), imported below.
from clubb_jax.src.Radiation.extended_atmosphere_module import determine_extended_atmos_bounds


def build_rad_grid_setup(zm_grid, zm_grid_spacing, p_in_Pa_zm, radiation_top, ext):
    """Precompute the static (state-independent) radiation-grid layout + reversed std-atm extension columns.
    Returns a dict consumed by compute_bugsrad_radiation: the 4 bounds, `buffer`=lin_int_buffer+range_size,
    and the extended T/sp_hmdty/o3l/p_in_mb already sliced top→bottom for radiation indices 1..range_size."""
    bottom, top, rsize, linbuf = determine_extended_atmos_bounds(
        zm_grid, zm_grid_spacing, p_in_Pa_zm, radiation_top, ext)
    buffer = linbuf + rsize
    # extended_X( top_level : bottom_level : -1 ) → 0-based slice [top-1 down to bottom-1], length rsize
    sl = slice(top - 1, bottom - 2 if bottom >= 2 else None, -1)
    return dict(bottom_level=bottom, top_level=top, range_size=rsize, lin_int_buffer=linbuf, buffer=buffer,
                ext_T=np.asarray(ext["T_in_K"])[sl].copy(),
                ext_q=np.asarray(ext["sp_hmdty"])[sl].copy(),
                ext_o3=np.asarray(ext["o3l"])[sl].copy(),
                ext_p_mb=np.asarray(ext["p_in_mb"])[sl].copy())


def _assemble_zt(model_rev, ext_col, setup):
    """Build a radiation-grid zt-type array (ncol, nzt+buffer): [extended | lin-interp | model(reversed)].
    `model_rev`=(ncol, nzt) the model-grid field already reversed (top→surface). `ext_col`=(range_size,)
    static extended values. Linear-interp buffer blends model-top (rad idx buffer+1) and extended-bottom
    (rad idx range_size), per the Fortran z-loop weights."""
    ncol = model_rev.shape[0]
    rsize, linbuf, buf = setup["range_size"], setup["lin_int_buffer"], setup["buffer"]
    ext_part = jnp.broadcast_to(jnp.asarray(ext_col)[None, :], (ncol, rsize))      # rad idx 1..range_size
    col_top = model_rev[:, 0]                                                      # rad idx buffer+1 (model top)
    col_ext = ext_part[:, rsize - 1]                                              # rad idx range_size (ext bottom)
    z1, z2 = buf + 1, rsize                                                        # Fortran indices of the endpoints
    if linbuf > 0:
        fz = np.arange(rsize + 1, buf + 1)                                         # 0-based idx range_size..buffer-1
        w1 = (z2 - (fz + 1)) / (z2 - z1)                                           # coeff of X(z1)=col_top; fz+1=Fortran idx
        w2 = ((fz + 1) - z1) / (z2 - z1)                                           # coeff of X(z2)=col_ext
        lin_part = col_top[:, None] * jnp.asarray(w1)[None, :] + col_ext[:, None] * jnp.asarray(w2)[None, :]
        return jnp.concatenate([ext_part, lin_part, model_rev], axis=1)
    return jnp.concatenate([ext_part, model_rev], axis=1)


def compute_bugsrad_radiation(setup, T_in_K, rcm, rtm, rsm, rim, cloud_frac, ice_supersat_frac,
                              p_in_Pa, p_in_Pam, rho_zm, exner, amu0, slr,
                              alvdr, alvdf, alndr, alndf, sol_const=1367.0):
    """Map CLUBB state → BUGSrad grid, run bugs_rad, return CLUBB-grid heating rate + fluxes.
    `setup`=build_rad_grid_setup(...). State arrays (ncol,nzt): T_in_K, rcm, rtm, rsm, rim, cloud_frac,
    ice_supersat_frac, p_in_Pa, exner. (ncol,nzm): p_in_Pam, rho_zm. Scalars: amu0, slr, surface albedos,
    sol_const. Returns dict(radht, radht_SW, radht_LW, Frad, Frad_SW_up/down, Frad_LW_up/down) on the CLUBB
    zt/zm grid (surface→top)."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    T_in_K, rcm, rtm, rsm, rim = a(T_in_K), a(rcm), a(rtm), a(rsm), a(rim)
    cloud_frac, ice_supersat_frac = a(cloud_frac), a(ice_supersat_frac)
    p_in_Pa, p_in_Pam, rho_zm, exner = a(p_in_Pa), a(p_in_Pam), a(rho_zm), a(exner)
    ncol, nzt = T_in_K.shape
    nzm = p_in_Pam.shape[1]
    buf, rsize = setup["buffer"], setup["range_size"]

    # --- model-grid derived fields (bugsrad_driver.F90:252-262) ---
    rcm_in_cloud = rcm / jnp.maximum(cloud_frac, CLOUD_FRAC_MIN)
    p_in_mb_model = p_in_Pa / PASCAL_PER_MB
    o3l_model = (5.4e-5 / rho_zm[:, :nzt]) / G_PER_KG
    sp_hum_model = jnp.where(rtm < rcm, 0.0, (rtm - rcm) / (1.0 + rtm))
    playerinmb_model = p_in_Pam / PASCAL_PER_MB                    # (ncol, nzm)

    # --- assemble the top-down radiation zt arrays (flip model + buffer + std-atm extension) ---
    rev = lambda X: X[:, ::-1]
    T_rad   = _assemble_zt(rev(T_in_K),       setup["ext_T"],  setup)
    q_rad   = _assemble_zt(rev(sp_hum_model), setup["ext_q"],  setup)
    o3_rad  = _assemble_zt(rev(o3l_model),    setup["ext_o3"], setup)
    p_rad   = _assemble_zt(rev(p_in_mb_model),setup["ext_p_mb"],setup)
    # cloud fields: extended + lin-interp region is cloud-free (zeros), only the model region is filled
    # (Fortran rcil/rsm_rad/rcm_in_cloud_rad/cloud_frac_rad init to 0, only [buffer+1:] set; the
    #  lin-interp z-loop touches only T/q/o3/p, never the cloud fields).
    def _assemble_cloud(model_rev):
        ext_lin = jnp.zeros((model_rev.shape[0], buf))      # extended + lin-interp region, all zero
        return jnp.concatenate([ext_lin, model_rev], axis=1)
    rcil_rad  = _assemble_cloud(rev(rim))                   # cloud ice  (qcil)
    rcmic_rad = _assemble_cloud(rev(rcm_in_cloud))          # cloud water (qcwl)
    cf_rad    = _assemble_cloud(rev(cloud_frac))            # cloud fraction (acld)

    # --- level-pressure (playerinmb) grid, size nzm+buffer (bugsrad_driver.F90:274-322) ---
    play_model_rev = playerinmb_model[:, ::-1]                     # reversed model levels → rad idx buffer+1..nzm+buffer
    # buffer levels = average of adjacent layer pressures p_rad: playerinmb(2:buffer+1)=(p_rad(1:buffer)+p_rad(2:buffer+1))/2
    play_buf = (p_rad[:, 0:buf] + p_rad[:, 1:buf + 1]) / 2.0       # → rad idx 2..buffer+1 (0-based 1..buffer)
    # index 1 (0-based 0): extrapolate from indices 2,3 (0-based 1,2 of the final grid = play_buf[:,0], play_buf[:,1])
    p2 = play_buf[:, 0]; p3 = play_buf[:, 1]
    play_0 = jnp.where(2.0 * p2 - p3 > 0.0, 2.0 * p2 - p3, 0.5 * p2)
    # final playerinmb: [play_0 | play_buf (len buffer) | play_model_rev without its first entry]
    # play_buf covers rad idx 2..buffer+1; play_model_rev covers buffer+1..nzm+buffer, so drop its first (overlap at buffer+1)
    playerinmb = jnp.concatenate([play_0[:, None], play_buf, play_model_rev[:, 1:]], axis=1)
    diff_pres = playerinmb[:, 1:] - playerinmb[:, :-1]             # (ncol, nzm-1+buffer) = (ncol, nzt+buffer)
    ts = T_rad[:, nzt + buf - 1]                                  # surface temp = bottom radiation level

    # --- call bugs_rad per column-batch (it is already vectorised over ncol; it does the
    #     den/rmix/cwrho conversion internally, so pass raw mixing ratios + cloud fraction) ---
    fn = _jitted_bugs_rad(float(sol_const), _GRAV_CLUBB, _CP_CLUBB, 330.0, 1.6, 0.28)
    out = fn(playerinmb, p_rad, diff_pres, T_rad, q_rad, rcmic_rad, rcil_rad, o3_rad,
             ts, amu0, slr, alvdf, alndf, alvdr, alndr, cf_rad)

    # --- flip heating rates back to the CLUBB grid (bugsrad_driver.F90:363-367) ---
    # radht_*_rad(:, nzt+buffer : buffer+1 : -1) → take the model region [buffer:nzt+buffer], reverse
    asl = out["asl"][:, buf:buf + nzt][:, ::-1]
    atl = out["atl"][:, buf:buf + nzt][:, ::-1]
    radht_SW = asl * (1.0 / exner)
    radht_LW = atl * (1.0 / exner)
    radht = radht_SW + radht_LW

    # --- fluxes on the zm grid (bugsrad_driver.F90:369-375); flux arrays have nzm+buffer levels ---
    # Frad_*(:, :) = Frad_*( nzm+buffer : buffer+1 : -1 ) → model levels [buffer:nzm+buffer], reversed
    def _flux_model(F):
        return F[:, buf:buf + nzm][:, ::-1]
    Frad_SW_up   = _flux_model(out["fusw"]); Frad_SW_down = _flux_model(out["fdsw"])
    Frad_LW_up   = _flux_model(out["fulw"]); Frad_LW_down = _flux_model(out["fdlw"])
    Frad_SW = Frad_SW_up - Frad_SW_down
    Frad_LW = Frad_LW_up - Frad_LW_down
    Frad = Frad_SW + Frad_LW

    return dict(radht=radht, radht_SW=radht_SW, radht_LW=radht_LW, Frad=Frad,
                Frad_SW_up=Frad_SW_up, Frad_SW_down=Frad_SW_down,
                Frad_LW_up=Frad_LW_up, Frad_LW_down=Frad_LW_down)
